/************************************************************************
 ************************************************************************
 **                                                                    **
 **           Finite Mixture of Regressions Lasso CPP                  **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

// ---------------------------------------------------------------------
// Author: Luke Lloyd-Jones
// Date started:      04/09/2015
// Date last updated: 04/02/2016
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Includes and headers
// ---------------------------------------------------------------------

#include <cmath>
#include <cassert>
#include "fmr_lasso.hpp"

using namespace boost::math;
using namespace arma;
using namespace std;

// ---------------------------------------------------------------------
// Class methods
// ---------------------------------------------------------------------


// Constructor for matrix of a given size
// Allocates memory, and initialises entries
// to zero

FMRLasso::FMRLasso(double lambda,
                   double criterion,
                   double epsilon,
                   unsigned int N,
                   unsigned int P,
                   unsigned int G)
{
    assert(N > 0);
    assert(P > 0);
    assert(lambda >= 0);
    mlambda = lambda;
    mcriterion = criterion;
    mepsilon = epsilon;
    mN = N;
    mP = P;
    mG = G;
}

// Overwritten destructor

FMRLasso::~FMRLasso()
{
}

// Calculate lasso estimates of regression Betas by
// coordinate descent lasso
mat FMRLasso::FMRLassoRun(const mat& geno,
                          const mat& pheno,
                          double sig1,
                          double sig2,
                          double mu_1,
                          double mu_2,
                          double pi1,
                          double pi2,
                          double lambda_in_1,
                          double lambda_in_2,
                          unsigned int maxit,
                          unsigned int mode,
                          const mat& beta_str,
                          const mat& beta_ind,
                          double* bic,
                          int no_active,
                          int verbose,
                          string outpath)
{
    assert(geno.n_rows == mN);
    assert(geno.n_cols == mP);
    // Set the old obj for the obj criterion
    double old_obj = -10e20;
    // cout << "Lasso constraint parameter is " << mlambda << endl;
    // Initialise some elements
    maxit = maxit;
    mat beta_old(mP, mG, fill::zeros);
    mat beta = beta_str;
    // cout << X << endl;
    // cout << Y << endl;
    // Assume we have read in the genotype and phenotype matrices
    // Call them X and Y in the code
    // ---------------------------------------------------
    // Initialise the beta vector, lambda (lasso penalty),
    // alpha (intercept), variances, and pis
    // ---------------------------------------------------
    mat resid(mN, mG, fill::zeros);
    vec lambda;
    // cout << "Betas after initialisation" << endl << beta << endl;
    if (mode == 0)
    {
        cout << "You are in lasso mode" << endl;
        lambda = {lambda_in_1, lambda_in_2};
    }
    else
    {
        cout << "You are in non lasso mode" << endl;
        lambda = {0.0, 0.0};
    }
    vec alpha  = {mu_1, mu_2};
    vec sig_2  = {sig1, sig2};
    vec sig_2_old = {sig1, sig2};
    vec pi     = {pi1, pi2};
    // cout << "Lasso vector is " << lambda << endl;
    // ---------------------------------------------------
    // Printing some elements to see how we are going
    // ---------------------------------------------------
    if (verbose == 1)
    {
        cout << "Lasso vector is" << endl << lambda  << endl;
        cout << "Mean  vector is" << endl << alpha   << endl;
        cout << "Sigma vector is" << endl << sig_2   << endl;
        cout << "Pi    vector is" << endl << pi      << endl;
        // cout << beta << endl;
        cout << "Size of X is " << geno.n_rows  << " by " << geno.n_cols  << endl;
        cout << "Size of Y is " << pheno.n_rows << " by " << pheno.n_cols << endl;
        cout << "Size of X is "    << mN << " by " << mP << endl;
        cout << "Size of beta is " << mP << " by " << mG << endl;
    }
    // ---------------------------------------------------
    // Initialise the objective function
    // ---------------------------------------------------
    double obj = 0.0;
    double mean = 0.0;
    double stdev = 0.0;
    for (unsigned int i = 0; i < mN; i++)
    {
        double inner = 0.0;
        for (unsigned int k = 0; k < mG; k++)
        {
            // cout << beta.col(k) << endl;
            mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
            stdev = sqrt(sig_2(k));
            // cout << "sigma inside initial objective " << stdev << endl;
            boost::math::normal_distribution<double>myNormal(mean, stdev);
            //cout << "Mean: "  << myNormal.mean() << endl;
            //cout << "Stdev: " << myNormal.standard_deviation() << endl;
            //cout << "Yi: "    << pheno(i, 0) << endl;
            //cout << "pdf: "   << boost::math::pdf(myNormal, pheno(i, 0)) << endl;
            inner += pi(k) * boost::math::pdf(myNormal, pheno(i, 0));
            //cout << "inner "  << inner << endl;
        }
        // obj += (1 / (double)(mN)) * log(inner);
        obj += log(inner);
        //cout << "inner obj "  << obj << endl;
    }
    if (verbose == 1)
    {
        cout << "Initial objective 1 "  <<  obj << endl;
    }
    // ---------------------------------------------------
    // Subtract the lasso penalisation
    // ---------------------------------------------------
    if (mode == 0)
    {
      for (unsigned int k = 0; k < mG; k++)
      {
          // cout << "abs beta" << sum(abs(beta.row(k))) << endl;
          obj -= pi(k) * lambda(k) * sum(abs(beta.col(k)));
      }
    }
    else
    {
        // If we are not in lasso mode
        obj = obj;
    }
    if (verbose == 1)
    {
        cout << "Initial objective 2 "  <<  obj << endl;
    }
    // ---------------------------------------------------
    // Initialise the taus and the weights
    // ---------------------------------------------------
    mat tau(mG, mN, fill::ones);
    tau = (1 / mG) * tau;
    mat weight(mP, mG, fill::ones);
    //cout << "tau "  <<  tau << endl;
    //cout << "weight "  <<  weight << endl;
    mat geno_sqr = geno % geno;
    // ---------------------------------------------------
    // Main loop
    // ---------------------------------------------------
     for (unsigned int iteration = 0; iteration < maxit; iteration++)
     {
        if (iteration == (maxit - 1))
        {
            cout << "You have reached maxim iterations and have failed to converge" << endl;
            cout << "Your beta estimates at this stage" << "\n" << beta << endl;
            exit(1);
        }
        // Update the tau scores
        // ---------------------
        for (unsigned int i = 0; i < mN; i++)
        {
            double inner = 0.0;
            for (unsigned int k = 0; k < mG; k++)
            {
                // cout << beta(k, i) << endl;
                mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
                stdev = sqrt(sig_2(k));
                //cout << "sigma inside tau 11 " << stdev << endl;
                boost::math::normal_distribution<double>myNormal(mean, stdev);
                inner += pi(k) * boost::math::pdf(myNormal, pheno(i));
            }
            for (unsigned int k = 0; k < mG; k++)
            {
                mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
                stdev = sqrt(sig_2(k));
                //cout << "sigma inside tau 12 " << stdev << endl;
                boost::math::normal_distribution<double>myNormal(mean, stdev);
//                if (inner == 0.0)
//                {
//                    tau(k, i) = pi(k);
//                } else
//                {
                    tau(k, i) = pi(k) * boost::math::pdf(myNormal, pheno(i)) / inner;
//                }
            }
        }
        // cout << "inner inside "  <<  inner << endl;
        // cout << "tau inside "  << endl << tau << endl;
        // Calculate the roots of the Lagrangian
        // -------------------------------------
        vec rho(mG);
        vec psi(mG);
        double gamma_c1;
        double gamma_c2;
        double gamma_1;
        double gamma_2;
        for (unsigned int k = 0; k < mG; k++)
        {
            rho(k) = sum(tau.row(k));
            psi(k) = lambda(k) * sum(abs(beta.col(k)));
        }
        //cout << "rho "  <<  rho << endl;
        //cout << "psi "  <<  psi << endl;
        gamma_c1 = 0.5 * (sum(rho) + sum(psi));
        gamma_c2 = 0.5 * sqrt(sum(rho % rho) + sum(psi % psi) +
                              2 * rho(0) * rho(1) -
                              2 * rho(0) * psi(1) +
                              2 * rho(0) * psi(0) +
                              2 * rho(1) * psi(1) -
                              2 * rho(1) * psi(0) -
                              2 * psi(0) * psi(1));
        gamma_1 = gamma_c1 + gamma_c2;
        gamma_2 = gamma_c1 - gamma_c2;
        //cout << "gamma 1 " << gamma_1 << endl;
        //cout << "gamma 2 " << gamma_2 << endl;
        // Update the pis with the correct gamma solution
        // i.e., the one that gives the solution between (0, 1)
        // ----------------------------------------------------
        vec pi_1(mG);
        vec pi_2(mG);
        for (unsigned int k = 0; k < mG; k++)
        {
            pi_1(k) = rho(k) / (gamma_1 - psi(k));
            pi_2(k) = rho(k) / (gamma_2 - psi(k));
        }
        // cout << "pi1 and pi2 "  << pi_1 << "\t" << pi_2 << endl;
        if ((max(pi_1) < 1) * (min(pi_1) > 0) == 1)
        {
            pi = pi_1;
        } else
        {
            pi = pi_2;
        }
        // cout << "pi chosen " << pi << endl;
        // Update the intercepts
        // ---------------------
        for (unsigned int k = 0; k < mG; k++)
        {
            double top = 0.0;
            double bottom = sum(tau.row(k));
            for (unsigned int i = 0; i < mN; i++)
            {
                top += tau(k, i) * (pheno(i) - sum(beta.col(k) % trans(geno.row(i))));
            }
            alpha(k) = top / bottom;
            
        }
        //cout << "alphas" << alpha << endl;
        // Update the variances
        // ---------------------
        sig_2_old = sig_2;
        for (unsigned int k = 0; k < mG; k++)
        {
            double top = 0.0;
            double bottom = sum(tau.row(k));
            for (unsigned int i = 0; i < mN; i++)
            {
                top += tau(k, i) * pow((pheno(i) - alpha(k) - sum(beta.col(k) % trans(geno.row(i)))), 2);
            }
            // Experimental with the degredation of sigma
            //sig_2[k] = 0.1 * max((top / bottom), 0.01) + 0.8 * sig_2_old[k];
            sig_2[k] = max((top / bottom), 0.1);
            //cout << "sigma bottom" << bottom << endl;
            //cout << "sigma top"    << top << endl;
        }
        
        //cout << "sigmas " << sig_2  << endl;
        // Update the taus given pis and sigmas
        // ------------------------------------
        for (unsigned int i = 0; i < mN; i++)
        {
            double inner = 0.0;
            for (unsigned int k = 0; k < mG; k++)
            {
                // cout << beta(k, i) << endl;
                mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
                stdev = sqrt(sig_2(k));
                //cout << "sigma inside tau 12 " << stdev << endl;
                boost::math::normal_distribution<double>myNormal(mean, stdev);
                inner += pi(k) * boost::math::pdf(myNormal, pheno(i));
            }
            for (unsigned int k = 0; k < mG; k++)
            {
                mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
                stdev = sqrt(sig_2(k));
                //cout << "sigma inside tau 22 " << stdev << endl;
                boost::math::normal_distribution<double>myNormal(mean, stdev);
//                if (inner == 0.0)
//                {
//                    tau(k, i) = pi(k);
//                } else
//                {
                    tau(k, i) = pi(k) * boost::math::pdf(myNormal, pheno(i)) / inner;
//                }
                //cout << "inner inside "  <<  inner << endl;
            }
        }
        // cout << "inner inside "  <<  inner << endl;
        // cout << "tau inside "  <<  tau << endl;
        // Update the weights
        // ------------------
        for (unsigned int k = 0; k < mG; k++)
        {
            for (unsigned int j = 0; j < mP; j++)
            {
                weight(j, k) = 1 / (abs(beta(j, k)) + mepsilon);
            }
            
        }
        // cout << "weights updated" << weight << endl;
        // Update the betas. Most critical step
        // ------------------------------------
        beta_old = beta;
        for (unsigned int k = 0; k < mG; k++)
        {
            resid.col(k) = pheno - geno * beta_old.col(k) - alpha(k);
            for (unsigned int j = 0; j < mP; j++)
            {
                if ((mode != 0) & (beta_ind(j, k) == 1))
                {
                    beta(j, k) = 0.0;
                }
                else
                {
                double c1 = 0.0;
                double c2 = 0.0;
                double c3 = 0.0;
                c1 = mP * as_scalar(tau.row(k) *  (geno.col(j) % geno.col(j)));
                c2 = sum(trans(tau.row(k)) % geno.col(j) % resid.col(k));
                c3 = sig_2(k) * pi(k) * lambda(k) * weight(j, k);
                beta(j, k) = (beta_old(j, k) * c1 + c2) / (c1 + c3);
                }
            }
        }
        // cout << "beta updated " << beta << endl;
        // Update the objective function based on all new updates
        // ------------------------------------------------------
        obj = 0.0;
        for (unsigned int i = 0; i < mN; i++)
        {
            double inner = 0.0;
            for (unsigned int k = 0; k < mG; k++)
            {
                // cout << beta(k, i) << endl;
                mean = alpha(k) + sum(beta.col(k) % trans(geno.row(i)));
                stdev = sqrt(sig_2(k));
                // cout << "sigma inside likelihood " << stdev << endl;
                boost::math::normal_distribution<double>myNormal(mean, stdev);
                inner += pi(k) * boost::math::pdf(myNormal, pheno(i));
            }
            obj +=  log(inner);
        }
        // ---------------------------------------------------
        // Subtract the lasso penalisation
        // ---------------------------------------------------
        for (unsigned int k = 0; k < mG; k++)
        {
              // cout << "abs beta" << sum(abs(beta.row(k))) << endl;
              obj -= pi(k) * lambda(k) * sum(abs(beta.col(k)));
        }
        // ---------------------------------------------------
        // Manage the objective function
        // ---------------------------------------------------
        // Check for descent failure or convergence. If neither occurs,
        // record the new value of the objective function
         if ((iteration % 100 == 0))
         {
            if (verbose == 1)
            {
            cout << "Objective difference at iteration " << iteration <<
                    "\t" << (obj - old_obj) << endl;
            }
            if (mode != 0)
            {
                // ---------------------------------------------------
                // Calculate the BIC and print to std output
                // ---------------------------------------------------
                *bic = -2 * obj + (3 * mG + no_active) * log(mN);
                cout << *bic << endl;
            }
         }
        // cout << "Old objective " << old_obj << endl;
        // cout << "Objective "     << obj << endl;
        // cout << "Objective difference " << obj - old_obj << endl;
        if ((obj - old_obj) < 0)
        {
            if (verbose == 1)
            {
                cout << "Old objective " << old_obj << endl;
                cout << "Objective "     << obj << endl;
                cout << "Betas"  << endl << beta  << endl;
                cout << "Sigmas" << endl << sig_2  << endl;
                cout << "Alphas" << endl << alpha  << endl;
                cout << "Pis"    << endl << pi     << endl;
            }
            cout << "Objective difference " << obj - old_obj << endl;
            cout << "*** ERROR *** OBJECTIVE FUNCTION INCREASE" << endl;
            //exit(1);
            break;
        }
        if (abs(obj - old_obj) < mcriterion)
        {
            if (verbose == 1)
            {
                // cout << estimate << endl;
                cout << "We have convergence" << endl;
                cout << "Objective at convergence " << obj << endl;
            }
            // Break from the function
            break;
            
        }
        else
        {
            old_obj = obj;
        }

    } // For loop finished
    if (verbose == 1)
    {
        // Print to screen and write to outpath the results
        string beta_out;
        string sigma_out;
        string alpha_out;
        string pi_out;
        // Sigmas
        cout << "Betas"  << endl << beta  << endl;
        beta_out = outpath + "beta_estimates.txt";
        beta.save(beta_out, csv_ascii);
        // Sigmas
        cout << "Sigmas" << endl << sig_2  << endl;
        sigma_out = outpath + "sigma_estimates.txt";
        sig_2.save(sigma_out, csv_ascii);
        // Alphas
        cout << "Alphas" << endl << alpha  << endl;
        alpha_out = outpath + "alpha_estimates.txt";
        alpha.save(alpha_out, csv_ascii);
        // Pis
        cout << "Pis"    << endl << pi    << endl;
        pi_out = outpath + "pi_estimates.txt";
        pi.save(pi_out, csv_ascii);
    }
    //cout << "Estimate in the loop are " << estimate << endl;
    return beta;
}

// Standardise the genotypes to have mean 0 and sd 1 by allele frequencies

mat FMRLasso::GenoStd(mat& geno)
{
    assert(geno.n_rows == mN);
    assert(geno.n_cols == mP);
    
    double q;
    for ( unsigned int j = 0; j < mP; j++)
    {
        q = sum(geno.col(j)) / (2.0 * mN); // BE CAREFUL FOR A VECTOR 0 INDEX
        //cout << q << endl;
        geno.col(j) = (geno.col(j) - 2.0 * q) / sqrt(2.0 * q * (1.0 - q));
    }
    
    return geno;
}

// Function to run the lasso run function for a given lambda
// and return the BIC criterions

double FMRLasso::BICRun(const mat& geno,
                        const mat& pheno,
                        double sig_1,
                        double sig_2,
                        double mu_1,
                        double mu_2,
                        double pi_1,
                        double pi_2,
                        double lambda_in_1,
                        double lambda_in_2,
                        double dists,
                        unsigned int maxit,
                        int verbose,
                        string outpath)
{
    // Declare the est matrix to be returned
    mat est;
    double bic = 0.0;
    int no_not_active = 0;
    int no_active = geno.n_cols * dists;
    // Decalre the initial values for beta to be passed to the function
    // as warm sarts. Initially all betas are 0.01
    mat beta_str(geno.n_cols, dists, fill::ones);
    // Initiate betas by a read in from the outpath
    string beta_str_in;
    beta_str_in = outpath + "betas_str.txt";
    beta_str.load(beta_str_in, csv_ascii);
    // cout << "Loaded betas" << endl << beta_str << endl;
    mat est_ind(geno.n_cols, dists, fill::zeros);
    // cout << sig_1 << endl;
    // ------------------------------------------
    // Run the penalised lasso
    // ------------------------------------------
    est = FMRLassoRun(geno, pheno, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2,
                      lambda_in_1, lambda_in_2, maxit, 0, beta_str, est_ind,
                      &bic, no_active, 0, outpath);
    // Find which values of est that are greater than a threshold
    ActiveSet(est, est_ind, &no_not_active);
    // Print some items for checking
    // cout << "Estimates " << est << endl;
    // cout << "Estimates indicator matrix " << est_ind << endl;
    // Caclulate the number of active elements to
    no_active = no_active - no_not_active;
    // ------------------------------------------
    // Run the unpenalised FMR
    // ------------------------------------------
    if (verbose == 1)
    {
        est = FMRLassoRun(geno, pheno, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2,
                          lambda_in_1, lambda_in_2, maxit, 1, est, est_ind,
                          &bic, no_active, 1, outpath);
    } else
    {
        FMRLassoRun(geno, pheno, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2,
                    lambda_in_1, lambda_in_2, maxit, 1, est, est_ind,
                    &bic, no_active, 0, outpath);
    }
    return bic;
}

// Function to set elements of the parameter matrix to 0 and alter
// an index matrix with 1s for those elements that have been set to
// zero. Also alters the no of non active elements

void FMRLasso::ActiveSet(mat& est, mat& est_ind, int* no_not_active)
{
    double delta = 1e-10;
    for (unsigned int i = 0; i < est.n_rows; i++)
    {
        for (unsigned int j = 0; j < est.n_cols; j++)
        {
            if (abs(est(i, j)) < delta)
            {
                est(i, j) = 0;
                est_ind(i, j) = 1;
                *no_not_active += 1;
            }
        }
    }
}


//// Nelder and Mead's simplex function
//// ----------------------------------
//
std::vector<double> FMRLasso::Simplex(std::vector<double> init,    //initial guess of the parameters
                                 double tol, //termination criteria
                                 int iterations,
                                 const mat& geno,
                                 const mat& pheno,
                                 double sig_1,
                                 double sig_2,
                                 double mu_1,
                                 double mu_2,
                                 double pi_1,
                                 double pi_2,
                                 double lambda_in_1,
                                 double lambda_in_2,
                                 double dists,
                                 unsigned int maxit,
                                 int verbose,
                                 string outpath)
{
    //iteration step number
    unsigned int N = init.size();                         //space dimension
    const double a = 1.0, b = 1.0, g = 0.5, h = 0.5;   //coefficients
    std::vector<std::vector<double> > x =  std::vector<std::vector<double> >(); //x: The Simplex
    //a: reflection  -> xr
    //b: expansion   -> xe
    //g: contraction -> xc
    //h: full contraction to x1
    std::vector<double> xcentroid_old(N, 0);   //simplex center * (N+1)
    std::vector<double> xcentroid_new(N, 0);   //simplex center * (N+1)
    std::vector<double> vf(N + 1,0);           //f evaluated at simplex vertexes
    unsigned int x1 = 0, xn = 0, xnp1 = 0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
    // xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
    //   xn: f(xn) < f(xnp1) && f(xn) > all other f(x_i)
    int cnt = 0; //iteration step number
    
    if(x.size() == 0) //if no initial simplex is specified
    {   // construct the trial simplex
        // based upon the initial guess parameters
        std::vector<double> del( init );
        std::transform(del.begin(), del.end(), del.begin(),
                       std::bind2nd( std::divides<double>() , 20) );//'20' is picked
        //assuming initial trial close to true
        
        for(unsigned int i = 0; i < N; ++i){
            std::vector<double> tmp( init );
            tmp[i] +=  del[i];
            //cout << "The constrcuted simplex " << tmp[i] << endl;
            x.push_back( tmp );
        }
        x.push_back(init);//x.size()=N+1, x[i].size()=N
        
        //xcentriod
        std::transform(init.begin(), init.end(),
                       xcentroid_old.begin(), std::bind2nd(std::multiplies<double>(), N + 1) );
    }//constructing the simplex finished
    cout << "Initial simplex is " << endl;
    cout << x[0].at(0) << endl;
    cout << x[0].at(1) << endl;
    cout << x[1].at(0) << endl;
    cout << x[1].at(1) << endl;
    cout << x[2].at(0) << endl;
    cout << x[2].at(1) << endl;
    //optimization begins
    for (cnt = 0; cnt < iterations; ++cnt)
    {
        for (unsigned int i = 0; i < N + 1; ++i)
        {
            // vf[i] = f(x[i]);
            vf[i] = BICRun(geno,
                           pheno,
                           sig_1,
                           sig_2,
                           mu_1,
                           mu_2,
                           pi_1,
                           pi_2,
                           x[i].at(0),
                           x[i].at(1),
                           dists,
                           maxit,
                           verbose,
                           outpath);
        }
        // cout << "poo" << vf.at(0) << endl;
        x1 = 0; xn = 0; xnp1 = 0; //find index of max, second max, min of vf.
        for(unsigned int i = 0; i < vf.size(); ++i)
        {
            if(vf[i] < vf[x1]){
                x1 = i;
            }
            if(vf[i] > vf[xnp1]){
                xnp1 = i;
            }
        }
        // cout << xnp1 << endl;
        xn = x1;
        for(unsigned int i = 0; i < vf.size(); ++i)
        {
            if(vf[i] < vf[xnp1] && vf[i] > vf[xn])
            {
                xn = i;
            }
        }
        //x1, xn, xnp1 are found
        // cout << "The xs are " << x1 << xn << xnp1 << endl;
        std::vector<double> xg(N, 0);//xg: centroid of the N best vertexes
        for(unsigned int i = 0; i < x.size(); ++i){
            if(i != xnp1)
            {
                std::transform(xg.begin(), xg.end(), x[i].begin(), xg.begin(), std::plus<double>() );
            }
        }
        std::transform(xg.begin(),
                       xg.end(),
                       x[xnp1].begin(),
                       xcentroid_new.begin(),
                       std::plus<double>());
        std::transform(xg.begin(),
                       xg.end(),
                       xg.begin(),
                       std::bind2nd(std::divides<double>(), N) );
        //xg found, xcentroid_new updated
        //cout << "xg found" << xg.at(1) << endl;
        
        //termination condition
        double diff = 0;          //calculate the difference of the simplex centers
        //see if the difference is less than the termination criteria
        for(unsigned int i = 0; i < N; ++i)
        {
            diff += fabs(xcentroid_old[i]-xcentroid_new[i]);
        }
        if (diff / N < tol)
        {
            break;              //terminate the optimizer
        }
        else
        {
            xcentroid_old.swap(xcentroid_new); //update simplex center
        }
        cout << "Difference in objective at current iteration is " << diff << endl;
        //reflection:
        std::vector<double> xr(N, 0);
        for(unsigned int i = 0; i < N; ++i)
        {
            xr[i] = xg[i] + a * (xg[i] - x[xnp1][i]);
        }
        //reflection, xr found
        
        double fxr = BICRun(geno,
                            pheno,
                            sig_1,
                            sig_2,
                            mu_1,
                            mu_2,
                            pi_1,
                            pi_2,
                            xr.at(0),
                            xr.at(1),
                            dists,
                            maxit,
                            verbose,
                            outpath);//record function at xr
        // cout << "The fxr is? " << fxr << endl;
        if(vf[x1] <= fxr && fxr <= vf[xn])
        {
            std::copy(xr.begin(), xr.end(), x[xnp1].begin());
        }
        else if (fxr < vf[x1]) //expansion:
        {
            std::vector<double> xe(N, 0);
            for (unsigned int i = 0; i < N; ++i)
            {
                xe[i] = xr[i] + b * (xr[i] - xg[i]);
            }
            double fxe = BICRun(geno,
                                pheno,
                                sig_1,
                                sig_2,
                                mu_1,
                                mu_2,
                                pi_1,
                                pi_2,
                                xe.at(0),
                                xe.at(1),
                                dists,
                                maxit,
                                verbose,
                                outpath);
            // cout << "XE1" << xe.at(0) << " " << "XE2" << xe.at(1) << endl;
            // cout << "fxe" << fxe << endl;
            if (fxe < fxr)
            {
                std::copy(xe.begin(), xe.end(), x[xnp1].begin());
            }
            else
            {
                std::copy(xr.begin(), xr.end(), x[xnp1].begin());
            }
        }//expansion finished,  xe is not used outside the scope
        else if (fxr > vf[xn]) //contraction:
        {
            std::vector<double> xc(N, 0);
            for (unsigned int i = 0; i < N; ++i)
            {
                xc[i] = xg[i] + g * (x[xnp1][i] - xg[i]);
            }
            double fxc = BICRun(geno,
                                pheno,
                                sig_1,
                                sig_2,
                                mu_1,
                                mu_2,
                                pi_1,
                                pi_2,
                                xc.at(0),
                                xc.at(1),
                                dists,
                                maxit,
                                verbose,
                                outpath);
            if (fxc < vf[xnp1])
            {
                std::copy(xc.begin(), xc.end(), x[xnp1].begin());
            }
            else
            {
                for (unsigned int i = 0; i < x.size(); ++i)
                {
                    if (i != x1)
                    {
                        for(unsigned int j = 0; j < N; ++j)
                        {
                            x[i][j] = x[x1][j] + h * (x[i][j] - x[x1][j]);
                        }
                    }
                }
                
            }
        }//contraction finished, xc is not used outside the scope

    //optimization is finished

    if (cnt == (iterations - 1))
    {//max number of iteration achieves before tol is satisfied
        std::cout << "Iteration limit achieves, result may not be optimal" << std::endl;
    }
        cout << "COMPLETED RUNNING OF NELDER-MEAD SIMPLEX ALGOIRTHM " << cnt << endl;
    }
    return x[x1];
}




















