/************************************************************************
 ************************************************************************
 **                                                                    **
 **           Finite Mixture of Regressions Lasso Main                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

// ---------------------------------------------------------------------
// Author: Luke Lloyd-Jones
// Date started:   13/06/2015
// Date last updated: 10/06/2016
// ---------------------------------------------------------------------

// Includes

#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <fstream>
#include <sstream>
#include <ostream>
#include <array>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <iomanip>  // setprecision()
#include <string>
#include <vector>
#include <cerrno>
#include <armadillo>
#include <time.h>
#include "fmr_lasso.cpp"
// #include "simplex.hpp"
#include <typeinfo>
#include <boost/math/distributions/normal.hpp>

// Set the namespace std

using namespace std;
using namespace arma;
using namespace boost::math;

// MAIN COMPUTATION BLOCK
// ----------------------
// ----------------------
// ----------------------

int main(int argc, char* argv[])
{
    // ----------------------------------------------------------
    // Print the arguments that have been passed
    // ----------------------------------------------------------
    std::cout << "Number of command line arguments " << argc << "\n";
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Argument " << i << " = " << argv[i];
        std::cout << "\n";
    }
    // ----------------------------------------------------------
    // Manage the argument passing
    // ----------------------------------------------------------
    double lasso_b1    = atof(argv[3]);
    double lasso_b2    = atof(argv[4]);
    double con_crit    = atof(argv[5]);
    double maxit       = atoi(argv[6]);
    double pterb       = atof(argv[7]);
    static const std::string outpath  = argv[8];
    static const std::string beta_str_in  = argv[9];
    unsigned int dists = (argc - 9) / 3;
    cout << "Number of dists " << dists << " seconds" << endl;
    // Fill the initial vectors for the mus, pis, and sigs
    vec sig_srt(dists);
    vec mu_srt(dists);
    vec pi_srt(dists);
    for (unsigned int i = 0; i < dists; i++)
    {
        sig_srt(i) = atof(argv[(10 + i)]);
        mu_srt(i)  = atof(argv[(10 + dists + i)]);
        pi_srt(i)  = atof(argv[(10 + 2 * dists + i)]);
    }
//     cout << "sig start " << sig_srt << endl;
//     cout << "mu  start " << mu_srt << endl;
//     cout << "pi  start " << pi_srt << endl;
    clock_t t1, t2, t3, t4;
    // ----------------------------------------------------------
    // Reading data in .csv format to eliminate possible data
    // read problems
    // ----------------------------------------------------------
    t1 = clock();
    mat X;
    mat Y;
    X.load(argv[1], csv_ascii);
    Y.load(argv[2], csv_ascii);
    t2 = clock();
    float diff1 ((float)t2 - (float)t1);
    float seconds1 = diff1 / CLOCKS_PER_SEC;
    cout << "File input took " << seconds1 << " seconds" << endl;
    // ----------------------------------------------------------
    // Load up the starting betas
    // ----------------------------------------------------------
    mat beta_str(X.n_cols, dists, fill::ones);
    // Initiate betas by a read in from the outpath
    beta_str.load(beta_str_in, csv_ascii);
    // ----------------------------------------------------------
    // Declare an instance of the class
    // ----------------------------------------------------------
    t3 = clock();
    FMRLasso l1(lasso_b1, con_crit, pterb, X.n_rows, X.n_cols, dists);
    // -------------------------------------------------------------
    // Golden search over the lasso paramter
    // -------------------------------------------------------------
    const double golden_ratio = 2 / (sqrt(5) + 1);
    //double lb = gs_str(1);
    //double ub = gs_str(0);
    double lb = lasso_b1;
    double ub = lasso_b2;
    double f1 = 0.0;
    double f2 = 0.0;
    double tolerance = 0.1;
    int iteration = 0;
    // Set the initial points with the golden ratio
    double x1 = ub - golden_ratio * (ub - lb);
    double x2 = lb + golden_ratio * (ub - lb);
    // Evaluate the function at these test points
    vec lambdas_str_1(dists);
    vec lambdas_str_2(dists);
    for (unsigned int i = 0; i < dists; i++)
    {
        lambdas_str_1(i) = x1;
        lambdas_str_2(i) = x2;
    }
    f1 = l1.BICRun(X, Y, beta_str, sig_srt, mu_srt, pi_srt, lambdas_str_1, dists, maxit, 1, outpath);
    f2 = l1.BICRun(X, Y, beta_str, sig_srt, mu_srt, pi_srt, lambdas_str_2, dists, maxit, 1, outpath);
    cout << "BIC at initial end points " << x1 << "," << f1 << " " << x2 << "," << f2 << endl;
    // Use golden search to find best lambda
    while (abs(ub - lb) > tolerance)
    {
        iteration++;
        if (f2 > f1)
        {
            // Minimum is to the left of x2
            // Set x2 to be the new upper bound
            ub = x2;
            // Set new upper test point
            x2 = x1;
            f2 = f1;
            x1 = ub - golden_ratio * (ub - lb);
            for (unsigned int i = 0; i < dists; i++)
            {
                lambdas_str_1(i) = x1;
            }
            f1 = l1.BICRun(X, Y, beta_str, sig_srt, mu_srt, pi_srt, lambdas_str_1, dists, maxit, 0, outpath);
        }
        else
        {
            // Minimum is to the right of x1
            lb = x1;
            x1 = x2;
            f1 = f2;
            x2 = lb + golden_ratio * (ub - lb);
            for (unsigned int i = 0; i < dists; i++)
            {
                lambdas_str_2(i) = x2;
            }
            f2 = l1.BICRun(X, Y, beta_str, sig_srt, mu_srt, pi_srt, lambdas_str_2, dists, maxit, 0, outpath);
        }
        cout << "Bounds at iteration "       << iteration << " are (" << lb << "," << ub << ")" << endl;
        cout << "BIC bounds at iteration "   << iteration << " are (" << f1 << "," << f2 << ")" << endl;
    }
    double gs_min;
    gs_min = (ub + lb) / 2;
    cout << "Estimate of minimum Golden search lambda 1 " << gs_min << endl;
    double min_bic1 = 0.0;
    vec lambdas_min(dists);
    for (unsigned int i = 0; i < dists; i++)
    {
        lambdas_min(i) = gs_min;
    }
    min_bic1 = l1.BICRun(X, Y, beta_str, sig_srt, mu_srt, pi_srt, lambdas_min, dists, maxit, 1, outpath);
    // ---------------------------------------------------------------------------
    // Update the betas so we are not starting from the a bad position for simplex
    // ---------------------------------------------------------------------------
    mat beta_update;
    int no_active = X.n_cols * dists;
    mat est_ind(X.n_cols, dists, fill::zeros);
    beta_update = l1.FMRLassoRun(X, Y, sig_srt, mu_srt, pi_srt,
                                 lambdas_min, maxit, 0, beta_str, est_ind,
                                 &min_bic1, no_active, 0, outpath);
    // cout << "What's at beta " << beta_update << endl;
    cout << "Minimum Golden search BIC " << min_bic1 << endl;
    cout << "Estimate of minimum Golden search lambda " << gs_min << endl;
    // ------------------------------------------------------------------------------
    // Simplex algoirthm search over two parameters starting at Golden search minimum
    // ------------------------------------------------------------------------------
    vector<double> simpl_res;
    vector<double> init;
    int a = 1;
    double arnd = 0.0;
    for (unsigned int i = 0; i < dists; i++)
    {
        a = a * -1 ;
        arnd = (gs_min + 0.5 * a);
        init.push_back(arnd);
    }
    simpl_res = l1.Simplex(init, 1E-3, 100, X, Y, beta_update, sig_srt, mu_srt, pi_srt, lambdas_min,
                           dists, maxit, 0, outpath);
    // ------------------------------------------------------------------------------
    // Run final set with this lambda and report
    // ------------------------------------------------------------------------------
    double min_bic2;
    vec simpl_res_arma(dists);
    for (unsigned int i = 0; i < dists; i++)
    {
        simpl_res_arma(i) = simpl_res.at(i);
    }
    min_bic2 = l1.BICRun(X, Y, beta_update, sig_srt, mu_srt, pi_srt, simpl_res_arma, dists, maxit, 1, outpath);
    cout << "Minimum simplex BIC " << min_bic2 << endl;
    cout << "Best simplex lambda vector is " << simpl_res_arma << endl;
    // ------------------------------------------------------------------------------
    // Print out the time to estimate parameters
    // ------------------------------------------------------------------------------
    t4 = clock();
    float diff2 ((float)t4 - (float)t3);
    float seconds2 = diff2 / CLOCKS_PER_SEC;
    cout << "Time to estimate lasso parameters was " << seconds2 << " seconds" << endl;
    return 0;
}

