/************************************************************************
 ************************************************************************
 **                                                                    **
 **           Finite Mixture of Regressions Lasso HPP                  **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

// ---------------------------------------------------------------------
// Author: Luke Lloyd-Jones
// Date started:      04/09/2015
// Date last updated: 05/08/2015
// ---------------------------------------------------------------------

#ifndef FMRLASSOHEADERDEF
#define FMRLASSOHEADERDEF

#include <cmath>
#include <algorithm>
#include <armadillo>
#include <time.h>
#include <boost/math/distributions/normal.hpp>

using namespace arma;
using namespace std;
using namespace boost::math;

class FMRLasso
{
private:
    double mlambda; // lasso parameter
    double mcriterion;
    double mepsilon;
    unsigned int mN;
    unsigned int mP;
    unsigned int mG;
public:
    FMRLasso(double lambda,
             double criterion,
             double epsilon,
             unsigned int N,
             unsigned int P,
             unsigned int G);
    ~FMRLasso();
    // Declare the variables to be used in algorithm
    unsigned int maxit;
    unsigned int nt_zero;
    double obj;
    mat beta_old;
    mat beta;
    const mat beta_str;
    mat resid;
    vec lambda;
    vec alpha;
    vec sig_2;
    vec pi;
    vec pi_brent;
    double mean;
    double stdev;
    double inner;
    double old_obj;
    mat tau;
    mat weight;
    mat geno_sqr;
    vec rho;
    vec psi;
    double gamma_c1;
    double gamma_c2;
    double gamma_1;
    double gamma_2;
    vec pi_1;
    vec pi_2;
    double c1 = 0.0;
    double c2 = 0.0;
    double c3 = 0.0;
    // Main coordinate descent method
    mat FMRLassoRun(const mat& geno,
                    const mat& pheno,
                    vec sigs,
                    vec mus,
                    vec pis,
                    vec lambdas,
                    unsigned int maxit,
                    unsigned int mode,
                    const mat& beta_str,
                    const mat& beta_ind,
                    double* bic,
                    int no_active,
                    int verbose,
                    string outpath);
    // Method to standardise the genotype matrix
    mat GenoStd(mat& X);
    void ActiveSet(mat& est, mat& est_ind, int* no_not_active);
    double BICRun(const mat& geno,
                  const mat& pheno,
                  const mat& beta_str,
                  vec sigs,
                  vec mus,
                  vec pis,
                  vec lambdas,
                  double dists,
                  unsigned int maxit,
                  int verbose,
                  string outpath);
    // Simplex method
    vector<double> Simplex(vector<double> init,
                           double tol,
                           int iterations,
                           const mat& geno,
                           const mat& pheno,
                           const mat& beta_update,
                           vec sigs,
                           vec mus,
                           vec pis,
                           vec lambdas,
                           double dists,
                           unsigned int maxit,
                           int verbose,
                           string outpath);
};



#endif
