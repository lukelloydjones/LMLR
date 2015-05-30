// ---------------------------------------------------------------------
// Author: Luke Lloyd-Jones
// Date: 30/05/2015
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Includes and headers
// ---------------------------------------------------------------------


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
#include <algorithm>
#include <iomanip>  // setprecision()
#include <string>
#include <vector>
#include <cerrno> // something to help with reading in buffers
#include <armadillo>
#include <time.h>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>

// Set the namespace std

using namespace std;
using namespace arma;

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                    Finite Mixture of Regressions                   **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main(int argc, char* argv[])
{
    // Preparing some arguments
    
    std::cout << "Number of command line arguments " << argc << "\n";
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Argument " << i << " = " << argv[i];
        std::cout << "\n";
    }
    //clock_t t1, t2, t3, t4;
    //t1 = clock();
    int ncol = 15;
    int nrow = 20;
    int g = 2; // Number of groups in 1, ..., k
    int p = ncol; // Number of covariates in 1, ..., p
    int n = nrow; // Number of individuals in 1, ..., n
    double eff_size = 5;
    int nt_zero = 5;
    double lambda_str = 1 / 25.0;
    mat X(nrow, ncol);
    X.load(argv[1], csv_ascii);
    Y.load(argv[2], csv_ascii);
    // Assume we have read in the genotype and phenotype matrices
    // Call them X and Y in the code
    // ---------------------------------------------------
    // Initialise the beta vector, lambda (lasso penalty),
    // alpha (intercept), variances, and pis
    // ---------------------------------------------------
    mat beta(g, p, fill::zeros);
    for (int j = 0; j < nt_zero; j++)
    {
        beta(0, j) = eff_size;
        beta(1, (j + nt_zero)) = eff_size;
        
    }
    vec lambda(g);
    vec alpha = {20, -20};
    vec sig_2(g);
    vec pi(g);
    for (int k = 0; k < g; k++)
    {
        lambda(k) = lambda_str;
        sig_2(k)  = 1;
        pi(k)     = 1 / (double)(g);
    }
    // ---------------------------------------------------
    // Printing some elements to see how we are going
    // ---------------------------------------------------
    cout << lambda << endl;
    cout << alpha << endl;
    cout << sig_2 << endl;
    cout << pi << endl;
    cout << beta << endl;
    // ---------------------------------------------------
    // Initialise the objective function
    // ---------------------------------------------------
    double obj = 0.0;
    double mean;
    double stdev;
    double Yi = 0.25
    for (int i = 0; i < n; i++)
    {
        double inner = 0.0;
        for (int k = 0; k < g; k++)
        {
            cout << beta.row(k)(i) << endl;
            mean = alpha(k) + sum(beta.row(k) % X.row(i));
            stdev = sqrt(sig_2(k));
            Yi = 0.25;
            normal_distribution<> myNormal(mean, stdev);
            cout << "pdf: " << pdf(myNormal, Y(i)) << endl;
            // inner += Y(i) alpha(k) beta.row(k)(i) X.row(i) sig_2(k)
            
            
        }
        obj = obj + (1 / n) * log(inner);
    }
    
    
    
    
    
    
    
    
    return 0;
}

//normal_distribution<> myNormal(1.0, 10.0);
//cout << "Mean: " << myNormal.mean() << ", standard deviation: " << myNormal.standard_deviation() << endl;
//
//// Distributional properties
//double x = 10.25;
//
//cout << "pdf: " << pdf(myNormal, x) << endl;
//cout << "cdf: " << cdf(myNormal, x) << endl;
//
//// Choose another data type and now a N(0,1) variate
//normal_distribution<float> myNormal2;
//cout << "Mean: " << myNormal2.mean() << ", standard deviation: " << myNormal2.standard_deviation() << endl;
//
//cout << "pdf: " << pdf(myNormal2, x) << endl;
//cout << "cdf: " << cdf(myNormal2, x) << endl;
//cout << "cdf complement:" << cdf(complement(myNormal2, x));
//
//// Choose precision
//cout.precision(10); // Number of values behind the comma
//
//// Other properties
//cout << "n***normal distribution: n";
//cout << "mean: " << mean(myNormal) << endl;
//cout << "variance: " << variance(myNormal) << endl;
//cout << "median: " << median(myNormal) << endl;
//cout << "mode: " << mode(myNormal) << endl;
//cout << "kurtosis excess: " << kurtosis_excess(myNormal; cout << "kurtosis: " << kurtosis(myNormal) << endl;
//                                               cout << "skewness: " << skewness(myNormal) << endl;
//                                               cout << "characteristic function: " << chf(myNormal, x);
//                                               cout << "hazard: " << hazard(myNormal, x) << endl;
//                                               cout << "cumulative hazard: " << chf(myNormal, x) << endl;
//                                               
//                                               // Return the value of random variable
//                                               // s.t. dist(myNormal, x) = p
//                                               double p = 0.3;
//                                               cout << "quantile: " << quantile(myNormal, p) << endl;
//                                               
//                                               // Other properties; these functions return a pair
//                                               cout << "range: (" << range(myNormal).first << ","
//                                               << range(myNormal).second << ")" << endl;
//                                               cout << "support: (" << support(myNormal).first << ","
//                                               << support(myNormal).second << ")" << endl;
//








