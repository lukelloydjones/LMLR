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
// Date last updated: 22/02/2016
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
    // Print the header
    // ----------------------------------------------------------
    cout << "          +--------------------------------------------------+" << endl;
    cout << "          |                                                  |" << endl;
    cout << "          | LMLR, v1.0                                       |" << endl;
    cout << "          | (C) 2016 Luke Lloyd-Jones and Hien Nguyen        |" << endl;
    cout << "          | GNU General Public License v3                    |" << endl;
    cout << "          |                                                  |" << endl;
    cout << "          +--------------------------------------------------+" << endl;
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
    double sig_1       = atof(argv[8]);
    double sig_2       = atof(argv[9]);
    double mu_1        = atof(argv[10]);
    double mu_2        = atof(argv[11]);
    double pi_1        = atof(argv[12]);
    double pi_2        = atof(argv[13]);
    static const std::string outpath  = argv[14];
    unsigned int dists = 2;
    // cout << "lasso value is " << lasso << endl;
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
    // Declare an instance of the class
    // ----------------------------------------------------------
    t3 = clock();
    FMRLasso l1(lasso_b1, con_crit, pterb, X.n_rows, X.n_cols, dists);
    // -------------------------------------------------------------
    // Golden search over the lasso paramter
    // -------------------------------------------------------------
    cout << "INITIALISING GOLDEN SECTION SEARCH ALGORITHM" << endl;
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
    f1 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, x1, x1, dists, maxit, 0, outpath);
    f2 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, x2, x2, dists, maxit, 0, outpath);
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
            f1 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, x1, x1, dists, maxit, 0, outpath);
        }
        else
        {
            // Minimum is to the right of x1
            lb = x1;
            x1 = x2;
            f1 = f2;
            x2 = lb + golden_ratio * (ub - lb);
            f2 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, x2, x2, dists, maxit, 0, outpath);
        }
        cout << "Bounds at iteration "     << iteration << " are (" << lb << "," << ub << ")" << endl;
        cout << "BIC bounds at iteration "   << iteration << " are (" << f1 << "," << f2 << ")" << endl;
    }
    double gs_min;
    gs_min = (ub + lb) / 2;
    cout << "Estimate of minimum Golden search lambda 1 " << gs_min << endl;
    double min_bic1 = 0.0;
    min_bic1 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, gs_min, gs_min, dists, maxit, 0, outpath);
    cout << "Minimum Golden search BIC " << min_bic1 << endl;
    cout << "Estimate of minimum Golden search lambda " << gs_min << endl;
    // ------------------------------------------------------------------------------
    // Simplex algoirthm search over two parameters starting at Golden search minimum
    // ------------------------------------------------------------------------------
    vector<double> simpl_res;
    vector<double> init;
    init.push_back(gs_min - 0.5);
    init.push_back(gs_min + 0.5);
    simpl_res = l1.Simplex(init, 1E-3, 50, X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, lasso_b1, lasso_b2,
                           dists, maxit, 0, outpath);
    // ------------------------------------------------------------------------------
    // Run final set with this lambda and report
    // ------------------------------------------------------------------------------
    double min_bic2;
    min_bic2 = l1.BICRun(X, Y, sig_1, sig_2, mu_1, mu_2, pi_1, pi_2, simpl_res.at(0), simpl_res.at(1), dists, maxit, 1, outpath);
    cout << "Minimum simplex BIC " << min_bic2 << endl;
    cout << "Best simplex lambda vector is (" << simpl_res.at(0) << ", " << simpl_res.at(1) << ")" << endl;
    // ------------------------------------------------------------------------------
    // Print out the time to estimate parameters
    // ------------------------------------------------------------------------------
    t4 = clock();
    float diff2 ((float)t4 - (float)t3);
    float seconds2 = diff2 / CLOCKS_PER_SEC;
    cout << "Time to estimate lasso parameters was " << seconds2 << " seconds" << endl;
    return 0;
}

