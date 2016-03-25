// ---------------------------------------------------------------------
// Author: Luke Lloyd-Jones
// Date: 22/05/2015
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// Includes and headers
// ---------------------------------------------------------------------


#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <functional>
#include <stdio.h>
#include <gsl/gsl_poly.h>

// ---------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------

int main(int argc, char* argv[])
{
    int i;
    /* coefficients of P(x) = -1 + x^5 */
    double a[3] = {4, 4, 1};
    double z[10];
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(3);
    gsl_poly_complex_solve (a, 3, w, z);
    gsl_poly_complex_workspace_free (w);
    for (i = 0; i < 5; i++)
    {
        printf ("z%d = %+.18f %+.18f\n", i, z[2 * i], z[2 * i + 1]);
    }

    return 0;

}











