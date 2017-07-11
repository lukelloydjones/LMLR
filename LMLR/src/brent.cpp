# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iostream>
#include <armadillo>
using namespace arma;

using namespace std;
const double check_tolerance(1e-12);
# include "brent.hpp"

namespace brent{

//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
  const int TIME_SIZE(40);

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
}

//****************************************************************************80

void zero_rc ( double a, double b, double t, double &arg, int &status,
  double value )

//****************************************************************************80
//
//  Purpose:
//
//    ZERO_RC seeks the root of a function F(X) using reverse communication.
//
//  Discussion:
//
//    The interval [A,B] must be a change of sign interval for F.
//    That is, F(A) and F(B) must be of opposite signs.  Then
//    assuming that F is continuous implies the existence of at least
//    one value C between A and B for which F(C) = 0.
//
//    The location of the zero is determined to within an accuracy
//    of 6 * MACHEPS * abs ( C ) + 2 * T.
//
//    The routine is a revised version of the Brent zero finder
//    algorithm, using reverse communication.
//
//    Thanks to Thomas Secretin for pointing out a transcription error in the
//    setting of the value of P, 11 February 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the change of sign interval.
//
//    Input, double T, a positive error tolerance.
//
//    Output, double &ARG, the currently considered point.  The user
//    does not need to initialize this value.  On return with STATUS positive,
//    the user is requested to evaluate the function at ARG, and return
//    the value in VALUE.  On return with STATUS zero, ARG is the routine's
//    estimate for the function's zero.
//
//    Input/output, int &STATUS, used to communicate between
//    the user and the routine.  The user only sets STATUS to zero on the first
//    call, to indicate that this is a startup call.  The routine returns STATUS
//    positive to request that the function be evaluated at ARG, or returns
//    STATUS as 0, to indicate that the iteration is complete and that
//    ARG is the estimated zero
//
//    Input, double VALUE, the function value at ARG, as requested
//    by the routine on the previous call.
//
{
  static double c;
  static double d;
  static double e;
  static double fa;
  static double fb;
  static double fc;
  double m;
  static double macheps;
  double p;
  double q;
  double r;
  double s;
  static double sa;
  static double sb;
  double tol;
//
//  Input STATUS = 0.
//  Initialize, request F(A).
//
  if ( status == 0 )
  {
    macheps = r8_epsilon ( );

    sa = a;
    sb = b;
    e = sb - sa;
    d = e;

    status = 1;
    arg = a;
    return;
  }
//
//  Input STATUS = 1.
//  Receive F(A), request F(B).
//
  else if ( status == 1 )
  {
    fa = value;
    status = 2;
    arg = sb;
    return;
  }
//
//  Input STATUS = 2
//  Receive F(B).
//
  else if ( status == 2 )
  {
    fb = value;

    if ( 0.0 < fa * fb )
    {
      status = -1;
      return;
    }
    c = sa;
    fc = fa;
  }
  else
  {
    fb = value;

    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
//
//  Compute the next point at which a function value is requested.
//
  if ( fabs ( fc ) < fabs ( fb ) )
  {
    sa = sb;
    sb = c;
    c = sa;
    fa = fb;
    fb = fc;
    fc = fa;
  }

  tol = 2.0 * macheps * fabs ( sb ) + t;
  m = 0.5 * ( c - sb );

  if ( fabs ( m ) <= tol || fb == 0.0 )
  {
    status = 0;
    arg = sb;
    return;
  }

  if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
  {
    e = m;
    d = e;
  }
  else
  {
    s = fb / fa;

    if ( sa == c )
    {
      p = 2.0 * m * s;
      q = 1.0 - s;
    }
    else
    {
      q = fa / fc;
      r = fb / fc;
      p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
      q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
    }

    if ( 0.0 < p )
    {
      q = - q;
    }
    else
    {
      p = - p;
    }
    s = e;
    e = d;

    if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
         p < fabs ( 0.5 * s * q ) )
    {
      d = p / q;
    }
    else
    {
      e = m;
      d = e;
    }
  }

  sa = sb;
  fa = fb;

  if ( tol < fabs ( d ) )
  {
    sb = sb + d;
  }
  else if ( 0.0 < m )
  {
    sb = sb + tol;
  }
  else
  {
    sb = sb - tol;
  }

  arg = sb;
  status = status + 1;

  return;
}

// ======================================================================
// === Simple wrapper functions
// === for convenience and/or compatibility.
//
// === The three functions are the same as above,
// === except that they take a plain function F
// === instead of a c++ functor.  In all cases, the
// === input and output of F() are of type double.

typedef double DoubleOfDouble (double);

class func_wrapper : public func_base {
  DoubleOfDouble* func;
public:
  func_wrapper(DoubleOfDouble* f) {
    func = f;
  }
  virtual double operator() (double x){
    return func(x);
  }
};

//****************************************************************************80
    
double zero_rc_pi_root ( double a, double b, double t, double f (double x, unsigned int g,
                         const vec& rho, const vec& psi),  unsigned int g,
                         const vec& rho, const vec& psi, string title )
    
//****************************************************************************80
    //
    //  Purpose:
    //
    //    TEST_ZERO_RC_ONE tests ZERO_RC on one test function.
    //
    //  Licensing:
    //s
    //    This code is distributed under the GNU LGPL license.
    //
    //  Modified:
    //
    //    10 June 2017
    //
    //  Author:
    //
    //    Luke Lloyd-Jones adapted from John Burkardt's code
    //
    //  Parameters:
    //
    //    Input, double A, B, the two endpoints of the change of sign
    //    interval.
    //
    //    Input, double MACHEP, an estimate for the relative machine
    //    precision.
    //
    //    Input, double T, a positive error tolerance.
    //
    //    Input, double F ( double x ), the name of a user-supplied
    //    function which evaluates the function whose zero is being sought.
    //
    //    Input, string TITLE, a title for the problem.
    //
    {
        double arg;
        int status;
        double value = 0.0;
        // cout << "\n";
        //cout << "  " << title << "\n";
        //cout << "\n";
        //cout << "    STATUS      X               F(X)\n";
        //cout << "\n";
        
        status = 0;
        
        for ( ; ; )
        {
            zero_rc ( a, b, t, arg, status, value );
            
            if ( status < 0 )
            {
                cout << "\n";
                cout << "  ZERO_RC returned an error flag!\n";
                break;
            }
            
            value = f ( arg, g, rho, psi );
            
            //cout << "  " << setw(8) << status
            //     << "  " << setw(14) << arg
            //     << "  " << setw(14) << value << "\n";
            
            if ( status == 0 )
            {
                break;
            }
        }
        if (abs(value) > check_tolerance) {
            cerr << "*** error ***" << endl;
            cerr << "final value " << value
            << " exceeds check_tolerance " << check_tolerance << endl;
            exit(1);
        }
        
        return arg;
}

//****************************************************************************
double pi_root ( double x, unsigned int g, const vec& rho, const vec& psi )
//****************************************************************************8
    //
    //  Purpose:
    //
    //    Evaluate the .
    //
    //  Licensing:
    //
    //    This code is distributed under the GNU LGPL license.
    //
    //  Author:
    //
    //    Luke Lloyd-Jones
    //
    //  Parameters:
    //
    //    Input, double X, the point at which pi_root is to be evaluated.
    //           g, the number of mixtures
    //           rho, a row vector of rho elements
    //           psi, a row vector of psi elements
    //
    //    Output, double value which corresponds to the solution of the Lagrangian
    //
{
    double value = 0.0;
    double eval = 0.0;
    for (unsigned int i = 0; i < g; ++i){
         eval = rho(i) / (x + psi(i));
         value = value + eval;
    }
    value = value - 1;
    return value;
}
    
// ======================================================================
// Generally useful functor to evaluate a monic polynomial.
// For details, see class definition in brent.hpp

double monicPoly::operator()(double x){
  double rslt(1);
  for (int ii = coeff.size()-1; ii >= 0; ii--){
    rslt *= x;
    rslt += coeff[ii];
  }
  return rslt;
}

// Similarly, evaluate a general polynomial (not necessarily monic):
double Poly::operator()(double x){
  double rslt(0);
  for (int ii = coeff.size()-1; ii >= 0; ii--){
    rslt *= x;
    rslt += coeff[ii];
  }
  return rslt;
}

} // end namespace brent
