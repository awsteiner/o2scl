/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#include <o2scl/mmin_bfgs2.h>

using namespace std;
using namespace o2scl;

double mmin_linmin_gsl::interp_quad(double f0, double fp0, double f1, 
				     double zl, double zh) {
      
  double fl = f0 + zl*(fp0 + zl*(f1 - f0 -fp0));
  double fh = f0 + zh*(fp0 + zh*(f1 - f0 -fp0));
  /* curvature */
  double c = 2 * (f1 - f0 - fp0);       

  double zmin = zl, fmin = fl;

  if (fh < fmin) { zmin = zh; fmin = fh; }

  /* positive curvature required for a minimum */
  if (c > 0)  {

    /* location of minimum */
    double z = -fp0 / c;      
    if (z > zl && z < zh) {
      double f = f0 + z*(fp0 + z*(f1 - f0 -fp0));
      if (f < fmin) { zmin = z; fmin = f; };
    }
  }
  
  return zmin;
}

double mmin_linmin_gsl::cubic(double c0, double c1, double c2, 
			      double c3, double z) {
  return c0 + z * (c1 + z * (c2 + z * c3));
}
    
void mmin_linmin_gsl::check_extremum(double c0, double c1, double c2, 
				     double c3, double z,
				     double *zmin, double *fmin) {
  /* could make an early return by testing curvature >0 for minimum */
  
  double y = cubic (c0, c1, c2, c3, z);
  
  if (y < *fmin) {

    /* accepted new point*/
    *zmin = z;  
    *fmin = y;
  }

}

double mmin_linmin_gsl::interp_cubic(double f0, double fp0, double f1, 
				     double fp1, double zl, double zh) {
  double eta = 3 * (f1 - f0) - 2 * fp0 - fp1;
  double xi = fp0 + fp1 - 2 * (f1 - f0);
  double c0 = f0, c1 = fp0, c2 = eta, c3 = xi;
  double zmin, fmin;
  double z0, z1;

  zmin = zl; 
  fmin = cubic(c0, c1, c2, c3, zl);
  check_extremum (c0, c1, c2, c3, zh, &zmin, &fmin);

  {
    int n = gsl_poly_solve_quadratic (3 * c3, 2 * c2, c1, &z0, &z1);
    
    if (n == 2) { 
      /* found 2 roots */
      if (z0 > zl && z0 < zh) {
	check_extremum (c0, c1, c2, c3, z0, &zmin, &fmin);
      }
      if (z1 > zl && z1 < zh) {
	check_extremum (c0, c1, c2, c3, z1, &zmin, &fmin);
      }
    } else if (n == 1)  {
      /* found 1 root */
      if (z0 > zl && z0 < zh) {
	check_extremum (c0, c1, c2, c3, z0, &zmin, &fmin);
      }
    }
  }

  return zmin;
}

double mmin_linmin_gsl::interpolate(double a, double fa, double fpa, 
				     double b, double fb, double fpb, 
				     double xmin, double xmax, int order) {
      
  /* Map [a,b] to [0,1] */
  double z, alpha, zmin, zmax;

  zmin = (xmin - a) / (b - a);
  zmax = (xmax - a) / (b - a);

  if (zmin > zmax) {
    double tmp = zmin;
    zmin = zmax;
    zmax = tmp;
  }

  if (order > 2 && GSL_IS_REAL(fpb)) {
    z = interp_cubic (fa, fpa * (b - a), fb, fpb * (b - a), zmin, zmax);
  } else {
    z = interp_quad (fa, fpa * (b - a), fb, zmin, zmax);
  }

  alpha = a + z * (b - a);

  return alpha;
}

int mmin_linmin_gsl::minimize(mmin_wrap_gsl &wrap, double rho, 
			      double sigma, double tau1, double tau2, 
			      double tau3, int order, double alpha1, 
			      double *alpha_new) {
  
  double f0, fp0, falpha, falpha_prev, fpalpha=0.0;
  double fpalpha_prev, delta, alpha_next;
  double alpha = alpha1, alpha_prev = 0.0;
  double a, b, fa, fb, fpa, fpb;
  const size_t bracket_iters = 100, section_iters = 100;
  size_t i = 0;
  
  wrap.wrap_fdf(0.0,&f0,&fp0);
  falpha_prev = f0;
  fpalpha_prev = fp0;

  /* Avoid uninitialized variables morning */
  a = 0.0; b = alpha;
  fa = f0; fb = 0.0;
  fpa = fp0; fpb = 0.0;

  /* Begin bracketing */

  while (i++ < bracket_iters) {

    falpha=wrap.wrap_f(alpha);
	
    /* Fletcher's rho test */

    if (falpha > f0 + alpha * rho * fp0 || falpha >= falpha_prev) {
      a = alpha_prev; fa = falpha_prev; fpa = fpalpha_prev;
      b = alpha; fb = falpha; fpb = GSL_NAN;
      break;                /* goto sectioning */
    }

    fpalpha=wrap.wrap_df(alpha);
	
    /* Fletcher's sigma test */

    if (fabs (fpalpha) <= -sigma * fp0) {
      *alpha_new = alpha;
      return success;
    }
	  
    if (fpalpha >= 0) {
      a = alpha; fa = falpha; fpa = fpalpha;
      b = alpha_prev; fb = falpha_prev; fpb = fpalpha_prev;
      break;                /* goto sectioning */
    }

    delta = alpha - alpha_prev;

    {
      double lower = alpha + delta;
      double upper = alpha + tau1 * delta;
	    
      alpha_next = interpolate (alpha_prev, falpha_prev, fpalpha_prev,
				alpha, falpha, fpalpha, lower, 
				upper, order);

    }

    alpha_prev = alpha;
    falpha_prev = falpha;
    fpalpha_prev = fpalpha;
    alpha = alpha_next;
	
    double temp = wrap.wrap_f(alpha);

  }

  /*  Sectioning of bracket [a,b] */

  while (i++ < section_iters) {
    delta = b - a;
	
    {
      double lower = a + tau2 * delta;
      double upper = b - tau3 * delta;
	  
      alpha = interpolate (a, fa, fpa, b, fb, fpb, lower, upper, order);
    }
	
    falpha=wrap.wrap_f(alpha);
	
#ifdef O2SCL_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif

    if ((a-alpha)*fpa <= dbl_eps) {
      /* roundoff prevents progress */
      return exc_enoprog;
    };
	
    if (falpha > f0 + rho * alpha * fp0 || falpha >= fa) {

      /*  a_next = a; */
      b = alpha; fb = falpha; fpb = GSL_NAN;

    } else {

      fpalpha=wrap.wrap_df(alpha);
	  
      if (fabs(fpalpha) <= -sigma * fp0) {
	*alpha_new = alpha;
	return success;  /* terminate */
      }
	  
      if ( ((b-a) >= 0 && fpalpha >= 0) || ((b-a) <=0 && fpalpha <= 0))  {
	b = a; fb = fa; fpb = fpa;
	a = alpha; fa = falpha; fpa = fpalpha;
      } else {
	a = alpha; fa = falpha; fpa = fpalpha;
      }
    }
  }
      
  return success;
}
