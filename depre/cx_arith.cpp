/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

#include <o2scl/cx_arith.h>

using namespace std;
using namespace o2scl;

gsl_complex o2scl::complex_to_gsl(std::complex<double> &d) {
  gsl_complex g;
  g.dat[0]=d.real();
  g.dat[1]=d.imag();
  return g;
}

std::complex<double> o2scl::gsl_to_complex(gsl_complex &g) {
  std::complex<double> d(g.dat[0],g.dat[1]);
  return d;
}

gsl_complex o2scl::operator+(gsl_complex x, gsl_complex y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]+y.dat[0];
  r.dat[1]=x.dat[1]+y.dat[1];
  return r;
}

gsl_complex o2scl::operator-(gsl_complex x, gsl_complex y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]-y.dat[0];
  r.dat[1]=x.dat[1]-y.dat[1];
  return r;
}

gsl_complex o2scl::operator*(gsl_complex x, gsl_complex y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]*y.dat[0]-x.dat[1]*y.dat[1];
  r.dat[1]=x.dat[0]*y.dat[1]+x.dat[1]*y.dat[0];
  return r;
}

gsl_complex o2scl::operator/(gsl_complex x, gsl_complex y) {
  gsl_complex r;
  double s=1.0/hypot(y.dat[0],y.dat[1]);
  double sbr=s*y.dat[0];
  double sbi=s*y.dat[1];
  r.dat[0]=(x.dat[0]*sbr+x.dat[1]*sbi)*s;
  r.dat[1]=(x.dat[1]*sbr-x.dat[0]*sbi)*s;
  return r;
}

gsl_complex o2scl::operator+=(gsl_complex &x, gsl_complex y) {
  x.dat[0]=x.dat[0]+y.dat[0];
  x.dat[1]=x.dat[1]+y.dat[1];
  return x;
}

gsl_complex o2scl::operator-=(gsl_complex &x, gsl_complex y) {
  x.dat[0]=x.dat[0]-y.dat[0];
  x.dat[1]=x.dat[1]-y.dat[1];
  return x;
}

gsl_complex o2scl::operator*=(gsl_complex &x, gsl_complex y) {
  x.dat[0]=x.dat[0]*y.dat[0]-x.dat[1]*y.dat[1];
  x.dat[1]=x.dat[0]*y.dat[1]+x.dat[1]*y.dat[0];
  return x;
}

gsl_complex o2scl::operator/=(gsl_complex &x, gsl_complex y) {
  double s=1.0/hypot(y.dat[0],y.dat[1]);
  double sbr=s*y.dat[0];
  double sbi=s*y.dat[1];
  double tmp=(x.dat[0]*sbr+x.dat[1]*sbi)*s;
  x.dat[1]=(x.dat[1]*sbr-x.dat[0]*sbi)*s;
  x.dat[0]=tmp;
  return x;
}

gsl_complex o2scl::operator+(gsl_complex x, double y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]+y;
  r.dat[1]=x.dat[1];
  return r;
}

gsl_complex o2scl::operator+(double y, gsl_complex x) {
  gsl_complex r;
  r.dat[0]=y+x.dat[0];
  r.dat[1]=x.dat[1];
  return r;
}

gsl_complex o2scl::operator-(gsl_complex x, double y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]-y;
  r.dat[1]=x.dat[1];
  return r;
}

gsl_complex o2scl::operator-(double y, gsl_complex x) {
  gsl_complex r;
  r.dat[0]=y-x.dat[0];
  r.dat[1]=-x.dat[1];
  return r;
}

gsl_complex o2scl::operator*(gsl_complex x, double y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]*y;
  r.dat[1]=x.dat[1]*y;
  return r;
}

gsl_complex o2scl::operator*(double y, gsl_complex x) {
  gsl_complex r;
  r.dat[0]=y*x.dat[0];
  r.dat[1]=y*x.dat[1];
  return r;
}

gsl_complex o2scl::operator/(gsl_complex x, double y) {
  gsl_complex r;
  r.dat[0]=x.dat[0]/y;
  r.dat[1]=x.dat[1]/y;
  return r;
}

double o2scl::arg(gsl_complex x) {
  return atan2(x.dat[1],x.dat[0]);
}

double o2scl::abs(gsl_complex x) {
  return hypot(x.dat[0],x.dat[1]);
}

double o2scl::abs2(gsl_complex z) {
  double x=z.dat[0];
  double y=z.dat[1];
  return x*x+y*y;
}

double o2scl::logabs(gsl_complex z) {
  double xabs = fabs (z.dat[0]);
  double yabs = fabs (z.dat[1]);
  double max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }
  
  // Handle underflow when u is close to 0
  return std::log (max) + 0.5 * log1p (u * u);
}

gsl_complex o2scl::conjugate(gsl_complex a) {
  gsl_complex z={{a.dat[0],-a.dat[1]}};
  return z;
}

gsl_complex o2scl::sqrt(gsl_complex a) {
  
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double x = fabs (GSL_REAL (a));
      double y = fabs (GSL_IMAG (a));
      double w;

      if (x >= y)
	{
	  double t = y / x;
	  w = std::sqrt(x)*std::sqrt (0.5*(1.0+std::sqrt(1.0+t*t)));
	}
      else
	{
	  double t = x / y;
	  w = std::sqrt(y)*std::sqrt(0.5*(t+std::sqrt(1.0+t*t)));
	}

      if (GSL_REAL (a) >= 0.0)
	{
	  double ai = GSL_IMAG (a);
	  GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
	}
      else
	{
	  double ai = GSL_IMAG (a);
	  double vi = (ai >= 0) ? w : -w;
	  GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
	}
    }

  return z;
}

gsl_complex o2scl::sqrt_real(double x) {

  gsl_complex z;

  if (x >= 0)
    {
      GSL_SET_COMPLEX (&z, std::sqrt (x), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 0.0, std::sqrt (-x));
    }

  return z;
}

gsl_complex o2scl::exp(gsl_complex a) {           
  double rho = std::exp (GSL_REAL (a));
  double theta = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, rho * std::cos (theta), rho * std::sin (theta));
  return z;
}

gsl_complex o2scl::pow(gsl_complex a, gsl_complex b) {           
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0.0, 0.0);
    }
  else
    {
      double logr = logabs (a);
      double theta = arg (a);

      double br = GSL_REAL (b), bi = GSL_IMAG (b);

      double rho = std::exp (logr * br - bi * theta);
      double beta = theta * br + bi * logr;
      
      GSL_SET_COMPLEX (&z, rho * std::cos (beta), rho * std::sin (beta));
    }

  return z;
}

gsl_complex o2scl::pow_real(gsl_complex a, double b) {
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double logr = logabs (a);
      double theta = arg (a);
      double rho = std::exp (logr * b);
      double beta = theta * b;
      GSL_SET_COMPLEX (&z, rho * std::cos (beta), rho * std::sin (beta));
    }

  return z;
}

gsl_complex o2scl::log(gsl_complex a) {
  double logr = logabs (a);
  double theta = arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, logr, theta);
  return z;
}

gsl_complex o2scl::log10(gsl_complex a) {
  return gsl_complex_mul_real (log (a), 1 / std::log (10.));
}

gsl_complex o2scl::log_b(gsl_complex a, gsl_complex b) {
  return gsl_complex_div (log (a), log (b));
}

gsl_complex o2scl::sin(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0) 
    {
      // avoid returing negative zero (-0.0) for the imaginary part  

      GSL_SET_COMPLEX (&z, std::sin (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, std::sin (R) * std::cosh (I), 
		       std::cos (R) * std::sinh (I));
    }

  return z;
}

gsl_complex o2scl::cos(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0) 
    {
      // avoid returing negative zero (-0.0) for the imaginary part  

      GSL_SET_COMPLEX (&z, std::cos (R), 0.0);  
    } 
  else 
    {
      GSL_SET_COMPLEX (&z, std::cos (R) * std::cosh (I), 
		       std::sin (R) * std::sinh (-I));
    }

  return z;
}

gsl_complex o2scl::tan(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs (I) < 1)
    {
      double D = std::pow (std::cos (R), 2.0) + 
	std::pow (std::sinh (I), 2.0);

      GSL_SET_COMPLEX (&z, 0.5 * std::sin (2 * R) / D, 
		       0.5 * std::sinh (2 * I) / D);
    }
  else
    {
      double u = std::exp (-I);
      double C = 2 * u / (1 - std::pow (u, 2.0));
      double D = 1 + std::pow (std::cos (R), 2.0) * std::pow (C, 2.0);

      double S = std::pow (C, 2.0);
      double T = 1.0 / std::tanh (I);

      GSL_SET_COMPLEX (&z, 0.5 * std::sin (2 * R) * S / D, T / D);
    }

  return z;
}

gsl_complex o2scl::sec(gsl_complex a) {
  gsl_complex z = gsl_complex_cos (a);
  return gsl_complex_inverse (z);
}

gsl_complex o2scl::csc(gsl_complex a) {
  gsl_complex z = gsl_complex_sin (a);
  return gsl_complex_inverse(z);
}


gsl_complex o2scl::cot(gsl_complex a) {
  gsl_complex z = gsl_complex_tan (a);
  return gsl_complex_inverse (z);
}

gsl_complex o2scl::asin(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arcsin_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
	{
	  real = std::asin (B);
	}
      else
	{
	  if (x <= 1)
	    {
	      double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
	      real = std::atan (x / std::sqrt (D));
	    }
	  else
	    {
	      double Apx = A + x;
	      double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
	      real = std::atan (x / (y * std::sqrt (D)));
	    }
	}

      if (A <= A_crossover)
	{
	  double Am1;

	  if (x < 1)
	    {
	      Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
	    }
	  else
	    {
	      Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
	    }

	  imag = log1p (Am1 + std::sqrt (Am1 * (A + 1)));
	}
      else
	{
	  imag = std::log (A + std::sqrt (A * A - 1));
	}

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
    }

  return z;
}

gsl_complex o2scl::asin_real(double a) {
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, std::asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (&z, -M_PI_2, ::acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, M_PI_2, -::acosh (a));
	}
    }

  return z;
}

gsl_complex o2scl::acos(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arccos_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
	{
	  real = ::acos (B);
	}
      else
	{
	  if (x <= 1)
	    {
	      double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
	      real = ::atan (std::sqrt (D) / x);
	    }
	  else
	    {
	      double Apx = A + x;
	      double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
	      real = ::atan ((y * std::sqrt (D)) / x);
	    }
	}

      if (A <= A_crossover)
	{
	  double Am1;

	  if (x < 1)
	    {
	      Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
	    }
	  else
	    {
	      Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
	    }

	  imag = log1p (Am1 + std::sqrt (Am1 * (A + 1)));
	}
      else
	{
	  imag = std::log (A + std::sqrt (A * A - 1));
	}

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : M_PI - real, 
		       (I >= 0) ? -imag : imag);
    }

  return z;
}

gsl_complex o2scl::acos_real(double a) {
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, ::acos (a), 0);
    }
  else
    {
      if (a < 0.0)
	{
	  GSL_SET_COMPLEX (&z, M_PI, -::acosh (-a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, 0, ::acosh (a));
	}
    }

  return z;
}

gsl_complex o2scl::atan(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      GSL_SET_COMPLEX (&z, ::atan (R), 0);
    }
  else
    {
      // FIXME: This is a naive implementation which does not fully
      //take into account cancellation errors, overflow, underflow
      //etc.  It would benefit from the Hull et al treatment.

      double r = hypot (R, I);

      double imag;

      double u = 2 * I / (1 + r * r);

      // FIXME: the following cross-over should be optimized but 0.1
      //seems to work ok 

      if (fabs (u) < 0.1)
	{
	  imag = 0.25 * (log1p (u) - log1p (-u));
	}
      else
	{
	  double A = hypot (R, I + 1);
	  double B = hypot (R, I - 1);
	  imag = 0.5 * std::log (A / B);
	}

      if (R == 0)
	{
	  if (I > 1)
	    {
	      GSL_SET_COMPLEX (&z, M_PI_2, imag);
	    }
	  else if (I < -1)
	    {
	      GSL_SET_COMPLEX (&z, -M_PI_2, imag);
	    }
	  else
	    {
	      GSL_SET_COMPLEX (&z, 0, imag);
	    };
	}
      else
	{
	  GSL_SET_COMPLEX (&z, 0.5 * atan2 (2 * R, 
					    ((1 + r) * (1 - r))), imag);
	}
    }

  return z;
}

gsl_complex o2scl::asec(gsl_complex a) {
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arccos (z);
}

gsl_complex o2scl::asec_real(double a) {
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, ::acos (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (&z, 0, ::acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, M_PI, -::acosh (-1 / a));
	}
    }

  return z;
}

gsl_complex o2scl::acsc(gsl_complex a) {
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arcsin (z);
}

gsl_complex o2scl::acsc_real(double a) {
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, ::asin (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
	{
	  GSL_SET_COMPLEX (&z, M_PI_2, -::acosh (1 / a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, -M_PI_2, ::acosh (-1 / a));
	}
    }

  return z;
}

gsl_complex o2scl::acot(gsl_complex a) {
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, M_PI_2, 0);
    }
  else
    {
      z = gsl_complex_inverse (a);
      z = gsl_complex_arctan (z);
    }

  return z;
}

gsl_complex o2scl::sinh(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, std::sinh (R) * std::cos (I), 
		   std::cosh (R) * std::sin (I));
  return z;
}

gsl_complex o2scl::cosh(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, std::cosh (R) * std::cos (I), 
		   std::sinh (R) * std::sin (I));
  return z;
}

gsl_complex o2scl::tanh(gsl_complex a) {
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs(R) < 1.0) 
    {
      double D = std::pow (std::cos (I), 2.0) + std::pow (std::sinh (R), 2.0);
      
      GSL_SET_COMPLEX (&z, std::sinh (R) * std::cosh (R) / D, 
		       0.5 * std::sin (2 * I) / D);
    }
  else
    {
      double D = std::pow (std::cos (I), 2.0) + std::pow (std::sinh (R), 2.0);
      double F = 1 + std::pow (std::cos (I) / std::sinh (R), 2.0);

      GSL_SET_COMPLEX (&z, 1.0 / (std::tanh (R) * F), 
		       0.5 * std::sin (2 * I) / D);
    }

  return z;
}

gsl_complex o2scl::sech(gsl_complex a) {
  gsl_complex z = gsl_complex_cosh (a);
  return gsl_complex_inverse (z);
}

gsl_complex o2scl::csch(gsl_complex a) {
  gsl_complex z = gsl_complex_sinh (a);
  return gsl_complex_inverse (z);
}

gsl_complex o2scl::coth(gsl_complex a) {
  gsl_complex z = gsl_complex_tanh (a);
  return gsl_complex_inverse (z);
}

gsl_complex o2scl::asinh(gsl_complex a) {
  gsl_complex z = gsl_complex_mul_imag(a, 1.0);
  z = gsl_complex_arcsin (z);
  z = gsl_complex_mul_imag (z, -1.0);
  return z;
}

gsl_complex o2scl::acosh(gsl_complex a) {
  gsl_complex z = gsl_complex_arccos (a);
  z = gsl_complex_mul_imag (z, GSL_IMAG(z) > 0 ? -1.0 : 1.0);
  return z;
}

gsl_complex o2scl::acosh_real(double a) {
  gsl_complex z;

  if (a >= 1)
    {
      GSL_SET_COMPLEX (&z, ::acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
	{
	  GSL_SET_COMPLEX (&z, 0, ::acos (a));
	}
      else
	{
	  GSL_SET_COMPLEX (&z, ::acosh (-a), M_PI);
	}
    }

  return z;
}

gsl_complex o2scl::atanh(gsl_complex a) {
  if (GSL_IMAG (a) == 0.0)
    {
      return gsl_complex_arctanh_real (GSL_REAL (a));
    }
  else
    {
      gsl_complex z = gsl_complex_mul_imag(a, 1.0);
      z = gsl_complex_arctan (z);
      z = gsl_complex_mul_imag (z, -1.0);
      return z;
    }
}

gsl_complex o2scl::atanh_real(double a) {
  gsl_complex z;

  if (a > -1.0 && a < 1.0)
    {
      GSL_SET_COMPLEX (&z, ::atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, ::atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

gsl_complex o2scl::asech(gsl_complex a) {
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arccosh (t);
}

gsl_complex o2scl::acsch(gsl_complex a) {
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arcsinh (t);
}

gsl_complex o2scl::acoth(gsl_complex a) {
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arctanh (t);
}
