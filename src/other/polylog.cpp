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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/polylog.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

polylog::polylog() {
  li2neg1=-1.0/12.0*pi2;
  li4neg1=-7.0/720.0*pi2*pi2;
  li6neg1=-31.0/30240.0*pi2*pi2*pi2;
  ifstream fp;
  int i;
  string fn=o2scl_settings.get_data_dir();
  fn+="litable";
  arg=new double[102];
  two=new double[102];
  three=new double[102];
  four=new double[102];
  five=new double[102];
  six=new double[102];
  fp.open(fn.c_str());
  for(i=1;i<=101;i++) {
    fp >> arg[i] >> two[i] >> three[i] >> four[i] >> five[i] >> six[i];
  }
  fp.close();
  return;
}

polylog::~polylog() {
  delete[] arg;
  delete[] two;
  delete[] three;
  delete[] four;
  delete[] five;
  delete[] six;
}

double polylog::li0(double x) {
  if (x==1.0) {
    O2SCL_ERR("Infinite result (argument==1) for Li0.",exc_efailed);
    return 0.0;
  }
  return x/(1.0-x);
}

double polylog::li1(double x) {
  if (x>1.0) {
    O2SCL_ERR("Infinite result (argument>1) for Li1.",exc_efailed);
    return 0.0;
  }
  return -1.0*log(1.0-x);
}

double polylog::li2(double x) {
  return gsl_sf_dilog(x);
}

double polylog::li3(double x) {
  int i;
  if (x>0.0) {
    O2SCL_ERR("Infinite result for Li3.",exc_efailed);
  } else if (x>=-1.0) {
    i=((int)(x*(-100.0)))+1;
    if (i>100) i=100;
    return(three[i]+(three[i+1]-three[i])*(x-arg[i])/(arg[i+1]-arg[i]));
  } else {
    return(-1.0/6.0*pow(log(-1.0*x),3.0)+2.0*log(-1.0*x)*li2neg1+li3(1.0/x));
  }
  return 0.0;
}

double polylog::li4(double x) {
  int i;
  if (x>0.0) {
    O2SCL_ERR("Infinite result for Li4.",exc_efailed);
  } else if (x>=-1.0) {
    i=((int)(x*(-100.0)))+1;
    if (i>100) i=100;
    return(four[i]+(four[i+1]-four[i])*(x-arg[i])/(arg[i+1]-arg[i]));
  } else {
    return(-1.0/24.0*pow(log(-1.0*x),4.0)+pow(log(-1.0*x),2.0)*li2neg1+
	   2.0*li4neg1-li4(1.0/x));
  }
  return 0.0;
}

double polylog::li5(double x) {
  int i;
  if (x>0.0) {
    O2SCL_ERR("Infinite result for Li5.",exc_efailed);
  } else if (x>=-1.0) {
    i=((int)(x*(-100.0)))+1;
    if (i>100) i=100;
    return(five[i]+(five[i+1]-five[i])*(x-arg[i])/(arg[i+1]-arg[i]));
  } else {
    return(-1.0/120.0*pow(log(-1.0*x),5.0)+pow(log(-1.0*x),3.0)*li2neg1/3.0+
	   2.0*log(-1.0*x)*li4neg1+li5(1.0/x));
  }
  return 0.0;
}

double polylog::li6(double x) {
  int i;
  if (x>0.0) {
    O2SCL_ERR("Infinite result for Li6.",exc_efailed);
  } else if (x>=-1.0) {
    i=((int)(x*(-100.0)))+1;
    if (i>100) i=100;
    return(six[i]+(six[i+1]-six[i])*(x-arg[i])/(arg[i+1]-arg[i]));
  } else {
    return(-1.0/720.0*pow(log(-1.0*x),6.0)+pow(log(-1.0*x),4.0)*li2neg1/12.0+
	   pow(log(-1.0*x),2.0)*li4neg1+2.0*li6neg1-li6(1.0/x));
  }
  return 0.0;
}

