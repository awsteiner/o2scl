/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2021, Andrew W. Steiner
  
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
#include <o2scl/prob_dens_func.h>

using namespace std;
using namespace o2scl;

prob_dens_hist::prob_dens_hist() {
  n=0;
}

prob_dens_hist::~prob_dens_hist() {
  if (n>0) {
    sum.clear();
    range.clear();
    n=0;
  }
}

double prob_dens_hist::lower_limit() const {
  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::lower_limit().",
	      exc_efailed);
  }
  return range[0];
}

double prob_dens_hist::upper_limit() const {
  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::lower_limit().",
	      exc_efailed);
  }
  return range[n];
}

void prob_dens_hist::init(hist &h) {
  
  n=h.size();
  sum.resize(n+1);
  range.resize(n+1);

  for(size_t i=0;i<n;i++) {
    if (h[i]<0.0) {
      O2SCL_ERR("Bins negative in prob_dens_hist::init().",exc_efailed);
    }
    range[i]=h.get_bin_low_i(i);
  }
  range[n]=h.get_bin_high_i(n-1);

  double mean=0.0;
  double dsum=0.0;
  for(size_t i=0;i<n;i++) {
    mean+=(h[i]-mean)/((double)(i+1));
  }
  sum[0]=0.0;
  for(size_t i=0;i<n;i++) {
    dsum+=(h[i]/mean)/n;
    sum[i+1]=dsum;
  }
  sv.set_vec(n+1,sum);
  
  return;
}

double prob_dens_hist::operator()() const {
  
  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::operator().",
	      exc_efailed);
  }
  
  double r=rng.random();
  if (r==1.0) r=0.0;
  size_t cache=0;
  size_t ix=sv.find_inc_const(r,cache);
  
  double delta=(r-sum[ix])/(sum[ix+1]-sum[ix]);
  double x=range[ix]+delta*(range[ix+1]-range[ix]);
  
  return x;
}

double prob_dens_hist::log_pdf(double r) const {
  return log(pdf(r));
}

double prob_dens_hist::pdf(double r) const {

  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::pdf().",exc_efailed);
  }
  
  if (r<range[0]) return 0.0;
  if (r>=range[n]) return 0.0;

  return 0.0;
}
    
double prob_dens_hist::cdf(double r) const {

  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::cdf().",exc_efailed);
  }
  
  if (r<range[0]) return 0.0;
  if (r>=range[n]) return 1.0;

  return 0.0;
}

double prob_dens_hist::invert_cdf(double r) const {

  if (n==0) {
    O2SCL_ERR("No data specified in prob_dens_hist::cdf().",exc_efailed);
  }
  
  if (r<range[0]) return 0.0;
  if (r>=range[n]) return 1.0;

  return 0.0;
}

