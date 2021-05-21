/*
  -------------------------------------------------------------------
  
  Copyright (C) 2021-2021, Andrew W. Steiner and Jesse Farr
  
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

#include <o2scl/vec_stats.h>

using namespace std;
using namespace o2scl;

double o2scl::kl_div_gaussian(double mean_prior, double mean_post,
                          double covar_prior, double covar_post) {
  
  double covar_prior_inv=1/covar_prior;
  
  double prod1=covar_prior_inv*covar_post;
  
  double diff=mean_prior-mean_post;
  
  double prod2=diff*covar_prior_inv*diff;
  
  double div=0.5*(prod1+prod2-1.0+log(covar_prior/covar_post));
  
  return div;
}
