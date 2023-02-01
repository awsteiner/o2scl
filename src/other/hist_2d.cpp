/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2023, Andrew W. Steiner
  
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

#include <o2scl/hist_2d.h>

using namespace std;
using namespace o2scl;

hist_2d::hist_2d() {
  xrmode=rmode_avg;
  yrmode=rmode_avg;
  extend_rhs=false;
  extend_lhs=false;
  hsize_x=0;
  hsize_y=0;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

hist_2d::~hist_2d() {
  // Don't call is_valid() here, because we want
  // destructors not to throw exceptions
}

hist_2d::hist_2d(const hist_2d &h) {
  xrmode=h.xrmode;
  yrmode=h.yrmode;
  extend_rhs=h.extend_rhs;
  hsize_x=h.hsize_x;
  hsize_y=h.hsize_y;
  xa=h.xa;
  ya=h.ya;
  xrep=h.xrep;
  yrep=h.yrep;
  user_xrep=h.user_xrep;
  user_yrep=h.user_yrep;
  wgt=h.wgt;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

hist_2d &hist_2d::operator=(const hist_2d &h) {
  if (&h!=this) {
    xrmode=h.xrmode;
    yrmode=h.yrmode;
    extend_rhs=h.extend_rhs;
    hsize_x=h.hsize_x;
    hsize_y=h.hsize_y;
    xa=h.xa;
    ya=h.ya;
    xrep=h.xrep;
    yrep=h.yrep;
    user_xrep=h.user_xrep;
    user_yrep=h.user_yrep;
    wgt=h.wgt;
  }
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  return *this;
}

void hist_2d::set_bin_edges(uniform_grid<double> gx, uniform_grid<double> gy) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (gx.get_nbins()==0 || gy.get_nbins()==0) {
    O2SCL_ERR2("Requested zero bins in ",
	       "hist_2d::set_bin_edges().",exc_einval);
  }
  // Allocate if necessary
  if (gx.get_nbins()!=hsize_x || gy.get_nbins()!=hsize_y) {
    if (hsize_x!=0 || hsize_y!=0) {
      O2SCL_ERR2("Requested binning change in non-empty ",
		 "histogram in hist_2d::set_bin_edges().",exc_efailed);
    }
    allocate(gx.get_nbins(),gy.get_nbins());
  }
  // Set the rep modes from the uniform_grid objects
  if (xrmode!=rmode_user) {
    if (gx.is_log()) {
      xrmode=rmode_gmean;
    } else {
      xrmode=rmode_avg;
    }
  }
  if (yrmode!=rmode_user) {
    if (gy.is_log()) {
      yrmode=rmode_gmean;
    } else {
      yrmode=rmode_avg;
    }
  }
  // Set the uniform_grid
  gx.vector(xa);
  gy.vector(ya);

  // Reset internal reps
  if (xrep.size()>0) xrep.resize(0);
  if (yrep.size()>0) yrep.resize(0);

  return;
}

void hist_2d::allocate(size_t nx, size_t ny) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (nx==0 || ny==0) {
    O2SCL_ERR("Tried to allocate zero bins in hist_2d::allocate().",
	      exc_efailed);
  }
  xa.resize(nx+1);
  ya.resize(ny+1);
  xrep.resize(nx);
  yrep.resize(ny);
  wgt.resize(nx,ny);
  hsize_x=nx;
  hsize_y=ny;

  // Set all weights to zero
  for(size_t i=0;i<nx;i++) {
    for(size_t j=0;j<ny;j++) {
      wgt(i,j)=0.0;
    }
  }

  xrmode=rmode_avg;
  yrmode=rmode_avg;
  return;
}

double hist_2d::sum_wgts() {
  double sum=0.0;
  for(size_t i=0;i<hsize_x;i++) {
    for(size_t j=0;j<hsize_y;j++) {
      sum+=wgt(i,j);
    }
  }
  return sum;
}

double hist_2d::integ_wgts() {
  double sum=0.0;
  for(size_t i=0;i<hsize_x;i++) {
    for(size_t j=0;j<hsize_y;j++) {
      sum+=wgt(i,j)*(xa[i+1]-xa[i])*(ya[j+1]-ya[j]);
    }
  }
  return sum;
}

void hist_2d::set_reps_auto() {
  if (xrmode!=rmode_user) {
    xrep.resize(hsize_x);
    for(size_t i=0;i<hsize_x;i++) {
      if (xrmode==rmode_avg) xrep[i]=(xa[i]+xa[i+1])/2.0;
      else if (xrmode==rmode_low) xrep[i]=xa[i];
      else if (xrmode==rmode_high) xrep[i]=ya[i];
      else if (xrmode==rmode_gmean) xrep[i]=sqrt((xa[i]*xa[i+1]));
      else {
	O2SCL_ERR("Invalid xrmode in hist_2d::set_reps_auto().",
		  exc_efailed);
      }
    }
  }
  if (yrmode!=rmode_user) {
    yrep.resize(hsize_y);
    for(size_t i=0;i<hsize_y;i++) {
      if (yrmode==rmode_avg) yrep[i]=(ya[i]+ya[i+1])/2.0;
      else if (yrmode==rmode_low) yrep[i]=ya[i];
      else if (yrmode==rmode_high) yrep[i]=ya[i];
      else if (yrmode==rmode_gmean) yrep[i]=sqrt((ya[i]*ya[i+1]));
      else {
	O2SCL_ERR("Invalid yrmode in hist_2d::set_reps_auto().",
		  exc_efailed);
      }
    }
  }
  return;
}

void hist_2d::clear_wgts() {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (hsize_x>0 || hsize_y>0) {
    for(size_t i=0;i<wgt.size1();i++) {
      for(size_t j=0;j<wgt.size2();j++) {
	wgt(i,j)=0.0;
      }
    }
  }
  return;
}

void hist_2d::clear() {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (hsize_x>0 || hsize_y>0) {
    xa.resize(0);
    ya.resize(0);
    user_xrep.resize(0);
    xrep.resize(0);
    user_yrep.resize(0);
    yrep.resize(0);
    wgt.resize(0,0);
    hsize_x=0;
    hsize_y=0;
  }
  return;
}

void hist_2d::get_bin_indices(double x, double y, 
			      size_t &i, size_t &j) const {
  i=get_x_bin_index(x);
  j=get_y_bin_index(y);
  return;
}

size_t hist_2d::get_x_bin_index(double x) const {
  size_t i;

  if (hsize_x==0) {
    O2SCL_ERR2("Histogram has zero size in ",
	      "hist::get_x_bin_index().",exc_einval);
  }

  // Compute x index
  search_vec<const ubvector> sv(xa.size(),xa);
  if (xa[0]<xa[hsize_x]) {
    if (x<xa[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(x)+"' smaller than smallest "+
	  "bin '"+dtos(xa[0])+"' in hist_2d::get_x_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (x>xa[hsize_x]) {
      if (extend_rhs) {
	return hsize_x-1;
      } else {
	std::string s="Value '"+dtos(x)+"' larger than largest "+
	  "bin '"+dtos(xa[hsize_x])+"' in hist_2d::get_x_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    i=sv.find_inc(x);
  } else {
    if (x>xa[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(x)+"' larger than largest "+
	  "bin '"+dtos(xa[0])+"' in hist_2d::get_x_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (x<xa[hsize_x]) {
      if (extend_rhs) {
	return hsize_x-1;
      } else {
	std::string s="Value '"+dtos(x)+"' smaller than smallest "+
	  "bin '"+dtos(xa[hsize_x])+"' in hist_2d::get_x_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    i=sv.find_dec(x);
  }

  return i;
}

size_t hist_2d::get_y_bin_index(double y) const {
  size_t j;

  if (hsize_y==0) {
    O2SCL_ERR2("Histogram has zero size in ",
	      "hist::get_y_bin_index().",exc_einval);
  }

  // Compute y index
  search_vec<const ubvector> sv2(ya.size(),ya);
  if (ya[0]<ya[hsize_y]) {
    if (y<ya[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(y)+"' smaller than smallest "+
	  "bin '"+dtos(ya[0])+"' in hist_2d::get_y_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (y>ya[hsize_y]) {
      if (extend_rhs) {
	return hsize_y-1;
      } else {
	std::string s="Value '"+dtos(y)+"' larger than largest "+
	  "bin '"+dtos(ya[hsize_y])+"' in hist_2d::get_y_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    j=sv2.find_inc(y);
  } else {
    if (y>ya[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(y)+"' larger than largest "+
	  "bin '"+dtos(ya[0])+"' in hist_2d::get_y_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (y<ya[hsize_y]) {
      if (extend_rhs) {
	return hsize_y-1;
      } else {
	std::string s="Value '"+dtos(y)+"' smaller than smallest "+
	  "bin '"+dtos(ya[hsize_y])+"' in hist_2d::get_y_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    j=sv2.find_dec(y);
  }

  return j;
}

double &hist_2d::get_x_low_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize_x) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize_x)+
      " in hist_2d::get_x_low().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return xa[0];
  }
  if (xrep.size()>0) xrep.resize(0);
  return xa[i];
}

const double &hist_2d::get_x_low_i(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize_x) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize_x)+
      " in hist_2d::get_x_low().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return xa[0];
  }
  return xa[i];
}

double &hist_2d::get_x_high_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize_x) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize_x)+
      " in hist_2d::get_x_high().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return xa[0];
  }
  if (xrep.size()>0) xrep.resize(0);
  return xa[i+1];
}

const double &hist_2d::get_x_high_i(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize_x) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize_x)+
      " in hist_2d::get_x_high().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return xa[0];
  }
  return xa[i+1];
}

double &hist_2d::get_y_low_i(size_t j) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (j>=hsize_y) {
    std::string s="Index "+itos(j)+" must be smaller than "+
      "the histogram size "+itos(hsize_y)+
      " in hist_2d::get_y_low().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ya[0];
  }
  if (yrep.size()>0) yrep.resize(0);
  return ya[j];
}

const double &hist_2d::get_y_low_i(size_t j) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (j>=hsize_y) {
    std::string s="Index "+itos(j)+" must be smaller than "+
      "the histogram size "+itos(hsize_y)+
      " in hist_2d::get_y_low().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ya[0];
  }
  return ya[j];
}

double &hist_2d::get_y_high_i(size_t j) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (j>=hsize_y) {
    std::string s="Index "+itos(j)+" must be smaller than "+
      "the histogram size "+itos(hsize_y)+
      " in hist_2d::get_y_high().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ya[0];
  }
  if (yrep.size()>0) yrep.resize(0);
  return ya[j+1];
}

const double &hist_2d::get_y_high_i(size_t j) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (j>=hsize_y) {
    std::string s="Index "+itos(j)+" must be smaller than "+
      "the histogram size "+itos(hsize_y)+
      " in hist_2d::get_y_high().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ya[0];
  }
  return ya[j+1];
}

void hist_2d::set_rep_mode(size_t x_mode, size_t y_mode) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (x_mode>4 || y_mode>4) {
    O2SCL_ERR2("Valid rep mode values are 0-4 in ",
	       "hist_2d::set_rep_mode().",exc_einval);
  }
  if (x_mode==rmode_user) {
    if (xrmode==rmode_user) {
      return;
    }
    O2SCL_ERR2("Need to set user rep. mode for x with set_reps() in ",
	       "hist::set_rep_mode().",exc_efailed);
  }
  if (y_mode==rmode_user) {
    if (yrmode==rmode_user) {
      return;
    }
    O2SCL_ERR2("Need to set user rep. mode for y with set_reps() in ",
	       "hist::set_rep_mode().",exc_efailed);
  }
  xrmode=x_mode;
  yrmode=y_mode;
  if (xrep.size()>0) xrep.resize(0);
  if (yrep.size()>0) yrep.resize(0);
  return;
}

double hist_2d::get_x_rep_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize_x) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize_x)+
      " in hist_2d::get_x().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return xrep[0];
  }
  if (xrmode==rmode_user) return user_xrep[i];
  if (xrep.size()==0) set_reps_auto();
  return xrep[i];
}

double hist_2d::get_y_rep_i(size_t j) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (j>=hsize_y) {
    std::string s="Index "+itos(j)+" must be smaller than "+
      "the histogram size "+itos(hsize_y)+
      " in hist_2d::get_y().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return yrep[0];
  }
  if (yrmode==rmode_user) return user_yrep[j];
  if (yrep.size()==0) set_reps_auto();
  return yrep[j];
}

const double &hist_2d::get_wgt_i(size_t i, size_t j) const {
  if (i>=hsize_x || j>=hsize_y) {
    std::string s="Indices "+itos(i)+" and "+itos(j)+" must be smaller "+
      "than the histogram sizes "+itos(hsize_x)+" and "+itos(hsize_y)+
      " in const double &hist_2d::get_wgt().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return wgt(0,0);
  }
  return wgt(i,j);
}

double &hist_2d::get_wgt_i(size_t i, size_t j) {
  if (i>=hsize_x || j>=hsize_y) {
    std::string s="Indices "+itos(i)+" and "+itos(j)+" must be smaller "+
      "than the histogram sizes "+itos(hsize_x)+" and "+itos(hsize_y)+
      " in double &hist_2d::get_wgt_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return wgt(0,0);
  }
  return wgt(i,j);
}

void hist_2d::set_wgt_i(size_t i, size_t j, double val) {
  if (i>=hsize_x || j>=hsize_y) {
    std::string s="Indices "+itos(i)+" and "+itos(j)+" must be smaller "+
      "than the histogram sizes "+itos(hsize_x)+" and "+itos(hsize_y)+
      " in double &hist_2d::set_wgt_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
  }
  wgt(i,j)=val;
  return;
}

void hist_2d::is_valid() const {
  if ((hsize_x==0 && hsize_y>0) || 
      (hsize_y==0 && hsize_x>0)) {
    std::string str=((std::string)"Semi-empty histogram (size_x=")+
      o2scl::szttos(hsize_x)+", size_y="+o2scl::szttos(hsize_x)+
      ") in hist_2d::is_valid().";
    O2SCL_ERR(str.c_str(),exc_efailed);
  }
  if (hsize_x==0) {
    if (xa.size()>0 || ya.size()>0 || wgt.size1()>0 || wgt.size2()>0 || 
	xrep.size()>0 || yrep.size()>0 || user_xrep.size()>0 || 
	user_yrep.size()>0) {
      std::string str=((std::string)"Histogram size is zero but ")+
	"vectors are not empty:"+
	" xa.size()="+o2scl::szttos(xa.size())+
	", ya.size()="+o2scl::szttos(ya.size())+
	", wgt.size1()="+o2scl::szttos(wgt.size1())+
	", wgt.size2()="+o2scl::szttos(wgt.size2())+
	", xrep.size()="+o2scl::szttos(xrep.size())+
	", yrep.size()="+o2scl::szttos(yrep.size())+
	" in hist_2d::is_valid().";
      O2SCL_ERR(str.c_str(),exc_efailed);
    }
  } else {
    if (xa.size()!=hsize_x+1 || ya.size()!=hsize_y+1 || 
	wgt.size1()!=hsize_x || wgt.size2()!=hsize_y ||
	(xrep.size()>0 && xrep.size()!=hsize_x) ||
	(yrep.size()>0 && yrep.size()!=hsize_y)) {
      std::string str=((std::string)"Vector/matrix sizes do not match ")+
	"histogram sizes:"+
	" hsize_x="+o2scl::szttos(hsize_x)+
	", hsize_y="+o2scl::szttos(hsize_y)+
	", xa.size()="+o2scl::szttos(xa.size())+
	", ya.size()="+o2scl::szttos(ya.size())+
	", wgt.size1()="+o2scl::szttos(wgt.size1())+
	", wgt.size2()="+o2scl::szttos(wgt.size2())+
	", xrep.size()="+o2scl::szttos(xrep.size())+
	", yrep.size()="+o2scl::szttos(yrep.size())+
	" in hist_2d::is_valid().";
      O2SCL_ERR(str.c_str(),exc_efailed);
    }
    if (xrmode!=rmode_user && user_xrep.size()>0) {
      O2SCL_ERR2("Rep. mode for x is not user mode but user rep vector ",
		 "is not empty in hist_2d::is_valid().",exc_efailed);
    }
    if (yrmode!=rmode_user && user_yrep.size()>0) {
      O2SCL_ERR2("Rep. mode for y is not user mode but user rep vector ",
		 "is not empty in hist_2d::is_valid().",exc_efailed);
    }
    if (xrmode==rmode_user && user_xrep.size()!=hsize_x) {
      O2SCL_ERR2("Rep. mode for x is user mode but user rep vector ",
		 "size doesn't match in hist_2d::is_valid().",exc_efailed);
    }
    if (yrmode==rmode_user && user_yrep.size()!=hsize_y) {
      O2SCL_ERR2("Rep. mode for y is user mode but user rep vector ",
		 "size doesn't match in hist_2d::is_valid().",exc_efailed);
    }
  }
}

void hist_2d::copy_to_table3d(table3d &t, std::string xreps_name,
			      std::string yreps_name, std::string weights) {
  
  t.clear();

  // Set the grid
  ubvector xreps(size_x()), yreps(size_y());
  for(size_t i=0;i<size_x();i++) {
    xreps[i]=get_x_rep_i(i);
  }
  for(size_t i=0;i<size_y();i++) {
    yreps[i]=get_y_rep_i(i);
  }
  t.set_xy(xreps_name,size_x(),xreps,yreps_name,size_y(),yreps);

  t.new_slice(weights);

  for(size_t i=0;i<size_x();i++) {
    for(size_t j=0;j<size_y();j++) {
      t.set(i,j,weights,get_wgt_i(i,j));
    }
  }
  
  return;
}
