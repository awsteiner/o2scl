/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2016, Andrew W. Steiner
  
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

#include <o2scl/hist.h>

using namespace std;
using namespace o2scl;

hist::hist() {
  itype=1;
  rmode=rmode_avg;
  extend_lhs=false;
  extend_rhs=false;
  hsize=0;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

hist::~hist() {
  // Don't call is_valid() here, because we want
  // destructors not to throw exceptions
}

hist::hist(const hist &h) {
  itype=h.itype;
  rmode=h.rmode;
  extend_rhs=h.extend_rhs;
  extend_lhs=h.extend_lhs;
  hsize=h.hsize;
  ubin=h.ubin;
  urep=h.urep;
  uwgt=h.uwgt;
  user_rep=h.user_rep;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

hist &hist::operator=(const hist &h) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (&h!=this) {
    itype=h.itype;
    rmode=h.rmode;
    extend_rhs=h.extend_rhs;
    extend_lhs=h.extend_lhs;
    hsize=h.hsize;
    ubin=h.ubin;
    urep=h.urep;
    uwgt=h.uwgt;
    user_rep=h.user_rep;
  }
  return *this;
}

void hist::set_reps_auto() {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    O2SCL_ERR2("Function set_reps_auto() was called in ",
	       "user rep mode.",exc_esanity);
  }
  urep.resize(hsize);
  for(size_t i=0;i<hsize;i++) {
    if (rmode==rmode_avg) urep[i]=(ubin[i]+ubin[i+1])/2.0;
    else if (rmode==rmode_low) urep[i]=ubin[i];
    else if (rmode==rmode_high) urep[i]=uwgt[i];
    else if (rmode==rmode_gmean) urep[i]=sqrt((ubin[i]*ubin[i+1]));
    else {
      O2SCL_ERR("Invalid rmode in hist::set_reps().",
		exc_esanity);
    }
  }
  return;
}

void hist::normalize(double new_sum) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (hsize==0) {
    O2SCL_ERR("Empty histogram in hist::normalize().",exc_einval);
  }
  double fac=new_sum/integ(ubin[0],ubin[hsize]);
  for(size_t i=0;i<size();i++) {
    uwgt[i]*=fac;
  }
  return;
}

void hist::allocate(size_t n) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (n==0) {
    O2SCL_ERR("Tried to allocate zero bins in hist::allocate().",
	      exc_efailed);
  }
  ubin.resize(n+1);
  uwgt.resize(n);
  hsize=n;

  // Set all weights to zero
  for(size_t i=0;i<n;i++) uwgt[i]=0.0;

  rmode=rmode_avg;
}

void hist::set_bin_edges(uniform_grid<double> g) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (g.get_nbins()==0) {
    O2SCL_ERR2("Requested zero bins in ",
		   "hist::set_bin_edges().",exc_einval);
  }
  // Allocate if necessary
  if (g.get_nbins()!=hsize) {
    if (hsize!=0) {
      O2SCL_ERR2("Requested binning change in non-empty ",
		     "histogram in hist::set_bin_edges().",exc_efailed);
    }
    allocate(g.get_nbins());
  }
  // Set the rep mode from the uniform_grid if we're not in user mode
  if (rmode!=rmode_user) {
    if (g.is_log()) {
      rmode=rmode_gmean;
    } else {
      rmode=rmode_avg;
    }
  }
  // Set the bins from the uniform grid
  g.vector(ubin);
  // Reset internal reps
  if (urep.size()>0) urep.resize(0);
  return;
}

void hist::clear_wgts() {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (hsize>0) {
    for(size_t i=0;i<uwgt.size();i++) uwgt[i]=0.0;
  }
  return;
}

void hist::clear() {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (hsize>0) {
    ubin.resize(0);
    uwgt.resize(0);
    if (urep.size()>0) urep.resize(0);
    if (user_rep.size()>0) user_rep.resize(0);
    hsize=0;
  }
  return;
}

size_t hist::get_bin_index(double x) const {
  if (hsize==0) {
    O2SCL_ERR2("Histogram has zero size in ",
	      "hist::get_bin_index().",exc_einval);
  }
  search_vec<const ubvector> sv(ubin.size(),ubin);
  // Increasing case
  if (ubin[0]<ubin[hsize]) {
    if (x<ubin[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(x)+"' smaller than smallest "+
	  "bin '"+dtos(ubin[0])+"' (increasing) in hist::get_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (x>ubin[hsize]) {
      if (extend_rhs) {
	return hsize-1;
      } else {
	std::string s="Value '"+dtos(x)+"' larger than largest "+
	  "bin '"+dtos(ubin[hsize])+"' (increasing) in hist::get_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    return sv.find_inc(x);
  } else {
    // Decreasing case
    if (x>ubin[0]) {
      if (extend_lhs) {
	return 0;
      } else {
	std::string s="Value '"+dtos(x)+"' larger than largest "+
	  "bin '"+dtos(ubin[0])+"' (decreasing) in hist::get_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    if (x<ubin[hsize]) {
      if (extend_rhs) {
	return hsize-1;
      } else {
	std::string s="Value '"+dtos(x)+"' smaller than smallest "+
	  "bin '"+dtos(ubin[hsize])+"' (decreasing) in hist::get_bin_index().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
    }
    return sv.find_dec(x);
  }
}

double &hist::get_bin_low_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_bin_low_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ubin[0];
  }
  if (urep.size()>0) urep.resize(0);
  return ubin[i];
}

const double &hist::get_bin_low_i(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_bin_low_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ubin[0];
  }
  return ubin[i];
}

double &hist::get_bin_high_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_bin_high_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ubin[0];
  }
  if (urep.size()>0) urep.resize(0);
  return ubin[i+1];
}

const double &hist::get_bin_high_i(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_bin_high_i().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return ubin[0];
  }
  return ubin[i+1];
}

void hist::set_rep_mode(size_t mode) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (mode>4) {
    O2SCL_ERR2("Valid rep mode values are 0-4 in ",
		   "hist::set_rep_mode().",exc_einval);
  }
  if (mode==rmode_user) {
    if (rmode==rmode_user) {
      return;
    }
    O2SCL_ERR2("Need to set user rep. mode with set_reps() in ",
	       "hist::set_rep_mode().",exc_efailed);
  }
  if (rmode!=mode) {
    rmode=mode;
    if (urep.size()>0) urep.resize(0);
  }
  return;
}

double hist::get_rep_i(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+" in hist::get_bin().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return urep[0];
  }
  // If we're in user mode, just return the user value
  if (rmode==rmode_user) return user_rep[i];
  // Otherwise, we have to compute the internal reps
  if (urep.size()==0) set_reps_auto();
  return urep[i];
}

void hist::update(double x, double val) {
  size_t loc=get_bin_index(x);
  uwgt[loc]+=val;
  return;
}

const double &hist::get_wgt_i(size_t i) const {
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_wgt().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return uwgt[0];
  }
  return uwgt[i];
}

double &hist::get_wgt_i(size_t i) {
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::get_wgt().";
    O2SCL_ERR(s.c_str(),exc_einval);
    return uwgt[0];
  }
  return uwgt[i];
}

void hist::set_wgt_i(size_t i, double val) {
  if (i>=hsize) {
    std::string s="Index "+itos(i)+" must be smaller than "+
      "the histogram size "+itos(hsize)+
      " in hist::set().";
    O2SCL_ERR(s.c_str(),exc_einval);
  }
  uwgt[i]=val;
  return;
}

double hist::operator()(double x) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    interp_t si(hsize,user_rep,uwgt,itype);
    return si.eval(x);
  } 
  if (urep.size()==0) set_reps_auto();
  interp_t si(hsize,urep,uwgt,itype);
  return si.eval(x);
}

double hist::interp(double x) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    interp_t si(hsize,user_rep,uwgt,itype);
    return si.eval(x);
  } 
  if (urep.size()==0) set_reps_auto();
  interp_t si(hsize,urep,uwgt,itype);
  return si.eval(x);
}

double hist::deriv(double x) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    interp_t si(hsize,user_rep,uwgt,itype);
    return si.deriv(x);
  } 
  if (urep.size()==0) set_reps_auto();
  interp_t si(hsize,urep,uwgt,itype);
  return si.deriv(x);
}

double hist::deriv2(double x) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    interp_t si(hsize,user_rep,uwgt,itype);
    return si.deriv2(x);
  } 
  if (urep.size()==0) set_reps_auto();
  interp_t si(hsize,urep,uwgt,itype);
  return si.deriv2(x);
}

double hist::integ(double x, double y) {
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  if (rmode==rmode_user) {
    interp_t si(hsize,user_rep,uwgt,itype);
    return si.integ(x,y);
  } 
  if (urep.size()==0) set_reps_auto();
  interp_t si(hsize,urep,uwgt,itype);
  return si.integ(x,y);
}


void hist::set_interp_type(size_t interp_type) {
  itype=interp_type;
  return;
}

double hist::get_max_wgt() const {
  double max=uwgt[0];
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]>max) max=uwgt[i];
  }
  return max;
}

size_t hist::get_max_index() const {
  double max=uwgt[0];
  size_t max_ix=0;
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]>max) {
      max=uwgt[i];
      max_ix=i;
    }
  }
  return max_ix;
}

double hist::get_max_rep() {
  double max=uwgt[0];
  size_t max_ix=0;
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]>max) {
      max=uwgt[i];
      max_ix=i;
    }
  }
  if (urep.size()==0) set_reps_auto();
  return urep[max_ix];
}

double hist::get_min_wgt() const {
  double min=uwgt[0];
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]<min) min=uwgt[i];
  }
  return min;
}

size_t hist::get_min_index() const {
  double min=uwgt[0];
  size_t min_ix=0;
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]<min) {
      min=uwgt[i];
      min_ix=i;
    }
  }
  return min_ix;
}

double hist::get_min_rep() {
  double min=uwgt[0];
  size_t min_ix=0;
  for(size_t i=1;i<hsize;i++) {
    if (uwgt[i]<min) {
      min=uwgt[i];
      min_ix=i;
    }
  }
  if (urep.size()==0) set_reps_auto();
  return urep[min_ix];
}

void hist::is_valid() const {
  if (hsize==0) {
    if (ubin.size()>0 || uwgt.size()>0 || urep.size()>0 || 
	user_rep.size()>0) {
      O2SCL_ERR2("Histogram size is zero but vectors are not ",
		 "empty in hist::is_valid().",exc_efailed);
    }
  } else {
    if (ubin.size()!=hsize+1 || uwgt.size()!=hsize ||
	(urep.size()>0 && urep.size()!=hsize)) {
      O2SCL_ERR2("Vector sizes don't match histogram size ",
		 "in hist::is_valid().",exc_efailed);
    }
    if (rmode!=rmode_user && user_rep.size()>0) {
      O2SCL_ERR2("Rep. mode is not user mode but user rep vector ",
		 "is not empty in hist::is_valid().",exc_efailed);
    }
    if (rmode==rmode_user && user_rep.size()!=hsize) {
      O2SCL_ERR2("Rep. mode is user mode but user rep vector ",
		 "size doesn't match in hist::is_valid().",exc_efailed);
    }
  }
}

void hist::copy_to_table(table<> &t, std::string reps, std::string lower_edges, 
			 std::string upper_edges, std::string weights) {
  
  if (t.get_nlines()<hsize) {
    size_t old_lines=t.get_nlines();
    t.set_nlines(hsize);
    // Clear out new rows for current columns
    for(size_t i=old_lines;i<hsize;i++) {
      for(size_t j=0;j<t.get_ncolumns();j++) {
	t.set(j,i,0.0);
      }
    }
  }
  t.new_column(reps);
  t.new_column(lower_edges);
  t.new_column(upper_edges);
  t.new_column(weights);
  for(size_t i=0;i<hsize;i++) {
    t.set(reps,i,get_rep_i(i));
    t.set(lower_edges,i,get_bin_low_i(i));
    t.set(upper_edges,i,get_bin_high_i(i));
    t.set(weights,i,get_wgt_i(i));
  }
  
  return;
}
