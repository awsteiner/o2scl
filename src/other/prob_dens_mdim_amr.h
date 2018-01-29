/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018, Andrew W. Steiner
  
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
/** \file prob_dens_mdim_amr.h
    \brief File defining \ref o2scl::matrix_view, 
    \ref o2scl::matrix_view_table, and \ref o2scl::prob_dens_mdim_amr
*/
#ifndef O2SCL_PROB_DENS_MDIM_AMR_H
#define O2SCL_PROB_DENS_MDIM_AMR_H

#include <o2scl/table.h>
#include <o2scl/err_hnd.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/rng_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A simple matrix view object
   */
  class matrix_view {
  
  public:
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    const double &operator()(size_t row, size_t col) const;
    /** \brief Return the number of rows
     */
    size_t size1();
    /** \brief Return the number of columns
     */
    size_t size2();
  
  };

  /** \brief View a o2scl::table object as a matrix

      \note This stores a pointer to the table and the user must ensure
      that the pointer is valid with the matrix view is accessed.
  */
  template<class vec_t> 
    class matrix_view_table : public matrix_view {
  
  protected:
  
    /// The number of columns
    size_t nc;
    /// Pointers to each column
    std::vector<const vec_t *> col_ptrs;
    /// Pointer to the table
    o2scl::table<vec_t> *tp;

  public:

    /** \brief Create a matrix view object from the specified 
	table and list of columns
    */
    matrix_view_table(o2scl::table<vec_t> &t,
		      std::vector<std::string> cols) {
      nc=cols.size();
      col_ptrs.resize(nc);
      for(size_t i=0;i<nc;i++) {
	col_ptrs[i]=&t[cols[i]];
      }
      tp=&t;
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() {
      return tp->get_nlines();
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() {
      return nc;
    }
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    const double &operator()(size_t row, size_t col) const {
      const vec_t *cp=col_ptrs[col];
      return (*cp)[row];
    }
    
  };

  /** \brief Probability distribution from an adaptive mesh
      created using a matrix of points

      \note This class is experimental.
  */
  template<class vec_t, class mat_t=matrix_view_table<vec_t> >
    class prob_dens_mdim_amr : public o2scl::prob_dens_mdim<vec_t> {

  protected:

  /** \brief A hypercube class for \ref o2scl::prob_dens_mdim_amr
   */
  class hypercube {

  public:

    /** \brief The number of dimensions
     */
    size_t ndim;
    /** \brief The corner of smallest values 
     */
    std::vector<double> low;
    /** \brief The corner of largest values 
     */
    std::vector<double> high;
    /** \brief The list of indices inside
     */
    std::vector<size_t> inside;
    /** \brief The fractional volume enclosed
     */
    double frac_vol;
    /** \brief The weight 
     */
    double weight;

    /** \brief Create an empty hypercube
     */
    hypercube() {
      ndim=0;
    }

    /** \brief Set the hypercube information
     */
    template<class vec2_t>
      void set(vec2_t &l, vec2_t &h, size_t in, double fvol, double wgt) {
      ndim=l.size();
      low.resize(l.size());
      high.resize(h.size());
      inside.resize(1);
      inside[0]=in;
      for(size_t i=0;i<l.size();i++) {
	if (low[i]>high[i]) {
	  low[i]=h[i];
	  high[i]=l[i];
	} else {
	  low[i]=l[i];
	  high[i]=h[i];
	}
      }
      frac_vol=fvol;
      weight=wgt;
      return;
    }
  
    /** \brief Copy constructor
     */
    hypercube(const hypercube &h) {
      ndim=h.ndim;
      low=h.low;
      high=h.high;
      inside=h.inside;
      frac_vol=h.frac_vol;
      weight=h.weight;
      return;
    }

    /** \brief Copy constructor through <tt>operator=()</tt>
     */
    hypercube &operator=(const hypercube &h) {
      if (this!=&h) {
	ndim=h.ndim;
	low=h.low;
	high=h.high;
	inside=h.inside;
	frac_vol=h.frac_vol;
	weight=h.weight;
      }
      return *this;
    }
  
    /** \brief Test if point \c v is inside this hypercube
     */
    template<class vec2_t> bool is_inside(const vec2_t &v) const {
      for(size_t i=0;i<ndim;i++) {
	if (high[i]<v[i] || v[i]<low[i]) {
	  return false;
	}
      }
      return true;
    }
  
  };

  public:

  /** \brief Internal random number generator
   */
  mutable o2scl::rng_gsl rg;
 
  /** \brief Number of dimensions
   */
  size_t ndim;

  /** \brief Corner of smallest values
   */
  vec_t low;
  
  /** \brief Corner of largest values
   */
  vec_t high;

  /** \brief Vector of length scales
   */
  vec_t scale;

  /** \brief Mesh stored as an array of hypercubes
   */
  std::vector<hypercube> mesh;

  /** \brief Verbosity parameter
   */
  int verbose;

  prob_dens_mdim_amr() {
    ndim=0;
  }
  
  /** \brief Initialize a probability distribution from the corners
   */
  prob_dens_mdim_amr(vec_t &l, vec_t &h, vec_t &s) {
    set(l,h,s);
  }

  /** \brief Set the mesh limits

      This function is called by the constructor
   */
  void set(vec_t &l, vec_t &h, vec_t &s) {
    mesh.clear();
    if (h.size()<l.size() || s.size()<l.size()) {
      O2SCL_ERR2("Vector sizes not correct in ",
		"prob_dens_mdim_amr::set().",o2scl::exc_einval);
    }
    low.resize(l.size());
    high.resize(h.size());
    scale.resize(s.size());
    for(size_t i=0;i<l.size();i++) {
      low[i]=l[i];
      high[i]=h[i];
      scale[i]=s[i];
      if (s[i]<=0.0) {
	O2SCL_ERR2("Scale parameter zero or negative in ",
		   "prob_dens_mdim_amr::set().",o2scl::exc_einval);
      }
    }
    ndim=low.size();
    verbose=1;
    return;
  }
 
  /** \brief Insert point at row \c ir, creating a new hypercube 
      for the new point
   */
  void insert(size_t ir, mat_t &m) {
    if (ndim==0) {
      O2SCL_ERR2("Region limits and scales not set in ",
		 "prob_dens_mdim_amr::insert().",o2scl::exc_einval);
    }

    if (mesh.size()==0) {
      if (verbose>1) {
	std::cout << "Creating cube with point ";
	for(size_t k=0;k<ndim;k++) {
	  std::cout << m(0,k) << " ";
	}
	std::cout << std::endl;
      }
      
      // Initialize the mesh with the first point
      mesh.resize(1);
      mesh[0].set(low,high,0,1.0,m(0,ndim));
      return;
    }
   
    // Convert the row to a vector
    std::vector<double> v;
    for(size_t k=0;k<ndim;k++) {
      v.push_back(m(ir,k));
    }
    if (verbose>1) {
      std::cout << "Finding cube with point ";
      for(size_t k=0;k<ndim;k++) {
	std::cout << v[k] << " ";
      }
      std::cout << std::endl;
    }
   
    // Find the right hypercube
    bool found=false;
    size_t jm=0;
    for(size_t j=0;j<mesh.size() && found==false;j++) {
      if (mesh[j].is_inside(v)) {
	found=true;
	jm=j;
      }
    }
    if (found==false) {
      O2SCL_ERR2("Couldn't find point inside mesh in ",
		 "prob_dens_mdim_amr::insert().",o2scl::exc_efailed);
    }
    hypercube &h=mesh[jm];
    if (verbose>1) {
      std::cout << "Found cube " << jm << std::endl;
    }
   
    // Find coordinate with largest relative variation
    size_t max_ip=0;
    double max_var=fabs(v[0]-m(h.inside[0],0))/scale[0];
    for(size_t ip=1;ip<ndim;ip++) {
      double var=fabs(v[ip]-m(h.inside[0],ip))/scale[ip];
	scale[i];
      if (var>max_var) {
	max_ip=ip;
	max_var=var;
      }
    }
    if (verbose>1) {
      std::cout << "Found coordinate " << max_ip << " with variance "
      << max_var << std::endl;
      std::cout << v[max_ip] << " " << m(h.inside[0],max_ip) << " "
      << h.high[max_ip] << " " << h.low[max_ip] << std::endl;
    }
   
    // Slice the mesh in coordinate max_ip
    double loc=(v[max_ip]+m(h.inside[0],max_ip))/2.0;
    double old_vol=h.frac_vol;
    double old_low=h.low[max_ip];
    double old_high=h.high[max_ip];

    size_t old_inside=h.inside[0];

    if (verbose>1) {
      std::cout << "Old limits:" << std::endl;
      for(size_t i=0;i<ndim;i++) {
	std::cout << h.low[i] << " " << h.high[i] << std::endl;
      }
    }
    
    // Set values for hypercube currently in mesh
    h.low[max_ip]=loc;
    h.high[max_ip]=old_high;
    h.frac_vol=old_vol*(old_high-loc)/(old_high-old_low);
   
    // Set values for new hypercube
    hypercube h_new;
    std::vector<double> low_new, high_new;
    o2scl::vector_copy(h.low,low_new);
    o2scl::vector_copy(h.high,high_new);
    low_new[max_ip]=old_low;
    high_new[max_ip]=loc;
    double new_vol=old_vol*(loc-old_low)/(old_high-old_low);
    h_new.set(low_new,high_new,ir,new_vol,m(ir,ndim));

    // --------------------------------------------------------------
    // Todo: this test is unnecessarily slow, and can be replaced by a
    // simple comparison between v[max_ip], old_low, old_high, and
    // m(h.inside[0],max_ip)
    
    if (h.is_inside(v)) {
      h.inside[0]=ir;
      h_new.inside[0]=old_inside;
    } else {
      h.inside[0]=old_inside;
      h_new.inside[0]=ir;
    }

    // --------------------------------------------------------------
    
    if (verbose>1) {
      std::cout << "New limits:" << std::endl;
      for(size_t i=0;i<ndim;i++) {
	std::cout << h.low[i] << " " << h.high[i] << std::endl;
      }
      std::cout << "New cube " << mesh.size() << std::endl;
      for(size_t i=0;i<ndim;i++) {
	std::cout << h_new.low[i] << " " << h_new.high[i] << std::endl;
      }
    }

    // Add new hypercube to mesh
    mesh.push_back(h_new);
   
    return;
  }
 
  /** \brief Parse the matrix \c m, creating a new hypercube
      for every point 
   */
  void initial_parse(mat_t &m) {
   
    for(size_t ir=0;ir<m.size1();ir++) {
      insert(ir,m);
    }
   
    return;
  }

  /** \brief Check the total volume by adding up the fractional
      part of the volume in each hypercube
   */
  double total_volume() {
    if (mesh.size()==0) {
      O2SCL_ERR2("Mesh empty in ",
		 "prob_dens_mdim_amr::total_volume().",o2scl::exc_einval);
    }
    double ret=0.0;
    for(size_t i=0;i<mesh.size();i++) {
      ret+=mesh[i].frac_vol;
    }
    return ret;
  }
 
  /// The normalized density 
  virtual double pdf(const vec_t &x) const {

    if (mesh.size()==0) {
      O2SCL_ERR2("Mesh empty in ",
		 "prob_dens_mdim_amr::pdf().",o2scl::exc_einval);
    }

    // Find the right hypercube
    bool found=false;
    size_t jm=0;
    for(size_t j=0;j<mesh.size() && found==false;j++) {
      if (mesh[j].is_inside(x)) {
	found=true;
	jm=j;
      }
    }
    if (found==false) {
      O2SCL_ERR("Error 2.",o2scl::exc_esanity);
    }
    return mesh[jm].weight;
  }

  /// Select a random point in the largest weighted box
  virtual void select_in_largest(vec_t &x) const {
   
    if (mesh.size()==0) {
      O2SCL_ERR2("Mesh empty in ",
		 "prob_dens_mdim_amr::operator()().",o2scl::exc_einval);
    }

    size_t im=0;
    double wgt=mesh[0].frac_vol*mesh[0].weight;
    for(size_t i=1;i<mesh.size();i++) {
      if (mesh[i].frac_vol*mesh[i].weight>wgt) {
	im=i;
	wgt=mesh[i].frac_vol*mesh[i].weight;
      }
    }
    for(size_t j=0;j<ndim;j++) {
      x[j]=rg.random()*(mesh[im].high[j]-mesh[im].low[j])+mesh[im].low[j];
    }

    return;
  }

  /// Sample the distribution
  virtual void operator()(vec_t &x) const {
   
    if (mesh.size()==0) {
      O2SCL_ERR2("Mesh empty in ",
		 "prob_dens_mdim_amr::operator()().",o2scl::exc_einval);
    }

    double total_weight=0.0;
    for(size_t i=0;i<mesh.size();i++) {
      total_weight+=mesh[i].weight*mesh[i].frac_vol;
    }
   
    double this_weight=rg.random()*total_weight;
    double cml_wgt=0.0;
    for(size_t j=0;j<mesh.size();j++) {
      cml_wgt+=mesh[j].frac_vol*mesh[j].weight;
      if (this_weight<cml_wgt || j==mesh.size()-1) {
	for(size_t i=0;i<ndim;i++) {
	  x[i]=mesh[j].low[i]+rg.random()*
	    (mesh[j].high[i]-mesh[j].low[i]);
	}
	return;
      }
    }

    return;
  }
 
  };
 
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
