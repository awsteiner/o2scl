/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2025, Andrew W. Steiner

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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_HDF_IO_H
#define O2SCL_HDF_IO_H

/** \file hdf_io.h
    \brief File defining HDF I/O for selected \o2 objects
*/
#include <fnmatch.h>

#include <boost/numeric/ublas/vector.hpp>

#include <regex>

#include <o2scl/hdf_file.h>
#include <o2scl/table.h>
#include <o2scl/table_units.h>
#include <o2scl/hist.h>
#include <o2scl/hist_2d.h>
#include <o2scl/table3d.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/expval.h>
#include <o2scl/contour.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/prob_dens_mdim_amr.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/exp_max.h>
#include <o2scl/funct_to_fp.h>

/** \brief The \o2 namespace for I/O with HDF
 */
namespace o2scl_hdf {

  /** \brief Input a \ref o2scl::prob_dens_mdim_gaussian object from a 
      \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_gaussian,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_gaussian, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_gaussian, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
  */
  template<class vec_t, class mat_t> 
  void hdf_input_n(hdf_file &hf,
                   o2scl::prob_dens_mdim_gaussian<vec_t,mat_t> &p,
                   std::string &name) {

    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_object_by_type("prob_dens_mdim_gaussian",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type prob_dens_mdim_gaussian found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    
    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    vec_t peak;
    mat_t chol, covar_inv;
    double norm;

    hf.getd("norm",norm);
    hf.getd_vec_copy("peak",peak);
    hf.getd_mat_copy("chol",chol);
    hf.getd_mat_copy("covar_inv",covar_inv);

    p.set_alt(peak.size(),peak,chol,covar_inv,norm);

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    return;
  }

  /** \brief Input a \ref o2scl::prob_dens_mdim_gaussian object from a \ref
      hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_gaussian,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_gaussian, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_gaussian, then the
      error handler is called.
  */
  template<class vec_t, class mat_t> 
  void hdf_input(hdf_file &hf,
                 o2scl::prob_dens_mdim_gaussian<vec_t,mat_t> &p,
                 std::string name="") {
    hdf_input_n<vec_t,mat_t>(hf,p,name);
    return;
  }
  
  /** \brief Output a \ref o2scl::prob_dens_mdim_gaussian 
      object to a \ref hdf_file
  */
  template<class vec_t, class mat_t> 
  void hdf_output(hdf_file &hf,
                  o2scl::prob_dens_mdim_gaussian<vec_t,mat_t> &p,
                  std::string name) {
    
    if (hf.has_write_access()==false) {
      O2SCL_ERR2("File not opened with write access in hdf_output",
		 "(hdf_file,prob_dens_mdim_gaussian<>,string).",
		 o2scl::exc_efailed);
    }
    
    // Start group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    
    // Add typename
    hf.sets_fixed("o2scl_type","prob_dens_mdim_gaussian");
    
    const vec_t &peak=p.get_peak();
    const mat_t &chol=p.get_chol();
    const mat_t &covar_inv=p.get_covar_inv();
    double norm=p.get_norm();

    hf.setd("norm",norm);
    hf.setd_vec_copy("peak",peak);
    hf.setd_mat_copy("chol",chol);
    hf.setd_mat_copy("covar_inv",covar_inv);

    // Close prob_dens_mdim_gaussian group
    hf.close_group(group);
    
    // Return location to previous value
    hf.set_current_id(top);
    
    return;
  }
  
  /** \brief Input a \ref o2scl::prob_dens_mdim_gmm object from a 
      \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_gmm,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_gmm, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_gmm, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
  */
  template<class mat_t, class vecp_t> 
  void hdf_input_n(hdf_file &hf,
                   o2scl::prob_dens_mdim_gmm<mat_t,vecp_t> &p,
                   std::string &name) {

    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_object_by_type("prob_dens_mdim_gmm",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type prob_dens_mdim_gmm found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    
    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    size_t n;
    hf.get_szt("n",n);
    p.weights.resize(n);
    hf.getd_vec_copy("weights",p.weights);

    p.pdmg.resize(n);
    for(size_t i=0;i<n;i++) {
      std::string pdmg_name=((std::string)"pdmg")+o2scl::szttos(i);
      hdf_input_n(hf,p.pdmg[i],pdmg_name);
    }

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    return;
  }
  
  /** \brief Input a \ref o2scl::prob_dens_mdim_gmm object from a \ref
      hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_gmm,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_gmm, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_gmm, then the
      error handler is called.
  */
  template<class mat_t, class vecp_t> 
  void hdf_input(hdf_file &hf,
                 o2scl::prob_dens_mdim_gmm<mat_t,vecp_t> &p,
                 std::string name="") {
    hdf_input_n<mat_t,vecp_t>(hf,p,name);
    return;
  }
  
  /** \brief Output a \ref o2scl::prob_dens_mdim_gmm 
      object to a \ref hdf_file
  */
  template<class mat_t, class vecp_t> 
  void hdf_output(hdf_file &hf,
                  o2scl::prob_dens_mdim_gmm<mat_t,vecp_t> &p,
                  std::string name) {
    
    if (hf.has_write_access()==false) {
      O2SCL_ERR2("File not opened with write access in hdf_output",
		 "(hdf_file,prob_dens_mdim_gmm<>,string).",
		 o2scl::exc_efailed);
    }
    
    // Start group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    
    // Add typename
    hf.sets_fixed("o2scl_type","prob_dens_mdim_gmm");
    
    hf.set_szt("n",p.pdmg.size());
    hf.setd_vec_copy("weights",p.weights);
    for(size_t i=0;i<p.pdmg.size();i++) {
      std::string pdmg_name=((std::string)"pdmg")+o2scl::szttos(i);
      hdf_output(hf,p.pdmg[i],pdmg_name);
    }

    // Close prob_dens_mdim_gmm group
    hf.close_group(group);
    
    // Return location to previous value
    hf.set_current_id(top);
    
    return;
  }
  
  /** \brief Input a \ref o2scl::prob_dens_mdim_amr object from a 
      \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_amr,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_amr, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_amr, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
  */
  template<class vec_t, class mat_t> 
  void hdf_input_n(hdf_file &hf,
                   o2scl::prob_dens_mdim_amr<vec_t,mat_t> &p,
                   std::string &name) {

    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_object_by_type("prob_dens_mdim_amr",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type prob_dens_mdim_amr found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    
    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    size_t nd, dc, ms;
    std::vector<double> data;
    std::vector<size_t> insides;
    
    hf.get_szt("n_dim",nd);
    hf.get_szt("dim_choice",dc);
    hf.get_szt("mesh_size",ms);
    hf.getd_vec("data",data);
    hf.get_szt_vec("insides",insides);

    p.set_from_vectors(nd,dc,ms,data,insides);

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    return;
  }

  /** \brief Input a \ref o2scl::prob_dens_mdim_amr object from a \ref
      hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::prob_dens_mdim_amr,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::prob_dens_mdim_amr, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::prob_dens_mdim_amr, then the
      error handler is called. 
  */
  template<class vec_t, class mat_t> 
  void hdf_input(hdf_file &hf,
                 o2scl::prob_dens_mdim_amr<vec_t,mat_t> &p,
                 std::string name="") {
    hdf_input_n<vec_t,mat_t>(hf,p,name);
    return;
  }
  
  /** \brief Output a \ref o2scl::prob_dens_mdim_amr 
      object to a \ref hdf_file
  */
  template<class vec_t, class mat_t> 
  void hdf_output(hdf_file &hf,
                  o2scl::prob_dens_mdim_amr<vec_t,mat_t> &p,
                  std::string name) {
    
    if (hf.has_write_access()==false) {
      O2SCL_ERR2("File not opened with write access in hdf_output",
		 "(hdf_file,prob_dens_mdim_amr<>,string).",
		 o2scl::exc_efailed);
    }
    
    // Start group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    
    // Add typename
    hf.sets_fixed("o2scl_type","prob_dens_mdim_amr");
    
    size_t nd, dc, ms;
    std::vector<double> data;
    std::vector<size_t> insides;
    
    p.copy_to_vectors(nd,dc,ms,data,insides);

    hf.set_szt("n_dim",nd);
    hf.set_szt("dim_choice",dc);
    hf.set_szt("mesh_size",ms);
    hf.setd_vec("data",data);
    hf.set_szt_vec("insides",insides);

    // Close prob_dens_mdim_amr group
    hf.close_group(group);
    
    // Return location to previous value
    hf.set_current_id(top);
    
    return;
  }
  
  /** \brief Output a \ref o2scl::table object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table<> &t, std::string name);

  /** \brief Input a \ref o2scl::table object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::table, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::table, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
  */
  template<class vec_t, class fp_t> 
  void hdf_input_n(hdf_file &hf, o2scl::table<vec_t,fp_t> &t,
                   std::string &name) {

    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_object_by_type("table",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type table found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }

    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    std::string type;
    hf.gets_fixed("o2scl_type",type);
    if (type!="table") {
      O2SCL_ERR("Object is not of type table in hdf_input_n().",
                o2scl::exc_einval);
    }
    
    // Input the table data
    hdf_input_data(hf,t);

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    t.is_valid();

    return;
  }

  /** \brief Input a \ref o2scl::table object from a \ref
      hdf_file

      \note This function is declared in <tt>table.h</tt>, and 
      thus the default argument to the \c name parameter 
      appears there instead of here.

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::table, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::table, then the
      error handler is called. 
  */
  template<class vec_t, class fp_t> 
  void hdf_input(hdf_file &hf, o2scl::table<vec_t,fp_t> &t, std::string name) {
    hdf_input_n<vec_t>(hf,t,name);
    return;
  }
  /** \brief Internal function for outputting a \ref o2scl::table object
   */
  void hdf_output_data(hdf_file &hf, o2scl::table<> &t);

  /** \brief Internal function for inputting a \ref o2scl::table object
   */
  template<class vec_t, class fp_t> 
  void hdf_input_data(hdf_file &hf, o2scl::table<vec_t,fp_t> &t) {
    hid_t group=hf.get_current_id();

    // Clear previous data
    t.clear_table();
    t.clear_constants();

    // Check typename
    std::string type2;
    hf.gets_fixed("o2scl_type",type2);
    if (type2!="table") {
      O2SCL_ERR2("Typename in HDF group does not match ",
		 "class in o2scl_hdf::hdf_input_data().",o2scl::exc_einval);
    }

    // Storage
    std::vector<std::string> cnames, cols;
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector cvalues;
      
    // Get constants
    hf.gets_vec_copy("con_names",cnames);
    hf.getd_vec_copy("con_values",cvalues);
    if (cnames.size()!=cvalues.size()) {
      O2SCL_ERR2("Size mismatch between constant names and values ",
		 "in o2scl_hdf::hdf_input_data().",o2scl::exc_einval);
    }
    for(size_t i=0;i<cnames.size();i++) {
      t.add_constant(cnames[i],cvalues[i]);
    }

    // Get column names
    hf.gets_vec_copy("col_names",cols);
    for(size_t i=0;i<cols.size();i++) {
      t.new_column(cols[i]);
    }

    // Get number of lines
    int nlines2;
    hf.geti("nlines",nlines2);
    t.set_nlines(nlines2);
    
    // Output the interpolation type
    hf.get_szt_def("itype",o2scl::itp_cspline,t.itype);

    // Open data group
    hid_t group2=hf.open_group("data");
    hf.set_current_id(group2);

    if (nlines2>0) {
    
      // Get data
      for(size_t i=0;i<t.get_ncolumns();i++) {
	ubvector vtmp(nlines2);
	hf.getd_vec_copy(t.get_column_name(i),vtmp);
	for(int j=0;j<nlines2;j++) {
	  t.set(t.get_column_name(i),j,vtmp[j]);
	}
      }

    }

    // Close groups
    hf.close_group(group2);

    hf.set_current_id(group);

    // Check that input created a valid table
    t.is_valid();

    return;
  }
  
  /** \brief Output a \ref o2scl::table_units object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table_units<> &t, 
		  std::string name);

  /** \brief Input a \ref o2scl::table_units object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table_units, otherwise, it reads
      the object of the specified name. If \c name is specified and no
      object is found, or if the object with the specified name is not
      of type \ref o2scl::table_units, then the error handler is called. If
      \c name is unspecified and there is no object with type \ref
      o2scl::table_units, then the error handler is called. Upon exit, \c
      name contains the name of the object which was read.
  */
  template<class vec_t> 
  void hdf_input_n(hdf_file &hf, o2scl::table_units<vec_t> &t, 
                   std::string &name) {

    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_object_by_type("table",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type table found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    
    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    // Input the table_units data
    hdf_input_data(hf,t);

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    t.is_valid();
    
    return;
  }

  /** \brief Input a \ref o2scl::table_units object from a \ref
      hdf_file

      \note This function is declared in <tt>table_units.h</tt>, and 
      thus the default argument to the \c name parameter 
      appears there instead of here.

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table_units, otherwise, it reads
      the object of the specified name. If \c name is specified and no
      object is found, or if the object with the specified name is not
      of type \ref o2scl::table_units, then the error handler is called. If
      \c name is unspecified and there is no object with type \ref
      o2scl::table_units, then the error handler is called.
  */
  template<class vec_t> 
  void hdf_input(hdf_file &hf, o2scl::table_units<vec_t> &t, 
                 std::string name) {
    hdf_input_n<vec_t>(hf,t,name);
    return;
  }
  
  /** \brief Internal function for outputting a \ref o2scl::table_units object
   */
  void hdf_output_data(hdf_file &hf, o2scl::table_units<> &t);

  /** \brief Internal function for inputting a \ref o2scl::table_units object
   */
  template<class vec_t> 
  void hdf_input_data(hdf_file &hf, o2scl::table_units<vec_t> &t) {
    
    // Input base table object
    o2scl::table<vec_t,double> *tbase=
      dynamic_cast<o2scl::table_units<vec_t> *>(&t);
    if (tbase==0) {
      O2SCL_ERR2("Cast failed in hdf_input_data",
		 "(hdf_file &, table_units &).",o2scl::exc_efailed);
    }
    hdf_input_data(hf,*tbase);
  
    // Get unit flag
    int uf;
    hf.geti("unit_flag",uf);

    // If present, get units
    if (uf>0) {
      std::vector<std::string> units;
      hf.gets_vec_copy("units",units);
      for(size_t i=0;i<units.size();i++) {
	t.set_unit(t.get_column_name(i),units[i]);
      }
    }

    return;
  }
  
  /** \brief Output a \ref o2scl::hist object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::hist &h, std::string name);

  /** \brief Input a \ref o2scl::hist object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::hist,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::hist, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::hist, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::hist &h, std::string &name);

  /** \brief Input a \ref o2scl::hist object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::hist,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::hist, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::hist, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::hist &h, std::string name="");

  /** \brief Output a \ref o2scl::hist_2d object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, const o2scl::hist_2d &h, std::string name);

  /** \brief Input a \ref o2scl::hist_2d object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::hist_2d,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::hist_2d, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::hist_2d, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::hist_2d &h, std::string &name);

  /** \brief Input a \ref o2scl::hist_2d object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::hist_2d,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::hist_2d, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::hist_2d, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::hist_2d &h, std::string name="");

  /** \brief Output a \ref o2scl::table3d object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, const o2scl::table3d &h, std::string name);

  /** \brief Input a \ref o2scl::table3d object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table3d,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::table3d, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::table3d, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::table3d &h, std::string &name);

  /** \brief Input a \ref o2scl::table3d object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::table3d,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::table3d, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::table3d, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::table3d &h, std::string name="");

  /** \brief Output a \ref o2scl::expval_scalar object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::expval_scalar &h,
		  std::string name);

  /** \brief Input a \ref o2scl::expval_scalar object from a \ref
      hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_scalar,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_scalar, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_scalar, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::expval_scalar &h,
                   std::string &name);

  /** \brief Input a \ref o2scl::expval_scalar object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_scalar,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_scalar, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_scalar, then the
      error handler is called.
   */
  void hdf_input(hdf_file &hf, o2scl::expval_scalar &h,
		 std::string name="");

  /** \brief Output a \ref o2scl::expval_vector object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::expval_vector &h,
		  std::string name);

  /** \brief Input a \ref o2scl::expval_vector object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_vector,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_vector, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_vector, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::expval_vector &h, std::string &name);

  /** \brief Input a \ref o2scl::expval_vector object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_vector,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_vector, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_vector, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::expval_vector &h, std::string name="");

  /** \brief Output a \ref o2scl::expval_matrix object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::expval_matrix &h,
		  std::string name);

  /** \brief Input a \ref o2scl::expval_matrix object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_matrix,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_matrix, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_matrix, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::expval_matrix &h, std::string &name);

  /** \brief Input a \ref o2scl::expval_matrix object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::expval_matrix,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::expval_matrix, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::expval_matrix, then the
      error handler is called.
   */
  void hdf_input(hdf_file &hf, o2scl::expval_matrix &h, std::string name="");

  /** \brief Output a \ref o2scl::uniform_grid object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		  std::string name);

  /** \brief Input a \ref o2scl::uniform_grid object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::uniform_grid,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::uniform_grid, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::uniform_grid, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::uniform_grid<double> &h, 
                   std::string &name);
  
  /** \brief Input a \ref o2scl::uniform_grid object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::uniform_grid,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::uniform_grid, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::uniform_grid, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		 std::string name="");
  
  /** \brief Output a vector of \ref o2scl::contour_line objects to a
      \ref hdf_file
   */
  void hdf_output(hdf_file &hf, const std::vector<o2scl::contour_line> &cl, 
		  std::string name);
  
  /** \brief Input a vector of \ref o2scl::contour_line objects from a
      \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::contour_line,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::contour_line, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::contour_line, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, std::vector<o2scl::contour_line> &cl, 
                   std::string &name);
  
  /** \brief Input a vector of \ref o2scl::contour_line objects from a
      \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::contour_line,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::contour_line, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::contour_line, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, std::vector<o2scl::contour_line> &cl, 
		 std::string name="");
  
  /** \brief Output a vector of \ref o2scl::edge_crossings objects to
      a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, const std::vector<o2scl::edge_crossings> &ec, 
		  std::string name);
  
  /** \brief Input a vector of \ref o2scl::edge_crossings objects from 
      a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::edge_crossings,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::edge_crossings, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::edge_crossings, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
  */
  void hdf_input_n(hdf_file &hf, std::vector<o2scl::edge_crossings> &ec, 
                   std::string &name);
  
  /** \brief Input a vector of \ref o2scl::edge_crossings objects from
      a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::edge_crossings,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::edge_crossings, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::edge_crossings, then the
      error handler is called.
   */
  void hdf_input(hdf_file &hf, std::vector<o2scl::edge_crossings> &ec, 
		 std::string name="");
  
  /** \brief Output a \ref o2scl::tensor_grid object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::tensor_grid<std::vector<double>,
		  std::vector<size_t> > &t, std::string name);
  
  /** \brief Input a \ref o2scl::tensor_grid object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::tensor_grid,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::tensor_grid, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::tensor_grid, then the
      error handler is called. Upon exit, \c name contains the name
      of the object which was read.
   */
  void hdf_input_n(hdf_file &hf, o2scl::tensor_grid<std::vector<double>,
                   std::vector<size_t> > &t, std::string &name);
  
  /** \brief Input a \ref o2scl::tensor_grid object from a \ref hdf_file

      If \c name has a non-zero length, then this function first reads
      the first object of type \ref o2scl::tensor_grid,
      otherwise, it reads the object of the specified name. If \c name
      is specified and no object is found, or if the object with the
      specified name is not of type \ref o2scl::tensor_grid, then
      the error handler is called. If \c name is unspecified and there
      is no object with type \ref o2scl::tensor_grid, then the
      error handler is called. 
   */
  void hdf_input(hdf_file &hf, o2scl::tensor_grid<std::vector<double>,
		 std::vector<size_t> > &t, std::string name="");

  /** \brief Write a \ref o2scl::table_units object to an HDF5 file 
      with a given filename
   */
  void hdf5_write_file(o2scl::table_units<> &t, std::string fn,
                        std::string name="table");
  
  /** \brief A value specified by a string
      
      The first part of the specification is a "type" followed by a
      colon, followed by arguments which depend on the type. If no
      colon is present, then a "func:" prefix is assumed. The
      different types for a value specification are:
      
      1. <numeric value or function> - Value equal to the result of
         <function>, e.g. "7.6" or "sin(0.5)". See the \c functions
         help topic for a list of functions that can be used.
      
         For example:
      
         <tt>acol -create double "sqrt(5)" -output</tt>
      
      2. hdf5:<file>:<object name>:[addl. spec.] - Read an HDF5 value
         and obtain the value from object named <object name>. For
         some object types, additional specifications are required to
         specify which value should be used. A list of object types
         and additional specifications and more detail is given below.
      
         - <tt>double</tt>: (no addl. spec.)
         - <tt>int</tt>: (no addl. spec.)
         - <tt>size_t</tt>: (no addl. spec.)
         - <tt>double[]</tt>: index
         - <tt>int[]</tt>: index
         - <tt>size_t[]</tt>: index
         - <tt>uniform_grid<double></tt>: index
         - <tt>table</tt>: column name,row index
      
         For example:
      
         <tt>acol -create double hdf5:data/o2scl/apr98.o2:apr:rho,0 
         -output</tt>
      
      3. shell:<shell command> - Set the value equal to the first
         result obtained using the specified shell command. For
         example (using bash):
      
         <tt>acol -create double shell:"ls | wc | awk '{print $1}'" 
         -output</tt>
      
      4. python:<python code> - Set the value equal to the result
         obtained using the specified python code. For example (using
         bash):
      
         <tt>acol -create double 
         $'python:\"import numpy\nprint(numpy.sin(4))\"' -output</tt>
  */
  int value_spec(std::string spec, double &d,
		 int verbose=0, bool err_on_fail=true);
  
  /** \brief A vector specified by a string
      
      Some acol commands take arguments which are "vector
      specifications", i.e. an array specified as a string. The
      different parts of the string are separated by a colon, and the
      first part specifes the type of vector specification. The
      different types are:

      1. val:<number or function> - Create a vector with one element,
      which may be a number or a simple function, e.g.
      <tt>val:sin(0.5)</tt>.

      2. list:<entry 0>,<entry 1>, ..., <entry n-1> - Create a vector
      with a simple list of numbers or functions, e.g.
      <tt>list:3.0,1.0e-3,sqrt(2.0)</tt>.

      3. func:<N>:<function of i> - Create a vector by specifying the
      length of the vector and a function used to fill the elements.
      For example: <tt>func:41:sin(i/20.0*acos(-1))</tt>.

      4. grid:<begin>,<end>,<width>,["log"] - Create a vector equal to
      a uniform grid, e.g. use 'grid:1.0,10.0,1.0' for a 10-element
      vector filled with the numbers 1 to 10. The grid arguments can
      be values or mathematical expressions. Note that, in O2scl, grids
      must have at least 2 elements.

      5. text:<filename>:<column index> - Read a text file and extract
      a vector of numbers from a column of the text file (starting
      with zero for the first column), ignoring any header rows which
      contain non-numeric values. For example
      <tt>text:~/temp.dat:2</tt> will construct a vector from the
      third column of the file "temp.dat" in the user's home
      directory.

      6. hdf5:<file name>:<object name>:[addtional spec.] - Read an
      HDF5 file and obtain a vector from the object with the specified
      name. The remaining parts of the string contain additional
      information which may be needed depending on the type of object
      stored in the HDF5 file. A list of object types and additional
      specifications and more detail is given below.

      - double: (no addl. spec.) Implies vector of size 1
      - double[]: (no addl. spec.)
      - hist: (no addl. spec.) Vector of histogram weights
      - int: (no addl. spec.) Implies vector of size 1
      - int[]: (no addl. spec.)
      - size_t: (no addl. spec.) Implies vector of size 1
      - size_t[]: (no addl. spec.)
      - table: <column> Selected column from table
      - table: <row>:<col pat> Selected row and columns  
      - uniform_grid<double>: (no addl. spec.)

      For table <row>:<col pat>, the first additional specification is
      a row number, which can be negative to refer to counting from the end
      of the table. The second additional specification is a pattern of
      column names using either '*' or '?'.

      End of runtime documentation.

      \note Any data in the vector \c before the function is 
      called will be lost.

      \warning Experimental.
  */
  template<class vec_t> int vector_spec
  (std::string spec, vec_t &v, bool use_regex=false,
   int verbose=0, bool err_on_fail=true) {

    if (verbose>2) {
      std::cout << "vector_spec(): parsing: " << spec << std::endl;
      std::cout << "vector_spec(): verbose is " << verbose << std::endl;
    }
      
    if (spec.find("val:")==0) {

      std::string temp=spec.substr(4,spec.length()-4);
      if (verbose>1) {
	std::cout << "vector_spec(): single value " << temp
		  << std::endl;
      }
      v.resize(1);
      v[0]=o2scl::function_to_double(temp);
	
    } else if (spec.find("list:")==0) {
	
      // List
      std::string list=spec.substr(5,spec.length()-5);
      std::vector<std::string> sv;
      o2scl::split_string_delim(list,sv,',');
      size_t n=sv.size();
      if (n==0) {
	if (err_on_fail) {
	  O2SCL_ERR2("String split failed, spec empty? ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 10;
	}
      }
      if (verbose>1) {
	std::cout << "vector_spec(): List " << list << std::endl;
	std::cout << "  " << n << " " << sv[0] << " ... "
                  << sv[n-1] << std::endl;
      }
      v.resize(n);
      for(size_t i=0;i<n;i++) {
	v[i]=o2scl::function_to_double(sv[i]);
      }
	
    } else if (spec.find("func:")==0) {
	
      // Function
      if (verbose>1) {
	std::cout << "vector_spec(): Function " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function specified but no array length specified ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
      }
      size_t n;
      int cret=o2scl::stoszt_nothrow(temp.substr(0,ncolon),n);
      if (cret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Conversion to size_t failed ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 2;
	}
      }
      if (verbose>1) {
	std::cout << "Size: " << n << std::endl;
      }
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent function specified ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 3;
	}
      }
      std::string func=temp.substr(ncolon+1,temp.length()-ncolon-1);
      if (verbose>1) {
	std::cout << "vector_spec(): Function: " << func << std::endl;
      }
      o2scl::calc_utf8<> calc;
      std::map<std::string,double> vars;
      calc.compile(func.c_str(),&vars);
      v.resize(n);
      for(size_t i=0;i<n;i++) {
	vars["i"]=((double)i);
	v[i]=calc.eval(&vars);
      }
	
    } else if (spec.find("grid:")==0) {
	
      // Grid
      if (verbose>1) {
	std::cout << "vector_spec(): Grid " << spec << std::endl;
      }

      std::string temp=spec.substr(5,spec.length()-5);
      
      std::vector<std::string> sv;
      o2scl::split_string_delim(temp,sv,',');
      if (sv.size()<3) {
        if (verbose>0) {
          std::cerr << "vector_spec(): Grid spec only given " << sv.size()
                    << " arguments: " << temp << std::endl;
        }
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough information for grid ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 7;
	}
      }
      if (verbose>1) {
	std::cout << "vector_spec(): Begin,end,width "
		  << sv[0] << " " << sv[1] << " " << sv[2] << std::endl;
      }
      if (sv.size()>=4 && sv[3]=="log") {
	o2scl::uniform_grid_log_end_width<double>
	  ug(o2scl::function_to_double(sv[0]),
	     o2scl::function_to_double(sv[1]),
	     o2scl::function_to_double(sv[2]));
	ug.vector(v);
      } else {
	o2scl::uniform_grid_end_width<double>
	  ug(o2scl::function_to_double(sv[0]),
	     o2scl::function_to_double(sv[1]),
	     o2scl::function_to_double(sv[2]));
	ug.vector(v);
      }
	
    } else if (spec.find("text:")==0) {
	
      // Text
      if (verbose>1) {
	std::cout << "vector_spec(): Text " << spec << std::endl;
      }

      std::vector<std::string> sv;
      o2scl::split_string_delim(spec,sv,':');
      
      if (sv.size()<3) {
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough information for text file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 12;
	}
      }

      // The column index is 'col'
      int col=o2scl::stoi(sv[2]);
      if (verbose>1) {
	std::cout << "vector_spec(): Filename,column "
                  << sv[1] << " " << col << std::endl;
      }
      if (col<0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Column is negative for text file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 13;
	}
      }

      // Open the file
      std::ifstream fin;
      fin.open(sv[1].c_str());

      // Read the header, if there is any
      std::string line, word;
      bool in_header=true;
      do {
	getline(fin,line);
	std::istringstream is(line);
	int i=0;
	for(;i<=col && (is >> word) && in_header;i++) {
	  if (i==col && o2scl::is_number(word)) {
	    in_header=false;
	  }
	  if (i==col && verbose>2) {
	    if (in_header) {
	      std::cout << "vector_spec(): Word: " << word
                        << " header." << std::endl;
	    } else {
	      std::cout << "vector_spec(): Word: " << word
                        << " start of data." << std::endl;
	    }
	  }
	}
      } while (in_header && !fin.eof());

      // If we got to the end of the file and there wasn't
      // any data then call the error handler.
      if (in_header) {
	if (err_on_fail) {
	  O2SCL_ERR2("Couldn't find a number in text file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 13;
	}
      }

      // Now store the data in a temporary vector
      std::vector<double> tempv;
      bool end_of_data=false;
      do {
	tempv.push_back(o2scl::stod(word));
	
	getline(fin,line);
	std::istringstream is(line);
	for(int i=0;i<=col;i++) {
	  if (!(is >> word)) {
	    i=col;
	    end_of_data=true;
	  }
	}
	if (!o2scl::is_number(word)) end_of_data=true;
	
      } while (!end_of_data && !fin.eof());

      // Close the file
      fin.close();

      // Copy the data back to the user-specified vector
      v.resize(tempv.size());
      for(size_t i=0;i<tempv.size();i++) {
	v[i]=tempv[i];
      }
	
    } else if (spec.find("hdf5:")==0) {
	
      // HDF5 object in a file
      if (verbose>1) {
	std::cout << "vector_spec(): HDF5 file " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 4;
	}
      }
      std::string fname=temp.substr(0,ncolon);
      if (verbose>1) {
	std::cout << "vector_spec(): Filename " << fname << std::endl;
      }
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 5;
	}
      }
      
      std::string obj_name=temp.substr(ncolon+1,temp.length()-ncolon-1);
      std::string addl_spec;
      ncolon=obj_name.find(':');
      if (ncolon!=std::string::npos) {
	addl_spec=obj_name.substr(ncolon+1,obj_name.length()-ncolon-1);
	obj_name=obj_name.substr(0,ncolon);
      } 
      if (verbose>1) {
        if (obj_name[0]=='_') {
          std::cout << "vector_spec(): Object type: "
                    << obj_name.substr(1,obj_name.length()) << std::endl;
        }
	std::cout << "vector_spec(): Object name: " << obj_name << std::endl;
	std::cout << "vector_spec(): Additional specification: "
                  << addl_spec << std::endl;
      }
      o2scl_hdf::hdf_file hf;
	
      std::string fname_old=fname;
      std::vector<std::string> matches;
      int wret=o2scl::wordexp_wrapper(fname_old,matches);
      if (matches.size()>1 || matches.size()==0 || wret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function wordexp_wrapper() failed ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 9;
	}
      }
      fname=matches[0];
      if (verbose>1) {
	std::cout << "vector_spec(): Filename after wordexp_wrapper(): "
                  << fname << std::endl;
      }
	
      hf.open(fname);
      std::string type;
      int find_ret;
      // If the object name begins with an underscore,
      // interpret it as a type and read the first object
      // of that time
      if (obj_name[0]=='_') {
        type=obj_name.substr(1,obj_name.length()-1);
        find_ret=hf.find_object_by_type(type,obj_name);
      } else {
        find_ret=hf.find_object_by_name(obj_name,type);
      }
      if (find_ret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Object not found in file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 11;
	}
      }
      if (verbose>1) {
        std::cout << "vector_spec(): Object type and name: "
                  << type << " " << obj_name << std::endl;
      }
	
      if (type=="table") {
	
	if (addl_spec.length()==0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("No table column name specified ",
		       "in vector_spec().",o2scl::exc_einval);
	  } else {
	    return 6;
	  }
	}

	if (addl_spec.find(':')!=std::string::npos) {

	  if (verbose>2) {
	    std::cout << "vector_spec(): Found colon, so assuming "
                      << "row:column patterns." << std::endl;
	  }
	  std::vector<std::string> sv2;
	  o2scl::split_string_delim(addl_spec,sv2,':');
	  
          o2scl::calc_utf8<> calc;
	  calc.compile(sv2[0].c_str(),0);
	  int row=(int)calc.eval(0);
	  if (verbose>2) {
	    std::cout << "vector_spec(): Row is: " << row << std::endl;
	  }
	  
	  o2scl::table_units<> t;
	  o2scl_hdf::hdf_input(hf,t,obj_name);
	  
	  if (row<0) row+=t.get_nlines();
	  if (verbose>2) {
	    std::cout << "vector_spec(): Row+nlines is: " << row << std::endl;
	  }
	  
	  if (row<0 || row>=((int)t.get_nlines())) {
	    if (err_on_fail) {
	      O2SCL_ERR2("Row specification incorrect in ",
			 "vector_spec().",o2scl::exc_einval);
	    } else {
	      return 6;
	    }
	  }
	  
	  std::vector<std::string> col_list;
	  std::vector<std::string> col_patterns;
	  o2scl::split_string_delim(sv2[1],col_patterns,',');
	  if (verbose>2) {
	    std::cout << "vector_spec(): Column patterns: ";
	    o2scl::vector_out(std::cout,col_patterns,true);
	  }
	  
	  for(size_t k=0;k<col_patterns.size();k++) {
	    if (verbose>2) {
	      std::cout << "vector_spec(): Column pattern "
                        << col_patterns[k] << std::endl;
	    }

            if (use_regex) {
              std::regex r(col_patterns[k]);
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (std::regex_search(t.get_column_name(j),r)) {
                  col_list.push_back(t.get_column_name(j));
                  if (verbose>2) {
                    std::cout << "vector_spec(): Found match (using regex): "
                              << t.get_column_name(j)
                              << std::endl;
                  }
                }
              }
            } else {
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (fnmatch(col_patterns[k].c_str(),
                            t.get_column_name(j).c_str(),0)==0) {
                  col_list.push_back(t.get_column_name(j));
                  if (verbose>2) {
                    std::cout << "vector_spec(): Found match (using fnmatch): "
                              << t.get_column_name(j)
                              << std::endl;
                  }
                }
              }
            }
            
	  }
	  
	  v.resize(col_list.size());
	  for(size_t i=0;i<col_list.size();i++) {
	    v[i]=t.get(col_list[i],row);
	    if (verbose>2) {
	      std::cout << "vector_spec(): Getting entry at: "
                        << col_list[i] << " "
			<< row << " " << v[i] << std::endl;
	    }
	  }

	} else {
	  
	  if (verbose>2) {
	    std::cout << "vector_spec(): No colon, so assuming "
                      << "column name " << addl_spec << std::endl;
	  }
	  
	  o2scl::table_units<> t;
	  o2scl_hdf::hdf_input(hf,t,obj_name);
	  v.resize(t.get_nlines());
	  for(size_t i=0;i<t.get_nlines();i++) {
	    v[i]=t.get(addl_spec,i);
	  }
	}
	
      } else if (type=="double[]") {
	
	std::vector<double> vtemp;
	hf.getd_vec(obj_name,vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	
      } else if (type=="hist") {
	
	o2scl::hist ht;
	hdf_input(hf,ht,obj_name);
	typedef boost::numeric::ublas::vector<double> ubvector;
	const ubvector &wgts=ht.get_wgts();
	v.resize(wgts.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=wgts[i];
	}
	
      } else if (type=="int[]") {
	
	std::vector<int> vtemp;
	hf.geti_vec(obj_name,vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	
      } else if (type=="size_t[]") {
	
	std::vector<size_t> vtemp;
	hf.get_szt_vec(obj_name,vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	
      } else if (type=="uniform_grid<double>") {
	
	o2scl::uniform_grid<double> ug;
	hdf_input(hf,ug,obj_name);
	std::vector<double> vtemp;
	ug.vector(vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	
      } else if (type=="int") {
	
	int itemp;
	hf.geti(obj_name,itemp);
	v.resize(1);
	v[0]=itemp;
	
      } else if (type=="double") {
	
	double dtemp;
	hf.getd(obj_name,dtemp);
	v.resize(1);
	v[0]=dtemp;
	
      } else if (type=="size_t") {
	
	size_t szttemp;
	hf.get_szt(obj_name,szttemp);
	v.resize(1);
	v[0]=szttemp;
	
      } else {

	if (err_on_fail) {
	  O2SCL_ERR2("Cannot handle type ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
	
      }
      
      hf.close();
	
    } else {

      if (verbose>0) {
	std::cout << "Could not understand prefix in vector_spec()."
		  << std::endl;
      }
	
      if (err_on_fail) {
	O2SCL_ERR2("Could not parse specification ",
		   "in vector_spec().",o2scl::exc_einval);
      } else {
	return 8;
      }
	
    }

    return 0;
  }

  /** \brief Specification of several strings

      Some acol commands take arguments which are 'string list
      specifications'. The different parts of the string are separated
      by a colon, and the first part specifes the type of vector
      specification. The different types are:

      1. list:<comma-separated list> - A list of strings

      2. shell:<command> - The lines obtained from the result of a
      shell command, with a maximum of 256 characters per line.

      3. pattern:N:x[0][a][A] - The N strings obtained from a pattern.
      Occurrences of [0] are replaced with the integer 'i' where i
      runs from 0 to N-1. Occurrences of [a] are replaced with 'a'
      through 'z' from 0 through 25, and 'aa' through 'zz' for i from
      26 to 701. Occurrences of [A] are replaced with 'A' through 'Z'
      from 0 through 25, and 'AA' through 'ZZ' for i from 26 to 701.

      4. text:<filename> - The lines in the text file.

      5. hdf5: - Unfinished.

      This function is used for the acol slack command.

      End of runtime documentation.

      \note Previous strings in \c v are left unchanged and
      new strings are just added to the end.

      \warning Experimental.
  */
  template<class vec_t> int strings_spec(std::string spec, vec_t &v,
					 int verbose=0,
					 bool err_on_fail=true) {

    if (verbose>2) {
      std::cout << "Function strings_spec is parsing: " << spec << std::endl;
    }
      
    if (spec.find("list:")==0) {

      // List
      std::string list=spec.substr(5,spec.length()-5);
      std::vector<std::string> sv;
      o2scl::split_string_delim(list,sv,',');

      size_t n=sv.size();
      if (n==0) {
	if (err_on_fail) {
	  O2SCL_ERR2("String split failed, spec empty? ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 10;
	}
      }
      if (verbose>1) {
	std::cout << "strings_spec(): List " << list << std::endl;
	std::cout << "  N,first,last: "
                  << n << " " << sv[0] << " " << sv[n-1] << std::endl;
      }
      for(size_t i=0;i<n;i++) {
	v.push_back(sv[i]);
      }
	
    } else if (spec.find("pattern:")==0) {

      std::vector<std::string> sv;
      o2scl::split_string_delim(spec,sv,':');
      if (sv.size()<3) {
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough arguments for pattern specification ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
      }

      size_t np=o2scl::stoszt(sv[1]);
      std::string pat=sv[2];
      if (verbose>2) {
	std::cout << "pattern: np,pat: " << np << " " << pat
		  << std::endl;
      }
      
      for(size_t i=0;i<np;i++) {
	std::string cs=pat;
	while (cs.find("[0]")!=std::string::npos) {
	  cs=cs.replace(cs.find("[0]"),3,o2scl::szttos(i));
	}
	if (i<26) {
	  int ai='a';
	  int Ai='A';
	  char aic=ai+i;
	  char Aic=Ai+i;
	  std::string ais;
	  std::string Ais;
	  ais+=aic;
	  Ais+=Aic;
	  while (cs.find("[a]")!=std::string::npos) {
	    cs=cs.replace(cs.find("[a]"),3,ais);
	  }
	  while (cs.find("[A]")!=std::string::npos) {
	    cs=cs.replace(cs.find("[A]"),3,Ais);
	  }
	} else if (i<26*26+26) {
	  int ai='a';
	  int Ai='A';
	  size_t i2=i-26;
	  char aic=ai+i2/26;
	  char Aic=Ai+i2/26;
	  std::string ais;
	  std::string Ais;
	  ais+=aic;
	  Ais+=Aic;
	  char aic2=ai+i2%26;
	  char Aic2=Ai+i2%26;
	  ais+=aic2;
	  Ais+=Aic2;
	  while (cs.find("[a]")!=std::string::npos) {
	    cs=cs.replace(cs.find("[a]"),3,ais);
	  }
	  while (cs.find("[A]")!=std::string::npos) {
	    cs=cs.replace(cs.find("[A]"),3,Ais);
	  }
	} else {
	  if (cs.find("[a]")!=std::string::npos ||
	      cs.find("[A]")!=std::string::npos) {
	    if (err_on_fail) {
	      O2SCL_ERR2("Ran out of alphabet in ",
			 "strings_spec().",o2scl::exc_einval);
	    }
	    std::cerr << "Ran out of alphabet." << std::endl;
	    return 1;
	  }
	}
	if (verbose>2) {
	  std::cout << "String (" << i << "/" << np << ") from pattern "
		    << pat << " is " << cs << std::endl;
	}
	v.push_back(cs);
      }
      
    } else if (spec.find("shell:")==0) {
      
      // Result from shell command
      
#ifdef HAVE_POPEN
      
      std::string cmd=spec.substr(6,spec.length()-6);
      std::cout << "Using shell command: " << cmd << std::endl;
      FILE *ps_pipe=popen(cmd.c_str(),"r");
      if (!ps_pipe) {
	if (err_on_fail) {
	  O2SCL_ERR2("Pipe could not be opened in ",
		     "strings_spec().",
		     o2scl::exc_efilenotfound);
	}
	return o2scl::exc_efilenotfound;
      }
      
      static const size_t cl_max=256;
      char line1[cl_max];
      char *cret=fgets(line1,cl_max,ps_pipe);
      while (cret!=0) {
	std::string sline1=line1;
	if (sline1[sline1.length()-1]=='\n') {
	  sline1=sline1.substr(0,sline1.length()-1);
	}
	if (verbose>0) {
	  std::cout << "Read line: "
		    << sline1 << std::endl;
	}
	v.push_back(sline1);
	cret=fgets(line1,cl_max,ps_pipe);
      }
      
      if (pclose(ps_pipe)!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Pipe could not be closed in ",
		     "strings_spec().",o2scl::exc_efailed);
	}
	return o2scl::exc_efailed;
      }

#endif
      
      return 0;

    } else if (spec.find("text:")==0) {
      
      std::string file=spec.substr(5,spec.length()-5);
      std::vector<std::string> matches;
      int wret=o2scl::wordexp_wrapper(file,matches);
      if (matches.size()>1 || matches.size()==0 || wret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function wordexp_wrapper() failed ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 9;
	}
      }
      file=matches[0];
      
      std::cout << "Reading from text file: " << file << std::endl;
      std::ifstream fin(file.c_str());
      std::string stmp;
      while (getline(fin,stmp)) {
	v.push_back(stmp);
      }
      fin.close();
      
      return 0;

    } else if (spec.find("python:")==0) {
      
      // Result from shell command
      
#ifdef HAVE_POPEN
      
      std::string cmd=((std::string)"python -c ")+
	spec.substr(7,spec.length()-7);
      std::cout << "Using shell command: " << cmd << std::endl;
      FILE *ps_pipe=popen(cmd.c_str(),"r");
      if (!ps_pipe) {
	if (err_on_fail) {
	  O2SCL_ERR2("Pipe could not be opened in ",
		     "convert_units::convert_gnu_units().",
		     o2scl::exc_efilenotfound);
	}
	return o2scl::exc_efilenotfound;
      }
      
      static const size_t cl_max=256;
      char line1[cl_max];
      char *cret=fgets(line1,cl_max,ps_pipe);
      while (cret!=0) {
	std::string sline1=line1;
	if (sline1[sline1.length()-1]=='\n') {
	  sline1=sline1.substr(0,sline1.length()-1);
	}
	if (verbose>0) {
	  std::cout << "Read line "
		    << sline1 << std::endl;
	}
	v.push_back(sline1);
	cret=fgets(line1,cl_max,ps_pipe);
      }
      
      if (pclose(ps_pipe)!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Pipe could not be closed in ",
		     "strings_spec().",o2scl::exc_efailed);
	}
	return o2scl::exc_efailed;
      }

#endif
      
      return o2scl::exc_efailed;

    } else if (spec.find("hdf5:")==0) {
	
      // HDF5 object in a file
      if (verbose>1) {
	std::cout << "strings_spec(): HDF5 file " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 4;
	}
      }
      std::string fname=temp.substr(0,ncolon);
      if (verbose>1) {
	std::cout << "Filename " << fname << std::endl;
      }
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 5;
	}
      }
      std::string obj_name=temp.substr(ncolon+1,temp.length()-ncolon-1);
      std::string addl_spec;
      ncolon=obj_name.find(':');
      if (ncolon!=std::string::npos) {
	addl_spec=obj_name.substr(ncolon+1,obj_name.length()-ncolon-1);
	obj_name=obj_name.substr(0,ncolon);
      } 
      if (verbose>1) {
        if (obj_name[0]=='_') {
          std::cout << "Object type: "
                    << obj_name.substr(1,obj_name.length()) << std::endl;
        } else {
          std::cout << "Object name: " << obj_name << std::endl;
        }
	std::cout << "Additional specification: " << addl_spec << std::endl;
      }
      o2scl_hdf::hdf_file hf;
      
      std::string fname_old=fname;
      std::vector<std::string> matches;
      int wret=o2scl::wordexp_wrapper(fname_old,matches);
      if (matches.size()>1 || matches.size()==0 || wret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function wordexp_wrapper() failed ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 9;
	}
      }
      fname=matches[0];
      if (verbose>1) {
	std::cout << "Filename after wordexp() " << fname << std::endl;
      }
	
      hf.open(fname);
      std::string type;
      int find_ret;
      // If the object name begins with an underscore,
      // interpret it as a type and read the first object
      // of that time
      if (obj_name[0]=='_') {
        type=obj_name.substr(1,obj_name.length()-1);
        find_ret=hf.find_object_by_type(type,obj_name);
      } else {
        find_ret=hf.find_object_by_name(obj_name,type);
      }
      if (find_ret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Object not found in file ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 11;
	}
      }
      if (verbose>1) {
        std::cout << "Object type and name: " << type << " " << obj_name
                  << std::endl;
      }

      if (type=="string") {
	
	std::string stmp;
	hf.gets(obj_name,stmp);
	v.push_back(stmp);
	
      } else if (type=="string[]") {

	std::vector<std::string> vtmp;
	hf.gets_vec_copy(obj_name,vtmp);
	for(size_t k=0;k<vtmp.size();k++) {
	  v.push_back(vtmp[k]);
	}
	
      } else {

	if (err_on_fail) {
	  O2SCL_ERR2("Cannot handle type ",
		     "in strings_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
	
      }
      
    } else {

      v.push_back(spec);
      
    }
    
    return 0;
  }

  /** \brief Convert a vector specification to a 
      \c std::vector
  */
  std::vector<double> vector_spec(std::string spec);
  
  /** \brief A list of vectors specified by a string

      Some acol commands take arguments which are 'multiple vector
      specifications', i.e. a set of arrays specified as a string. The
      different parts of the string are separated by a colon, and the
      first part specifes the type of multiple vector specification.
      The different types are:

      1. func:<N>:<function of i>:<function of i and j> - Specify the
      number of vectors, a function of "i" which determines the length
      of the ith vector, and a function of "i" and "j" which specifies
      the jth element of the ith vector.

      2. text:<filename pattern>:<numeric column list> - Read one or
      more text files and extract vectors of numbers from columns of
      the text file, ignoring any header rows which contain
      non-numeric values. For example 'text:~/temp.dat:2-4' will
      construct vectors from the third, fourth, and fifth columns of
      the file 'temp.dat' in the user's home directory.

      3. hdf5:<filename pattern>:<object name>:[additional spec.] -
      Read one or more HDF5 files and obtain a vector from the object with
      the specified name. The remaining parts of the string contain
      additional information which may be needed depending on the type of
      object stored in the HDF5 file. 

      - type: addl. spec. 
      - table: <column pattern> table:<row list>:<column pattern>
      - table3d: <slice>:<index>

      Also, many normal vector specifications (from 'acol -help
      vector-spec') also work as multiple vector specifications. These
      include specifications which begin with 'val:', 'list:',
      'grid:', and 'table-row:'. Also included are 'hdf5:'
      specifications which refer to objects of type double, double[],
      hist, int, int[], size_t, size_t[], and uniform_grid<double>.

      End of runtime documentation.
      
      \note The vector \c v is not cleared and new vector specified
      in \c spec are added to the end of \c v.

      \warning Experimental.

      Used in acol create table-mv.

      \future When hdf5 is a single vector spec, it has to
      close the file so that vector_spec() can open it again.
      This should be fixed. Maybe the way to improve this is to
      break it up into several functions.
  */
  template<class vec_t> int mult_vector_spec
  (std::string spec, std::vector<vec_t> &v, bool use_regex=false,
   int verbose=0, bool err_on_fail=true) {

    if (spec.find("list:")==0 || spec.find("grid:")==0 ||
	spec.find("val:")==0) {
      
      // If the user specifies a list, grid, or value, then just use
      // the vector_spec() code
      vec_t v2;
      vector_spec(spec,v2,verbose,err_on_fail);
      v.push_back(v2);
      
    } else if (spec.find("func:")==0) {
	
      // Function
      if (verbose>1) {
	std::cout << "mult_vector_spec(): Function " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      std::vector<std::string> sv;
      o2scl::split_string_delim(temp,sv,':');

      if (sv.size()<2) {
	if (err_on_fail) {
	  O2SCL_ERR2("Less than three parts for function type ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
      }

      // There were only two arguments to the func specification,
      // so presume a normal vector spec
      if (sv.size()==2) {

	if (verbose>1) {
	  std::cout << "mult_vector_spec(): Only two arguments "
		    << "to \"func:\", presuming vector spec."
		    << std::endl;
	}
	
	vec_t v2;
	vector_spec(spec,v2,verbose,err_on_fail);
	v.push_back(v2);
	
      } else {

	// --------------------------------------------------------
	// Otherwise, there were three arguments to the function
	// specification, so proceed as normal

	// First determine the number of functions
	
	size_t n;
        if (verbose>1) {
          std::cout << "mult_vector_spec(): string for size: "
                    << sv[0] << std::endl;
        }
	int cret=o2scl::stoszt_nothrow(sv[0],n);
	if (cret!=0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Conversion to size_t failed ",
		       "in mult_vector_spec().",o2scl::exc_einval);
	  } else {
	    return 2;
	  }
	}
	if (verbose>1) {
	  std::cout << "Size " << n << std::endl;
	}

	// Now, loop through each function
	
	for(size_t i=0;i<n;i++) {

	  // Compile the function for the length of the ith vector
	  
          o2scl::calc_utf8<> calc2;
	  std::map<std::string,double> vars2;
	  vars2["i"]=((double)i);
	  
	  int cret2=calc2.compile_nothrow(sv[1].c_str(),&vars2);
	  if (cret2!=0) {
	    if (err_on_fail) {
	      O2SCL_ERR2("Function to get vector size failed ",
			 "in mult_vector_spec().",o2scl::exc_einval);
	    } else {
	      return 2;
	    }
	  }

	  // Evaulate the function for the length of the ith vector

	  double ce;
	  int cret3=calc2.eval_nothrow(&vars2,ce);
	  if (cret3!=0) {
	    if (err_on_fail) {
	      O2SCL_ERR2("Function to get vector size failed ",
			 "in mult_vector_spec().",o2scl::exc_einval);
	    } else {
	      return 2;
	    }
	  }
	  size_t n2=(size_t)ce;
	  
	  if (verbose>1) {
	    std::cout << "Size of vector " << n << " is " << n2 << std::endl;
	    std::cout << "Function " << sv[2] << std::endl;
	  }
	  
	  // Compile the function for the ith vector

          o2scl::calc_utf8<> calc;
	  std::map<std::string,double> vars;
	  vars["i"]=((double)i);
	  
	  int cret4=calc.compile_nothrow(sv[2].c_str(),&vars);
	  if (cret4!=0) {
	    if (err_on_fail) {
	      O2SCL_ERR2("Function to get vector size failed ",
			 "in mult_vector_spec().",o2scl::exc_einval);
	    } else {
	      return 2;
	    }
	  }
	  
	  // Evaluate the function for the ith vector

	  std::vector<double> vtemp;
	  vtemp.resize(n2);
	  for(size_t j=0;j<n2;j++) {
	    vars["j"]=((double)j);
	    vtemp[j]=calc.eval(&vars);
	  }

	  // Add the temporary vector to the list
	  v.push_back(vtemp);
	  
	}

      }
	
    } else if (spec.find("text:")==0) {

      // Text file specification
      
      if (verbose>1) {
	std::cout << "mult_vector_spec(): Text " << spec << std::endl;
      }

      std::vector<std::string> sv;
      o2scl::split_string_delim(spec,sv,':');
      
      if (sv.size()<3) {
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough information for text file ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 12;
	}
      }

      std::string col_list=sv[2];
      if (verbose>1) {
	std::cout << "Filename,column list: " << sv[1] << " "
		  << col_list << std::endl;
      }

      std::vector<size_t> col_list2;
      o2scl::string_to_uint_list(col_list,col_list2);
      o2scl::vector_sort<std::vector<size_t>,size_t>(col_list2.size(),
						     col_list2);
      if (verbose>1) {
	std::cout << "Column list: ";
	o2scl::vector_out(std::cout,col_list2,true);
      }
      
      if (col_list2.size()<1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No column list in ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 20;
	}
      }

      // Open the file
      std::ifstream fin;
      fin.open(sv[1].c_str());

      // Read the header, if there is any
      std::string line, word;
      bool in_header;
      do {
	in_header=false;
        getline(fin,line);
	if (verbose>2) {
	  std::cout << "line: " << line << std::endl;
	}
        std::istringstream is(line);
	while ((is >> word) && (in_header==false)) {
	  if (verbose>2) {
	    std::cout << "word: " << word << " is_number: "
		      << o2scl::is_number(word) << " eof: " << fin.eof()
		      << std::endl;
	  }
          if (!o2scl::is_number(word)) {
            in_header=true;
          }
        }
	if (verbose>2) {
	  std::cout << "in_header: " << in_header << " eof: "
		    << fin.eof() << std::endl;
	}
      } while (in_header && !fin.eof());
      
      // If we got to the end of the file and there wasn't
      // any data then call the error handler.
      if (in_header) {
        if (err_on_fail) {
          O2SCL_ERR2("Couldn't find a number in text file ",
                     "in mult_vector_spec().",o2scl::exc_einval);
        } else {
          return 13;
        }
      }

      // Compute the last column from the column list so we can skip
      // the latter entries in the line which we won't use
      size_t last_col=o2scl::vector_max_value<std::vector<size_t>,size_t>
	(col_list2);
      if (verbose>1) {
	std::cout << "Last column is: " << last_col << std::endl;
      }
      
      // Store the data in a temporary object
      std::vector<std::vector<double> > vtemp;
      vtemp.resize(col_list2.size());
      
      bool end_of_data=false;
      
      do {
	
        std::istringstream is(line);

	for(size_t i=0;i<=last_col && end_of_data==false;i++) {
          if (!(is >> word)) {
	    end_of_data=true;
	  } else {
	    size_t ix;
	    if (o2scl::vector_search(col_list2,i,ix)==true) {
	      vtemp[ix].push_back(o2scl::stod(word));
	    }
	  }
	}

	if (end_of_data==false) {
	  getline(fin,line);
	}
        
      } while (!end_of_data && !fin.eof());

      // Add the temporary data to the vector list
      for(size_t i=0;i<vtemp.size();i++) {
	v.push_back(vtemp[i]);
      }
      
      // Close the file
      fin.close();
      
    } else if (spec.find("hdf5:")==0) {

      // HDF5 object in a file
      if (verbose>1) {
	std::cout << "mult_vector_spec(), HDF5 file, spec: "
                  << spec << std::endl;
      }

      // Count the number of colons after "hdf5:"
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent file name specified ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 4;
	}
      }

      // Extract the filename
      std::string fname=temp.substr(0,ncolon);
      if (verbose>1) {
	std::cout << "Filename: " << fname << std::endl;
      }
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 5;
	}
      }

      // Extract the object name and any additional specifications
      std::string obj_name=temp.substr(ncolon+1,temp.length()-ncolon-1);
      std::string addl_spec;
      ncolon=obj_name.find(':');
      if (ncolon!=std::string::npos) {
	addl_spec=obj_name.substr(ncolon+1,obj_name.length()-ncolon-1);
	obj_name=obj_name.substr(0,ncolon);
      } 
      if (verbose>1) {
        if (obj_name[0]=='_') {
          std::cout << "Object type: "
                    << obj_name.substr(1,obj_name.length()-1) << std::endl;
        } else {
          std::cout << "Object name: " << obj_name << std::endl;
        }
	std::cout << "Additional specification: " << addl_spec << std::endl;
      }

      // Use wordexp_wrapper() to determine the filelist
      o2scl_hdf::hdf_file hf;
      std::string fname_old=fname;
      std::vector<std::string> matches;
      int wret=o2scl::wordexp_wrapper(fname_old,matches);

      if (matches.size()==0 || wret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function wordexp_wrapper() failed ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 9;
	}
      }

      size_t nfiles=matches.size();
      if (verbose>1) {
	std::cout << "Wordexp() found " << nfiles << " files."
		  << std::endl;
      }

      // Loop over each file
      for(size_t ifile=0;ifile<nfiles;ifile++) {
        
	fname=matches[ifile];
	if (verbose>1) {
	  std::cout << "Filename for index " << ifile << " is " << fname
		    << std::endl;
	}

        // Read object by type or name
	hf.open(fname);
	std::string type;
        int find_ret;
        // If the object name begins with an underscore,
        // interpret it as a type and read the first object
        // of that time
        if (obj_name[0]=='_') {
          type=obj_name.substr(1,obj_name.length()-1);
          find_ret=hf.find_object_by_type(type,obj_name);
        } else {
          find_ret=hf.find_object_by_name(obj_name,type);
        }
	if (find_ret!=0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Object not found in file ",
		       "in mult_vector_spec().",o2scl::exc_einval);
	  } else {
	    return 11;
	  }
	}
	if (verbose>1) {
          std::cout << "Object type and name: " << type << " " << obj_name
                    << std::endl;
	}

        // Handle each object type separately
	if (type=="table") {
	  
	  if (addl_spec.length()==0) {
	    if (err_on_fail) {
	      O2SCL_ERR2("No table column name specified ",
			 "in mult_vector_spec().",o2scl::exc_einval);
	    } else {
	      return 6;
	    }
	  }

          o2scl::table_units<> t;
          o2scl_hdf::hdf_input(hf,t,obj_name);

          // The row:column pattern case for a table object
	  if (addl_spec.find(':')!=std::string::npos) {
	    if (verbose>1) {
	      std::cout << "Found row:column pattern in mult_vector_spec()."
			<< std::endl;
	    }

            int pos=addl_spec.find(':');
            if (pos==0 || pos==((int)addl_spec.length())-1) {
              std::cout << "Colon at beginning or end of additional "
                        << "specification in mult_vector_spec()."
                        << std::endl;
              return 2;
            }
            std::string row_list=addl_spec.substr(0,pos);
            std::string col_spec=addl_spec.substr(pos+1,
                                                  addl_spec.length()-pos-1);
            if (verbose>1) {
              std::cout << "row_list: " << row_list << std::endl;
              std::cout << "col_spec: " << col_spec << std::endl;
            }

            std::vector<size_t> uint_list;
            if (row_list!="*") {
              o2scl::string_to_uint_list(row_list,uint_list);
              if (uint_list.size()<=0) {
                if (err_on_fail) {
                  O2SCL_ERR2("Could not parse row list ",
                             "in mult_vector_spec().",o2scl::exc_einval);
                } else {
                  return 16;
                }
              } else if (verbose>1) {
                std::cout << "rows: " << uint_list[0];
                if (uint_list.size()>=2) {
                  std::cout << "," << uint_list[1];
                }
                if (uint_list.size()>=3) {
                  std::cout << "," << uint_list[2];
                }
                if (uint_list.size()==4) {
                  std::cout << "," << uint_list[3] << std::endl;
                } else if (uint_list.size()>4) {
                  std::cout << ",...," << uint_list[uint_list.size()-1]
                            << std::endl;
                } else {
                  std::cout << std::endl;
                }
              }
            }
            
            std::vector<size_t> col_ix;

            if (use_regex) {
              std::regex r(col_spec);
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (std::regex_search(t.get_column_name(j),r)) {
                  col_ix.push_back(j);
                  if (verbose>1) {
                    std::cout << "Found match (using regex): "
                              << t.get_column_name(j)
                              << std::endl;
                  }
                }
              }
            } else {
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (fnmatch(col_spec.c_str(),
                            t.get_column_name(j).c_str(),0)==0) {
                  col_ix.push_back(j);
                  if (verbose>1) {
                    std::cout << "Found match (using fnmatch): "
                              << t.get_column_name(j)
                              << std::endl;
                  }
                }
              }
            }

            if (row_list=="*") {
              for(size_t j=0;j<t.get_nlines();j++) {
                vec_t vtemp(col_ix.size());
                for(size_t k=0;k<col_ix.size();k++) {
                  vtemp[k]=t.get(col_ix[k],j);
                }
                v.push_back(vtemp);
              }
            } else {
              for(size_t j=0;j<uint_list.size();j++) {
                vec_t vtemp(col_ix.size());
                for(size_t k=0;k<col_ix.size();k++) {
                  vtemp[k]=t.get(col_ix[k],uint_list[j]);
                }
                v.push_back(vtemp);
              }
            }

            // End of the row:column pattern case for a table object
            
	  } else {

            // The set of columns case for a table object
            
            if (use_regex) {
              std::regex r(addl_spec);
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (std::regex_search(t.get_column_name(j),r)) {
                  if (verbose>1) {
                    std::cout << "Column " << t.get_column_name(j)
                              << " matches regex " << addl_spec
                              << "." << std::endl;
                  }
                  vec_t vtemp(t.get_nlines());
                  for(size_t k=0;k<t.get_nlines();k++) {
                    vtemp[k]=t.get(j,k);
                  }
                  v.push_back(vtemp);
                }
              }
            } else {
              for(size_t j=0;j<t.get_ncolumns();j++) {
                if (fnmatch(addl_spec.c_str(),
                            t.get_column_name(j).c_str(),0)==0) {
                  if (verbose>1) {
                    std::cout << "Column " << t.get_column_name(j)
                              << " matches pattern " << addl_spec
                              << "." << std::endl;
                  }
                  vec_t vtemp(t.get_nlines());
                  for(size_t k=0;k<t.get_nlines();k++) {
                    vtemp[k]=t.get(j,k);
                  }
                  v.push_back(vtemp);
                }
              }
            }

            // End of the set of columns case for a table object
          }
            
          hf.close();

          // End of the table object case
          
        } else if (type=="table3d") {

	  if (addl_spec.length()==0) {
	    if (err_on_fail) {
	      O2SCL_ERR2("No table3d slice or index specified ",
			 "in mult_vector_spec().",o2scl::exc_einval);
	    } else {
	      return 6;
	    }
	  }

          o2scl::table3d t3d;
          o2scl_hdf::hdf_input(hf,t3d,obj_name);

	  if (addl_spec.find(':')==std::string::npos) {
            if (err_on_fail) {
              O2SCL_ERR2("No additional spec. ",
			 "in mult_vector_spec().",o2scl::exc_einval);
            } else {
              return 21;
            }
          }
          
          int pos=addl_spec.find(':');
          if (pos==0 || pos==((int)addl_spec.length())-1) {
            std::cout << "Colon at beginning or end of additional "
                      << "specification in mult_vector_spec()."
                      << std::endl;
            return 22;
          }
          std::string slice_name=addl_spec.substr(0,pos);
          int index=o2scl::stoi(addl_spec.substr(pos+1,
                                                 addl_spec.length()-pos-1));

          if (verbose>2) {
            std::cout << "Using slice name " << slice_name
                      << " and index " << index << std::endl;
          }
          if (index==0) {
            for(size_t i=0;i<t3d.get_nx();i++) {
              std::vector<double> row;
              for(size_t j=0;j<t3d.get_ny();j++) {
                row.push_back(t3d.get(i,j,slice_name));
              }
              v.push_back(row);
            }
          } else if (index==1) {
            for(size_t i=0;i<t3d.get_ny();i++) {
              std::vector<double> col;
              for(size_t j=0;j<t3d.get_nx();j++) {
                col.push_back(t3d.get(j,i,slice_name));
              }
              v.push_back(col);
            }
          }

          std::cout << "v.size(): " << v.size() << std::endl;
          
          hf.close();
          
	} else if (type=="double" || type=="double[]" || type=="hist" ||
		   type=="int" || type=="int[]" || type=="size_t" ||
		   type=="size_t[]" || type=="uniform_grid<double>") {

          // If the spec can be interpreted as a single vector spec,
	  // then just use the vector_spec() function. Close the
	  // file so vector_spec() can reopen it.
	  hf.close();
	  
	  std::vector<double> vtemp;
	  int iret=vector_spec(spec,vtemp,verbose,err_on_fail);
	  v.push_back(vtemp);

	  if (iret!=0) return iret;
	  
	} else {
	  
	  if (err_on_fail) {
	    O2SCL_ERR2("Cannot handle type ",
		       "in mult_vector_spec().",o2scl::exc_einval);
	  } else {
	    return 1;
	  }
	  
	}

        // End of the loop over files
      }

    } else {
	
      if (err_on_fail) {
	O2SCL_ERR2("Could not parse specification ",
		   "in mult_vector_spec().",o2scl::exc_einval);
      } else {
	return 8;
      }
	
    }

    return 0;
  }

  /** \brief Convert a vector specification to a 
      \c std::vector
  */
  std::vector<std::vector<double>> mult_vector_spec(std::string spec);

  /** \brief Functions for acol
      
      Functions can be created using the operators and functions
      listed below. Examples are <tt>x==5 && y<1</tt>,
      <tt>acos(-1)</tt>, and <tt>sin(x>5)</tt>. Comparison operators
      result in either 1.0 (true) or 0.0 (false).

      There are two additional limitations. First, in order to avoid
      confusion between arguments and functions, use parenthesis and
      quotes, i.e. <tt>"(-x*2)"</tt> instead of <tt>-x*2</tt>. Also,
      do not use a unary minus next to a binary operator, i.e. use
      <tt>a>(-1)</tt> instead of <tt>a>-1</tt>.

      Operators:

      () ^ * / % + - == != < > && || << >> >= <=

      Power functions:

      sqrt(x), cbrt(x), pow(x,y), and hypot(x,y)

      Exponential functions:

      exp(x), log(x), log10(x), log1p(x), and expm1(x)

      Trigonometric functions:

      asin(x), acos(x), atan(x), sinh(x), cosh(x), tanh(x), asinh(x),
      acosh(x), atanh(x), and atan2(y,x)

      Exponential functions:

      erf(x) [2/sqrt(pi) ∫_0^{x} exp(-t^2) dt],
      erfc(x) [2/sqrt(pi) ∫_x^∞ exp(-t^2) dt = 1-erf(x)],
      tgamma(x) [∫ t_0^∞ t^{x-1} exp(-t) dt = (x-1)! for ints],
      and lgamma(x) [ln |tgamma(x)|]

      Bessel functions:

      cyl_bessel_i(nu,x) [i.e. Iᵥ(x)], cyl_bessel_j(nu,x) [i.e. Jᵥ(x)], 
      cyl_bessel_k(nu,x) [i.e. Kᵥ(x)], and cyl_neumann(nu,x) [i.e. Yᵥ(x)]
      [all defined only for x>0],

      Additional integrals (experimental):

      fdint(nu,x) [Fermi-Dirac integ., ∫_0^∞ t^nu/(exp(t-x)+1)],
      beint(nu,x) [Bose-Einstein integ. ∫_0^∞ t^nu/(exp(t-x)-1)],
      and polylog(nu,x) [for x<-1 only].

      Other functions:

      abs(x), ceil(x), floor(x), if(t,x,y) [If t>0.5 then x,
      otherwise y.], isfinite(x), isinf(x), isnan(x),
      max(x,y), min(x,y), and sqrt1pm1(x) [√(1+x)-1].

      Special values:

      false = 0, true = 1, rand (experimental) = random number

      Use <tt>acol -help function</tt> to get more information on the
      type-specific command called \c function .
  */
  // sph_bessel(nu,x) [i.e. jᵥ(x)], and sph_neumann(nu,x) [i.e. yᵥ(x)].
  // For sph_bessel and sph_neumann, the first argument must be an
  // integer (fractional parts are truncated). 
  void functions();

  /** \brief Index specifications for acol

      The \c tensor rearrange commands use index specifications to
      specify how the tensor should be rearranged. Index
      specifications may be specified as separate arguments or
      multiple index specifications may be given in a single argument
      separated by spaces or commas. The indices begin with 0, the
      first index so that index 1 is the second index. The list of
      index specifications is:

      - index(ix): Retain index ix in the new tensor.

      - fixed(ix): Fix the value of index ix.

      - sum(ix): Sum over the value of index ix

      - trace(ix1,ix2): Trace (sum) over indices ix and ix2. If the
        number of entries in either index is smaller than the other,
        then the remaining entries are ignored in the sum.

      - reverse(ix): Retain index ix but reverse the order.

      - range(ix,start,end): Retain index ix but modify range. Ranges
        include both of their endpoints.

      - interp(ix,value) (for \c tensor_grid ): fix index ix by
        interpolating 'value' into the grid for index ix.

      - grid(ix,begin,end,n_bins,log) (for \c tensor_grid ):
        interpolate the specified index on a grid to create a new
        index. If the value of 'log' is 1, then the grid is
        logarithmic.

      - gridw(ix,begin,end,bin_width,log) (for \c tensor_grid ):
        interpolate the specified index on a grid with a fixed bin
        width to create a new index. If the value of 'log' is 1, then
        the grid is logarithmic and 'bin_width' is the
        multiplicative factor between bin edges.

      Note that the index specifications which result in a tensor
      index (all except 'fixed', 'sum', 'trace' and 'interp') must be
      given in the order they should appear in the tensor which
      results. Also, the \c rearrange commands require that the result
      of the rearrangement must have at least one index left.

      Examples:

      <tt>acol -create tensor_grid 3 9 10 11 -set-grid 0 "i^2"
      -set-grid 1 "2*i" -set-grid 2 "(2-i)" -function "sin(x0+x1*x2)"
      -entry 4 2 5 -rearrange "index(1) index(0)" "fixed(2,5)"
      -entry 2 4</tt>

      This example creates a rank three tensor of size 9 by 10 by 11
      and then flips the first two indices while fixing the third
      index to select the 6th value to create a rank 2 tensor of 
      size 10 by 9. Additionally, this example demonstrates that
      one can put multiple index specifications into the same
      argument (as in the first argument to rearrange) or place
      additional index specifications in additional arguments.
      The \c entry command shows an entry in the original rank 3 
      tensor and the \c output command shows that that entry 
      has now moved to the to the (2,4) location in the final
      rank 2 tensor.

      <tt>acol -create tensor_grid 3 8 10 6 -set-grid 0 "i*m" \\
      -set-grid 1 "m/5*i+x" -set-grid 2 "(2-i)" -function \\
      "sin(x0+x1*x2)" -get-grid -entry 5 1 3 -entry 5 2 3 -rearrange \\
      "reverse(2),range(0,7,3),interp(1,4.5)" -get-grid -entry 2 2</tt>

      This example creates a rank three tensor of size 8 by 10 by 6
      and then creates a rank 2 tensor from it. In this example, the
      <tt>set-grid</tt> commands use variables 'm' and 'x'. In the
      first call to <tt>set-grid</tt>, the variable m refers to the
      value 8, while in the second call, the variable m refers to the
      value 10. The value 'x' always refers to the previous grid value
      which is set by the <tt>create</tt> command to be equal to the
      index. The first index in the rank 2 tensor is a reversed
      version of the last index in the rank 3 tensor. The second index
      of the rank 2 tensor is constructed from a subset of the second
      index of the rank 3 tensor and also reversed (since 7 and 3 are
      given in the opposite order). The rank two tensor is constructed
      by linearly interpolating the value 4.5 into the second index of
      the rank three tensor, thus halfway between values when the
      index takes values 1 and 2. The final tensor thus has size 6x5.
      As an example, the (2,2) entry in the final tensor comes from
      the average of the (5,1,3) and (5,2,3) entries in the original
      rank 3 tensor.

      <tt>acol -create tensor_grid 3 6 8 4 -set-grid 0 "erf(i*m/20)"
      -set-grid 1 "m/(i+1)" -set-grid 2 "exp(i)" -function
      "sin(x0+x1*x2+i0+i2)" -get-grid -rearrange "grid(2,5,20,5)"
      "gridw(1,7,1.2,0.9,1)" "interp(0,(erf(6/20)+erf(12/20))/2)" 
      -get-grid</tt>

      This example creates a rank three tensor of size 6 by 8 by 4 and
      then creates a rank 2 tensor from it. In this example, the
      <tt>set-grid</tt> commands uses the error and exponential
      functions to set the initial grid. The first index of the final
      tensor comes from the last index of the initial tensor, and the
      values are obtained by interpolating the original exponential
      grid to a linear grid beginning at 5 and ending at 20 with 5
      bins (i.e. 6 grid points). The second index of the final tensor
      is obtained by interpolating the second index of the initial
      tensor, creating a logarithmic grid from 7 to 1.2. The final
      tensor is also created by interpolating the result of the
      mathematical expression "(erf(6/20)+erf(12/20))/2" into the
      first index of the initial tensor.
  */
  void index_spec();
  
}

#endif
