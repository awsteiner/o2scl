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
#ifndef O2SCL_HDF_IO_H
#define O2SCL_HDF_IO_H

/** \file hdf_io.h
    \brief File defining HDF I/O for selected \o2 objects
*/
#include <boost/numeric/ublas/vector.hpp>

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

/** \brief The \o2 namespace for I/O with HDF
 */
namespace o2scl_hdf {
  
  /** \brief Output a \ref o2scl::table object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table<> &t, std::string name);

  /** \brief Input a \ref o2scl::table object from a \ref hdf_file

      \todo Removed default value for \c name for compiling at nersc
  */
  template<class vec_t> 
    void hdf_input(hdf_file &hf, o2scl::table<vec_t> &t, std::string name) {
      
    // If no name specified, find name of first group of specified type
    std::cout << "In hit1." << std::endl;
    if (name.length()==0) {
      hf.find_group_by_type(hf,"table",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type table found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    std::cout << "In hit2." << std::endl;

    // Open main group
    hid_t top=hf.get_current_id();
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);

    // Input the table data
    std::cout << "In hit3." << std::endl;
    hdf_input_data(hf,t);
    std::cout << "In hit4." << std::endl;

    // Close group
    hf.close_group(group);
    std::cout << "In hit5." << std::endl;

    // Return location to previous value
    hf.set_current_id(top);

    return;
  }

  /** \brief Internal function for outputting a \ref o2scl::table object
   */
  void hdf_output_data(hdf_file &hf, o2scl::table<> &t);

  /** \brief Internal function for inputting a \ref o2scl::table object
   */
  template<class vec_t> 
    void hdf_input_data(hdf_file &hf, o2scl::table<vec_t> &t) {
    hid_t group=hf.get_current_id();

    // Clear previous data
    std::cout << "In hitd1." << std::endl;
    t.clear_table();
    t.clear_constants();

    // Check typename
    std::cout << "In hitd2." << std::endl;
    std::string type2;
    hf.gets_fixed("o2scl_type",type2);
    if (type2!="table") {
      O2SCL_ERR2("Typename in HDF group does not match ",
		 "class in o2scl_hdf::hdf_input().",o2scl::exc_einval);
    }
    std::cout << "In hitd3." << std::endl;

    // Storage
    std::vector<std::string> cnames, cols;
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector cvalues;
      
    // Get constants
    hf.gets_vec("con_names",cnames);
    hf.getd_vec_copy("con_values",cvalues);
    if (cnames.size()!=cvalues.size()) {
      O2SCL_ERR2("Size mismatch between constant names and values ",
		 "in o2scl_hdf::hdf_input().",o2scl::exc_einval);
    }
    for(size_t i=0;i<cnames.size();i++) {
      t.add_constant(cnames[i],cvalues[i]);
    }

    // Get column names
    hf.gets_vec("col_names",cols);
    for(size_t i=0;i<cols.size();i++) {
      t.new_column(cols[i]);
    }

    std::cout << "In hitd4." << std::endl;
    // Get number of lines
    int nlines2;
    hf.geti("nlines",nlines2);
    t.set_nlines(nlines2);

    // Output the interpolation type
    hf.get_szt_def("itype",o2scl::itp_cspline,t.itype);

    // Open data group
    hid_t group2=hf.open_group("data");
    std::cout << "In hitd4b: " << group2 << std::endl;
    hf.set_current_id(group2);
    std::cout << "In hitd4c: " << group2 << std::endl;

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

    std::cout << "In hitd5: " << group2 << std::endl;
    // Close groups
    hf.close_group(group2);
    std::cout << "In hitd6." << std::endl;

    hf.set_current_id(group);

    // Check that input created a valid table
    std::cout << "In hitd7." << std::endl;
    t.check_synchro();

    return;
  }

  /** \brief Output a \ref o2scl::table_units object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table_units<> &t, 
		  std::string name);

  /** \brief Input a \ref o2scl::table_units object from a \ref hdf_file

      \todo Removed default value for \c name for compiling at nersc
  */
  template<class vec_t> 
    void hdf_input(hdf_file &hf, o2scl::table_units<vec_t> &t, 
		   std::string name) {
      
    std::cout << "In hitu1." << std::endl;
    // If no name specified, find name of first group of specified type
    if (name.length()==0) {
      hf.find_group_by_type(hf,"table",name);
      if (name.length()==0) {
	O2SCL_ERR2("No object of type table found in ",
		   "o2scl_hdf::hdf_input().",o2scl::exc_efailed);
      }
    }
    std::cout << "In hitu2." << std::endl;

    // Open main group
    hid_t top=hf.get_current_id();
    std::cout << "In hitu3." << std::endl;
    hid_t group=hf.open_group(name);
    std::cout << "In hitu4." << std::endl;
    hf.set_current_id(group);

    // Input the table_units data
    hdf_input_data(hf,t);

    // Close group
    std::cout << "In hitu5." << std::endl;
    hf.close_group(group);
    std::cout << "In hitu6." << std::endl;

    // Return location to previous value
    hf.set_current_id(top);

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
    o2scl::table<vec_t> *tbase=dynamic_cast<o2scl::table_units<vec_t> *>(&t);
    if (tbase==0) {
      O2SCL_ERR2("Cast failed in hdf_output_data",
		 "(hdf_file &, table_units &).",o2scl::exc_efailed);
    }
    std::cout << "In hidtu1." << std::endl;
    hdf_input_data(hf,*tbase);
    std::cout << "In hidtu2." << std::endl;
  
    // Get unit flag
    int uf;
    hf.geti("unit_flag",uf);

    // If present, get units
    std::cout << "In hidtu3." << std::endl;
    if (uf>0) {
      std::vector<std::string> units;
      hf.gets_vec("units",units);
      for(size_t i=0;i<units.size();i++) {
	t.set_unit(t.get_column_name(i),units[i]);
      }
    }
    std::cout << "In hidtu4." << std::endl;

    return;
  }

  /// Output a \ref o2scl::hist object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::hist &h, std::string name);
  /// Input a \ref o2scl::hist object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::hist &h, std::string name="");
  /// Output a \ref o2scl::hist_2d object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::hist_2d &h, std::string name);
  /// Input a \ref o2scl::hist_2d object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::hist_2d &h, std::string name="");
  /// Output a \ref o2scl::table3d object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::table3d &h, std::string name);
  /// Input a \ref o2scl::table3d object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::table3d &h, std::string name="");
  /// Output a \ref o2scl::tensor_grid object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::tensor_grid &h, std::string name);
  /// Input a \ref o2scl::tensor_grid object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::tensor_grid &h, std::string name="");
  /// Output a \ref o2scl::expval_scalar object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_scalar &h, std::string name);
  /// Input a \ref o2scl::expval_scalar object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_scalar &h, std::string name="");
  /// Output a \ref o2scl::expval_vector object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_vector &h, std::string name);
  /// Input a \ref o2scl::expval_vector object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_vector &h, std::string name="");
  /// Output a \ref o2scl::expval_matrix object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_matrix &h, std::string name);
  /// Input a \ref o2scl::expval_matrix object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_matrix &h, std::string name="");
  /// Output a \ref o2scl::uniform_grid object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		  std::string name);
  /// Input a \ref o2scl::uniform_grid object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		 std::string name="");
  /// Output a vector of \ref o2scl::contour_line objects to a \ref hdf_file
  void hdf_output(hdf_file &hf, std::vector<o2scl::contour_line> &cl, 
		  std::string name);
  /// Input a vector of \ref o2scl::contour_line objects from a \ref hdf_file
  void hdf_input(hdf_file &hf, std::vector<o2scl::contour_line> &cl, 
		 std::string name="");
  /// Output a vector of \ref o2scl::edge_crossings objects to a \ref hdf_file
  void hdf_output(hdf_file &hf, std::vector<o2scl::edge_crossings> &ec, 
		  std::string name);
  /// Input a vector of \ref o2scl::edge_crossings objects from a \ref hdf_file
  void hdf_input(hdf_file &hf, std::vector<o2scl::edge_crossings> &ec, 
		 std::string name="");

}

#endif
