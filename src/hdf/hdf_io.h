/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2019, Andrew W. Steiner

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

#include <fnmatch.h>

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

/** \brief The \o2 namespace for I/O with HDF
 */
namespace o2scl_hdf {

  /** \brief Input a \ref o2scl::prob_dens_mdim_amr object from a 
      \ref hdf_file
  */
  template<class vec_t, class mat_t> 
    void hdf_input(hdf_file &hf,
		   o2scl::prob_dens_mdim_amr<vec_t,mat_t> &p,
		   std::string name) {
      
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

    // Close table_units group
    hf.close_group(group);
    
    // Return location to previous value
    hf.set_current_id(top);
    
    return;
  }
  
  /** \brief Output a \ref o2scl::table object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table<> &t, std::string name);

#ifndef O2SCL_NO_HDF_INPUT  
  /** \brief Input a \ref o2scl::table object from a \ref hdf_file

      \comment
      Note that a default value is not allowed here because this
      is a template function
      \endcomment
  */
  template<class vec_t> 
    void hdf_input(hdf_file &hf, o2scl::table<vec_t> &t, std::string name) {
      
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

    // Input the table data
    hdf_input_data(hf,t);

    // Close group
    hf.close_group(group);

    // Return location to previous value
    hf.set_current_id(top);

    t.check_synchro();

    return;
  }
#endif

  /** \brief Internal function for outputting a \ref o2scl::table object
   */
  void hdf_output_data(hdf_file &hf, o2scl::table<> &t);

  /** \brief Internal function for inputting a \ref o2scl::table object
   */
  template<class vec_t> 
    void hdf_input_data(hdf_file &hf, o2scl::table<vec_t> &t) {
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
    hf.gets_vec("con_names",cnames);
    hf.getd_vec_copy("con_values",cvalues);
    if (cnames.size()!=cvalues.size()) {
      O2SCL_ERR2("Size mismatch between constant names and values ",
		 "in o2scl_hdf::hdf_input_data().",o2scl::exc_einval);
    }
    for(size_t i=0;i<cnames.size();i++) {
      t.add_constant(cnames[i],cvalues[i]);
    }

    // Get column names
    hf.gets_vec("col_names",cols);
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
    t.check_synchro();

    return;
  }
  
  /** \brief Output a \ref o2scl::table_units object to a \ref hdf_file
   */
  void hdf_output(hdf_file &hf, o2scl::table_units<> &t, 
		  std::string name);

  /** \brief Input a \ref o2scl::table_units object from a \ref hdf_file

      \comment
      Note that a default value is not allowed here because this
      is a template function
      \endcomment
  */
  template<class vec_t> 
    void hdf_input(hdf_file &hf, o2scl::table_units<vec_t> &t, 
		   std::string name) {
      
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

    t.check_synchro();
    
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
      hf.gets_vec("units",units);
      for(size_t i=0;i<units.size();i++) {
	t.set_unit(t.get_column_name(i),units[i]);
      }
    }

    return;
  }
  
  /// Output a \ref o2scl::hist object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::hist &h, std::string name);
  /// Input a \ref o2scl::hist object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::hist &h, std::string name="");
  /// Output a \ref o2scl::hist_2d object to a \ref hdf_file
  void hdf_output(hdf_file &hf, const o2scl::hist_2d &h, std::string name);
  /// Input a \ref o2scl::hist_2d object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::hist_2d &h, std::string name="");
  /// Output a \ref o2scl::table3d object to a \ref hdf_file
  void hdf_output(hdf_file &hf, const o2scl::table3d &h, std::string name);
  /// Input a \ref o2scl::table3d object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::table3d &h, std::string name="");
  /// Output a \ref o2scl::expval_scalar object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_scalar &h,
		  std::string name);
  /// Input a \ref o2scl::expval_scalar object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_scalar &h,
		 std::string name="");
  /// Output a \ref o2scl::expval_vector object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_vector &h,
		  std::string name);
  /// Input a \ref o2scl::expval_vector object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_vector &h, std::string name="");
  /// Output a \ref o2scl::expval_matrix object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::expval_matrix &h,
		  std::string name);
  /// Input a \ref o2scl::expval_matrix object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::expval_matrix &h, std::string name="");
  /// Output a \ref o2scl::uniform_grid object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		  std::string name);
  /// Input a \ref o2scl::uniform_grid object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::uniform_grid<double> &h, 
		 std::string name="");
  /// Output a vector of \ref o2scl::contour_line objects to a \ref hdf_file
  void hdf_output(hdf_file &hf, const std::vector<o2scl::contour_line> &cl, 
		  std::string name);
  /// Input a vector of \ref o2scl::contour_line objects from a \ref hdf_file
  void hdf_input(hdf_file &hf, std::vector<o2scl::contour_line> &cl, 
		 std::string name="");
  /// Output a vector of \ref o2scl::edge_crossings objects to a \ref hdf_file
  void hdf_output(hdf_file &hf, const std::vector<o2scl::edge_crossings> &ec, 
		  std::string name);
  /// Input a vector of \ref o2scl::edge_crossings objects from a \ref hdf_file
  void hdf_input(hdf_file &hf, std::vector<o2scl::edge_crossings> &ec, 
		 std::string name="");
  /// Output a \ref o2scl::tensor_grid object to a \ref hdf_file
  void hdf_output(hdf_file &hf, o2scl::tensor_grid<std::vector<double>,
		  std::vector<size_t> > &t, std::string name);
  /// Input a \ref o2scl::tensor_grid object from a \ref hdf_file
  void hdf_input(hdf_file &hf, o2scl::tensor_grid<std::vector<double>,
		 std::vector<size_t> > &t, std::string name="");

  /** \brief A value specified by a string
      
      Formats:
      - func: \<function\>
      - HDF5 object in file: 
      hdf5:\<file name\>:\<object name\>:[additional specification]

      \note unfinished.
  */
  template<class vec_t> int value_spec(std::string spec,
					  double &d,
					  int verbose=0,
					  bool err_on_fail=true) {
    if (verbose>2) {
      std::cout << "Function vector_spec is parsing: " << spec << std::endl;
    }
      
    if (spec.find("func:")==0) {

      std::string temp=spec.substr(5,spec.length()-5);
      if (verbose>1) {
	std::cout << "vector_spec(): single value " << temp
		  << std::endl;
      }
      
      o2scl::calculator calc;
      std::map<std::string,double> vars;
      calc.compile(temp.c_str(),&vars);

      d=calc.eval(&vars);
      return 0;

    } else if (spec.find("hdf5:")==0) {
	
      // HDF5 object in a file
      if (verbose>1) {
	std::cout << "value_spec(): HDF5 file " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in value_spec().",o2scl::exc_einval);
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
		     "in value_spec().",o2scl::exc_einval);
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
	std::cout << "Object name " << obj_name << std::endl;
	std::cout << "Additional specification " << addl_spec << std::endl;
      }
      o2scl_hdf::hdf_file hf;
	
      std::string fname_old=fname;
      std::vector<std::string> matches;
      int wret=o2scl::wordexp_wrapper(fname_old,matches);
      if (matches.size()>1 || matches.size()==0 || wret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function wordexp_wrapper() failed ",
		     "in value_spec().",o2scl::exc_einval);
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
      int find_ret=hf.find_object_by_name(obj_name,type);
      if (find_ret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Object not found in file ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 11;
	}
      }
      if (verbose>1) {
	std::cout << "Object type from file: " << type << std::endl;
      }
	
      if (type=="table") {
	if (addl_spec.length()==0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("No table column name specified ",
		       "in value_spec().",o2scl::exc_einval);
	  } else {
	    return 6;
	  }
	}
	/*
	  o2scl::table_units<> t;
	  o2scl_hdf::hdf_input(hf,t,obj_name);
	  v.resize(t.get_nlines());
	  for(size_t i=0;i<t.get_nlines();i++) {
	  v[i]=t.get(addl_spec,i);
	  }
	*/
      } else if (type=="double[]") {
	/*
	  std::vector<double> vtemp;
	  hf.getd_vec(obj_name,vtemp);
	  v.resize(vtemp.size());
	  for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	  }
	*/
      } else if (type=="hist") {
	/*
	  o2scl::hist ht;
	  hdf_input(hf,ht,obj_name);
	  typedef boost::numeric::ublas::vector<double> ubvector;
	  const ubvector &wgts=ht.get_wgts();
	  v.resize(wgts.size());
	  for(size_t i=0;i<v.size();i++) {
	  v[i]=wgts[i];
	  }
	*/
      } else if (type=="int[]") {
	/*
	  std::vector<int> vtemp;
	  hf.geti_vec(obj_name,vtemp);
	  v.resize(vtemp.size());
	  for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	  }
	*/
      } else if (type=="size_t[]") {
	/*
	std::vector<size_t> vtemp;
	hf.get_szt_vec(obj_name,vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	*/
      } else if (type=="uniform_grid<double>") {
	/*
	o2scl::uniform_grid<double> ug;
	hdf_input(hf,ug,obj_name);
	std::vector<double> vtemp;
	ug.vector(vtemp);
	v.resize(vtemp.size());
	for(size_t i=0;i<v.size();i++) {
	  v[i]=vtemp[i];
	}
	*/
      } else if (type=="int") {
	int itemp;
	hf.geti(obj_name,itemp);
	d=itemp;
      } else if (type=="double") {
	double dtemp;
	hf.getd(obj_name,dtemp);
	d=dtemp;
      } else if (type=="size_t") {
	size_t szttemp;
	hf.get_szt(obj_name,szttemp);
	d=szttemp;
      }
      hf.close();
	
    }
    
    if (verbose>0) {
      std::cout << "Could not understand prefix in value_spec()."
		<< std::endl;
    }
    
    if (err_on_fail) {
      O2SCL_ERR2("Could not parse specification ",
		 "in value_spec().",o2scl::exc_einval);
    }
    
    return 1;
  }
  
  /** \brief A vector specified by a string
      
      Formats:
      - single value: val:\<value\>
      - list of values: list:\<entry 0\>,\<entry 1\>, ...,\<entry n-1\>
      - function: func:\<N\>:\<function of i\>
      - grid: grid:\<begin\>:\<end\>:\<width\>:["log"]
      - column in text file: text:\<filename\>:\<column\> 
      - HDF5 object in file: 
      hdf5:\<file name\>:\<object name\>:[additional specification]
      
      Additional specifications
      - table: \<column\>
  */
  template<class vec_t> int vector_spec(std::string spec, vec_t &v,
					int verbose=0,
					bool err_on_fail=true) {

    if (verbose>2) {
      std::cout << "Function vector_spec is parsing: " << spec << std::endl;
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
	std::cout << n << " " << sv[0] << " " << sv[n-1] << std::endl;
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
	std::cout << "Size " << n << std::endl;
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
	std::cout << "Function " << func << std::endl;
      }
      o2scl::calculator calc;
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
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough information for grid ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 7;
	}
      }
      if (verbose>1) {
	std::cout << "Begin,end,width "
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
      size_t col=o2scl::stoszt(sv[2]);
      if (verbose>1) {
	std::cout << "Filename,column " << sv[1] << " " << col << std::endl;
      }
      if (col==0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Column is zero for text file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 13;
	}
      }

      bool in_header=true;
      std::ifstream fin;
      std::string line, word;
      fin.open(sv[1].c_str());

      do {
	getline(fin,line);
	std::istringstream is(line);
	size_t i=0;
	for(;i<col && (is >> word) && in_header;i++) {
	  if (i==col-1 && o2scl::is_number(word)) {
	    in_header=false;
	  }
	  if (i==col-1 && verbose>2) {
	    if (in_header) {
	      std::cout << "Word: " << word << " header." << std::endl;
	    } else {
	      std::cout << "Word: " << word << " start of data."
			<< std::endl;
	    }
	  }
	}
      } while (in_header && !fin.eof());

      if (in_header) {
	if (err_on_fail) {
	  O2SCL_ERR2("Couldn't find a number in text file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 13;
	}
      }

      std::vector<double> tempv;
      
      bool end_of_data=false;
      do {
	tempv.push_back(o2scl::stod(word));
	
	getline(fin,line);
	std::istringstream is(line);
	for(size_t i=0;i<col;i++) {
	  if (!(is >> word)) {
	    i=col;
	    end_of_data=true;
	  }
	}
	if (!o2scl::is_number(word)) end_of_data=true;
	
      } while (!end_of_data && !fin.eof());
      
      fin.close();

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
	std::cout << "Filename " << fname << std::endl;
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
	std::cout << "Object name " << obj_name << std::endl;
	std::cout << "Additional specification " << addl_spec << std::endl;
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
	std::cout << "Filename after wordexp() " << fname << std::endl;
      }
	
      hf.open(fname);
      std::string type;
      int find_ret=hf.find_object_by_name(obj_name,type);
      if (find_ret!=0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Object not found in file ",
		     "in vector_spec().",o2scl::exc_einval);
	} else {
	  return 11;
	}
      }
      if (verbose>1) {
	std::cout << "Object type from file: " << type << std::endl;
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
	o2scl::table_units<> t;
	o2scl_hdf::hdf_input(hf,t,obj_name);
	v.resize(t.get_nlines());
	for(size_t i=0;i<t.get_nlines();i++) {
	  v[i]=t.get(addl_spec,i);
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

  /** \brief Return a std vector specified by a string
   */
  std::vector<double> vector_spec(std::string spec);
  
  /** \brief A list of vectors specified by a string

      Formats:
      - function: func:\<num\>:\<len(i)\>:\<function of j\>
      - columns in text file: text:\<column list\> 
      - HDF5 object(s) in file(s): 
      hdf5:\<file name(s)\>:\<object name(s)\>:[additional specification]

      \note unfinished.
  */
  template<class vec_t> int mult_vector_spec(std::string spec,
					     std::vector<vec_t> &v,
					     int verbose=0,
					     bool err_on_fail=true) {
    
    if (spec.find("func:")==0) {
	
      // Function
      if (verbose>1) {
	std::cout << "mult_vector_spec(): Function " << spec << std::endl;
      }
      std::string temp=spec.substr(5,spec.length()-5);
      size_t ncolon=temp.find(':');
      if (ncolon==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("Function specified but no array length specified ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 1;
	}
      }
      size_t n;
      int cret=o2scl::stoszt_nothrow(temp.substr(0,ncolon),n);
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
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent function specified ",
		     "in mult_vector_spec().",o2scl::exc_einval);
	} else {
	  return 3;
	}
      }
      std::string func=temp.substr(ncolon+1,temp.length()-ncolon-1);
      if (verbose>1) {
	std::cout << "Function " << func << std::endl;
      }
      o2scl::calculator calc;
      std::map<std::string,double> vars;
      calc.compile(func.c_str(),&vars);
      v.resize(n);
      for(size_t i=0;i<n;i++) {
	vars["i"]=((double)i);
	v[i]=calc.eval(&vars);
      }
	
    } else if (spec.find("hdf5:")==0) {
	
      // HDF5 object in a file
      if (verbose>1) {
	std::cout << "mult_vector_spec(): HDF5 file " << spec << std::endl;
      }
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
      std::string fname=temp.substr(0,ncolon);
      if (verbose>1) {
	std::cout << "Filename " << fname << std::endl;
      }
      if (temp.length()<ncolon+1) {
	if (err_on_fail) {
	  O2SCL_ERR2("No apparent object name specified ",
		     "in mult_vector_spec().",o2scl::exc_einval);
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
	std::cout << "Object name " << obj_name << std::endl;
	std::cout << "Additional specification " << addl_spec << std::endl;
      }
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
      
      for(size_t ifile=0;ifile<nfiles;ifile++) {
	fname=matches[0];
	if (verbose>1) {
	  std::cout << "Filename for index " << ifile << " is " << fname
		    << std::endl;
	}
	
	hf.open(fname);
	std::string type;
	int find_ret=hf.find_object_by_name(obj_name,type);
	if (find_ret!=0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Object not found in file ",
		       "in mult_vector_spec().",o2scl::exc_einval);
	  } else {
	    return 11;
	  }
	}
	if (verbose>1) {
	  std::cout << "Object type from file: " << type << std::endl;
	}
	
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
	  
	  for(size_t j=0;j<t.get_ncolumns();j++) {
	    if (fnmatch(addl_spec.c_str(),
			t.get_column_name(j).c_str(),0)==0) {
	      vec_t vtemp(t.get_nlines());
	      for(size_t k=0;k<t.get_nlines();k++) {
		vtemp[k]=t.get(j,k);
	      }
	      v.push_back(vtemp);
	    }
	  }
	} else if (type=="double[]") {
	  std::vector<double> vtemp;
	  hf.getd_vec(obj_name,vtemp);
	  v.resize(vtemp.size());
	  for(size_t i=0;i<v.size();i++) {
	    v[i]=vtemp[i];
	    std::cout << "X: " << i << " " << v[i] << std::endl;
	  }
	}
	hf.close();

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
  
}

#endif
