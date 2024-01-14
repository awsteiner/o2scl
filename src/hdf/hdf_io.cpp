/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2024, Andrew W. Steiner

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

void o2scl_hdf::index_spec() {
  return;
}

void o2scl_hdf::functions() {
  return;
}

void o2scl_hdf::hdf5_write_file(o2scl::table_units<> &t, std::string fn,
                                std::string name) {
  hdf_file hf;
  hf.open_or_create(fn);
  hdf_output(hf,t,name);
  hf.close();
  return;
}

int o2scl_hdf::value_spec(std::string spec, double &d,
			  int verbose, bool err_on_fail) {
  
  if (verbose>2) {
    std::cout << "Function value_spec is parsing: " << spec << std::endl;
  }
  
  if (spec.find(':')==std::string::npos) {

    if (verbose>1) {
      std::cout << "value_spec(): simple function or value "
		<< spec << std::endl;
    }
      
    d=o2scl::function_to_double(spec);
    
    return 0;
      
  } else if (spec.find("shell:")==0) {
	
    // Result from shell command
    
#ifdef HAVE_POPEN
    
    std::string cmd=spec.substr(6,spec.length()-6);
    cout << "Using shell command: " << cmd << endl;
    FILE *ps_pipe=popen(cmd.c_str(),"r");
    if (!ps_pipe) {
      if (err_on_fail) {
	O2SCL_ERR2("Pipe could not be opened in ",
		   "value_spec().",exc_efilenotfound);
      }
      return exc_efilenotfound;
    }
    
    // Variable 'cret' is unused, but put here to avoid
    // unused return value errors
    char line1[80];
    char *cret=fgets(line1,80,ps_pipe);
    if (verbose>0) {
      std::cout << "Shell command output is "
		<< line1 << std::endl;
    }

    if (pclose(ps_pipe)!=0) {
      if (err_on_fail) {
	O2SCL_ERR2("Pipe could not be closed in ",
		   "value_spec().",exc_efailed);
      }
      return 1;
    }
  
    // Read the output from the shell command
    std::string s=line1, t1;
    std::istringstream *ins=new std::istringstream(s);
    if (!((*ins) >> t1)) {
      return exc_efailed;
    }
    delete ins;

    // Finally, take the result string and convert to a double
    d=o2scl::stod(t1);

    return 0;
    
#else
    
    if (err_on_fail) {
      O2SCL_ERR2("Popen not supported in ",
		 "value_spec().",exc_efailed);
    }
    return 11;
    
#endif

  } else if (spec.find("python:")==0) {
	
    std::string cmd=spec.substr(7,spec.length()-7);
    cout << "Using python string: " << cmd << endl;
    string result;
    int ret=python_cmd_string(cmd,result,err_on_fail,200);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Function python_cmd_string() failed in ",
                   "value_spec().",exc_efailed);
      }
      return exc_efailed;
    }
    
    // Finally, take the result string and convert to a double
    d=o2scl::stod(result);

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
	  O2SCL_ERR2("No table column name or row index specified ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 6;
	}
      }
      size_t ncomma=addl_spec.find(',');
      if (ncomma==std::string::npos) {
	if (err_on_fail) {
	  O2SCL_ERR2("No comma found for hdf5:table ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 8;
	}
      }
      std::string col_name=addl_spec.substr(0,ncomma);
      std::string row_index=addl_spec.substr(ncomma+1,
					     addl_spec.length()-ncomma-1);
      if (col_name.length()==0 || row_index.length()==0) {
	if (err_on_fail) {
	  O2SCL_ERR2("Missing column name or row index in hdf5:table ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 9;
	}
      }
      o2scl::table_units<> t;
      o2scl_hdf::hdf_input(hf,t,obj_name);
      if (t.is_column(col_name)==false) {
	if (err_on_fail) {
	  O2SCL_ERR2("Can't find column hdf5:table ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 10;
	}
      }
      if (t.get_nlines()<=o2scl::stoszt(row_index)) {
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough lines in table in hdf5:table ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 11;
	}
      }
      d=t.get(col_name,o2scl::stoszt(row_index));
      
    } else if (type=="double[]") {
      
      size_t index=o2scl::stoszt(addl_spec);
      std::vector<double> vtemp;
      hf.getd_vec(obj_name,vtemp);
      d=vtemp[index];
      
    } else if (type=="int[]") {
      
      size_t index=o2scl::stoszt(addl_spec);
      std::vector<int> vtemp;
      hf.geti_vec(obj_name,vtemp);
      d=vtemp[index];
      
    } else if (type=="size_t[]") {
      
      size_t index=o2scl::stoszt(addl_spec);
      std::vector<size_t> vtemp;
      hf.get_szt_vec(obj_name,vtemp);
      d=vtemp[index];
      
    } else if (type=="uniform_grid<double>") {
      
      if (addl_spec.length()==0) {
	if (err_on_fail) {
	  O2SCL_ERR2("No table column name specified ",
		     "in value_spec().",o2scl::exc_einval);
	} else {
	  return 7;
	}
      }
      size_t index=o2scl::stoszt(addl_spec);
      o2scl::uniform_grid<double> ug;
      hdf_input(hf,ug,obj_name);
      std::vector<double> vtemp;
      ug.vector(vtemp);
      d=vtemp[index];
      
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

    return 0;
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

void o2scl_hdf::hdf_output(hdf_file &hf, o2scl::table<> &t, std::string name) {
  
  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,table<>,string).",exc_efailed);
  }
  
  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
  
  // Add typename
  hf.sets_fixed("o2scl_type","table");
  
  // Output table data
  hdf_output_data(hf,t);
  
  // Close table group
  hf.close_group(group);
  
  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_output(hdf_file &hf, o2scl::table_units<> &t, 
			   std::string name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,table_units<>,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
      
  // Add typename
  hf.sets_fixed("o2scl_type","table");

  // Output table_units data
  hdf_output_data(hf,t);
     
  // Close table_units group
  hf.close_group(group);
      
  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_output_data(hdf_file &hf, o2scl::table_units<> &t) {
      
  // Output base table object
  o2scl::table<> *tbase=dynamic_cast<o2scl::table_units<> *>(&t);
  if (tbase==0) {
    O2SCL_ERR2("Cast failed in hdf_output_data",
	       "(hdf_file &, table_units &).",o2scl::exc_efailed);
  }
  hdf_output_data(hf,*tbase);
      
  // Output unit info
  hf.seti("unit_flag",1);
      
  std::vector<std::string> units;
      
  // Restructure units
  for(size_t i=0;i<t.get_ncolumns();i++) {
    units.push_back(t.get_unit(t.get_column_name(i)));
  }
      
  // Output units
  hf.sets_vec_copy("units",units);

  return;
}

void o2scl_hdf::hdf_output_data(hdf_file &hf, o2scl::table<> &t) {

  // We're going to open a subgroup for data, so 
  // get the current ID to return to the main group later
  hid_t group=hf.get_current_id();

  // Restructure constants
  std::vector<std::string> cnames, cols;
  std::vector<double> cvalues;
      
  for(size_t i=0;i<t.get_nconsts();i++) {
    std::string name;
    double val;
    t.get_constant(i,name,val);
    cnames.push_back(name);
    cvalues.push_back(val);
  }
      
  // Output constants
  hf.sets_vec_copy("con_names",cnames);
  hf.setd_vec("con_values",cvalues);
      
  // Restructure column names
  for(size_t i=0;i<t.get_ncolumns();i++) {
    cols.push_back(t.get_column_name(i));
  }
      
  // Output column names
  hf.sets_vec_copy("col_names",cols);
      
  // Output unit info
  hf.seti("unit_flag",0);
      
  // Output number of lines
  hf.seti("nlines",((int)t.get_nlines()));

  // Output the interpolation type
  hf.set_szt("itype",t.get_interp_type());
      
  // Create data group
  hid_t group2=hf.open_group("data");
  hf.set_current_id(group2);
      
  if (t.get_nlines()>0) {
	
    // Output data
    for(size_t i=0;i<t.get_ncolumns();i++) {
      const std::vector<double> &col=t.get_column(t.get_column_name(i));
      // The actual vector is of size "maxlines", but we
      // only want to output the first "nlines" elements
      hf.setd_arr(t.get_column_name(i),t.get_nlines(),&(col[0]));
    }
	
  }

  // Close data group
  hf.close_group(group2);

  // Return ID to the main group
  hf.set_current_id(group);

  return;
}

void o2scl_hdf::hdf_output(hdf_file &hf, hist &h, std::string name) {
  
  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,hist,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","hist");

  // Add histogram
  hf.set_szt("size",h.size());
  hf.set_szt("rmode",h.get_rep_mode());
  hf.seti("extend_rhs",h.extend_rhs);
  hf.seti("extend_lhs",h.extend_lhs);
  hf.setd_vec_copy("bins",h.get_bins());
  hf.setd_vec_copy("weights",h.get_wgts());
  if (h.get_rep_mode()==hist::rmode_user) {
    hf.setd_vec_copy("ureps",h.user_rep);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(hdf_file &hf, hist &h, std::string &name) {
  
  // If no name specified, find name of first group of specified type
  if (name.length()==0) {
    hf.find_object_by_type("hist",name);
    if (name.length()==0) {
      O2SCL_ERR2("No object of type hist found in o2scl_hdf::hdf_",
		 "input(hdf_file &,hist &,string &).",exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Clear previous data
  h.clear();

  // Check typename
  string type;
  hf.gets_fixed("o2scl_type",type);
  if (type!="hist") {
    O2SCL_ERR2("Typename in HDF group does not match class in o2scl_hdf::",
	       "hdf_input(hdf_file &,hist &,string &).",exc_efailed);
  }

  size_t size2, mode;
  typedef std::vector<double> ubvector;
  ubvector bins, weights, ur;

  // Get histogram
  hf.get_szt("size",size2);
  hf.getd_vec("bins",bins);

  int itmp;

  // Get extend_lhs
  hf.geti_def("extend_lhs",0,itmp);
  if (itmp>1 || itmp<0) {
    O2SCL_ERR2("Failed to read extend_lhs in o2scl_hdf::",
	       "hdf_input(hdf_file &,hist &,string &).",exc_efailed);
  }
  if (itmp==1) h.extend_lhs=true;
  else h.extend_lhs=false;

  // Get extend_rhs
  hf.geti_def("extend_rhs",0,itmp);
  if (itmp>1 || itmp<0) {
    O2SCL_ERR2("Failed to read extend_rhs in o2scl_hdf::",
	       "hdf_input(hdf_file &,hist &,string &).",exc_efailed);
  }
  if (itmp==1) h.extend_rhs=true;
  else h.extend_rhs=false;

  // Set bins
  h.set_bin_edges(size2+1,bins);

  // Get weights
  hf.getd_vec("weights",weights);
  for(size_t i=0;i<size2;i++) {
    h.set_wgt_i(i,weights[i]);
  }

  // Get center mode
  hf.get_szt("rmode",mode);
  h.set_rep_mode(mode);

  if (h.get_rep_mode()==hist::rmode_user) {
    // Get reps
    hf.getd_vec("ureps",ur);
    h.set_reps(size2,ur);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  // Check that the histogram is valid
#if !O2SCL_NO_RANGE_CHECK
  h.is_valid();
#endif

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, hist &h, std::string name) {
  hdf_input_n(hf,h,name);
  return;
}


void o2scl_hdf::hdf_output(hdf_file &hf, const hist_2d &h, std::string name) {
  
  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,hist_2d,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","hist_2d");

  // Add histogram
  hf.set_szt("size_x",h.size_x());
  hf.set_szt("size_y",h.size_y());
  hf.set_szt("xrmode",h.get_x_rep_mode());
  hf.set_szt("yrmode",h.get_y_rep_mode());
  hf.setd_vec_copy("x_bins",h.get_x_bins());
  hf.setd_vec_copy("y_bins",h.get_y_bins());
  hf.setd_mat_copy("weights",h.get_wgts());
  hf.seti("extend_rhs",h.extend_rhs);
  hf.seti("extend_lhs",h.extend_lhs);

  if (h.get_x_rep_mode()==hist_2d::rmode_user) {
    hf.setd_vec_copy("uxreps",h.get_user_reps_x());
  }

  if (h.get_y_rep_mode()==hist_2d::rmode_user) {
    hf.setd_vec_copy("uyreps",h.get_user_reps_y());
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(hdf_file &hf, hist_2d &h, std::string &name) {

  // If no name specified, find name of first group of specified type
  if (name.length()==0) {
    hf.find_object_by_type("hist_2d",name);
    if (name.length()==0) {
      O2SCL_ERR2("No object of type hist_2d found in o2scl_hdf::hdf_",
		 "input(hdf_file &,hist_2d &,string &).",exc_efailed);
    }
  }
  
  h.clear();

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Check typename
  string type;
  hf.gets_fixed("o2scl_type",type);
  if (type!="hist_2d") {
    O2SCL_ERR2("Typename in HDF group does not match ",
	       "class in hdf_input().",exc_einval);
  }

  size_t size_x2, size_y2, mode_x, mode_y;
  typedef std::vector<double> ubvector;
  ubvector bins_x, bins_y, weights;

  // Get sizes
  hf.get_szt("size_x",size_x2);
  hf.get_szt("size_y",size_y2);

  // Get bins
  hf.getd_vec("x_bins",bins_x);
  hf.getd_vec("y_bins",bins_y);

  // Set bins
  h.set_bin_edges(bins_x.size(),bins_x,bins_y.size(),bins_y);

  // Get modes
  hf.get_szt("xrmode",mode_x);
  hf.get_szt("yrmode",mode_y);

  // Set modes
  h.set_rep_mode(mode_x,mode_y);
  
  // Get weights
  hf.getd_mat_copy("weights",h.get_wgts());

  typedef std::vector<double> ubvector;

  // Get x reps
  if (mode_x==hist_2d::rmode_user) {
    ubvector urx;
    hf.getd_vec("uxreps",urx);
    h.set_x_reps(urx.size(),urx);
  }

  // Get y reps
  if (mode_y==hist_2d::rmode_user) {
    ubvector ury;
    hf.getd_vec("uyreps",ury);
    h.set_y_reps(ury.size(),ury);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  // Check that the histogram is valid
#if !O2SCL_NO_RANGE_CHECK
  h.is_valid();
#endif

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, hist_2d &h, std::string name) {
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, const table3d &t, 
			   std::string name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,table3d,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
  
  // Add typename
  hf.sets_fixed("o2scl_type","table3d");

  string sx;
  hf.gets_fixed("o2scl_type",sx);

  // Restructure constants
  std::vector<std::string> cnames;
  std::vector<double> cvalues;
  
  for(size_t i=0;i<t.get_nconsts();i++) {
    string n;
    double val;
    t.get_constant(i,n,val);
    cnames.push_back(n);
    cvalues.push_back(val);
  }

  // Output constants
  hf.sets_vec_copy("con_names",cnames);
  hf.setd_vec("con_values",cvalues);
      
  // Restructure slice names
  std::vector<std::string> slnames;
  for(size_t i=0;i<t.get_nslices();i++) {
    slnames.push_back(t.get_slice_name(i));
  }

  // Output slice names
  hf.sets_vec_copy("slice_names",slnames);

  // Output unit info
  hf.seti("unit_flag",0);
  
  hf.seti("size_set",t.is_size_set());

  if (t.is_size_set()) {

    // Output grid
    hf.seti("numx",t.get_nx());
    hf.seti("numy",t.get_ny());
    
    hf.seti("xy_set",t.is_xy_set());
    
    if (t.is_xy_set()) {
      hf.sets("xname",t.get_x_name());
      hf.sets("yname",t.get_y_name());
      hf.setd_vec_copy("xval",t.get_x_data());
      hf.setd_vec_copy("yval",t.get_y_data());
    }
    
    // Create data group
    hid_t group2=hf.open_group("data");
    hf.set_current_id(group2);
    
    // Output data
    for(size_t i=0;i<t.get_nslices();i++) {
      hf.setd_mat_copy(t.get_slice_name(i),t.get_slice(i));
    }
    
    hf.close_group(group2);

  }

  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, table3d &t, 
                            std::string &name) {

  // If no name specified, find name of first group of specified type
  if (name.length()==0) {
    hf.find_object_by_type("table3d",name);
    if (name.length()==0) {
      O2SCL_ERR2("No object of type table3d found in o2scl_hdf::hdf_",
		 "input(hdf_file &,table3d &,string &).",exc_efailed);
    }
  }
  
  typedef std::vector<double> ubvector;

  // Clear previous data
  t.clear();

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Check typename
  string type2;
  hf.gets_fixed("o2scl_type",type2);
  if (type2!="table3d") {
    O2SCL_ERR2("Typename in HDF group does not match ",
	       "class in table::hdf_input().",exc_einval);
  }

  // Storage
  std::vector<std::string> cnames;
  ubvector cvalues;

  // Get constants
  hf.gets_vec_copy("con_names",cnames);
  hf.getd_vec("con_values",cvalues);
  for(size_t i=0;i<cnames.size();i++) {
    t.add_constant(cnames[i],cvalues[i]);
  }

  int size_set2, xy_set2;

  hf.geti("size_set",size_set2);

  if (size_set2>0) {
    
    int nx, ny;

    // Get grid
    hf.geti("numx",nx);
    hf.geti("numy",ny);
    
    hf.geti("xy_set",xy_set2);
    
    if (xy_set2>0) {
      
      string xname2, yname2;
      ubvector xv, yv;
      hf.gets("xname",xname2);
      hf.gets("yname",yname2);
      hf.getd_vec("xval",xv);
      hf.getd_vec("yval",yv);
      t.set_xy(xname2,nx,xv,yname2,ny,yv);

    } else {

      t.set_size(nx,ny);

    }

    // Get slice names
    std::vector<std::string> slnames;
    hf.gets_vec_copy("slice_names",slnames);
    for(size_t i=0;i<slnames.size();i++) {
      t.new_slice(slnames[i]);
    }

    // Open data group
    hid_t group2=hf.open_group("data");
    hf.set_current_id(group2);
    
    // Get data
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    double *d=new double[nx*ny];
    for(size_t i=0;i<t.get_nslices();i++) {
      hf.getd_mat_prealloc(t.get_slice_name(i),nx,ny,d);
      ubmatrix &m=t.get_slice(i);
      for(int ii=0;ii<nx;ii++) {
	for(int jj=0;jj<ny;jj++) {
	  m(ii,jj)=d[ii*ny+jj];
	}
      }
    }
    delete[] d;

    hf.close_group(group2);

  }

  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, table3d &h, std::string name) {
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, expval_scalar &sev,
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,expval_scalar,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","expval_scalar");

  // Add data
  hf.set_szt("iblock",sev.iblock);
  hf.set_szt("i",sev.i);
  hf.set_szt("nblocks",sev.nblocks);
  hf.set_szt("nperblock",sev.nperblock);
  hf.setd("current",sev.current);
  hf.setd_vec_copy("vals",sev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, expval_scalar &sev,
                            std::string &hdf_name) {
  
  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("expval_scalar",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type expval_scalar found in o2scl_hdf::hdf_",
		 "input(hdf_file &,expval_scalar &,string &).",exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get data
  hf.get_szt("iblock",sev.iblock);
  hf.get_szt("i",sev.i);
  hf.get_szt("nblocks",sev.nblocks);
  hf.get_szt("nperblock",sev.nperblock);
  hf.getd("current",sev.current);
  hf.getd_vec_copy("vals",sev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, expval_scalar &h,
                            std::string name) {
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, expval_vector &vev,
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,expval_vector,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","expval_vector");

  // Add data
  hf.set_szt("nvec",vev.nvec);
  hf.set_szt("iblock",vev.iblock);
  hf.set_szt("i",vev.i);
  hf.set_szt("nblocks",vev.nblocks);
  hf.set_szt("nperblock",vev.nperblock);
  hf.setd_vec_copy("current",vev.current);
  hf.setd_mat_copy("vals",vev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, expval_vector &vev,
                            std::string &hdf_name) {
  
  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("expval_vector",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type expval_vector found in o2scl_hdf::hdf_",
		 "input(hdf_file &,expval_vector &,string &).",exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get data
  hf.get_szt("nvec",vev.nvec);
  hf.get_szt("iblock",vev.iblock);
  hf.get_szt("i",vev.i);
  hf.get_szt("nblocks",vev.nblocks);
  hf.get_szt("nperblock",vev.nperblock);
  hf.getd_vec_copy("current",vev.current);
  hf.getd_mat_copy("vals",vev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, expval_vector &h,
                            std::string name) {
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, expval_matrix &mev, 
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,expval_matrix,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","expval_matrix");

  // Add data
  hf.set_szt("nr",mev.nr);
  hf.set_szt("nc",mev.nc);
  hf.set_szt("iblock",mev.iblock);
  hf.set_szt("i",mev.i);
  hf.set_szt("nblocks",mev.nblocks);
  hf.set_szt("nperblock",mev.nperblock);
  hf.setd_mat_copy("current",mev.current);
  hf.setd_ten("vals",mev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, expval_matrix &mev,
                            std::string &hdf_name) {

  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("expval_matrix",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type expval_matrix found in o2scl_hdf::hdf_",
		 "input(hdf_file &,expval_matrix &,string &).",exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get data
  hf.get_szt("nr",mev.nr);
  hf.get_szt("nc",mev.nc);
  hf.get_szt("iblock",mev.iblock);
  hf.get_szt("i",mev.i);
  hf.get_szt("nblocks",mev.nblocks);
  hf.get_szt("nperblock",mev.nperblock);
  hf.getd_mat_copy("current",mev.current);
  hf.getd_ten("vals",mev.vals);

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, expval_matrix &h,
                            std::string name) {
  
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf,
			   uniform_grid<double> &ug, 
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,uniform_grid<double>,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","uniform_grid<double>");

  // Set data
  hf.setd("start",ug.g_start);
  hf.setd("end",ug.g_end);
  hf.setd("width",ug.g_width);
  hf.set_szt("n_bins",ug.g_n_bins);
  if (ug.g_log) {
    hf.seti("log",1);
  } else {
    hf.seti("log",0);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, uniform_grid<double> &ug,
                            std::string &hdf_name) {

  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("uniform_grid<double>",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type uniform_grid found in o2scl_hdf::hdf_",
		 "input(hdf_file &,uniform_grid<double> &,string &).",
                 exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get data
  hf.getd("start",ug.g_start);
  hf.getd("end",ug.g_end);
  hf.getd("width",ug.g_width);
  hf.get_szt("n_bins",ug.g_n_bins);
  int tmp;
  hf.geti("log",tmp);
  if (tmp<=0) {
    ug.g_log=false;
  } else {
    ug.g_log=true;
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, uniform_grid<double> &h,
                            std::string name) {
  
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf,
			   const vector<contour_line> &cl, 
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,vector<contour_line>,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","vector<contour_line>");

  // Number of contour lines
  hf.set_szt("n_lines",cl.size());

  // Loop through each line
  for(size_t i=0;i<cl.size();i++) {

    // Create line group
    hid_t line_group=hf.open_group(((string)"line_")+itos(i));
    hf.set_current_id(line_group);

    // Output level and points
    hf.setd("level",cl[i].level);
    hf.setd_vec("x",cl[i].x);
    hf.setd_vec("y",cl[i].y);

    // Close line group
    hf.close_group(line_group);
    hf.set_current_id(group);

  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, vector<contour_line> &cl,
                            std::string &hdf_name) {
  
  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("vector<contour_line>",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type vector<contour_line> found in ",
		 "o2scl_hdf::hdf_input_n().",exc_efailed);
    }
  }

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get number of contour lines
  size_t n_lines;
  hf.get_szt("n_lines",n_lines);

  // Empty object
  contour_line empty;

  // Loop through each line
  for(size_t i=0;i<n_lines;i++) {

    // Open line group
    hid_t line_group=hf.open_group(((string)"line_")+itos(i));
    hf.set_current_id(line_group);

    // Push back a blank contour line
    cl.push_back(empty);
    
    // Input level and points
    hf.getd("level",cl[i].level);
    hf.getd_vec("x",cl[i].x);
    hf.getd_vec("y",cl[i].y);

    // Close line group
    hf.close_group(line_group);
    hf.set_current_id(group);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, vector<contour_line> &h,
                            std::string name) {
  
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf,
			   const vector<edge_crossings> &ec, 
			   std::string hdf_name) {

  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output",
	       "(hdf_file,vector<edge_crossings>,string).",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","vector<edge_crossings>");

  // Number of contour lines
  hf.set_szt("n",ec.size());

  // Loop through each line
  for(size_t i=0;i<ec.size();i++) {

    // Create line group
    hid_t ec_group=hf.open_group(((string)"ec_")+itos(i));
    hf.set_current_id(ec_group);

    // Output level and points
    hf.seti_mat_copy("status",ec[i].status);
    hf.setd_mat_copy("values",ec[i].values);

    // Close line group
    hf.close_group(ec_group);
    hf.set_current_id(group);

  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, vector<edge_crossings> &ec,
                            std::string &hdf_name) {

  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_object_by_type("vector<edge_crossings>",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type vector<edge_crossings> found in ",
		 "o2scl_hdf::hdf_input().",
                 exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Get number of contour lines
  size_t n;
  hf.get_szt("n",n);

  // Empty object
  edge_crossings empty;

  // Loop through each line
  for(size_t i=0;i<n;i++) {

    // Open line group
    hid_t ec_group=hf.open_group(((string)"ec_")+itos(i));
    hf.set_current_id(ec_group);

    // Push back a blank contour line
    ec.push_back(empty);
    
    // Input level and points
    hf.geti_mat_copy("status",ec[i].status);
    hf.getd_mat_copy("values",ec[i].values);

    // Close line group
    hf.close_group(ec_group);
    hf.set_current_id(group);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, vector<edge_crossings> &h,
                            std::string name) {
  
  hdf_input_n(hf,h,name);
  return;
}

void o2scl_hdf::hdf_output(hdf_file &hf,
			   o2scl::tensor_grid<std::vector<double>,
			   std::vector<size_t>> &t, std::string name) {

  // Check that tensor_grid object is valid
  t.is_valid();
  
  if (hf.has_write_access()==false) {
    O2SCL_ERR2("File not opened with write access in hdf_output(hdf_file,",
	       "tensor_grid<vector<double>,vector<size_t> >,string).",
	       exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
      
  // Add typename
  hf.sets_fixed("o2scl_type","tensor_grid");
      
  // Add rank
  hf.seti("rank",t.get_rank());
      
  // Add dimensions
  std::vector<int> size_arr;
  for(size_t i=0;i<t.get_rank();i++) {
    size_arr.push_back(t.get_size(i));
  }
  hf.seti_vec("size",size_arr);
      
  // Add data
  const std::vector<double> &d=t.get_data();
  const double *ptr=&d[0];
  hf.setd_arr("data",d.size(),ptr);
      
  // Add grid
  if (t.is_grid_set()) hf.seti("grid_set",1);
  else hf.seti("grid_set",0);
      
  if (t.is_grid_set()) {
    std::vector<double> grid2;
    for(size_t j=0;j<t.get_rank();j++) {
      for(size_t k=0;k<t.get_size(j);k++) {
	grid2.push_back(t.get_grid(j,k));
      }
    }
    hf.setd_vec("grid",grid2);
  }
      
  // Close group
  hf.close_group(group);
      
  // Return location to previous value
  hf.set_current_id(top);
      
  return;
}
  
void o2scl_hdf::hdf_input_n(hdf_file &hf,
                            o2scl::tensor_grid<std::vector<double>,
                            std::vector<size_t>> &t, std::string &name) {
    
  // If no name specified, find name of first group of specified type
  if (name.length()==0) {
    hf.find_object_by_type("tensor_grid",name);
    if (name.length()==0) {
      O2SCL_ERR2("No object of type tensor_grid found in o2scl_hdf::hdf_",
		 "input(hdf_file &,tensor_grid<vector<double> > &,string &).",
                 exc_efailed);
    }
  }
  
  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);
      
  // Check typename
  std::string type;
  hf.gets_fixed("o2scl_type",type);
  if (type!="tensor_grid") {
    O2SCL_ERR2("Typename in HDF group does not match ",
	       "class in hdf_input().",o2scl::exc_einval);
  }
      
  // Get rank
  int rank;
  hf.geti("rank",rank);
  int *size_i=new int[rank];
  hf.geti_vec_prealloc("size",rank,size_i);
      
  size_t *size_s=new size_t[rank];
  for(int k=0;k<rank;k++) size_s[k]=size_i[k];
  t.resize(rank,size_s);
      
  delete[] size_i;
  delete[] size_s;

  // Get non-const pointer to first element
  vector<size_t> zero(rank);
  for(int k=0;k<rank;k++) zero[k]=0;
  double *start=&t.get(zero);
  
  // The tensor_grid class really doesn't allow one to obtain a
  // non-const vector, so we just use a const pointer instead
  hf.getd_arr("data",t.total_size(),start);
  //hf.getd_vec("data",t.get_data());
  
  // Get grid
  bool grid_set2=false;
  int igrid_set;
  hf.geti("grid_set",igrid_set);
  if (igrid_set>0) grid_set2=true;
  if (grid_set2) {
    std::vector<double> ogrid;
    hf.getd_vec("grid",ogrid);
    size_t ix=0;
    double **grid2=new double *[rank];
    for(size_t j=0;j<((size_t)rank);j++) {
      grid2[j]=new double[t.get_size(j)];
      for(size_t k=0;k<t.get_size(j);k++) {
	grid2[j][k]=ogrid[ix];
	ix++;
      }
    }
    t.set_grid_packed(ogrid);
    for(size_t j=0;j<((size_t)rank);j++) {
      delete[] grid2[j];
    }
    delete[] grid2;
  }
      
  // Close group
  hf.close_group(group);
      
  // Return location to previous value
  hf.set_current_id(top);

  // Check that tensor_grid object is valid
  t.is_valid();
  
  return;
}

void o2scl_hdf::hdf_input(hdf_file &hf, tensor_grid<vector<double> > &h,
                          std::string name) {
  
  hdf_input_n(hf,h,name);
  
  return;
}

std::vector<double> o2scl_hdf::vector_spec(std::string spec) {
  std::vector<double> v;
  vector_spec<std::vector<double> >(spec,v);
  return v;
}

std::vector<std::vector<double>>
o2scl_hdf::mult_vector_spec(std::string spec) {
  std::vector<std::vector<double>> v;
  mult_vector_spec<std::vector<double> >(spec,v);
  return v;
}
