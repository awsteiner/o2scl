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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>
#include <o2scl/cursesw.h>
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/interp_krige.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interpm_python.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/set_python.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

/*
  int acol_manager::comm_generic()
  int acol_manager::comm_get_grid()
  int acol_manager::comm_get_row()
  int acol_manager::comm_get_unit()
  int acol_manager::comm_help()
  int acol_manager::comm_h5_copy()
  int acol_manager::comm_index()
  int acol_manager::comm_insert()
  int acol_manager::comm_insert_full()
  int acol_manager::comm_integ()
  int acol_manager::comm_internal()
  int acol_manager::comm_interp()
  int acol_manager::comm_interp_type()
  int acol_manager::comm_interactive()
*/

int acol_manager::comm_generic(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype;

  // Delete previous object
  command_del(type);
  clear_obj();
  
  int ret=get_input_one(sv,"Enter type of object to create",ctype,"create",
			itive_com);
  if (ret!=0) return ret;

  vector<string> sv2=sv;
  vector<string>::iterator it=sv2.begin();
  sv2.erase(it+1);

  string fname;
  if (sv2[1]!=((std::string)"cin")) {
    string fname_old=sv2[1];
    std::vector<std::string> matches;
    int wret=wordexp_wrapper(fname_old,matches);
    if (matches.size()>1 || matches.size()==0 || wret!=0) {
      cerr << "Function wordexp_wrapper() returned non-zero value. "
	   << "Bad filename?" << endl;
      return 1;
    }
    fname=matches[0];
    if (verbose>1) {
      cout << "Function wordexp() converted "
	   << fname_old << " to " << fname << endl;
    }
  } else {
    fname="cin";
  }
  
  // Open the file 
  ifstream ifs;
  if (fname!=((std::string)"cin")) {
    ifs.open(fname.c_str());
    if (!(ifs)) {
      cerr << "Read of file named '" << fname
	   << "' failed. Non-existent file?" << endl;
      return exc_efailed;
    }
  }

  istream *istr;
  
  //if (fname!=((std::string)"cin")) {
  if (fname!=((std::string)"cin")) {
    istr=&ifs;
  } else {
    istr=&cin;
  }
  
  if (ctype=="table") {
    
    table_obj.read_generic(*istr,verbose);

  } else if (ctype=="table3d") {
    
    table3d_obj.read_gen3_list(*istr,verbose);
    
  } else if (ctype=="uniform_grid<double>") {

    double start, end, width;
    size_t n_bins;
    std::string slog;
    (*istr) >> n_bins >> start >> end >> width >> slog;
    bool is_log=o2scl::stob(slog);
    if (is_log) {
      uniform_grid_end_width<double> ug(start,end,width);
      ug_obj=ug;
    } else {
      uniform_grid_log_end_width<double> ug(start,end,width);
      ug_obj=ug;
    }
    
  } else if (ctype=="hist") {

    hist_obj.clear();
    vector<double> edges, reps;
    double dtemp;
    bool edge_is_next=true;
    while ((*istr) >> dtemp) {
      if (edge_is_next) {
        edges.push_back(dtemp);
      } else {
        reps.push_back(dtemp);
      }
      edge_is_next=!edge_is_next;
    }
    if (edges.size()!=reps.size()+1) {
      cerr << "Command generic failed to read type hist." << endl;
      return 1;
    }
    hist_obj.set_bin_edges(edges.size(),edges);
    hist_obj.set_reps(reps.size(),reps);
    
  } else if (ctype=="prob_dens_mdim_gmm") {
    
    pgmm_obj.read_generic(*istr);
    
  } else if (ctype=="prob_dens_mdim_gaussian") {
    
    pdmg_obj.read_generic(*istr);
    
  } else if (ctype=="int") {

    (*istr) >> int_obj;
    
  } else if (ctype=="char") {

    (*istr) >> char_obj;
    
  } else if (ctype=="double") {

    (*istr) >> double_obj;
    
  } else if (ctype=="size_t") {

    (*istr) >> size_t_obj;
    
  } else if (ctype=="string") {

    getline(*istr,string_obj);
    
  } else if (ctype=="int[]") {

    int itmp;
    intv_obj.clear();
    while ((*istr) >> itmp) {
      intv_obj.push_back(itmp);
    }
    
  } else if (ctype=="int[]-line") {

    std::string row;
    getline((*istr),string_obj);
    int d;
    istringstream ins(row);
    intv_obj.clear();
    while (ins >> d) {
      intv_obj.push_back(d);
      return 0;
    }
    
  } else if (ctype=="int[]-n") {

    intv_obj.clear();
    size_t n;
    (*istr) >> n;
    intv_obj.resize(n);
    for(size_t i=0;i<n;i++) {
      (*istr) >> intv_obj[i];
    }
    
  } else if (ctype=="double[]") {
    
    double dtmp;
    doublev_obj.clear();
    while ((*istr) >> dtmp) {
      doublev_obj.push_back(dtmp);
    }
    
  } else if (ctype=="double[]-line") {

    std::string row;
    getline(*istr,string_obj);
    double d;
    istringstream ins(row);
    doublev_obj.clear();
    while (ins >> d) {
      doublev_obj.push_back(d);
      return 0;
    }
    
  } else if (ctype=="hist_2d") {

    std::string row;
    double d;
    vector<double> x_grid, y_grid;
    
    getline(*istr,row);
    istringstream ins(row);
    while (ins >> d) {
      x_grid.push_back(d);
    }
    
    getline(*istr,row);
    istringstream ins2(row);
    while (ins2 >> d) {
      y_grid.push_back(d);
    }
    
    hist_2d_obj.set_bin_edges(x_grid.size(),x_grid,y_grid.size(),y_grid);
    
    for(size_t i=0;i<hist_2d_obj.size_x();i++) {
      for(size_t j=0;j<hist_2d_obj.size_y();j++) {
        (*istr) >> d;
        hist_2d_obj.set_wgt_i(i,j,d);
      }
    }
    
  } else if (ctype=="vector<contour_line>") {

    size_t n;
    (*istr) >> n;
    cont_obj.clear();
    for(size_t k=0;k<n;k++) {
      contour_line cline;
      size_t m;
      (*istr) >> cline.level >> m;
      cline.x.resize(m);
      cline.y.resize(m);
      for(size_t i=0;i<m;k++) {
        (*istr) >> cline.x[i] >> cline.y[i];
      }
      cont_obj.push_back(cline);
    }

  } else if (ctype=="vec_vec_string") {

    size_t m, n;
    (*istr) >> m;
    vvstring_obj.resize(m);
    for(size_t i=0;i<m;i++) {
      (*istr) >> n;
      vvstring_obj[i].resize(n);
      for(size_t j=0;j<n;j++) {
        getline((*istr),vvstring_obj[i][j]);
      }
    }
    
  } else if (ctype=="vec_vec_double") {

    size_t m, n;
    (*istr) >> m;
    vvdouble_obj.resize(m);
    for(size_t i=0;i<m;i++) {
      (*istr) >> n;
      vvdouble_obj[i].resize(n);
      for(size_t j=0;j<n;j++) {
        (*istr) >> vvdouble_obj[i][j];
      }
    }
    
  } else if (ctype=="vvd-equal-cols") {

    vvdouble_obj.clear();
    std::string row;
    
    while (getline((*istr),row)) {
      stringv_obj.push_back(row);
      double d;
      istringstream ins(row);
      vector<double> vd;
      while (ins >> d) {
        vd.push_back(d);
      }
      vvdouble_obj.push_back(vd);
    }
    
  } else if (ctype=="double[]-n") {

    doublev_obj.clear();
    size_t n;
    (*istr) >> n;
    doublev_obj.resize(n);
    for(size_t i=0;i<n;i++) {
      (*istr) >> doublev_obj[i];
    }
    
  } else if (ctype=="size_t[]") {
    
      size_t sttmp;
      size_tv_obj.clear();
      while ((*istr) >> sttmp) {
	size_tv_obj.push_back(sttmp);
      }
    
  } else if (ctype=="size_t[]-line") {

    std::string row;
    getline((*istr),string_obj);
    size_t d;
    istringstream ins(row);
    size_tv_obj.clear();
    while (ins >> d) {
      size_tv_obj.push_back(d);
      return 0;
    }
    
  } else if (ctype=="size_t[]-n") {

    size_tv_obj.clear();
    size_t n;
    (*istr) >> n;
    size_tv_obj.resize(n);
    for(size_t i=0;i<n;i++) {
      (*istr) >> size_tv_obj[i];
    }
    
  } else if (ctype=="string[]") {

    std::string stmp;
    while (getline((*istr),stmp)) {
      stringv_obj.push_back(stmp);
    }
    
  } else if (ctype=="string[]-line") {

    std::string row;
    getline((*istr),string_obj);
    string d;
    istringstream ins(row);
    stringv_obj.clear();
    while (ins >> d) {
      stringv_obj.push_back(d);
      return 0;
    }
    
  } else if (ctype=="string[]-n") {

    stringv_obj.clear();
    size_t n;
    string tmp;
    getline((*istr),tmp);
    n=o2scl::stoi(tmp);
    stringv_obj.resize(n);
    for(size_t i=0;i<n;i++) {
      getline((*istr),stringv_obj[i]);
    }
    
  } else {

    cerr << "Cannot read generic text file for object of type "
	 << ctype << endl;
    return 1;
    
  }

  command_add(ctype);
  type=ctype;
  
  if (fname!=((std::string)"cin")) {
    ifs.close();
  }

  return 0;
}

int acol_manager::comm_get_grid(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    vector<vector<string> > string_mat(3);
    vector<int> align_spec(3);
    
    size_t max_size=table3d_obj.get_nx();
    if (table3d_obj.get_ny()>max_size) {
      max_size=table3d_obj.get_ny();
    }

    // The first column which enumerates the grid points
    align_spec[0]=columnify::align_left;
    string_mat[0].resize(max_size+1);
    for(size_t ell=0;ell<max_size;ell++) {
      string_mat[0][ell+1]=o2scl::szttos(ell)+".";
    }

    // The first row which labels the grids
    string_mat[1].resize(max_size+1);
    align_spec[1]=columnify::align_right;
    string_mat[1][0]=table3d_obj.get_x_name();
    string_mat[2].resize(max_size+1);
    align_spec[2]=columnify::align_right;
    string_mat[2][0]=table3d_obj.get_y_name();
    
    // Now the grid data
    for(size_t ell=0;ell<max_size;ell++) {
      if (ell<table3d_obj.get_nx()) {
	string_mat[1][ell+1]=o2scl::dtos(table3d_obj.get_grid_x(ell),
                                         precision);
      }
      if (ell<table3d_obj.get_ny()) {
	string_mat[2][ell+1]=o2scl::dtos(table3d_obj.get_grid_y(ell),
                                         precision);
      }
    }

    columnify col;
    vector<string> aligned(max_size+1);
    col.align(string_mat,3,max_size+1,aligned,align_spec);
    for(size_t i=0;i<aligned.size();i++) {
      cout << aligned[i] << endl;
    }

  } else if (type=="tensor_grid") {

    size_t rank=tensor_grid_obj.get_rank();

    size_t max_size=0;
    for(size_t k=0;k<rank;k++) {
      if (tensor_grid_obj.get_size(k)>max_size) {
	max_size=tensor_grid_obj.get_size(k);
      }
    }

    vector<vector<string> > string_mat(rank+1);
    vector<int> align_spec(rank+1);

    // The first column which enumerates the grid points
    align_spec[0]=columnify::align_left;
    string_mat[0].resize(max_size+1);
    for(size_t ell=0;ell<max_size;ell++) {
      string_mat[0][ell+1]=o2scl::szttos(ell)+".";
    }

    // The first row which labels the grids
    for(size_t k=0;k<rank;k++) {
      string_mat[k+1].resize(max_size+1);
      align_spec[k+1]=columnify::align_right;
      string_mat[k+1][0]=((string)"Grid ")+o2scl::szttos(k);
    }

    // Now the grid data
    for(size_t ell=0;ell<max_size;ell++) {
      for(size_t k=0;k<rank;k++) {
	if (ell<tensor_grid_obj.get_size(k)) {
	  string_mat[k+1][ell+1]=
	    o2scl::dtos(tensor_grid_obj.get_grid(k,ell),precision);
	}
      }
    }

    columnify col;
    vector<string> aligned(max_size+1);
    col.align(string_mat,rank+1,max_size+1,aligned,align_spec);
    for(size_t i=0;i<aligned.size();i++) {
      cout << aligned[i] << endl;
    }

  } else {

    cout << "Not implemented for type " << type << endl;
  }
  
  return 0;
}

int acol_manager::comm_get_row(std::vector<std::string> &sv,
                               bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  if (table_obj.get_nlines()==0 || table_obj.get_nlines()==0) {
    cerr << "No table or empty table in get-row." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter row number to get",
			i1,"get-row",itive_com);
  if (ret!=0) return ret;

  int ix=((int)(o2scl::function_to_double(i1)));
  
  // If negative, view as distance from end
  if (ix<0) ix+=table_obj.get_nlines();

  if (ix<0) {
    cerr << "Requested negative row in 'get-row'." << endl;
    return exc_efailed;
  }
  if (ix>((int)table_obj.get_nlines())-1) {
    cerr << "Requested row beyond end of table in get-row." << endl;
    return exc_efailed;
  }
  
  //--------------------------------------------------------------------
  // Compute number of screen columns

  int ncols_loc;
  if (ncols<=0) {
    int srow, scol;
    int iret=get_screen_size_ioctl(srow,scol);
    if (scol>10 || iret!=0) ncols_loc=80;
    else ncols_loc=scol;
  } else {
    ncols_loc=ncols;
  }
  cout << "Number of columns: " << ncols_loc << endl;

  //--------------------------------------------------------------------
  // Temporary storage strings for names and data

  vector<string> row_names, row_data;

  //--------------------------------------------------------------------
  // Process and/or output names

  if (names_out==true) {

    if (pretty) {

      size_t running_width=0;
      ostringstream str;
      // Clear ostringstream with str.str(""); and str.clear();
      str.setf(ios::scientific);
      str.precision(precision);
      
      for(size_t i=0;i<table_obj.get_ncolumns();i++) {

	// Count for space between columns and sign
	size_t this_col=2;
	// Count column name
	this_col+=table_obj.get_column_name(i).size();
	// Count extra spaces to format number
	int num_spaces=precision+6-
          ((int)(table_obj.get_column_name(i).size()));
	if (num_spaces>0) this_col+=num_spaces;
	// See if there will be space
	if (running_width>0 && ((int)(running_width+this_col))>=ncols_loc) {
	  row_names.push_back(str.str());
	  str.str("");
	  str.clear();
	  str.setf(ios::scientific);
	  str.precision(precision);
	  running_width=0;
	}
	// Output this column name
	str << ' ' << table_obj.get_column_name(i) << ' ';
	for(int j=0;j<num_spaces;j++) {
	  str << ' ';
	}
	running_width+=this_col;
      }
      row_names.push_back(str.str());
      str.str("");
      str.clear();
      
    } else {
      
      cout.precision(precision);
  
      for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	cout << table_obj.get_column_name(i) << ' ';
      }
      cout << endl;

    }
  }
  
  //--------------------------------------------------------------------
  // Process and/or output data
  
  if (pretty) {
    
    size_t running_width=0;
    ostringstream str;
    str.setf(ios::scientific);
    str.precision(precision);
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      
      // Count space for number
      size_t this_col=precision+8;
      // Count extra spaces if necessary
      int num_spaces=((int)(table_obj.get_column_name(i).size())-
                      precision-6);
      if (num_spaces>0) this_col+=num_spaces;
      // See if there will be space
      if (running_width>0 && ((int)(running_width+this_col))>=ncols_loc) {
	row_data.push_back(str.str());
	str.str("");
	str.clear();
	str.setf(ios::scientific);
	str.precision(precision);
	running_width=0;
      }
      // Output the data
      if (table_obj.get(i,ix)>=0.0) {
	str << ' ' << table_obj.get(i,ix) << ' ';
      } else {
	str << table_obj.get(i,ix) << ' ';
      }
      for(int j=0;j<num_spaces;j++) {
	str << ' ';
      }
      running_width+=this_col;
    }
    row_data.push_back(str.str());
    str.str("");
    str.clear();
    
    //--------------------------------------------------------------------
    // Now output both names and data to cout

    if (row_names.size()!=row_data.size()) {
      O2SCL_ERR("Names and data size don't match in get-row.",
		exc_esanity);
    }
    for(size_t k=0;k<row_names.size();k++) {
      cout << row_names[k] << endl;
      cout << row_data[k] << endl;
    }
    
  } else {
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      cout << table_obj.get(i,ix) << ' ';
    }
    cout << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_get_unit(std::vector<std::string> &sv,
                                bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to get units for." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to get units of");
  int ret=get_input(sv,pr,in,"get-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  cout << "Units of column " << in[0] << " are: " 
       << table_obj.get_unit(in[0]) << endl;

  return 0;
}

int acol_manager::comm_help(std::vector<std::string> &sv, bool itive_com) {

  terminal ter;

  int ncols_loc;
  if (ncols<=0) {
    int srow, scol;
    int iret=get_screen_size_ioctl(srow,scol);
    if (scol<10 || iret!=0) ncols_loc=80;
    else ncols_loc=scol;
  } else {
    ncols_loc=ncols;
  }
  
  // Create a line for separating help text sections
  string line=ter.hrule(ncols_loc-2);

  // Handle the 'help type command' case for type-specific commands
  if (sv.size()==3) {
    string temp_type=sv[1];
    string cur_type=type;

    bool valid_type=false;
    for(size_t j=0;j<type_list.size() && valid_type==false;j++) {
      if (temp_type==type_list[j]) {
	valid_type=true;
      }
    }

    if (valid_type) {
      
      command_del(type);
      command_add(temp_type);
      
      if (cl->is_valid_option(sv[2])) {

	string cmd=sv[2];
	
	std::vector<std::string>::iterator it=sv.begin();
	it++;
	sv.erase(it);
	
	cout << "Help for command " << command_color
	     << cmd << default_color << " given object of type "
             << type_color << temp_type << default_color
             << ".\n" << endl;
        
	int ret=cl->comm_option_help(sv,itive_com);
	
	command_del(temp_type);
	command_add(cur_type);
	
	return ret;

      } else {
	
	command_del(temp_type);
	command_add(cur_type);
	
	cout << "Command " << command_color
	     << sv[2] << default_color << " not valid for type "
             << type_color << sv[1] << default_color 
	     << ".\n" << endl;
	
	vector<string> sv2={"help","help"};
	cl->comm_option_help(sv2,itive_com);
	return 1;
	
      }
      
    } else {
      cout << "Type " << type_color
	   << sv[1] << default_color << " not found.\n" << endl;

      vector<string> sv2={"help","help"};
      cl->comm_option_help(sv2,itive_com);
      return 2;
    }
    
  }

  // Handle the special case 'help types'
  if (sv.size()==2 && sv[1]=="types") {
    cout << "Documentation for help topic: " << help_color
	 << "types" << default_color << endl;
    cout << line << "\n" << endl;

    string str="The O2scl types which can be handled by "+cl->cmd_name;
    str+=" are ";
    for(size_t i=0;i<type_list.size()-1;i++) {
      str+=type_color+type_list[i]+default_color+", ";
    }
    str+="and "+type_color+
      type_list[type_list.size()-1]+default_color+'.';
    str+="\n\n";
    str+="For more information on some of the types:\n";
    str+=type_color+"hist"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/hist.html"+
      default_color+"\n";
    str+=type_color+"hist_2d"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/hist_2d.html"+
      default_color+"\n";
    str+=type_color+"table"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/table.html"+
      default_color+"\n";
    str+=type_color+"table3d"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/table3d.html"+
      default_color+"\n";
    str+=type_color+"tensor"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/tensor.html"+
      default_color+"\n";
    str+=type_color+"tensor_grid"+default_color+
      ": "+url_color+
      "https://awsteiner.org/code/o2scl/html/class/tensor_grid.html"+
      default_color+"\n";

    std::vector<std::string> svx;
    o2scl::rewrap_ignore_vt100(str,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
    cout << endl;
      
    return 0;
  }
  
  // Handle the special case 'help functions'
  if (sv.size()==2 && (sv[1]=="functions" || sv[1]=="functions")) {
    cout << "Documentation for help topic: " << help_color
	 << "functions" << default_color << endl;
    cout << line << "\n" << endl;

    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="functions") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            if (help_doc_strings[j][kk].length()>0) {
              help+="\n\n"+help_doc_strings[j][kk];
            }
          }
        }
      }
    }

    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help index-spec'
  if (sv.size()==2 && (sv[1]=="index-spec" || sv[1]=="index_spec")) {
    cout << "Documentation for help topic: " << help_color
	 << "index-spec" << default_color << endl;
    cout << line << "\n" << endl;

    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="index_spec") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            if (help_doc_strings[j][kk].length()>0) {
              help+="\n\n"+help_doc_strings[j][kk];
            }
          }
        }
      }
    }

    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help value-spec'
  if (sv.size()==2 && (sv[1]=="value-spec" || sv[1]=="value_spec")) {
    cout << "Documentation for help topic: " << help_color
	 << "value-spec" << default_color << endl;
    cout << line << "\n" << endl;

    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="value_spec") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            if (help_doc_strings[j][kk].length()>0) {
              help+="\n\n"+help_doc_strings[j][kk];
            }
          }
        }
      }
    }

    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help vector-spec'
  if (sv.size()==2 && (sv[1]=="vector-spec" || sv[1]=="vector_spec")) {
    cout << "Documentation for help topic: " << help_color
	 << "vector-spec" << default_color << endl;
    cout << line << "\n" << endl;

    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="vector_spec") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            help+="\n\n"+help_doc_strings[j][kk];
          }
        }
      }
    }
    
    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help strings-spec'
  if (sv.size()==2 && (sv[1]=="strings-spec" || sv[1]=="strings_spec")) {
    cout << "Documentation for help topic: " << help_color
	 << "strings-spec" << default_color << endl;
    cout << line << "\n" << endl;
    
    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="strings_spec") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            if (help_doc_strings[j][kk].length()>0) {
              help+="\n\n"+help_doc_strings[j][kk];
            }
          }
        }
      }
    }
    
    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help mult-vector-spec'
  if (sv.size()==2 && (sv[1]=="mult-vector-spec" ||
		       sv[1]=="mult_vector_spec")) {
    cout << "Documentation for help topic: " << help_color
	 << "mult-vector-spec" << default_color << endl;
    cout << line << "\n" << endl;
    
    std::string help;
    for(size_t j=0;j<help_doc_strings.size();j++) {
      if (help_doc_strings[j][0]=="mult_vector_spec") {
        for(size_t kk=1;kk<help_doc_strings[j].size();kk++) {
          if (kk==1) {
            help=help_doc_strings[j][kk];
          } else {
            if (help_doc_strings[j][kk].length()>0) {
              help+="\n\n"+help_doc_strings[j][kk];
            }
          }
        }
      }
    }

    color_replacements(help);
    std::vector<std::string> svx;
    o2scl::rewrap_keep_endlines(help,svx,ncols_loc-1);
    for(size_t i=0;i<svx.size();i++) {
      cout << svx[i] << endl;
    }
      
    return 0;
  }

  // Handle the case 'help <type>' where <type> is an acol o2scl type
  for(size_t i=0;i<type_list.size();i++) {
    if (sv.size()>=2 && sv[1]==type_list[i]) {
      cout << "Documentation for type: " << type_color
	   << type_list[i] << default_color << endl;
      cout << line << "\n" << endl;
      std::string str="Objects of type "+type_color+
	type_list[i]+default_color+
	" can be read, written, or modified by "+cl->cmd_name+".";
      std::vector<std::string> svx;
      o2scl::rewrap_keep_endlines(str,svx,ncols_loc-1);
      for(size_t j=0;j<svx.size();j++) {
	cout << svx[j] << endl;
      }
      cout << endl;
      svx.clear();
      std::map<std::string,std::vector<std::string> >::iterator it;
      for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
	if (it->first==type_list[i]) {
	  std::vector<std::string> &clist=it->second;
	  if (clist.size()>1) {
	    str=((string)"The type-specific commands for objects ")+
	      "of type "+type_list[i]+" are: ";
	    for(size_t j=0;j<clist.size()-1;j++) {
	      str+=command_color+clist[j]+default_color+", ";
	    }
	    str+=" and "+command_color+
	      clist[clist.size()-1]+default_color+".";
	    o2scl::rewrap_ignore_vt100(str,svx,ncols_loc-1);
	    for(size_t j=0;j<svx.size();j++) {
	      cout << svx[j] << endl;
	    }
	  } else {
	    cout << "The only type-specific command for objects of type "
		 << type_color
		 << type_list[i] << default_color
		 << " is " << command_color
		 << clist[0] << default_color << "." << endl;
	  }
	}
      }
      
      return 0;
    }
  }

  // Handle the case 'help <command>' where <command> is a type-specific
  // command
  if (sv.size()>=2 && cl->is_valid_option(sv[1])==false) {

    bool found=false;

    std::map<std::string,std::vector<std::string> >::iterator it;
    for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
      std::vector<std::string> &clist=it->second;
      for(size_t j=0;j<clist.size();j++) {
	if (clist[j]==sv[1]) {

	  if (found==false) {
	    cout << "Command " << command_color << sv[1]
                 << default_color << " is a type specific "
		 << "command. Below are "
		 << "the various\ndescriptions of its operation with "
		 << "the relevant types." << endl;
	    if (sv[1]=="function") {
	      cout << "Type \"" << cl->cmd_name << " -help functions\" "
		   << "for help on specifying functions as\ncommand "
		   << "arguments." << endl;
	    }
	    cout << endl;
	    found=true;
	  }
	  
	  cout << line << endl;
	  
	  string cur_type=type;
	  
	  vector<string> sv2;
	  
	  sv2.push_back("help");
	  sv2.push_back(sv[1]);
	  
	  command_del(type);
	  command_add(it->first);
	  type=it->first;

	  cout << "Type " << type_color << it->first
	       << default_color << ":" << endl;
	  int ret=cl->comm_option_help(sv2,itive_com);
	  cout << endl;

	  command_del(type);
	  command_add(cur_type);
	  type=cur_type;
	}
      }
    }

    if (found) return 0;
  }

  // ----------------------------------------------------------------
  // Handle the generic help case using cli::comm_option_help(). We
  // have to define the additional help text here because we need to
  // know if the output is being redirected to decide whether or not
  // to print colors and lines.
  
  string stemp;
  string dsc="\n"+line+"\nNotes:\n\n";
  vector<std::string> sv2;
  
  stemp="1. Help for commands which apply to the current object ";
  stemp+="may be obtained with '-help ";
  stemp+="[command]'. Help for type-specific commands can be obtained ";
  stemp+="by '-help [type] [command]'. A list of commands for each type ";
  stemp+="can be obtained with 'commands [type]', or for all commands ";
  stemp+="use '-commands all'. Required arguments ";
  stemp+="are surrounded by ";
  stemp+="<>'s and optional arguments are surrounded by []'s.\n";
  rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
  dsc+=sv2[0]+"\n";
  for(size_t j=1;j<sv2.size();j++) {
    dsc+="   "+sv2[j]+"\n";
  }

  stemp="2. Options may also be specified in the environment variable ";
  stemp+=env_var_name+".\n";
  rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
  dsc+=sv2[0]+"\n";
  for(size_t j=1;j<sv2.size();j++) {
    dsc+="   "+sv2[j]+"\n";
  }

  stemp="3. Long options may be preceeded by two dashes.\n";
  rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
  dsc+=sv2[0]+"\n";
  for(size_t j=1;j<sv2.size();j++) {
    dsc+="   "+sv2[j]+"\n";
  }

  stemp="4. In order to avoid confusion between arguments and functions, ";
  stemp+="use parenthesis and quotes, i.e. '(-x*2)' instead of -x*2.\n";
  rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
  dsc+=sv2[0]+"\n";
  for(size_t j=1;j<sv2.size();j++) {
    dsc+="   "+sv2[j]+"\n";
  }

  stemp="5. Also, do not use a unary minus next to a binary operator, ";
  stemp+="i.e. use 'a>(-1)' instead of 'a>-1'.\n\n";
  rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
  dsc+=sv2[0]+"\n";
  for(size_t j=1;j<sv2.size();j++) {
    dsc+="   "+sv2[j]+"\n";
  }

  if (ter.is_redirected()==false) {
    
    stemp="6. Types are denoted as "+type_color+"char";
    stemp+=default_color+", commands as "+command_color;
    stemp+="function"+default_color+",\n   get/set parameters as ";
    stemp+=param_color+"verbose"+default_color;
    stemp+=", help topics as\n   ";
    stemp+=help_color+"functions"+default_color+", command-line ";
    stemp+="scripts as "+exec_color+"acol -help"+default_color;
    stemp+=", and URLs as "+url_color+"https://arxiv.org";
    stemp+=default_color+".\n";
    
    rewrap_ignore_vt100(stemp,sv2,ncols_loc-4);
    
    dsc+=sv2[0]+"\n";
    for(size_t j=1;j<sv2.size();j++) {
      dsc+="   "+sv2[j]+"\n";
    }
    
  }
  
  dsc+="\n"+line+"\n";

  cl->addl_help_cmd=dsc;
  cl->addl_help_cli=dsc;
  
  int ret=cl->comm_option_help(sv,itive_com);

  if (sv.size()>=2 && cl->is_valid_option(sv[1])==false &&
      cl->is_parameter(sv[1])==false) {
    cout << "\n" << line << endl;
  }
  
  if (sv.size()<2 || (cl->is_valid_option(sv[1])==false &&
		      cl->is_parameter(sv[1])==false)) {
    
    cout << "List of additional help topics (e.g. \"acol -help [topic]\"): ";
    cout << help_color << "functions" << default_color
	 << "," << endl;
    cout << help_color << "index-spec"
	 << default_color << ", ";
    cout << help_color << "mult-vector-spec"
	 << default_color << ", ";
    cout << help_color << "strings-spec" << default_color
	 << ", ";
    cout << help_color << "types" << default_color
	 << ", ";
    cout << help_color << "value-spec" << default_color
	 << ", and\n";
    cout << help_color << "vector-spec" << default_color
	 << ".\n" << endl;

    cout << line << endl;
  
#ifndef O2SCL_UBUNTU_PKG
    cout << ((string)"Compiled at ")+((string)__TIME__)+" on "+
      ((string)__DATE__)+" for "+exec_color+"O₂scl"+
      default_color+", version "+((string)VERSION)+".\n" << endl;
#else
    cout << ((string)"Compiled for ")+exec_color+"O₂scl"+
      default_color+", version "+((string)VERSION)+".\n" << endl;
#endif

  }
  
  return ret;
}

int acol_manager::comm_h5_copy(std::vector<std::string> &sv, 
			       bool itive_com) {

  vector<string> in, pr;
  
  pr.push_back("Source file");
  pr.push_back("Destination file");
  int ret=get_input(sv,pr,in,"h5-copy",itive_com);
  if (ret!=0) return ret;

  if (in[0]==in[1]) {
    cerr << "Command 'h5-copy' will not copy a file onto itself." << endl;
    return 2;
  }

  // Use hdf_file to open the file
  hdf_file hf, hf2;
  int hfret=hf.open(in[0].c_str(),false,false);
  if (hfret!=0) {
    cerr << "Failed to read file named " << in[0].c_str() << endl;
    return exc_efailed;
  }
  hf2.open_or_create(in[1].c_str());

  hf.copy(verbose,hf2);
  
  return 0;
}

int acol_manager::comm_index(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to add line numbers to." << endl;
    return exc_efailed;
  }

  std::string i1="N";
  if (sv.size()>1) i1=sv[1];
  int j;
  for (j=0;j<10 && table_obj.is_column(i1);j++) {
    i1+="_";
  }
  if (j==10) {
    cerr << "Could not find unique column name in 'index'." << endl;
    return 1;
  }
  table_obj.new_column(i1);
  if (verbose>0) {
    cout << "Created new column named " << i1 << "." << endl;
  }
  for(size_t i=0;i<table_obj.get_nlines();i++) table_obj.set(i1,i,((double)i));

  return 0;
}

int acol_manager::comm_insert(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table" && type!="table3d") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (type=="table3d") {

    std::string in[4], pr[4]=
      {"Enter filename of external table (or blank to stop): ",
       "Enter name of table in file (or blank for first table): ",
       "Enter slice in external table (or blank to stop): ",
       "Enter name of new slice in present table (blank to keep old name): "};

    if (sv.size()>3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
      if (sv.size()>4) {
	in[3]=sv[4];
      } else {
	in[3]="";
      }
    } else {
      if (itive_com) {
	for(size_t is=0;is<4;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (is!=3 && is!=1 && in[is].length()==0) {
	    cout << "Command 'insert' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to 'insert'" << endl;
	return exc_efailed;
      }
    }

    if (true || verbose>2) {
      cout << "Read table3d named " << in[1] << " from file "
           << in[0] << endl;
      cout << "old slice, new slice: " << in[2] << " " << in[3] << endl;
    }

    // Read table from file
    hdf_file hf;
    table3d tmp;
    int hfret=hf.open(in[0],false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << in[0] << endl;
      return exc_efailed;
    }
    
    std::string tmp_name;
    if (in[1].length()>0) tmp_name=in[1];
    hdf_input(hf,tmp,tmp_name);
    hf.close();

    if (verbose>2) {
      cout << "Original table3d object numx,numy: "
           << table3d_obj.get_nx() << " "
           << table3d_obj.get_ny() << endl;
      cout << "New table3d object numx,numy: "
           << tmp.get_nx() << " "
           << tmp.get_ny() << endl;
    }

    table3d_obj.add_slice_from_table(tmp,in[2],in[3],3);

    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into." << endl;
    return exc_efailed;
  }

  std::string in[6], pr[6]=
    {"Enter filename of external table (or blank to stop): ",
     "Enter name of table in file (or blank for first table): ",
     "Enter index column in external table (or blank to stop): ",
     "Enter data column in external table (or blank to stop): ",
     "Enter index column in present table (or blank to stop): ",
     "Enter name of new column in present table (or blank to keep old name): "};
  if (sv.size()>5) {
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
    in[3]=sv[4];
    in[4]=sv[5];
    if (sv.size()>6) {
      in[5]=sv[6];
    } else {
      in[5]="";
    }
  } else {
    if (itive_com) {
      for(size_t is=0;is<6;is++) {
	in[is]=cl->cli_gets(pr[is].c_str());
	if (is!=5 && is!=1 && in[is].length()==0) {
	  cout << "Command 'insert' cancelled." << endl;
	  return 0;
	}
      }
    } else {
      cerr << "Not enough arguments to 'insert'" << endl;
      return exc_efailed;
    }
  }

  if (true || verbose>2) {
    cout << "Read table named " << in[1] << " from file " << in[0] << endl;
    cout << "oldx,oldy,newx,newy: " << in[2] << " " << in[3] << " " 
	 << in[4] << " " << in[5] << endl;
    cout << endl;
  }

  // Read table from file
  hdf_file hf;
  table_units<> tmp;
  int hfret=hf.open(in[0],false,false);
  if (hfret!=0) {
    cerr << "Failed to read file named " << in[0] << endl;
    return exc_efailed;
  }
  std::string tmp_name;
  if (in[1].length()>0) tmp_name=in[1];
  hdf_input(hf,tmp,tmp_name);
  hf.close();

  table_obj.add_col_from_table(tmp,in[2],in[3],in[4],in[5]);

  return 0;

}

int acol_manager::comm_insert_full(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into in command 'insert-full'."
	 << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
    
  if (sv.size()==2) {
    
    cout << "No index columns or table name specified. "
	 << "Trying without interpolation." << endl;
    in.push_back(sv[1]);
    
  } else if (sv.size()==3) {
    
    cout << "No index columns specified. Trying without interpolation."
	 << endl;
    in.push_back(sv[1]);
    in.push_back(sv[2]);
    
  } else if (sv.size()==4) {
    
    cout << "Presuming index column in source and destination are "
	 << "the same." << endl;
    
    in.push_back(sv[1]);
    in.push_back(sv[2]);
    in.push_back(sv[3]);
    
  } else {
    
    pr.push_back("Enter filename of external table");
    pr.push_back("Enter name of table in file");
    pr.push_back("Enter index column in current table");
    pr.push_back("Enter index column in external table");
    int ret=get_input(sv,pr,in,"insert-full",itive_com);
    if (ret!=0) return ret;

  }

  // Read table from file
  hdf_file hf;
  table_units<> tmp;
  int hfret=hf.open(in[0],false,false);
  if (hfret!=0) {
    cerr << "Failed to read file named " << in[0] << endl;
    return exc_efailed;
  }

  std::string tmp_name;
  if (in.size()>=2) {
    if (in[1].length()>0) tmp_name=in[1];
  }
  hdf_input(hf,tmp,tmp_name);
  hf.close();

  if (in.size()==1 || in.size()==2) {

    // No index columns were specified, so see if we can avoid
    // interpolation

    if (table_obj.get_nlines()!=tmp.get_nlines()) {
      cout << "Number of lines in current table "
	   << table_obj.get_nlines()
	   << "\n\tdoes not match number of lines in external table "
	   << tmp.get_nlines() << endl;
      return 3;
    }

    for (size_t j=0;j<tmp.get_ncolumns();j++) {
      string ty=tmp.get_column_name(j);
      if (!table_obj.is_column(ty)) table_obj.new_column(ty);
      for(size_t k=0;k<table_obj.get_nlines();k++) {
	table_obj.set(ty,k,tmp.get(ty,k));
      }
    }
    
  } else {
  
    for (size_t j=0;j<tmp.get_ncolumns();j++) {
      string ty=tmp.get_column_name(j);
      if (ty!=in[2]) {
	if (in.size()==3) {
	  table_obj.add_col_from_table(tmp,in[2],ty);
	} else {
	  table_obj.add_col_from_table(tmp,in[2],ty,in[3]);
	}
      }
    }
  }

  return 0;
}

int acol_manager::comm_integ(std::vector<std::string> &sv, bool itive_com) {

  ff.set_sig_figs(precision+1);
  
  if (type=="table") {

    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to integrate." << endl;
      return exc_efailed;
    }
    vector<string> pr, in;
    pr.push_back("Enter 'x' column");
    pr.push_back("Enter 'y' column");
    pr.push_back("Enter name of new column");
    int ret=get_input(sv,pr,in,"integ",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    if (table_obj.is_column(in[1])==false) {
      cerr << "Could not find column named '" << in[1] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.integ(in[0],in[1],in[2]);
  
  } else if (type=="hist") {

    double d=hist_obj.integ_wgts();
    cout << "Integral over all "
         << hist_obj.size() << " bins in the histogram is "
         << d << " (" << ff.convert(d) << ")." << endl;
  
  } else if (type=="hist_2d") {

    double d=hist_2d_obj.integ_wgts();
    cout << "Integral over all "
         << hist_2d_obj.size_x() << " by "
         << hist_2d_obj.size_y() << " bins in the histogram is "
         << d << " (" << ff.convert(d) << ")." << endl;
  
  } else {
    cerr << "Command 'integ' not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }


  return 0;
}

int acol_manager::comm_internal(std::vector<std::string> &sv, bool itive_com) {
  
  std::string i1;
  int ret=get_input_one(sv,"Enter filename",i1,"internal",itive_com);
  if (ret!=0) return ret;
  
  if (verbose>2) {
    cout << "Opening O2scl file: " << i1 << endl;
  }
  
  hdf_file hf;
  hf.compr_type=compress;
  hf.open_or_create(i1);
  
  if (type=="int") {

    hf.seti(obj_name,int_obj);
    
  } else if (type=="double") {

    hf.setd(obj_name,double_obj);
    
  } else if (type=="char") {

    hf.setc(obj_name,char_obj);
    
  } else if (type=="string") {

    hf.sets(obj_name,string_obj);
    
  } else if (type=="size_t") {

    hf.set_szt(obj_name,size_t_obj);
    
  } else if (type=="double[]") {

    hf.setd_vec(obj_name,doublev_obj);
    
  } else if (type=="tensor") {

    hf.setd_ten(obj_name,tensor_obj);
    
  } else if (type=="int[]") {

    hf.seti_vec(obj_name,intv_obj);
    
  } else if (type=="size_t[]") {

    hf.set_szt_vec(obj_name,size_tv_obj);
    
  } else if (type=="string[]") {

    hf.sets_vec_copy(obj_name,stringv_obj);
    
  } else if (type=="vec_vec_string") {

    hf.sets_vec_vec_copy(obj_name,vvstring_obj);
    
  } else if (type=="vec_vec_double") {

    hf.setd_vec_vec_copy(obj_name,vvdouble_obj);
    
  } else if (type=="table3d") {
    
    hdf_output(hf,((const table3d &)(table3d_obj)),obj_name);
    
  } else if (type=="tensor_grid") {
    
    hdf_output(hf,tensor_grid_obj,obj_name);
    
  } else if (type=="table") {
    
    hdf_output(hf,table_obj,obj_name);

  } else if (type=="prob_dens_mdim_amr") {
    
    hdf_output(hf,pdma_obj,obj_name);

  } else if (type=="prob_dens_mdim_gaussian") {
    
    hdf_output(hf,pdmg_obj,obj_name);

  } else if (type=="prob_dens_mdim_gmm") {
    
    hdf_output(hf,pgmm_obj,obj_name);

  } else if (type=="hist") {

    hdf_output(hf,hist_obj,obj_name);
    
  } else if (type=="hist_2d") {

    hdf_output(hf,((const hist_2d &)(hist_2d_obj)),obj_name);
    
  } else if (type=="vector<contour_line>") {

    hdf_output(hf,cont_obj,obj_name);
    
  } else if (type=="uniform_grid<double>") {

    hdf_output(hf,ug_obj,obj_name);
    
  } else if (type=="") {

    cerr << "There is no current object, so the 'internal' command "
         << "cannot write to a file." << endl;

  } else {
    
    O2SCL_ERR("Canno output object of this type.",
              o2scl::exc_esanity);
    
  }

  hf.close();
    
  return 0;
}

int acol_manager::comm_interp(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    // --------------------------------------------------------------
    // 3d table interpolation
    
    if (type!="table3d") {
      cerr << "No table to interpolate into." << endl;
      return exc_efailed;
    }

    std::string in[3], pr[3]=
      {"Enter slice name (or blank to stop): ",
       "Enter x value (or blank to stop): ",
       "Enter y value (or blank to stop): "};
    if (sv.size()>=3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
    } else {
      if (itive_com) {
	for(size_t is=0;is<3;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (in[is].length()==0) {
	    cout << "Command 'interp' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to command 'interp'" << endl;
	return exc_efailed;
      }
    }
  
    double ret=table3d_obj.interp(function_to_double(in[1]),
				  function_to_double(in[2]),in[0]);
    if (err_hnd->get_errno()!=0) {
      cerr << "Interpolation failed." << endl;
      return exc_efailed;
    } else {
      cout << "Interpolation result: " << ret << endl;
    }

    return 0;

  } else if (type=="table") {

    // --------------------------------------------------------------
    // 2d table interpolation
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to interpolate into." << endl;
      return exc_efailed;
    }
    
    std::string in[3], pr[3]=
      {"Enter column name of independent variable (or blank to stop): ",
       "Enter value of independent variable (or blank to stop): ",
       "Enter column name of dependent variable (or blank to stop): "};
    if (sv.size()>=3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
    } else {
      if (itive_com) {
	for(size_t is=0;is<3;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (in[is].length()==0) {
	    cout << "Command 'interp' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to 'interp'" << endl;
	return exc_efailed;
      }
    }
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    if (table_obj.is_column(in[2])==false) {
      cerr << "Could not find column named '" << in[2] << "'." << endl;
      return exc_efailed;
    }

    double ret=table_obj.interp(in[0],function_to_double(in[1]),in[2]);
    
    if (err_hnd->get_errno()!=0) {
      cerr << "Interpolation failed." << endl;
      return exc_efailed;
    } else {
      cout << "Interpolation result: " << ret << endl;
    }

    // --------------------------------------------------------------
    
  } else if (type=="double[]") {

    size_t n=doublev_obj.size();
    std::vector<double> index(n);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp_vec<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    it.set(n,index,doublev_obj);
    cout << "Interpolation result: " << it.eval(x) << endl;

  } else if (type=="tensor_grid") {
    
    vector<string> in, pr;
    for(size_t i=0;i<tensor_grid_obj.get_rank();i++) {
      pr.push_back(((std::string)"Value for index ")+
		   o2scl::szttos(i));
    }
    int ret=get_input(sv,pr,in,"interp",itive_com);
    if (ret!=0) return ret;

    vector<double> vals;
    for(size_t i=0;i<tensor_grid_obj.get_rank();i++) {
      vals.push_back(o2scl::stod(in[i]));
    }
    
    double res=tensor_grid_obj.interp_linear(vals);
    cout << "Interpolation result: " << res << endl;

  } else if (type=="int[]") {

    size_t n=intv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(intv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp_vec<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    it.set(n,index,value);
    cout << "Interpolation result: " << it.eval(x) << endl;
	 

  } else if (type=="size_t[]") {

    size_t n=size_tv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(size_tv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp_vec<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    it.set(n,index,value);
    cout << "Interpolation result: " << it.eval(x) << endl;

  } else {
    cout << "Not implemented for type " << type << endl;
    return 1;
  }    
  
  return 0;
}

int acol_manager::comm_interp_table3d(std::vector<std::string> &sv,
                                      bool itive_com) {

  if (type=="table") {
    
    // --------------------------------------------------------------
    // 2d table interpolation
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to interpolate into." << endl;
      return exc_efailed;
    }
    
    vector<std::string> in;
    
    if (sv.size()<8) {

      in.resize(7);
      vector<std::string> pr=
        {"Enter method (or blank to stop): ",
         "Enter options (or blank to stop): ",
         "Enter x column (or blank to stop): ",
         "Enter y column (or blank to stop): ",
         "Enter x grid, \"auto\", (or blank to stop): ",
         "Enter y grid, \"auto\", (or blank to stop): ",
         "Enter the z column (or blank to stop): "};
      int ret=get_input(sv,pr,in,"interp-table3d",itive_com);
      if (ret!=0) return ret;
      
    } else {
      for(size_t j=1;j<sv.size();j++) {
        in.push_back(sv[j]);
      }
    }

    // Determine method and options
    std::string method=in[0];
    if (verbose>1) {
      std::cout << "In acol command 'interp-table3d':" << endl;
      std::cout << "  Using method: " << method << endl;
    }
    std::string options=in[1];

    // Set up x and y-grids
    vector<double> gx, gy;
    
    double n_bins_d=pow(table_obj.get_nlines(),1.0/3.0);
    size_t n_bins;
    if (n_bins_d<5.0) {
      n_bins=5;
    } else {
      n_bins=(size_t)n_bins_d;
    }
    
    if (in[4]=="auto") {
      double min_x=vector_min_value<vector<double>,double>
        (table_obj[in[2]]);
      double max_x=vector_max_value<vector<double>,double>
        (table_obj[in[2]]);
      uniform_grid_end<double> ugx(min_x,max_x,n_bins);
      ugx.vector(gx);
    } else {
      vector_spec(in[4],gx,false,2);
      
    }
    if (in[5]=="auto") {
      double min_y=vector_min_value<vector<double>,double>
        (table_obj[in[3]]);
      double max_y=vector_max_value<vector<double>,double>
        (table_obj[in[3]]);
      uniform_grid_end<double> ugy(min_y,max_y,n_bins);
      ugy.vector(gy);
    } else {
      vector_spec(in[5],gy,false,2);
    }

    // Clear table3d, set up grids and slices
    table3d_obj.clear();
    
    table3d_obj.set_xy(in[2],gx.size(),gx,in[3],gy.size(),gy);
      
    for(size_t k=6;k<in.size();k++) {
      table3d_obj.new_slice(in[k]);
    }

    // Perform the final interpolation
    if (method=="idw") {

      vector<string> col_list={in[2],in[3]};
      for (size_t j=6;j<in.size();j++) {
        col_list.push_back(in[j]);
      }
      cout << "  Column list: ";
      vector_out(cout,col_list,true);
      
      const_matrix_view_table_transpose<> cmvt(table_obj,col_list);
      interpm_idw<ubvector,const_matrix_view_table_transpose<>,
                  const_matrix_view_table_transpose<>> imi;

      imi.set_data(2,col_list.size()-2,table_obj.get_nlines(),cmvt,
                   cmvt);
      
      for(size_t i=0;i<table3d_obj.get_nx();i++) {
        for(size_t j=0;j<table3d_obj.get_ny();j++) {
          vector<double> px={table3d_obj.get_grid_x(i),
            table3d_obj.get_grid_y(j)};
          vector<double> py(in.size()-6);
          imi.eval_tl(px,py);
          for(size_t k=6;k<in.size();k++) {
            table3d_obj.set(i,j,k-6,py[k-6]);
          }
        }
      }
      
      command_del(type);
      clear_obj();
      command_add("table3d");
      type="table3d";
      
    } else if (method=="gp") {
      
      vector<string> col_list_x={in[2],in[3]}, col_list_y;
      for (size_t j=6;j<in.size();j++) {
        col_list_y.push_back(in[j]);
      }
      cout << "  Column list for x: ";
      vector_out(cout,col_list_x,true);
      cout << "  Column list for y: ";
      vector_out(cout,col_list_y,true);

      typedef const const_matrix_row_gen<o2scl::const_matrix_view_table<>>
        mat_x_row_t;
      
      vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
      vmfrn.resize(col_list_y.size());
      
      const_matrix_view_table<> mvt_x(table_obj,col_list_x);
      matrix_view_table<> mvt_y(table_obj,col_list_y);
      
      interpm_krige_optim<> iko;

      vector<vector<vector<double>>> param_lists;
      
      for(size_t kk=0;kk<col_list_y.size();kk++) {
        std::shared_ptr<mcovar_funct_rbf_noise<
          ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                      mat_x_row_t>);
        vmfrn[kk]=mfrn;
        for(size_t i=0;i<col_list_y.size();i++) {
          mfrn->len.resize(2);
        }
      }
      
      for(size_t i=0;i<col_list_y.size();i++) {
        std::vector<std::vector<double>> ptemp;
        
        for(size_t j=0;j<2;j++) {
          std::vector<double> diffs;
          o2scl::vector_diffs<std::vector<double>,std::vector<double>>
            (table_obj[in[j+2]],diffs);
          double min=vector_min_value<std::vector<double>,double>
            (diffs.size(),diffs);
          double max=vector_max_value<std::vector<double>,double>
            (diffs.size(),diffs);
          std::vector<double> len_list={min/10.0,max*10.0};
          ptemp.push_back(len_list);
        }
        vector<double> l10_list={-15,-9};
        ptemp.push_back(l10_list);
        param_lists.push_back(ptemp);
      }
      
      iko.full_min=true;
      iko.def_mmin.verbose=1;
      iko.set_covar(vmfrn,param_lists);
      
      iko.set_data(2,col_list_y.size(),table_obj.get_nlines(),
                   mvt_x,mvt_y);
      
      for(size_t i=0;i<table3d_obj.get_nx();i++) {
        for(size_t j=0;j<table3d_obj.get_ny();j++) {
          ubvector px(2);
          px[0]=table3d_obj.get_grid_x(i);
          px[1]=table3d_obj.get_grid_y(j);
          ubvector py(in.size()-6);
          iko.eval(px,py);
          for(size_t k=6;k<in.size();k++) {
            table3d_obj.set(i,j,k-6,py[k-6]);
          }
        }
      }
      
      command_del(type);
      clear_obj();
      command_add("table3d");
      type="table3d";
      
    } else if (method=="py_gp" || method=="py_dnn") {
      
#ifdef O2SCL_SET_PYTHON
      
      tensor<> tin, tout;
      vector<size_t> in_size={table_obj.get_nlines(),2};
      vector<size_t> out_size={table_obj.get_nlines(),in.size()-6};
      tin.resize(2,in_size);
      tout.resize(2,out_size);
      for(size_t j=0;j<table_obj.get_nlines();j++) {
        vector<size_t> ix={j,0};
        tin.get(ix)=table_obj.get(in[2],j);
        ix={j,1};
        tin.get(ix)=table_obj.get(in[3],j);
        for(size_t k=6;k<in.size();k++) {
          ix={j,k-6};
          tout.get(ix)=table_obj.get(in[k],j);
        }
      }
      o2scl_settings.py_init();
      interpm_python<> ip;
      if (method=="py_gp") {
        ip.set_function("interpm_sklearn_gp","verbose=2",2);
        ip.set_data_tensor(2,table_obj.get_nlines(),in.size()-6,
                           tin,tout);
      } else {
        ip.set_function("interpm_tf_dnn","verbose=2",2);
        ip.set_data_tensor(2,table_obj.get_nlines(),in.size()-6,
                           tin,tout);
      }
        
      for(size_t i=0;i<table3d_obj.get_nx();i++) {
        for(size_t j=0;j<table3d_obj.get_ny();j++) {
          vector<double> px={table3d_obj.get_grid_x(i),
            table3d_obj.get_grid_y(j)};
          vector<double> py(in.size()-6);
          ip.eval_std_vec(px,py);
          for(size_t k=6;k<in.size();k++) {
            table3d_obj.set(i,j,k-6,py[k-6]);
          }
        }
      }
      
      command_del(type);
      clear_obj();
      command_add("table3d");
      type="table3d";
#else
      cerr << "Python support not included." << endl;
      return 3;
#endif
      
    } else {
      cerr << "Method " << method << " not understood in interp-table3d."
           << endl;
      return 2;
    }
    
  } else {
    cerr << "Command interp-table3d implemented for type " << type << endl;
    return 1;
  }    
  
  return 0;
}

int acol_manager::comm_interactive(std::vector<std::string> &sv, 
				   bool itive_com) {

  post_interactive=!post_interactive;
  if (verbose>0) {
    if (post_interactive) {
      cout << "Interactive mode will run after command-line is parsed." 
	   << endl;
    } else {
      cout << "Interactive mode will not run after command-line is parsed." 
	   << endl;
    }
  }
  return 0;
}

