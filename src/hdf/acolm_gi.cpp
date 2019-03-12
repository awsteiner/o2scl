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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int acol_manager::comm_get_conv
(std::vector<std::string> &sv, bool itive_com) {
  
  vector<string> in, pr;
  pr.push_back("Old unit");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"get-conv",itive_com);
  if (ret!=0) return ret;
  
  if (unit_fname.length()>0) {
    cng.units_cmd_string=((string)"units -f ")+unit_fname;
    if (verbose>=2) {
      cout << "Units command string: " << cng.units_cmd_string
	   << endl;
    }
  }
  
  // Set the proper output precision and mode
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);

  if (verbose>=2) {
    cng.verbose=1;
  } else {
    cng.verbose=0;
  }

  // If cng.verbose is non-zero, then cng.convert may output
  // verbose information to cout
  double val;
  int cret=cng.convert_ret(in[0],in[1],1.0,val);
  if (cret!=0) {
    cerr << "Conversion failed." << endl;
    return 1;
  }
  
  cout << "Conversion factor is: " << val << endl;
  
  return 0;
}

int acol_manager::comm_get_grid(std::vector<std::string> &sv, bool itive_com) {

  if (type=="tensor_grid") {

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
	    o2scl::dtos(tensor_grid_obj.get_grid(k,ell),prec);
	}
      }
    }

    columnify col;
    vector<string> aligned(max_size+1);
    col.align(string_mat,rank+1,max_size+1,aligned,align_spec);
    for(size_t i=0;i<aligned.size();i++) {
      cout << aligned[i] << endl;
    }

    return 0;
  }

  cout << "Not implemented for type " << type << endl;
  return 1;
}

int acol_manager::comm_get_unit(std::vector<std::string> &sv, bool itive_com) {

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

int acol_manager::comm_get_row(std::vector<std::string> &sv, bool itive_com) {

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

  if (user_ncols<=0) {
    char *ncstring=getenv("COLUMNS");
    if (ncstring) {
      int nc2;
      int sret=o2scl::stoi_nothrow(ncstring,nc2);
      if (sret==0 && nc2>0) {
	ncols=nc2;
      } else {
	cerr << "Failed to interpret COLUMNS value " << ncstring
	     << " as a positive number of columns." << endl;
      }
    }
  } else {
    ncols=user_ncols;
  }

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
      str.precision(prec);
      
      for(size_t i=0;i<table_obj.get_ncolumns();i++) {

	// Count for space between columns and sign
	size_t this_col=2;
	// Count column name
	this_col+=table_obj.get_column_name(i).size();
	// Count extra spaces to format number
	int num_spaces=prec+6-((int)(table_obj.get_column_name(i).size()));
	if (num_spaces>0) this_col+=num_spaces;
	// See if there will be space
	if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	  row_names.push_back(str.str());
	  str.str("");
	  str.clear();
	  str.setf(ios::scientific);
	  str.precision(prec);
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
      
      cout.precision(prec);
  
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
    str.precision(prec);
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      
      // Count space for number
      size_t this_col=prec+8;
      // Count extra spaces if necessary
      int num_spaces=((int)(table_obj.get_column_name(i).size())-prec-6);
      if (num_spaces>0) this_col+=num_spaces;
      // See if there will be space
      if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	row_data.push_back(str.str());
	str.str("");
	str.clear();
	str.setf(ios::scientific);
	str.precision(prec);
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

int acol_manager::comm_generic(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype;

  // Delete previous object
  command_del();
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
  
  if (ctype=="table") {
    
    if (fname!=((std::string)"cin")) {
      table_obj.read_generic(ifs,verbose);
    } else {
      table_obj.read_generic(std::cin,verbose);
    }

  } else if (ctype=="table3d") {
    
    if (fname!=((std::string)"cin")) {
      table3d_obj.read_gen3_list(ifs,verbose);
    } else {
      table3d_obj.read_gen3_list(std::cin,verbose);
    }
    
  } else if (ctype=="int") {

    if (fname!=((std::string)"cin")) {
      ifs >> int_obj;
    } else {
      cin >> int_obj;
    }
    
  } else if (ctype=="char") {

    if (fname!=((std::string)"cin")) {
      ifs >> char_obj;
    } else {
      cin >> char_obj;
    }
    
  } else if (ctype=="double") {

    if (fname!=((std::string)"cin")) {
      ifs >> double_obj;
    } else {
      cin >> double_obj;
    }
    
  } else if (ctype=="size_t") {

    if (fname!=((std::string)"cin")) {
      ifs >> size_t_obj;
    } else {
      cin >> size_t_obj;
    }
    
  } else if (ctype=="string") {

    if (fname!=((std::string)"cin")) {
      getline(ifs,string_obj);
    } else {
      getline(cin,string_obj);
    }
    
  } else if (ctype=="int[]") {

    if (fname!=((std::string)"cin")) {
      int itmp;
      intv_obj.clear();
      while (ifs >> itmp) {
	intv_obj.push_back(itmp);
      }
    } else {
      int itmp;
      intv_obj.clear();
      while (cin >> itmp) {
	intv_obj.push_back(itmp);
      }
    }
    
  } else if (ctype=="double[]") {

    if (fname!=((std::string)"cin")) {
      double dtmp;
      doublev_obj.clear();
      while (ifs >> dtmp) {
	doublev_obj.push_back(dtmp);
      }
    } else {
      double dtmp;
      doublev_obj.clear();
      while (cin >> dtmp) {
	doublev_obj.push_back(dtmp);
      }
    }
    
  } else if (ctype=="size_t[]") {
    
    if (fname!=((std::string)"cin")) {
      size_t sttmp;
      size_tv_obj.clear();
      while (ifs >> sttmp) {
	size_tv_obj.push_back(sttmp);
      }
    } else {
      size_t sttmp;
      size_tv_obj.clear();
      while (cin >> sttmp) {
	size_tv_obj.push_back(sttmp);
      }
    }
    
  } else if (ctype=="string[]") {

    if (fname!=((std::string)"cin")) {
      std::string stmp;
      while (getline(ifs,stmp)) {
	stringv_obj.push_back(stmp);
      }
    } else {
      std::string stmp;
      while (getline(cin,stmp)) {
	stringv_obj.push_back(stmp);
      }
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

int acol_manager::comm_help(std::vector<std::string> &sv, bool itive_com) {

  // Handle the 'help type command' case for type-specific commands
  if (sv.size()==3) {
    string temp_type=sv[1];
    string cur_type=type;

    command_del();
    command_add(temp_type);
    
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
    
    int ret=cl->comm_option_help(sv,itive_com);

    command_del();
    command_add(cur_type);
    return ret;
  }

  // Handle the special case 'help functions'
  if (sv.size()==2 && sv[1]=="functions") {
    string str=((std::string)"Functions can be created using the ");
    str+="operators and functions listed below. Examples are ";
    str+="\"x==5 && y<1\", \"acos(-1)\", and \"sin(x>5)\".\n\n";
    str+="Operators:\n\n() ^ * / % + - == != < > && || << >> >= <=\n\n";
    str+="Functions:\n\n";
    str+="exp(x) log(x) log10(x) sin(x) cos(x) tan(x) sqrt(x) abs(x) ";
    str+="asin(x) acos(x) atan(x) sinh(x) cosh(x) tanh(x) ";
    str+="asinh(x) acosh(x) atanh(x)\n\n";
    /*
      dsc+="atan2(x,y) if(x,y,z)\n";
      dsc+="cot(x) csc(x) sec(x)\n";
      dsc+="ceil(x) floor(x) int(x) max(x,y) min(x,y)\n";
    */
    std::vector<std::string> sv;
    o2scl::rewrap_keep_endlines(str,sv);
    for(size_t i=0;i<sv.size();i++) {
      cout << sv[i] << endl;
    }
      
    return 0;
  }
  
  // Handle the special case 'help vector-spec'
  if (sv.size()==2 && sv[1]=="vector-spec") {
    
    std::string str=((std::string)"Vector specification ")+
      "description:\n\nSome acol commands take arguments "+
      "which are 'vector specifications', i.e. an "+
      "array specified as a string. The different parts of the string "+
      "are separated by a colon, and the first part specifes the type "+
      "of vector specification. The different types are:\n\n"+
      "  1. val:<value> - Create a vector with one element equal to "+
      "<value>, which\nmay be a number or a simple function, e.g. "+
      "'val:sin(0.5)'.\n\n"+
      "  2. list:<entry 0>,<entry 1>, ..., <entry n-1> - Create a "+
      "vector with a\nsimple list of numbers or functions, e.g. "+
      "'list:3.0,1.0e-3,sqrt(2.0)'.\n\n"+
      "  3. func:<N>:<function of i> - Create a vector by specifying "+
      "the length of\nthe vector and a function used to fill the "+
      "elements. For example: 'func:41:sin(i/20.0*acos(-1))'.\n\n"+
      "  4. grid:<begin>:<end>:<width>:[\"log\"] - Create a "+
      "vector equal to a\nuniform grid, e.g. use 'grid:1.0,10.0,1.0' "+
      "for a 10-element vector filled with the numbers 1 to 10.\n\n"+
      "  5. text:<filename>:<column> - Read a text file and extract "+
      "a vector of\nnumbers from a column of the text file, ignoring "+
      "any header rows which contain non-numeric values. For "+
      "example 'text:~/temp.dat:2' will construct a vector from the "+
      "second row of the file 'temp.dat' in the user's home directory.\n\n"+
      "  6. hdf5:<file name>:<object name>:[additional spec.] - "+
      "Read an HDF5 file\nand obtain a vector from the object with "+
      "the specified name. The remaining parts of the string contain "+
      "additional information which may be needed depending on the "+
      "type of object stored in the HDF5 file.";

    std::vector<std::string> sv;
    o2scl::rewrap_keep_endlines(str,sv);
    for(size_t i=0;i<sv.size();i++) {
      cout << sv[i] << endl;
    }
      
    return 0;
  }

  // Handle the special case 'help mult-vector-spec'
  if (sv.size()==2 && sv[1]=="mult-vector-spec") {
    
    std::string str=((std::string)"Multiple vector specification ")+
      "description:\n\nSome acol commands take arguments "+
      "which are 'multiple vector specifications', i.e. a set of "+
      "arrays specified as a string. The different parts of the string "+
      "are separated by a colon, and the first part specifes the type "+
      "of multiple vector specification. The different types are:\n\n"+
      "  1. func:<N>:<function of i>:<function of i and j> - "+
      "Specify the\nnumber of vectors, a function of \"i\" which "+
      "determines the length of the ith vector, and a function of "+
      "\"i\" and \"j\" which specifies the jth element of the ith "+
      "vector.\n\n"
      "  2. text:<filename pattern>:<column list> - Read one or\nmore "+
      "text files and extract "+
      "vectors of numbers from columns of the text file, ignoring "+
      "any header rows which contain non-numeric values. For "+
      "example 'text:~/temp.dat:2' will construct a vector from the "+
      "second row of the file 'temp.dat' in the user's home directory.\n\n"+
      "  3. hdf5:<filename pattern>:<object name>:[additional spec.] - "+
      "Read one or\nmore HDF5 files and obtain a vector from the object "+
      "with "+
      "the specified name. The remaining parts of the string contain "+
      "additional information which may be needed depending on the "+
      "type of object stored in the HDF5 file.";

    std::vector<std::string> sv;
    o2scl::rewrap_keep_endlines(str,sv);
    for(size_t i=0;i<sv.size();i++) {
      cout << sv[i] << endl;
    }
      
    return 0;
  }

  // Handle the case 'help <type>' where <type> is an acol o2scl type
  for(size_t i=0;i<type_list.size();i++) {
    if (sv[1]==type_list[i]) {
      std::string str="Type "+type_list[i]+" is one of the types of "+
	"objects which can be read, written, or modified by "+
	cl->cmd_name+".";
      std::vector<std::string> svx;
      o2scl::rewrap_keep_endlines(str,svx);
      for(size_t i=0;i<svx.size();i++) {
	cout << svx[i] << endl;
      }
      cout << endl;
      svx.clear();
      std::map<std::string,std::vector<std::string> >::iterator it;
      for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
	if (it->first==type_list[i]) {
	  std::vector<std::string> &clist=it->second;
	  if (clist.size()>1) {
	    str=((string)"The type-specific commands for \"")+
	      type_list[i]+"\" are: ";
	    for(size_t j=0;j<clist.size()-1;j++) {
	      str+=clist[j]+", ";
	    }
	    str+=" and "+clist[clist.size()-1]+".";
	    o2scl::rewrap_keep_endlines(str,svx);
	    for(size_t j=0;j<svx.size();j++) {
	      cout << svx[j] << endl;
	    }
	  } else {
	    cout << "The only type-specific command for \""
		 << type_list[i] << " is " << clist[0] << "." << endl;
	  }
	}
      }
      
      return 0;
    }
  }

  // Handle the case 'help <command>' where <command> is a type-specific
  // command
  if (sv.size()>2 && cl->is_valid_option(sv[1])==false) {

    bool found=false;
    
    std::map<std::string,std::vector<std::string> >::iterator it;
    for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
      std::vector<std::string> &clist=it->second;
      for(size_t j=0;j<clist.size();j++) {
	if (clist[j]==sv[1]) {

	  if (found==false) {
	    cout << "Command \"" << sv[1] << "\" is a type specific "
		 << "command. Below are "
		 << "the various descriptions\nof its operation with "
		 << "the relevant types.\n" << endl;
	    found=true;
	  }
	  
	  {
	    // Draw a line
	    ostringstream oss;
	    oss << ((char)27) << '(' << '0';
	    for(size_t i=0;i<78;i++) oss << 'q';
	    oss << ((char)27) << '(' << 'B';
	    string s=oss.str();
	    cout << s << endl;
	  }
	  
	  string cur_type=type;
	  
	  vector<string> sv2;
	  
	  sv2.push_back("help");
	  sv2.push_back(sv[1]);
	  
	  command_del();
	  command_add(it->first);
	  type=it->first;

	  cout << "Type " << it->first << ":" << endl;
	  int ret=cl->comm_option_help(sv2,itive_com);
	  cout << endl;

	  command_del();
	  command_add(cur_type);
	  type=cur_type;
	}
      }
    }

    if (found) return 0;
  }

  // Handle the usual cli case
  return cl->comm_option_help(sv,itive_com);
}

/*
  int acol_manager::comm_interp_type(std::vector<std::string> &sv, 
  bool itive_com) {

  if (type=="table3d") {
    
  if (type!="table3d") {
  cout << "No table to get interpolation type of in 'interp-type'." << endl;
  return exc_efailed;
  }

  if (sv.size()>1) {
  if (o2scl::stoi(sv[1])>7 || o2scl::stoi(sv[1])<0) {
  cout << "Invalid interpolation type in 'interp-type'." << endl;
  return exc_efailed;
  }
  table3d_obj.set_interp_type(o2scl::stoi(sv[1]));
  }
    
  if (sv.size()==1 || verbose>0) {
  size_t itype=table3d_obj.get_interp_type();
  cout << "Current interpolation type is " << itype;
      
  if (itype<4) {
  if (itype<3) {
  if (itype==1) {
  cout << " (linear)." << endl;
  } else {
  cout << " (cubic spline)." << endl;
  }
  } else {
  cout << " (cubic spline, periodic)." << endl;
  }
  } else {
  if (itype<6) {
  if (itype<5) {
  cout << " (Akima spline)." << endl;
  } else {
  cout << " (Akima spline, periodic)." << endl;
  }
  } else {
  if (itype<7) {
  cout << " (Monotonicity-preserving)." << endl;
  } else {
  cout << " (Steffen's monotonic)." << endl;
  }
  }
  }
  }
    
  return 0;
  }
  
  if (table_obj.get_nlines()==0) {
  cout << "No table to get interpolation type of." << endl;
  return exc_efailed;
  }

  if (sv.size()>1) {
  if (o2scl::stoi(sv[1])>7 || o2scl::stoi(sv[1])<0) {
  cout << "Invalid interpolation type in interp-type." << endl;
  return exc_efailed;
  }
  table_obj.set_interp_type(o2scl::stoi(sv[1]));
  }

  if (sv.size()==1 || verbose>0) {
  size_t itype=table_obj.get_interp_type();
  cout << "Current interpolation type is " << itype;
    
  if (itype<4) {
  if (itype<3) {
  if (itype==1) {
  cout << " (linear)." << endl;
  } else {
  cout << " (cubic spline)." << endl;
  }
  } else {
  cout << " (cubic spline, periodic)." << endl;
  }
  } else {
  if (itype<6) {
  if (itype<5) {
  cout << " (Akima spline)." << endl;
  } else {
  cout << " (Akima spline, periodic)." << endl;
  }
  } else {
  if (itype<7) {
  cout << " (Monotonicity-preserving)." << endl;
  } else {
  cout << " (Steffen's monotonic)." << endl;
  }
  }
  }
  }
  
  return 0;
  }
*/

int acol_manager::comm_integ(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

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

    hf.sets_vec(obj_name,stringv_obj);
    
  } else if (type=="table3d") {
    
    hdf_output(hf,((const table3d &)(table3d_obj)),obj_name);
    
  } else if (type=="tensor_grid") {
    
    hdf_output(hf,tensor_grid_obj,obj_name);
    
  } else if (type=="table") {
    
    hdf_output(hf,table_obj,obj_name);

  } else if (type=="hist") {

    hdf_output(hf,hist_obj,obj_name);
    
  } else if (type=="hist_2d") {

    hdf_output(hf,((const hist_2d &)(hist_2d_obj)),obj_name);
    
  } else if (type=="vector<contour_line>") {

    hdf_output(hf,cont_obj,obj_name);
    
  } else if (type=="uniform_grid<double>") {

    hdf_output(hf,ug_obj,obj_name);
    
  }

  hf.close();
    
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
  table_obj.new_column(i1);
  for(size_t i=0;i<table_obj.get_nlines();i++) table_obj.set(i1,i,((double)i));

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
      cout << "Read table3d  named " << in[1] << " from file " << in[0] << endl;
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

    table3d_obj.add_slice_from_table(tmp,in[2],in[3]);

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

  std::string in[3], pr[3]=
    {"Enter filename of external table (or blank to stop): ",
     "Enter index column in present table (or blank to stop): ",
     "Enter index column in external table (or blank to stop): "};
  if (sv.size()>=3) {
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
  } else {
    if (itive_com) {
      for(size_t is=0;is<3;is++) {
	in[is]=cl->cli_gets(pr[is].c_str());
	if (in[is].length()==0) {
	  cout << "Command 'insert' cancelled." << endl;
	  return 0;
	}
      }
    } else {
      cerr << "Not enough arguments to command 'insert'" << endl;
      return exc_efailed;
    }
  }

  cout << "Unimplemented." << endl;

  return 1;

#ifdef O2SCL_NEVER_DEFINED

  for (size_t j=0;j<tmp->get_ncolumns();j++) {
    string ty=tmp->get_column_name(j);
    int tret=table_obj.add_col_from_table(in[1],*tmp,in[2],ty,ty);
    if (tret!=0) {
      cerr << "Adding column " << ty << " failed." << endl;
      ret=tret;
      // We don't return here so that "delete tmp;" is called below
    }
  }
  
#endif
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

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,doublev_obj) << endl;

  } else if (type=="int[]") {

    size_t n=intv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(intv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,value) << endl;

  } else if (type=="size_t[]") {

    size_t n=size_tv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(size_tv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,value) << endl;

  } else {
    cout << "Not implemented for type " << type << endl;
    return 1;
  }    
  
  return 0;
}

