/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

int acol_manager::comm_list(std::vector<std::string> &sv, bool itive_com) {
  cout.precision(prec);
  if (type=="table3d") {
    cout << "table3d name: " << obj_name << endl;
    table3d_obj.summary(&cout,ncols);
  } else if (type=="table") {
    cout << "table name: " << obj_name << endl;
    if (table_obj.get_nunits()>0) {
      table_obj.summary(&cout,ncols);
    } else {
      table_obj.table<std::vector<double> >::summary(&cout,ncols);
    }
  } else if (type=="tensor") {
    cout << "tensor name: " << obj_name << endl;
    size_t rk=tensor_obj.get_rank();
    cout << "Rank: " << rk << endl;
    const std::vector<size_t> &sarr=tensor_obj.get_size_arr();
    for(size_t j=0;j<rk;j++) {
      cout << "Size of rank " << j << " is " << sarr[j] << endl;
    }
  } else if (type=="tensor<int>") {
    cout << "tensor<int> name: " << obj_name << endl;
    size_t rk=tensor_int_obj.get_rank();
    cout << "Rank: " << rk << endl;
    const std::vector<size_t> &sarr=tensor_int_obj.get_size_arr();
    for(size_t j=0;j<rk;j++) {
      cout << "Size of rank " << j << " is " << sarr[j] << endl;
    }
  } else if (type=="tensor<size_t>") {
    cout << "tensor<size_t> name: " << obj_name << endl;
    size_t rk=tensor_size_t_obj.get_rank();
    cout << "Rank: " << rk << endl;
    const std::vector<size_t> &sarr=tensor_size_t_obj.get_size_arr();
    for(size_t j=0;j<rk;j++) {
      cout << "Size of rank " << j << " is " << sarr[j] << endl;
    }
  } else if (type=="tensor_grid") {
    cout << "tensor_grid name: " << obj_name << endl;
    size_t rk=tensor_grid_obj.get_rank();
    cout << "Rank: " << rk << endl;
    const std::vector<size_t> &sarr=tensor_grid_obj.get_size_arr();
    if (tensor_grid_obj.is_grid_set()==false) {
      for(size_t j=0;j<rk;j++) {
	cout << "Size of rank " << j << " is " << sarr[j] << endl;
      }
      cout << "Grid not set." << endl;
    } else {
      for(size_t j=0;j<rk;j++) {
	cout << "Grid " << j << " (" << sarr[j] << "): ";
	cout << tensor_grid_obj.get_grid(j,0) << " ";
	if (sarr[j]>1) {
	  cout << tensor_grid_obj.get_grid(j,1) << " ";
	}
	if (sarr[j]>2) {
	  cout << "... " << tensor_grid_obj.get_grid(j,sarr[j]-1);
	}
	cout << endl;
      }
    }
  } else {
    cerr << "Cannot 'list' for type " << type << "." << endl;
    return exc_efailed;
  }
  return 0;
}

int acol_manager::comm_nlines(std::vector<std::string> &sv, 
			      bool itive_com) {
  if (type!="table") {
    cerr << "No table in 'nlines'." << endl;
    return 1;
  }

  if (table_obj.is_constant("nlines")) {
    cerr << "Constant 'nlines' already exists." << endl;
    return 2;
  }

  table_obj.add_constant("nlines",table_obj.get_nlines());
  
  return 0;
}

int acol_manager::comm_nrf(std::vector<std::string> &sv, 
			   bool itive_com) {

  table3d_obj.clear();

  if (sv.size()<12) {
    cerr << "Not enough arguments for nrf()." << endl;
    return 2;
  }
  
  double x_min=o2scl::function_to_double(sv[1]);
  double x_max=o2scl::function_to_double(sv[2]);
  double y_min=o2scl::function_to_double(sv[3]);
  double y_max=o2scl::function_to_double(sv[4]);
  size_t n_grid=o2scl::stoszt(sv[5]);
  string fx=sv[6];
  string fy=sv[7];
  string dfxdx=sv[8];
  string dfxdy=sv[9];
  string dfydx=sv[10];
  string dfydy=sv[11];

  uniform_grid<double> grid_x=uniform_grid_end<double>(x_min,x_max,n_grid-1);
  uniform_grid<double> grid_y=uniform_grid_end<double>(y_min,y_max,n_grid-1);

  table3d_obj.set_xy("x",grid_x,"y",grid_y);
  table3d_obj.new_slice("it");
  table3d_obj.new_slice("root");
  table3d_obj.set_slice_all("it",0.0);
  table3d_obj.set_slice_all("root",0.0);

  vector<size_t> max_count, min_count;

  vector<double> x_roots, y_roots;
  calculator c_fx, c_fy, c_dfxdx, c_dfxdy, c_dfydx, c_dfydy;
  c_fx.compile(fx.c_str());
  c_fy.compile(fy.c_str());
  c_dfxdx.compile(dfxdx.c_str());
  c_dfxdy.compile(dfxdy.c_str());
  c_dfydx.compile(dfydx.c_str());
  c_dfydy.compile(dfydy.c_str());
  std::map<std::string, double> vars;
  
  static const size_t kmax=1000;

  for(size_t i=0;i<n_grid;i++) {
    cout << "nrf progress: " << i+1 << "/" << n_grid << endl;
    for(size_t j=5;j<n_grid;j++) {
      double x=grid_x[i];
      double y=grid_y[j];
      bool found=false;
      for(size_t k=0;k<kmax && found==false;k++) {
	vars["x"]=x;
	vars["y"]=y;
	// x_{n+1} = x_n - J(x_n)^{-1} f(x_n)
	double fxd=c_fx.eval(&vars);
	double fyd=c_fy.eval(&vars);
	cout << fxd << " " << fyd << endl;
	// The Jacobian
	double a=c_dfxdx.eval(&vars);
	double b=c_dfxdy.eval(&vars);
	double c=c_dfydx.eval(&vars);
	double d=c_dfydy.eval(&vars);
	std::cout << a << " " << b << " " << c << " " << d << std::endl;
	// The inverse of the Jacobian
	double det=a*d-b*c;
	double ai=d/det;
	double bi=-b/det;
	double ci=-c/det;
	double di=a/det;
	std::cout << ai << " " << bi << " "
		  << ci << " " << di << std::endl;
	// The NR step
	double x1=x-ai*fxd-bi*fyd;
	double y1=y-ci*fxd-di*fyd;
	std::cout << x << " " << y << std::endl;
	std::cout << x1 << " " << y1 << std::endl;
	exit(-1);
	if (fabs(x1-x)+fabs(y1-y)<1.0e-14) {
	  found=true;
	  bool root_found=false;
	  for(size_t m=0;m<x_roots.size();m++) {
	    if (fabs(x1-x_roots[m])<1.0e-13 &&
		fabs(y1-y_roots[m])<1.0e-13) {
	      root_found=true;
	      table3d_obj.set(i,j,"it",k+1);
	      table3d_obj.set(i,j,"root",m+1);
	      if (k+1>max_count[m]) {
		max_count[m]=k+1;
	      }
	      if (k+1<min_count[m]) {
		min_count[m]=k+1;
	      }
	    }
	  }
	  if (root_found==false) {
	    size_t m=x_roots.size();
	    x_roots.push_back(x1);
	    y_roots.push_back(y1);
	    max_count.push_back(k+1);
	    min_count.push_back(k+1);
	    table3d_obj.set(i,j,"it",k+1);
	    table3d_obj.set(i,j,"root",m+1);
	  }
	}
	x=x1;
	y=y1;
      }
      if (found==false) {
	table3d_obj.set(i,j,"it",kmax);
      }
    }
  }

  for(size_t m=0;m<max_count.size();m++) {
    table3d_obj.add_constant(((string)"max_")+o2scl::szttos(m+1),
			     max_count[m]);
    table3d_obj.add_constant(((string)"min_")+o2scl::szttos(m+1),
			     min_count[m]);
  }

  table3d_obj.add_constant("k_max",kmax);
  table3d_obj.add_constant("n_roots",x_roots.size());
  
  // Delete previous object
  command_del(type);
  clear_obj();

  // Set new type to table3d
  command_add("table3d");
  type="table3d";
  
  return 0;
}

int acol_manager::comm_output(std::vector<std::string> &sv, bool itive_com) {

  if (type.length()==0) {
    cerr << "No object to output." << endl;
    return 3;
  }

  //--------------------------------------------------------------------
  // Output formatting
  
  ostream *fout;
  ofstream ffout;
  
  if (sv.size()==1) {
    fout=&cout;
  } else {
    ffout.open(sv[1].c_str());
    fout=&ffout;
  }
  
  if (scientific) fout->setf(ios::scientific);
  else fout->unsetf(ios::scientific);
  fout->precision(prec);

  if (type=="table3d") {
    
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    if (nx!=0 && ny!=0) {

      (*fout) << "Grid x: ";
      for(size_t i=0;i<((size_t)nx);i++) {
	(*fout) << table3d_obj.get_grid_x(i) << " ";
      }
      (*fout) << endl;
      (*fout) << "Grid y: ";
      for(size_t i=0;i<((size_t)ny);i++) {
	(*fout) << table3d_obj.get_grid_y(i) << " ";
      }
      (*fout) << endl;
      
      size_t nt=table3d_obj.get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  (*fout) << "Slice " << k << ": ";
	  (*fout) << table3d_obj.get_slice_name(k) << endl;
	  if (k==0) {
	    (*fout) << "Outer loops over x grid, inner loop over y grid." 
		    << endl;
	  }
	  for(size_t i=0;i<nx;i++) {
	    for(size_t j=0;j<ny;j++) {
	      (*fout) << table3d_obj.get(i,j,k) << " ";
	    }
	    (*fout) << endl;
	  }
	  fout->unsetf(ios::showpos);
	}
      }
    }

    return 0;

  } else if (type=="table") {
    
    if (table_obj.get_ncolumns()>0) {

      //--------------------------------------------------------------------
      // Count column widths

      vector<size_t> col_wids(table_obj.get_ncolumns());

      for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	col_wids[i]=prec+6;
      }

      if (names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	  if (table_obj.get_column_name(i).size()>col_wids[i]) {
	    col_wids[i]=table_obj.get_column_name(i).size();
	  }
	  std::string tunit=table_obj.get_unit(table_obj.get_column_name(i));
	  if (tunit.size()+2>col_wids[i]) {
	    col_wids[i]=tunit.size()+2;
	  }
	}
      }

      //--------------------------------------------------------------------
      // Output column names
      
      if (names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	  
	  // Preceeding space
	  if (pretty==true) {
	    (*fout) << ' ';
	  }
	  
	  // Column name
	  (*fout) << table_obj.get_column_name(i) << " ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=col_wids[i]-table_obj.get_column_name(i).size();
	    if (nsp<0) {
	      O2SCL_ERR("Column size anomaly (col names).",exc_efailed);
	    }
	    for(int j=0;j<nsp;j++) (*fout) << ' ';
	  }
	
	}
	(*fout) << endl;
      }
      
      //--------------------------------------------------------------------
      // Output units
    
      if (table_obj.get_nunits()>0 && names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {

	  // Preceeding space
	  if (pretty==true) {
	    (*fout) << ' ';
	  }
	
	  // Unit name
	  (*fout) << '['
		  << table_obj.get_unit(table_obj.get_column_name(i)) << "] ";
	
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=col_wids[i]-
	      table_obj.get_unit(table_obj.get_column_name(i)).size()-2;
	    if (nsp<0) {
	      O2SCL_ERR("Column size anomaly (units).",exc_efailed);
	    }
	    for(int j=0;j<nsp;j++) (*fout) << ' ';
	  }
	
	}
	(*fout) << endl;
      }
      
      //--------------------------------------------------------------------
      // Output data
      
      for(int i=0;i<((int)table_obj.get_nlines());i++) {
	
	for(size_t j=0;j<table_obj.get_ncolumns();j++) {
	  
	  // Otherwise, for normal output
	  if (pretty==true && table_obj.get(j,i)>=0.0) {
	    (*fout) << ' ';
	  }
	  (*fout) << table_obj.get(j,i) << ' ';
	  if (pretty==true) {
	    int nsp=((int)(table_obj.get_column_name(j).size()-prec-6));
	    for(int kk=0;kk<nsp;kk++) {
	      (*fout) << ' ';
	    }
	  }
	
	}
	(*fout) << endl;
      
	//--------------------------------------------------------------------
	// Continue to next row
      }
    
    }

  } else if (type=="hist") {
    
    for(size_t k=0;k<hist_obj.size();k++) {
      (*fout) << hist_obj.get_bin_low_i(k) << " ";
    }
    (*fout) << hist_obj.get_bin_high_i(hist_obj.size()-1) << endl;
    for(size_t k=0;k<hist_obj.size();k++) {
      (*fout) << hist_obj.get_wgt_i(k) << endl;
    }
    
  } else if (type=="int") {
    
    (*fout) << int_obj << endl;
    
  } else if (type=="char") {
    
    (*fout) << char_obj << endl;
    
  } else if (type=="double") {
    
    (*fout) << double_obj << endl;
    
  } else if (type=="size_t") {

    (*fout) << size_t_obj << endl;

  } else if (type=="string") {

    (*fout) << string_obj << endl;

  } else if (type=="int[]") {

    if (pretty) {
      vector<string> sv, sv_out;
      for(size_t k=0;k<intv_obj.size();k++) {
	sv.push_back(o2scl::itos(intv_obj[k])+' ');
      }
      screenify(intv_obj.size(),sv,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),intv_obj,true);
    }

  } else if (type=="double[]") {

    if (pretty) {
      vector<string> sv, sv_out;
      for(size_t k=0;k<doublev_obj.size();k++) {
	if (has_minus_sign(&doublev_obj[k])) {
	  sv.push_back(o2scl::dtos(doublev_obj[k])+' ');
	} else {
	  sv.push_back(" "+o2scl::dtos(doublev_obj[k])+' ');
	}
      }
      screenify(doublev_obj.size(),sv,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),doublev_obj,true);
    }

  } else if (type=="size_t[]") {

    if (pretty) {
      vector<string> sv, sv_out;
      for(size_t k=0;k<size_tv_obj.size();k++) {
	sv.push_back(o2scl::szttos(size_tv_obj[k])+' ');
      }
      screenify(size_tv_obj.size(),sv,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),size_tv_obj,true);
    }
    
  } else if (type=="string[]") {
    
    for(size_t k=0;k<stringv_obj.size();k++) {
      (*fout) << stringv_obj[k] << endl;
    }
    (*fout) << endl;

  } else if (type=="hist_2d") {
    
    for(size_t k=0;k<hist_2d_obj.size_x();k++) {
      (*fout) << hist_2d_obj.get_x_low_i(k) << " ";
    }
    (*fout) << hist_2d_obj.get_x_high_i(hist_2d_obj.size_x()-1) << endl;

    for(size_t k=0;k<hist_2d_obj.size_y();k++) {
      (*fout) << hist_2d_obj.get_y_low_i(k) << " ";
    }
    (*fout) << hist_2d_obj.get_y_high_i(hist_2d_obj.size_y()-1) << endl;

    for(size_t ki=0;ki<hist_2d_obj.size_x();ki++) {
      for(size_t kj=0;kj<hist_2d_obj.size_y();kj++) {
	(*fout) << hist_2d_obj.get_wgt_i(ki,kj) << " ";
      }
      (*fout) << endl;
    }

  } else if (type=="vector<contour_line>") {

    (*fout) << cont_obj.size() << endl;
    for(size_t k=0;k<cont_obj.size();k++) {
      (*fout) << cont_obj[k].level << " " << cont_obj[k].x.size() << endl;
      for(size_t kk=0;kk<cont_obj[k].x.size();kk++) {
	(*fout) << cont_obj[k].x[kk] << " ";
	(*fout) << cont_obj[k].y[kk] << endl;
      }
    }

  } else if (type=="tensor") {

    tensor_out(*fout,tensor_obj,pretty);
    
  } else if (type=="tensor<int>") {

    tensor_out(*fout,tensor_int_obj,pretty);
    
  } else if (type=="tensor<size_t>") {

    tensor_out(*fout,tensor_size_t_obj,pretty);
    
  } else if (type=="uniform_grid<double>") {

    (*fout) << ug_obj.get_nbins() << " ";
    (*fout) << ug_obj.get_start() << " ";
    (*fout) << ug_obj.get_end() << " ";
    (*fout) << ug_obj.get_width() << endl;

  } else {

    cerr << "Cannot output type " << type << "." << endl;
    return 2;
    
  }

  if (sv.size()!=1) {
    ffout.close();
  }
  
  return 0;
}

int acol_manager::comm_max(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (table3d_obj.get_nslices()==0) {
      cerr << "No slices to find the maximum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find maximum value of",
			  i1,"max",itive_com);
    if (ret!=0) return ret;

    size_t ix;
    if (!table3d_obj.is_slice(i1,ix)) {
      cerr << "No slice named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const ubmatrix &mat=table3d_obj.get_slice(ix);
    size_t i, j;
    double max;
    matrix_max_index(mat,i,j,max);

    cout << "Maximum value of slice '" << i1 << "' is: " 
	 << max << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="hist_2d") {
    
    const ubmatrix &mat=hist_2d_obj.get_wgts();
    size_t i, j;
    double max;
    matrix_max_index(mat,i,j,max);
    
    cout << "Maximum weight is: "
	 << max << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to find the maximum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column to find maximum value of",
			  i1,"max",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    double max;
    size_t ix;
    vector_max(table_obj.get_nlines(),(table_obj)[i1],ix,max);
    cout << "Maximum value of column '" << i1 << "' is: " 
	 << max << " at row with index " << ix << "." << endl;

  } else if (type=="double[]") {

    double val;
    size_t loc;
    o2scl::vector_max<vector<double>,double>(doublev_obj.size(),
					     doublev_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="tensor") {

    double val;
    size_t loc;
    const vector<double> &v=tensor_obj.get_data();
    o2scl::vector_max<vector<double>,double>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_obj.get_rank());
    tensor_obj.unpack_index(loc,ix);
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor<size_t>") {

    size_t val;
    size_t loc;
    const vector<size_t> &v=tensor_size_t_obj.get_data();
    o2scl::vector_max<vector<size_t>,size_t>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_size_t_obj.get_rank());
    tensor_size_t_obj.unpack_index(loc,ix);
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor<int>") {

    int val;
    size_t loc;
    const vector<int> &v=tensor_int_obj.get_data();
    o2scl::vector_max<vector<int>,int>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_int_obj.get_rank());
    tensor_int_obj.unpack_index(loc,ix);
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor_grid") {

    double val;
    size_t loc;
    const vector<double> &v=tensor_grid_obj.get_data();
    o2scl::vector_max<vector<double>,double>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_grid_obj.get_rank());
    vector<double> dx(tensor_grid_obj.get_rank());
    tensor_grid_obj.unpack_index(loc,ix);
    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
      dx[j]=tensor_grid_obj.get_grid(j,ix[j]);
    }
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    cout << "  and grid point ";
    vector_out(cout,dx,true);
    
  } else if (type=="int[]") {
    
    int val;
    size_t loc;
    o2scl::vector_max<vector<int>,int>(intv_obj.size(),
				       intv_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="size_t[]") {
    
    size_t val;
    size_t loc;
    o2scl::vector_max<vector<size_t>,size_t>(size_tv_obj.size(),
					     size_tv_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_min(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d" || table3d_obj.get_nslices()==0) {
      cerr << "No table3d with slices to find the minimum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find minimum of",
			  i1,"min",itive_com);
    if (ret!=0) return ret;

    size_t ix;
    if (!table3d_obj.is_slice(i1,ix)) {
      cerr << "No slice named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const ubmatrix &mat=table3d_obj.get_slice(ix);
    size_t i, j;
    double min;
    matrix_min_index(mat,i,j,min);

    cout << "Minimum value of slice '" << i1 << "' is: " 
	 << min << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="hist_2d") {
    
    const ubmatrix &mat=hist_2d_obj.get_wgts();
    size_t i, j;
    double min;
    matrix_min_index(mat,i,j,min);
    
    cout << "Minimum weight is: "
	 << min << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to find the minimum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column to find minimum of",
			  i1,"min",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    double min;
    size_t ix;
    vector_min(table_obj.get_nlines(),(table_obj)[i1],ix,min);
    cout << "Minimum value of column '" << i1 << "' is: " 
	 << min << " at row with index " << ix << "." << endl;
    
  } else if (type=="double[]") {

    double val;
    size_t loc;
    o2scl::vector_min<vector<double>,double>(doublev_obj.size(),
					     doublev_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="tensor") {

    double val;
    size_t loc;
    const vector<double> &v=tensor_obj.get_data();
    o2scl::vector_min<vector<double>,double>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_obj.get_rank());
    tensor_obj.unpack_index(loc,ix);
    cout << "Minimum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor<size_t>") {

    size_t val;
    size_t loc;
    const vector<size_t> &v=tensor_size_t_obj.get_data();
    o2scl::vector_max<vector<size_t>,size_t>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_size_t_obj.get_rank());
    tensor_size_t_obj.unpack_index(loc,ix);
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor<int>") {

    int val;
    size_t loc;
    const vector<int> &v=tensor_int_obj.get_data();
    o2scl::vector_max<vector<int>,int>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_int_obj.get_rank());
    tensor_int_obj.unpack_index(loc,ix);
    cout << "Maximum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    
  } else if (type=="tensor_grid") {

    double val;
    size_t loc;
    const vector<double> &v=tensor_grid_obj.get_data();
    o2scl::vector_min<vector<double>,double>(v.size(),v,loc,val);
    vector<size_t> ix(tensor_grid_obj.get_rank());
    vector<double> dx(tensor_grid_obj.get_rank());
    tensor_grid_obj.unpack_index(loc,ix);
    for(size_t j=0;j<tensor_grid_obj.get_rank();j++) {
      dx[j]=tensor_grid_obj.get_grid(j,ix[j]);
    }
    cout << "Minimum value is " << val << " at indices ";
    vector_out(cout,ix,true);
    cout << "  and grid point ";
    vector_out(cout,dx,true);
    
  } else if (type=="int[]") {
    
    int val;
    size_t loc;
    o2scl::vector_min<vector<int>,int>(intv_obj.size(),
				       intv_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="size_t[]") {
    
    size_t val;
    size_t loc;
    o2scl::vector_min<vector<size_t>,size_t>(size_tv_obj.size(),
					     size_tv_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  }
  
  return 0;
}

