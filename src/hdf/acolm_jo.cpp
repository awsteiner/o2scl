/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/cursesw.h>
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/inte_double_exp_boost.h>
#include <o2scl/inte_adapt_cern.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/*
  int acol_manager::comm_list()
  int acol_manager::comm_max()
  int acol_manager::comm_min()
  int acol_manager::comm_ninteg()
  int acol_manager::comm_nlines()
  int acol_manager::comm_output()
*/

int acol_manager::comm_list(std::vector<std::string> &sv, bool itive_com) {

  cout.precision(precision);

  int ncols_loc;
  if (ncols<=0) {
    int srow, scol;
    int iret=get_screen_size_ioctl(srow,scol);
    //std::cout << "iret,srow,scol: " << iret << " " << srow << " "
    //<< scol << std::endl;
    if (scol>10 || iret!=0) ncols_loc=scol;
    else ncols_loc=80;
  } else {
    ncols_loc=ncols;
  }
  if (verbose>1) {
    cout << "Number of screen columns: " << ncols_loc << endl;
  }
  
  if (type=="table3d") {
    cout << "table3d name: " << obj_name << endl;
    table3d_obj.summary(&cout,ncols_loc);
  } else if (type=="table") {
    cout << "table name: " << obj_name << endl;
    if (table_obj.get_nunits()>0) {
      table_obj.summary(&cout,ncols_loc);
    } else {
      table_obj.table<std::vector<double> >::summary(&cout,ncols_loc);
    }
  } else if (type=="hist_2d") {
    cout << "hist_2d name: " << obj_name << endl;
    cout << "x y" << endl;
    size_t max=hist_2d_obj.size_x()+1;
    if (hist_2d_obj.size_y()+1>max) {
      max=hist_2d_obj.size_y()+1;
    }
    for(size_t i=0;i<max;i++) {
      if (max>=10000) cout.width(8);
      else if (max>=1000) cout.width(4);
      else if (max>=100) cout.width(3);
      else if (max>=10) cout.width(2);
      cout << i << " ";
      if (i<hist_2d_obj.size_x()) {
	cout << hist_2d_obj.get_x_low_i(i) << " ";
      } else if (i==hist_2d_obj.size_x()) {
	cout << hist_2d_obj.get_x_high_i(i-1) << " ";
      } else {
	cout << "       ";
        for(int jj=0;jj<precision;jj++) cout << ' ';
          
      }
      if (i<hist_2d_obj.size_y()) {
	cout << hist_2d_obj.get_y_low_i(i) << endl;
      } else if (i==hist_2d_obj.size_y()) {
	cout << hist_2d_obj.get_y_high_i(i-1) << endl;
      } else {
	cout << "       ";
        for(int jj=0;jj<precision;jj++) cout << ' ';
	cout << endl;
      }
    }
  } else if (type=="tensor") {
    cout << "tensor name: " << obj_name << endl;
    size_t rk=tensor_obj.get_rank();
    cout << "Rank: " << rk << endl;
    const std::vector<size_t> &sarr=tensor_obj.get_size_arr();
    for(size_t j=0;j<rk;j++) {
      cout << "Size of rank " << j << " is " << sarr[j] << endl;
    }
  } else if (type=="hist") {
    cout << "hist name: " << obj_name << endl;
    cout << hist_obj.size() << " bins" << endl;
    cout << "Bin edges: " << endl;
    vector_out(cout,hist_obj.get_bins(),true);
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
  } else if (type=="vector<contour_line>") {
    cout << "vector<contour_line> name: " << obj_name << endl;
    cout << cont_obj.size() << " contour lines." << endl;
    for(size_t j=0;j<cont_obj.size();j++) {
      cout << "Contour line " << j << " has level "
           << cont_obj[j].level << " and has "
           << cont_obj[j].x.size() << " points." << endl;
    }
  } else if (type=="vec_vec_double") {
    cout << "vector<vector<double>> name: " << obj_name << endl;
    cout << vvdouble_obj.size() << " entries." << endl;
    for(size_t j=0;j<vvdouble_obj.size();j++) {
      cout << "Entry " << j << " is an array of size "
           << vvdouble_obj[j].size() << endl;
    }
  } else if (type=="vec_vec_string") {
    cout << "vector<vector<string>> name: " << obj_name << endl;
    cout << vvstring_obj.size() << " entries." << endl;
    for(size_t j=0;j<vvstring_obj.size();j++) {
      size_t count=0;
      for(size_t k=0;k<vvstring_obj[j].size();k++) {
        count+=vvstring_obj[j][k].length();
      }
      cout << "Entry " << j << " is an array of size "
           << vvstring_obj[j].size() << " with " << count
           << " total characters." << endl;
    }
  } else if (type=="prob_dens_mdim_gmm") {
    cout << "prob_dens_mdim_gmm name: " << obj_name << endl;
    size_t n=pgmm_obj.pdmg.size();
    for(size_t j=0;j<n;j++) {
      const ubvector &peak=pgmm_obj.pdmg[j].get_peak();
      cout << "Gaussian " << j << " has mean ";
      vector_out(std::cout,peak,true);
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

int acol_manager::comm_ninteg(std::vector<std::string> &sv, bool itive_com) {

  vector<string> in, pr;
  std::string kw;
  
  if (sv.size()>5) {

    in.resize(4);
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
    in[3]=sv[4];
    kw=sv[5];

  } else if (sv.size()>4) {

    in.resize(4);
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
    in[3]=sv[4];
    
  } else {
    
    pr.push_back("Function");
    pr.push_back("Integration variable");
    pr.push_back("Lower limit");
    pr.push_back("Upper limit");
    pr.push_back("Additional arguments");
    int ret=get_input(sv,pr,in,"ninteg",itive_com);
    if (ret!=0) return ret;

    kw=in[4];
    
  }

  bool multiprecision=false;
  string method="kb";
  if (kw.length()>0) {
    kwargs kwa(kw);
    multiprecision=kwa.get_bool("multip",false);
    method=kwa.get_string("method","kb");
    if (verbose>1) {
      cout << "Using method: " << method << endl;
    }
  }
  
  if (sv.size()<5) {
    cerr << "Not enough arguments for ninteg." << endl;
    return 1;
  }
  std::string func=in[0];
  std::string var=in[1];

  inte_kronrod_boost<61> ikb;
  inte_multip_double_exp_boost ideb;
  inte_adapt_cern iac;

  if (multiprecision) {
    
#ifndef O2SCL_OSX

    std::cerr << "Multiprecision for ninteg only works for OSX "
              << "at the moment." << std::endl;
    return 5;

#else

    funct_multip_string fms;
    fms.set_function(func,var);
    funct_multip_string *fmsp=&fms;
    
    funct_multip fm2;
    
    // C++ prints out precision+1 significant figures so we add one to
    // 'precision' to construct the integration tolerance.
    if (method=="kb") {
      ikb.tol_rel_multip=pow(10.0,-precision-1);
      ikb.verbose=verbose;
    } else if (method=="deb") {
      ideb.tol_rel_multip=pow(10.0,-precision-1);
      ideb.verbose=verbose;
    } else if (method=="ac") {
      iac.tol_rel_multip=pow(10.0,-precision-1);
      if (verbose>0) {
        iac.verbose=verbose-1;
      } else {
        iac.verbose=0;
      }
    }
    
    if (precision>49) {
      
      cerr << "Requested precision too large for the ninteg "
           << "command (maximum is 49)." << endl;
      return 2;
      
    } else if (precision>34) {

      cpp_dec_float_50 d=0, err, lower_lim, upper_lim;
      convert_units<cpp_dec_float_50> cu;
      if (in[2]=="-infty") {
        lower_lim=-std::numeric_limits<cpp_dec_float_50>::infinity();
      } else {
        function_to_fp_nothrow(in[2],lower_lim,cu);
      }
      if (in[3]=="infty") {
        upper_lim=std::numeric_limits<cpp_dec_float_50>::infinity();
      } else {
        function_to_fp_nothrow(in[3],upper_lim,cu);
      }
      int retx;
      if (method=="kb") {
        retx=ikb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else if (method=="deb") {
        retx=ideb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else {
        retx=iac.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      }
    
    if (retx!=0) {
        cerr << "Integrating " << func << " failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_50): ";
      cout << dtos(d,precision) << endl;
      return 0;
      
    } else if (precision>24) {
      
      cpp_dec_float_35 d=0, err, lower_lim, upper_lim;
      convert_units<cpp_dec_float_35> cu;
      if (in[2]=="-infty") {
        lower_lim=-std::numeric_limits<cpp_dec_float_35>::infinity();
      } else {
        function_to_fp_nothrow(in[2],lower_lim,cu);
      }
      if (in[3]=="infty") {
        upper_lim=std::numeric_limits<cpp_dec_float_35>::infinity();
      } else {
        function_to_fp_nothrow(in[3],upper_lim,cu);
      }
      int retx;
      if (method=="kb") {
        retx=ikb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else if (method=="deb") {
        retx=ideb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else {
        retx=iac.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      }
        
      if (retx!=0) {
        cerr << "Integrating " << func << " failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_35): ";
      cout << dtos(d,precision) << endl;
      return 0;
      
    } else if (precision>17) {
      
      cpp_dec_float_25 d=0, err, lower_lim, upper_lim;
      convert_units<cpp_dec_float_25> cu;
      if (in[2]=="-infty") {
        lower_lim=-std::numeric_limits<cpp_dec_float_25>::infinity();
      } else {
        function_to_fp_nothrow(in[2],lower_lim,cu);
      }
      if (in[3]=="infty") {
        upper_lim=std::numeric_limits<cpp_dec_float_25>::infinity();
      } else {
        function_to_fp_nothrow(in[3],upper_lim,cu);
      }
      int retx;
      if (method=="kb") {
        retx=ikb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else if (method=="deb") {
        retx=ideb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else {
        retx=iac.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      }
        
      if (retx!=0) {
        cerr << "Integrating " << func << " failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_25): ";
      cout << dtos(d,precision) << endl;
      
      return 0;
      
    } else if (precision>14) {
      
      long double d=0, err, lower_lim, upper_lim;
      convert_units<long double> cu;
      if (in[2]=="-infty") {
        lower_lim=-std::numeric_limits<long double>::infinity();
      } else {
        function_to_fp_nothrow(in[2],lower_lim,cu);
      }
      if (in[3]=="infty") {
        upper_lim=std::numeric_limits<long double>::infinity();
      } else {
        function_to_fp_nothrow(in[3],upper_lim,cu);
      }
      int retx;
      if (method=="kb") {
        retx=ikb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else if (method=="deb") {
        retx=ideb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else {
        retx=iac.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      }
        
      if (retx!=0) {
        cerr << "Integrating " << func << " failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (long double): ";
      cout << dtos(d,precision) << endl;
      
      return 0;
      
    } else {
      
      double d=0, err, lower_lim, upper_lim;
      convert_units<double> cu;
      if (in[2]=="-infty") {
        lower_lim=-std::numeric_limits<double>::infinity();
      } else {
        function_to_fp_nothrow(in[2],lower_lim,cu);
      }
      if (in[3]=="infty") {
        upper_lim=std::numeric_limits<double>::infinity();
      } else {
        function_to_fp_nothrow(in[3],upper_lim,cu);
      }
      int retx;
      if (method=="kb") {
        retx=ikb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else if (method=="deb") {
        retx=ideb.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      } else {
        retx=iac.integ_err_multip([fmsp](auto &&t) mutable
        { return (*fmsp)(t); },lower_lim,upper_lim,d,err);
      }
        
      if (retx!=0) {
        cerr << "Integrating " << func << " failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (double): ";
      cout << dtos(d,precision) << endl;
      
    }

#endif
    
  } else {
    
    // Normal double-precision integration

    if (precision>16) {
      std::cerr << "Warning: multiprecision is required to numerically "
                << "integrate to the\n requested precision."
                << std::endl;
    }
    
    double d=0, err, lower_lim, upper_lim;
    convert_units<double> cu;
    if (in[2]=="-infty") {
      lower_lim=-std::numeric_limits<double>::infinity();
    } else {
      function_to_fp_nothrow(in[2],lower_lim,cu);
    }
    if (in[3]=="infty") {
      upper_lim=std::numeric_limits<double>::infinity();
    } else {
      function_to_fp_nothrow(in[3],upper_lim,cu);
    }
    funct_string fs(func,var);
    funct f=std::bind(std::mem_fn<double(double) const>
                      (&funct_string::operator()),&fs,
                      std::placeholders::_1);
    int retx;
    if (method=="kb") {
      retx=ikb.integ_err(f,lower_lim,upper_lim,d,err);
    } else if (method=="deb") {
      retx=ideb.integ_err(f,lower_lim,upper_lim,d,err);
    } else {
      retx=iac.integ_err(f,lower_lim,upper_lim,d,err);
    }
    if (retx!=0) {
      cerr << "Integrating " << func << " failed." << endl;
      return 1;
    }
    
    if (scientific) cout.setf(ios::scientific);
    else cout.unsetf(ios::scientific);
    cout.precision(precision);
    if (verbose>0) cout << "Result: ";
    cout << d << " ± " << err << endl;
    std::string us;
    if (verbose>1) {
      us=unc_to_string(d,err,1);
    } else {
      us=unc_to_string(d,err);
    }
    cout << us << endl;

  }

  return 0;
}

int acol_manager::comm_nlines(std::vector<std::string> &sv, 
			      bool itive_com) {
  if (type!="table") {
    cerr << "No table in 'nlines'." << endl;
    return 1;
  }

  table_obj.add_constant("nlines",table_obj.get_nlines());
  cout << "The table has " << table_obj.get_nlines() << " lines." << endl;
  
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
  fout->precision(precision);

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
	col_wids[i]=precision+6;
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
	    int nsp=((int)(table_obj.get_column_name(j).size()-precision-6));
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
      vector<string> svx, sv_out;
      for(size_t k=0;k<intv_obj.size();k++) {
	svx.push_back(o2scl::itos(intv_obj[k])+' ');
      }
      screenify(intv_obj.size(),svx,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),intv_obj,true);
    }

  } else if (type=="double[]") {

    if (pretty) {
      vector<string> svx, sv_out;
      for(size_t k=0;k<doublev_obj.size();k++) {
	if (has_minus_sign(&doublev_obj[k])) {
	  svx.push_back(o2scl::dtos(doublev_obj[k])+' ');
	} else {
	  svx.push_back(" "+o2scl::dtos(doublev_obj[k])+' ');
	}
      }
      screenify(doublev_obj.size(),svx,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),doublev_obj,true);
    }

  } else if (type=="size_t[]") {

    if (pretty) {
      vector<string> svx, sv_out;
      for(size_t k=0;k<size_tv_obj.size();k++) {
	svx.push_back(o2scl::szttos(size_tv_obj[k])+' ');
      }
      screenify(size_tv_obj.size(),svx,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	(*fout) << sv_out[k] << endl;
      }
    } else {
      vector_out((*fout),size_tv_obj,true);
    }
    
  } else if (type=="string[]") {
    
    for(size_t k=0;k<stringv_obj.size();k++) {
      if (stringv_obj[k].length()==0) {
        (*fout) << "(empty string)" << endl;
      } else {
        (*fout) << stringv_obj[k] << endl;
      }
    }
    (*fout) << endl;

  } else if (type=="vec_vec_string") {
    
    (*fout) << vvstring_obj.size() << endl;
    for(size_t k=0;k<vvstring_obj.size();k++) {
      (*fout) << vvstring_obj[k].size() << endl;
      for(size_t kk=0;kk<vvstring_obj[k].size();kk++) {
        (*fout) << '\"' << vvstring_obj[k][kk] << '\"' << endl;
      }
    }

  } else if (type=="vec_vec_double") {
    
    (*fout) << vvdouble_obj.size() << endl;
    for(size_t k=0;k<vvdouble_obj.size();k++) {
      for(size_t kk=0;kk<vvdouble_obj[k].size();kk++) {
        (*fout) << vvdouble_obj[k][kk] << " ";
      }
      (*fout) << endl;
    }

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

  } else if (type=="prob_dens_mdim_gmm") {

    pgmm_obj.write_generic(*fout);

  } else if (type=="tensor") {

    tensor_out(*fout,tensor_obj,pretty);
    
  } else if (type=="tensor_grid") {

    tensor_grid_out(*fout,tensor_grid_obj,pretty);
    
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

