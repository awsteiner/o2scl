/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>
#include <o2scl/xml.h>
#include <o2scl/gmm_python.h>

#include <o2scl/set_mpfr.h>

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
int acol_manager::comm_to_gaussian
int acol_manager::comm_to_hist
int acol_manager::comm_to_hist_2d
int acol_manager::comm_to_table
int acol_manager::comm_to_table3d
int acol_manager::comm_to_table3d_sum
int acol_manager::comm_to_tensor
int acol_manager::comm_to_tensor_grid
int acol_manager::comm_to_tg_fermi
int acol_manager::comm_type
int acol_manager::comm_value
int acol_manager::comm_version
int acol_manager::comm_wdocs
int acol_manager::comm_wstats
void acol_manager::xml_replacements
int acol_manager::comm_xml_to_o2
int acol_manager::comm_x_name
int acol_manager::comm_y_name
 */

int acol_manager::comm_thin_mcmc(std::vector<std::string> &sv,
                                 bool itive_com) {
  if (type=="table") {
    if (sv.size()<2) {
      cerr << "Not enough arguments in command 'thin-mcmc'." << endl;
      return 2;
    }
    size_t window=stoszt(sv[1]);
    std::string mult_col;
    if (sv.size()>=3) {
      mult_col=sv[2];
    }
    size_t running_sum=0;
    size_t count=0;
    table_units<> tnew;

    for(size_t i=0;i<table_obj.get_nconsts();i++) {
      string tnam;
      double tval;
      table_obj.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      tnew.add_constant(tnam,tval);
    }
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      std::string col=table_obj.get_column_name(i);
      std::string unit=table_obj.get_unit(col);
      tnew.new_column(col);
      tnew.set_unit(col,unit);
    }
    
    for(size_t j=0;j<table_obj.get_nlines();j++) {
      if (mult_col.length()==0 || ((size_t)table_obj.get(mult_col,j))>0) {
        while (count<=running_sum) {
          tnew.copy_row(table_obj,j);
          count+=window;
        }
        if (mult_col.length()>0) {
          running_sum+=((size_t)(table_obj.get(mult_col,j)));
        } else {
          running_sum++;
        }
      }
    }
    table_obj=tnew;
  } else {
    cerr << "Command 'thin-mcmc' not supported for objects of "
         << "type " << type << endl;
    return 1;
  }
  
  return 0;
}

int acol_manager::comm_to_gaussian(std::vector<std::string> &sv,
                                   bool itive_com) {
  if (type=="table") {

    if (sv.size()<3) {
      cerr << "Not enough arguments for to-gaussian." << endl;
    }

    int n_dim=sv.size()-1;
    vector<string> col_names;
    cout << "X columns: " << endl;
    for(size_t i=1;i<sv.size();i++) {
      col_names.push_back(sv[i]);
      cout << i-1 << ": " << sv[i] << endl;
    }

    matrix_view_table<> mvt(table_obj,col_names);
    
    pdmg_obj.set<matrix_view_table<> >(n_dim,table_obj.get_nlines(),mvt);
    
    command_del(type);
    clear_obj();
    command_add("prob_dens_mdim_gaussian");
    type="prob_dens_mdim_gaussian";

    /*

      This doesn't work because we need to create a 
      matrix_view_vec_vec temporary which can't go out 
      of scope
      
      } else if (type=="vec_vec_double") {
      
      if (sv.size()<3) {
      cerr << "Not enough arguments for to-gaussian." << endl;
      }
      
      int n_dim=sv.size()-1;
      vector<size_t> vec_list;
      string_to_uint_list(sv[1],vec_list);
      cout << "X columns: ";
      vector_out(cout,vec_list,true);
      
      vector<vector<double>> vvd2;
      for(size_t i=0;i<vec_list.size();i++) {
      vvd2.push_back(vvdouble_obj[vec_list[i]]);
      }
      
      pdmg_obj.set<matrix_view_vec_vec>(vvd2.size(),vvd2[0].size(),vvd2);
      
      command_del(type);
      clear_obj();
      command_add("prob_dens_mdim_gaussian");
      type="prob_dens_mdim_gaussian";

    */
    
  }
    
  return 0;
}

int acol_manager::comm_to_gmm(std::vector<std::string> &sv,
                                bool itive_com) {
  if (type=="table") {

#ifdef O2SCL_PYTHON
    if (sv.size()<3) {
      cerr << "Not enough arguments for to-gmm." << endl;
    }

    size_t n_gauss=o2scl::stoszt(sv[1]);
    vector<string> col_names;
    cout << "X columns: " << endl;
    for(size_t i=2;i<sv.size();i++) {
      col_names.push_back(sv[i]);
      cout << i-2 << ": " << sv[i] << endl;
    }
    
    pgmm_obj.verbose=verbose;

    // Copy the table data to a tensor for use in gmm_python
    tensor<> tin;
    vector<size_t> in_size={table_obj.get_nlines(),col_names.size()};
    tin.resize(2,in_size);

    for(size_t i=0;i<table_obj.get_nlines();i++) {
      for(size_t j=0;j<col_names.size();j++) {
        vector<size_t> ix;
        ix={i,j};
        tin.get(ix)=table_obj.get(col_names[j],i);
      }
    }
    
    gmm_python gp("o2sclpy",n_gauss,tin,
                  ((string)"verbose=")+o2scl::itos(verbose),
                  "gmm_sklearn",verbose);
                  
    gp.get_python();

    pgmm_obj=gp.get_gmm();
    
    command_del(type);
    clear_obj();
    command_add("prob_dens_mdim_gmm");
    type="prob_dens_mdim_gmm";

#else
      cerr << "Python support not included." << endl;
      return 3;
#endif
    
  }
    
  return 0;
}

int acol_manager::comm_to_kde(std::vector<std::string> &sv,
                              bool itive_com) {
  if (type=="table") {

#ifdef O2SCL_PYTHON
    
    if (sv.size()<3) {
      cerr << "Not enough arguments for to-kde." << endl;
    }

    string options=sv[1];
    
    vector<string> col_names;
    cout << "Columns used for KDE: " << endl;
    for(size_t i=2;i<sv.size();i++) {
      col_names.push_back(sv[i]);
      cout << i-2 << ": " << sv[i] << endl;
    }

    int kde_verbose=verbose-1;
    if (kde_verbose<0) kde_verbose=0;
    pkde_obj.verbose=kde_verbose;
    kwargs kw;
    if (options!="none" && options!="None") {
      kw.set(options);
    }
    
    // Copy the table data to a tensor for use in kde_python
    tensor<> ttemp;
    vector<size_t> in_size={table_obj.get_nlines(),col_names.size()};
    ttemp.resize(2,in_size);

    for(size_t i=0;i<table_obj.get_nlines();i++) {
      for(size_t j=0;j<col_names.size();j++) {
        vector<size_t> ix;
        ix={i,j};
        ttemp.get(ix)=table_obj.get(col_names[j],i);
      }
    }

    if (kw.get_string("method")=="scipy") {
      vector<double> weights;
      if (kw.is_set("weights")) {
        std::string wcol=kw.get_string("weights");
        cout << "Using weights from column: " << wcol << endl;
        for(size_t i=0;i<table_obj.get_nlines();i++) {
          weights.push_back(table_obj.get(wcol,i));
        }
      }
      cout << "Herez: " << col_names.size() << " "
           << table_obj.get_nlines() << endl;
      pkde_obj.set_function("o2sclpy",ttemp,
                            weights,((string)"verbose=")+
                            o2scl::itos(kde_verbose),
                            "kde_scipy",kde_verbose);
    } else {
      uniform_grid_log_end<double> ug(1.0e-3,1.0e3,99);
      vector<double> bw_array;
      ug.vector(bw_array);
      pkde_obj.set_function("o2sclpy",ttemp,
                            bw_array,((string)"verbose=")+
                            o2scl::itos(kde_verbose),"kde_sklearn",
                            kde_verbose);
    }
    
    command_del(type);
    clear_obj();
    command_add("prob_dens_mdim_kde");
    type="prob_dens_mdim_kde";

#else
      cerr << "Python support not included." << endl;
      return 3;
#endif    
    
  }
    
  return 0;
}

int acol_manager::comm_to_pdma(std::vector<std::string> &sv,
                               bool itive_com) {
  if (type=="table") {

    if (sv.size()<3) {
      cerr << "Not enough arguments for to-pdma." << endl;
    }

    int n_dim=sv.size()-1;
    vector<string> col_names;
    cout << "X columns: " << endl;
    vector<double> max, min;
    for(size_t i=1;i<sv.size();i++) {
      col_names.push_back(sv[i]);
      if (i==sv.size()-1) {
        cout << "Y column: " << endl;
        cout << sv[i] << endl;
      } else {
        cout << i-1 << ": " << sv[i] << endl;
      }
      if (i!=sv.size()-1) {
        max.push_back(table_obj.max(sv[i]));
        min.push_back(table_obj.min(sv[i]));
      }
    }

    const_matrix_view_table<> mvt(table_obj,col_names);

    pdma_obj.set(min,max);
    pdma_obj.initial_parse(mvt);
    
    command_del(type);
    clear_obj();
    command_add("prob_dens_mdim_amr");
    type="prob_dens_mdim_amr";
    
  }
    
  return 0;
}

int acol_manager::comm_to_hist(std::vector<std::string> &sv, 
			       bool itive_com) {

  if (type=="table") {

    // AWS, 11/22/22, This looks like it could be old, so I'm removing
    //it for now.
    // std::string i1;
    //int ret=get_input_one(sv,((string)"Enter \"2d\" for
    //2d histogram ")+ +"and \"1d\" for 1d histogram",i1,"to-hist",
    //itive_com); if (ret!=0) return ret;
    
    vector<string> in, pr;
    pr.push_back("Column name");
    pr.push_back("Number of bins");
    int ret=get_input(sv,pr,in,"to-hist",itive_com);
    if (ret!=0) return ret;
      
    std::string col2;
    if (sv.size()>3) {
      col2=sv[3];
    } else if (itive_com) {
      col2=cl->cli_gets("Column for weights (or blank for none): ");
    }

    size_t nbins;
    int sret=o2scl::stoszt_nothrow(in[1],nbins);
    if (sret!=0 || nbins==0) {
      cerr << "Failed to interpret " << in[1]
	   << " as a positive number of bins." << endl;
      return exc_einval;
    }
    if (col2.length()==0) {
      hist_obj.from_table(table_obj,in[0],nbins);
    } else {
      hist_obj.from_table(table_obj,in[0],col2,nbins);
    }

    command_del(type);
    clear_obj();
    command_add("hist");
    type="hist";

    return 0;
  } 

  if (type=="table3d") {

    vector<string> in, pr;
    pr.push_back("Slice name");
    pr.push_back("Number of bins");
    int ret=get_input(sv,pr,in,"to-hist",itive_com);
    if (ret!=0) return ret;

    size_t ix;
    if (!table3d_obj.is_slice(in[0],ix)) {
      cerr << "No slice named " << in[0] << " in table3d." << endl;
      return 11;
    }

    hist_obj=table3d_obj.to_hist(in[0],o2scl::stoszt(in[1]),verbose);
    
    command_del(type);
    clear_obj();
    command_add("hist");
    type="hist";

    return 0;
  } 

  cerr << "Cannot convert object of type " << type << " to histogram."
       << endl;
  
  return 1;
}

int acol_manager::comm_to_hist_2d(std::vector<std::string> &sv, 
				  bool itive_com) {

  if (type=="table") {
    
    vector<string> in, pr;
    pr.push_back("Column name for x-axis");
    pr.push_back("Column name for y-axis");
    pr.push_back("Number of bins in x direction");
    pr.push_back("Number of bins in y direction");
    int ret=get_input(sv,pr,in,"to-hist-2d",itive_com);
    if (ret!=0) return ret;
    
    std::string col2;
    if (sv.size()>5) {
      col2=sv[5];
      // We don't want to prompt for weights if the user has given
      // enough arguments to proceed, so we test for sv.size()<5
    } else if (itive_com && sv.size()<5) {
      col2=cl->cli_gets("Column for weights (or blank for none): ");
    }
    
    size_t nbinsx, nbinsy;
    int sret=o2scl::stoszt_nothrow(in[2],nbinsx);
    if (sret!=0 || nbinsx==0) {
      cerr << "Failed to interpret " << in[2]
	   << " as a positive number of bins." << endl;
      return exc_einval;
    }
    sret=o2scl::stoszt_nothrow(in[3],nbinsy);
    if (sret!=0 || nbinsy==0) {
      cerr << "Failed to interpret " << in[2]
	   << " as a positive number of bins." << endl;
      return exc_einval;
    }
    
    if (col2.length()==0) {
      hist_2d_obj.from_table(table_obj,in[0],in[1],nbinsx,nbinsy);
    } else {
      hist_2d_obj.from_table(table_obj,in[0],in[1],col2,
			     nbinsx,nbinsy);
    }

    command_del(type);
    clear_obj();
    command_add("hist_2d");
    type="hist_2d";
    
    return 0;
    
  } else if (type=="table3d") {

    std::string i1;
    int ret=get_input_one(sv,"Enter slice name",i1,"to-hist-2d",itive_com);
    if (ret!=0) return ret;

    hist_2d_obj=table3d_obj.to_hist_2d(i1);
    
    command_del(type);
    clear_obj();
    command_add("hist_2d");
    type="hist_2d";

    return 0;
  } 

  cerr << "Cannot convert object of type " << type << " to histogram."
       << endl;
  
  return 1;
}

int acol_manager::comm_to_table(std::vector<std::string> &sv, bool itive_com) {

  if (type=="double[]") {

    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(doublev_obj.size());
    for(size_t i=0;i<doublev_obj.size();i++) {
      table_obj.set(i1,i,doublev_obj[i]);
    }

    command_del(type);
    clear_obj();
    command_add("table");
    type="table";
    
  } else if (type=="int[]") {

    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(intv_obj.size());
    for(size_t i=0;i<intv_obj.size();i++) {
      table_obj.set(i1,i,intv_obj[i]);
    }

    command_del(type);
    clear_obj();
    command_add("table");
    type="table";

  } else if (type=="size_t[]") {
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(size_tv_obj.size());
    for(size_t i=0;i<size_tv_obj.size();i++) {
      table_obj.set(i1,i,size_tv_obj[i]);
    }

    command_del(type);
    clear_obj();
    command_add("table");
    type="table";
    
  } else if (type=="hist") {
    
    table_obj.clear();
    hist_obj.copy_to_table(table_obj,"rep","low","high","wgt");
    command_del(type);
    clear_obj();
    command_add("table");
    type="table";
    
  } else if (type=="prob_dens_mdim_kde") {

#ifdef O2SCL_PYTHON
    
    if (pkde_obj.dim()!=1) {
      cerr << "Command to-table only works on a 1-dimensional KDE" << endl;
      return 1;
    }
    
    table_obj.clear();
    
    const o2scl::tensor<> t=pkde_obj.get_data();

    // Get x_min and x_max
    double x_min, x_max;
    vector<size_t> ix(1);
    
    vector<double> v;
    for(size_t i=0;i<t.get_size(0);i++) {
      ix[0]=i;
      v.push_back(t.get(ix));
    }
    vector_minmax_value(v.size(),v,x_min,x_max);

    if (verbose>2) {
      cout << "  to-table: x_min,x_max: " << x_min << " " << x_max << endl;
    }

    // Approximate y_min and y_max
    cout << pkde_obj.dim() << endl;
    vector<double> x_arr(1);
    x_arr[0]=x_min;
    double y_min=pkde_obj.pdf(x_arr);
    x_arr[0]=x_max;
    double y_max=pkde_obj.pdf(x_arr);
    double dx=(x_max-x_min)/10;
    x_arr[0]=x_min+dx;
    for(size_t i=0;i<10;i++) {
      double y=pkde_obj.pdf(x_arr);
      if (y<y_min) y_min=y;
      if (y>y_max) y_max=y;
      x_arr[0]+=dx;
    }

    if (verbose>2) {
      cout << "  to-table: y_min,y_max: " << y_min << " " << y_max << endl;
    }
    
    // Determine the x and y-values at the boundaries
    double x_left=x_min, x_right=x_max;
    x_arr[0]=x_left;
    double y_left=pkde_obj.pdf(x_arr);
    x_arr[0]=x_right;
    double y_right=pkde_obj.pdf(x_arr);

    if (verbose>2) {
      cout << "  to-table: x_left,x_right: " << x_left << " " << x_right
           << endl;
    }
    
    // There may be significant probability at the boundaries, so
    // try to extend them outwards if needed to capture the
    // full distribution
    size_t j=0;
    while (y_left>y_max/1.0e3 && j<10) {
      x_left-=dx;
      x_arr[0]=x_left;
      y_left=pkde_obj.pdf(x_arr);
      std::cout << y_left << " " << j << std::endl;
      j++;
    }
    j=0;
    while (y_right>y_max/1.0e3 && j<10) {
      x_right+=dx;
      x_arr[0]=x_right;
      y_right=pkde_obj.pdf(x_arr);
      std::cout << y_right << " " << j << std::endl;
      j++;
    }

    if (verbose>2) {
      cout << "  to-table: x_left,x_right: " << x_left << " "
           << x_right << endl;
      cout << "  to-table: y_left,y_right: " << y_left << " "
           << y_right << endl;
    }
    
    // Now create the table based on a uniform grid
    table_obj.line_of_names("x y");
    
    dx=(x_right-x_left)/200;
    for(x_arr[0]=x_left;x_arr[0]<x_right+dx/2.0;x_arr[0]+=dx) {
      double line[2]={x_arr[0],pkde_obj.pdf(x_arr)};
      table_obj.line_of_data(2,line);
    }
    
    command_del(type);
    clear_obj();
    command_add("table");
    type="table";

#else
      cerr << "Python support not included." << endl;
      return 3;
#endif
    
  } else if (type=="tensor_grid") {
    
    size_t rank=tensor_grid_obj.get_rank();
    
    vector<string> in, pr;
    pr.push_back("Index to vary");
    pr.push_back("Grid name");
    pr.push_back("Data name");
    int ret=get_input(sv,pr,in,"to-table",itive_com);
    if (ret!=0) return ret;

    size_t ix=o2scl::stoszt(in[0]);
    if (ix>=rank) {
      cerr << "Index larger than rank." << endl;
      return 1;
    }

    for(size_t i=0;i<3;i++) {
      std::vector<std::string>::iterator it=sv.begin();
      it++;
      sv.erase(it);
    }
    
    vector<string> in2, pr2;
    for(size_t i=0;i<rank;i++) {
      if (i!=ix) {
	pr2.push_back(((std::string)"Value for index ")+o2scl::szttos(i));
      }
    }
    int ret2=get_input(sv,pr2,in2,"to-table",itive_com);
    if (ret2!=0) return ret2;

    vector<double> values(rank);
    size_t i2=0;
    for(size_t i=0;i<rank;i++) {
      if (i!=ix) {
	values[i]=o2scl::stod(in2[i2]);
	if (verbose>0) {
	  cout << "Fixing value for index " << i << " to " << in2[i2] << endl;
	}
	i2++;
      }
    }

    if (verbose>0) {
      cout << "Index " << ix << " is free. "
	   << "New columns are: " << in[1] << " and " << in[2] << endl;
    }

    table_obj.clear();
    table_obj.new_column(in[1]);
    table_obj.new_column(in[2]);
    for(size_t i=0;i<tensor_grid_obj.get_size(ix);i++) {
      values[ix]=tensor_grid_obj.get_grid(ix,i);
      double line[2]={values[ix],tensor_grid_obj.interp_linear(values)};
      table_obj.line_of_data(2,line);
    }

    command_del(type);
    clear_obj();
    command_add("table");
    type="table";
    
  } else if (type=="table3d") {

    table_obj.clear();
    table_obj.new_column(table3d_obj.get_x_name());
    table_obj.new_column(table3d_obj.get_y_name());
    for(size_t j=0;j<table3d_obj.get_nslices();j++) {
      table_obj.new_column(table3d_obj.get_slice_name(j));
    }

    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
        vector<double> line;
        line.push_back(table3d_obj.get_grid_x(i));
        line.push_back(table3d_obj.get_grid_y(j));
        for(size_t k=0;k<table3d_obj.get_nslices();k++) {
          line.push_back(table3d_obj.get(table3d_obj.get_grid_x(i),
                                         table3d_obj.get_grid_y(j),
                                         table3d_obj.get_slice_name(k)));
        }
        table_obj.line_of_data(line.size(),line);
      }
    }
    
    command_del(type);
    clear_obj();
    command_add("table");
    type="table";
    
  }
  
  return 0;
}

int acol_manager::comm_to_table3d(std::vector<std::string> &sv,
				  bool itive_com) {

  if (type=="table") {

    vector<string> in, pr;
    pr.push_back("Column for x grid");
    pr.push_back("Column for y grid");

    int ret=get_input(sv,pr,in,"to-table3d",itive_com);
    if (ret!=0) return ret;

    std::string xname=in[0];
    std::string yname=in[1];

    double empty_value=0.0;
    double eps=1.0e-12;
    if (in.size()>2) {
      empty_value=o2scl::stod(in[2]);
    }
    if (in.size()>3) {
      eps=o2scl::stod(in[3]);
    }
    
    table3d_obj.clear();
    int cret=table3d_obj.read_table(table_obj,xname,yname,empty_value,
				    verbose,false,eps);
    if (cret!=0) {
      cerr << "Convert 'table' to 'table3d' failed." << endl;
      return 1;
    }

    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";
    
  } else if (type=="hist_2d") {
    
    table3d_obj.clear();

    vector<string> in, pr;
    pr.push_back("Name for x grid");
    pr.push_back("Name for y grid");
    pr.push_back("Name for weights");

    int ret=get_input(sv,pr,in,"to-table3d",itive_com);
    if (ret!=0) return ret;

    std::string xname=in[0];
    std::string yname=in[1];
    std::string wname=in[2];

    hist_2d_obj.copy_to_table3d(table3d_obj,xname,yname,wname);
    
    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";
    
  } else if (type=="prob_dens_mdim_kde") {
    
#ifdef O2SCL_PYTHON

    if (pkde_obj.dim()!=2) {
      cerr << "Command to-table only works on a 2-dimensional KDE" << endl;
      return 1;
    }
    
    table3d_obj.clear();

    const o2scl::tensor<> t=pkde_obj.get_data();

    // Get x_min and x_max for both coordinates
    double x_min[2], x_max[2];
    size_t N=t.get_size(0);
    vector<size_t> ix;
    ix[0]=0;
    ix[1]=0;
    x_min[0]=t.get(ix);
    x_max[0]=t.get(ix);
    ix[1]=1;
    x_min[1]=t.get(ix);
    x_max[1]=t.get(ix);
    for(size_t i=1;i<N;i++) {
      ix[0]=i;
      ix[1]=0;
      double val=t.get(ix);
      if (val<x_min[0]) x_min[0]=val;
      if (val>x_max[0]) x_max[0]=val;
      ix[1]=1;
      val=t.get(ix);
      if (val<x_min[1]) x_min[1]=val;
      if (val>x_max[1]) x_max[1]=val;
    }

    uniform_grid_end<double> ugx(x_min[0],x_max[0],49);
    uniform_grid_end<double> ugy(x_min[1],x_max[1],49);

    table3d_obj.set_xy("x",ugx,"y",ugy);
    table3d_obj.new_slice("z");
    
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
        vector<double> xx(2);
        xx[0]=table3d_obj.get_grid_x(i);
        xx[1]=table3d_obj.get_grid_y(j);
        table3d_obj.set(i,j,"z",pkde_obj.pdf(xx));
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
    
  } else if (type=="tensor" || type=="tensor<size_t>" || type=="tensor<int>") {

    size_t rank=tensor_obj.get_rank();

    vector<string> in, pr;
    pr.push_back("First index to vary");
    pr.push_back("Second index to vary");
    pr.push_back("Slice name");

    int ret=get_input(sv,pr,in,"to-table3d",itive_com);
    if (ret!=0) return ret;

    size_t ix_x=o2scl::stoszt(in[0]);
    size_t ix_y=o2scl::stoszt(in[1]);
    if (ix_x>=rank || ix_y>=rank) {
      cerr << "Index larger than rank." << endl;
      return 1;
    }
    
    vector<string> in2, pr2;
    if (rank>2) {
      for(size_t i=0;i<3;i++) {
	std::vector<std::string>::iterator it=sv.begin();
	it++;
	sv.erase(it);
      }
      
      for(size_t i=0;i<rank;i++) {
	if (i!=ix_x && i!=ix_y) {
	  pr2.push_back(((std::string)"Fixed index for rank ")+
			o2scl::szttos(i));
	}
      }
      
      int ret2=get_input(sv,pr2,in2,"to-table3d",itive_com);
      if (ret!=0) return ret2;
    }
    
    uniform_grid<double> ugx, ugy;
    if (type=="tensor") {
      ugx=uniform_grid_end<double>(0,tensor_obj.get_size(ix_x)-1,
				   tensor_obj.get_size(ix_x)-1);
      ugy=uniform_grid_end<double>(0,tensor_obj.get_size(ix_y)-1,
				   tensor_obj.get_size(ix_y)-1);
    } else if (type=="tensor<int>") {
      ugx=uniform_grid_end<double>(0,tensor_int_obj.get_size(ix_x)-1,
				   tensor_int_obj.get_size(ix_x)-1);
      ugy=uniform_grid_end<double>(0,tensor_int_obj.get_size(ix_y)-1,
				   tensor_int_obj.get_size(ix_y)-1);
    } else {
      ugx=uniform_grid_end<double>(0,tensor_size_t_obj.get_size(ix_x)-1,
				   tensor_size_t_obj.get_size(ix_x)-1);
      ugy=uniform_grid_end<double>(0,tensor_size_t_obj.get_size(ix_y)-1,
				   tensor_size_t_obj.get_size(ix_y)-1);
    }
    table3d_obj.clear();
    table3d_obj.set_xy("x",ugx,"y",ugy);
    table3d_obj.new_slice(in[2]);
    vector<size_t> ix(rank);
    size_t j=0;
    for(size_t i=0;i<rank;i++) {
      if (i!=ix_x && i!=ix_y) {
	ix[i]=o2scl::stoszt(in2[j]);
	j++;
      }
    }
    if (type=="tensor") {
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t jx=0;jx<table3d_obj.get_ny();jx++) {
	ix[ix_x]=i;
	ix[ix_y]=jx;
	table3d_obj.set(i,jx,in[2],tensor_obj.get(ix));
      }
    }
    } else if (type=="tensor<int>") {
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t jx=0;jx<table3d_obj.get_ny();jx++) {
	ix[ix_x]=i;
	ix[ix_y]=jx;
	table3d_obj.set(i,jx,in[2],tensor_int_obj.get(ix));
      }
    }
    } else {
      for(size_t i=0;i<table3d_obj.get_nx();i++) {
	for(size_t jx=0;jx<table3d_obj.get_ny();jx++) {
	  ix[ix_x]=i;
	  ix[ix_y]=jx;
	  table3d_obj.set(i,jx,in[2],tensor_size_t_obj.get(ix));
	}
      }
    }

    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";

  } else if (type=="tensor_grid") {

    size_t rank=tensor_grid_obj.get_rank();

    vector<string> in, pr;
    pr.push_back("First index to vary");
    pr.push_back("Second index to vary");
    pr.push_back("Slice name");
    int ret=get_input(sv,pr,in,"to-table3d",itive_com);
    if (ret!=0) return ret;

    size_t ix_x=o2scl::stoszt(in[0]);
    size_t ix_y=o2scl::stoszt(in[1]);
    if (ix_x>=rank || ix_y>=rank) {
      cerr << "Index larger than rank." << endl;
      return 1;
    }

    for(size_t i=0;i<3;i++) {
      std::vector<std::string>::iterator it=sv.begin();
      it++;
      sv.erase(it);
    }
    
    vector<string> in2, pr2;
    for(size_t i=0;i<rank;i++) {
      if (i!=ix_x && i!=ix_y) {
	pr2.push_back(((std::string)"Value for index ")+o2scl::szttos(i));
      }
    }
    int ret2=get_input(sv,pr2,in2,"to-table3d",itive_com);
    if (ret2!=0) return ret2;

    vector<double> values(rank);
    size_t i2=0;
    for(size_t i=0;i<rank;i++) {
      if (i!=ix_x && i!=ix_y) {
	if (verbose>0) {
	  cout << "Fixing value for index " << i << " to " << in2[i2] << endl;
	}
	values[i]=o2scl::function_to_double(in2[i2]);
	i2++;
      }
    }

    if (verbose>0) {
      cout << "Indices " << ix_x << " and " << ix_y << " are free. "
	   << "New slice name is: " << in[2] << endl;
    }

    table3d_obj.clear();
    tensor_grid_obj.copy_table3d_interp_values_setxy<vector<double> >
      (ix_x,ix_y,values,table3d_obj,"x","y",in[2]);
    
    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";

  } else if (type=="prob_dens_mdim_amr") {

    vector<string> in, pr;
    pr.push_back("First index to vary");
    pr.push_back("Second index to vary");
    pr.push_back("Number of x grid points");
    pr.push_back("Number of y grid points");
    pr.push_back("Slice name");
    int ret=get_input(sv,pr,in,"to-table3d",itive_com);
    if (ret!=0) return ret;
    
    table3d_obj.clear();
    size_t i, j, ni, nj;
    o2scl::stoszt_nothrow(in[0],i);
    o2scl::stoszt_nothrow(in[1],j);
    o2scl::stoszt_nothrow(in[2],ni);
    o2scl::stoszt_nothrow(in[3],nj);
    
    table3d_obj.set_xy("x",uniform_grid_end<double>(pdma_obj.low[i],
						    pdma_obj.high[i],ni-1),
		       "y",uniform_grid_end<double>(pdma_obj.low[j],
						    pdma_obj.high[j],nj-1));
    
    cout << "Converting pdma to table3d, using index "
	 << in[0] << ", index " << in[1] << ", and slice "
	 << in[4] << endl;
    pdma_obj.two_indices_to_density(i,j,table3d_obj,in[4]);
    
    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";

  } else {
    
    cerr << "Cannot use command 'to-table3d' for type "
	 << type << "." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_to_table3d_sum(std::vector<std::string> &sv,
				      bool itive_com) {

  if (type=="tensor") {

    size_t rank=tensor_obj.get_rank();

    vector<string> in, pr;
    pr.push_back("First index name");
    pr.push_back("First index to vary");
    pr.push_back("Second index name");
    pr.push_back("Second index to vary");
    pr.push_back("Slice name");

    int ret=get_input(sv,pr,in,"to-table3d-sum",itive_com);
    if (ret!=0) return ret;

    size_t ix_x=o2scl::stoszt(in[1]);
    size_t ix_y=o2scl::stoszt(in[3]);
    if (ix_x>=rank || ix_y>=rank) {
      cerr << "Index larger than rank." << endl;
      return 1;
    }

    vector<string> in2, pr2;
    if (rank>2) {
      for(size_t i=0;i<5;i++) {
	std::vector<std::string>::iterator it=sv.begin();
	it++;
	sv.erase(it);
      }
      
      for(size_t i=0;i<rank;i++) {
	if (i!=ix_x && i!=ix_y) {
	  pr2.push_back(((std::string)"Fixed index for rank ")+
			o2scl::szttos(i));
	}
      }
      
      int ret2=get_input(sv,pr,in,"to-table3d-sum",itive_com);
      if (ret!=0) return ret2;
    }

    table3d_obj.clear();
    tensor_obj.copy_table3d_sum(ix_x,ix_y,table3d_obj,in[0],
				   in[2],in[4]);

    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";

  } else {
    
    cerr << "Cannot use command 'to-table3d-sum' for type "
	 << type << "." << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_to_tensor(std::vector<std::string> &sv,
				 bool itive_com) {
  
  if (type=="tensor_grid") {

    tensor_obj=tensor_grid_obj;

    command_del(type);
    clear_obj();
    command_add("tensor");
    type="tensor";

  } else {
    
    cerr << "Cannot use command 'to-tensor' for type "
	 << type << "." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_to_tensor_grid(std::vector<std::string> &sv,
				      bool itive_com) {

  if (type=="table3d") {

    if (sv.size()<2) {
      cerr << "Need slice name." << endl;
      return 1;
    }
    string name=sv[1];

    vector<size_t> sz;
    sz.push_back(table3d_obj.get_nx());
    sz.push_back(table3d_obj.get_ny());
    tensor_grid_obj.resize(2,sz);
    vector<double> grid;
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      grid.push_back(table3d_obj.get_grid_x(i));
    }
    for(size_t i=0;i<table3d_obj.get_ny();i++) {
      grid.push_back(table3d_obj.get_grid_y(i));
    }
    tensor_grid_obj.set_grid_packed(grid);
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
	vector<size_t> ix={i,j};
	tensor_grid_obj.set(ix,table3d_obj.get(i,j,name));
      }
    }
    
    command_del(type);
    clear_obj();
    command_add("tensor_grid");
    type="tensor_grid";

  } else if (type=="tensor") {

    // Get rank
    size_t rank=tensor_obj.get_rank();

    // Get sizes and functions
    vector<string> funcs(rank);
    vector<size_t> sz(rank);
    for(size_t j=0;j<rank;j++) {
      if (sv.size()>j+1) {
	funcs[j]=sv[j+1];
      } else {
	funcs[j]=((string)"func:")+
	  o2scl::szttos(tensor_obj.get_size(j))+":i";
      }
      sz[j]=tensor_obj.get_size(j);
    }

    // Resize tensor_grid object
    tensor_grid_obj.resize(rank,sz);

    // Create grids
    std::vector<std::vector<double> > grid(rank);

    for(size_t j=0;j<rank;j++) {

      // If it's a function, automatically put in the size parameter
      if (funcs[j].find("func:")==0) {
	std::vector<std::string> svx;
	split_string_delim(funcs[j],svx,':');
	if (svx.size()==1) {
	  cerr << "Function specification incomplete." << endl;
	  return 2;
	} else if (svx.size()==2) {
	  svx.push_back(svx[1]);
	}
	svx[1]=o2scl::szttos(sz[j]);
	funcs[j]=svx[0]+':'+svx[1]+':'+svx[2];
	if (verbose>1) {
	  cout << "Added size to function specification: "
	       << funcs[j] << endl;
	}
      }

      int vs_ret=vector_spec(funcs[j],grid[j],false,verbose,false);
      if (vs_ret!=0) {
	cerr << "Function vector_spec() failed." << endl;
	return 1;
      }
    }

    tensor_grid_obj.set_grid(grid);
    
    // Swap data from tensor into tensor_grid
    vector<double> d(tensor_obj.total_size());
    tensor_obj.swap_data(d);
    tensor_grid_obj.swap_data(d);

    command_del(type);
    clear_obj();
    command_add("tensor_grid");
    type="tensor_grid";

  } else {
    
    cerr << "Cannot use command 'to-tensor-grid' for type "
	 << type << "." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_to_tg_fermi(std::vector<std::string> &sv,
				      bool itive_com) {

  if (type=="table3d") {

    if (sv.size()<2) {
      cerr << "Need slice name." << endl;
      return 1;
    }
    string name=sv[1];

    if (sv.size()<3) {
      cerr << "Need number of points." << endl;
      return 2;
    }
    size_t n_points=o2scl::stoszt(sv[2]);

    double low=0.0, high=0.0, width=0.0;
    if (sv.size()>=4) low=o2scl::stod(sv[3]);
    if (sv.size()>=5) high=o2scl::stod(sv[4]);
    if (sv.size()>=6) width=o2scl::stod(sv[5]);
    
    tensor_grid_obj.from_table3d_fermi(table3d_obj,name,n_points,
                                       low,high,width);
    
    command_del(type);
    clear_obj();
    command_add("tensor_grid");
    type="tensor_grid";

  } else {
    
    cerr << "Cannot use command 'to-tg-fermi' for type "
	 << type << "." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_to_vector(std::vector<std::string> &sv,
				      bool itive_com) {

  if (type=="uniform_grid<double>") {

    ug_obj.vector(doublev_obj);
    
    command_del(type);
    clear_obj();
    command_add("double[]");
    type="double[]";

  } else {
    
    cerr << "Cannot use command 'to-tg-fermi' for type "
	 << type << "." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_type(std::vector<std::string> &sv, 
			    bool itive_com) {
  if (type.length()==0) {
    cerr << "No current object to display type of." << endl;
    return 1;
  }
  cout << "The current object has type " << type << "." << endl;
  return 0;
}

int acol_manager::comm_value(std::vector<std::string> &sv, bool itive_com) {

  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(precision);
  
  if (type=="table") {

    // If we have no 'value' entry and we're not in interactive mode,
    // then just presume that the user is asking us to just output
    // the current value
    vector<string> pr, in;
    if (itive_com==true || sv.size()<3) {
    
      pr.push_back("Enter column name");
      pr.push_back("Enter row index");
      pr.push_back("Enter new value (or \"none\") to keep original value");
      int ret=get_input(sv,pr,in,"entry",itive_com);
      if (ret!=0) return ret;

    } else if (sv.size()>=4) {
      in.resize(3);
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
    } else {
      in.resize(2);
      in[0]=sv[1];
      in[1]=sv[2];
    }

    int row;
    int ret2=o2scl::stoi_nothrow(in[1],row);
    if (ret2!=0) {
      std::cerr << "Failed to convert " << in[1]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (row<0) {
      std::cerr << "Conversion of " << in[1]
		<< " resulted in, " << row << ", a negative number." << endl;
      return exc_efailed;
    }

    if (in.size()>=3) {
      // Convert in[2] to lower case
      std::transform(in[2].begin(),in[2].end(),in[2].begin(),::tolower);
    }

    if (in.size()<=2 || in[2]=="none") {
      cout << "Entry for column " << in[0] << " at row " << in[1] << " is "
	   << table_obj.get(in[0],row) << endl;
    } else {
      cout << "Entry for column " << in[0] << " at row " << in[1]
	   << " is has been changed from "
	   << table_obj.get(in[0],row);
      table_obj.set(in[0],row,o2scl::function_to_double(in[2]));
      cout << " to " << table_obj.get(in[0],row) << endl;
    }
    
  } else if (type=="table3d") {

    vector<string> pr, in;
    if (sv.size()==4) {
      // If there are three arguments, then presume a new value
      // wasn't specified, so we have enough information to proceed
      in.push_back(sv[1]);
      in.push_back(sv[2]);
      in.push_back(sv[3]);
    } else {
      pr.push_back("Enter slice name");
      pr.push_back("Enter first index");
      pr.push_back("Enter second index");
      pr.push_back("Enter new value (or \"none\") to keep original value");
      int ret=get_input(sv,pr,in,"entry",itive_com);
      if (ret!=0) return ret;
    }

    int xix;
    int ret2=o2scl::stoi_nothrow(in[1],xix);
    if (ret2!=0) {
      std::cerr << "Failed to convert " << in[1]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (xix<0) {
      std::cerr << "Conversion of " << in[1]
		<< " resulted in, " << xix
		<< ", a negative number." << endl;
      return exc_efailed;
    }

    int yix;
    int ret3=o2scl::stoi_nothrow(in[2],yix);
    if (ret3!=0) {
      std::cerr << "Failed to convert " << in[2]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (yix<0) {
      std::cerr << "Conversion of " << in[2]
		<< " resulted in, " << yix
		<< ", a negative number." << endl;
      return exc_efailed;
    }

    if (in.size()>3) {
      // Convert in[3] to lower case
      std::transform(in[3].begin(),in[3].end(),in[3].begin(),::tolower);
    }
    
    if (in.size()<=3 || in[3]=="none") {
      cout << "Entry for slice " << in[0] << " at (" << in[1] << ","
	   << in[2] << ") is "
	   << table3d_obj.get(xix,yix,in[0]) << endl;
    } else {
      cout << "Entry for slice " << in[0] << " at (" << in[1] << ","
	   << in[2] << ") has been changed from "
	   << table3d_obj.get(xix,yix,in[0]);
      table3d_obj.set(xix,yix,in[0],o2scl::function_to_double(in[3]));
      cout << " to " << table3d_obj.get(xix,yix,in[0]) << endl;
    }
    
  } else if (type=="tensor") {

    size_t rk=tensor_obj.get_rank();
    
    // Handle arguments
    vector<string> in;
    if (sv.size()<rk+1) {
      vector<string> pr;
      for(size_t i=0;i<rk;i++) {
	pr.push_back(((std::string)"Index ")+
		     o2scl::szttos(i));
      }
      pr.push_back("Enter new value (or \"none\") to keep original value");
      int ret=get_input(sv,pr,in,"entry",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t i=0;i<sv.size()-1;i++) {
	in.push_back(sv[i+1]);
      }
    }

    // Parse to array
    vector<size_t> ix;
    for(size_t i=0;i<rk;i++) {
      ix.push_back(o2scl::stoszt(in[i]));
    }

    // Set value if necessary
    if (in.size()>tensor_grid_obj.get_rank()) {
      // convert to lower case
      std::transform(in[rk].begin(),
		     in[rk].end(),
		     in[rk].begin(),::tolower);
      if (in[rk]!="none") {
	tensor_obj.set
	  (ix,o2scl::stod(in[rk]));
      }
    }

    // Output indices and value
    double value=tensor_obj.get(ix);
    cout << "Indices, value: ";
    vector_out(cout,ix,false);
    cout << " " << value << endl;
    
  } else if (type=="tensor_grid") {

    size_t rk=tensor_grid_obj.get_rank();

    // Handle arguments if they weren't specified
    vector<string> in;
    if (sv.size()<rk+1) {
      vector<string> pr;
      for(size_t i=0;i<rk;i++) {
	pr.push_back(((std::string)"Index ")+
		     o2scl::szttos(i));
      }
      pr.push_back("Enter new value (or \"none\") to keep original value");
      int ret=get_input(sv,pr,in,"entry",itive_com);
      if (ret!=0) return ret;
    } else {
      // If the user specified enough indices, then interpret the
      // command as a 'get' request and proceed without prompting for
      // more information.
      for(size_t i=0;i<sv.size()-1;i++) {
	in.push_back(sv[i+1]);
      }
    }

    // Parse to array
    vector<size_t> ix;
    for(size_t i=0;i<rk;i++) {
      ix.push_back(o2scl::stoszt(in[i]));
    }

    // Set value if necessary
    if (in.size()>rk && in[rk]!="none") {
      tensor_grid_obj.set(ix,o2scl::stod(in[rk]));
    }

    // Get associated grid points
    std::vector<double> vals(rk);
    for(size_t i=0;i<rk;i++) {
      vals[i]=tensor_grid_obj.get_grid(i,ix[i]);
    }
    
    // Output indices, grid point, value
    cout << "Indices, grid point, value: ";
    vector_out(cout,ix,false);
    cout << ", ";
    vector_out(cout,vals,false);
    cout << ", " << tensor_grid_obj.get(ix) << endl;

  } else if (type=="int") {
    
    if (sv.size()>1) {
      int_obj=o2scl::stoi(sv[1]);
    }
    cout << "Value of " << obj_name << " is " << int_obj << endl;
    
  } else if (type=="double") {
    
    if (sv.size()>1) {
      int vsret=value_spec(sv[1],double_obj,verbose,false);
      if (vsret!=0) {
	cerr << "Function value_spec() failed." << endl;
	return 1;
      }
    }
    cout << "Value of " << obj_name << " is " << double_obj << endl;
    
  } else if (type=="char") {
    
    if (sv.size()>1) {
      char_obj=sv[1][0];
    }
    cout << "Value of " << obj_name << " is " << char_obj << endl;
    
  } else if (type=="size_t") {
    
    if (sv.size()>1) {
      size_t_obj=o2scl::stoszt(sv[1]);
    }
    cout << "Value of " << obj_name << " is " << size_t_obj << endl;
    
  } else if (type=="string") {
    
    if (sv.size()>1) {
      string_obj=sv[1];
    }
    cout << "Value of " << obj_name << " is " << string_obj << endl;
    
  } else if (type=="double[]") {
    
    if (sv.size()>1) {
      size_t ix=o2scl::stoszt(sv[1]);
      if (ix>=doublev_obj.size()) {
        cerr << "The double[] object only has "
             << doublev_obj.size() << " elements. Use the 'resize' "
             << "command to resize." << endl;
        return 1;
      }
      if (sv.size()>2) {
        doublev_obj[ix]=o2scl::function_to_double(sv[2]);
      }
      cout << "The value at index " << ix << " is " << doublev_obj[ix]
           << "." << endl;
      return 0;
    }

  } else if (type=="string[]") {
    
    if (sv.size()>1) {
      size_t ix=o2scl::stoszt(sv[1]);
      if (ix>=stringv_obj.size()) {
        cerr << "The string[] object only has "
             << stringv_obj.size() << " elements. Use the 'resize' "
             << "command to resize." << endl;
        return 1;
      }
      if (sv.size()>2) {
        stringv_obj[ix]=sv[2];
      }
      cout << "The value at index " << ix << " is " << stringv_obj[ix]
           << "." << endl;
      return 0;
    }
    
  } else if (type=="int[]") {
    
    if (sv.size()>1) {
      size_t ix=o2scl::stoszt(sv[1]);
      if (ix>=intv_obj.size()) {
        cerr << "The double[] object only has "
             << intv_obj.size() << " elements. Use the 'resize' "
             << "command to resize." << endl;
        return 1;
      }
      if (sv.size()>2) {
        intv_obj[ix]=o2scl::stoi(sv[2]);
      }
      cout << "The value at index " << ix << " is " << intv_obj[ix]
           << "." << endl;
      return 0;
    }
    
  } else if (type=="size_t[]") {
    
    if (sv.size()>1) {
      size_t ix=o2scl::stoszt(sv[1]);
      if (ix>=size_tv_obj.size()) {
        cerr << "The double[] object only has "
             << size_tv_obj.size() << " elements. Use the 'resize' "
             << "command to resize." << endl;
        return 1;
      }
      if (sv.size()>2) {
        size_tv_obj[ix]=o2scl::stoi(sv[2]);
      }
      cout << "The value at index " << ix << " is " << size_tv_obj[ix]
           << "." << endl;
      return 0;
    }
    
  } else {
    
    cerr << "Command 'value' not implemented for type " << type << " ."
	 << endl;
    
    return exc_efailed;
  }
  
  return 0;
}
  
int acol_manager::comm_value_grid(std::vector<std::string> &sv,
				  bool itive_com) {

  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(precision);
  
  if (type=="table") {

    vector<string> pr, in;
    pr.push_back("Enter index column name");
    pr.push_back("Enter index column value");
    pr.push_back("Enter target column name");
    pr.push_back("Enter new value (or \"none\") to keep original value");
    int ret=get_input(sv,pr,in,"value-grid",itive_com);
    if (ret!=0) return ret;

    double val=o2scl::stod(in[1]);
    int row=table_obj.lookup(in[0],val);
    cout << "Looking up value " << val << " in column " << in[0]
	 << " results in row " << row << " with value "
	 << table_obj.get(in[0],row) << endl;
    
    if (in.size()>=4) {
      // Convert in[3] to lower case
      std::transform(in[3].begin(),in[3].end(),in[3].begin(),::tolower);
    }
    
    if (in.size()<=2 || in[3]=="none") {
      cout << "Entry for column " << in[2] << " at row " << row << " is "
	   << table_obj.get(in[2],row) << endl;
    } else {
      cout << "Entry for column " << in[2] << " at row " << row
	   << " is has been changed from "
	   << table_obj.get(in[2],row);
      table_obj.set(in[0],row,o2scl::function_to_double(in[3]));
      cout << " to " << table_obj.get(in[2],row) << endl;
    }
    
  } else if (type=="table3d") {

    vector<string> pr, in;
    pr.push_back("Enter slice name");
    pr.push_back("Enter first grid point");
    pr.push_back("Enter second grid point");
    pr.push_back("Enter new value (or \"none\") to keep original value");
    int ret=get_input(sv,pr,in,"value-grid",itive_com);
    if (ret!=0) return ret;

    int xix;
    int ret2=o2scl::stoi_nothrow(in[1],xix);
    if (ret2!=0) {
      std::cerr << "Failed to convert " << in[1]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (xix<0) {
      std::cerr << "Conversion of " << in[1]
		<< " resulted in, " << xix
		<< ", a negative number." << endl;
      return exc_efailed;
    }

    int yix;
    int ret3=o2scl::stoi_nothrow(in[2],yix);
    if (ret3!=0) {
      std::cerr << "Failed to convert " << in[2]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (yix<0) {
      std::cerr << "Conversion of " << in[2]
		<< " resulted in, " << yix
		<< ", a negative number." << endl;
      return exc_efailed;
    }

    if (in.size()>3) {
      // Convert in[3] to lower case
      std::transform(in[3].begin(),in[3].end(),in[3].begin(),::tolower);
    }
    
    if (in.size()<=3 || in[3]=="none") {
      cout << "Entry for slice " << in[0] << " at (" << in[1] << ","
	   << in[2] << ") is "
	   << table3d_obj.get(xix,yix,in[0]) << endl;
    } else {
      cout << "Entry for slice " << in[0] << " at (" << in[1] << ","
	   << in[2] << ") has been changed from "
	   << table3d_obj.get(xix,yix,in[0]);
      table3d_obj.set(xix,yix,in[0],o2scl::function_to_double(in[3]));
      cout << " to " << table3d_obj.get(xix,yix,in[0]) << endl;
    }
    
  } else if (type=="tensor_grid") {

    size_t rk=tensor_grid_obj.get_rank();
    
    vector<string> in;
    if (sv.size()<rk+1) {
      // Handle arguments
      vector<string> pr;
      for(size_t i=0;i<rk;i++) {
	pr.push_back(((std::string)"Value for index ")+
		     o2scl::szttos(i));
      }
      pr.push_back("Enter new value (or \"none\") to keep original value");
      int ret=get_input(sv,pr,in,"value-grid",itive_com);
      if (ret!=0) return ret;
      
    } else {
      // If the user specified enough indices, then interpret the
      // command as a 'get' request and proceed without prompting for
      // more information.
      for(size_t i=0;i<sv.size()-1;i++) {
	in.push_back(sv[i+1]);
      }
    }

    // Parse to array
    vector<double> vals;
    for(size_t i=0;i<rk;i++) {
      vals.push_back(o2scl::stod(in[i]));
    }

    // Set value if necessary
    if (in.size()>rk && in[rk]!="none") {
      tensor_grid_obj.set_val
	(vals,o2scl::stod(in[rk]));
    }

    // Lookup closest grid point
    double value=tensor_grid_obj.get_val(vals,vals);
    vector<size_t> ix(rk);
    tensor_grid_obj.lookup_grid_vec(vals,ix);

    if (scientific) cout.setf(ios::scientific);
    else cout.unsetf(ios::scientific);
    cout.precision(precision);
    
    // Output indices, grid point, value
    cout << "Indices, grid point, value: ";
    vector_out(cout,ix,false);
    cout << " ";
    vector_out(cout,vals,false);
    cout << " " << value << endl;
    
  } else {
    cerr << "Command 'value-grid' not implemented for type " << type << " ."
	 << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_version(std::vector<std::string> &sv, bool itive_com) {

  cout << "\n" << cl->desc << endl;

  terminal ter;
  if (cl->cmd_name==((string)"acol")) {
    cout << ((string)"Compiled at ")+((string)__TIME__)+" on "+
      ((string)__DATE__)+" for "+exec_color+"O₂scl"+
      default_color+", version "+((string)VERSION)+".\n" << endl;
  }

  cout << exec_color << "O₂scl" << default_color
       << " version: " << o2scl_settings.o2scl_version() << endl;
  cout << "Range checking: " << o2scl_settings.range_check() << endl;
  if (true) {
    unsigned maj, min, rel;
    o2scl_settings.hdf5_header_version(maj,min,rel);
    cout << "  HDF5 version numbers when O₂scl was compiled: "
	 << maj << " " << min << " " << rel << endl;
    o2scl_settings.hdf5_lib_version(maj,min,rel);
    cout << "  HDF5 version numbers in libraries currently linked: "
	 << maj << " " << min << " " << rel << endl;
  }
  cout << "  HDF5 compression support: "
       << o2scl_settings.hdf5_compression_support() << endl;
  cout << "Python support: " << o2scl_settings.python_support() << endl;
  if (o2scl_settings.python_support()) {
    cout << "Python version: " << o2scl_settings.py_version() << endl;
  }
  cout << "Armadillo support: "
       << o2scl_settings.armadillo_support() << endl;
  cout << "Eigen support: " << o2scl_settings.eigen_support() << endl;
  cout << "FFTW support: " << o2scl_settings.fftw_support() << endl;
  cout << "Cubature support: " << o2scl_settings.cubature_support() << endl;
  cout << "OpenMP support: " << o2scl_settings.openmp_support() << endl;
  cout << "Readline support: " << o2scl_settings.readline_support() << endl;
  cout << "Module support: " << o2scl_settings.module_support() << endl;
  cout << "MPFR support: " << o2scl_settings.mpfr_support() << endl;
  cout << "Ncurses support: " << o2scl_settings.ncurses_support() << endl;
  cout << "Data directory: " << o2scl_settings.get_data_dir() << endl;

  cout << "Documentation directory: "
       << o2scl_settings.get_doc_dir() << endl;
  cout << "Local documentation URL:\n  file://"
       << o2scl_settings.get_doc_dir() << "html/index.html" << endl;
  cout << "Online documentation URL:\n  http://neutronstars.utk.edu/code/o2scl"
       << "/html/index.html" << endl;
  cout << "System type: " << o2scl_settings.system_type() << endl;
  cout << endl;
  cout << "o2scl_name: " << o2scl_settings.o2scl_name() << endl;
  cout << "o2scl_package: " << o2scl_settings.o2scl_package() << endl;
  cout << "o2scl_bugreport: " << o2scl_settings.o2scl_bugreport() << endl;
  cout << "o2scl_string: " << o2scl_settings.o2scl_string() << endl;
  cout << "o2scl_tarname: " << o2scl_settings.o2scl_tarname() << endl;
  cout << endl;

  // Typically,
  // type              digits10 max_digits10 max          log_prec
  // --------------------------------------------------------------
  // double            15       17           1.8e308       39.1
  // long double       18       21           1.2e4932      48.4
  // cpp_dec_float_35  35       64           1.0e67108864 147.4
  // cpp_dec_float_50  50       80           1.0e67108864 184.2
  // cpp_dec_float_100 100      128          1.0e67108864 294.7

  cout.precision(4);
  cout.setf(ios::left);

  cout.width(18);
  cout << "type";
  cout << "d10 ";
  cout << "md10 ";
  cout.width(17);
  cout << "max" << " ";
  cout.width(15);
  cout << "log10(max_d10)" << " ";
  cout << "epsilon" << endl;
  cout << "------------------------------------"
       << "-----------------------------------"
       << endl;
  
  cout.width(18);
  cout << "double";
  cout.width(3);
  cout << std::numeric_limits<double>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<double>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<double>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,std::numeric_limits<double>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<double>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "long double";
  cout.width(3);
  cout << std::numeric_limits<long double>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<long double>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<long double>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,std::numeric_limits<long double>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<long double>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "cpp_dec_float_25";
  cout.width(3);
  cout << std::numeric_limits<cpp_dec_float_25>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<cpp_dec_float_25>::max_digits10 << " "; 
  cout.width(17);
  cout << std::numeric_limits<cpp_dec_float_25>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<cpp_dec_float_25>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<cpp_dec_float_25>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "cpp_dec_float_35";
  cout.width(3);
  cout << std::numeric_limits<cpp_dec_float_35>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<cpp_dec_float_35>::max_digits10 << " "; 
  cout.width(17);
  cout << std::numeric_limits<cpp_dec_float_35>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<cpp_dec_float_35>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<cpp_dec_float_35>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "cpp_dec_float_50";
  cout.width(3);
  cout << std::numeric_limits<cpp_dec_float_50>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<cpp_dec_float_50>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<cpp_dec_float_50>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<cpp_dec_float_50>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<cpp_dec_float_50>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "cpp_dec_float_100";
  cout.width(3);
  cout << std::numeric_limits<cpp_dec_float_100>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<cpp_dec_float_100>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<cpp_dec_float_100>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<cpp_dec_float_100>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<cpp_dec_float_100>::epsilon()
       << std::endl;

#ifdef O2SCL_SET_MPFR

  cout.width(18);
  cout << "mpfr_25";
  cout.width(3);
  cout << std::numeric_limits<mpfr_25>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<mpfr_25>::max_digits10 << " "; 
  cout.width(17);
  cout << std::numeric_limits<mpfr_25>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<mpfr_25>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<mpfr_25>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "mpfr_35";
  cout.width(3);
  cout << std::numeric_limits<mpfr_35>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<mpfr_35>::max_digits10 << " "; 
  cout.width(17);
  cout << std::numeric_limits<mpfr_35>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<mpfr_35>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<mpfr_35>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "mpfr_50";
  cout.width(3);
  cout << std::numeric_limits<mpfr_50>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<mpfr_50>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<mpfr_50>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<mpfr_50>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<mpfr_50>::epsilon()
       << std::endl;
  
  cout.width(18);
  cout << "mpfr_100";
  cout.width(3);
  cout << std::numeric_limits<mpfr_100>::digits10 << " ";
  cout.width(4);
  cout << std::numeric_limits<mpfr_100>::max_digits10 << " ";
  cout.width(17);
  cout << std::numeric_limits<mpfr_100>::max() << " ";
  cout.width(15);
  cout << log(pow(10.0,
                  std::numeric_limits<mpfr_100>::max_digits10));
  cout << " ";
  cout << std::numeric_limits<mpfr_100>::epsilon()
       << std::endl;

#endif
  
  cout.unsetf(ios::left);
  cout << endl;
  
  cout << "config.h: " << endl;
  o2scl_settings.config_h_report();
  return 0;
}

int acol_manager::comm_wdocs(std::vector<std::string> &sv, bool itive_com) {

  string cmd;

#ifdef O2SCL_LINUX
  cmd="xdg-open ";
#else
#ifdef O2SCL_OSX
  cmd="open "; 
#else
  cmd="xdg-open ";
#endif
#endif
  
  if (sv.size()>=3 || (sv.size()==2 && sv[1]!="dev")) {
    bool dev=false;
    string term=sv[1];
    string section;
    
    if (sv.size()>=3 && sv[1]==((string)"dev")) {
      term=sv[2];
      dev=true;
    }
    
    if (term.length()>40) {
      term=term.substr(0,40);
    }
    
    for(size_t i=0;i<term.length();i++) {
      // If there is a space, then replace it with "%20"
      if (term[i]==' ') {
	term.replace(term.begin()+i,term.begin()+i+1,"%20");
	i=0;
      } else if (!isalnum(term[i]) && term[i]!='%' && term[i]!='_') {
	// If there is some other non-alphanumeric, remove it
	term.replace(term.begin()+i,term.begin()+i+1,"");
	i=0;
      }
    }
    if (section=="part") {
      if (dev) {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl-dev/part/html/search.html?q="+term+" &";
      } else {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl/part/html/search.html?q="+term+" &";
      }
    } else if (section=="eos") {
      if (dev) {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl-dev/eos/html/search.html?q="+term+" &";
      } else {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl/eos/html/search.html?q="+term+" &";
      }
    } else {
      if (dev) {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl-dev/html/search.html?q="+term+" &";
      } else {
        cmd+=((string)"https://neutronstars.utk.edu/code/")+
          "o2scl/html/search.html?q="+term+" &";
      }
    }
  } else if (sv[1]=="dev") {
    cmd+="https://neutronstars.utk.edu/code/o2scl-dev/html/acol.html &";
  } else {
    cmd+="https://neutronstars.utk.edu/code/o2scl/html/acol.html &";
  }
  
  cout << "Using command: " << cmd << endl;

  int xret=system(cmd.c_str());
  
  return 0;
}

int acol_manager::comm_wstats(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }
  
  if (table_obj.get_nlines()==0) {
    cerr << "No table to analyze." << endl;
    return exc_efailed;
  }
  
  vector<std::string> in, pr;
  pr.push_back("Enter column to get info on");
  pr.push_back("Enter column for weights");
  int ret=get_input(sv,pr,in,"wstats",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0]
	 << "' in acol command wstats." << endl;
    return exc_efailed;
  }
  if (table_obj.is_column(in[1])==false) {
    cerr << "Could not find column named '" << in[1]
	 << "' in acol command wstats." << endl;
    return exc_efailed;
  }

  const vector<double> &cref=table_obj.get_column(in[0]);
  const vector<double> &wref=table_obj.get_column(in[1]);
  cout << "N        : " << table_obj.get_nlines() << endl;
  cout << "Mean     : "
       << wvector_mean(table_obj.get_nlines(),cref,wref) << endl;
  cout << "Std. dev.: "
       << wvector_stddev(table_obj.get_nlines(),cref,wref) << endl;
  
  return 0;
}

void acol_manager::xml_replacements(std::string &s,
                                    std::vector<std::string> &clist) {

  terminal ter;

  /*
    We can't make the full color replacements here because
    the actual colors aren't known until runtime, so we
    use [c], [d], [e], [h], [p], [t], and [u]
   */
  
  vector<string> subs={((string)"<computeroutput> Vector ")+
    "specifications </computeroutput>",
    "'acol -help [h]vector-spec[d]",
    "<computeroutput> Multiple vector specifications </computeroutput>",
    "'acol -help [h]mult-vector-spec[d]",
    ((string)"\"acol -help <computeroutput> Function ")+
    "specifications </computeroutput> \"",
    "'acol -help [h]functions[d]",
    "`String list specifications`","[e]acol -help strings-spec[d]",
    "`Vector specifications`","[e]acol -help vector-spec[d]",
    "`Multiple vector specifications`","[e]acol -help mult-vector-spec[d]",
    "`Value specifications`","[e]acol -help value-spec[d]",
    "`Function specifications`","[e]acol -help functions[d]",
    "cpp:func:`o2scl_hdf::functions()`","[e]acol -help functions[d]",
    "cpp:func:`o2scl_hdf::value_spec()`","[e]acol -help value-spec[d]",
    "cpp:func:`o2scl_hdf::vector_spec()`","[e]acol -help vector-spec[d]",
    "cpp:func:`o2scl_hdf::mult_vector_spec()`",
    "[e]acol -help mult_vector-spec[d]",
    "cpp:func:`o2scl_hdf::index_spec()`","[e]acol -help index-spec[d]",
    "cpp:func:`o2scl_hdf::strings_spec()`","[e]acol -help strings-spec[d]"};

  string_replace(s,"\n"," ");
  
  // Make the manual replacements from the 'subs' list
  // above
  for(size_t i=0;i<subs.size();i+=2) {
    string_replace(s,subs[i],subs[i+1]);
  }
  
  // Make all of the type replacements
  for(size_t i=0;i<type_list.size();i++) {
    string_replace(s,"<computeroutput> <ref> "+type_list[i]+
                   " </ref> </computeroutput>",
                   "[t]"+type_list[i]+"[d]");
    string_replace(s,"<computeroutput> "+type_list[i]+
                   " </computeroutput>",
                   "[t]"+type_list[i]+"[d]");
    string_replace(s,"<ref> "+type_list[i]+
                   " </ref>",
                   "[t]"+type_list[i]+"[d]");
  }
                
  // Make the command replacements
  for(size_t i=0;i<clist.size();i++) {
    string_replace(s,"<computeroutput> "+clist[i]+
                   " </computeroutput>",
                   "[c]"+clist[i]+"[d]");
  }
  
  // Make the command replacements from the current parameter list
  for(cli::par_t it=cl->par_list.begin();it!=cl->par_list.end();it++) {
    string_replace(s,"<computeroutput> "+it->first+
                   " </computeroutput>",
                   "[c]"+it->first+"[d]");
  }
  
  // Make the help topic replacements
  string_replace(s,"<computeroutput> functions </computeroutput>",
                 "[h]functions[d]");
  string_replace(s,"<computeroutput> types </computeroutput>",
                 "[h]types[d]");
  string_replace(s,"<computeroutput> value-spec </computeroutput>",
                 "[h]value-spec[d]");
  string_replace(s,"<computeroutput> vector-spec </computeroutput>",
                 "[h]vector-spec[d]");
  string_replace(s,"<computeroutput> mult-vector-spec </computeroutput>",
                 "[h]mult-vector-spec[d]");
  string_replace(s,"<computeroutput> strings-spec </computeroutput>",
                 "[h]strings-spec[d]");
  string_replace(s,"<computeroutput> index-spec </computeroutput>",
                 "[h]index-spec[d]");

  // Other miscellaneous replacements
  string_replace(s,"<itemizedlist> <listitem>","★");
  string_replace(s,"</listitem> <listitem>","★");
  string_replace(s,"</listitem> </itemizedlist>","");
  string_replace(s,"<orderedlist> <listitem>","★");
  string_replace(s,"</listitem> </orderedlist>","");
  string_replace(s,"<simplesect>","");
  string_replace(s,"</simplesect>","");
  string_replace(s," ★","★");

  /*
  if (s.find("<computeroutput>")!=std::string::npos &&
      s.find("</computeroutput>",s.find("<computeroutput>"))!=
      std::string::npos) {
    size_t s1=s.find("<computeroutput>");
    size_t s2=s.find("</computeroutput>",s.find("<computeroutput>"));
    if (s2-s1>16 && s2+17<s.length()) {
      s.replace(s2,17,default_color);
      s.replace(s1,16,exec_color);
    }
  }
  */
  if (s.find("Arguments:")==0) {
    string_replace(s,"<computeroutput> ","");
    string_replace(s," </computeroutput>","");
  } else {
    string_replace(s,"<computeroutput> ","[e]");
    string_replace(s," </computeroutput>","[d]");
  }
  string_replace(s,"<linebreak> ","");
  string_replace(s," </linebreak>","");

  string_replace(s,"<verbatim> embed:rst","");
  string_replace(s," </verbatim>","");
  
  string_replace(s,"See:ref:"," ");
  string_replace(s,"See:cpp:func:`","See ");
  string_replace(s," See :[","See [");
  string_replace(s," and :["," and [");
  string_replace(s,"See[","See [");

  string_replace(s,"  "," ");
  string_replace(s," )",")");
  string_replace(s," ,",",");
  string_replace(s," .",".");
  string_replace(s," :",":");
  string_replace(s,"`","");
  /*
  string_replace(s,"</itemizedlist>","");
  string_replace(s,"</listitem>","");
  */

  string_replace(s,"<formula> $ 10^{-\\mathrm{precision}-1} $ </formula>",
                 "10^{-precision-1}");
  string_replace(s,"<formula> $ y(x) $ </formula>","y(x)");
                  
  return;
}

int acol_manager::comm_xml_to_o2(std::vector<std::string> &sv,
                                 bool itive_com) {

#ifdef O2SCL_PUGIXML

  terminal ter;
  
  // XML walkers
  ostream_walker ow;
  vec_string_walker vsw;
  vsw.indent=false;
  
  std::string stmp;
  
  vector<vector<std::string>> cmd_doc_strings_loc, help_doc_strings_loc,
    param_doc_strings_loc;

  // Create a list of the current options
  vector<string> clist=cl->get_option_list();

  // Add the remaining type-specific options
  for (std::map<std::string,std::vector<std::string> >::iterator
         it=type_comm_list.begin();it!=type_comm_list.end();
       it++) {
    for(size_t ii=0;ii<it->second.size();ii++) {
      if (std::find(clist.begin(),clist.end(),it->second[ii])==clist.end()) {
        clist.push_back(it->second[ii]);
      }
    }
  }

  if (verbose>2) {
    cout << "clist: ";
    vector_out(cout,clist,true);
  }

  // Loop over every command
  for(size_t j=0;j<clist.size();j++) {
    
    // The command name, the brief description, the parameter
    // description, and then all of the paragraphs in the long-form
    // documentation, in that order.
    vector<std::string> vs_tmp;

    pugi::xml_document doc;
    pugi::xml_document doc2;
    
    std::string cmd_name=clist[j];
    std::string fn_name="comm_";
    for(size_t k=0;k<cmd_name.length();k++) {
      if (cmd_name[k]=='-') {
        fn_name+='_';
      } else {
        fn_name+=cmd_name[k];
      }
    }
    
    std::string fn="doc/o2scl/xml/classo2scl__acol_1_1acol__manager.xml";

    pugi::xml_node n3=doxygen_xml_member_get
      (fn,"acol_manager",fn_name,"briefdescription",doc);

    if (n3!=0) {
      
      // We found a brief description, so add the command name
      // to the list
      vs_tmp.push_back(cmd_name);
      
      if (verbose>2) {
        cout << "brief_desc name,value: "
             << n3.name() << " " << n3.value() << endl;
        n3.traverse(ow);
      }
      
      // Store the XML for the brief description in vsw.output
      n3.traverse(vsw);
      
      // Combine with spaces, and remove the outer paragraph
      stmp="";
      for(size_t k=0;k<vsw.output.size();k++) {
        if (stmp.length()!=0) {
          stmp+=' ';
        }
        if (vsw.output[k]!=((string)"<para>") &&
            vsw.output[k]!=((string)"</para>")) {
          stmp+=vsw.output[k];
        }
      }
      
      xml_replacements(stmp,clist);

      // Add brief description to stmp
      vs_tmp.push_back(stmp);

      pugi::xml_node n4=doxygen_xml_member_get
        (fn,"acol_manager",fn_name,"detaileddescription",doc2);
      
      if (n4!=0) {
        
        if (verbose>2) {
          cout << "desc name,value: "
               << n4.name() << " " << n4.value() << endl;
          n4.traverse(ow);
        }
        
        bool found=false;
        
        // Store the XML for the detailed description in vsw.output
        n4.traverse(vsw);
        
        bool done=false;
        stmp="";
        
        for(size_t k=0;k<vsw.output.size() && done==false;k++) {
          
          if (vsw.output[k].find("End of runtime documentation.")!=
              string::npos) {
            done=true;
          } else {
            if (vsw.output[k]!=((string)"<para>")) {
              if (vsw.output[k]==((string)"</para>")) {
                found=true;
                
                if (verbose>1) {
                  cout << "stmp before: " << stmp << endl;
                }
      
                xml_replacements(stmp,clist);
                
                if (verbose>1) {
                  cout << "stmp after: " << stmp << endl;
                }
      
                if (stmp.length()>0) vs_tmp.push_back(stmp);
                stmp.clear();
              } else {
                if (stmp.length()==0) stmp=vsw.output[k];
                else stmp+=' '+vsw.output[k];
              }
            }
          }
        }
        
        // End of if (n4!=0) 
      }
      
      if (vs_tmp.size()>=2) {
        if (verbose>0) {
          cout << "Adding documention for command: " << vs_tmp[0]
               << endl;
          for(size_t jj=0;jj<vs_tmp.size();jj++) {
            cout << jj << ": \"" << vs_tmp[jj] << "\"" << endl;
          }
          cout << endl;
        }
        cmd_doc_strings_loc.push_back(vs_tmp);
      }
      
      // End of if (n3!=0) 
    }
  }

  // Go through all the parameters
  for(cli::par_t itp=cl->par_list.begin();itp!=cl->par_list.end();itp++) {

    // This parameter name and then all of the paragraphs in the
    // long-form description, in that order.
    vector<std::string> vs_tmp;
    
    pugi::xml_document doc;
    pugi::xml_document doc2;
    
    if (verbose>1) {
      cout << "parameter name, doc_name: " << itp->first << " "
           << itp->second->doc_name << endl;
    }
    
    std::string fn="doc/o2scl/xml/classo2scl__acol_1_1acol__manager.xml";
    
    pugi::xml_node n3=doxygen_xml_member_get
      (fn,itp->first,itp->first,"briefdescription",doc);
    if (n3==0) {
      cout << "n3 0 " << itp->first << " "
           << itp->second->doc_name << " " << fn << endl;
    }
    
    if (n3!=0) {
      
      // We found a brief description, so add the parameter name
      // to the list
      vs_tmp.push_back(itp->first);
      
      if (verbose>2) {
        cout << "brief desc. name, value: "
             << n3.name() << " " << n3.value() << endl;
        n3.traverse(ow);
      }
      
      // Store the XML for the brief description in vsw.output
      n3.traverse(vsw);
      
      // Combine with spaces, and remove the outer paragraph
      stmp="";
      for(size_t k=0;k<vsw.output.size();k++) {
        if (stmp.length()!=0) {
          stmp+=' ';
        }
        if (vsw.output[k]!=((string)"<para>") &&
            vsw.output[k]!=((string)"</para>")) {
          stmp+=vsw.output[k];
        }
      }

      xml_replacements(stmp,clist);

      // Remove trailing spaces
      while (stmp.length()>2 && stmp[stmp.length()-1]==' ') {
        stmp=stmp.substr(0,stmp.length()-1);
      }
      
      // Add a period at the end, because this brief description
      // is combined with the detailed description by the
      // cli 'help' command above.
      if (stmp.length()>2 && stmp[stmp.length()-1]!='.') {
        stmp+='.';
      }
      
      // Add brief description to stmp
      vs_tmp.push_back(stmp);

      // Look for a detailed description
      pugi::xml_node n4=doxygen_xml_member_get
        (fn,itp->first,itp->first,
         "detaileddescription",doc2);
      
      if (n4!=0) {
        
        if (verbose>2) {
          cout << "detailed desc. name, value: "
               << n4.name() << " " << n4.value() << endl;
          n4.traverse(ow);
        }
        
        // Store the XML for the brief description in vsw.output
        n4.traverse(vsw);
        
        bool done=false;
        stmp="";
        
        for(size_t k=0;k<vsw.output.size() && done==false;k++) {
          
          if (vsw.output[k].find("End of runtime documentation.")!=
              string::npos) {
            done=true;
          } else {
            if (vsw.output[k]!=((string)"<para>")) {
              if (vsw.output[k]==((string)"</para>")) {
                
                xml_replacements(stmp,clist);
                
                if (stmp.length()>0) vs_tmp.push_back(stmp);
                stmp.clear();
              } else {
                if (stmp.length()==0) stmp=vsw.output[k];
                else stmp+=' '+vsw.output[k];
              }
            }
          }
        }

        // End of if (n4!=0)
      }
      
      if (vs_tmp.size()>=2) {
        if (verbose>0) {
          cout << "Adding documentation for parameter " << vs_tmp[0] << endl;
          for(size_t jj=0;jj<vs_tmp.size();jj++) {
            cout << jj << ": \"" << vs_tmp[jj] << "\"" << endl;
          }
          cout << endl;
        }
        param_doc_strings_loc.push_back(vs_tmp);
      }

      // End of if (n3!=0)
    }
    
  }
  
  // Help topic list
  vector<string> flist={"value_spec","vector_spec","mult_vector_spec",
    "strings_spec","index_spec","functions"};
  
  for(size_t j=0;j<flist.size();j++) {

    vector<std::string> vs_tmp;

    pugi::xml_document doc;
    pugi::xml_document doc2;

    std::string fn_name=flist[j];
    
    std::string fn="doc/o2scl/xml/namespaceo2scl__hdf.xml";
    
    pugi::xml_node n3=doxygen_xml_get
      (fn,fn_name,"briefdescription",doc);

    pugi::xml_node n4=doxygen_xml_get
      (fn,fn_name,"detaileddescription",doc2);
    
    if (n3!=0 && n4!=0) {

      if (verbose>2) {
        cout << "dxmg: " << n3.name() << " " << n3.value() << endl;
        n3.traverse(ow);
        
        cout << "dxmg: " << n4.name() << " " << n4.value() << endl;
        n4.traverse(ow);
      }
      
      pugi::xml_node_iterator it=n4.begin();

      vs_tmp.push_back(fn_name);
      vs_tmp.push_back(n3.child_value("para"));

      bool found=false;
      
      n4.traverse(vsw);
      //cout << vsw.output.size() << endl;
      bool done=false;
      std::string stmpx;
      for(size_t k=0;k<vsw.output.size() && done==false;k++) {
        if (vsw.output[k].find("End of runtime documentation.")!=
            string::npos) {
          done=true;
        } else {
          if (vsw.output[k]!=((string)"<para>")) {
            if (vsw.output[k]==((string)"</para>")) {
              if (verbose>1) {
                cout << "stmp before: " << stmpx << endl;
              }
              found=true;

              xml_replacements(stmpx,clist);
              if (verbose>1) {
                cout << "stmp after: " << stmpx << endl;
              }
              
              vs_tmp.push_back(stmpx);
              stmpx.clear();
            } else {
              if (stmpx.length()==0) stmpx=vsw.output[k];
              else stmpx+=' '+vsw.output[k];
            }
          }
        }
      }

      if (found) {
        if (vs_tmp.size()>=4 && verbose>0) {
          cout << "Adding doccumentation for help topic " << vs_tmp[0] << endl;
          for(size_t jj=0;jj<vs_tmp.size();jj++) {
            cout << jj << ": \"" << vs_tmp[jj] << "\"" << endl;
          }
          cout << endl;
        }
        help_doc_strings_loc.push_back(vs_tmp);
        
      }
      
    }

  }

  hdf_file hf;
  hf.open_or_create("data/o2scl/acol_docs.o2");
  hf.sets_vec_vec_copy("cmd_doc_strings",cmd_doc_strings_loc);
  hf.sets_vec_vec_copy("help_doc_strings",help_doc_strings_loc);
  hf.sets_vec_vec_copy("param_doc_strings",param_doc_strings_loc);
  hf.close();
  cout << "Created file data/o2scl/acol_docs.o2." << endl;

#else

  cout << "Pugixml must be enabled to create the runtime documentation "
       << "from the doxygen\n XML output." << endl;
  
#endif

  return 0;
}

int acol_manager::comm_x_name(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d" || table3d_obj.is_xy_set()==false) {
      cerr << "No table3d object or no grid." << endl;
      return exc_efailed;
    }

    if (sv.size()==1) {
      cout << "X grid is named " << table3d_obj.get_x_name() << endl;
    } else {
      table3d_obj.set_x_name(sv[1]);
      cout << "X grid is now named " << table3d_obj.get_x_name() << endl;
    }

  } else {
    
    cerr << "Command 'x-name' not implemented for " << type
	 << " objects." << endl;
    return exc_efailed;
  }
  
  return 0;
}

int acol_manager::comm_y_name(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d" || table3d_obj.is_xy_set()==false) {
      cerr << "No table3d object or no grid." << endl;
      return exc_efailed;
    }

    if (sv.size()==1) {
      cout << "Y grid is named " << table3d_obj.get_y_name() << endl;
    } else {
      table3d_obj.set_y_name(sv[1]);
      cout << "Y grid is now named " << table3d_obj.get_y_name() << endl;
    }

  } else {
    
    cerr << "Command 'y-name' not implemented for " << type
	 << " objects." << endl;
    return exc_efailed;
  }
  
  return 0;
}

