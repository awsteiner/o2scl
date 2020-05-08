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
    tensor_obj.convert_table3d_sum(ix_x,ix_y,table3d_obj,in[0],
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
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
	ix[ix_x]=i;
	ix[ix_y]=j;
	table3d_obj.set(i,j,in[2],tensor_obj.get(ix));
      }
    }
    } else if (type=="tensor<int>") {
    for(size_t i=0;i<table3d_obj.get_nx();i++) {
      for(size_t j=0;j<table3d_obj.get_ny();j++) {
	ix[ix_x]=i;
	ix[ix_y]=j;
	table3d_obj.set(i,j,in[2],tensor_int_obj.get(ix));
      }
    }
    } else {
      for(size_t i=0;i<table3d_obj.get_nx();i++) {
	for(size_t j=0;j<table3d_obj.get_ny();j++) {
	  ix[ix_x]=i;
	  ix[ix_y]=j;
	  table3d_obj.set(i,j,in[2],tensor_size_t_obj.get(ix));
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

int acol_manager::comm_to_tensor_grid(std::vector<std::string> &sv,
				      bool itive_com) {

  if (type=="tensor") {

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
	std::vector<std::string> sv;
	split_string_delim(funcs[j],sv,':');
	if (sv.size()==1) {
	  cerr << "Function specification incomplete." << endl;
	  return 2;
	} else if (sv.size()==2) {
	  sv.push_back(sv[1]);
	}
	sv[1]=o2scl::szttos(sz[j]);
	funcs[j]=sv[0]+':'+sv[1]+':'+sv[2];
	if (verbose>1) {
	  cout << "Added size to function specification: "
	       << funcs[j] << endl;
	}
      }

      int vs_ret=vector_spec(funcs[j],grid[j],verbose,false);
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

int acol_manager::comm_value(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()>1) {
    if (type=="int") {
      int_obj=o2scl::stoi(sv[1]);
    } else if (type=="double") {
      int vsret=value_spec(sv[1],double_obj,verbose,false);
      if (vsret!=0) {
	cerr << "Function value_spec() failed." << endl;
	return 1;
      }
    } else if (type=="char") {
      char_obj=sv[1][0];
    } else if (type=="size_t") {
      size_t_obj=o2scl::stoszt(sv[1]);
    } else if (type=="string") {
      string_obj=sv[1];
    }
  }
  
  if (type=="int") {
    cout << "Value of " << obj_name << " is " << int_obj << endl;
  } else if (type=="double") {
    cout << "Value of " << obj_name << " is " << double_obj << endl;
  } else if (type=="char") {
    cout << "Value of " << obj_name << " is " << char_obj << endl;
  } else if (type=="size_t") {
    cout << "Value of " << obj_name << " is " << size_t_obj << endl;
  } else if (type=="string") {
    cout << "Value of " << obj_name << " is " << string_obj << endl;
  }
  
  return 0;
}
  
#ifdef O2SCL_NEVER_DEFINED
int acol_manager::make_unique_name(string &col, std::vector<string> &cnames) {
  bool done;
  do {
    done=true;
    for(size_t i=0;i<cnames.size();i++) {
      if (col==cnames[i]) {
        done=false;
        i=cnames.size();
      }
    }
    if (done==false) {
      col+='_';
    }
  } while (done==false);
  return 0;
}
#endif

int acol_manager::comm_version(std::vector<std::string> &sv, bool itive_com) {

  cout << "\nacol: A data table viewing and processing program for O2scl."
       << endl;
  cout << (((string)" Compiled at ")+((string)__TIME__)+" on "+
	   ((string)__DATE__)+" for "+((string)PACKAGE)+", version "+
	   ((string)VERSION)+".\n") << endl;

  cout << "O2scl version: " << o2scl_settings.o2scl_version() << endl;
  cout << "Range checking: " << o2scl_settings.range_check() << endl;
  cout << "EOS library: " << o2scl_settings.eos_installed() << endl;
  cout << "Particle library: " << o2scl_settings.part_installed() << endl;
  cout << "HDF5 support: " << o2scl_settings.hdf_support() << endl;
  if (o2scl_settings.hdf_support()) {
    unsigned maj, min, rel;
    o2scl_settings.hdf5_header_version(maj,min,rel);
    cout << "  HDF5 version numbers during O2scl compilation: "
	 << maj << " " << min << " " << rel << endl;
    o2scl_settings.hdf5_lib_version(maj,min,rel);
    cout << "  HDF5 version numbers in libraries currently linked: "
	 << maj << " " << min << " " << rel << endl;
  }
  cout << "  HDF5 compression support: "
       << o2scl_settings.hdf5_compression_support() << endl;
  cout << "Armadillo support: " << o2scl_settings.armadillo_support() << endl;
  cout << "Eigen support: " << o2scl_settings.eigen_support() << endl;
  cout << "FFTW support: " << o2scl_settings.fftw_support() << endl;
  cout << "Python support: " << o2scl_settings.python_support() << endl;
  cout << "GSL V2.0+ support: " << o2scl_settings.gsl2_support() << endl;
  cout << "OpenMP support: " << o2scl_settings.openmp_support() << endl;
  cout << "Data directory: " << o2scl_settings.get_data_dir() << endl;
  cout << "System type: " << o2scl_settings.system_type() << endl;
  cout << endl;
  cout << "o2scl_name: " << o2scl_settings.o2scl_name() << endl;
  cout << "o2scl_package: " << o2scl_settings.o2scl_package() << endl;
  cout << "o2scl_bugreport: " << o2scl_settings.o2scl_bugreport() << endl;
  cout << "o2scl_string: " << o2scl_settings.o2scl_string() << endl;
  cout << "o2scl_tarname: " << o2scl_settings.o2scl_tarname() << endl;
  cout << endl;
  cout << "config.h: " << endl;
  o2scl_settings.config_h_report();
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

int acol_manager::comm_type(std::vector<std::string> &sv, 
			    bool itive_com) {
  if (type.length()==0) {
    cerr << "No current object to display type of." << endl;
    return 1;
  }
  cout << "Type is " << type << " ." << endl;
  return 0;
}

int acol_manager::comm_to_hist_2d(std::vector<std::string> &sv, 
				  bool itive_com) {

  std::string i1;

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
  } 

  cerr << "Cannot convert object of type " << type << " to histogram."
       << endl;
  
  return 1;
}

int acol_manager::comm_to_hist(std::vector<std::string> &sv, 
			       bool itive_com) {

  std::string i1;

  if (type=="table") {

    int ret=get_input_one(sv,((string)"Enter \"2d\" for 2d histogram ")+
			  +"and \"1d\" for 1d histogram",i1,"to-hist",
			  itive_com);
    if (ret!=0) return ret;
    
    vector<string> in, pr;
    pr.push_back("Column name");
    pr.push_back("Number of bins");
    ret=get_input(sv,pr,in,"to-hist",itive_com);
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

  cerr << "Cannot convert object of type " << type << " to histogram."
       << endl;
  
  return 1;
}

