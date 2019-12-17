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

void *o2scl_create_acol_manager() {
  o2scl_acol::acol_manager *amp=new o2scl_acol::acol_manager;
  // We call acol_manager::run() in o2scl_acol_set_names() later
  return amp;
}
  
void o2scl_free_acol_manager(void *vp) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  delete amp;
  return;
}

void o2scl_acol_set_names(void *vp, int n1, char *cmd_name,
			  int n2, char *short_desc, int n3,
			  char *env_var) {

  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  string str;
  for(int i=0;i<n3;i++) str+=env_var[i];
  amp->env_var_name=str;
  str.clear();

  // We call run() here because the environment variable
  // name is now properly set
  amp->run(0,0,false);

  // And now that run() has been called, the 'cl' pointer is valid
  for(int i=0;i<n1;i++) str+=cmd_name[i];
  amp->cl->cmd_name=str;
  str.clear();
  for(int i=0;i<n2;i++) str+=short_desc[i];
  amp->cl->desc=str;
  str.clear();
  
  return;
}


std::vector<std::string> o2scl_acol_parse_arrays
(int n_entries, int *sizes, char *str) {
  std::vector<std::string> list;
  size_t ix=0;
  for(int i=0;i<n_entries;i++) {
    std::string tmp;
    if (sizes[i]<0) {
      O2SCL_ERR("Size cannot be negative in o2scl_acol_parse_arrays().",
		o2scl::exc_einval);
    }
    for(int j=0;j<sizes[i];j++) {
      tmp+=str[ix];
      ix++;
    }
    list.push_back(tmp);
  }
  return list;
}

void o2scl_acol_parse(void *vp, int n_entries, int *sizes, 
		      char *str) {
  std::vector<std::string> args=o2scl_acol_parse_arrays(n_entries,sizes,str);
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  std::vector<o2scl::cmd_line_arg> ca;
  amp->cl->process_args(args,ca,0);
  amp->cl->call_args(ca);
  return;
}

void o2scl_acol_alias_counts(void *vp, int n_entries, int *sizes, 
			     char *str, int &n_new, int &s_new) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  /*
    cout << "Before: " << n_entries << endl;
    for(int k=0;k<n_entries;k++) {
    cout << k << " " << sizes[k] << endl;
    }
  */
  std::vector<std::string> args=o2scl_acol_parse_arrays(n_entries,sizes,str);
  /*
    cout << "Before: " << endl;
    for(size_t k=0;k<args.size();k++) {
    cout << "args[k]: " << k << " " << args[k] << endl;
    }
  */
  amp->cl->parse_for_aliases(args);
  amp->cl->apply_aliases(args,0);
  /*
    cout << "After: " << endl;
    for(size_t k=0;k<args.size();k++) {
    cout << "args[k]: " << k << " " << args[k] << endl;
    }
  */
  n_new=args.size();
  s_new=0;
  for(size_t k=0;k<args.size();k++) s_new+=args[k].length();
  //cout << "n_new, s_new: " << n_new << " " << s_new << endl;
  return;
}

void o2scl_acol_apply_aliases(void *vp, int n_entries, int *sizes, 
			      char *str, int *sizes_new,
			      char *str_new) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  std::vector<std::string> args=o2scl_acol_parse_arrays(n_entries,sizes,str);
  amp->cl->apply_aliases(args,0);
  int cnt=0;
  for(size_t i=0;i<args.size();i++) {
    sizes_new[i]=args[i].length();
    for(size_t j=0;j<args[i].length();j++) {
      str_new[cnt]=args[i][j];
      cnt++;
    }
  }
  return;
}

int o2scl_acol_get_column(void *vp, char *col_name,
			  int &n, double *&ptr) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="table") {
    return 1;
  }
  n=amp->table_obj.get_nlines();
  std::string stmp=col_name;
  if (amp->table_obj.is_column(stmp)==false) {
    return 2;
  }
  const std::vector<double> &col=amp->table_obj.get_column(stmp);
  ptr=(double *)&col[0];
  return 0;
}

int o2scl_acol_pdma_get_base(void *vp, int &ndim, int &n,
			     double *&low, double *&high) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="prob_dens_mdim_amr") {
    return 1;
  }
  prob_dens_mdim_amr<> &pdma=amp->pdma_obj;
  ndim=pdma.n_dim;
  n=pdma.mesh.size();
  low=&(pdma.low[0]);
  high=&(pdma.high[0]);
  return 0;
}

int o2scl_acol_pdma_get_cube(void *vp, int ix, 
			     double *&low, double *&high,
			     double &frac_vol, double &weight) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="prob_dens_mdim_amr") {
    return 1;
  }
  prob_dens_mdim_amr<> &pdma=amp->pdma_obj;
  low=&(pdma.mesh[ix].low[0]);
  high=&(pdma.mesh[ix].high[0]);
  frac_vol=pdma.mesh[ix].frac_vol;
  weight=pdma.mesh[ix].weight;
  return 0;
}

int o2scl_acol_get_row_ser(void *vp, char *pattern, int row_index,
			   int &n, double *&ptr) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="table") {
    return 1;
  }
  amp->doublev_obj.clear();
  for(size_t j=0;j<amp->table_obj.get_ncolumns();j++) {
    if (fnmatch(pattern,amp->table_obj.get_column_name(j).c_str(),0)==0) {
      amp->doublev_obj.push_back(amp->table_obj.get(j,row_index));
    }
  }
  n=amp->doublev_obj.size();
  ptr=(double *)&amp->doublev_obj[0];
  return 0;
}

int o2scl_acol_get_double_arr(void *vp, int &n, double *&ptr) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type=="double[]") {
    n=amp->doublev_obj.size();
  } else if (amp->type=="int[]") {
    n=amp->intv_obj.size();
    amp->doublev_obj.resize(n);
    for(int i=0;i<n;i++) {
      amp->doublev_obj[i]=amp->intv_obj[i];
    }
  } else if (amp->type=="size_t[]") {
    n=amp->size_tv_obj.size();
    amp->doublev_obj.resize(n);
    for(int i=0;i<n;i++) {
      amp->doublev_obj[i]=amp->size_tv_obj[i];
    }
  } else {
    return 1;
  }
  ptr=&(amp->doublev_obj[0]);
  return 0;
}

int o2scl_acol_get_hist_reps(void *vp, int &n, double *&ptr) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="hist") return 1;
  n=amp->hist_obj.size();
  amp->xtemp.resize(n);
  for(int i=0;i<n;i++) amp->xtemp[i]=amp->hist_obj.get_rep_i(i);
  ptr=&(amp->xtemp[0]);
  return 0;
}

int o2scl_acol_get_tensor_grid3(void *vp, int &nx, int &ny,
				int &nz, const double *&xg, const double *&yg,
				const double *&zg, const double *&data) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type=="tensor_grid") {
    if (amp->tensor_grid_obj.get_rank()!=3) return 2;
    nx=amp->tensor_grid_obj.get_size(0);
    ny=amp->tensor_grid_obj.get_size(1);
    nz=amp->tensor_grid_obj.get_size(2);
    const std::vector<double> &grid=amp->tensor_grid_obj.get_grid();
    xg=&(grid[0]);
    yg=&(grid[nx]);
    zg=&(grid[nx+ny]);
    data=&(amp->tensor_grid_obj.get_data()[0]);
    return 0;
  }
  return 1;
}

int o2scl_acol_mult_vectors_to_conts(void *vp, char *str1,
				     char *str2) {
  
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  
  vector<vector<double> > v1, v2;
  string s1, s2;
  if (str1!=0) {
    s1=str1;
    int ret1=mult_vector_spec(s1,v1,amp->get_verbose(),false);
    if (ret1!=0) return ret1;
  }
  
  s2=str2;
  int ret2=mult_vector_spec(s2,v2,amp->get_verbose(),false);
  if (ret2!=0) return ret2;
  
  amp->cont_obj.clear();

  if (str1!=0) {
    for(size_t i=0;i<v1.size();i++) {
      for(size_t j=0;j<v2.size();j++) {
	size_t ntemp=v1[i].size();
	if (v2[j].size()<ntemp) ntemp=v2[j].size();
	contour_line cl;
	cl.level=0.0;
	cl.x.resize(ntemp);
	cl.y.resize(ntemp);
	for(size_t k=0;k<ntemp;k++) {
	  cl.x[k]=v1[i][k];
	  cl.y[k]=v2[j][k];
	}
	amp->cont_obj.push_back(cl);
      }
    }
  } else {
    for(size_t j=0;j<v2.size();j++) {
      size_t ntemp=v2[j].size();
      contour_line cl;
      cl.level=0.0;
      cl.x.resize(ntemp);
      cl.y.resize(ntemp);
      for(size_t k=0;k<ntemp;k++) {
	cl.x[k]=k;
	cl.y[k]=v2[j][k];
      }
      amp->cont_obj.push_back(cl);
    }
  }
  return 0;
}

int o2scl_acol_get_hist_wgts(void *vp, int &n, double *&ptr) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="hist") return 1;
  n=amp->hist_obj.size();
  amp->ytemp.resize(n);
  for(int i=0;i<n;i++) amp->ytemp[i]=amp->hist_obj.get_wgt_i(i);
  ptr=&(amp->ytemp[0]);
  return 0;
}

int o2scl_acol_contours_n(void *vp) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  return amp->cont_obj.size();
}

double o2scl_acol_contours_line(void *vp, int i, int &n, double *&ptrx,
				double *&ptry) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  double lev=amp->cont_obj[i].level;
  n=amp->cont_obj[i].x.size();
  ptrx=&(amp->cont_obj[i].x[0]);
  ptry=&(amp->cont_obj[i].y[0]);
  return lev;
}

void o2scl_acol_get_type(void *vp, int &n, char *&str) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  n=amp->type.length();
  if (n>0) {
    str=(char *)(amp->type.c_str());
  } else {
    str=0;
  }
  return;
}

int o2scl_acol_tensor_to_table3d(void *vp, int i1, int i2) {

  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  
  if (amp->type!="tensor" && amp->type!="tensor_grid" &&
      amp->type!="tensor<size_t>" && amp->type!="tensor<int>") {
    return 2;
  }

  if (i1*i2!=0 || i1+i2!=1) {
    cerr << "Indices malformed in o2scl_acol_tensor_to_table3d." << endl;
    return 3;
  }
  
  amp->table3d_obj.clear();
  
  if (amp->type=="tensor") {
    if (amp->tensor_obj.get_rank()!=2) {
      cerr << "Can only autoconvert to table3d if rank is 2." << endl;
      return 1;
    }
    vector<double> grid_x, grid_y;
    for(size_t i=0;i<amp->tensor_obj.get_size(i1);i++) {
      grid_x.push_back((double)i);
    }
    for(size_t i=0;i<amp->tensor_obj.get_size(i2);i++) {
      grid_y.push_back((double)i);
    }
    amp->table3d_obj.set_xy("x",grid_x.size(),grid_x,
		       "y",grid_y.size(),grid_y);
    amp->table3d_obj.new_slice("tensor");
    vector<size_t> ix(2);
    for(size_t i=0;i<amp->table3d_obj.get_nx();i++) {
      for(size_t j=0;j<amp->table3d_obj.get_nx();j++) {
	ix[i1]=i;
	ix[i2]=j;
	amp->table3d_obj.set(i,j,"tensor",amp->tensor_obj.get(ix));
      }
    }
  } else if (amp->type=="tensor<size_t>") {
    if (amp->tensor_size_t_obj.get_rank()!=2) {
      cerr << "Can only autoconvert to table3d if rank is 2." << endl;
      return 1;
    }
    vector<double> grid_x, grid_y;
    for(size_t i=0;i<amp->tensor_size_t_obj.get_size(i1);i++) {
      grid_x.push_back((double)i);
    }
    for(size_t i=0;i<amp->tensor_size_t_obj.get_size(i2);i++) {
      grid_y.push_back((double)i);
    }
    amp->table3d_obj.set_xy("x",grid_x.size(),grid_x,
		       "y",grid_y.size(),grid_y);
    amp->table3d_obj.new_slice("tensor");
    vector<size_t> ix(2);
    for(size_t i=0;i<amp->table3d_obj.get_nx();i++) {
      for(size_t j=0;j<amp->table3d_obj.get_nx();j++) {
	ix[i1]=i;
	ix[i2]=j;
	amp->table3d_obj.set(i,j,"tensor",amp->tensor_size_t_obj.get(ix));
      }
    }
  } else if (amp->type=="tensor<int>") {
    if (amp->tensor_int_obj.get_rank()!=2) {
      cerr << "Can only autoconvert to table3d if rank is 2." << endl;
      return 1;
    }
    vector<double> grid_x, grid_y;
    for(size_t i=0;i<amp->tensor_int_obj.get_size(i1);i++) {
      grid_x.push_back((double)i);
    }
    for(size_t i=0;i<amp->tensor_int_obj.get_size(i2);i++) {
      grid_y.push_back((double)i);
    }
    amp->table3d_obj.set_xy("x",grid_x.size(),grid_x,
		       "y",grid_y.size(),grid_y);
    amp->table3d_obj.new_slice("tensor");
    vector<size_t> ix(2);
    for(size_t i=0;i<amp->table3d_obj.get_nx();i++) {
      for(size_t j=0;j<amp->table3d_obj.get_nx();j++) {
	ix[i1]=i;
	ix[i2]=j;
	amp->table3d_obj.set(i,j,"tensor",amp->tensor_int_obj.get(ix));
      }
    }
  } else if (amp->type=="tensor_grid") {
    if (amp->tensor_grid_obj.get_rank()!=2) {
      cerr << "Can only autoconvert to table3d if rank is 2." << endl;
      return 1;
    }
    vector<double> grid_x, grid_y;
    for(size_t i=0;i<amp->tensor_grid_obj.get_size(i1);i++) {
      grid_x.push_back(amp->tensor_grid_obj.get_grid(i1,i));
    }
    for(size_t i=0;i<amp->tensor_grid_obj.get_size(i2);i++) {
      grid_y.push_back(amp->tensor_grid_obj.get_grid(i2,i));
    }
    amp->table3d_obj.set_xy("x",grid_x.size(),grid_x,
			    "y",grid_y.size(),grid_y);
    amp->table3d_obj.new_slice("tensor");
    vector<size_t> ix(2);
    for(size_t i=0;i<amp->table3d_obj.get_nx();i++) {
      for(size_t j=0;j<amp->table3d_obj.get_ny();j++) {
	ix[i1]=i;
	ix[i2]=j;
	amp->table3d_obj.set(i,j,"tensor",amp->tensor_grid_obj.get(ix));
      }
    }
  }

  return 0;
}

int o2scl_acol_get_slice(void *vp, char *slice_name,
			 int &nx, double *&xptr,
			 int &ny, double *&yptr,
			 double *&data) {

  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;

  // Comment this out because this function is used in
  // den-plot to plot tensor types and table3d_obj is
  // used for empty storage
  /*
    if (amp->type!="table3d") {
    return 1;
    }
  */

  nx=amp->table3d_obj.get_nx();
  amp->xtemp.resize(nx);
  o2scl::vector_copy(amp->table3d_obj.get_x_data(),amp->xtemp);
  xptr=(double *)&amp->xtemp[0];

  ny=amp->table3d_obj.get_ny();
  amp->ytemp.resize(ny);
  o2scl::vector_copy(amp->table3d_obj.get_y_data(),amp->ytemp);
  yptr=(double *)&amp->ytemp[0];

  amp->stemp.resize(nx*ny);
  std::string stmp=slice_name;
  size_t itmp;
  if (!amp->table3d_obj.is_slice(stmp,itmp)) {
    return 2;
  }
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  const ubmatrix &m=amp->table3d_obj.get_slice(stmp);
  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      amp->stemp[i*ny+j]=m(i,j);
    }
  }
  data=(double *)&amp->stemp[0];
  
  return 0;
}
  
int o2scl_acol_get_hist_2d(void *vp, 
			   int &nx, double *&xptr,
			   int &ny, double *&yptr,
			   double *&data) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="hist_2d") {
    return 1;
  }

  nx=amp->hist_2d_obj.size_x();
  amp->xtemp.resize(nx);
  for(int i=0;i<nx;i++) {
    amp->xtemp[i]=amp->hist_2d_obj.get_x_rep_i(i);
  }
  xptr=(double *)&amp->xtemp[0];

  ny=amp->hist_2d_obj.size_y();
  amp->ytemp.resize(ny);
  for(int i=0;i<ny;i++) {
    amp->ytemp[i]=amp->hist_2d_obj.get_y_rep_i(i);
  }
  yptr=(double *)&amp->ytemp[0];

  amp->stemp.resize(nx*ny);
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  const ubmatrix &m=amp->hist_2d_obj.get_wgts();
  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      amp->stemp[i*ny+j]=m(i,j);
    }
  }
  data=(double *)&amp->stemp[0];
  return 0;
}

