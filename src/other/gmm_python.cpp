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
#include <o2scl/gmm_python.h>

using namespace std;
using namespace o2scl;

gmm_python::gmm_python() {
  p_set_func=0;
  p_components_func=0;
  p_set_args=0;
  p_components_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_points=0;
  n_components=0;
}
    
/** \brief Specify the Python module and function
 */
gmm_python::gmm_python(std::string module, size_t n_comp,
                       const o2scl::tensor<> &params,
                       std::string options, 
                       std::string class_name, int v) {
                    
  verbose=v;

  if (o2scl_settings.py_initialized==false) {
    if (verbose>1) {
      std::cout << "Running py_init()." << std::endl;
    }
    o2scl_settings.py_init();
  }
  p_set_func=0;
  p_components_func=0;
  p_get_func=0;
  p_set_args=0;
  p_components_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_points=0;
  n_components=0;

  if (module.length()>0) {
    set_function(module,n_comp,params,options,class_name,v);
  }
}      
    
void gmm_python::free() {
  if (verbose>1) {
    std::cout << "Starting gmm_python::free()." << std::endl;
  }
  if (p_set_func!=0) {
    if (verbose>1) {
      std::cout << "Decref set_func." << std::endl;
    }
    Py_DECREF(p_set_func);
  }
  if (p_components_func!=0) {
    if (verbose>1) {
      std::cout << "Decref components_func." << std::endl;
    }
    Py_DECREF(p_components_func);
  }
  if (p_get_func!=0) {
    if (verbose>1) {
      std::cout << "Decref get_func." << std::endl;
    }
    Py_DECREF(p_get_func);
  }
  if (p_set_args!=0) {
    if (verbose>1) {
      std::cout << "Decref set_args." << std::endl;
    }
    Py_DECREF(p_set_args);
  }
  if (p_components_args!=0) {
    if (verbose>1) {
      std::cout << "Decref components_args." << std::endl;
    }
    Py_DECREF(p_components_args);
  }
  if (p_instance!=0) {
    if (verbose>1) {
      std::cout << "Decref instance." << std::endl;
    }
    Py_DECREF(p_instance);
  }
  if (p_class!=0) {
    if (verbose>1) {
      std::cout << "Decref class." << std::endl;
    }
    Py_DECREF(p_class);
  }
  if (p_module!=0) {
    if (verbose>1) {
      std::cout << "Decref module." << std::endl;
    }
    Py_DECREF(p_module);
  }
  if (p_name!=0) {
    if (verbose>1) {
      std::cout << "Decref name." << std::endl;
    }
    Py_DECREF(p_name);
  }

  p_set_func=0;
  p_components_func=0;
  p_components_args=0;
  p_instance=0;
  p_class=0;
  p_module=0;
  p_name=0;
      
  n_params=0;
  n_points=0;
      
  if (verbose>1) {
    std::cout << "Done in gmm_python::free()." << std::endl;
  }

  set_func="set_data_str";
  get_func="get_data";
  components_func="components";
  log_pdf_func="log_pdf";
}      
    
gmm_python::~gmm_python() {
  free();
}      

int gmm_python::set_function(std::string module, 
                             size_t n_comp, const o2scl::tensor<> &params,
                             std::string options, std::string class_name, 
                             int v) {
  int ret;
  void *vp=set_function_internal(module,n_comp,params,options,
                                 class_name,v,ret);
                                 
  return ret;
}

void *gmm_python::set_function_internal
(std::string module, size_t n_comp, const o2scl::tensor<> &params,
 std::string options, std::string class_name, int v, int &ret) {
  
  ret=0;
  
  if (params.get_rank()!=2) {
    O2SCL_ERR2("Invalid rank for input tensors in ",
               "gmm_python().",o2scl::exc_einval);
  }
      
  free();

  n_params=params.get_size(1);
  n_points=params.get_size(0);
  n_components=n_comp;
      
  if (options.length()>0) {
    options+=",n_components="+o2scl::szttos(n_comp);
  } else {
    options="n_components="+o2scl::szttos(n_comp);
  }
      
  // Get the Unicode name of the user-specified module
  if (verbose>1) {
    std::cout << "Python version: "
              << o2scl_settings.py_version() << std::endl;
    std::cout << "Staring gmm_python::set_function()."
              << std::endl;
    std::cout << "  Getting unicode for module named "
              << module << std::endl;
  }
  p_name=PyUnicode_FromString(module.c_str());
  if (p_name==0) {
    O2SCL_ERR2("Create module name failed in ",
               "gmm_python::set_function().",
               o2scl::exc_efailed);
  }
      
  // Import the user-specified module
  if (verbose>1) {
    std::cout << "  Importing module." << std::endl;
  }
  p_module=PyImport_Import(p_name);
  if (p_module==0) {
    O2SCL_ERR2("Load module failed in ",
               "gmm_python::set_function().",
               o2scl::exc_efailed);
  }

  if (class_name.length()>0) {
    if (verbose>1) {
      std::cout << "  Obtaining python class." << std::endl;
    }
    p_class=PyObject_GetAttrString(p_module,class_name.c_str());
    if (p_class==0) {
      O2SCL_ERR2("Get class failed in ",
                 "emulator_python::set().",o2scl::exc_efailed);
    }
        
    // Create an instance of the class
    if (verbose>1) {
      std::cout << "  Loading python class." << std::endl;
    }
    if (PyCallable_Check(p_class)==false) {
      O2SCL_ERR2("Check class callable failed in ",
                 "funct_python_method::set_function().",
                 o2scl::exc_efailed);
    }
        
    if (verbose>1) {
      std::cout << "  Loading python class instance." << std::endl;
    }
    p_instance=PyObject_CallObject(p_class,0);
    if (p_instance==0) {
      O2SCL_ERR2("Instantiate class failed in ",
                 "funct_python_method::set_function().",
                 o2scl::exc_efailed);
    }
  }
      
  // Setup the arguments to the python function
  if (verbose>1) {
    std::cout << "  Making argument object for set function."
              << std::endl;
  }
  p_set_args=PyTuple_New(2);
  if (p_set_args==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "gmm_python::set_function().",
               o2scl::exc_efailed);
  }

  // Setup the arguments to the python function
  if (verbose>1) {
    std::cout << "  Making argument object for components function."
              << std::endl;
  }
  p_components_args=PyTuple_New(1);
  if (p_components_args==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "gmm_python::set_function().",
               o2scl::exc_efailed);
  }

  if (class_name.length()>0) {

    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python member function components: "
                << components_func<< std::endl;
    }
    p_components_func=PyObject_GetAttrString(p_instance,
                                             components_func.c_str());
    if (p_components_func==0) {
      O2SCL_ERR2("Get components function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }
        
    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python member function get: "
                << get_func<< std::endl;
    }
    p_get_func=PyObject_GetAttrString(p_instance,get_func.c_str());
    if (p_get_func==0) {
      O2SCL_ERR2("Get get function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }
        
    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python member function set: "
                << set_func << std::endl;
    }
    p_set_func=PyObject_GetAttrString(p_instance,set_func.c_str());
    if (p_set_func==0) {
      O2SCL_ERR2("Get set function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }

  } else {
    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python function set." << std::endl;
    }
    p_set_func=PyObject_GetAttrString(p_module,set_func.c_str());
    if (p_set_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }

    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python function components." << std::endl;
    }
    p_components_func=PyObject_GetAttrString(p_module,components_func.c_str());
    if (p_components_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }

    // Load the python function
    if (verbose>1) {
      std::cout << "  Loading python function get." << std::endl;
    }
    p_get_func=PyObject_GetAttrString(p_module,get_func.c_str());
    if (p_get_func==0) {
      O2SCL_ERR2("Get function failed in ",
                 "gmm_python::set_function().",
                 o2scl::exc_efailed);
    }
  }

  // AWS, 2/21/23: I'm not sure why it has to be done here and not in
  // a different function, but if I don't do it here I get a seg fault.
  //void *vp=o2scl_settings.py_import_array();
  import_array();

  npy_intp params_dims[]={(npy_intp)params.get_size(0),
    (npy_intp)params.get_size(1)};
  if (verbose>1) {
    std::cout << "gmm_python::operator():" << std::endl;
  }
  PyObject *array_in=PyArray_SimpleNewFromData
    (2,params_dims,NPY_DOUBLE,(void *)(&(params.get_data()[0])));
         
  int pret=PyTuple_SetItem(p_set_args,0,array_in);
  if (pret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }
      
  if (verbose>1) {
    std::cout << "Creating python unicode for string: "
              << options.length() << " " << options << std::endl;
  }
  PyObject *p_options=PyUnicode_FromString(options.c_str());
  if (p_options==0) {
    O2SCL_ERR2("String creation failed in ",
               "emulator_python::set().",o2scl::exc_efailed);
  }
      
  int ret3=PyTuple_SetItem(p_set_args,1,p_options);
  if (ret3!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }

  // Call the python function
  if (verbose>1) {
    std::cout << "  Calling python set function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_set_func,p_set_args);
  if (result==0) {
    O2SCL_ERR2("Function set call failed in ",
               "gmm_python::operator().",o2scl::exc_efailed);
  }

  if (verbose>1) {
    std::cout << p_set_func << " " << p_components_func << " "
              << p_get_func << std::endl;
    std::cout << "Done with gmm_python::set_function()."
              << std::endl;
  }
      
  return 0;
}

int gmm_python::components(const std::vector<double> &x,
                           std::vector<double> &y) const {

  if (x.size()!=n_params) {
    O2SCL_ERR("Input vector does not have correct size.",
              o2scl::exc_einval);
  }
  
  if (p_set_func==0 || p_components_func==0) {
    O2SCL_ERR2("No functions found in ",
               "gmm_python::operator().",
               o2scl::exc_efailed);
  }

  npy_intp x_dims[]={(npy_intp)x.size()};
  if (verbose>1) {
    std::cout << "gmm_python::operator():" << std::endl;
    std::cout << "  Array x: " << x.size() << std::endl;
  }
  PyObject *array_x=PyArray_SimpleNewFromData
    (1,x_dims,NPY_DOUBLE,(void *)(&(x[0])));
      
  int ret=PyTuple_SetItem(p_components_args,0,array_x);
  if (ret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "mm_funct_python::operator().",o2scl::exc_efailed);
  }
      
  // Call the python function
  if (verbose>1) {
    std::cout << "  Calling python components function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_components_func,p_components_args);
  if (result==0) {
    O2SCL_ERR2("Function components call failed in ",
               "gmm_python::operator().",o2scl::exc_efailed);
  }

  if (PyArray_Check(result)==0) {
    O2SCL_ERR2("Function call did not return a numpy array in ",
               "gmm_python::operator().",o2scl::exc_efailed);
  }
      
  if (verbose>1) {
    std::cout << "  Obtaining output 1." << std::endl;
  }
  void *vp=PyArray_DATA((PyArrayObject *)result);
  double *dp=(double *)vp;
  for(size_t i=0;i<1;i++) {
    y[i]=dp[i];
    std::cout << "  i,y[i]: " << i << " " << y[i] << std::endl;
  }
      
  if (verbose>1) {
    std::cout << "  Decref result." << std::endl;
  }
  Py_DECREF(result);
  
  if (verbose>1) {
    std::cout << "Done in gmm_python::operator()."
              << std::endl;
  }

  return 0;
}      

int gmm_python::get_python() {
  int ret;
  void *vp=get_python_internal(ret);
  return ret;
}

void *gmm_python::get_python_internal(int &ret) {

  ret=0;
  
  import_array();
  
  if (p_set_func==0 || p_components_func==0 || p_get_func==0) {
    O2SCL_ERR2("No functions found in ",
               "gmm_python::operator().",
               o2scl::exc_efailed);
  }

  // Call the python function
  if (verbose>1) {
    std::cout << "  Calling python get function." << std::endl;
  }
  PyObject *result=PyObject_CallObject(p_get_func,0);
  if (result==0) {
    O2SCL_ERR2("Function get call failed in ",
               "gmm_python::operator().",o2scl::exc_efailed);
  }

  if (PyTuple_Check(result)==0) {
    O2SCL_ERR2("Function call did not return a tuple in ",
               "gmm_python::operator().",o2scl::exc_efailed);
  }
      
  if (verbose>1) {
    std::cout << "  Obtaining output 2." << std::endl;
  }

  PyObject *get_w, *get_m, *get_c, *get_p, *get_pc;
  
  get_w=PyTuple_GetItem(result,0);
  get_m=PyTuple_GetItem(result,1);
  get_c=PyTuple_GetItem(result,2);
  get_p=PyTuple_GetItem(result,3);
  get_pc=PyTuple_GetItem(result,4);
  if (get_w==0 || get_m==0 || get_c==0 || get_p==0 || get_pc==0) {
    O2SCL_ERR("Get tuple failed",o2scl::exc_einval);
  }
  
  if (PyArray_Check(get_w)==0) {
    O2SCL_ERR2("Function call did not return a numpy array 1 in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  if (PyArray_Check(get_m)==0) {
    O2SCL_ERR2("Function call did not return a numpy array 2 in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  if (PyArray_Check(get_c)==0) {
    O2SCL_ERR2("Function call did not return a numpy array 3 in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  if (PyArray_Check(get_p)==0) {
    O2SCL_ERR2("Function call did not return a numpy array 4 in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  if (PyArray_Check(get_pc)==0) {
    O2SCL_ERR2("Function call did not return a numpy array 5 in ",
               "interpm_python::operator().",o2scl::exc_efailed);
  }
  
  if (verbose>1) {
    std::cout << "  Obtaining output." << std::endl;
  }
  double *ptr_w=(double *)PyArray_DATA((PyArrayObject *)get_w);
  double *ptr_m=(double *)PyArray_DATA((PyArrayObject *)get_m);
  double *ptr_c=(double *)PyArray_DATA((PyArrayObject *)get_c);
  double *ptr_p=(double *)PyArray_DATA((PyArrayObject *)get_p);
  double *ptr_pc=(double *)PyArray_DATA((PyArrayObject *)get_pc);
  
  pdm_gmm.weights.resize(n_components);
  pdm_gmm.pdmg.resize(n_components);
  
  for(size_t i=0;i<n_components;i++) {
    pdm_gmm.weights[i]=ptr_w[i];
    if (verbose>1) {
      std::cout << "Component " << i << " with weight: "
                << pdm_gmm.weights[i] << std::endl;
    }
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    ubvector mean(n_params);
    ubmatrix covar(n_params,n_params);
    for(size_t j=0;j<n_params;j++) {
      mean[j]=ptr_m[n_params*i+j];
      if (verbose>1) {
        std::cout << "mean[" << j << "]: " << mean[j] << std::endl;
      }
      for(size_t k=0;k<n_params;k++) {
        covar(j,k)=ptr_c[n_params*n_params*i+n_params*j+k];
        if (verbose>1) {
          std::cout << "covar(" << j << "," << k << "): " << covar(j,k)
                    << std::endl;
        }
      }
    }
    //pdm_gmm.pdmg[i].verbose=2;
    
    // AWS, 2/25/23: Unfortunately sklearn handles it's Gaussian
    // distributions a bit differently than O2scl. O2scl uses the
    // Cholesky decomposition of the covariance matrix, while sklearn
    // uses the Cholesky decomposition of the "precisions" matrix,
    // which is the inverse of the covariance matrix. We recompute the
    // Gaussians here, but there is probably a faster way.
    if (verbose>1) {
      pdm_gmm.pdmg[i].verbose=1;
    }
    pdm_gmm.pdmg[i].set_covar(n_params,mean,covar);

    if (verbose>1) {
      cout << endl;
    }
  }
  
  if (verbose>1) {
    std::cout << "  Decref result." << std::endl;
  }
  Py_DECREF(result);
  
  if (verbose>1) {
    std::cout << "Done in gmm_python::get_python()."
              << std::endl;
  }

  return 0;
}      

#ifdef O2SCL_PYTHON

int o2scl::compare_methast_gmm_kde(size_t n_comp_start, size_t n_comp_end,
                                   const o2scl::table<> &table,
                                   std::vector<std::string> param_cols,
                                   std::string lw_col,
                                   std::vector<gmm_python> &gp,
                                   kde_python<> &kp,
                                   std::vector<double> bw_array,
                                   double test_size, size_t n_tests,
                                   int verbose, std::string gp_options,
                                   std::string kp_options) {

  typedef boost::numeric::ublas::vector<double> ubvector;
  
  // Collect sizes
  size_t n_params=param_cols.size();
  size_t n_dat=table.get_nlines();
  size_t n_test=n_dat*test_size;
  size_t n_train=n_dat-n_test;
  if (n_tests==0) {
    n_tests=n_test*(n_test-1)/2;
  }

  // Allocate space for train and test tensor
  o2scl::tensor<> dat_train, dat_test;
  size_t sz[2];
  sz[0]=n_train;
  sz[1]=n_params;
  dat_train.resize(2,sz);
  sz[0]=n_test;
  dat_test.resize(2,sz);

  // Shuffle 
  rng<> r;
  r.clock_seed();
  std::vector<size_t> index(n_dat);
  for(size_t i=0;i<index.size();i++) index[i]=i;
  vector_shuffle<std::vector<size_t>,size_t>(r,index.size(),index);

  // Copy data from table to train and test tensors
  for(size_t i=0;i<index.size();i++) {
    size_t ix_src[2], ix[2];
    if (i<n_train) {
      for(size_t j=0;j<n_params;j++) {
        ix[0]=i;
        ix[1]=j;
        dat_train.get(ix)=table.get(param_cols[j],i);
      }
    } else {
      for(size_t j=0;j<n_params;j++) {
        ix[0]=i-n_train;
        ix[1]=j;
        dat_test.get(ix)=table.get(param_cols[j],i);
      }
    }
  }

  // Number of models
  size_t n_models=n_comp_end-n_comp_start+2;
    
  // Store the average deviation in the log likelihood for all of the
  // GMM models and the KDE model
  std::vector<double> avgs(n_models);

  // Train all of the models
  for(size_t i=0;i<n_models;i++) {
    if (i==0) {
      kp.set_function("o2sclpy",dat_train,bw_array,
                      kp_options,"kde_sklearn"); 
    } else {
      size_t n_comp=n_comp_start+i-1;
      gp[i-1].set_function("o2sclpy",n_comp,
                           dat_train,gp_options,"gmm_sklearn");
      gp[i-1].get_python();
    }
  }
    
  // Loop over the requested number of tests
  for(size_t i=0;i<n_tests;i++) {
    size_t j_test=r.random_int(n_test);
    size_t j_test2=index[j_test]+n_train;
    size_t k_test=r.random_int(n_test);
    size_t k_test2=index[k_test]+n_train;

    // Loop over each model
    for(size_t j=0;j<n_models;j++) {

      std::vector<double> x_j(n_params), x_k(n_params);
      ubvector x_j2(n_params), x_k2(n_params);
      for(size_t k=0;k<n_params;k++) {
        x_j[k]=table.get(param_cols[k],j_test2);
        x_k[k]=table.get(param_cols[k],k_test2);
        x_j2[k]=table.get(param_cols[k],j_test2);
        x_k2[k]=table.get(param_cols[k],k_test2);
      }

      // The exact difference in the log weight
      double log_wgt_j=table.get(lw_col,j_test2);
      double log_wgt_k=table.get(lw_col,k_test2);
      double diff_exact=log_wgt_j-log_wgt_k;
        
        
      if (j==0) {
        double lp_j=kp.log_pdf(x_j);
        double lp_k=kp.log_pdf(x_k);
        // The estimated difference in the log weight from the KDE
        double diff_est=lp_j-lp_k;
        avgs[j]+=pow(diff_est-diff_exact,2.0)/pow(diff_exact,2.0);
      } else {
        double lp_j=gp[j-1].get_gmm().log_pdf(x_j2);
        double lp_k=gp[j-1].get_gmm().log_pdf(x_k2);
        // The estimated difference in the log weight from the GMM
        double diff_est=lp_j-lp_k;
        avgs[j]+=pow(diff_est-diff_exact,2.0)/pow(diff_exact,2.0);
      }
    }
  }

  return vector_min_index<std::vector<double>,double>
    (avgs.size(),avgs);
}

#endif
