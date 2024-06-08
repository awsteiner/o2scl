 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2022-2024, Mahamudul Hasan Anik, Satyajit Roy, and
  Andrew W. Steiner
  
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
/** \file emulator.h
 */
#ifndef O2SCL_EMULATOR_H
#define O2SCL_EMULATOR_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/set_python.h>

#ifdef O2SCL_SET_PYTHON
#include <Python.h>
#endif

#include <o2scl/interpm_idw.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/vec_stats.h>

namespace o2scl {
  
  /** \brief Emulator base class
   */
  template<class data_t,
           class vec_t=boost::numeric::ublas::vector<double>>
  class emulator_base {
    
  protected:
    
  public:
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat
    */
    virtual int eval(size_t n, const vec_t &p, double &log_wgt,
                     data_t &dat)=0;
    
  };

  /** \brief Emulator with uncertainty base class
   */
  template<class data_t, class data_unc_t,
           class vec_t=boost::numeric::ublas::vector<double>>
  class emulator_unc : public emulator_base<data_t,vec_t> {
    
  public:
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                 double &log_wgt_unc, data_t &dat, data_unc_t &dat_unc)=0;
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat
    */
    virtual int eval(size_t n, const vec_t &p, double &log_wgt,
             data_t &dat) {
      double log_wgt_unc;
      // This assignment effectively allocates memory for
      // a new dat_unc object using a copy constructor
      data_unc_t dat_unc=dat;
      return eval_unc(n,p,log_wgt,log_wgt_unc,dat,dat_unc);
    }
    
  };
  
  /** \brief Emulate data stored in a table object with interpm_idw

      \note Currently, this emulator only works if the data object
      from the MCMC class is a vector type so that it knows how
      to find the "log weight" output
      
      This class is experimental.
  */
  template<class vec2_t=boost::numeric::ublas::vector<double>,
           class vec_t=boost::numeric::ublas::vector<double>>
  class emulator_interpm_idw_table :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:
    
    /// The view of the user-specified table
    const_matrix_view_table_transpose<> cmvtt;
    
    /// Index of the "log weight" in the MCMC data vector
    size_t ix;
    
  public:
    
    /// The internal interpolation object
    o2scl::interpm_idw<o2scl::const_matrix_view_table_transpose<>> ii;
    
    /** \brief Create an emulator
     */
    emulator_interpm_idw_table() {
    }
    
    /** \brief Set the emulator

        Set the emulator using a table containing \c np parameters and
        \c n_out output quantities. The variable \c ix_log_wgt should
        be the index of the log_weight among all of the output
        variables, from 0 to <tt>n_out-1</tt>. The list, \c list,
        should include the column names of the parameters and then the
        output quantities (including the log weight column), in order.
     */
    void set(size_t np, size_t n_out, size_t ix_log_wgt,
             table<> &t, std::vector<std::string> list) {
      if (ix_log_wgt>list.size()-np) {
        O2SCL_ERR2("Index of log_wgt is too large in ",
                   "emulator_interpm_idw_table::set().",o2scl::exc_einval);
      }
      cmvtt.set(t,list);
      ix=ix_log_wgt;
      ii.set_data(np,n_out,t.get_nlines(),cmvtt);
      return;
    }
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                 double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
      
      ii.eval_err<vec_t,vec2_t,vec2_t>(p,dat,dat_unc);
      log_wgt=dat[ix];
      log_wgt_unc=dat_unc[ix];
      return 0;
    }
    
  };

  /** \brief Emulate data stored in a table object with interpm_krige

      \note Currently, this emulator only works if the data object
      from the MCMC class is a vector type so that it knows how
      to find the "log weight" output

      This class is experimental.
   */
  template<class vec2_t=boost::numeric::ublas::vector<double>,
           class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::matrix<double>,
           class mat_inv_t=o2scl_linalg::matrix_invert_det_cholesky<
             boost::numeric::ublas::matrix<double> >>
  class emulator_interpm_krige_table :
    public emulator_unc<vec2_t,vec2_t,vec_t> {
    
  protected:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef o2scl::matrix_view_table<> mat_x_t;
    typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
    typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
    typedef o2scl::matrix_view_table_transpose<> mat_y_t;
    typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
    
    /// Index of the "log weight" in the MCMC data vector
    size_t ix;

  public:
    
    /// The internal interpolation object
    interpm_krige_optim
    <std::vector<mcovar_funct_rbf_noise>,ubvector,mat_x_t,
     mat_x_row_t,mat_y_t,mat_y_row_t,mat_t,mat_inv_t> iko;

    /// The covariance functions
    std::vector<mcovar_funct_rbf_noise> mfrn;
    
    /** \brief Create an emulator
     */
    emulator_interpm_krige_table() {
    }
    
    /** \brief Set the emulator

        Set the emulator using a table containing \c np parameters and
        \c n_out output quantities. The variable \c ix_log_wgt should
        be the index of the log_weight among all of the output
        variables, from 0 to <tt>n_out-1</tt>. The list, \c list,
        should include the column names of the parameters and then the
        output quantities (including the log weight column), in order.
    */
    void set(size_t np, size_t n_out, size_t ix_log_wgt,
             table<> &t, std::vector<std::string> list) {

      if (ix_log_wgt>list.size()-np) {
        O2SCL_ERR2("Index of log_wgt is too large in ",
                   "emulator_interpm_idw_table::set().",o2scl::exc_einval);
      }

      iko.full_min=true;
      
      mfrn.resize(n_out);
      
      ix=ix_log_wgt;

      std::vector<std::vector<std::vector<double>>> param_lists;

      std::vector<std::string> col_list_x;
      std::vector<std::string> col_list_y;

      for(size_t j=0;j<list.size();j++) {
        if (j<np) col_list_x.push_back(list[j]);
        else col_list_y.push_back(list[j]);
      }
        
      matrix_view_table<> mvt_x;
      matrix_view_table_transpose<> mvt_y;
      
      mvt_x.set(t,col_list_x);
      mvt_y.set(t,col_list_y);
      
      std::vector<std::vector<double>> ptemp;
      
      for(size_t j=0;j<np;j++) {

        mat_x_col_t xj(mvt_x,j);
        std::vector<double> xj_sort(t.get_nlines());
        vector_copy(t.get_nlines(),xj,xj_sort);
        vector_sort_double(t.get_nlines(),xj_sort);
        
        std::vector<double> diffs;
        o2scl::vector_diffs<std::vector<double>,std::vector<double>>
          (xj_sort,diffs);
        double min=vector_min_value<std::vector<double>,double>
          (diffs.size(),diffs);
        double max=vector_max_value<std::vector<double>,double>
          (diffs.size(),diffs);
        
        std::vector<double> len_list={min/10.0,max*10.0};
        ptemp.push_back(len_list);
      }

      for(size_t iout=0;iout<n_out;iout++) {
        
        mat_y_row_t yiout(mvt_y,iout);
        double min=vector_min_value<mat_y_row_t,double>
          (t.get_nlines(),yiout);
        double max=vector_max_value<mat_y_row_t,double>
          (t.get_nlines(),yiout);

        std::vector<double> l10_list;
        if (max*min<0.0 || min==0.0) {
          l10_list={-15,-13,-11,-9};
        } else {
          min=fabs(min);
          max=fabs(max);
          l10_list.push_back(log10(min/(max/min*1.0e15)));
          l10_list.push_back(log10(min/(max/min*1.0e11)));
          l10_list.push_back(log10(min/(max/min*1.0e7)));
          l10_list.push_back(log10(min/(max/min*1.0e3)));
          vector_out(std::cout,l10_list,true);
        }
        
        mfrn[iout].len.resize(np);
        
        std::vector<std::vector<double>> ptemp2=ptemp;
        
        ptemp2.push_back(l10_list);
        param_lists.push_back(ptemp2);
        
      }

      iko.set_covar(mfrn,param_lists);
      iko.mode=iko.mode_loo_cv_bf;
      iko.set_data(np,n_out,t.get_nlines(),mvt_x,mvt_y,true);

      return;
    }
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                         double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
      
      iko.eval(p,dat);
      iko.sigma(p,dat_unc);
      log_wgt=dat[ix];
      log_wgt_unc=abs(dat_unc[ix]);
      return 0;
    }
    
  };

#ifdef O2SCL_SET_PYTHON

  /** \brief A semi-generic interface for a python emulator
   */
  template<class vec2_t, class vec_t> class emulator_python :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:

    /// The name of the module
    PyObject *p_modname;

    /// The module
    PyObject *p_module;

    /// The class
    PyObject *p_class;

    /// An instance of the class
    PyObject *p_instance;

    /// The code to compute a point using the emulator
    PyObject *p_point_func;

    /// Python object for the number of parameters
    PyObject *p_np;

    /// Python object for the index of log weight
    PyObject *p_ix;

    /// The number of parameters
    size_t num_param;

    /// The number of output data points (not including log_wgt)
    size_t num_out;

    /// If true, then the emulator provides uncertainties (default true)
    bool has_unc;

    /// Index of the "log weight" in the MCMC data vector
    size_t ix;
    
  public:
    
    /** \brief Create an emulator
     */
    emulator_python() {
      p_modname=0;
      p_module=0;
      p_class=0;
      p_instance=0;
      p_point_func=0;
      p_np=0;
      num_param=0;
      num_out=0;
      has_unc=true;
      verbose=0;
    }
    
    /** \brief Clear the memory for all of the python objects
     */
    void decref() {
      if (p_modname!=0) {
        if (verbose>0) {
          std::cout << "Executing decrefs." << std::endl;
          std::cout << "p_modname " << p_modname << std::endl;
        }
        Py_DECREF(p_modname);
        if (verbose>0) std::cout << "p_module." << p_module << std::endl;
        Py_DECREF(p_module);
        if (verbose>0) std::cout << "p_class " << p_class << std::endl;
        Py_DECREF(p_class);
        if (verbose>0) std::cout << "p_instance " << p_instance << std::endl;
        Py_DECREF(p_instance);
        if (verbose>0) std::cout << "p_point_func "
                                 << p_point_func << std::endl;
        Py_DECREF(p_point_func);
        if (verbose>0) std::cout << "p_np " << p_np << std::endl;
        Py_DECREF(p_np);
        if (verbose>0) std::cout << "p_ix " << p_ix << std::endl;
        Py_DECREF(p_ix);
        p_modname=0;
        p_module=0;
        p_class=0;
        p_instance=0;
        p_point_func=0;
        p_np=0;
        p_ix=0;
        if (verbose>0) {
          std::cout << "Finished decrefs." << std::endl;
        }
      }
      return;
    }

    /// Verbosity parameter
    int verbose;
    
    virtual ~emulator_python() {
      decref();
    }
    
    /** \brief Set the emulator

        Set the emulator presuming table containing \c np parameters
        and several output quantities. The string \c log_wgt should be
        the column of the log_weight in the table. The list, \c list,
        should include the column names of the parameters and then the
        output quantities (including the log weight column), in order.
     */
    void set(std::string module, std::string class_name,
             std::string train_func, std::string point_func,
             size_t np, std::string file,
             size_t ix_log_wgt, std::vector<std::string> list,
             bool has_uncerts=true) {
      
      decref();

      ix=ix_log_wgt;
      
      int ret;
      
      num_param=np;
      num_out=list.size()-np;
      has_unc=has_uncerts;
      if (verbose>0) {
        std::cout << "emulator_python::set():" << std::endl;
        std::cout << "  " << np << " parameters." << std::endl;
        std::cout << "  " << num_out << " output quantities." << std::endl;
        if (has_unc) {
          std::cout << "  Uncertainties provided." << std::endl;
        } else {
          std::cout << "  Uncertainties not provided." << std::endl;
        }
      }

      if (verbose>0) {
        std::cout << "Getting module name in unicode." << std::endl;
      }
      p_modname=PyUnicode_FromString(module.c_str());
      if (p_modname==0) {
        O2SCL_ERR2("Module name string creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Obtaining python module." << std::endl;
      }
      p_module=PyImport_Import(p_modname);
      if (p_module==0) {
        O2SCL_ERR2("Module import failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Obtaining python class." << std::endl;
      }
      p_class=PyObject_GetAttrString(p_module,class_name.c_str());
      if (p_class==0) {
        O2SCL_ERR2("Get class failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      // Create an instance of the class
      if (verbose>0) {
        std::cout << "Loading python class." << std::endl;
      }
      if (PyCallable_Check(p_class)==false) {
        O2SCL_ERR2("Check class callable failed in ",
                   "funct_python_method::set_function().",
                   o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Loading python class instance." << std::endl;
      }
      p_instance=PyObject_CallObject(p_class,0);
      if (p_instance==0) {
        O2SCL_ERR2("Instantiate class failed in ",
                   "funct_python_method::set_function().",
                   o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Loading python class training function "
                  << train_func << std::endl;
      }
      PyObject *p_train_func=
        PyObject_GetAttrString(p_instance,train_func.c_str());
      if (p_train_func==0) {
        O2SCL_ERR2("Get training function failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Loading python class point function "
                  << point_func << std::endl;
      }
      p_point_func=PyObject_GetAttrString(p_instance,point_func.c_str());
      if (p_point_func==0) {
        O2SCL_ERR2("Get point function failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Loading python number of parameters." << std::endl;
      }
      p_np=PyLong_FromSsize_t(np);
      if (p_np==0) {
        O2SCL_ERR2("Parameter number object failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Loading python filename " << file << std::endl;
      }
      PyObject *p_file=PyUnicode_FromString(file.c_str());
      if (p_file==0) {
        O2SCL_ERR2("File name string creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Loading python index of log weight." << std::endl;
      }
      p_ix=PyLong_FromSsize_t(ix);
      if (p_ix==0) {
        O2SCL_ERR2("Parameter number object failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Reformatting list." << std::endl;
      }
      std::string list_reformat=list[0];
      for(size_t k=1;k<list.size();k++) {
        list_reformat+=","+list[k];
      }

      if (verbose>0) {
        std::cout << "Creating python list." << std::endl;
      }
      PyObject *p_list=PyUnicode_FromString(list_reformat.c_str());
      if (p_list==0) {
        O2SCL_ERR2("String creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Creating python tuple." << std::endl;
      }
      PyObject *p_args=PyTuple_New(5);
      
      if (verbose>0) {
        std::cout << "Setting first item." << std::endl;
      }
      ret=PyTuple_SetItem(p_args,0,p_np);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 0 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Setting second item." << std::endl;
      }
      ret=PyTuple_SetItem(p_args,1,p_file);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 1 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Setting third item." << std::endl;
      }
      ret=PyTuple_SetItem(p_args,2,p_ix);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 2 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Setting fourth item." << std::endl;
      }
      ret=PyTuple_SetItem(p_args,3,p_list);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 3 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      if (verbose>0) {
        std::cout << "Setting verbose parameter." << std::endl;
      }
      PyObject *p_verbose=PyLong_FromSsize_t(verbose);
      if (p_verbose==0) {
        O2SCL_ERR2("Verbose object failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      ret=PyTuple_SetItem(p_args,4,p_verbose);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 4 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Calling training function." << std::endl;
      }
      PyObject *p_result=PyObject_CallObject(p_train_func,p_args);
      if (p_result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      if (verbose>0) {
        std::cout << "Decref value and result." << std::endl;
      }
      
      Py_DECREF(p_train_func);

      // AWS 3/7/23: once we put the file and list in the args
      // tuple, then they take over the memory so no need to
      // call Py_DECREF() for them.

      Py_DECREF(p_args);
      Py_DECREF(p_result);
      
      if (verbose>0) {
        std::cout << "Done in emulator_python::set()." << std::endl;
      }
      
      return;
    }
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                         double &log_wgt_unc, vec2_t &dat,
                         vec2_t &dat_unc) {

      if (p_modname==0) {
        O2SCL_ERR2("Emulator was not set in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      
      // Create the list object
      PyObject *p_list=PyList_New(n);
      if (p_list==0) {
        O2SCL_ERR2("List creation failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }
      
      // Create a python object from the vector
      std::vector<PyObject *> p_values(n);
      if (verbose>0) {
        std::cout << "Creating python object from vector." << std::endl;
      }
      
      for(size_t i=0;i<n;i++) {
        p_values[i]=PyFloat_FromDouble(p[i]);
        std::cout << "p[i]: " << p[i] << std::endl;
        if (p_values[i]==0) {
          O2SCL_ERR2("Value creation failed in ",
                     "emulator_python::eval_unc().",o2scl::exc_efailed);
        }
        
        // Set the python function arguments
        int iret=PyList_SetItem(p_list,i,p_values[i]);
        if (iret!=0) {
          O2SCL_ERR2("Item set failed in ",
                     "emulator_python::eval_unc().",o2scl::exc_efailed);
        }
      }

      PyObject *p_args=PyTuple_New(1);
      if (p_args==0) {
        O2SCL_ERR2("Create arg tuple failed in ",
                   "multi_funct_python::set_function().",o2scl::exc_efailed);
      }
      
      int ret=PyTuple_SetItem(p_args,0,p_list);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      // Call the python function
      if (verbose>0) {
        std::cout << "Call python point function " << p_args << std::endl;
      }
      PyObject *p_result=PyObject_CallObject(p_point_func,p_args);
      if (p_result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      if (PyList_Check(p_result)==0) {
        O2SCL_ERR2("Function call did not return a list ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }
      
      for(size_t i=0;i<num_out;i++) {
        PyObject *p_y_val=PyList_GetItem(p_result,i);
        if (p_y_val==0) {
          O2SCL_ERR2("Failed to get y list value in ",
                     "mm_funct_python::operator().",o2scl::exc_efailed);
        }
        if (i==0) {
          log_wgt=PyFloat_AsDouble(p_y_val);
        } else {
          dat[i]=PyFloat_AsDouble(p_y_val);
        }
      }

      if (has_unc) {
        size_t n2=num_out;
        for(size_t i=n2;i<2*n2;i++) {
          PyObject *p_y_val=PyList_GetItem(p_result,i);
          if (p_y_val==0) {
            O2SCL_ERR2("Failed to get y list value in ",
                       "mm_funct_python::operator().",o2scl::exc_efailed);
          }
          if (i==n2) {
            log_wgt=PyFloat_AsDouble(p_y_val);
          } else {
            dat_unc[i-n2-1]=PyFloat_AsDouble(p_y_val);
          }
        }
      }

      if (verbose>0) {
        std::cout << "Decref list." << std::endl;
      }
      Py_DECREF(p_list);
      
      if (verbose>0) {
        std::cout << "Decref result." << std::endl;
      }
      Py_DECREF(p_result);
  
      return 0;
    }
    
  };

#endif
  
  /** \brief Placeholder for an adaptive emulator
   */
  template<class emu_t, class exact_t,
    class vec2_t, class vec_t> class emulator_adapt :
    public emulator_base<vec2_t,vec_t> {

  public:

    /** \brief Desc
     */
    emu_t emu_base;
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat
    */
    virtual int eval(size_t n, const vec_t &p, double &log_wgt,
             vec2_t &dat) {
      
      return 0;
    }
    
    
  };
  
  // End of namespace
}

#endif
