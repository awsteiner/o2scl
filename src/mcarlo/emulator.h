/*
  -------------------------------------------------------------------
  
  Copyright (C) 2022, Andrew W. Steiner
  
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
/** \file emulator.h
 */
#ifndef O2SCL_EMULATOR_H
#define O2SCL_EMULATOR_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

#include <o2scl/interpm_idw.h>
#include <o2scl/interpm_krige.h>

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  /** \brief Emulator base class
   */
  template<class data_t, class vec_t=ubvector> class emulator_base {
    
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
           class vec_t=ubvector> class emulator_unc : public
  emulator_base<data_t,vec_t> {
    
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
      // dat_unc object using a copy constructor
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
  template<class vec2_t, class vec_t> class emulator_interpm_idw_table :
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
  template<class vec2_t, class vec_t> class emulator_interpm_krige_table :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef o2scl::matrix_view_table<> mat_x_t;
    typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
    typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
    typedef o2scl::matrix_view_table_transpose<> mat_y_t;
    typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
    
    /// The view of the user-specified table
    matrix_view_table_transpose<> cmvtt;

    size_t ix;

    std::vector<std::string> col_list_x;
    std::vector<std::string> col_list_y;
    matrix_view_table<> mvt_x;
    matrix_view_table_transpose<> mvt_y;
    
  public:
    
    /// The internal interpolation object
    interpm_krige_optim
    <ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
     mat_y_t,mat_y_row_t,ubmatrix,
     o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > iko;

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

      ix=ix_log_wgt;

      col_list_x.clear();
      col_list_y.clear();
      for(size_t j=0;j<list.size();j++) {
        if (j<np) col_list_x.push_back(list[j]);
        else col_list_y.push_back(list[j]);
      }
      mvt_x.set(t,col_list_x);
      mvt_y.set(t,col_list_y);
      
      ubvector len_precompute;
      iko.set_data<>(2,1,8,mvt_x,mvt_y,len_precompute);

      return;
    }
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                 double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
      
      //iko.eval_err<vec_t,vec2_t,vec2_t>(p,dat,dat_unc);
      log_wgt=dat[ix];
      log_wgt_unc=dat_unc[ix];
      return 0;
    }
    
  };

#ifdef O2SCL_PYTHON
  
  template<class vec2_t, class vec_t> class emulator_python :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:

    PyObject *p_modname;
    PyObject *p_module;
    PyObject *p_point_func;
    PyObject *p_np;
    
  public:
    
    /** \brief Create an emulator
     */
    emulator_python() {
      p_modname=0;
      p_module=0;
      p_point_func=0;
      p_np=0;
    }
    
    /** \brief Clear the memory for all of the python objects
     */
    void decref() {
      if (p_modname!=0) {
        Py_DECREF(p_modname);
        Py_DECREF(p_module);
        Py_DECREF(p_point_func);
        Py_DECREF(p_np);
        p_modname=0;
        p_module=0;
        p_point_func=0;
        p_np=0;
      }
      return;
    }

    virtual ~emulator_python() {
      decref();
    }
    
    /** \brief Set the emulator

        Set the emulator using a table containing \c np parameters and
        \c n_out output quantities. The variable \c ix_log_wgt should
        be the index of the log_weight among all of the output
        variables, from 0 to <tt>n_out-1</tt>. The list, \c list,
        should include the column names of the parameters and then the
        output quantities (including the log weight column), in order.
     */
    void set(std::string module, std::string train_func,
             std::string point_func, size_t np, std::string file,
             std::string log_wgt, std::vector<std::string> list) {

      decref();
      
      int ret;

      p_modname=PyUnicode_FromString(module.c_str());
      if (p_modname==0) {
        O2SCL_ERR2("Module name string creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      p_module=PyImport_Import(p_modname);
      if (p_module==0) {
        O2SCL_ERR2("Module import failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      PyObject *p_train_func=
        PyObject_GetAttrString(p_module,train_func.c_str());
      if (p_train_func==0) {
        O2SCL_ERR2("Get training function failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      p_point_func=PyObject_GetAttrString(p_module,point_func.c_str());
      if (p_point_func==0) {
        O2SCL_ERR2("Get point function failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      p_np=PyLong_FromSsize_t(np);
      if (p_np==0) {
        O2SCL_ERR2("Parameter number object failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      PyObject *p_file=PyUnicode_FromString(file.c_str());
      if (p_file==0) {
        O2SCL_ERR2("File name string creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      std::string list_reformat=log_wgt;
      for(size_t k=0;k<list_list.size();k++) {
        list_reformat+=","+list_list[k];
      }
      PyObject *p_list=PyUnicode_FromString(list_reformat.c_str());
      if (p_list==0) {
        O2SCL_ERR2("String creation failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      PyObject *p_args=PyTuple_New(3);
      ret=PyTuple_SetItem(p_args,0,p_np);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 0 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      ret=PyTuple_SetItem(p_args,1,p_file);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 1 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      ret=PyTuple_SetItem(p_args,2,p_out);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set 2 failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }

      PyObject *p_result=PyObject_CallObject(p_train_func,p_args);
      if (p_result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "emulator_python::set().",o2scl::exc_efailed);
      }
      
      //if (verbose>0) {
      //cout << "Decref value and result." << endl;
      //}
      Py_DECREF(p_train_func);
      Py_DECREF(p_file);
      Py_DECREF(p_list);
      Py_DECREF(p_args);
      Py_DECREF(p_result);
      
      return;
    }
    
    /** \brief Evaluate the emulator at the point \c p returning
        \c log_wgt and \c dat and their uncertainties
     */
    virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                         double &log_wgt_unc, vec2_t &dat,
                         vec2_t &dat_unc) {

      // Create the list object
      PyObject *p_list=PyList_New(n);
      if (p_list==0) {
        O2SCL_ERR2("List creation failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      double y;
      
      // Create a python object from the vector
      std::vector<PyObject *> p_values(n);
      if (verbose>0) {
        std::cout << "Creating python object from vector." << std::endl;
      }
      
      for(size_t i=0;i<n;i++) {
        p_values[i]=PyFloat_FromDouble(p[i]);
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
      
      int ret=PyTuple_SetItem(pArgs,0,p_list);
      if (ret!=0) {
        O2SCL_ERR2("Tuple set failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      // Call the python function
      if (verbose>0) {
        std::cout << "Call python function." << std::endl;
      }
      PyObject *p_result=PyObject_CallObject(p_point_func,p_args);
      if (result==0) {
        O2SCL_ERR2("Function call failed in ",
                   "emulator_python::eval_unc().",o2scl::exc_efailed);
      }

      size_t n=10;
      for(size_t i=0;i<n;i++) {
        PyObject *p_y_val=PyList_GetItem(p_result,i);
        if (yval==0) {
          O2SCL_ERR2("Failed to get y list value in ",
                     "mm_funct_python::operator().",o2scl::exc_efailed);
        }
        if (i==0) {
          log_wgt=PyFloat_AsDouble(p_y_val);
        } else {
          dat[i-1]=PyFloat_AsDouble(p_y_val);
        }
        if (verbose>0) {
          std::cout << "Decref yval " << i << " of " << pValues.size()
                    << std::endl;
        }
        Py_DECREF(yval);
      }
      
      //log_wgt=dat[ix];
      //log_wgt_unc=dat_unc[ix];
      
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
