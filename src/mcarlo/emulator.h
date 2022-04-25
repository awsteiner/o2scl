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

#include <o2scl/interpm_idw.h>

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
      // dat_unc
      data_unc_t dat_unc=dat;
      return eval_unc(n,p,log_wgt,log_wgt_unc,dat,dat_unc);
    }
    
  };
  
  /** \brief Emulate data stored in a table object with interpm_idw
   */
  template<class vec2_t, class vec_t> class emulator_interpm_idw_table :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:
    
    /// The view of the user-specified table
    const_matrix_view_table_transpose<> cmvtt;

    size_t ix;
    
  public:
    
    /// The internal interpolation objec
    o2scl::interpm_idw<o2scl::const_matrix_view_table_transpose<>> ii;
    
    /** \brief Create an emulator
     */
    emulator_interpm_idw_table() {
    }
    
    /** \brief Set the emulator

        Set the emulator using a table containing \c np parameters and
        \c n_out output quantities using \c ix_log_wgt as the column
        index of the column which contains log_wgt and \c list as
        the the column names of the parameters and output quantities,
        in order.
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

      //ii.verbose=3;
      //std::cout << "E: " << dat.size() << " " << dat_unc.size()
      //<< std::endl;
      ii.eval_err<vec_t,vec2_t,vec2_t>(p,dat,dat_unc);
      log_wgt=dat[ix];
      log_wgt_unc=dat_unc[ix];
      return 0;
    }
    
  };

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
