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

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  /** \brief Emulator base class
   */
  template<class data_t, class vec_t=ubvector> class emulator_base {
    
  protected:
    
  public:

  /** \brief Desc
   */
  int eval(size_t n, const vec_t &p, double &log_wgt,
           data_t &dat);
  
  };
  
  /** \brief Emulator with uncertainty base class
   */
  template<class data_t, class data_unc_t,
    class vec_t=ubvector> class emulator_unc : public
    emulator_base<data_t,vec_t> {
    
  public:

  /** \brief Desc
   */
  int eval_unc(size_t n, const vec_t &p, double &log_wgt,
               double &log_wgt_unc, data_t &dat, data_unc_t &dat_unc);
  
  /** \brief Desc
   */
  int eval(size_t n, const vec_t &p, double &log_wgt,
           data_t &dat) {
    double log_wgt_unc,
    data_unc_t dat_unc;
    return eval_unc(n,p,log_wgt,log_wgt_unc,dat,dat_unc);
  }
    
  };
  
  /** \brief Emulate data stored in a table object with interpm_idw
   */
  template<class vec2_t, class vec_t> class emulator_interpm_idw_table :
    public emulator_unc<vec2_t,vec2_t,vec_t> {

  protected:
    
    /// Desc
    o2scl::interpm_idw<> ii;

    /// Desc
    const_matrix_view_table_transpose<> cmvtt;
    
  public:
    
    /** \brief Desc
     */
    emulator_interpm_idw_table() {
    }
    
    /** \brief Desc
     */
    void set(size_t np, size_t n_out, size_t ix_log_wgt,
             table<> &t, std::vector<std::string> list) {
      cmvtt.set(t,list);
      set_data(np,list.size()-np,t.get_nlines(),cmvtt);
    }
    
    /** \brief Desc
     */
    int eval_unc(size_t n, const vec_t &p, double &log_wgt,
                 double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
      
      ii.eval_err<vec2_t>(p,dat,dat_unc);
      log_wgt=dat[ix_log_wgt];
      log_wgt_unc=dat_unc[ix_log_wgt];
      return 0;
    }
    
  };

  /** \brief Desc
   */
  template<class emu_t, class exact_t,
    class vec2_t, class vec_t> class emulator_adapt :
    public emulator_base<vec2_t,vec_t> {

  public:

    /** \brief Desc
     */
    emu_t emu_base;
    
    /** \brief Desc
     */
    int eval(size_t n, const vec_t &p, double &log_wgt,
             vec2_t &dat) {
      
      return 0;
    }
    
    
  };
  
  // End of namespace
}

#endif
