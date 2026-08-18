/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_MULTIP_H
#define O2SCL_FERMION_MULTIP_H

/** \file fermion.h
    \brief File defining \ref o2scl::fermion_tl
*/
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include <o2scl/set_multip.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>
#include <o2scl/fermion_nonrel.h>

namespace o2scl {

#ifdef O2SCL_SET_MULTIP
  
  /** \brief 25-digit precision version of \ref
      o2scl::part_deriv_press_tl
   */
  typedef part_deriv_press_tl<o2fp_25> part_deriv_press_fp25;
  
#endif
  
#ifdef O2SCL_SET_MULTIP
  /** \brief 25-digit floating point version of \ref fermion_deriv_tl
   */
  typedef fermion_deriv_tl<o2fp_25> fermion_deriv_fp25;
#endif

#ifdef O2SCL_SET_MULTIP
  /** \brief 25-digit floating point version of \ref boson_deriv_tl
   */
  typedef boson_deriv_tl<o2fp_25> boson_deriv_fp25;
#endif

#ifdef O2SCL_SET_MULTIP
  /** \brief 25-digit floating point version of \ref part_deriv_tl
   */
  typedef part_deriv_tl<o2fp_25> part_deriv_fp25;
#endif

#ifdef O2SCL_SET_MULTIP
  
  /** \brief 25-digit precision thermodynamics object
   */
  typedef thermo_tl<o2fp_25> thermo_fp25;
  
#endif
  
#ifdef O2SCL_SET_MULTIP
  
  /// Boson type for 25-digit floating points
  typedef boson_tl<o2fp_25> boson_fp25;
  
#endif
  
#ifdef O2SCL_SET_MULTIP
  
  /** \brief Long double version of 
      \ref o2scl::boson_rel_tl 
  */
  typedef boson_rel_tl
  <bessel_K_exp_integ_boost<long double,
                            o2fp_25>,long double> boson_rel_ld;
  
  /** \brief 25-digit version of 
      \ref o2scl::boson_rel_tl 
  */
  typedef boson_rel_tl
  <bessel_K_exp_integ_boost<o2fp_25,
                            o2fp_35>,o2fp_25> boson_rel_fp25;
   
  
#endif
  
#ifdef O2SCL_SET_MULTIP
  
  /** \brief Long double version of 
      \ref o2scl::fermion_nonrel_tl 
  */
  typedef fermion_nonrel_tl
  <fermion_tl<long double>,
   fermi_dirac_integ_direct<long double,funct_fp25,
                            o2fp_25>,
   bessel_K_exp_integ_boost<long double,
                            o2fp_25>,
   root_brent_gsl<funct_ld,long double>,
   funct_ld,long double> fermion_nonrel_ld;
  
#endif

#ifdef O2SCL_SET_MPFR
#ifdef O2SCL_SET_MULTIP
  
  typedef fermion_tl<o2fp_25> fermion_fp25;
  typedef fermion_tl<o2fp_35> fermion_fp35;
  typedef fermion_tl<o2fp_50> fermion_fp50;
  typedef fermion_tl<o2fp_100> fermion_fp100;
  
#endif
#endif
  
#if defined (O2SCL_SET_MULTIP) || defined (DOXYGEN)
  
  /** \brief Long double version of 
      \ref o2scl::fermion_rel_tl 
  */
  class fermion_rel_ld : public
  fermion_rel_tl<
    fermion_tl<long double>,
    fermi_dirac_integ_direct<long double,funct_fp25,
                             o2fp_25>,
    bessel_K_exp_integ_boost<long double,
                             o2fp_25>,
    inte_double_exp_boost<>,
    inte_double_exp_boost<>,
    root_cern<funct_ld,long double>,
    funct_ld,long double> {
    
  public:
    
    fermion_rel_ld() {
      //density_root.test_form=2;

      // AWS, 2/19/25: I haven't yet optimized the value of
      // upper_limit_fac for this type
      upper_limit_fac=40;
    }
    
  };
  
  /** \brief 25-digit version of 
      \ref o2scl::fermion_rel_tl 
  */
  class fermion_rel_fp25 : public
  fermion_rel_tl<fermion_tl<o2fp_25>,
                 fermi_dirac_integ_direct<
                   o2fp_25,funct_fp35,
                   o2fp_35>,
                 bessel_K_exp_integ_boost<o2fp_25,
                                          o2fp_35>,
                 inte_double_exp_boost<>,
                 inte_double_exp_boost<>,
                 root_cern<funct_fp25,o2fp_25>,
                 funct_fp25,
                 o2fp_25> {

  public:
    
    fermion_rel_fp25() {
      //density_root.test_form=2;

      // AWS, 2/19/25: I haven't yet optimized the value of
      // upper_limit_fac for this type
      upper_limit_fac=60;
    }
    
  };
  
#endif  

#ifdef O2SCL_SET_MULTIP
  
  /** \brief Long double version of 
      \ref o2scl::fermion_deriv_rel_tl 
  */
  typedef fermion_deriv_rel_tl<fermion_deriv_tl<long double>,
                               fermion_rel_ld,
			       inte_double_exp_boost<>,
			       inte_double_exp_boost<>,
			       long double>
  fermion_deriv_rel_ld;

  /** \brief 25-digit version of 
      \ref o2scl::fermion_deriv_rel_tl 
  */
  typedef fermion_deriv_rel_tl<fermion_deriv_tl<o2fp_25>,
                               fermion_rel_fp25,
			       inte_double_exp_boost<>,
			       inte_double_exp_boost<>,
			       o2fp_25>
  fermion_deriv_rel_fp25;
  
#endif  
  
}

#endif
