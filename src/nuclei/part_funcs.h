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
#ifndef O2SCL_PART_FUNCS_H
#define O2SCL_PART_FUNCS_H

/** \file part_funcs.h
    \brief File defining \ref o2scl::nucmass_dz_table and other classes
*/

#include <o2scl/table_units.h>
#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/inte_qag_gsl.h>

namespace o2scl {

  /** \brief Desc
   */
  class part_funcs {
    
  public:
    
    part_funcs();

    /** \brief Desc
     */
    int few78(int Z, int N, double T, double &pf, double &dpfdT);
    
    /** \brief Desc
     */
    int rtk97(int Z, int N, double T, double &pf, double &dpfdT);
    
    /** \brief Desc
     */
    int rt00(int Z, int N, double T, double &pf, double &dpfdT);
    
    /** \brief Desc
     */
    int r03(int Z, int N, double T, double &pf, double &dpfdT);

    /** \brief Desc
     */
    int load_rt00(std::string fname="", bool external=false);
    
    /** \brief Desc
     */
    int load_r03(std::string fname="", bool external=false);

    /** \brief Desc
     */
    void compare_spin_deg();

    /** \brief Desc
     */
    double get_spin_deg(int Z, int N);
    
    /** \brief Load all of the data necessary
     */
    void load(std::string dir="") {
      if (dir.length()>0) {
        load_rt00(dir+"/pf_frdm_low.o2");
        load_r03(dir+"/pf_frdm_high.o2");
        o2scl_hdf::ame_load_ext(ame,dir+"/ame12.o2","ame12.o2");
        o2scl_hdf::mnmsk_load(mnmsk,"msis16",
                              dir+"/msis16.o2");
        o2scl_hdf::hfb_sp_load(hfb,27,dir);
      } else {
        load_rt00();
        load_r03();
        o2scl_hdf::ame_load(ame);
        o2scl_hdf::mnmsk_load(mnmsk);
        o2scl_hdf::hfb_sp_load(hfb,27);
      }
      std::cout << "rt00: " << tab_rt00.get_nlines() << std::endl;
      std::cout << "r03: " << tab_r03.get_nlines() << std::endl;
      std::cout << "ame: " << ame.is_included(28,28) << std::endl;
      std::cout << "mnmsk: " << mnmsk.is_included(28,28) << std::endl;
      std::cout << "hfb: " << hfb.is_included(28,28) << std::endl;
      return;
    }
    
  protected:

    /// For unit conversions (set in constructor)
    convert_units<double> &cu;
    
    /// Desc
    table_units<> tab_rt00;

    /// Desc
    table_units<> tab_r03;

    /// Desc
    nucmass_mnmsk mnmsk;

    /// Desc
    nucmass_ame ame;

    /// Desc
    nucmass_hfb_sp hfb;
    
    /// Integrator
    o2scl::inte_qag_gsl<> iqg;

    double rtk_alpha;
    double rtk_beta;
    double rtk_gamma;

    /** \brief Desc
     */
    int shen10(int Z, int N, double T, double &pf, double &dpfdT,
               int a_delta=0);

    /** \brief Integrand for partition function when \f$ \delta \f$
        is smaller than \f$ E_d \f$

        From eq. (25) & (27) in Shen 10.
    */
    double delta_small_iand(double E, double T_MeV, double delta,
                            double a);
    
    /** \brief Integrand for derivative of partition function
        with respect to temperature when \f$ \delta \f$
        is smaller than \f$ E_d \f$

        From eq. (26) & (27) in Shen 10.
    */
    double delta_small_iand_prime(double E, double T_MeV, double delta,
                                  double a);
    
    /** \brief Integrand when \f$ \delta \f$ is greater than \f$ E_d \f$

        From eq. (25) and (30) in Shen 10.
    */
    double delta_large_iand(double E, double T_MeV, double delta,
                            double Tc, double C);
    
    /** \brief Integrand for temperature derivative when \f$ \delta
        \f$ is greater than \f$ E_d \f$

        From eq. (26) and (30) in Shen 10.
    */
    double delta_large_iand_prime(double E, double T_MeV, double delta,
                                  double Tc, double C);
    
  };

}

#endif
