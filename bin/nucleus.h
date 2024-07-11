/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2021-2024, Andrew W. Steiner
  
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
#include <iostream>
#include <string>
#include <o2scl/cli.h>
//#include <vector>

#include <o2scl/nucmass.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_ktuy.h>
#include <o2scl/nucmass_dglg.h>
#include <o2scl/nucmass_wlw.h>
#include <o2scl/nucmass_sdnp.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/nucmass_gen.h>
#include <o2scl/nucmass_ldrop_shell.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/fermion.h>
#include <o2scl/hdf_nucmass_io.h>

/** \brief Desc
 */
class nucleus_class {
  
 protected:

  /** \brief Desc
   */
  o2scl::nucmass_info nmi;

  /// \name Experimental nuclear mass tables
  //@{
  /** \brief Desc
   */
  o2scl::nucmass_ame ame16;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame20exp;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame20round;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame95rmd;

  /** \brief Desc
   */
  o2scl::nucmass_gen ddme2;

  /** \brief Desc
   */
  o2scl::nucmass_gen ddmed;

  /** \brief Desc
   */
  o2scl::nucmass_gen ddpc1;

  /** \brief Desc
   */
  o2scl::nucmass_gen nl3s;

  /** \brief Desc
   */
  o2scl::nucmass_gen sly4;

  /** \brief Desc
   */
  o2scl::nucmass_gen skms;

  /** \brief Desc
   */
  o2scl::nucmass_gen skp;

  /** \brief Desc
   */
  o2scl::nucmass_gen sv_min;

  /** \brief Desc
   */
  o2scl::nucmass_gen unedf0;

  /** \brief Desc
   */
  o2scl::nucmass_gen unedf1;
  
  /** \brief Desc
   */
  o2scl::nucmass_ame ame03round;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame03;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame95exp;

  /** \brief Desc
   */
  o2scl::nucmass_ame ame12;
  //@}

  /// \name Theoretical nuclear mass tables
  //@{
  /// Number of tables (set in constructor)
  size_t n_tables;
  
  /** \brief Desc
   */
  o2scl::nucmass_mnmsk m95;

  /** \brief Desc
   */
  o2scl::nucmass_mnmsk m16;

  /** \brief Desc
   */
  o2scl::nucmass_ktuy kt;

  /** \brief Desc
   */
  o2scl::nucmass_ktuy kt2;

  /** \brief Desc
   */
  o2scl::nucmass_wlw wlw1;

  /** \brief Desc
   */
  o2scl::nucmass_wlw wlw2;

  /** \brief Desc
   */
  o2scl::nucmass_wlw wlw3;

  /** \brief Desc
   */
  o2scl::nucmass_sdnp sdnp1;

  /** \brief Desc
   */
  o2scl::nucmass_sdnp sdnp2;

  /** \brief Desc
   */
  o2scl::nucmass_sdnp sdnp3;

  /** \brief Desc
   */
  o2scl::nucmass_hfb hfb2;

  /** \brief Desc
   */
  o2scl::nucmass_hfb hfb8;

  /** \brief Desc
   */
  o2scl::nucmass_hfb hfb14;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb hfb14_v0;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb17;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb21;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb22;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb23;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb24;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb25;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb26;
  
  /** \brief Desc
   */
  o2scl::nucmass_hfb_sp hfb27;

  /** \brief Desc
   */
  o2scl::nucmass_dz_table dz;

  /** \brief Desc
   */
  o2scl::nucmass_gen ng[10];
  //@}

  /// \name Nuclear mass fits
  //@{
  /// Number of fits (set in constructor)
  size_t n_fits;
  
  /** \brief Desc
   */
  o2scl::nucmass_semi_empirical se;

  /** \brief Desc
   */
  o2scl::nucmass_frdm frdm;
  
  /** \brief Desc
   */
  o2scl::nucmass_dz_fit dzf;
  
  /** \brief Desc
   */
  o2scl::nucmass_dz_fit_33 dzf33;

  /** \brief Desc
   */
  o2scl::nucmass_frdm_shell frdm_shell;
  
  /** \brief Desc
   */
  o2scl::nucmass_ldrop_shell ldrop_shell;
  //@}

  /** \brief Desc
   */
  o2scl::eos_had_skyrme sk;

  /** \brief Desc
   */
  o2scl::fermion nrn, nrp;
  
  /** \brief Desc
   */
  o2scl::nucmass_fit fitter;
  
  /** \brief Verbosity parameter (default 1)
   */
  int verbose;

  /// \name Parameter objects
  //@{
  /// Verbosity parameter
  o2scl::cli::parameter_int p_verbose;
  //@}

  /// Table list
  std::vector<o2scl::nucmass_table *> nmd;

  /// Table list
  std::vector<o2scl::nucmass_fit_base *> nmfd;

  /// List
  std::vector<std::string> table_names;

  /// List of fits
  std::vector<std::string> fit_names;

  /// Mass of the neutron
  double m_neut;

  /// Mass of the proton
  double m_prot;

  /// Mass of the electron
  double m_elec;

  /// Atomic mass unit
  double m_amu;

  /// Desc
  int drip_lines(std::vector<int> &v_Z, std::vector<int> &v_Npd,
                 std::vector<int> &v_Nnd, int file_index=-1);
  
public:
  
  nucleus_class();
  
  virtual ~nucleus_class() {
  }

  /** \brief Desc
   */
  o2scl::table_units<> erler;
  
  /** \brief Desc
   */
  std::vector<o2scl::nucleus> common_dist;

  /** \brief Desc
   */
  std::vector<o2scl::nucleus> exp_dist;

  /** \brief Desc
   */
  std::vector<o2scl::nucleus> moller_dist;

  /** \brief Desc
   */
  int get(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Desc
   */
  int isotope(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Desc
   */
  int tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int fits(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int refs(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int cdist(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int info_matrix(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int mass_fit(std::vector<std::string> &sv, bool itive_com);
    
  /** \brief Compute the neutron and proton drip lines
   */
  int drip_lines_esym(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  virtual void setup_cli(o2scl::cli &cl);
  
};
  
