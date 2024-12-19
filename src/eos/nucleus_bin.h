/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2021-2025, Andrew W. Steiner
  
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
#include <o2scl/nucmass_ldrop_ext.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/fermion.h>
#include <o2scl/hdf_nucmass_io.h>

/** \brief Desc
 */
class nucleus_bin {
  
public:
  
  /// \name Atomic mass evaluations
  //@{
  o2scl::nucmass_ame ame16;
  o2scl::nucmass_ame ame20exp;
  o2scl::nucmass_ame ame20round;
  o2scl::nucmass_ame ame95rmd;
  o2scl::nucmass_ame ame03round;
  o2scl::nucmass_ame ame03;
  o2scl::nucmass_ame ame95exp;
  o2scl::nucmass_ame ame12;
  //@}

  /// \name Nuclear mass tables from FRIB Mass Explorer
  //@{
  o2scl::nucmass_gen ddme2;
  o2scl::nucmass_gen ddmed;
  o2scl::nucmass_gen ddpc1;
  o2scl::nucmass_gen nl3s;
  o2scl::nucmass_gen sly4;
  o2scl::nucmass_gen skms;
  o2scl::nucmass_gen skp;
  o2scl::nucmass_gen sv_min;
  o2scl::nucmass_gen unedf0;
  o2scl::nucmass_gen unedf1;
  //@}

  /// \name Other theoretical nuclear mass tables
  //@{
  o2scl::nucmass_mnmsk m95;
  o2scl::nucmass_mnmsk m16;
  o2scl::nucmass_ktuy kt;
  o2scl::nucmass_ktuy kt2;
  o2scl::nucmass_wlw wlw1;
  o2scl::nucmass_wlw wlw2;
  o2scl::nucmass_wlw wlw3;
  o2scl::nucmass_wlw wlw4;
  o2scl::nucmass_wlw wlw5;
  o2scl::nucmass_sdnp sdnp1;
  o2scl::nucmass_sdnp sdnp2;
  o2scl::nucmass_sdnp sdnp3;
  o2scl::nucmass_dz_table dz;
  //@}

  /// \name HFB mass tables from the Brussels group
  //@{
  o2scl::nucmass_hfb hfb2;
  o2scl::nucmass_hfb hfb8;
  o2scl::nucmass_hfb hfb14;
  o2scl::nucmass_hfb hfb14_v0;
  o2scl::nucmass_hfb_sp hfb17;
  o2scl::nucmass_hfb_sp hfb21;
  o2scl::nucmass_hfb_sp hfb22;
  o2scl::nucmass_hfb_sp hfb23;
  o2scl::nucmass_hfb_sp hfb24;
  o2scl::nucmass_hfb_sp hfb25;
  o2scl::nucmass_hfb_sp hfb26;
  o2scl::nucmass_hfb_sp hfb27;
  //@}

  /// \name Nuclear mass fits
  //@{
  o2scl::nucmass_semi_empirical se;
  o2scl::nucmass_frdm frdm;
  o2scl::nucmass_dz_fit dzf;
  o2scl::nucmass_dz_fit_33 dzf33;
  o2scl::nucmass_frdm_shell frdm_shell;
  o2scl::nucmass_ldrop_shell ldrop_shell;
  o2scl::nucmass_ldrop_ext ldrop_ext;
  //@}

  /// \name Other objects
  //@{
  o2scl::eos_had_skyrme sk;
  o2scl::fermion nrn;
  o2scl::fermion nrp;
  o2scl::nucmass_fit fitter;
  
  /// Mass of the neutron
  double m_neut;

  /// Mass of the proton
  double m_prot;

  /// Mass of the electron
  double m_elec;

  /// Atomic mass unit
  double m_amu;

  /** \brief Verbosity parameter (default 1)
   */
  int verbose;
  //@}
  
protected:
  
  /// \name Other nuclear mass objects
  //@}
  /** \brief Information object to look up element name given Z
   */
  o2scl::nucmass_info nmi;

  /// Number of tables (theory and experiment; set in constructor)
  size_t n_tables;
  
  /// Number of fits (set in constructor)
  size_t n_fits;
  
  /// List of pointers to tables
  std::vector<o2scl::nucmass_table *> nmd;

  /// List of pointers to fits
  std::vector<o2scl::nucmass_fit_base *> nmfd;

  /// List of table names
  std::vector<std::string> table_names;

  /// List of fit names
  std::vector<std::string> fit_names;
  //@}
  
  /// \name Parameter objects
  //@{
  /// Verbosity parameter
  o2scl::cli::parameter_int p_verbose;
  //@}

  // Compute the drip lines
  //int drip_lines(std::vector<int> &v_Z, std::vector<int> &v_Npd,
  //std::vector<int> &v_Nnd, int file_index=-1);
  
public:
  
  nucleus_bin();
  
  virtual ~nucleus_bin() {
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

  /** \brief Get a nucleus by Z and N
   */
  int get(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Desc
   */
  int isotope(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Output all table names and number of nuclei
   */
  int tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Improve the fits using an optimizer
   */
  int fits(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compare the quality of the fits and tables to experiment
   */
  int compare(std::vector<std::string> &sv, bool itive_com);

  /** \brief Output table and fit references
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

  /** \brief Setup the command-line interface
   */
  virtual void setup_cli(o2scl::cli &cl);
  
};
  
