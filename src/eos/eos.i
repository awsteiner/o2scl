# Interface file for o2scl eos classes
#
namespace o2scl
py_class_doc |
| Python interface for O2scl class ``%name%``.
| See
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _eos:
|
| EOS classes
| ===========
|
| :ref:`O2sclpy <o2sclpy>`
# 
# Include statements for C++ header file
# 
h_include <o2scl/eos_base.h>
h_include <o2scl/eos_leptons.h>
h_include <o2scl/eos_had_base.h>
h_include <o2scl/eos_had_skyrme.h>
h_include <o2scl/eos_had_apr.h>
h_include <o2scl/eos_had_rmf.h>
h_include <o2scl/eos_quark_bag.h>
h_include <o2scl/eos_quark_njl.h>
h_include <o2scl/part.h>
h_include <o2scl/fermion_nonrel.h>
h_include <o2scl/fermion_deriv_nr.h>
h_include <o2scl/hdf_eos_io.h>
h_include <o2scl/nstar_cold.h>
h_include <o2scl/nstar_rot.h>
h_include <o2scl/tov_love.h>
h_include <o2scl/eos_tov.h>
h_include <o2scl/nucleus_rmf.h>
h_include <o2scl/nucmass_ldrop.h>
h_include <o2scl/nucmass_ldrop_shell.h>
h_include <o2scl/nucleus_bin.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/eos_python.h>
cpp_include <o2scl/part.h>
cpp_include <o2scl/fermion_nonrel.h>
cpp_include <o2scl/fermion_deriv_nr.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
cpp_using o2scl_hdf
#
# Additional python headers
#
# We need shared_ptr_table_units from base
#
py_header from o2sclpy.base import *
py_header from o2sclpy.part import *
py_header from o2sclpy.nuclei import *
#
# ------------------------------------------------------
# 
# Class eos_base
#
class eos_base abstract
- o2scl::thermo def_thermo
# 1/17/2021 removing set_thermo and get_thermo for now, as yanic does
# not yet handle these kind of functions.
# - function set_thermo
#  - void
#  - o2scl::thermo &th
# - function get_thermo
#   - const o2scl::thermo &
#
# ------------------------------------------------------
# 
# Class eos_leptons
#
class eos_leptons
- o2scl::thermo th
- o2scl::fermion e
- o2scl::fermion mu  
- o2scl::boson ph
- o2scl::part_deriv_press ed  
- o2scl::part_deriv_press mud  
- o2scl::part_deriv_press phd
- bool include_muons  
- bool include_deriv
- bool pde_from_density
- int verbose
- bool err_nonconv
- function default_acc
  - void
- function improved_acc
  - void
- function ld_acc
  - void
- function fp_25_acc
  - void
- function pair_mu
  - int
  - double T
- function pair_mu_eq
  - int
  - double T
- function pair_density
  - int
  - double T
- function pair_density_eq
  - int
  - double nq
  - double T
- o2scl::fermion_rel frel
#
# ------------------------------------------------------
# 
# Class eos_had_base
#
class eos_had_base abstract
- parent eos_base
- double eoa
- double msom
- double comp
- double n0
- double esym
- double kprime
- bool err_nonconv
- o2scl::fermion def_neutron
- o2scl::fermion def_proton
- function calc_e
  - int
  - o2scl::fermion &n
  - o2scl::fermion &p
  - o2scl::thermo &th
- function calc_p
  - int
  - o2scl::fermion &n
  - o2scl::fermion &p
  - o2scl::thermo &th
- function fcomp
  - double
  - double nb
  - double delta [0.0]
- function fcomp_err
  - double
  - double nb
  - double delta
  - out double &unc
- function feoa
  - double
  - double nb
  - double delta [0.0]
- function fesym
  - double
  - double nb
  - double delta [0.0]
- function fesym_err
  - double
  - double nb
  - double delta
  - out double &unc
- function fesym_slope
  - double
  - double nb
  - double delta [0.0]
- function fesym_curve
  - double
  - double nb
  - double delta [0.0]
- function fesym_skew
  - double
  - double nb
  - double delta [0.0]
- function fesym_diff
  - double
  - double nb
- function feta
  - double
  - double nb
- function feta_prime
  - double
  - double nb
- function fkprime
  - double
  - double nb
  - double delta [0.0]
- function fmsom
  - double
  - double nb
  - double delta [0.0]
- function f_effm_neut
  - double
  - double nb
  - double delta [0.0]
- function f_effm_prot
  - double
  - double nb
  - double delta [0.0]
- function f_effm_scalar
  - double
  - double nb
  - double delta
- function f_effm_vector
  - double
  - double nb
  - double delta
- function fn0
  - int
  - double delta
  - out double &nb                             
  - out double &leoa
- function f_number_suscept
  - void
  - double mun
  - double mup
  - out double &dPdnn
  - out double &dPdnp
  - out double &dPdpp
- function f_inv_number_suscept
  - void
  - double mun
  - double mup
  - out double &dednn
  - out double &dednp
  - out double &dedpp
- function saturation
  - int
- function calc_mun_e
  - double
  - double nn
  - double np
- function calc_mup_e
  - double
  - double nn
  - double np
- function calc_ed
  - double
  - double nn
  - double np
- function calc_pr
  - double
  - double nn
  - double np
- function calc_nn_p
  - double
  - double mun
  - double mup
- function calc_np_p
  - double
  - double nn
  - double mup
- function calc_dmu_delta
  - double
  - double nb
  - double delta
- function calc_musum_delta
  - double
  - double nb
  - double delta
- function calc_pressure_nb
  - double
  - double nb
  - double delta
- function calc_edensity_nb
  - double
  - double nb
  - double delta
# 
# Class eos_had_eden_base
#
class eos_had_eden_base abstract
- parent eos_had_base
# 
# Class eos_had_pres_base
#
class eos_had_pres_base abstract
- parent eos_had_base
# 
# Class eos_had_temp_base 
#
class eos_had_temp_base abstract
- parent eos_had_base
- function calc_temp_e
  - int
  - io fermion &n
  - io fermion &p
  - double T
  - io thermo &th
- function calc_temp_p
  - int
  - io fermion &n
  - io fermion &p
  - double T
  - io thermo &th
# 
# Class eos_had_temp_eden_base
#
class eos_had_temp_eden_base abstract
- parent eos_had_temp_base
# 
# Class eos_had_temp_pres_base
#
class eos_had_temp_pres_base abstract
- parent eos_had_temp_base
# 
# Class eos_had_skyrme
#
class eos_had_skyrme
- parent eos_had_temp_eden_base
- double t0
- double t1
- double t2
- double t3
- double x0
- double x1
- double x2
- double x3
- double alpha
- double a
- double b
- double W0
- double b4
- double b4p
- bool parent_method
- std::string reference
- o2scl::fermion_deriv_nr nrfd
# 
# Class eos_had_apr
#
class eos_had_apr
- parent eos_had_temp_eden_base
- int pion
- bool parent_method
# 
# Class eos_had_rmf
#
class eos_had_rmf
- parent eos_had_temp_pres_base
- size_t calc_e_steps
- bool calc_e_relative
- bool zm_mode
- int verbose
- double mnuc
- double ms
- double mw
- double mr
- double cs
- double cw
- double cr
- double b
- double c
- double zeta
- double xi
- double a1
- double a2
- double a3
- double a4
- double a5
- double a6
- double b1
- double b2
- double b3
- function get_fields
  - int
  - out double &sig
  - out double &ome
  - out double &rho
- function set_fields
  - int
  - io double &sig
  - io double &ome
  - io double &rho
# 
# Class eos_quark
#
class eos_quark abstract
- parent eos_base
# 
# Class eos_quark_bag
#
class eos_quark_bag
- parent eos_quark
- double bag_constant
# 
# Class eos_quark_njl
#
class eos_quark_njl
- parent eos_quark
- double B0
- double L
- double G
- double K
- double limit
- bool from_qq
# - o2scl::quark def_up
# - o2scl::quark def_down
# - o2scl::quark def_strange
# 
# Class eos_tov
#
class eos_tov abstract
- int verbose
- function has_baryons
  - bool
- function ed_from_pr
  - double
  - double pr
- function pr_from_ed
  - double
  - double ed
- function nb_from_ed
  - double
  - double ed
- function nb_from_pr
  - double
  - double pr
- function ed_from_nb
  - double
  - double nb
- function pr_from_nb
  - double
  - double nb
- function ed_nb_from_pr
  - void
  - double pr
  - out double &ed
  - out double &nb
# 
# Class eos_tov_buchdahl
#
class eos_tov_buchdahl
- parent eos_tov
- double Pstar
- double G_km_Msun
- function set_baryon_density
  - void
  - double nb
  - double ed
- function rad_from_gm
  - double
  - double gm
- function ed_from_r_gm
  - double
  - double r
  - double beta
- function pr_from_r_gm
  - double
  - double r
  - double beta
- function exp2lam_from_r_gm
  - double
  - double r
  - double beta
- function exp2phi_from_r_gm
  - double
  - double r
  - double beta
# 
# Class eos_tov_polytrope
#
class eos_tov_polytrope
- parent eos_tov
- function set_coeff_index
  - void
  - double coeff
  - double index
- function set_baryon_density
  - void
  - double nb
  - double ed
# 
# Class eos_tov_linear
#
class eos_tov_linear
- parent eos_tov
- function set_cs2_eps0
  - void
  - double cs2
  - double eps0
# 
# Class eos_tov_interp
#
class eos_tov_interp
- parent eos_tov
- bool err_nonconv
- std::vector<double> full_vece
- std::vector<double> full_vecp
- std::vector<double> full_vecnb
- function read_table
  - void
  - table_units<> &eos
  - std::string s_cole
  - std::string s_colp
  - std::string s_colnb [""]
- function default_low_dens_eos
  - void
- function sho11_low_dens_eos
  - void
- function s12_low_dens_eos
  - void
  - std::string model ["SLy4"]
  - bool external [false]
- function gcp10_low_dens_eos
  - void
  - std::string model ["BSk20"]
  - bool external [false]
- function ngl13_low_dens_eos
  - void
  - double L
  - std::string model ["PNM"]
  - bool external [false]
- function ngl13_low_dens_eos2
  - void
  - double S
  - double L
  - double nt
  - std::string fname [""]
- function no_low_dens_eos
  - void
# 
# Class tov_solve
#
class tov_solve
- no_cc  
- size_t buffer_size
- size_t max_table_size
- double mass
- double rad
- double bmass
- double gpot
- double last_rjw
- double last_f
- double domega_rat
- double pcent_max
- bool reformat_results
- double baryon_mass
- bool ang_vel
- bool gen_rel
- bool calc_gpot
- double step_min
- double step_max
- double step_start
- int verbose
- size_t max_integ_steps
- bool err_nonconv
- double pmax_default
- double prbegin
- double prend
- double princ
- double fixed_pr_guess
- double max_begin
- double max_end
- double max_inc
- function set_eos
  - void
  - eos_tov &eos
- function mvsr
  - int
- function fixed
  - int
  - double mass
  - double pmax [1.0e20]
- function fixed_pr
  - int
  - double pcent
  - double pmax [1.0e20]
- function max
  - int
- function get_results
  - shared_ptr table_units<>
# 
# Class tov_love
#
class tov_love
- int show_ode
- bool addl_testing
- bool err_nonconv
- table_units<> results
- double delta
- double eps
- shared_ptr table_units<> tab
- function calc_y
  - int
  - out double &yR
  - out double &beta
  - out double &k2
  - out double &lambda_km5
  - out double &lambda_cgs
  - bool tabulate
- function add_disc
  - void
  - double rd
- function clear_discs
  - void
- function calc_H
  - int
  - out double &yR
  - out double &beta
  - out double &k2
  - out double &lambda_km5
  - out double &lambda_cgs
# 
# ------------------------------------------------------
# 
# Class nstar_cold
#
class nstar_cold
- function set_eos
  - void
  - io eos_had_base &eos
- function calc_eos
  - int
  - double np_0 [0.0]
- function calc_nstar
  - int
- function fixed
  - int
  - double target_mass
- double pressure_dec_nb
- double allow_urca_nb
- double deny_urca_nb
- double acausal_nb
- double acausal_ed
- double acausal_pr
- tov_solve def_tov                             
- bool eos_neg
- int verbose
- function get_eos_results
  - shared_ptr table_units<>
- function get_tov_results
  - shared_ptr table_units<>
- double nb_start
- double nb_end
- double dnb
- size_t max_row  
- bool remove_rows  
- bool include_muons
- bool err_nonconv
# 
# ------------------------------------------------------
# 
# Class eos_nstar_rot
#
class eos_nstar_rot abstract
- parent eos_tov
- function enth_from_pr
  - double
  - double pr
- function pr_from_enth
  - double
  - double pr
- function enth_from_nb
  - double
  - double pr
# 
# ------------------------------------------------------
# 
# Class eos_nstar_rot_interp
#
class eos_nstar_rot_interp
- parent eos_nstar_rot
- function set_eos_fm
  - void
  - size_t n
  - vector<double> &eden
  - vector<double> &pres
  - vector<double> &nb
# 
# ------------------------------------------------------
# 
# Class eos_nstar_rot_C
#
class eos_nstar_rot_C
- parent eos_nstar_rot_interp
- function set
  - void
  - bool rns_constants [false]
# 
# ------------------------------------------------------
# 
# Class eos_nstar_rot_L
#
class eos_nstar_rot_L
- parent eos_nstar_rot_interp
- function set
  - void
  - bool rns_constants [false]
# 
# ------------------------------------------------------
# 
# Class nstar_rot
#
class nstar_rot
- bool err_nonconv
- function output_table
  - void
  - io table3d &t
- function set_eos
  - void
  - io eos_nstar_rot &eos
- function polytrope_eos
  - void
  - double index
- function fix_cent_eden_axis_rat
  - void
  - double cent_eden
  - double axis_rat
  - bool use_guess [false]
- function fix_cent_eden_grav_mass
  - void
  - double cent_eden
  - double grav_mass
- function fix_cent_eden_bar_mass
  - void
  - double cent_eden
  - double bar_mass
- function fix_cent_eden_with_kepler
  - int
  - double cent_eden
- function fix_cent_eden_ang_vel
  - int
  - double cent_eden
  - double ang_vel
- function fix_cent_eden_ang_mom
  - int
  - double cent_eden
  - double ang_mom
- function fix_cent_eden_non_rot
  - int
  - double cent_eden
- function fix_cent_eden_with_kepler_alt
  - int
  - double cent_eden
  - bool use_guess [false]
- function fix_cent_eden_grav_mass_alt
  - int
  - double cent_eden
  - double grav_mass
  - bool use_guess [false]
- function fix_cent_eden_bar_mass_alt
  - int
  - double cent_eden
  - double bar_mass
  - bool use_guess [false]
- function fix_cent_eden_ang_vel_alt
  - int
  - double cent_eden
  - double ang_vel
  - bool use_guess [false]
- function fix_cent_eden_ang_mom_alt
  - int
  - double cent_eden
  - double ang_mom
  - bool use_guess [false]
- function constants_rns
  - void
- function constants_o2scl
  - void
- double eq_radius_tol_rel
- double alt_tol_rel
- int verbose
- double e_center
- double r_ratio
- double r_e
- double r_p
- double s_p
- double s_e
- double R_e
- double Mass_p
- double Mass
- double Mass_0
- double J
- double Omega
- double T
- double I
- double W
- double Z_p
- double Z_f
- double Z_b
- double Omega_K
- double eccentricity
# I'm not sure if ubvectors work yet
- ubvector v_plus  
- ubvector v_minus
- double vel_plus  
- double vel_minus
- double h_plus
- double h_minus
- double Omega_plus
- double om_over_Om
- double mass_quadrupole
- double cf
- double C
- double G
- double MSUN
- double KAPPA
- double MB
- double KSCALE
- double PI
# 
# ------------------------------------------------------
#
# Class nucleus_rmf
#
class nucleus_rmf
- function run_nucleus
  - int
  - int nucleus_Z
  - int nucleus_N
  - int unocc_Z
  - int unocc_N
- function get_profiles
  - shared_ptr table_units<>
- function get_chden
  - shared_ptr table_units<>
- double stens
- double rnrp
- double rnrms
- double rprms
- double etot
- double r_charge  
- double r_charge_cm
#- eos_had_rmf def_rmf
#- function set_eos
#  - int
#  - io eos_had_rmf &r
#
# ------------------------------------------------------
#
# Class nucmass_ldrop
#                              
class nucmass_ldrop
- parent nucmass_fit_base
- double n1
- double n0
- double surften
- double coul_coeff
- double nn
- double np
- double Rn
- double Rp
- double surf
- double bulk
- double coul
- bool large_vals_unphys
- function mass_excess_d
  - double
  - double Z
  - double N
- function mass_excess
  - double
  - int Z
  - int N
- function binding_energy_densmat
  - double
  - double Z
  - double N
  - double npout
  - double nnout
  - double neout
  - double T
- eos_had_skyrme def_had_eos
- fermion def_neutron
- fermion def_proton
- thermo th
- function set_n_and_p
  - void
  - io fermion &un
  - io fermion &up
- function set_eos_had_temp_base
  - int
  - io eos_had_temp_base &uhe
#
# ------------------------------------------------------
#
# Class nucmass_ldrop_skin
#                              
class nucmass_ldrop_skin
- parent nucmass_ldrop
- bool full_surface
- bool new_skin_mode
- double doi
- double ss
- double pp
- double a0
- double a2
- double a4
- bool rel_vacuum
- double Tchalf
#
# ------------------------------------------------------
#
# Class nucmass_ldrop_pair
#                              
class nucmass_ldrop_pair
- parent nucmass_ldrop_skin
- double Epair
- double pair
#
# ------------------------------------------------------
#
# Class nucmass_ldrop_shell
#                              
class nucmass_ldrop_shell
- parent nucmass_ldrop_pair
#
# ------------------------------------------------------
#
# Class nucmass_frdm_shell
#                              
class nucmass_frdm_shell
- parent nucmass_frdm
#
# ------------------------------------------------------
#
# Class nucleus_bin
#                              
class nucleus_bin
- nucmass_ame ame16
- nucmass_ame ame20exp
- nucmass_ame ame20round
- nucmass_ame ame95rmd
- nucmass_ame ame03round
- nucmass_ame ame03
- nucmass_ame ame95exp
- nucmass_ame ame12
- nucmass_gen ddme2 
- nucmass_gen ddmed 
- nucmass_gen ddpc1 
- nucmass_gen nl3s 
- nucmass_gen sly4 
- nucmass_gen skms 
- nucmass_gen skp 
- nucmass_gen sv_min 
- nucmass_gen unedf0 
- nucmass_gen unedf1 
- nucmass_mnmsk m95 
- nucmass_mnmsk m16 
- nucmass_ktuy kt 
- nucmass_ktuy kt2 
- nucmass_wlw wlw1 
- nucmass_wlw wlw2 
- nucmass_wlw wlw3 
- nucmass_wlw wlw4 
- nucmass_wlw wlw5 
- nucmass_sdnp sdnp1 
- nucmass_sdnp sdnp2 
- nucmass_sdnp sdnp3 
- nucmass_dz_table dz 
- nucmass_hfb hfb2 
- nucmass_hfb hfb8 
- nucmass_hfb hfb14 
- nucmass_hfb hfb14_v0 
- nucmass_hfb_sp hfb17 
- nucmass_hfb_sp hfb21 
- nucmass_hfb_sp hfb22 
- nucmass_hfb_sp hfb23 
- nucmass_hfb_sp hfb24 
- nucmass_hfb_sp hfb25 
- nucmass_hfb_sp hfb26 
- nucmass_hfb_sp hfb27 
# 
# HDF functions
#
function skyrme_load
- void
- io eos_had_skyrme &sk
- std::string model
- bool external [false]
- int verbose [0]
function rmf_load
- void
- io eos_had_rmf &rmf
- std::string model
- bool external [false]

