/*
  -------------------------------------------------------------------
  
  This file is part of O2scl. It has been adapted from RNS v1.1d
  written by N. Stergioulas and S. Morsink. The modifications made in
  this version from the original are copyright (C) 2015-2021, Andrew
  W. Steiner.
  
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
/*
  -------------------------------------------------------------------
  Relativistic models of rapidly rotating compact stars,
  using tabulated or polytropic equations of state.
  
  Author:  Nikolaos Stergioulas
  
  Current Address:
  
  Department of Physics
  University of Wisconsin-Milwaukee
  PO Box 413, Milwaukee, WI 53201, USA
  
  E-mail: niksterg@csd.uwm.edu, or
  niksterg@pauli.phys.uwm.edu
  
  Version: 1.1
  
  Date:    June, 1995
  
  Changes made to code by Sharon Morsink
   
  03-03-97: Corrected the units for polytropic stars
  10-28-98: Added the star's quadrupole moment to the output.
  
  References:
  KEH : H. Komatsu, Y. Eriguchi and I. Hachisu, Mon. Not. R. astr. Soc. 
  (1989) 237, 355-379.
  CST : G. Cook, S. Shapiro and S. Teukolsky, Ap. J (1992) 398, 203-223.
  
  -------------------------------------------------------------------
*/
/** \file nstar_rot.h
    \brief File defining \ref o2scl::nstar_rot
*/
#ifndef NSTAR_ROT_H
#define NSTAR_ROT_H

#include <cmath>
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/search_vec.h>
#include <o2scl/test_mgr.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/lib_settings.h>
#include <o2scl/interp.h>
#include <o2scl/eos_tov.h>
#include <o2scl/table3d.h>
#include <o2scl/tensor.h>
#include <o2scl/mroot_hybrids.h>

#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

namespace o2scl {
  
  /** \brief An EOS for \ref nstar_rot

      The \ref nstar_rot class uses the functions
      - \ref eos_tov::pr_from_ed()
      - \ref eos_tov::ed_from_pr()
      - \ref eos_tov::nb_from_ed()
      - \ref enth_from_pr()
      - \ref pr_from_enth()

      It is assumed that all of these functions operate with
      the following unit scheme
      - \c ed is the energy density in 
      units of \f$ 10^{15}~\mathrm{g}/\mathrm{cm}^3 \f$
      - \c pr is the pressure in 
      units of \f$ 10^{15}~\mathrm{g}/\mathrm{cm}^3 \f$
      - \c enth is the pseudo-enthalpy defined by
      \f[
      \hat{h}_{\mathrm{RNS}} = 
      10^{-8} \left[\frac{\mu}{\mu(P=0)}\right] 
      \f]
      - \c nb is the baryon density in 
      units of \f$ 1/\mathrm{cm}^3 \f$
   */
  class eos_nstar_rot : public eos_tov {
    
  public:

    /** \brief From the pressure, return the pseudo-enthalpy
     */
    virtual double enth_from_pr(double pr)=0;

    /** \brief From the pseudo-enthalpy, return the pressure
     */
    virtual double pr_from_enth(double enth)=0;

    /** \brief From the baryon density, return the pseudo-enthalpy
     */
    virtual double enth_from_nb(double nb)=0;
  };
  
  /** \brief Create a tabulated EOS for \ref nstar_rot using interpolation

      The tabulated EOS is stored in the following units
      - \ref log_e_tab is the log of the energy density in 
      units of \f$ 10^{15}~\mathrm{g}/\mathrm{cm}^3 \f$
      - \ref log_p_tab is the log of the pressure in 
      units of \f$ 10^{15}~\mathrm{g}/\mathrm{cm}^3 \f$
      - \ref log_h_tab is the log of the pseudo-enthalpy
      - \ref log_n0_tab is the log of the baryon density in 
      units of \f$ 1/\mathrm{cm}^3 \f$

      The user specified EOS must have a monotonically increasing
      baryon chemical potential, as guaranteed, for example, by 
      \ref o2scl::eos_nstar_rot_C.

      \future Replace arrays with vectors and fix unit-indexing.
   */
  class eos_nstar_rot_interp : public eos_nstar_rot {
    
  protected:

    /// Array search object
    o2scl::search_vec<double *> sv;
    
    /// Search in array \c x of length \c n for value \c val
    int new_search(int n, double *x, double val);
    
    /** \brief number of tabulated EOS points */
    int n_tab;                           
    /** \brief rho points in tabulated EOS */
    double log_e_tab[200];               
    /** \brief p points in tabulated EOS */
    double log_p_tab[200];               
    /** \brief h points in EOS file */
    double log_h_tab[200];               
    /** \brief number density in EOS file */  
    double log_n0_tab[200];              

    /// \name Constants
    //@{
    /** \brief Speed of light in vacuum (in CGS units) */ 
    double C;
    /** \brief Gravitational constant (in CGS units) */ 
    double G;
    /** \brief Square of length scale in CGS units, 
        \f$ \kappa \equiv 10^{-15} c^2/G \f$
    */
    double KAPPA;
    /** \brief The value \f$ 10^{-15}/c^2 \f$ */
    double KSCALE;
    //@}

    /** \brief Cache for interpolation
     */
    int n_nearest;
    
    /** \brief Driver for the interpolation routine. 
        
        First we find the tab. point nearest to xb, then we
        interpolate using four points around xb.

         Note that this version, since the EOS arrays are now 0 indexed,
         is different than the nstar_rot version. There is also a bit of
         extra arithmetic in this function below which is probably
         unnecessary which is left over from the shift from unit to
         zero-indexing. 
    */
    double interp(double xp[], double yp[], int np, double xb);

  public:
    
    eos_nstar_rot_interp();

    /** \brief Set the EOS from four vectors in the native unit system
     */
    template<class vec1_t, class vec2_t, class vec3_t, class vec4_t>
      void set_eos_native(vec1_t &eden, vec2_t &pres, vec3_t &enth,
                          vec4_t &nb) {
      
      double C=o2scl_cgs::speed_of_light;
      double G=o2scl_cgs::gravitational_constant;
      double KAPPA=1.0e-15*C*C/G;
      double KSCALE=KAPPA*G/(C*C*C*C);
  
      n_tab=eden.size();

      // Note that conv1*C*C*KSCALE is identical to conv2*KSCALE.
      
      for(int i=0;i<n_tab;i++) {
        // Convert from g/cm^3 to 1.0e15 grams/cm^3 and take the log
        log_e_tab[i]=log10(eden[i]*C*C*KSCALE);
        // Convert from dyne/cm^2 to 1.0e15 grams/cm^3 and take the log
        log_p_tab[i]=log10(pres[i]*KSCALE);
        // Convert from cm^2/s^2 to a unitless quantity and take the log
        log_h_tab[i]=log10(enth[i]/(C*C));
        // Take the log of a quantity in units of 1/cm^3
        log_n0_tab[i]=log10(nb[i]);
      }
      
      return;
    }

    /** \brief Set the EOS from energy density, pressure, and
        baryon density stored in powers of \f$ \mathrm{fm} \f$ .

        This function appends the default RNS crust and thus
        the \c eden, \c pres, and \c nb vectors should not
        contain points corresponding to the EOS for baryon
        densities smaller than \f$ 0.07463~\mathrm{fm}^{-3} \f$.
    */
    template<class vec1_t, class vec2_t, class vec3_t>
      void set_eos_fm(size_t n, vec1_t &eden, vec2_t &pres, vec3_t &nb) {

      // The number of points in the hard-coded crust EOS
      static const int n_crust=78;
      
      if (n>200-n_crust) {
        O2SCL_ERR2("Too many EOS points in ",
                   "nstar_rot::set_eos().",o2scl::exc_einval);
      }

      // Conversion factor for energy density
      double conv1=o2scl_settings.get_convert_units().convert
        ("1/fm^4","g/cm^3",1.0);
      // Conversion factor for pressure
      double conv2=o2scl_settings.get_convert_units().convert
        ("1/fm^4","dyne/cm^2",1.0);

      n_tab=n+n_crust;

      /* Use the original RNS crust from eosC, except for the enthalpy
         which is computed by hand below. This appears to work better
         than the default O2scl crust, and this may have to do with
         the fact that the default O2scl crust has decreasing mu with
         increasing density at low densities.

         These columns are:
         - energy density in g/cm^3
         - pressure in dyne/cm^2
         - baryon density in 1/cm^3

         The first few rows are close to the FMT EOS stored in fmt49.o2.
      */
      double nst_arr[n_crust][3]={
        {7.800e+00,1.010e+08,4.698795180722962e+24},
        {7.860e+00,1.010e+09,4.734939759036205e+24},
        {7.900e+00,1.010e+10,4.759036144578364e+24},
        {8.150e+00,1.010e+11,4.909638554215315e+24},
        {1.160e+01,1.210e+12,6.987951807098076e+24},
        {1.640e+01,1.400e+13,9.879518070489597e+24},
        {4.510e+01,1.700e+14,2.716867462904601e+25},
        {2.120e+02,5.820e+15,1.277108403508764e+26},
        {1.150e+03,1.900e+17,6.927709645088004e+26},
        {1.044e+04,9.744e+18,6.289148562640985e+27},
        {2.622e+04,4.968e+19,1.579513843816999e+28},
        {6.587e+04,2.431e+20,3.968050678245718e+28},
        {1.654e+05,1.151e+21,9.963748410271617e+28},
        {4.156e+05,5.266e+21,2.503563031417219e+29},
        {1.044e+06,2.318e+22,6.288917532113082e+29},
        {2.622e+06,9.755e+22,1.579410809416864e+30},
        {6.588e+06,3.911e+23,3.968207649843547e+30},
        {8.293e+06,5.259e+23,4.995116726219748e+30},
        {1.655e+07,1.435e+24,9.967984755458204e+30},
        {3.302e+07,3.833e+24,1.988624478073943e+31},
        {6.589e+07,1.006e+25,3.967807406359445e+31},
        {1.315e+08,2.604e+25,7.917691186982454e+31},
        {2.624e+08,6.676e+25,1.579648605894070e+32},
        {3.304e+08,8.738e+25,1.988876577393412e+32},
        {5.237e+08,1.629e+26,3.152005155076383e+32},
        {8.301e+08,3.029e+26,4.995278531652059e+32},
        {1.045e+09,4.129e+26,6.287859551784352e+32},
        {1.316e+09,5.036e+26,7.917701445937253e+32},
        {1.657e+09,6.860e+26,9.968319738044036e+32},
        {2.626e+09,1.272e+27,1.579408507997411e+33},
        {4.164e+09,2.356e+27,2.503766293549853e+33},
        {6.601e+09,4.362e+27,3.967852390467774e+33},
        {8.312e+09,5.662e+27,4.995474308724729e+33},
        {1.046e+10,7.702e+27,6.285277578607203e+33},
        {1.318e+10,1.048e+28,7.918132634568090e+33},
        {1.659e+10,1.425e+28,9.964646988214994e+33},
        {2.090e+10,1.938e+28,1.255052800774333e+34},
        {2.631e+10,2.503e+28,1.579545673652798e+34},
        {3.313e+10,3.404e+28,1.988488463504033e+34},
        {4.172e+10,4.628e+28,2.503379640977065e+34},
        {5.254e+10,5.949e+28,3.151720931652274e+34},
        {6.617e+10,8.089e+28,3.968151735612910e+34},
        {8.332e+10,1.100e+29,4.994995310195290e+34},
        {1.049e+11,1.495e+29,6.286498800006776e+34},
        {1.322e+11,2.033e+29,7.919521253825185e+34},
        {1.664e+11,2.597e+29,9.964341016667146e+34},
        {1.844e+11,2.892e+29,1.104024323001462e+35},
        {2.096e+11,3.290e+29,1.254619611126682e+35},
        {2.640e+11,4.473e+29,1.579588892045295e+35},
        {3.325e+11,5.816e+29,1.988565738933728e+35},
        {4.188e+11,7.538e+29,2.503561780689725e+35},
        {4.299e+11,7.805e+29,2.569780082714395e+35},
        {4.460e+11,7.890e+29,2.665824694449485e+35},
        {5.228e+11,8.352e+29,3.123946525953616e+35},
        {6.610e+11,9.098e+29,3.948222384313103e+35},
        {7.964e+11,9.831e+29,4.755697604312120e+35},
        {9.728e+11,1.083e+30,5.807556544067428e+35},
        {1.196e+12,1.218e+30,7.138304213736713e+35},
        {1.471e+12,1.399e+30,8.777653631971616e+35},
        {1.805e+12,1.683e+30,1.076837272716171e+36},
        {2.202e+12,1.950e+30,1.313417953138369e+36},
        {2.930e+12,2.592e+30,1.747157788902558e+36},
        {3.833e+12,3.506e+30,2.285004034820638e+36},
        {4.933e+12,4.771e+30,2.939983642627298e+36},
        {6.248e+12,6.481e+30,3.722722765704268e+36},
        {7.801e+12,8.748e+30,4.646805278760175e+36},
        {9.611e+12,1.170e+31,5.723413975645761e+36},
        {1.246e+13,1.695e+31,7.417258934884369e+36},
        {1.496e+13,2.209e+31,8.902909532230595e+36},
        {1.778e+13,2.848e+31,1.057801059193907e+37},
        {2.210e+13,3.931e+31,1.314278492046241e+37},
        {2.988e+13,6.178e+31,1.775810743961577e+37},
        {3.767e+13,8.774e+31,2.237518046976615e+37},
        {5.081e+13,1.386e+32,3.015480061626022e+37},
        {6.193e+13,1.882e+32,3.673108933334910e+37},
        {7.732e+13,2.662e+32,4.582250451016437e+37},
        {9.826e+13,3.897e+32,5.817514573447143e+37},
        {1.262e+14,5.861e+32,7.462854442694524e+37}};

      double mu_start;
      // Note that there is no c^2 needed in the computation of the
      // enthalpy as the original code removes it.
      for(size_t i=0;i<n_crust;i++) {

        // Convert from g/cm^3 to 1.0e15 grams/cm^3 and take the log
        log_e_tab[i]=log10(nst_arr[i][0]*C*C*KSCALE);

        // Convert from dyne/cm^2 to 1.0e15 grams/cm^3 and take the log
        log_p_tab[i]=log10(nst_arr[i][1]*KSCALE);
        
        // Convert the first term from g/cm^3 to 1/fm^4 and the second
        // term from dyne/cm^2 to 1/fm^4 and then divide by a quantity
        // in units of 1/fm^3 to obtain mu in units of 1/fm.
        double mu=(nst_arr[i][0]/conv1+nst_arr[i][1]/conv2)/
          nst_arr[i][2]*1.0e39;
        if (i==0) {
          mu_start=mu;
        } else {
          log_h_tab[i]=log10(log(mu/mu_start));
        }
        // Take the log of a quantity in units of 1/cm^3
        log_n0_tab[i]=log10(nst_arr[i][2]);
        //std::cout << "1." << log_e_tab[i] << " " 
        //<< log_p_tab[i] << " " << log_n0_tab[i] << " " << mu
        //<< std::endl;
      }

      // RNS chooses the psuedo-enthalpy at the smallest density
      // to be 1.0
      log_h_tab[0]=log10(1.0/(C*C));

      // Note that conv1*C*C*KSCALE is identical to conv2*KSCALE.

      for(size_t i=0;i<n;i++) {
        // Convert from 1/fm^4 to 1.0e15 grams/cm^3 and take the log
        log_e_tab[i+n_crust]=log10(eden[i]*conv1*C*C*KSCALE);
        // Convert from 1/fm^4 to 1.0e15 grams/cm^3 and take the log
        log_p_tab[i+n_crust]=log10(pres[i]*conv2*KSCALE);
        // Take the log of a unitless quantity
        log_h_tab[i+n_crust]=log10(log((eden[i]+pres[i])/nb[i]/
                                         mu_start));
        // Convert from 1/fm^3 to 1/cm^3 and take the log
        if (nb[i]<0.07643) {
          O2SCL_ERR2("EOS overlaps with crust in ",
                     "eos_nstar_rot_interp::set_eos_fm().",
                     o2scl::exc_einval);
        }
        log_n0_tab[i+n_crust]=log10(nb[i]*1.0e39);
        double mu=(eden[i]+pres[i])/nb[i];
        //std::cout << "1b." << log_e_tab[i+n_crust] << " " 
        //<< log_p_tab[i+n_crust] << " "
        //<< log_n0_tab[i+n_crust] << " " << mu
        //<< std::endl;
      }

      return;
    }

    /** \brief Set the EOS from energy density, pressure, and
        baryon density stored in powers of \f$ \mathrm{fm} \f$ .

        This function presumes that the user-specified vectors that
        contain the crust EOS. Note that the user-specified EOS must
        have a monotonically increasing baryon chemical potential.
    */
    template<class vec1_t, class vec2_t, class vec3_t>
      void set_eos_crust_fm(size_t n, vec1_t &eden, vec2_t &pres,
                            vec3_t &nb) {

      if (n>200) {
        O2SCL_ERR2("Too many EOS points in ",
                   "nstar_rot::set_eos().",o2scl::exc_einval);
      }
      
      n_tab=n;
      
      // Conversion factor for energy density
      double conv1=o2scl_settings.get_convert_units().convert
        ("1/fm^4","g/cm^3",1.0);
      // Conversion factor for pressure
      double conv2=o2scl_settings.get_convert_units().convert
        ("1/fm^4","dyne/cm^2",1.0);
      
      // Note that conv1*C*C*KSCALE is identical to conv2*KSCALE.
      
      double mu_start;
      for(size_t i=0;i<n;i++) {
        // Convert from 1/fm^4 to 1.0e15 grams/cm^3 and take the log
        log_e_tab[i]=log10(eden[i]*conv1*C*C*KSCALE);
        // Convert from 1/fm^4 to 1.0e15 grams/cm^3 and take the log
        log_p_tab[i]=log10(pres[i]*conv2*KSCALE);
        
        double mu=(eden[i]+pres[i])/nb[i];
        if (i==0) {
          mu_start=mu;
        } else {
          log_h_tab[i]=log10(log(mu/mu_start));
          if (!std::isfinite(log_h_tab[i])) {
            O2SCL_ERR2("Non-increasing chemical potential ",
                      "in set_eos_crust_fm().",o2scl::exc_einval);
          }
        }
        
        log_n0_tab[i]=log10(nb[i]*1.0e39);
        //std::cout << "2." << log_e_tab[i] << " " 
        //<< log_p_tab[i] << " " << log_n0_tab[i] << " " << mu
        //<< std::endl;
      }
      
      // RNS chooses the pseudo-enthalpy at the smallest density
      // to be 1.0
      log_h_tab[0]=log10(1.0/(C*C));
      
      return;
    }

    /** \brief From the pressure, return the energy density
     */
    virtual double ed_from_pr(double pr);

    /** \brief From the energy density, return the pressure
     */
    virtual double pr_from_ed(double ed);

    /** \brief From the pressure, return the baryon density
     */
    virtual double nb_from_pr(double pr);

    /** \brief From the baryon density, return the pressure
     */
    virtual double pr_from_nb(double nb);

    /** \brief From the baryon density, return the energy density
     */
    virtual double ed_from_nb(double nb);

    /** \brief From the energy density, return the baryon density
     */
    virtual double nb_from_ed(double ed);

    /** \brief From the pressure, return the enthalpy
     */
    virtual double enth_from_pr(double pr);

    /** \brief From the baryon density, return the enthalpy
     */
    virtual double enth_from_nb(double nb);

    /** \brief From the enthalpy, return the pressure
     */
    virtual double pr_from_enth(double enth);

    /** \brief Given the pressure, compute the energy and number densities
    */
    virtual void ed_nb_from_pr(double pr, double &ed, double &nb);

    /** \brief Output EOS table to screen

        \comment 
        This is mostly for testing and should be replaced once
        the arrays are properly replaced with vectors
        \endcomment
     */
    void output() {
      for(int i=n_tab-1;i>=0;i--) {
        std::cout << log_e_tab[i] << " " << log_p_tab[i] << " "
                  << log_h_tab[i] << " " << log_n0_tab[i] << std::endl;
      }
      std::cout << std::endl;
      return;
    }

    /** \brief Output EOS to an HDF5 file
     */
    void output_table(std::string fname) {
      o2scl_hdf::hdf_file hf;
      o2scl::table_units<> t;
      t.line_of_names("log_e log_p log_h log_n0");
      for(int i=0;i<n_tab;i++) {
        double line[4]={log_e_tab[i],log_p_tab[i],log_h_tab[i],
                        log_n0_tab[i]};
        t.line_of_data(4,line);
      }
      hf.open_or_create(fname);
      hdf_output(hf,t,"eos");
      hf.close();
      return;
    }
    
  };
  
  /** \brief Tabulated EOS for \ref nstar_rot from Bethe et al. 1974

      \verbatim embed:rst
      From [Bethe74]_.
      \endverbatim
   */
  class eos_nstar_rot_C : public eos_nstar_rot_interp {
  public:
    eos_nstar_rot_C(bool rns_constants=false);
  };
  
  /** \brief Tabulated EOS for \ref nstar_rot from Pandharipande75 et al. 1975

      \verbatim embed:rst
      From [Pandharipande75]_.
      \endverbatim
   */
  class eos_nstar_rot_L : public eos_nstar_rot_interp {
  public:
    eos_nstar_rot_L(bool rns_constants=false);
  };
  
  /** \brief Rotating neutron star class based on RNS v1.1d from
      N. Stergioulas et al.
      
      \note This class is still experimental.

      Several changes have been made to the original code. The code
      using Numerical Recipes has been removed and replaced with an
      equivalent based on GSL and \o2. The overall interface has
      been changed and some code has been updated with C++
      equivalents.

      <b>Initial guess</b>

      The original RNS code suggests that the initial guess is
      typically a star with a smaller angular momentum.

      <b>References</b> 

      The original RNS v1.1d can be obtained from
      http://www.gravity.phys.uwm.edu/rns/, and you may find Nick
      Stergioulas's web page http://www.astro.auth.gr/~niksterg/, or
      Sharon Morsink's page http://fermi.phys.ualberta.ca/~morsink/
      useful. 

      \verbatim embed:rst
      See [Bonazzola73]_, [Bonazzola94]_, [Cook92]_,
      [Cook94]_, [Friedman88]_, [Gourgoulhon94]_, 
      [Komatsu89]_, [Laarakkers99]_, [Nozawa98]_, [Stergioulas95]_,
      and [Stergioulas03]_.
      \endverbatim

      \verbatim embed:rst
      .. todo:: 

         In class nstar_rot:

         - Better documentation is needed everywhere.
         - Test the resize() function
         - It appears that KAPPA and KSCALE contant an arbitrary
           constant, try changing it and see if we get identical
           results. Try to ensure that the values are consistent between
           the eos_nstar_rot class and the nstar_rot class.
         - Variables r_is_gp, p_center, h_center, and others only
           occur in spherical_star(), integrate(), and make_center(),
           and can be moved to function parameters or otherwise
           reorganized.
         - Directly compare spherical_star() output with 
           tov_solve results

      \endverbatim

      \future Consider moving the int_z() algorithm to vector_derint.h
      \future Remove the unit-indexed arrays everywhere.
      \future Try moving some of the storage to the heap?
      \future Some of the arrays seem larger than necessary.
      \future The function \ref o2scl::nstar_rot::new_search() is
      inefficient because it has to handle the boundary conditions
      separately. This could be improved.
      \future Make the solvers more robust. The ang_vel() and ang_vel_alt()
      functions appear particularly unstable.

      <b>Equation of State</b> 

      The thermodynamic enthalpy is 
      \f[
      H = E + P V
      \f]
      and the specific enthalpy is often written 
      \f$ h = H/M \f$ 
      where \f$ M \f$ is some mass scale. In GR this is often chosen
      to be a generic "baryon mass", \f$ m_B \f$ so that 
      \f[
      h = \frac{\varepsilon + P}{\rho}
      \f]
      where \f$ \varepsilon \equiv E/V \f$ and \f$ \rho \equiv m_B/V
      \f$. The pseudo-enthalpy is typically defined by 
      \f[ 
      d\hat{h} = dP/(P+\varepsilon) 
      \f] 
      (sometimes denoted as \f$ H \f$ or \f$ h
      \f$ and called `h` in the code, but referred to here as \f$
      \hat{h} \f$ to avoid confusion with the enthalpy and
      specific enthalpy defined
      above). Additionally, this quantity is sometimes defined
      with an additional factor of \f$ c^2 \f$.
      At \f$ T=0 \f$ and presuming only one conserved charge, 
      one can use the Gibbs-Duhem relation, \f$
      dP = (\varepsilon + P)/\mu d \mu \f$ to write 
      \f[ 
      d \hat{h} = \frac{d \mu}{\mu}
      \f]
      \f[
      \hat{h} = \ln \mu + C
      \f]
      The constant \f$ C \f$ is often chosen to be 
      \f$ C \equiv - \ln \mu(P=0) \f$, so that 
      \f[
      \hat{h} = \ln \frac{\mu}{\mu(P=0)}
      \f]
      Typically, \f$ \mu(P=0) \f$ is around 931 MeV. The \ref
      eos_nstar_rot_interp class takes an EOS tabulated in powers of
      \f$ \mathrm{fm}^{-1} \f$ and recasts it into a form which can be
      used by \ref nstar_rot.

      <b>Coordinate system</b>

      The space is mapped to coordinates \f$ (s,\mu) \f$ where
      \f[
      r = r_{\mathrm{eq}} \left( \frac{s}{1-s} \right)
      \f]
      and \f$ r_{\mathrm{eq}} \f$ is the equatorial radius,
      thus 
      \f[
      s = \frac{r}{r+r_{\mathrm{eq}}}
      \f]
      (This is Eq. 48 in Stergioulas' thesis)

      <b>Draft documentation</b> 

      For spherical stars, the isotropic radius \f$ r_{\mathrm{is}}
      \f$ is defined by
      \f[
      \frac{d r}{d r_{\mathrm{is}}} = 
      \left(\frac{r}{r_{\mathrm{is}}}\right)
      \left(1 - 2 \frac{m}{r}\right)^{1/2}
      \f]

      <b>Quadrupole moments</b>

      \verbatim embed:rst
      Quadrupole moments computed using the method in [Laarakkers99]_. 
      \endverbatim

      <b>Axisymmetric Instability</b>

      \verbatim embed:rst
      [Friedman88]_ shows that a secular axisymmetric instability
      sets in when the mass becomes maximum along a sequence of
      constant angular momentum. Equivalently, [Cook92]_ shows that
      the instability occurs when the angular momentum becomes minimum
      along a sequence of constant rest mass.
      \endverbatim

      \verbatim embed:rst
      A GR virial theorem for a stationary and axisymmetric system was
      found in [Bonazzola73]_. A more general two-dimensional virial
      identity was found in [Bonazzola94]_. The three-dimensional
      virial identity found in [Gourgoulhon94]_ is a generalization
      of the Newtonial virial theorem.
      \endverbatim

      Using the stationary and axisymmetric metric ( \f$ G = c = 1 \f$
      )
      \f[
      ds^2 = - e^{\gamma+\rho} dt^2 + e^{2 \alpha} \left( dr^2 + 
      r^2 d\theta^2 \right) + e^{\gamma-\rho} r^2 \sin^2 \theta
      ( d \phi - \omega dt) ^2
      \f]
      one solves for the four metric functions \f$ \rho(r,\theta) \f$,
      \f$ \gamma(r,\theta) \f$, \f$ \alpha(r,\theta) \f$ and \f$
      \omega(r,\theta) \f$ .

      It is assumed that matter is a perfect fluid, and the 
      stress-energy tensor is
      \f[
      T^{\mu \nu} = \left( \rho_0 + \rho_i + P \right) u^{\mu} u^{\nu}
      + P g^{\mu \nu}
      \f]
      
      Einstein's field equations imply four field equations for
      a specified rotation law,
      \f[
      u^{t} u_{\phi} = F(\Omega) 
      \f]
      for some function \f$ F(\omega) \f$ .

      \verbatim embed:rst
      Using Eq. (27) in [Cook92]_, one can write
      \endverbatim
      \f[
      \rho(s,\mu) = - e^{-\gamma/2} \sum_{n=0}^{\infty}
      P_{2n}(\mu) \int_0^{1}~ds^{\prime} \int_0^1~d \mu 
      f_{\rho}(n,s,s^{\prime}) P_{2n}{\mu^{\prime}} 
      \tilde{S}(s^{\prime},\mu^{\prime})
      \f]
      where the function \f$ f_{\rho} \f$ is defined by
      \f[
      f_{\rho} \equiv \Theta(s^{\prime}-s)
      \left(\frac{1-s}{s}\right)^{2 n+1} \left[\frac{s^{\prime
      2n}}{(1-s^{\prime})^{2n+2}}\right] + \Theta(s^{\prime}-s)
      \left(\frac{1-s}{s}\right)^{2 n+1} \left[\frac{s^{\prime
      2n}}{(1-s^{\prime})^{2n+2}}\right]
      \f]
      This function is stored in \ref f_rho . Similar 
      definitions are made for \ref f_gamma and \ref f_omega .

      The Keplerial orbit at the equator is 
      \f[
      \Omega_K = \frac{\omega^{\prime}}{2 \psi^{\prime}} ...
      \f]
      \verbatim embed:rst
      (eq. 31 in [Stergioulas03]_)
      \endverbatim
      
  */
  class nstar_rot {
  
  public:    
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::range ub_range;
    typedef boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<double> > ubvector_range;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    /// The number of grid points in the \f$ \mu \f$ direction
    int MDIV;
    /// The number of grid points in the \f$ s \f$ direction
    int SDIV;
    /// The number of Legendre polynomials
    int LMAX;

    /// Resize the grid
    void resize(int MDIV_new, int SDIV_new, int LMAX_new,
                int RDIV_new);
    
    /** \brief Default solver
     */
    o2scl::mroot_hybrids<> def_mroot;

    /** \brief Set new solver
     */
    int set_solver(o2scl::mroot<> &m) {
      mrootp=&m;
      return 0;
    }
    
  protected:

    /// Solver
    o2scl::mroot<> *mrootp;

    /// Solve for the Keplerian velocity
    int solve_kepler(size_t nv, const ubvector &x, ubvector &y);

    /// Solve for the gravitational mass
    int solve_grav_mass(size_t nv, const ubvector &x, ubvector &y,
                        double grav_mass);
    
    /// Solve for the gravitational mass
    int solve_bar_mass(size_t nv, const ubvector &x, ubvector &y,
                        double bar_mass);
    
    /// Solve for the gravitational mass
    int solve_ang_vel(size_t nv, const ubvector &x, ubvector &y,
                        double ang_vel);
    
    /// Solve for the gravitational mass
    int solve_ang_mom(size_t nv, const ubvector &x, ubvector &y,
                        double ang_mom);
    
    /** \brief Subclass of \ref nstar_rot which specifies the function
        to invert a polytropic EOS
    */
    class polytrope_solve {

    protected:

      /** \brief The polytropic index
       */
      double _Gamma_P;

      /** \brief The energy density
       */
      double _ee;

    public:

      /** \brief Create a function object with specified 
          polytropic index and ?
      */
      polytrope_solve(double Gamma_P, double ee) {
        _Gamma_P=Gamma_P;
        _ee=ee;
      }
      
      /** \brief The function
       */
      double operator()(double rho0) {
        return pow(rho0,_Gamma_P)/(_Gamma_P-1.0)+rho0-_ee;
      }
      
    };

    /// The polytrope solver
    o2scl::root_bkt_cern<polytrope_solve> rbc;

    /// Array search object
    o2scl::search_vec<ubvector_range> sv_ub;

    /** \brief The number of grid points in integration of TOV equations
        for spherical stars
    */ 
    int RDIV;
  
    /** \brief Maximum value of s-coordinate (default 0.9999) */  
    double SMAX;
    /** \brief Spacing in \f$ s \f$ direction, 
        \f$ \mathrm{SMAX}/(\mathrm{SDIV}-1) \f$ 
    */
    double DS;
    /** \brief Spacing in \f$ \mu \f$ direction, \f$ 1/(\mathrm{MDIV}-1) \f$ 
     */ 
    double DM;

    /// Minimum radius for spherical stars (default \f$ 10^{-15} \f$)
    double RMIN;

    /// \name Grid quantities set in make_grid()
    //@{
    /** \brief The quantity \f$ s \f$, which runs from 0 to SMAX
        (which defaults to 0.9999) in steps of DS=SMAX/(SDIV-1)
    */
    ubvector s_gp;
    /// \f$ s (1-s) \f$
    ubvector s_1_s;
    /// \f$ 1-s \f$
    ubvector one_s;
    /** \brief The quantity \f$ \mu \f$ which runs from 0 to 1
        in steps of DM=1/(MDIV-1)
     */
    ubvector mu;
    /// \f$ 1-\mu^2 \f$
    ubvector one_m2;
    /// \f$ \theta \f$ defined by \f$ \mathrm{acos}~\mu \f$
    ubvector theta;
    /// \f$ \sin \theta = \sqrt{1-\mu^2}\f$
    ubvector sin_theta;
    //@}

    /// \name Grid values computed in integrate() for spherical_star()
    //@{
    /// Isotropic radius
    ubvector r_gp;
    /// Radial coordinate
    ubvector r_is_gp;
    /// Metric function \f$ \lambda \f$
    ubvector lambda_gp;
    /// Metric function \f$ \nu \f$
    ubvector nu_gp;
    /// Enclosed gravitational mass
    ubvector m_gp;
    /// Energy density
    ubvector e_d_gp;
    //@}

    /// Desc
    ubmatrix dgds;
    /// Desc
    ubmatrix dgdm;

    /// \name Metric functions
    //@{
    /** \brief potential \f$ \rho \f$ */ 
    ubmatrix rho;
    /** \brief potential \f$ \gamma \f$ */ 
    ubmatrix gamma;
    /** \brief potential \f$ \omega \f$ */ 
    ubmatrix omega;
    /** \brief potential \f$ \alpha \f$ */ 
    ubmatrix alpha;
    //@}

    /// \name Initial guess computed by the comp() function
    //@{
    /// Guess for the equatorial radius
    double r_e_guess;
    /** \brief Guess for \f$ \rho \f$ */ 
    ubmatrix rho_guess;
    /** \brief Guess for \f$ \gamma \f$ */
    ubmatrix gamma_guess;
    /** \brief Guess for \f$ \alpha \f$ */
    ubmatrix omega_guess;
    /** \brief Guess for \f$ \omega \f$ */
    ubmatrix alpha_guess;
    //@}

    /// \name EOS quantities
    //@{
    /** \brief Energy density \f$ \epsilon \f$ */
    ubmatrix energy;
    /** \brief Pressure */ 
    ubmatrix pressure;
    /** \brief Enthalpy */
    ubmatrix enthalpy;
    //@}

    /// \name Other quantities defined over the full two-dimensional grid
    //@{
    /** \brief Proper velocity squared */
    ubmatrix velocity_sq;
    /** \brief Derivative of \f$ \alpha \f$ with respect to \f$ \mu \f$ */
    ubmatrix da_dm;
    //@}

    /// \name Quantities defined for fixed values of mu
    //@{
    /** \brief \f$ \gamma(s) \f$ at \f$ \mu=1 \f$ */
    ubvector gamma_mu_1;
    /** \brief \f$ \gamma(s) \f$ at \f$ \mu=0 \f$ */
    ubvector gamma_mu_0;
    /** \brief \f$ \rho(s) \f$ at \f$ \mu=1 \f$ */
    ubvector rho_mu_1;
    /** \brief \f$ \rho(s) \f$ at \f$ \mu=0 \f$ */
    ubvector rho_mu_0;
    /** \brief \f$ \omega(s) \f$ at \f$ \mu=0 \f$ */
    ubvector omega_mu_0;
    //@}

    /** \brief The value of \f$ \hat{\gamma} \f$ at the pole */  
    double gamma_pole_h;                  
    /** \brief The value of \f$ \hat{\gamma} \f$ at the center */
    double gamma_center_h;                
    /** \brief The value of \f$ \hat{\gamma} \f$ at the equator */
    double gamma_equator_h;               
    /** \brief The value of \f$ \hat{\rho} \f$ at the pole */ 
    double rho_pole_h;                   
    /** \brief The value of \f$ \hat{\rho} \f$ at the center */
    double rho_center_h;                 
    /** \brief The value of \f$ \hat{\rho} \f$ at the equator */ 
    double rho_equator_h;                
    /** \brief The value of \f$ \hat{\omega} \f$ at the equator */
    double omega_equator_h;              
    /** \brief Angular velocity, \f$ \hat{\omega} \f$ */
    double Omega_h;                      
    /** \brief Central pressure */ 
    double p_center;                     
    /** \brief Central enthalpy */
    double h_center;                     

    /// \name Helper functions for Green's function expansion
    //@{
    /** \brief \f$ f_{\rho}(s,n,s') \f$ */
    tensor3<> f_rho;
    /** \brief \f$ f_{\gamma}(s,n,s') \f$ */
    tensor3<> f_gamma;
    /** \brief \f$ f_{\omega}(s,n,s') \f$ */
    tensor3<> f_omega;
    //@}
  
    /// \name Legendre polynomials
    //@{
    /** \brief Legendre polynomial \f$ P_{2n}(\mu) \f$ 
     */  
    ubmatrix P_2n;
    /** \brief Associated Legendre polynomial \f$ P^1_{2n-1}(\mu) \f$ 
     */ 
    ubmatrix P1_2n_1;
    //@}

    /** \brief Integrated term over m in eqn for \f$ \rho \f$ */
    ubmatrix D1_rho;
    /** \brief Integrated term over m in eqn for \f$ \gamma \f$ */
    ubmatrix D1_gamma;
    /** \brief Integ. term over m in eqn for \f$ \omega \f$ */
    ubmatrix D1_omega;
    /** \brief Integrated term over s in eqn for \f$ \rho \f$ */
    ubmatrix D2_rho;
    /** \brief Integrated term over s in eqn for \f$ \gamma \f$ */
    ubmatrix D2_gamma;
    /** \brief Integ. term over s in eqn for \f$ \omega \f$ */
    ubmatrix D2_omega;

    /** \brief source term in eqn for \f$ \gamma \f$ */
    ubmatrix S_gamma;
    /** \brief source term in eqn for \f$ \rho \f$ */
    ubmatrix S_rho;
    /** \brief source term in eqn for \f$ \omega \f$ */
    ubmatrix S_omega;

    /** \brief The tolerance for the functions with the prefix "fix" 
        (default \f$ 10^{-4} \f$ )
    */
    double tol_abs;

    /// \name Thermodyanmic quantities near the surface
    //@{
    /// Pressure at the surface
    double p_surface;
    /// Energy density at the surface
    double e_surface;
    /** \brief Minimum specific enthalpy
     */
    double enthalpy_min;                 
    //@}

    /// \name Polytrope parameters
    //@{
    /// Polytropic index
    double n_P;
    /// Polytropic exponent
    double Gamma_P;
    //@}

    /// \name Interpolation functions
    //@{
    /** \brief Cache for interpolation
     */
    int n_nearest;
  
    /// Search in array \c x of length \c n for value \c val
    int new_search(int n, ubvector &x, double val);

    /** \brief Driver for the interpolation routine. 
        
        First we find the tab. point nearest to xb, then we
        interpolate using four points around xb.
        
        Used by \ref int_z(), \ref e_at_p(), \ref p_at_e(), \ref
        p_at_h(), \ref h_at_p(), \ref n0_at_e(), \ref comp_omega(),
        \ref comp_M_J(), \ref comp(), \ref spherical_star(), \ref
        iterate().
    */  
    double interp(ubvector &xp, ubvector &yp, int np, double xb);

    /** \brief Driver for the interpolation routine.

        Four point interpolation at a given offset the index of the
        first point k.

        Used in \ref comp() .
    */
    double interp_4_k(ubvector &xp, ubvector &yp, int np, double xb, int k);
    //@}

    /** \brief Integrate f[mu] from m-1 to m. 

        This implements a 8-point closed Newton-Cotes formula.
        
        Used in \ref comp() .
    */
    double int_z(ubvector &f, int m);

    /// \name EOS functions
    //@{
    /** \brief Compute \f$ \varepsilon(P) \f$  
        
        Used in \ref dm_dr_is(), \ref dp_dr_is(), \ref integrate()
        and \ref iterate(). 
    */
    double e_at_p(double pp);

    /** \brief Compute \f$ P(\varepsilon) \f$  
        
        Used in \ref make_center() and \ref integrate().
    */
    double p_at_e(double ee);

    /** \brief Pressure at fixed enthalpy

        Used in \ref iterate().
    */
    double p_at_h(double hh);

    /** \brief Enthalpy at fixed pressure 

        Used in \ref make_center() and \ref integrate().
    */
    double h_at_p(double pp);
    
    /** \brief Baryon density at fixed energy density 

        Used in \ref comp_M_J() and \ref comp() .
    */
    double n0_at_e(double ee);
    //@}

    /// \name Derivatives on the grid
    //@{
    /** \brief Returns the derivative w.r.t. s of an array f[SDIV+1]. 
     */ 
    double s_deriv(ubvector &f, int s);

    /** \brief Returns the derivative w.r.t. mu of an array f[MDIV+1]. 
     */ 
    double m_deriv(ubvector &f, int m);

    /** \brief Returns the derivative w.r.t. s  
     */
    double deriv_s(ubmatrix &f, int s, int m);

    /** \brief Returns the derivative w.r.t. mu 
     */ 
    double deriv_m(ubmatrix &f, int s, int m);

    /** \brief Returns the derivative w.r.t. s and mu 
     */ 
    double deriv_sm(ubmatrix &f, int s, int m);
    //@}

    /// \name Initialization functions
    //@{
    /** \brief Returns the Legendre polynomial of degree n, evaluated at x. 

        This uses the recurrence relation and is used in \ref comp_f_P()
        which is called by the constructor.
    */
    double legendre(int n, double x);

    /** \brief Compute two-point functions
        
        This function computes the 2-point functions \f$
        f^m_{2n}(r,r') \f$ used to integrate the potentials \f$ \rho,
        \gamma \f$ and \f$ \omega \f$ (See Komatsu et al. 1989 for
        details). Since the grid points are fixed, we can compute the
        functions \ref f_rho, \ref f_gamma, \ref f_omega, \ref P_2n,
        and \ref P1_2n_1 once at the beginning.

        \verbatim embed:rst
        See Eqs. 27-29 of [Cook92]_ and Eqs. 33-35 of 
        [Komatsu89]_. This function is called by the constructor.
        \endverbatim
    */
    void comp_f_P();

    /** \brief Create computational mesh. 

        Create the computational mesh for \f$ s=r/(r+r_e) \f$
        (where \f$ r_e \f$ is the coordinate equatorial radius) 
        and \f$ \mu = \cos \theta \f$
        using 
        \f[
        s[i]=\mathrm{SMAX}\left(\frac{i-1}{\mathrm{SDIV}-1}\right)
        \f]
        \f[
        \mu[j]=\left(\frac{i-1}{\mathrm{MDIV}-1}\right)
        \f]
        When \f$ r=0 \f$, \f$ s=0 \f$, when \f$ r=r_e \f$, 
        \f$ s=1/2 \f$, and when \f$ r = \infty \f$, \f$ s=1 \f$ .
        Inverting the relationship between \f$ r \f$ and \f$ s \f$
        gives \f$ r = r_e s / (1-s) \f$ .
        \comment
        (Note that some versions of the manual have a typo,
        giving \f$ 1-i \f$ rather than \f$ i-1 \f$ above.)
        \endcomment
        
        Points in the mu-direction are stored in the array
        <tt>mu[i]</tt>. Points in the s-direction are stored in the
        array <tt>s_gp[j]</tt>.

        This function sets \ref s_gp, \ref s_1_s, \ref one_s,
        \ref mu, \ref one_m2, \ref theta and \ref sin_theta .
        All of these arrays are unit-indexed. It is called by
        the constructor.
    */
    void make_grid();
    //@}

    /** \brief Compute central pressure and enthalpy from central
        energy density

        For polytropic EOSs, this also computes <tt>rho0_center</tt> .
    */
    void make_center(double e_center);

    /// \name Post-processing functions
    //@{
    /** \brief Compute Omega and Omega_K. 
     */
    void comp_omega();
  
    /** \brief Compute rest mass and angular momentum. 
     */
    void comp_M_J();

    /** \brief Compute various quantities.

        The main post-processing function
    */
    void comp();
    //@}

    /// \name For computing spherical stars
    //@{
    /** \brief Computes a spherically symmetric star 
        
        The metric is 
        \f[
        ds^2 = -e^{2\nu}dt^2 + e^{2\lambda} dr^2 + r^2 d\theta^2 + 
        r^2 sin^2\theta d\phi^2
        \f]
        where \f$ r \f$ is an isotropic radial coordinate 
        (corresponding to <tt>r_is</tt> in the code).
      
        \todo AWS: 7/20/20: Better document out how this metric
        definition leads to \f$ \gamma=\nu+\lambda \f$ and \f$ \rho =
        \nu - \lambda \f$ and the relationship between r and r_is .

        This function computes \ref r_e_guess, \ref R_e, 
        \ref Mass, and \ref Z_p .
    */
    void spherical_star();

    /** \brief Derivative of gravitational mass with respect to
        isotropic radius */
    double dm_dr_is(double r_is, double r, double m, double p);
 
    /** \brief Derivative of pressure with respect to isotropic radius */
    double dp_dr_is(double r_is, double r, double m, double p);

    /** \brief Derivative of radius with respect to isotropic radius */
    double dr_dr_is(double r_is, double r, double m);
  
    /** \brief Integrate one of the differential equations for 
        spherical stars*/
    void integrate(int i_check, double &r_final, double &m_final,
                   double &r_is_final);
    //@}

    /// \name Desc
    //@{
    /** \brief Main iteration function
     */
    int iterate(double r_ratio, double tol_rel);
    //@}

    /// \name EOS member variables
    //@{ 
    /** \brief If true, then an EOS has been set
     */
    bool eos_set;
  
    /** \brief If true, then use a polytrope and rescale
     */
    bool scaled_polytrope;

    /** \brief Pointer to the user-specified EOS
     */
    eos_nstar_rot *eosp;
    //@}

    /** \brief Compute masses and angular momentum
     */
    void calc_masses_J(ubmatrix &rho_0);
    
  public:

    nstar_rot();

    /** \brief Relative accuracy for the equatorial radius,
        \f$ r_e \f$ (default \f$ 10^{-5} \f$) 

        Used in \ref iterate() .
    */
    double eq_radius_tol_rel;                    

    /** \brief Accuracy for equatorial radius using alternate
        solvers (default \f$ 10^{-9} \f$) 
     */
    double alt_tol_rel;                    

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Create an output table3d object from results

        This function creates slices named
        "ed pr h vsq rho gamma omega alpha" on the grid of
        \f$ s \f$ and \f$ \mu \f$.
     */
    void output_table(o2scl::table3d &t);
    
    /// \name Output
    //@{
    /** \brief Central energy density (in units of 
        \f$ 10^{15} \mathrm{g}/\mathrm{cm}^3 \f$) 
    */
    double e_center;                     
    /** \brief Ratio of polar to equatorial radius
     */ 
    double r_ratio;                      
    /** \brief Coordinate equatorial radius
     */ 
    double r_e;                          
    //@}

    /// \name Quantities computed by nstar_rot::comp() (in order)
    //@{
    /** \brief Radius at pole */      
    double r_p;                          
    /** \brief The value of the s-coordinate at the pole */
    double s_p;                          
    /** \brief The value of the s-coordinate at the equator */
    double s_e; 
    /// The velocity at the equator
    double velocity_equator;              
    /** \brief Circumferential radius in cm (i.e. the radius defined
        such that \f$ 2 \pi R_e \f$ is the proper circumference) */
    double R_e;                          
    /// Proper mass (in \f$ \mathrm{g} \f$ )
    double Mass_p;
    /// Gravitational mass (in \f$ \mathrm{g} \f$ )
    double Mass;
    /// Baryonic mass (in \f$ \mathrm{g} \f$ )
    double Mass_0;
    /// Angular momentum (in \f$ \mathrm{g}~\mathrm{cm}^2 / \mathrm{s} \f$ )
    double J;
    /// Angular velocity (in \f$ \mathrm{radians}/\mathrm{s} \f$ )
    double Omega;
    /// Total rotational kinetic energy
    double T;
    /// Moment of inertia
    double I;
    /// Gravitational binding energy
    double W;
    /// Polar redshift
    double Z_p;
    /// Forward equatorial redshift
    double Z_f;
    /// Backward equatorial redshift
    double Z_b;
    /** \brief Kepler rotation frequency (in 1/s) */  
    double Omega_K;                      
    /// The eccentricity
    double eccentricity;
    /// Desc
    ubvector v_plus;
    /// Desc
    ubvector v_minus;
    /// Desc
    double vel_plus;
    /// Desc
    double vel_minus;
    /** \brief Height from surface of last stable co-rotating circular 
        orbit in equatorial plane

        If this is zero then all orbits are stable.
    */
    double h_plus;
    /** \brief Height from surface of last stable counter-rotating circular 
        orbit in equatorial plane
        
        If this is zero then all orbits are stable.
    */
    double h_minus;
    /// Desc
    double Omega_plus;
    /// Desc
    double u_phi;
    /// Angular velocity of a particle in a circular orbit at the equator
    double Omega_p;
    /// Desc
    double grv2;
    /// Desc
    double grv2_new;
    /// Desc
    double grv3;
    /** \brief Ratio of potential \f$ \omega \f$ to angular 
        velocity \f$ \Omega \f$
    */
    double om_over_Om;
    /** \brief Mass quadrupole moment
     */
    double mass_quadrupole;
    //@}

    /// \name Settings
    //@{
    /// The convergence factor (default 1.0)
    double cf;
    //@}

    /// \name Internal constants
    //@{
    /** \brief Use the values of the constants from the original RNS
        code
    */
    void constants_rns();
    /** \brief Use the O2scl constants (the default)
     */
    void constants_o2scl();
    /** \brief Speed of light in vacuum (in CGS units) */ 
    double C;
    /** \brief Gravitational constant (in CGS units) */ 
    double G;
    /** \brief Mass of sun (in g) */
    double MSUN;
    /** \brief Square of length scale in CGS units, 
        \f$ \kappa \equiv 10^{-15} c^2/G \f$
    */
    double KAPPA;
    /** \brief The mass of one baryon (in g)
     */
    double MB;
    /** \brief The value \f$ \kappa G c^{-4} \f$ */
    double KSCALE;
    /// The constant \f$ \pi \f$
    double PI;
    //@}

    /// \name Basic Usage
    //@{
    /** \brief Set the EOS
     */
    void set_eos(eos_nstar_rot &eos) {
      eosp=&eos;
      eos_set=true;
      scaled_polytrope=false;
      return;
    }

    /** \brief Use a polytropic EOS with a specified index
     */
    void polytrope_eos(double index) {
      n_P=index;
      scaled_polytrope=true;
      eos_set=true;
      return;
    }
    
    /** \brief Construct a configuration with a fixed central 
        energy density and a fixed axis ratio

        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ and the axis ratio is unitless.
        This is fastest of the high-level interface functions as it
        doesn't require an additional solver.
    */
    int fix_cent_eden_axis_rat(double cent_eden, double axis_rat,
                               bool use_guess=false);
    
    /** \brief Construct a configuration with a fixed central 
        energy density and a fixed gravitational mass
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ and the gravitational 
        mass should be in solar masses. 
    */
    int fix_cent_eden_grav_mass(double cent_eden, double grav_mass);

    /** \brief Construct a configuration with a fixed central 
        energy density and a fixed baryonic mass
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ and the baryonic 
        mass should be in solar masses. 
    */
    int fix_cent_eden_bar_mass(double cent_eden, double bar_mass);

    /** \brief Construct a configuration with a fixed central 
        energy density and the Keplerian rotation rate
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ .
    */
    int fix_cent_eden_with_kepler(double cent_eden);

    /** \brief Experimental alternate form for
        \ref fix_cent_eden_with_kepler()
    */
    int fix_cent_eden_with_kepler_alt(double cent_eden,
                                      bool use_guess=false);

    /** \brief Experimental alternate form for
        \ref fix_cent_eden_grav_mass()
    */
    int fix_cent_eden_grav_mass_alt(double cent_eden, double grav_mass,
                                    bool use_guess=false);
    
    /** \brief Experimental alternate form for
        \ref fix_cent_eden_bar_mass()
    */
    int fix_cent_eden_bar_mass_alt(double cent_eden, double bar_mass,
                                   bool use_guess=false);
    
    /** \brief Experimental alternate form for
        \ref fix_cent_eden_ang_vel()
    */
    int fix_cent_eden_ang_vel_alt(double cent_eden, double ang_vel,
                                  bool use_guess=false);
    
    /** \brief Experimental alternate form for
        \ref fix_cent_eden_ang_mom()
    */
    int fix_cent_eden_ang_mom_alt(double cent_eden, double ang_mom,
                                  bool use_guess=false);
    
    /** \brief Construct a non-rotating configuration with a fixed central 
        energy density
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ .
    */
    int fix_cent_eden_non_rot(double cent_eden);

    /** \brief Construct a configuration with a fixed central 
        energy density and a fixed angular velocity.
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$ and the angular
        velocity should be in \f$ \mathrm{rad}/\mathrm{s} \f$.
        The final angular velocity (possibly slightly different
        than <tt>ang_vel</tt> is stored in \ref Omega .
        
        \note In the original RNS code, the <tt>ang_vel</tt> argument
        is different because it was rescaled by a factor of \f$ 10^{4}
        \f$.
    */
    int fix_cent_eden_ang_vel(double cent_eden, double ang_vel);

    /** \brief Construct a configuration with a fixed central 
        energy density and a fixed angular momentum.
        
        The central energy density should be in \f$
        \mathrm{g}/\mathrm{cm}^3 \f$. The angular momentum should be
        in units of \f$ G M_{\odot}^2/C \f$ .
    */
    int fix_cent_eden_ang_mom(double cent_eden, double ang_mom);
    //@}
    
    /** \name Testing functions

        These functions compare this class with hard-coded results
        obtained with the RNS code.
    */
    //@{
    /** \brief Test determining configuration with fixed central
        energy density and fixed radius ratio with EOS C
    */    
    void test1(o2scl::test_mgr &t);
    
    /** \brief Test configuration rotating and Keplerian frequency
        with a fixed central energy density and EOS C
    */    
    void test2(o2scl::test_mgr &t);
    
    /** \brief Test fixed central energy density and fixed 
        gravitational mass with EOS C
    */    
    void test3(o2scl::test_mgr &t);
    
    /** \brief Test fixed central energy density and fixed baryonic 
        mass with EOS C
    */    
    void test4(o2scl::test_mgr &t);
    
    /** \brief Test fixed central energy density and fixed angular
        velocity with EOS C
    */    
    void test5(o2scl::test_mgr &t);
    
    /** \brief Test fixed central energy density and fixed angular 
        momentum with EOS C
    */    
    void test6(o2scl::test_mgr &t);

    /** \brief Test a series of non-rotating stars on a energy density
        grid with EOS C
    */    
    void test7(o2scl::test_mgr &t);
    
    /** \brief Test Keplerian frequency for a polytrope
     */    
    void test8(o2scl::test_mgr &t);
    //@}

  };

}

#endif
