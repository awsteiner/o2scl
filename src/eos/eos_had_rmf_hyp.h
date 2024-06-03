/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

#ifndef O2SCL_EOS_HAD_RMF_H
#define O2SCL_EOS_HAD_RMF_H

#include <string>
#include <cmath>
#include <o2scl/lib_settings.h>
#include <o2scl/constants.h>
#include <o2scl/mm_funct.h>
#include <o2scl/part.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/fermion.h>

namespace o2scl {
  
  /** \brief Relativistic mean field theory EOS with hyperons
      
      \verbatim embed:rst

      Based on [Glendenning91ro]_, but generalized for higher-order
      couplings as in \ref eos_had_rmf .
      
      .. todo:: 

         In class eos_had_rmf_hyp:

         - The couplings in the test code match the table but the
           maximum masses appear much smaller than GM91. I need to
           check that muons are added correctly, and it might be good 
           to compare with a different reference. This also might be 
           due to a different crust EOS. 
         - The interpretation of the calc_e() function is a bit
           unclear, so I need to more clearly figure out what 
           that function ought to do. I don't think it's really used
           at the moment. 

      \endverbatim
  */
  class eos_had_rmf_hyp : public eos_had_rmf {
    
  protected:

    /// The neutron object
    fermion *lambda;

    /// The proton object
    fermion *sigma_p;

    /// The neutron object
    fermion *sigma_z;

    /// The proton object
    fermion *sigma_m;

    /// The neutron object
    fermion *cascade_z;

    /// The proton object
    fermion *cascade_m;

    /// The function for calc_e()
    virtual int calc_e_solve_fun(size_t nv, const ubvector &ex, 
                         ubvector &ey);

    /// The function for calc_e_hyp_nobeta()
    int calc_e_nobeta_fun(size_t nv, const ubvector &ex, 
                          ubvector &ey, double nB,
                          double Ye, double Ys);
    
    /// The function for calc_e_hyp_nobeta_np()
    int calc_e_nobeta_np_fun(size_t nv, const ubvector &ex, 
                             ubvector &ey, double nn,
                             double np, double Ys);
    
    /** \brief Equation for solving for beta-equilibrium at T=0
    */
    virtual int solve_beta_eq_T0(size_t nv, const ubvector &x,
                                 ubvector &y, const double &nB,
                                 eos_leptons &elep);

  public:

    eos_had_rmf_hyp();

    /// \name Hyperon objects
    //@{
    /// The default Lambda hyperon
    fermion def_lambda;

    /// The default Sigma plus hyperon
    fermion def_sigma_p;

    /// The default Sigma zero hyperon
    fermion def_sigma_z;

    /// The default Sigma minus hyperon
    fermion def_sigma_m;

    /// The default Xi zero hyperon
    fermion def_cascade_z;

    /// The default Xi minus hyperon
    fermion def_cascade_m;
    //@}

    /// \name Hyperon-meson couplings
    //@{
    double xs;
    double xw;
    double xr;
    //@}

    /// If true, include cascade hyperons (default true)
    bool inc_cascade;

    /** \brief Equation of state and meson field equations 
        as a function of chemical potentials
    */
    virtual int calc_eq_hyp_p
    (fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
     fermion &sigm, fermion &casz, fermion &casm, double sig, double ome, 
     double lrho, double &f1, double &f2, double &f3, thermo &lth);
    
    /** \brief Equation of state and meson field equations 
        as a function of chemical potentials at finite temperature
    */
    virtual int calc_eq_temp_hyp_p
    (fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
     fermion &sigm, fermion &casz, fermion &casm, double T,
     double sig, double ome, 
     double lrho, double &f1, double &f2, double &f3, thermo &lth);
    
    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
        binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
    */
    void calc_xs(double lam_be);

    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
        binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
     */
    void calc_xw(double lam_be);

    /** \brief Equation of state as a function of density at 
        fixed baryon and charge density presuming the hyperons 
        are in beta-equilibrium with the nucleons

        Initial guesses for the chemical potentials are taken
        from the user-given values. Initial guesses for the fields
        can be set by set_fields(), or default values will be used.
        After the call to calc_e(), the final values of the fields
        can be accessed through get_fields(). 
    */
    virtual int calc_hyp_e(double n_baryon, double n_charge,
                           fermion &ne, fermion &pr,
                           fermion &lam, fermion &sigp, fermion &sigz, 
                           fermion &sigm, fermion &casz, fermion &casm,
                           thermo &lth);
    
    /** \brief Equation of state as a function of density at 
        fixed baryon and charge density presuming the hyperons 
        are in beta-equilibrium with the nucleons

        Initial guesses for the chemical potentials are taken
        from the user-given values. Initial guesses for the fields
        can be set by set_fields(), or default values will be used.
        After the call to calc_e(), the final values of the fields
        can be accessed through get_fields(). 
    */
    virtual int calc_temp_hyp_e(double n_baryon, double n_charge,
                                fermion &ne, fermion &pr,
                                fermion &lam, fermion &sigp, fermion &sigz, 
                                fermion &sigm, fermion &casz, fermion &casm,
                                double T, thermo &lth);

    /** \brief Equation of state as a function of density
        out of beta equilibrium
    */
    virtual int calc_hyp_e_nobeta(double nB, double Ye, double Ys,
                                  fermion &ne, fermion &pr,
                                  fermion &lam, fermion &sigp, fermion &sigz, 
                                  fermion &sigm, fermion &casz, fermion &casm,
                                  thermo &lth);
    
    /** \brief Equation of state as a function of density
        out of beta equilibrium at fixed neutron and proton densities
    */
    virtual int calc_hyp_e_nobeta_np(double Ys, fermion &ne, fermion &pr,
                                     fermion &lam, fermion &sigp,
                                     fermion &sigz, fermion &sigm,
                                     fermion &casz, fermion &casm,
                                     thermo &lth);
    
    /** \brief Set the hyperon objects
     */
    virtual void set_hyp(fermion &lam, fermion &sigp, fermion &sigz, 
                         fermion &sigm, fermion &casz, fermion &casm);
    
    /** \brief Compute the EOS in beta-equilibrium at 
        zero temperature
    */
    virtual int beta_eq_T0(ubvector &nB_grid, ubvector &guess,
                           eos_leptons &elep,
                           std::shared_ptr<table_units<> > results);

#ifdef O2SCL_NEVER_DEFINED
    
    /** \brief Function to solve for calc_temp_f_gen()
     */
    virtual int solve_ctf(size_t nv, const ubvector &x, ubvector &y,
                          double T, thermo &th, double nB, double nQ,
                          double nS) {
      def_lambda.n=x[0];
      def_sigma_p.n=x[1];
      def_sigma_z.n=x[2];
      def_cascade_z.n=x[3];
      def_cascade_m.n=x[4];

      def_neutron.n=nB+2.0*def_cascade_m.n+3.0*def_cascade_z.n-
        def_lambda.n-nQ-2.0*nS+2.0*def_sigma_p.n+def_sigma_z.n;
      def_proton.n=-def_cascade_m.n-2.0*def_cascade_z.n+nQ+nS-
        2.0*def_sigma_p.n-def_sigma_z.n;
      def_sigma_m.n=-2.0*def_cascade_m.n-2.0*def_cascade_z.n-nS-
        def_sigma_p.n-def_sigma_z.n;

      int ret=calc_temp_hyp_e(def_neutron,def_proton,def_lambda,
                              def_sigma_p,def_sigma_z,def_sigma_m,
                              def_cascade_z,def_cascade_m,T,th);

      double muB=def_neutron.mu;
      double muQ=def_proton.mu-def_neutron.mu;
      double muS=def_lambda.mu-def_neutron.mu;
      
      y[0]=def_sigma_p.mu-muB-muQ-muS;
      y[1]=def_sigma_z.mu-muB-muS;
      y[2]=def_sigma_m.mu-muB+muQ-muS;
      y[3]=def_cascade_z.mu-muB-2.0*muS;
      y[4]=def_cascade_m.mu-muB+muQ-2.0*muS;
      
      return ret;
    }
    
    /** \brief Equation of state as a function of baryon, charge,
        and strangeness density at finite temperature
    */
    virtual int calc_temp_f_gen(double nB, double nQ, double nS,
                                double T, thermo &th) {
      
      mm_funct fmf=std::bind
        (std::mem_fn<int(size_t,const ubvector &,ubvector &,
                         double,thermo &,double,double,double)>
         (&eos_had_rmf_hyp::solve_ctf),
         this,std::placeholders::_1,std::placeholders::_2,
         std::placeholders::_3,T,std::ref(th),nB,nQ,nS);

      ubvector x(5), y(5);
      x[0]=def_lambda.n;
      x[1]=def_sigma_p.n;
      x[2]=def_sigma_z.n;
      x[3]=def_cascade_z.n;
      x[4]=def_cascade_m.n;

      int ret=eos_mroot->msolve(5,x,fmf);

      return ret;
    }
    
#endif    
    
  };

}

#endif
