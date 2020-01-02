/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

#ifndef DOXYGENP
namespace o2scl {
#endif
  
  /** \brief Relativistic mean field theory EOS with hyperons
      
      Based on \ref Glendenning91ro, but generalized for higher-order
      couplings as in \ref eos_had_rmf .
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
    virtual int calc_eq_p
      (fermion &ne, fermion &pr, fermion &lam, fermion &sigp, fermion &sigz, 
       fermion &sigm, fermion &casz, fermion &casm, double sig, double ome, 
       double lrho, double &f1, double &f2, double &f3, thermo &lth);

    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
	binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
     */
    void calc_xs(double lam_be);

    /** \brief Compute xs assuming a fixed value of the \f$ \Lambda \f$
	binding energy in nuclear matter in \f$ \mathrm{fm}^{-1} \f$
     */
    void calc_xw(double lam_be);

    /** \brief Equation of state as a function of density

	Initial guesses for the chemical potentials are taken
	from the user-given values. Initial guesses for the fields
	can be set by set_fields(), or default values will be used.
	After the call to calc_e(), the final values of the fields
	can be accessed through get_fields(). 
    */
    virtual int calc_e(fermion &ne, fermion &pr,
		       fermion &lam, fermion &sigp, fermion &sigz, 
		       fermion &sigm, fermion &casz, fermion &casm,
		       thermo &lth);

#ifdef O2SCL_NEVER_DEFINED

    /** \brief Set the hyperon objects
     */
    virtual void set_hyp(fermion &lam, fermion &sigp, fermion &sigz, 
			 fermion &sigm, fermion &casz, fermion &casm) {
      lambda=&lam;
      sigma_p=&sigp;
      sigma_z=&sigz;
      sigma_m=&sigm;
      cascade_z=&casz;
      cascade_m=&casm;
      return 0;
    }
    
    /** \brief Compute the EOS in beta-equilibrium at 
	zero temperature
    */
    virtual int beta_eq_T0(double nB, std::vector<double> &guess,
			   fermion &e, bool include_muons,
			   fermion &mu, fermion_rel &frel,
			   std::vector<double> &res) {
      
      if (guess[0]<=0.0 || guess[0]>=nB) guess[0]=nB/2.0;
      this->set_fields(guess[1],guess[2],guess[3]);
      
      mm_funct fmf=std::bind
	(std::mem_fn<int(size_t,const ubvector &, ubvector &, 
			 const double, fermion &, bool,
			 fermion &, fermion_rel &)>
	 (&eos_had_eden_base::solve_beta_eq_T0),
	 this,std::placeholders::_1,std::placeholders::_2,
	 std::placeholders::_3,nB,std::ref(e),bool,std::ref(mu),
	 std::ref(frel));
      
      beta_mroot.solve(1,guess,fmf);

      // Final function evaluation to make sure, e.g.
      // eos_thermo object is correct
      ubvector y(1);
      fmf(1,guess,y);

      if (inc_cascade) {
	if (res.size()<26) res.resize(26);
      } else {
	if (res.size()<20) res.resize(20);
      }
      
      res[0]=neutron->n;
      res[1]=proton->n;
      res[2]=lambda->n;
      res[3]=sigma_p->n;
      res[4]=sigma_z->n;
      res[5]=sigma_m->n;
      
      res[6]=eos_thermo->ed;
      res[7]=eos_thermo->pr;
      
      res[8]=neutron->mu;
      res[9]=proton->mu;
      res[10]=lambda->mu;
      res[11]=sigma_p->mu;
      res[12]=sigma_z->mu;
      res[13]=sigma_m->mu;
      
      res[14]=neutron->kf;
      res[15]=proton->kf;
      res[16]=lambda->kf;
      res[17]=sigma_p->kf;
      res[18]=sigma_z->kf;
      res[19]=sigma_m->kf;
      
      if (inc_cascade) {
	res[20]=cascade_z->n;
	res[21]=cascade_m->n;
	res[22]=cascade_z->mu;
	res[23]=cascade_m->mu;
	res[24]=cascade_z->kf;
	res[25]=cascade_m->kf;
      }
      
      return 0;
    }

#endif

  };

#ifndef DOXYGENP
}
#endif

#endif
