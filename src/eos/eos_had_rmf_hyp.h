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
    virtual int beta_eq_T0(ubvector &nB_grid, ubvector &guess,
			   fermion &e, bool include_muons,
			   fermion &mu, fermion_rel &frel,
			   std::shared_ptr<table_units<> > results) {
      
      if (guess[0]<=0.0 || guess[0]>=nB) guess[0]=nB/2.0;
      if (guess[1]!=0.0 && guess[2]!=0.0 && guess[3]!=0.0) {
	this->set_fields(guess[1],guess[2],guess[3]);
      }

      double nB_temp;
      
      mm_funct fmf=std::bind
	(std::mem_fn<int(size_t,const ubvector &, ubvector &, 
			 const double &, fermion &, bool,
			 fermion &, fermion_rel &)>
	 (&eos_had_eden_base::solve_beta_eq_T0),
	 this,std::placeholders::_1,std::placeholders::_2,
	 std::placeholders::_3,std::cref(nB_temp),std::ref(e),
	 include_muons,std::ref(mu),std::ref(frel));
      
      results->clear();
      results->line_of_names(((std::string)"ed pr nb nn np nlam ")+
			     "nsigp nsigz nsigm mun mup mulam musigp musigz "+
			     "musigm kfn kfp kflam kfsigp kfsigz kfsigm");
      results->line_of_units(((std::string)"1/fm^4 1/fm^4 1/fm^3 ")+
			     "1/fm^3 1/fm^3 1/fm^3 1/fm^3 1/fm^3 "+
			     "1/fm^3 1/fm^3 1/fm 1/fm 1/fm 1/fm 1/fm 1/fm "+
			     "1/fm 1/fm 1/fm 1/fm 1/fm 1/fm");
      if (inc_cascade) {
	results->line_of_names("ncasz ncasm mucasz mucasm kfcasz kfcasm");
	results->set_unit("ncasz","1/fm^3");
	results->set_unit("ncasz","1/fm^3");
	results->set_unit("mucasz","1/fm");
	results->set_unit("mucasm","1/fm");
	results->set_unit("kfcasz","1/fm");
	results->set_unit("kfcasm","1/fm");
      }
      
      for(size_t i=0;i<nB_grid.size();i++) {
	nB_temp=nB_grid[i];
	
	beta_mroot.solve(1,guess,fmf);
	
	// Final function evaluation to make sure, e.g.
	// eos_thermo object is correct
	ubvector y(1);
	fmf(1,guess,y);

	std::vector<double> line={eos_thermo->ed,eos_thermo->pr,nB_temp,
				  neutron->n,proton->n,lambda->n,
				  sigma_p->n,sigma_z->n,sigma_m->n,
				  neutron->mu,proton->mu,lambda->mu,
				  sigma_p->mu,sigma_z->mu,sigma_m->mu,
				  neutron->kf,proton->kf,lambda->kf,
				  sigma_p->kf,sigma_z->kf,sigma_m->kf};
	results->line_of_data(line);
	if (inc_cascade) {
	  row=results->get_nlines()-1;
	  results->set("ncasz",row,cascade_z->n);
	  results->set("ncasm",row,cascade_z->n);
	  results->set("mucasz",row,cascade_z->mu);
	  results->set("mucasm",row,cascade_z->mu;
	  results->set("kfcasz",row,cascade_z->kf);
	  results->set("kfcasm",row,cascade_z->kf);
	}
	  
      }
      
      return 0;
    }

#endif

  };

#ifndef DOXYGENP
}
#endif

#endif
