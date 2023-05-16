/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
/** \file eos_had_tabulated.h
    \brief File defining \ref o2scl::eos_had_tabulated
*/
#ifndef O2SCL_TABULATED_EOS_H
#define O2SCL_TABULATED_EOS_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>
#include <o2scl/eos_had_apr.h>

namespace o2scl {

  /** \brief Zero-temperature EOS from a table

      This assumes a symmetry energy which depends quadratically on
      the isospin asymmetry in order to construct an EOS from 
      a table of baryon density and energy per baryon for both
      nuclear and pure neutron matter. 

      Note: If using a tabulated EOS to compute derivatives (like the
      compressibility which effectively requires a second derivative),
      it is important to tabulated the EOS precisely enough to ensure
      that the derivatives are accurate. In the case of ensuring that
      the compressibility at saturation density is reproduced to
      within 3 significant figures, the EOS should be specified with
      at least 6 digits of precision on a grid at least as small as
      0.002 \f$ \mathrm{fm}^{-3} \f$.
      
      \future Storage in a \ref o2scl::table object isn't necessary, and
      storing in a vector object is more efficient.
      \future Allow the pressure to be specified to make the
      EOS more accurate?
  */
  class eos_had_tabulated : public eos_had_eden_base {
    
  protected:

    /// True if the table has been allocated
    bool table_alloc;

    /// \name The EOS tables
    //@{
    table<> *tnuc;
    table<> *tneut;
    //@}

    /// If true, then tnuc and tneut point to the same table
    bool one_table;

    /// \name Strings for the column names
    //@{
    std::string srho_nuc, srho_neut, snuc, sneut;
    //@}

    /// Free the table memory
    void free_table() {
      if (table_alloc) {
	delete tnuc;
	if (!one_table) delete tneut;
	table_alloc=false;
      }
      return;
    }
    
  public:

    eos_had_tabulated() {
      table_alloc=false;
      one_table=false;
    }

    virtual ~eos_had_tabulated() {
      if (table_alloc) {
	delete tnuc;
	if (!one_table) delete tneut;
      }
    }

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(fermion &ne, fermion &pr, thermo &th) {
      
      if (table_alloc==false) {
	O2SCL_ERR("No EOS specified in eos_had_tabulated::calc_e().",
		    exc_einval);
      }
      double barn=ne.n+pr.n;
      double xp=pr.n/barn;
      double delta=(1.0-2.0*xp);

      // The energy density of nuclear matter
      double ednuc=(tnuc->interp(srho_nuc,barn,snuc)/o2scl_const::hc_mev_fm+
		    ne.m)*barn;
      // The symmetry energy density
      double edsym=(tneut->interp(srho_neut,barn,sneut)-
		    tnuc->interp(srho_nuc,barn,snuc))/
	o2scl_const::hc_mev_fm*barn;
      // The total energy density
      th.ed=ednuc+delta*delta*edsym;
      
      // The derivatives of the energy densities wrt the baryon density
      double dednucdn=tnuc->deriv(srho_nuc,barn,snuc)/
	o2scl_const::hc_mev_fm*barn+ednuc/barn;
      double dedsymdn=barn*(tneut->deriv(srho_neut,barn,sneut)-
			    tnuc->deriv(srho_nuc,barn,snuc))/
	o2scl_const::hc_mev_fm+edsym/barn;
      
      // The chemical potentials
      ne.mu=(dednucdn+delta*delta*dedsymdn)+4.0*delta*edsym*xp/barn;
      pr.mu=(dednucdn+delta*delta*dedsymdn)+4.0*delta*edsym*(xp-1.0)/barn;

      // The pressure
      th.pr=-th.ed+ne.n*ne.mu+pr.n*pr.mu;

      return 0;
    }
    
    /** \brief Set the EOS through vectors specifying the densities and 
	energies
    */
    template <class vec_t>
      int set_eos(size_t n, vec_t &rho, vec_t &Enuc, vec_t &Eneut) {

      free_table();

      tnuc=new table<>(n);
      tnuc->line_of_names("rho nuc neut");
      srho_nuc="rho";
      srho_neut="rho";
      snuc="nuc";
      sneut="neut";
      for(size_t i=0;i<n;i++) {
	double line[3]={rho[i],Enuc[i],Eneut[i]};
	tnuc->line_of_data(3,line);
      }
      tneut=tnuc;
      table_alloc=true;
      one_table=true;

      return 0;
    }

    /** \brief Set the EOS through vectors specifying the densities and 
	energies
    */
    template<class vec_t> 
      int set_eos(size_t n_nuc, vec_t &rho_nuc, vec_t &E_nuc, 
		  size_t n_neut, vec_t &rho_neut, vec_t &E_neut) {
      
      free_table();

      tnuc=new table<>(n_nuc);
      tneut=new table<>(n_neut);
      tnuc->line_of_names("rho nuc");
      tneut->line_of_names("rho neut");
      srho_nuc="rho";
      srho_neut="rho";
      snuc="nuc";
      sneut="neut";
      for(size_t i=0;i<n_nuc;i++) {
	double line[2]={rho_nuc[i],E_nuc[i]};
	tnuc->line_of_data(2,line);
      }
      for(size_t i=0;i<n_neut;i++) {
	double line[2]={rho_neut[i],E_neut[i]};
	tneut->line_of_data(2,line);
      }
      table_alloc=true;
      return 0;
    }

    /// Return the internal table 
    table<> &get_nuc_table() {
      return *tnuc;
    }

    /// Return the internal table 
    table<> &get_neut_table() {
      return *tneut;
    }

  };

}

#endif
