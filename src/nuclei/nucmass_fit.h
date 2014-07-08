/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef NUCMASS_FIT_H
#define NUCMASS_FIT_H

/** \file nucmass_fit.h
    \brief File defining \ref o2scl::nucmass_fit
*/

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mmin.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucdist.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Fit a nuclear mass formula

      There is an example of the usage of this class given in 
      \ref ex_nucmass_fit_sect.

      \future Convert to a real fit with errors and covariance, etc.
   */
  class nucmass_fit {

  public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;

    nucmass_fit();
    
    virtual ~nucmass_fit() {};

    /// \name Fitting method
    //@{
    /// Current fitting method
    int fit_method;
    /// RMS deviation in mass excess
    static const int rms_mass_excess=0;
    /// RMS deviation in binding_energy
    static const int rms_binding_energy=1;
    /// Chi-squared for mass excess using specified uncertainties
    static const int chi_squared_me=2;
    /// Chi-squared for binding energy using specified uncertainties
    static const int chi_squared_be=3;
    //@}
    
    /// If true, then only fit doubly-even nuclei (default false)
    bool even_even;
   
    /// Minimum proton number to fit (default 8)
    int minZ;
    
    /// Minimum neutron number to fit (default 8)
    int minN;

    /// Fit the nuclear mass formula
    virtual void fit(nucmass_fit_base &n, double &res);
    
    /** \brief Evaluate quality without fitting
     */
    virtual void eval(nucmass &n, double &res);

    /** \brief The default minimizer

	The value of def_mmin::ntrial is automatically multiplied by
	10 in the constructor because the minimization frequently
	requires more trials than the default.
    */
    mmin_simp2<> def_mmin;
    
    /// Change the minimizer for use in the fit
    void set_mmin(mmin_base<> &umm) {
      mm=&umm;
      return;
    }
    
    /** \brief Select the experimental nuclei to fit
     */
    std::vector<nucleus> dist;

    /** \brief Set the fit uncertainties (in MeV)
     */
    template<class vec_t> void set_uncerts(vec_t &u) {
      size_t nv=u.size();
      set_uncerts(nv,u);
      return;
    }

    /** \brief Set the fit uncertainties (in MeV) from the first \c nv
	elements of \c u
     */
    template<class vec_t> void set_uncerts(size_t nv, vec_t &u) {
      if (nv==0) {
	O2SCL_ERR2("Tried to give zero uncertainties in nucmass_fit::",
		   "set_uncerts().",exc_efailed);
      }
      if (uncs.size()>0) uncs.clear();
      uncs.resize(nv);
      vector_copy(nv,u,uncs);
      return;
    }
    
    /** \brief Evaluate isospin dependence of fit quality

	\todo More documentation and compute uncertainty
     */
    void eval_isospin_beta(nucmass &n, ubvector_int &n_qual,
			   ubvector &qual, int max_iso=20);
    
    /** \brief Evaluate isospin dependence of fit quality
     */
    void eval_isospin(nucmass &n, ubvector_int &n_qual,
		      ubvector &qual, int min_iso=-8, int max_iso=60);

    /** \brief The function to minimize
     */
    virtual double min_fun(size_t nv, const ubvector &x);

  protected:

#ifndef DOXYGEN_NO_O2NS

    /// Uncertainties
    ubvector uncs;
    
    /// The pointer to the minimizer
    mmin_base<> *mm;
    
    /** \brief The nuclear mass formula to fit to

	This pointer is set by fit() and eval().
     */
    nucmass_fit_base *nmf;
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
