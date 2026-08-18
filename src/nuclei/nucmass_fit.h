/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2026, Andrew W. Steiner
  
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
#ifndef NUCMASS_FIT_H
#define NUCMASS_FIT_H

/** \file nucmass_fit.h
    \brief File defining \ref o2scl::nucmass_fit
*/

#include <memory>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mmin.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucdist.h>
#include <o2scl/fit_nonlin.h>
#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

namespace o2scl {

  /** \brief Fit a nuclear mass formula

      \verbatim embed:rst
      There is an example of the usage of this class given in 
      the :ref:`Nuclear mass fit example`.
      \endverbatim

      \future Convert to a real fit with errors and covariance, etc.
   */
  class nucmass_fit {

  public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;

    nucmass_fit();

    virtual ~nucmass_fit() {};

    /** \brief Copy constructor

        \ref def_mmin is an \ref mmin_simp2 object, which
        intentionally disallows copying (as do most O2scl
        minimizers) since a blind field-by-field copy of its
        internal scratch state (e.g. its simplex, mid-search) would
        not be meaningful. This constructor therefore doesn't copy
        \ref def_mmin itself; instead, as \ref mmin_base's own copy
        constructor does for its own copyable descendants, it copies
        just \ref def_mmin's tunable settings (verbose, ntrial,
        tol_rel, tol_abs, err_nonconv) onto a freshly-constructed
        \ref def_mmin, alongside \c nf's other settings (\ref
        fit_method, \ref even_even, \ref minZ, \ref minN, \ref dist,
        and the protected uncertainty vector and minimizer pointer
        set via \ref set_uncerts() and \ref set_mmin()). This is
        primarily intended to let callers such as \ref
        nucmass_fit_iso give each of several OpenMP threads its own
        independent copy of a tuned \ref nucmass_fit object.

        \ref mm defaults (in the ordinary constructor) to
        <tt>&def_mmin</tt>, i.e. "use my own \ref def_mmin". If \c
        nf's \ref mm is still pointing at \c nf's own \ref def_mmin
        (the common case, when \ref set_mmin() was never called),
        this copy points the new object's \ref mm at its own \ref
        def_mmin too, not \c nf's -- otherwise every copy made this
        way would still share (and race on) \c nf's single \ref
        def_mmin regardless of how many independent copies exist.
        If \c nf's \ref mm was redirected to some other, external
        minimizer via \ref set_mmin(), that pointer is copied as-is;
        note that doing so means several thread-local copies of \c
        nf would then share (and need to independently be safe for
        concurrent use of) that same external minimizer object --
        \ref n_threads greater than 1 in \ref nucmass_fit_iso is
        only guaranteed race-free through this copy when using the
        default \ref def_mmin.
    */
    nucmass_fit(const nucmass_fit &nf) {
      fit_method=nf.fit_method;
      even_even=nf.even_even;
      minZ=nf.minZ;
      minN=nf.minN;
      dist=nf.dist;
      uncs=nf.uncs;
      mm=(nf.mm==&nf.def_mmin) ? &def_mmin : nf.mm;
      def_mmin.verbose=nf.def_mmin.verbose;
      def_mmin.ntrial=nf.def_mmin.ntrial;
      def_mmin.tol_rel=nf.def_mmin.tol_rel;
      def_mmin.tol_abs=nf.def_mmin.tol_abs;
      def_mmin.err_nonconv=nf.def_mmin.err_nonconv;
    }

    /** \brief Copy from operator=

        See the copy constructor above for exactly what is (and
        isn't) copied, and the caveat about \ref mm when \ref
        set_mmin() has redirected it to an external minimizer.
    */
    nucmass_fit &operator=(const nucmass_fit &nf) {
      if (this!=&nf) {
        fit_method=nf.fit_method;
        even_even=nf.even_even;
        minZ=nf.minZ;
        minN=nf.minN;
        dist=nf.dist;
        uncs=nf.uncs;
        mm=(nf.mm==&nf.def_mmin) ? &def_mmin : nf.mm;
        def_mmin.verbose=nf.def_mmin.verbose;
        def_mmin.ntrial=nf.def_mmin.ntrial;
        def_mmin.tol_rel=nf.def_mmin.tol_rel;
        def_mmin.tol_abs=nf.def_mmin.tol_abs;
        def_mmin.err_nonconv=nf.def_mmin.err_nonconv;
      }
      return *this;
    }

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
    /// Fit the mass excess and the neutron separation energy
    static const int rms_me_Sn=4;
    /// Fit the mass excess and the one and two neutron separation energies
    static const int rms_me_Sn_S2n=5;
    /// Fit the mass excess and the one and two neutron separation energies
    static const int rms_me_S2n=6;
    //@}
    
    /// If true, then only fit doubly-even nuclei (default false)
    bool even_even;
   
    /// Minimum proton number to fit (default 8)
    int minZ;
    
    /// Minimum neutron number to fit (default 8)
    int minN;

    /** \brief Fit the nuclear mass formula

        If \ref set_mmin() has been used to install a minimizer
        which inherits from \ref mmin_parallel_base and reports
        (via \ref mmin_parallel_base::mmin_n_threads()) that it may
        call the function being minimized from more than one thread
        at once, this method gives each of those threads its own
        private clone of \c n (via \ref
        nucmass_fit_base::clone()), so that concurrent calls into
        \c n's \ref nucmass_fit_base::fit_fun() from different
        threads don't race by writing to and reading from \c n's
        shared internal parameter state. This requires \c n's
        concrete type to actually implement \ref
        nucmass_fit_base::clone() (the default implementation calls
        the error handler); with a serial minimizer (the default),
        no cloning happens and \c n's own \ref
        nucmass_fit_base::clone() is never called. Regardless of
        which path is taken, \c n itself always ends up holding the
        best-fit parameters found once this function returns.
    */
    virtual void fit(nucmass_fit_base &n, double &res);
    
    /** \brief Evaluate quality without fitting
     */
    virtual void eval(nucmass &n, double &res);

    /** \brief Evaluate quality without fitting
     */
    virtual void eval_max(nucmass &n, double &res, double &max_abs_dev);

    /** \brief Desc
     */
    void eval_table(nucmass &n, double &fmin, double &max_abs_dev,
                    bool make_table, table<> &tab);
    
    /** \brief Fit a nuclear mass formula using least squares
        and report the associated \f$ \chi^2 \f$ and 
        covariance matrix

        \note This function only works for \ref fit_method equal
        to \ref chi_squared_me or \ref chi_squared_be .
     */
    void fit_covar(nucmass_fit_base &n, 
                   double &chi2, ubmatrix &covar);

    /** \brief The form of the fitting function which is set
        for a fitting object of type \ref o2scl::fit_nonlin
     */
    double fit_covar_fun(size_t np, const ubvector &p,
                         double x, const std::vector<size_t> &Zlist,
                         const std::vector<size_t> &Nlist);
    
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
        
        \verbatim embed:rst
        .. todo:: 

           - In nucmass_fit::eval_isospin_beta(): 
             More documentation and compute uncertainty.

        \endverbatim
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

    /// Uncertainties
    ubvector uncs;
    
    /// The pointer to the minimizer
    mmin_base<> *mm;
    
    /** \brief The nuclear mass formula to fit to

        This pointer is set by fit() and eval().
     */
    nucmass_fit_base *nmf;

    /** \brief Per-thread clones of the formula being fit, used only
        while a call to \ref fit() is in progress with a parallel
        minimizer

        Empty at all other times, including while \ref eval(), \ref
        fit_covar(), and similar are running (those always operate
        serially through \ref nmf directly). See \ref fit() and
        \ref min_fun() for how this is used.
    */
    std::vector<std::shared_ptr<nucmass_fit_base> > thread_clones;

  };

}

#endif
