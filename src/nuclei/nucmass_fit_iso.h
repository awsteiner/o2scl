/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2026, Andrew W. Steiner

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
#ifndef O2SCL_NUCMASS_FIT_ISO_H
#define O2SCL_NUCMASS_FIT_ISO_H

/** \file nucmass_fit_iso.h
    \brief File defining \ref o2scl::nucmass_fit_iso
*/

#include <map>
#include <memory>

#include <o2scl/nucmass.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/nucdist.h>
#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

namespace o2scl {

  /** \brief Nuclear mass formula built from independent fits to
      each isotopic chain

      Rather than fitting a single set of parameters to the entire
      nuclear chart, this class fits an independent copy of a
      user-specified mass formula (a \ref nucmass_fit_base
      descendant which overrides \ref nucmass_fit_base::clone()) to
      each isotopic chain (i.e. to each proton number Z separately),
      using only the isotopes of the neighboring chains from Z-x to
      Z+x (inclusive) as calibration data for that chain's fit. This
      is motivated by the observation (see the O2sclpy example
      <tt>nucmass_fit_range.py</tt>) that many global mass formulas
      extrapolate poorly far from the region they were fit to, and
      that fitting each chain to only its local neighborhood often
      gives a better local description near that chain than a
      single set of parameters fit to the entire chart, at the
      price of losing global smoothness (there is generally a small
      discontinuity in the fitted mass surface between neighboring
      chains).

      Because each chain is fit independently, this is not a single
      smooth mass formula; \ref mass_excess() simply looks up the
      chain fit for the requested Z (if one was successfully fit by
      \ref fit()) and evaluates that chain's fitted formula at
      (Z,N). Chains for which \ref fit() found insufficient
      calibration data, or whose fit failed to converge, are
      omitted, and \ref is_included() returns false for any (Z,N)
      whose chain wasn't fit (or which lies outside the range that
      chain's underlying formula itself reports as included).

      \note Only mass formula classes descended from \ref
      nucmass_fit_base whose \ref nucmass_fit_base::clone()
      function overrides the (default, error-throwing) base class
      implementation can be used with this class. As of this
      writing, \ref nucmass_dz_fit_33 and \ref nucmass_two_interp
      (for any choice of its template parameters) provide one.
  */
  class nucmass_fit_iso : public nucmass_table {

  protected:

    /** \brief The fitted formula object for each isotopic chain,
        indexed by proton number Z
    */
    std::map<int,std::shared_ptr<nucmass_fit_base> > chains;

  public:

    nucmass_fit_iso();

    virtual ~nucmass_fit_iso() {
    }

    /// Number of neighboring chains fit on either side (default 4)
    int x;

    /** \brief Minimum number of training isotopes required per fit
        parameter before a chain is even attempted (default 3)
    */
    int min_data_per_param;

    /// Verbosity parameter (default 0)
    int verbose;

    /** \brief Number of OpenMP threads to use in \ref fit() and
        \ref fit_interp_pass() (default 1)

        Fitting one isotopic chain is independent of fitting any
        other, so both \ref fit() and \ref fit_interp_pass() can
        split their loop over chains across \ref n_threads OpenMP
        threads. Each thread gets its own copy of \ref mf (so
        \ref mf's own settings, like <tt>mf.def_mmin.ntrial</tt>,
        should be set before calling \ref fit(), not after), and
        each chain already has its own independently cloned formula
        object, so there is no shared mutable state between chains
        during the fit itself.

        \warning If the formula being fit uses a library with its
        own global, not-thread-local random state -- for example,
        \ref interpm_libtorch, whose training draws from LibTorch's
        default generator -- setting \ref n_threads greater than 1
        is only safe if that library's own state is not being
        concurrently mutated in a way that would let one chain's
        random draws interfere with another's. See \ref
        interpm_libtorch's documentation for how it protects
        against this specifically. If in doubt for some other
        formula, test with a small Z range first and confirm the
        fitted parameters don't depend on \ref n_threads.
    */
    size_t n_threads;

    /** \brief The internal fitting object used by \ref fit(),
        exposed so its minimizer (\ref nucmass_fit::def_mmin) can be
        tuned before calling \ref fit(), e.g.
        <tt>mf.def_mmin.ntrial=60000</tt> for a formula (like \ref
        nucmass_dz_fit_33) whose minimization needs a larger trial
        budget than the default. The constructor sets
        <tt>mf.def_mmin.err_nonconv=false</tt> so that a
        non-converging chain is skipped by \ref fit() rather than
        aborting the process; this should not be changed back to
        \c true.
    */
    nucmass_fit mf;

    /** \brief Fit a copy of \c proto independently to each
        isotopic chain with sufficient calibration data in \c src

        For each proton number Z from \c minZ to \c maxZ
        (inclusive), a fresh copy of \c proto (via \ref
        nucmass_fit_base::clone()) is fit (via \ref mf) to the
        isotopes of \c src with proton number in [Z-x,Z+x]. If that
        subset has fewer than <tt>min_data_per_param*proto.nfit</tt>
        isotopes, or if the fit doesn't converge within \ref mf's
        minimizer trial budget (<tt>mf.def_mmin.ntrial</tt>), that
        chain is skipped. Any chain successfully fit here replaces a
        previously-fit chain at the same Z, if one exists; call
        clear() first to discard all previously-fit chains rather
        than only those refit this time. Returns the number of
        chains successfully (re-)fit.

        If \ref n_threads is greater than 1 and O2scl was compiled
        with OpenMP support, the loop over Z is split across \ref
        n_threads threads; see \ref n_threads for what this
        requires of \c proto.
    */
    size_t fit(nucmass_fit_base &proto, nucmass &src,
               int minZ=8, int maxZ=118);

    /** \brief For every fitted chain with a secondary fit stage
        (see \ref nucmass_fit_base::fit_interp()), train that stage
        on the same Z-x to Z+x calibration window used by \ref fit()

        \ref fit() only fits chains through \ref mf's simplex
        minimizer, which touches only the parameters exposed via
        \ref nucmass_fit_base::fit_fun()/guess_fun() -- for a \ref
        nucmass_two_interp chain, that's just its base formula (see
        \ref nucmass_two_interp::fit_fun()), since an interpolator
        such as \ref interpm_libtorch is trained separately, via
        gradient descent rather than a derivative-free simplex
        search. Call this as a second pass, after \ref fit(), to
        train each chain's secondary stage (via \ref
        nucmass_fit_base::fit_interp(), overridden by classes like
        \ref nucmass_two_interp) on the residual between experiment
        and the chain's now-fit primary formula. This calls \ref
        nucmass_fit_base::fit_interp() on every chain regardless of
        concrete type -- the default implementation is a no-op that
        reports "not applicable", so chains without a secondary
        stage are automatically skipped without any \c dynamic_cast
        or other type inspection here. Returns the number of chains
        whose secondary stage was (re-)trained.

        As in \ref fit(), if \ref n_threads is greater than 1 and
        O2scl was compiled with OpenMP support, this loop over
        already-fitted chains is also split across \ref n_threads
        threads. See the \ref n_threads warning about formulas
        (like \ref interpm_libtorch) which rely on a library with
        its own global random state.
    */
    size_t fit_interp_pass(nucmass &src);

    /// Remove all fitted chains
    virtual void clear() {
      chains.clear();
      n=0;
      return;
    }

    /** \brief Return true if a chain was fit for \c Z and that
        chain's underlying formula reports (Z,N) as included
    */
    virtual bool is_included(int Z, int N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);

    /// Return the type, \c "nucmass_fit_iso".
    virtual const char *type() const { return "nucmass_fit_iso"; }

    /** \brief Store all fitted chains in a named HDF5 group

        Writes \c x, \c min_data_per_param, and the sorted list of
        fitted Z values as metadata, then calls each chain's own
        (possibly overridden) \ref nucmass_fit_base::hdf_output() to
        write that chain's data into its own nested subgroup, named
        <tt>chain_&lt;Z&gt;</tt>. Because each chain serializes
        itself polymorphically, this works unchanged even for a
        future nucmass_fit_base descendant whose fit state isn't a
        flat parameter vector (e.g. a per-isotope neural network) --
        nucmass_fit_iso never needs to know how a chain's data is
        actually represented.
    */
    int hdf_output(o2scl_hdf::hdf_file &hf, std::string name) const;

    /** \brief Load fitted chains from a named HDF5 group written by
        \ref hdf_output()

        A fresh copy of \c proto (via \ref
        nucmass_fit_base::clone()) is made for each Z listed in the
        group, and that copy's own \ref nucmass_fit_base::hdf_input()
        is used to load its nested subgroup. \c proto must be the
        same (concrete) type used when the file was written; there
        is no way to recover the concrete type from the file alone,
        since \ref nucmass_fit_base::clone() has no generic
        by-name factory. Any chains already present are discarded
        first (as if clear() had been called).
    */
    void hdf_input(o2scl_hdf::hdf_file &hf, nucmass_fit_base &proto,
                   std::string name);

    /** \brief Load the first \ref nucmass_fit_iso found in \c hf,
        and return its name (see \ref hdf_input())
    */
    int hdf_input_n(o2scl_hdf::hdf_file &hf, nucmass_fit_base &proto,
                    std::string &name);

  };

}

#endif
