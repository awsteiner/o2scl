/*
   ───────────────────────────────────────────────────────────────────

   Copyright (C) 2010-2026, Edwin van Leeuwen and Andrew W. Steiner

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
#ifndef O2SCL_DIFF_EVO_PARA_H
#define O2SCL_DIFF_EVO_PARA_H

/** \file diff_evo_para.h
    \brief File defining \ref o2scl::diff_evo_para
*/

#include <vector>
#include <limits>
#include <exception>

#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

#include <o2scl/diff_evo_adapt.h>

namespace o2scl {

  /** \brief Multidimensional minimization by differential
      evolution (OpenMP version)

      This class parallelizes a single differential-evolution
      minimization across OpenMP threads.

      \verbatim embed:rst
      Differential evolution's population already provides natural
      per-generation parallel work: in the classic, generational
      ("rand/1/bin") update rule of [Storn97]_ , every trial vector in
      a generation is constructed from and evaluated against a fixed
      snapshot of the population taken at the start of that
      generation, and all accepted replacements are applied only
      afterwards. Because the evaluation of one agent's trial vector
      never depends on another agent's trial vector within the same
      generation, that evaluation loop can be split across threads
      without any risk of a data race, as long as the (comparatively
      cheap) application of the accepted trials back into the
      population happens serially afterward.
      \endverbatim

      This class always operates in this generational mode (see \ref
      o2scl::diff_evo_adapt::generational), regardless of what the
      user sets that member to, and reuses essentially all of \ref
      o2scl::diff_evo_adapt unchanged: it only overrides \ref
      o2scl::diff_evo_adapt::compute_trials() (to evaluate the
      generation's trial vectors in an OpenMP-parallel loop, each
      thread using its own random number generator and its own copy of
      \c func_t) and \ref o2scl::diff_evo_adapt::select_pop_size() (to
      round the automatically-selected population size up to a
      multiple of \ref n_threads, so that no thread is left idle
      during a generation). The function \ref
      o2scl::diff_evo_adapt::apply_trials(), the self-adaptive f and
      cr selection, and the outer generation loop in \ref
      o2scl::diff_evo_adapt::mmin() are all inherited unchanged.
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double>,
    class init_funct_t=mm_funct>
    class diff_evo_para :
    public diff_evo_adapt<func_t, vec_t, init_funct_t>,
    public mmin_parallel_base {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /** \brief The number of OpenMP threads (default 1)

        After a call to \ref mmin(), this is updated to the actual
        number of threads OpenMP granted, which may be smaller than
        the requested value.
    */
    size_t n_threads;

    diff_evo_para() : diff_evo_adapt<func_t,vec_t,init_funct_t>() {
      n_threads=1;
      this->generational=true;
    }

    /** \brief Calculate the minimum \c fmin of \c func w.r.t the
        array \c x of size \c nvar, using \ref n_threads OpenMP
        threads

        A separate copy of \c func is made for each thread, so any
        state \c func_t captures by value is duplicated across
        threads, while state captured by reference or pointer is
        shared. If \c func_t is not internally thread-safe, then the
        user should ensure that any such state is either duplicated in
        \c func_t's copy constructor or made thread-local before
        calling this function.
    */
    virtual int mmin(size_t nvar, vec_t &x0, double &fmin,
                      func_t &func) {

      // Set the number of threads, updating n_threads to the
      // actual value OpenMP grants us
#ifdef O2SCL_SET_OPENMP
      omp_set_num_threads(n_threads);
#pragma omp parallel
      {
        n_threads=omp_get_num_threads();
      }
#else
      n_threads=1;
#endif

      // One copy of the function object per thread
      vfuncs.assign(n_threads,func);

      // One independent random number generator per thread, so
      // that trial-vector construction in different threads doesn't
      // draw from the same, non-thread-safe stream
      vrng.resize(n_threads);
      unsigned long int s=time(0);
      for(size_t it=0;it<n_threads;it++) {
        vrng[it].set_seed(s+it);
      }

      // Generational update is required for race-free parallel
      // evaluation of a generation
      this->generational=true;

      return diff_evo_adapt<func_t,vec_t,init_funct_t>::mmin
        (nvar,x0,fmin,func);
    }

    /// Return string denoting type ("diff_evo_para")
    virtual const char *type() { return "diff_evo_para"; }

    /** \brief Return \ref n_threads, the number of threads this
        minimizer may use to concurrently evaluate the function
        being minimized (implements \ref o2scl::mmin_parallel_base)
    */
    virtual size_t mmin_n_threads() const { return n_threads; }

  protected:

    /// Independent random number generators, one per thread
    std::vector<rng<> > vrng;

    /// Independent copies of the function object, one per thread
    std::vector<func_t> vfuncs;

    /** \brief Auto-select the population size based on
        dimensionality, as in the parent class, then round up to
        the next multiple of \ref n_threads so that no thread is
        idle during a generation's parallel evaluation loop

        As in the parent class, this only applies when the user has
        not set \ref o2scl::diff_evo::pop_size directly.
    */
    virtual size_t select_pop_size(size_t nvar) {
      if (this->pop_size!=0) {
        return this->pop_size;
      }
      size_t pop_size_loc=
        diff_evo_adapt<func_t,vec_t,init_funct_t>::select_pop_size(nvar);
      if (n_threads>1) {
        pop_size_loc=
          ((pop_size_loc+n_threads-1)/n_threads)*n_threads;
      }
      return pop_size_loc;
    }

    /** \brief Evaluate the trial vectors for the whole generation
        in an OpenMP-parallel loop, each thread using its own
        random number generator (\ref vrng) and its own copy of the
        function object (\ref vfuncs)

        Each loop iteration only reads the read-only snapshot \c
        pop_src and only writes to the entries of \c trial, \c
        fmin_trial, \c f_trial, and \c cr_trial indexed by its own
        agent number \c x, so no two threads ever write to the same
        location and this is safe to parallelize directly. The
        function object \c func passed in by \ref
        o2scl::diff_evo_adapt::run_generation() is ignored here in
        favor of the per-thread copies in \ref vfuncs, which were
        set up in \ref mmin().

        As in \ref o2scl::diff_evo_adapt::compute_trials(), any
        exception thrown while evaluating a trial vector is caught
        and treated as a very bad (but finite) fitness value rather
        than allowed to propagate -- this matters even more here
        than in the serial base class, since an exception that
        escapes this OpenMP parallel region uncaught is immediately
        fatal (<tt>std::terminate()</tt>), and several threads
        hitting an uncaught exception at once is exactly the
        "terminate called recursively" crash pattern.
    */
    virtual void compute_trials(size_t nvar, size_t pop_size_loc,
                                 const vec_t &pop_src, func_t &func,
                                 std::vector<vec_t> &trial,
                                 ubvector &fmin_trial,
                                 vec_t &f_trial, vec_t &cr_trial) {

#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
      {
#pragma omp for
#endif
        for (size_t x=0;x<pop_size_loc;++x) {

#ifdef O2SCL_SET_OPENMP
          size_t ithread=omp_get_thread_num();
#else
          size_t ithread=0;
#endif

          double f_x, cr_x;
          this->adapt_f_cr(x,f_x,cr_x,vrng[ithread]);
          f_trial[x]=f_x;
          cr_trial[x]=cr_x;
          trial[x].resize(nvar);
          this->build_trial_vector(nvar,x,pop_size_loc,pop_src,f_x,cr_x,
                                    trial[x],vrng[ithread]);
          try {
            fmin_trial[x]=vfuncs[ithread](nvar,trial[x]);
          } catch (const std::exception &e) {
            fmin_trial[x]=std::numeric_limits<double>::max();
          }
        }
#ifdef O2SCL_SET_OPENMP
      }
#endif

    }

  };

}

#endif
