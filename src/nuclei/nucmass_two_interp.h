/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2025-2026, Andrew W. Steiner

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
#ifndef NUCMASS_TWO_INTERP_H
#define NUCMASS_TWO_INTERP_H

/** \file nucmass_two_interp.h
    \brief File defining \ref o2scl::nucmass_two_interp
*/

#include <cmath>
#include <type_traits>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_file.h>

namespace o2scl {

  /** \brief Compile-time traits used by \ref nucmass_two_interp to
      decide, at compile time, whether its \c interp_t template
      parameter can be serialized to HDF5

      \ref nucmass_two_interp::hdf_output() and \ref
      nucmass_two_interp::hdf_input() only serialize \ref
      nucmass_two_interp::interp when \c interp_t provides matching
      <tt>hdf_output(hdf_file &, std::string) const</tt> and
      <tt>hdf_input(hdf_file &, std::string)</tt> methods (as \ref
      interpm_libtorch does); \c interp_t types without them (e.g.
      \ref interpm_idw) remain usable, they just can't have their
      state saved or loaded. \c has_hdf_output and \c has_hdf_input
      detect this at compile time (via the standard \c void_t
      SFINAE idiom) so that \ref nucmass_two_interp::hdf_output()/
      hdf_input() can guard the relevant code with <tt>if
      constexpr</tt>, which keeps the discarded branch from ever
      being instantiated (and so from needing to compile) for \c
      interp_t types that lack these methods. This matters because,
      unlike an ordinary member function template, a *virtual*
      member function of a class template (which \ref
      nucmass_two_interp::hdf_output()/hdf_input() must be, to
      override \ref nucmass_fit_base's virtual functions) is
      instantiated as soon as the enclosing class itself is, not
      only when actually called -- so an unconditional call to \c
      interp.hdf_output() here would otherwise fail to compile for
      any \c interp_t without one, even if that code path is never
      reached at runtime.
      \comment
      detail namespace, not documented individually via \\name since
      these are implementation plumbing rather than part of the
      public API.
      \endcomment
  */
  namespace detail {

    template<class T, class=void>
    struct has_hdf_output : std::false_type {};

    template<class T>
    struct has_hdf_output
      <T,std::void_t<decltype
                     (std::declval<const T &>().hdf_output
                      (std::declval<o2scl_hdf::hdf_file &>(),
                       std::declval<std::string>()))> > :
      std::true_type {};

    template<class T, class=void>
    struct has_hdf_input : std::false_type {};

    template<class T>
    struct has_hdf_input
      <T,std::void_t<decltype
                     (std::declval<T &>().hdf_input
                      (std::declval<o2scl_hdf::hdf_file &>(),
                       std::declval<std::string>()))> > :
      std::true_type {};

  }

  /** \brief Nuclear mass formula combining a base formula with an
      interpolator fit to the residual

      This class owns a base mass formula \ref base (of type \c
      fit_t), fit first, e.g. via \ref nucmass_fit::fit(), and an
      interpolator \ref interp (of type \c interp_t), fit second,
      via \ref fit_interp(), to the residual between experiment and
      \ref base's own prediction. Once \ref interp_trained is set
      (by \ref fit_interp()), \ref mass_excess_d() returns the sum
      of the two; before that, it returns just \ref base's own
      prediction, so that a stage-1-only fit of \ref base (e.g. via
      \ref nucmass_fit::fit() called directly on this object, as
      \ref nucmass_fit_iso::fit() does) doesn't require \ref interp
      to have any data yet.

      \comment
      Both \c fit_t and \c interp_t are owned as value (not pointer)
      members, so that the compiler-generated copy constructor --
      and thus \ref clone() -- performs a correct, independent deep
      copy of both. This is safe specifically because \ref clone()
      is only ever used (e.g. by \ref nucmass_fit_iso::fit(), before
      either \ref fit_fun() or \ref fit_interp() has been called) on
      an as-yet-untrained object; a shallow/pointer-based design (as
      this class used before) would instead leave every clone
      aliasing the same underlying base formula and interpolator,
      defeating the whole point of fitting each chain independently.
      \endcomment

      \c fit_t must be a \ref nucmass_fit_base descendant, providing
      \c nfit, \c fit_fun(), and \c guess_fun(). \c interp_t must
      provide <tt>eval(const vec_t &,vec_t &) const</tt> and
      <tt>set_data(size_t,size_t,size_t,mat_x_t &,mat_y_t &)</tt>,
      following the interface of \ref interpm_base (e.g. \ref
      interpm_idw or \ref interpm_libtorch).

      \comment
      The default template
      arguments, \ref nucmass_dz_fit and \ref interpm_idw, require no
      external dependencies (in particular, no LibTorch), so that
      this class and \ref nucmass_fit_iso can be exercised without
      them; a LibTorch-backed configuration would instead use, e.g.,
      <tt>nucmass_two_interp<nucmass_semi_empirical,
      interpm_libtorch<> ></tt>.
      \endcomment
  */
  template<class fit_t=nucmass_dz_fit, class interp_t=interpm_idw<> >
  class nucmass_two_interp : public nucmass_fit_base {

  public:

    nucmass_two_interp() {
      nfit=base.nfit;
      interp_trained=false;
    }

    virtual ~nucmass_two_interp() {
    }

    /// The base nuclear mass formula, fit first
    fit_t base;

    /// The interpolator, fit second, to the residual left by \ref base
    interp_t interp;

    /** \brief The training table backing \ref interp's data, kept
        alive for the lifetime of this object

        Many \c interp_t implementations (e.g. \ref interpm_idw, via
        its default \c mat_x_t/\c mat_y_t of \ref
        const_matrix_view_table) store their training data as a
        non-owning VIEW into a \ref table, rather than copying it;
        \ref interp_t::set_data() only swaps the (lightweight) view
        objects themselves. \ref fit_interp() therefore builds its
        training table here, as a member, rather than as a local
        variable that would be destroyed on return -- doing the
        latter would leave \ref interp holding dangling references
        as soon as \ref fit_interp() returned, silently corrupting
        every subsequent \ref eval() call.
    */
    table<> fit_tab;

    /** \brief True once \ref fit_interp() has successfully trained
        \ref interp (default false)

        Until then, \ref mass_excess_d() returns just \ref base's
        own prediction, ignoring \ref interp entirely -- most \c
        interp_t implementations (e.g. \ref interpm_idw) throw if
        asked to \c eval() before any data has been set via \ref
        interp_t::set_data(), which would otherwise happen during
        the base-formula-only stage-1 fit (\ref
        nucmass_fit::fit()'s chi-squared evaluation calls this
        object's own \ref mass_excess(), not \ref base's).
    */
    bool interp_trained;

    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x) {
      return base.fit_fun(nv,x);
    }

    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x) const {
      return base.guess_fun(nv,x);
    }

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess_d(double Z, double N) {
      double ret=base.mass_excess_d(Z,N);
      if (interp_trained) {
        ubvector x(2);
        x[0]=Z;
        x[1]=N;
        ubvector y(1);
        interp.eval(x,y);
        ret+=y[0];
      }
      return ret;
    }

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N) {
      return mass_excess_d(Z,N);
    }

    /// Return the type, \c "nucmass_two_interp".
    virtual const char *type() const { return "nucmass_two_interp"; }

    /** \brief Store the fit parameters (and, if trained, the
        interpolator) in a named HDF5 group

        Overrides \ref nucmass_fit_base::hdf_output() to store the
        same "o2scl_type"/"nfit"/"params" entries the default
        implementation would (\ref base's current parameters, via
        \ref guess_fun()), plus \ref interp_trained and, if true,
        \ref interp itself, in a nested "interp" subgroup, via \c
        interp_t's own <tt>hdf_output()</tt>. Only actually stores
        \ref interp if \c interp_t provides an <tt>hdf_output()</tt>
        method (see \ref detail::has_hdf_output); for \c interp_t
        types that don't (e.g. \ref interpm_idw), \ref interp_trained
        is still stored, but \ref interp's state is silently
        omitted, so a round trip through \ref hdf_input() would
        reproduce only \ref base's fit, not \ref interp's -- callers
        needing a full round trip should use an \c interp_t with
        HDF5 support (e.g. \ref interpm_libtorch).
    */
    virtual int hdf_output(o2scl_hdf::hdf_file &hf,
                            std::string name) const {

      if (!hf.has_write_access()) {
        O2SCL_ERR2("File not opened with write access in ",
                   "nucmass_two_interp::hdf_output().",exc_efailed);
      }

      // Start group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Add typename, fit-parameter count, and base's current
      // parameter values, exactly as the default
      // nucmass_fit_base::hdf_output() does
      hf.sets_fixed("o2scl_type",type());
      hf.set_szt("nfit",nfit);
      ubvector p(nfit);
      guess_fun(nfit,p);
      hf.setd_vec_copy("params",p);

      // Additionally store whether interp has been trained, and,
      // if interp_t supports it, interp itself
      hf.seti("interp_trained",interp_trained ? 1 : 0);
      if constexpr (detail::has_hdf_output<interp_t>::value) {
        if (interp_trained) {
          interp.hdf_output(hf,"interp");
        }
      }

      // Close group
      hf.close_group(group);

      // Return location to previous value
      hf.set_current_id(top);

      return 0;
    }

    /** \brief Load the fit parameters (and, if present, the
        interpolator) from a named HDF5 group written by \ref
        hdf_output()

        See \ref hdf_output() for the conditions under which \ref
        interp's state is actually stored (and so can be reloaded
        here); if it wasn't (because \c interp_t lacks an
        <tt>hdf_input()</tt> method, see \ref detail::has_hdf_input),
        \ref interp_trained is loaded as stored by \ref hdf_output(),
        but if that was true, \ref interp itself is left untrained
        (i.e. \ref interp_trained is forced back to false), since
        there is nothing to load it from.
    */
    virtual void hdf_input(o2scl_hdf::hdf_file &hf, std::string name) {

      // Open main group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Check typename
      std::string type2;
      hf.gets_fixed("o2scl_type",type2);
      if (type2!=((std::string)type())) {
        O2SCL_ERR2("Typename in HDF group does not match class in ",
                   "nucmass_two_interp::hdf_input().",exc_einval);
      }

      // Check the fit-parameter count
      size_t nfit2;
      hf.get_szt("nfit",nfit2);
      if (nfit2!=nfit) {
        O2SCL_ERR2("Parameter count in HDF group does not match nfit in ",
                   "nucmass_two_interp::hdf_input().",exc_einval);
      }

      // Load and set base's parameter values
      ubvector p;
      hf.getd_vec_copy("params",p);
      fit_fun(nfit,p);

      // Load whether interp was trained, and, if interp_t supports
      // it and it was, interp itself
      int it_int;
      hf.geti("interp_trained",it_int);
      interp_trained=(it_int!=0);
      if constexpr (detail::has_hdf_input<interp_t>::value) {
        if (interp_trained) {
          std::string iname="interp";
          interp.hdf_input(hf,iname);
        }
      } else {
        // interp_t can't be reloaded even if it was trained and
        // stored (see hdf_output()'s doc comment); leave interp
        // untrained rather than silently returning base-only
        // predictions while claiming interp_trained is true.
        interp_trained=false;
      }

      // Close group
      hf.close_group(group);

      // Return location to previous value
      hf.set_current_id(top);

      return;
    }

    /** \brief Return a new, independent copy of this object

        Safe because \ref base and \ref interp are owned as value
        members: the compiler-generated copy constructor used here
        deep-copies both, rather than aliasing shared external
        state as a raw-pointer design would.
    */
    virtual nucmass_two_interp<fit_t,interp_t> *clone() {
      return new nucmass_two_interp<fit_t,interp_t>(*this);
    }

    /** \brief Set up and fit \ref base and \ref interp, in that
        order, to the full AME 2020 table

        A demonstration/default configuration: loads AME 2020, fits
        \ref base to it with \ref nucmass_fit, then calls \ref
        fit_interp() to train \ref interp on the residual.
    */
    void set_default() {

      nucmass_ame ame;
      ame.load("20");

      nucmass_fit nf;
      nucdist_set(nf.dist,ame);
      double res;
      nf.def_mmin.verbose=2;
      nf.fit(*this,res);

      fit_interp(nf.dist);

      return;
    }

    /** \brief Fit the interpolator \ref interp to the residual
        between the experimental mass excess and \ref base, for the
        nuclei in \c dist

        \ref base should already be fit (e.g. via \ref
        nucmass_fit::fit(), which is what \ref fit_fun() delegates
        to) before this is called; this function does not fit \ref
        base itself, only \ref interp. The interpolator is trained
        on <tt>me_exp-me_th</tt> (the residual), not <tt>me_th</tt>
        directly, since \ref mass_excess_d() adds the interpolator's
        output on top of \ref base's own prediction rather than
        replacing it. Overrides \ref
        nucmass_fit_base::fit_interp() so that \ref
        nucmass_fit_iso::fit_interp_pass() can train this second
        stage for every chain without needing to know which chains
        are actually \c nucmass_two_interp objects.
    */
    virtual int fit_interp(std::vector<nucleus> &dist) {

      nucmass_fit nf;
      nf.dist=dist;
      double res;
      
      // Build directly into fit_tab (a member -- see its doc
      // comment for why this can't be a local variable here).
      fit_tab.clear();
      double max_abs_dev;
      nf.eval_table(base,res,max_abs_dev,true,fit_tab);

      // Train on the RESIDUAL (me_exp-me_th), not me_th itself: see
      // the note above and in mass_excess_d().
      size_t nl=fit_tab.get_nlines();
      fit_tab.new_column("me_resid");
      for(size_t i=0;i<nl;i++) {
        fit_tab.set("me_resid",i,fit_tab.get("me_exp",i)-
                     fit_tab.get("me_th",i));
      }

      const_matrix_view_table<> ix(fit_tab,{"Z","N"});
      std::vector<std::string> out_cols={"me_resid"};
      const_matrix_view_table<> iy(fit_tab,out_cols);
      interp.set_data(2,out_cols.size(),fit_tab.get_nlines(),ix,iy);

      interp_trained=true;

      return 0;
    }

  };

}

#endif
