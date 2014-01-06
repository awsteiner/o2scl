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
#ifndef CFL_NJL_EOS_H
#define CFL_NJL_EOS_H

#include <iostream>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_poly.h>

#include <o2scl/constants.h> 
#include <o2scl/err_hnd.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mm_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/poly.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>

#include <o2scl/part.h>
#include <o2scl/nambujl_eos.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nambu Jona-Lasinio model with a schematic CFL 
      di-quark interaction at finite temperature
      
      The variable B0 must be set before use.
      
      The original Lagrangian is
      
      \f[
      {\cal L} =
      {\cal L}_{\mathrm{Dirac}} +
      {\cal L}_{\mathrm{4-fermion}} +
      {\cal L}_{\mathrm{6-fermion}} +
      {\cal L}_{CSC1} +
      {\cal L}_{CSC2} 
      \f]
      
      \f[
      {\cal L}_{\mathrm{Dirac}} = \bar{q}_{i \alpha}
      \left( i {\partial} \delta_{i j} \delta_{\alpha \beta} - 
      m_{i j} \delta_{\alpha \beta} - \mu_{i j,~\alpha \beta} \gamma^0
      \right) q_{j \beta}  
      \f]
      
      \f[
      {\cal L}_{\mathrm{4-fermion}} =
      G_S \sum_{a=0}^8 \left[ 
      \left( \bar{q} \lambda^a_f q \right)^2 +
      \left( \bar{q} i \gamma_5 \lambda^a_f q \right)^2 
      \right] 
      \f]
      
      \f[
      {\cal L}_{\mathrm{6-fermion}} =
      - G_{D} \left[ 
      {\mathrm det}_{i j} \, \bar{q}_{i \alpha} \left( 1 + i \gamma_5 \right) 
      q_{j \beta} +
      {\mathrm det}_{i j} \, \bar{q}_{i \alpha} \left( 1 - i \gamma_5 \right) 
      q_{j \beta} \right] \delta_{\alpha \beta} 
      \f]
      
      \f[
      {\cal L}_{CSC1} = 
      G_{DIQ} \sum_k \sum_{\gamma} \left[
      \left(\bar{q}_{i \alpha} \epsilon_{i j k} 
      \epsilon_{\alpha \beta \gamma} q^C_{j \beta}\right)
      \left(\bar{q}_{i^{\prime} \alpha^{\prime}}^C 
      \epsilon_{i^{\prime} j^{\prime} k} \epsilon_{\alpha^{\prime} 
      \beta^{\prime} \gamma} q_{j^{\prime} \beta^{\prime}}\right)\right]
      \f]
      
      \f[
      {\cal L}_{CSC2} = 
      G_{DIQ} \sum_k \sum_{\gamma} \left[
      \left(\bar{q}_{i \alpha} i \gamma_5 \epsilon_{i j k} 
      \epsilon_{\alpha \beta \gamma} q^C_{j \beta}\right)
      \left(\bar{q}_{i^{\prime} \alpha^{\prime}}^C i \gamma_5 
      \epsilon_{i^{\prime} j^{\prime} k} \epsilon_{\alpha^{\prime} 
      \beta^{\prime} \gamma} q_{j^{\prime} \beta^{\prime}}\right) \right] \,, 
      \f]
      
      where \f$ \mu \f$ is the \quark number chemical potential.  
      couplings \f$ G_S \f$, \f$ G_D \f$, and \f$ G_{DIQ} \f$
      ultra-violet three-momentum cutoff, \f$ \Lambda \f$
      
      The thermodynamic potential is
      \f[
      \Omega(\mu_i,\left<\bar{q} q\right>_i,\left< q q\right>_i,T)
      =  \Omega_{\mathrm{vac}}+\Omega_{\mathrm{stat}} + \Omega_0
      \f]
      where \f$ i \f$ runs over all nine (three colors times three
      flavors) quarks. We assume that the condensates are independent
      of color and 
      that the \quark chemical potentials
      are of the form
      \f$ \mu_Q=\mu_{\mathrm{Flavor(Q)}}+\mu_{\mathrm{Color(Q)}} \f$
      with
      \f[
      \mu_{\mathrm{red}} = \mu_3 + \mu_8/\sqrt{3}
      \quad
      \mu_{\mathrm{green}} = -\mu_3 + \mu_8/\sqrt{3}
      \quad
      \mu_{\mathrm{blue}} = -2 \mu_8 /\sqrt{3}
      \f]
      With these assumptions, the thermodynamic potential as given
      by the function thd_potential(), is a function of 12 variables
      \f[
      \Omega(\mu_u, \mu_d, \mu_s, \mu_3, \mu_8, \left<\bar{u} u\right>,
      \left<\bar{d} d\right>, \left<\bar{s} s\right>,
      \left< u d\right>, \left< u s\right>, \left< d s\right>,
      T)
      \f]
      
      The individual terms are
      \f[
      \Omega_{\mathrm{stat}} = - \frac{1}{2} \int \frac{d^3
      p}{\left(2 \pi\right)^3} \, \sum_{i=1}^{72} \left[ \frac{\lambda_i}{2} +
      T \ln{\left(1 + e^{-\lambda_i/T} \right)} \right] 
      \f]
      
      \f[
      \Omega_{\mathrm{vac}} =  
      - 2 G_S \sum_{i=u,d,s} \langle {\bar q_i} q_i \rangle^2  
      +4 G_D \left<{\bar u} u\right> \left<{\bar d} d \right> 
      \left<{\bar s} s\right>
      + \sum_k \sum_{\gamma} \frac{\left|\Delta^{k \gamma}\right|^2}{4 G_{DIQ}}
      \f]
      
      where \f$ \lambda_i \f$ are the eigenvalues of the (72 by 72) matrix
      (calculated by the function eigenvalues())
      \f[
      D = \left[
      \begin{array}{cc}
      - \gamma^0 \vec{\gamma} \cdot \vec{p} - M_{i} \gamma^0 + \mu_{i \alpha}  
      & \Delta i \gamma^0 \gamma_5 C \\
      i \Delta \gamma^0 C \gamma_5 
      & - \gamma^0 \vec{\gamma}^T \cdot \vec{p} + M_{i} \gamma^0 
      - \mu_{i \alpha}\end{array}
      \right] 
      \f]
      and \f$ C \f$ is the charge conjugation matrix (in the Dirac
      representation).
      
      The values of the various condensates are usually determined by
      the condition
      \f[
      \frac{\partial \Omega}{\left<\bar{q} q\right>_i} = 0
      \quad
      \frac{\partial \Omega}{\left<q q\right>_i} = 0
      \f]
      
      Note that setting fixed_mass to \c true and setting all of the
      gaps to zero when \c gap_limit is less than zero will reproduce
      an analog of the bag model with a momentum cutoff.
      
      The variable nambujl_eos::fromqq is automatically set to true in
      the constructor, as computations with \c fromqq=false are not
      implemented.

      \future This class internally mixes ubvector, ubmatrix, gsl_vector
      and gsl_matrix objects in a confusing and non-optimal way. Fix this. 
      \future Allow user to change derivative object? This isn't possible
      right now because the stepsize parameter of the derivative
      object is used.

      \hline
      <b>References:</b>

      Created for \ref Steiner02.
  */
  class cfl_njl_eos : public nambujl_eos {
  public:
    
    cfl_njl_eos();
    
    virtual ~cfl_njl_eos();
    
    /** \brief Set the parameters and the bag constant 'B0'
	
	This function allows the user to specify the momentum cutoff,
	\c lambda, the four-fermion coupling \c fourferm, the
	six-fermion coupling from the 't Hooft interaction \c sixferm,
	and the color-superconducting coupling, \c fourgap. If 0.0 is
	given for any of the values, then the default is used (\f$
	\Lambda=602.3/(\hbar c), G=1.835/\Lambda^2, K=12.36/\Lambda^5
	\f$).
	
	If the four-fermion coupling that produces a gap is not
	specified, it is automatically set to 3/4 G, which is the
	value obtained from the Fierz transformation.
	
	The value of the shift in the bag constant nambujl_eos::B0 is
	automatically calculated to ensure that the vacuum has zero
	energy density and zero pressure. The functions set_quarks()
	and set_thermo() must be used before hand to specify the \ref
	quark and \ref thermo objects.
	
    */
    virtual int set_parameters(double lambda=0.0, double fourferm=0.0, 
			       double sixferm=0.0, double fourgap=0.0);
    
    /** \brief Calculate the EOS
	
	Calculate the EOS from the quark condensates in \c u.qq, \c
	d.qq and \c s.qq. Return the mass gap equations in \c qq1, \c
	qq2, \c qq3, and the normal gap equations in \c gap1, \c gap2,
	and \c gap3.
	
	Using \c fromqq=false as in nambujl_eos and nambujl_eos does not
	work here and will return an error. Also, the quarks must be
	set through quark_eos::quark_set() before use.
	
	If all of the gaps are less than gap_limit, then the
	nambujl_eos::calc_temp_p() is used, and \c gap1, \c gap2, and
	\c gap3 are set to equal \c u.del, \c d.del, and \c s.del,
	respectively.

	\todo It surprises me that n3 is not -res[11]. Is there 
	a sign error in the color densities?
    */
    virtual int calc_eq_temp_p(quark &u, quark &d, quark &s, double &qq1, 
			       double &qq2, double &qq3, double &gap1, 
			       double &gap2, double &gap3, double mu3, 
			       double mu8, double &n3, double &n8, 
			       thermo &qb, double temper);
    
    /// Check the derivatives specified by eigenvalues()
    virtual int test_derivatives(double lmom, double mu3, double mu8,
				 test_mgr &t);
    
    /** \brief Calculate the energy eigenvalues as a function of the
	momentum
	
	Given the momentum \c mom, and the chemical potentials
	associated with the third and eighth gluons (\c mu3 and \c
	mu8), the energy eigenvalues are computed in egv[0]
	... egv[35].
    */
    virtual int eigenvalues(double lmom, double mu3, double mu8, 
			    double egv[36], double dedmuu[36], 
			    double dedmud[36], double dedmus[36], 
			    double dedmu[36], double dedmd[36], 
			    double dedms[36], double dedu[36], 
			    double dedd[36], double deds[36], 
			    double dedmu3[36], double dedmu8[36]);
    
    /// Set the routine for solving quartics
    int set_quartic(quartic_real_coeff &q) { 
      quartic=&q; 
      return 0;
    }
    
    /// The equal mass threshold
    double eq_limit;
    
    /// Set to true to test the integration (default false)
    bool integ_test;
    
    /// Test the integration routines
    int test_integration(test_mgr &t);

    //
    //gsl_quadratic_real_coeff quad;

    /** \brief The default quartic routine
	
	Slightly better accuracy (with slower execution times) can be
	achieved using \ref poly_real_coeff_gsl which polishes the
	roots of the quartics. For example

	\code
	cfl_njl_eos cfl;
	poly_real_coeff_gsl gp;
	cfl.set_quartic(gp);
	\endcode
     */
    quartic_real_coeff_cern def_quartic;
    //@}
    
    /** \brief Test the routine to compute the eigenvalues of 
	non-superfluid fermions
    */
    int test_normal_eigenvalues(test_mgr &t);

    /** \brief Test the routine to compute the eigenvalues of 
	superfluid fermions
    */
    int test_gapped_eigenvalues(test_mgr &t);
    
    /** \brief Smallest allowable gap (default 0.0)
	
	If any of the gaps are below this value, then it is assumed
	that they are zero and the equation of state is simplified
	accordingly. If all of the gaps are less than gap_limit, then
	the results from nambujl_eos are used in
	calc_eq_temp_p(), calc_temp_p() and thd_potential().
    */
    double gap_limit;
    
    /** \brief If this is true, then finite temperature
	corrections are ignored (default false)
	
	This implements some simplifications in the 
	momentum integration that are not possible at finite temperature.
    */
    bool zerot;
    
    /** \brief Use a fixed quark mass and ignore the quark
	condensates
    */
    bool fixed_mass;

    /// If true, then ensure color neutrality
    bool color_neut;
    
    /** \brief Diquark coupling constant (default 3 G/4)
	
	The default value is the one derived from a Fierz
	transformation. (\ref Buballa04)
    */
    double GD;

    /** \brief The absolute precision for the integration
	(default \f$ 10^{-4} \f$ )
	
	This is analogous to gsl_inte::epsabs
    */
    double inte_epsabs;

    /** \brief The relative precision for the integration 
	(default \f$ 10^{-4} \f$ )

	This is analogous to gsl_inte::epsrel
    */
    double inte_epsrel;
    
    /** \brief The number of points used in the last integration
	(default 0)

	This returns 21, 43, or 87 depending on the number of function
	evaluations needed to obtain the desired precision. If it
	the routine failes to obtain the desired precision, then
	this variable is set to 88.
    */
    size_t inte_npoints;

    /// Return string denoting type ("cfl_njl_eos")
    virtual const char *type() { return "cfl_njl_eos"; };
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /** \brief The integrands
	
	- res[0] is the thermodynamic potential, \f$ \Omega \f$
	- res[1] is \f$ d -\Omega / d T \f$
	- res[2] is \f$ d \Omega / d \mu_u \f$
	- res[3] is \f$ d \Omega / d \mu_d \f$
	- res[4] is \f$ d \Omega / d \mu_s \f$
	- res[5] is \f$ d \Omega / d m_u \f$
	- res[6] is \f$ d \Omega / d m_d \f$
	- res[7] is \f$ d \Omega / d m_s \f$
	- res[8] is \f$ d \Omega / d \Delta_{ds} \f$
	- res[9] is \f$ d \Omega / d \Delta_{us} \f$
	- res[10] is \f$ d \Omega / d \Delta_{ud} \f$
	- res[11] is \f$ d \Omega / d \mu_3 \f$
	- res[12] is \f$ d \Omega / d \mu_8 \f$
    */
    virtual int integrands(double p, double res[]);
    
    /// Compute ungapped eigenvalues and the appropriate derivatives
    int normal_eigenvalues(double m, double lmom, double mu, 
			   double lam[2], double dldmu[2], 
			   double dldm[2]);

    /** \brief Treat the simply gapped quarks in all cases gracefully

	This function uses the quarks \c q1 and \c q2 to construct the
	eigenvalues of the inverse propagator, properly handling the
	either zero or finite quark mass and either zero or finite
	quark gaps. In the case of finite quark mass and finite quark
	gaps, the quartic solver is used.
	
	The chemical potentials are separated so we can add the 
	color chemical potentials to the quark chemical potentials
	if necessary.

	This function is used by eigenvalues(). It does not work for
	the "ur-dg-sb" set of quarks which are paired in a non-trivial
	way.

	\todo In the code, the equal mass case seems to be commented
	out. Why?
    */
    int gapped_eigenvalues(double m1, double m2, double lmom,
			   double mu1, double mu2, double tdelta,
			   double lam[4], double dldmu1[4], 
			   double dldmu2[4], double dldm1[4],
			   double dldm2[4], double dldg[4]);    
    
    /// Temperature
    double temper;
    
    /// 3rd gluon chemical potential
    double smu3;
    
    /// 8th gluon chemical potential
    double smu8;

    /// \name Numerical methods
    //@{
    /// The routine to solve quartics
    quartic_real_coeff *quartic;
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<std::complex<double> > 
      ubvector_complex;
    typedef boost::numeric::ublas::matrix<std::complex<double> > 
      ubmatrix_complex;

    /// \name For computing eigenvalues
    //@{
    /// Inverse propagator matrix
    gsl_matrix_complex *iprop;
    
    /// The eigenvectors
    gsl_matrix_complex *eivec;

    /// The derivative of the inverse propagator wrt the ds gap
    ubmatrix_complex dipdgapu;
    /// The derivative of the inverse propagator wrt the us gap
    ubmatrix_complex dipdgapd;
    /// The derivative of the inverse propagator wrt the ud gap
    ubmatrix_complex dipdgaps;
    
    /// The eigenvalues
    gsl_vector *eval;

    /// Workspace for eigenvalue computation
    gsl_eigen_hermv_workspace *w;
    //@}

    /// \name For the integration
    //@{
    /// The error scaling function for integ_err
    double rescale_error(double err, double result_abs, 
			 double result_asc);
    
    /** \brief A new version of inte_qng_gsl to integrate several
	functions at the same time
    */
    int integ_err(double a, double b, const size_t nr,
		  ubvector &res, double &err2);
    //@}

  private:
    
    cfl_njl_eos(const cfl_njl_eos &);
    cfl_njl_eos& operator=(const cfl_njl_eos&);
    
#endif
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
