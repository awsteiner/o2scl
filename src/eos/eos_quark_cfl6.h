/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file eos_quark_cfl6.h
    \brief File defining \ref o2scl::eos_quark_cfl6
*/
#ifndef CFL6_EOS_H
#define CFL6_EOS_H

#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/eos_quark_cfl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief An EOS like \ref eos_quark_cfl but 
      with a color-superconducting 't Hooft interaction

      Beginning with the Lagrangian:
      \f[
      {\cal L} = {\cal L}_{Dirac} + {\cal L}_{NJL} + 
      {\cal L}_{'t Hooft} + {\cal L}_{SC} + {\cal L}_{SC6} 
      \f]
      \f[
      {\cal L}_{Dirac} = {\bar q} \left( i \partial -m - 
      \mu \gamma^0 \right) q
      \f]
      \f[
      {\cal L}_{NJL} = G_S \sum_{a=0}^8 
      \left[ \left( {\bar q} \lambda^a q \right)^2
      - \left( {\bar q} \lambda^a \gamma^5 q \right)^2 \right]
      \f]
      \f[
      {\cal L}_{'t Hooft} = G_D \left[
      \mathrm{det}_f {\bar q} \left(1-\gamma^5 \right) q
      +\mathrm{det}_f {\bar q} \left(1+\gamma^5 \right) q
      \right]
      \f]
      \f[
      {\cal L}_{SC} = G_{DIQ}
      \left( {\bar q}_{i \alpha} i \gamma^5 
      \varepsilon^{i j k} \varepsilon^{\alpha \beta \gamma} 
      q^C_{j \beta} \right)
      \left( {\bar q}_{\ell \delta} i \gamma^5 
      \epsilon^{\ell m k} 
      \epsilon^{\delta \varepsilon \gamma}
      q^C_{m \varepsilon} \right)
      \f]
      \f[
      {\cal L}_{SC6} = K_D
      \left( {\bar q}_{i \alpha} i \gamma^5 
      \varepsilon^{i j k} \varepsilon^{\alpha \beta \gamma} 
      q^C_{j \beta} \right)
      \left( {\bar q}_{\ell \delta} i \gamma^5 
      \epsilon^{\ell m n} 
      \epsilon^{\delta \varepsilon \eta}
      q^C_{m \varepsilon} \right)
      \left( {\bar q}_{k \gamma} q_{n \eta} \right)
      \f]

      We can simplify the relevant terms in \f${\cal L}_{NJL}\f$:
      \f[
      {\cal L}_{NJL} = G_S \left[ 
      \left({\bar u} u\right)^2+
      \left({\bar d} d\right)^2+
      \left({\bar s} s\right)^2
      \right]
      \f]
      and in \f${\cal L}_{'t Hooft}\f$:
      \f[
      {\cal L}_{NJL} = G_D \left( 
      {\bar u} u {\bar d} d {\bar s} s
      \right)
      \f]

      Using the definition:
      \f[
      \Delta^{k \gamma} =  \left< {\bar q} i \gamma^5 
      \epsilon \epsilon q^C_{} \right>
      \f]
      and the ansatzes:
      \f[
      ({\bar q}_1 q_2) ({\bar q}_3 q_4) \rightarrow
      {\bar q}_1 q_2 \left< {\bar q}_3 q_4 \right>
      +{\bar q}_3 q_4 \left< {\bar q}_1 q_2 \right>
      -\left< {\bar q}_1 q_2 \right> \left< {\bar q}_3 q_4 \right>
      \f]
      \f[
      ({\bar q}_1 q_2) ({\bar q}_3 q_4) ({\bar q}_5 q_6) \rightarrow
      {\bar q}_1 q_2 \left< {\bar q}_3 q_4 \right> 
      \left< {\bar q}_5 q_6 \right>
      +{\bar q}_3 q_4 \left< {\bar q}_1 q_2 \right>
      \left< {\bar q}_5 q_6 \right>
      +{\bar q}_5 q_6 \left< {\bar q}_1 q_2 \right>
      \left< {\bar q}_3 q_4 \right>
      -2\left< {\bar q}_1 q_2 \right> \left< {\bar q}_3 q_4 \right>
      \left< {\bar q}_5 q_6 \right>
      \f]
      for the mean field approximation, we can rewrite the Lagrangian
      \f[
      {\cal L}_{NJL} = 2 G_S \left[ 
      \left( {\bar u} u \right) \left< {\bar u} u \right>
      +\left( {\bar d} d \right) \left< {\bar d} d \right>
      +\left( {\bar s} s \right) \left< {\bar s} s \right>
      - \left< {\bar u} u \right>^2
      - \left< {\bar d} d \right>^2
      - \left< {\bar s} s \right>^2
      \right]
      \f]
      \f[
      {\cal L}_{'t Hooft} = - 2 G_D \left[
      \left( {\bar u} u \right) \left< {\bar u} u \right>
      \left< {\bar s} s \right>
      + \left( {\bar d} d \right) \left< {\bar u} u \right>
      \left< {\bar s} s \right>
      + \left( {\bar s} s \right) \left< {\bar u} u \right>
      \left< {\bar d} d \right>
      - 2 \left< {\bar u} u \right>\left< {\bar d} d \right>
      \left< {\bar s} s \right>
      \right]
      \f]
      \f[
      {\cal L}_{SC} = G_{DIQ} \left[
      \Delta^{k \gamma}
      \left( {\bar q}_{\ell \delta} i \gamma^5 
      \epsilon^{\ell m k} 
      \epsilon^{\delta \varepsilon \gamma}
      q^C_{m \varepsilon} \right)
      + \left( {\bar q}_{i \alpha} i \gamma^5 
      \varepsilon^{i j k} \varepsilon^{\alpha \beta \gamma} 
      q^C_{j \beta} \right)
      \Delta^{k \gamma \dagger}
      - \Delta^{k \gamma}
      \Delta^{k \gamma \dagger}
      \right]
      \f]
      \f[
      {\cal L}_{SC6} = K_D \left[
      \left( {\bar q}_{m \varepsilon} q_{n \eta} \right)
      \Delta^{k \gamma} \Delta^{m \varepsilon \dagger}
      + \left( {\bar q}_{i \alpha} i \gamma^5 
      \varepsilon^{i j k} \varepsilon^{\alpha \beta \gamma} 
      q^C_{j \beta} \right)
      \Delta^{m \varepsilon \dagger} 
      \left< {\bar q}_{m \varepsilon} q_{n \eta} \right>
      \right]
      \f]
      \f[
      + K_D \left[\Delta^{k \gamma}
      \left( {\bar q}_{\ell \delta} i \gamma^5 
      \epsilon^{\ell m n} 
      \epsilon^{\delta \varepsilon \eta}
      q^C_{m \varepsilon} \right)
      \left< {\bar q}_{m \varepsilon} q_{n \eta} \right>
      -2 
      \Delta^{k \gamma} \Delta^{m \varepsilon \dagger}
      \left< {\bar q}_{m \varepsilon} q_{n \eta} \right>
      \right]
      \f]
      
      If we make the definition \f$ {\tilde \Delta} =
      2 G_{DIQ} \Delta \f$

      \hline
      <b>References:</b>

      Created for \ref Steiner05.
  */
  class eos_quark_cfl6 : public eos_quark_cfl {
  public:
  
    eos_quark_cfl6();

    virtual ~eos_quark_cfl6();
  
    /** \brief Calculate the EOS
	\nothing
      
	Calculate the EOS from the quark condensates. Return the mass
	gap equations in \c qq1, \c qq2, \c qq3, and the normal gap
	equations in \c gap1, \c gap2, and \c gap3.
	
	Using \c fromqq=true as in eos_quark_njl and nambujl_temp_eos
	does not work here and will return an error.
	
	If all of the gaps are less than gap_limit, then the
	nambujl_temp_eos::calc_temp_p() is used, and \c gap1, \c gap2,
	and \c gap3 are set to equal \c u.del, \c d.del, and \c s.del,
	respectively.
      
    */
    virtual int calc_eq_temp_p(quark &u, quark &d, quark &s, 
			    double &qq1, double &qq2, double &qq3, 
			    double &gap1, double &gap2, double &gap3, 
			    double mu3, double mu8,
			    double &n3, double &n8, thermo &qb, 
			    double temper);

    /// The momentum integrands
    virtual int integrands(double p, double res[]);
    
    /// Check the derivatives specified by eigenvalues()
    virtual int test_derivatives(double lmom, double mu3, double mu8,
				 test_mgr &t);

    /** \brief Calculate the energy eigenvalues and their derivatives
	\nothing

	Given the momentum \c mom, and the chemical potentials
	associated with the third and eighth gluons (\c mu3 and \c mu8),
	this computes the eigenvalues of the inverse propagator and
	the assocated derivatives.
	
	Note that this is not the same as eos_quark_cfl::eigenvalues()
	which returns \c dedmu rather \c dedqqu.
    */
    virtual int eigenvalues6(double lmom, double mu3, double mu8, 
			     double egv[36], double dedmuu[36], 
			     double dedmud[36], double dedmus[36], 
			     double dedmu[36], double dedmd[36], 
			     double dedms[36], double dedu[36], 
			     double dedd[36], double deds[36], 
			     double dedmu3[36], double dedmu8[36]);

    /** \brief Construct the matrices, but don't solve the eigenvalue 
	problem
	\nothing
	
	This is used by check_derivatives() to make sure that the derivative
	entries are right.
     */
    virtual int make_matrices(double lmom, double mu3, double mu8, 
			      double egv[36], double dedmuu[36], 
			      double dedmud[36], double dedmus[36], 
			      double dedmu[36], double dedmd[36], 
			      double dedms[36], double dedu[36], 
			      double dedd[36], double deds[36], 
			      double dedmu3[36], double dedmu8[36]);
    
    /// The color superconducting 't Hooft coupling (default 0)
    double KD;

    /// Return string denoting type ("eos_quark_cfl6")
    virtual const char *type() { return "eos_quark_cfl6"; };

    /** \brief The absolute value below which the CSC 't Hooft coupling 
	is ignored(default \f$ 10^{-6} \f$)
    */
    double kdlimit;

    typedef boost::numeric::ublas::vector<std::complex<double> > 
      ubvector_complex;
    typedef boost::numeric::ublas::matrix<std::complex<double> > 
      ubmatrix_complex;

  protected:

#ifndef DOXYGEN_INTERNAL

    /// Set the quark effective masses from the gaps and the condensates
    int set_masses();
    
    /// The size of the matrix to be diagonalized
    static const int mat_size=36;
    
    /// Storage for the inverse propagator
    gsl_matrix_complex *iprop6;

    /// The eigenvectors
    gsl_matrix_complex *eivec6;

    /// The derivative wrt the ds gap
    ubmatrix_complex dipdgapu;

    /// The derivative wrt the us gap
    ubmatrix_complex dipdgapd;

    /// The derivative wrt the ud gap
    ubmatrix_complex dipdgaps;

    /// The derivative wrt the up quark condensate
    ubmatrix_complex dipdqqu;

    /// The derivative wrt the down quark condensate
    ubmatrix_complex dipdqqd;

    /// The derivative wrt the strange quark condensate
    ubmatrix_complex dipdqqs;

    /// Storage for the eigenvalues
    gsl_vector *eval6;

    /// GSL workspace for the eigenvalue computation
    gsl_eigen_hermv_workspace *w6;

  private:

    eos_quark_cfl6(const eos_quark_cfl6 &);
    eos_quark_cfl6& operator=(const eos_quark_cfl6&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
