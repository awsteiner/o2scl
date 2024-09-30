/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2018-2024, Xingfu Du and Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ───────────────────────────────────────────────────────────────────
*/
#ifndef EOS_HAD_VIRIAL_H
#define EOS_HAD_VIRIAL_H

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/linear_solver.h>
#include <o2scl/poly.h>

namespace o2scl {

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  /** \brief Virial EOS with spin (experimental)
   */
  class eos_had_virial_spin {
    
  protected:
    
  public:

    virtual int calc_temp_p_spin(fermion &n_up, fermion &n_dn,
                                 fermion &p_up, fermion &p_dn,
                                 double T, thermo &th) {
      
      double lam_n_up=sqrt(2.0*o2scl_const::pi/n_up.m/T);
      double lam_n_dn=sqrt(2.0*o2scl_const::pi/n_dn.m/T);
      double lam_p_up=sqrt(2.0*o2scl_const::pi/p_up.m/T);
      double lam_p_dn=sqrt(2.0*o2scl_const::pi/p_dn.m/T);
      
      double lam_n_up3=lam_n_up*lam_n_up*lam_n_up;
      double lam_n_dn3=lam_n_dn*lam_n_dn*lam_n_dn;
      double lam_p_up3=lam_p_up*lam_p_up*lam_p_up;
      double lam_p_dn3=lam_p_dn*lam_p_dn*lam_p_dn;

      double z_n_up;
      if (n_up.inc_rest_mass) {
        z_n_up=exp((n_up.mu-n_up.m)/T);
      } else {
        z_n_up=exp(n_up.mu/T);
      }
      double z_n_dn;
      if (n_dn.inc_rest_mass) {
        z_n_dn=exp((n_dn.mu-n_dn.m)/T);
      } else {
        z_n_dn=exp(n_dn.mu/T);
      }
      double z_p_up;
      if (p_up.inc_rest_mass) {
        z_p_up=exp((p_up.mu-p_up.m)/T);
      } else {
        z_p_up=exp(p_up.mu/T);
      }
      double z_p_dn;
      if (p_dn.inc_rest_mass) {
        z_p_dn=exp((p_dn.mu-p_dn.m)/T);
      } else {
        z_p_dn=exp(p_dn.mu/T);
      }

      double bn1, bn0, bpn1, bpn0;
      
      th.pr=T*(z_n_up/lam_n_up3+z_n_dn/lam_n_dn3+
               z_p_up/lam_p_up3+z_p_dn/lam_p_dn3)+
        bn1*(z_n_up*z_n_up+z_n_dn*z_n_dn+z_p_up*z_p_up+z_p_dn*z_p_dn)+
        bn0*(z_n_up*z_n_dn+z_p_up*z_p_dn)+
        2.0*bpn1*(z_n_up*z_p_up+z_n_dn*z_p_dn)+
        2.0*bpn0*(z_n_up*z_p_dn+z_n_dn*z_p_up);

      n_up.n=z_n_up/lam_n_up3+
        bn1/T*2.0*z_n_up*z_n_up+
        bn0/T*z_n_up*z_n_dn+
        2.0/T*bpn1*z_n_up*z_p_up+
        2.0/T*bpn0*z_n_up*z_p_dn;

      n_dn.n=z_n_dn/lam_n_dn3+
        bn1/T*2.0*z_n_dn*z_n_dn+
        bn0/T*z_n_dn*z_n_up+
        2.0/T*bpn1*z_n_dn*z_p_dn+
        2.0/T*bpn0*z_n_dn*z_p_up;

      p_up.n=z_p_up/lam_p_up3+
        bn1/T*2.0*z_p_up*z_p_up+
        bn0/T*z_p_up*z_p_dn+
        2.0/T*bpn1*z_p_up*z_n_up+
        2.0/T*bpn0*z_p_up*z_n_dn;

      p_dn.n=z_p_dn/lam_p_dn3+
        bn1/T*2.0*z_p_dn*z_p_dn+
        bn0/T*z_p_dn*z_p_up+
        2.0/T*bpn1*z_p_dn*z_n_dn+
        2.0/T*bpn0*z_p_dn*z_n_up;
      
      th.en=(z_n_up/lam_n_up3+z_n_dn/lam_n_dn3+
             z_p_up/lam_p_up3+z_p_dn/lam_p_dn3);
      
      double fr=n_up.mu*n_up.n+n_dn.mu*n_dn.n+
        p_up.mu*p_up.n+p_dn.mu*p_dn.n-th.pr;

      th.ed=fr+T*th.en;

      return 0;
    }
    
  };
  
  /** \brief Virial solver with derivatives
   */
  class eos_had_virial {
  
  protected:
  
    // Generic polynomial solver
    o2scl::poly_real_coeff_gsl<> quart;
  
    /// Storage for the four roots
    std::complex<double> res[4];

    /// Classical thermodynamics for \f$ \nu_n \f$ and \f$ \nu_p \f$
    classical_thermo cl;

  public:

    /// \name First derivatives of the fugacities
    //@{
    double dzndnn;
    double dzndnp;
    double dzpdnn;
    double dzpdnp;
    double dzndT;
    double dzpdT;
    //@}
  
    /// \name Second derivatives of the fugacities
    //@{
    double d2zndnn2;
    double d2zndnndnp;
    double d2zndnp2;
    double d2zpdnn2;
    double d2zpdnndnp;
    double d2zpdnp2;
    //@}

    /// \name Main functions
    //@{
    /** \brief Solve for the fugacity given the density
     */
    virtual void solve_fugacity(double nn, double np,
                                double lam_n, double lam_p,
                                double b_n, double b_pn,
                                double &zn, double &zp) {

      double npt=pow(lam_n,3)/2.0*np;
      double nnt=pow(lam_p,3)/2.0*nn;

      // At high densities or very low densities, just use the
      // non-interacting result
      //
      // AWS: 9/13/2020: I added the "|| nnt>1.0e5 || npt>1.0e5" option
      // later, and this might not be the best method
      // 
      if (nnt<5.0e-6 || npt<5.0e-6 || nnt>1.0e5 || npt>1.0e5) {
        zn=nnt;
        zp=npt;
        return;
      }

      zn=0.0;
      zp=0.0;
    
      double a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
      double b=-1+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
      double c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-nnt+npt-b_n*b_n*npt*2.0/b_pn/b_pn;
      double d=-b_n*npt/b_pn/b_pn+npt/2.0/b_pn;
      double e=b_n*npt*npt/2.0/b_pn/b_pn;
    
      quart.solve_rc(a,b,c,d,e,res[0],res[1],res[2],res[3]);
    
      std::vector<double> zp_list, zn_list;
      for(size_t k=0;k<4;k++) {
        if (res[k].imag()==0.0 && res[k].real()>0.0 && res[k].real()<500.0) {
          double r0, r1;
          gsl_poly_solve_quadratic(2.0*b_n,2.0*res[k].real()*
                                   b_pn+1.0,-nnt,&r0,&r1);
          if (r0>0.0 && r0<500.0) {
            if (r1>0.0 && r1<500.0) {
              O2SCL_ERR2("Unexpected pair of roots in ",
                         "eos_had_virial::solve_fugacity().",
                         o2scl::exc_einval);
            }
            std::cout << res[k] << "," << r0 << " ";
            zp_list.push_back(res[k].real());
            zn_list.push_back(r0);
          }
          if (r1>0.0 && r1<500.0) {
            zp_list.push_back(res[k].real());
            zn_list.push_back(r1);
          }
        }
      }
      if (zp_list.size()==1) {
        zp=zp_list[0];
        zn=zn_list[0];
      } else if (zp_list.size()==2) {
        double norm_0=zp_list[0]*zp_list[0]+zn_list[0]*zn_list[0];
        double norm_1=zp_list[1]*zp_list[1]+zn_list[1]*zn_list[1];
        if (norm_0<norm_1) {
          zp=zp_list[0];
          zn=zn_list[0];
        } else {
          zp=zp_list[1];
          zn=zn_list[1];
        }
      } else {
        std::cout << "eos_had_virial::solve_fugacity "
                  << "multiplicity problem:\n\t"
                  << "res0,res1: " << res[0] << " " << res[1] << "\n\t"
                  << "res2,res3: " << res[2] << " " << res[3] << std::endl;
        std::cout << "\tnn,np,lam_n,lam_p: " << nn << " " << np << " "
                  << lam_n << " " << lam_p << "\n\t"
                  << "nnt,npt,zp_list.size(): " << nnt << " " << npt << " "
                  << zp_list.size() << std::endl;
        O2SCL_ERR2("Unexpected root multiplicity in ",
                   "eos_had_virial::solve_fugacity().",o2scl::exc_einval);
      }
      return;
    }

    /** \brief Compute \f$ \nu_n \f$ and \f$ \nu_p \f$

        \note This function ignores the value of
        \ref fermion_tl::non_interacting and just sets the
        effective mass in \ref part_tl::ms equal to the
        bare mass \ref part_tl::m .
     */
    virtual void calc_nun_nup(fermion &n, fermion &p,
                              double T) {

      n.ms=n.m;
      p.ms=p.m;
      cl.calc_density(n,T);
      cl.calc_density(p,T);
      
      return;
    }
    
    /** \brief Compute the derivatives given the densities and the
        fugacities (i.e. after a call to \ref solve_fugacity())
    */
    virtual void calc_deriv(double nn, double np,
                            double lam_n, double lam_p,
                            double b_n, double b_pn,
                            double zn, double zp,
                            double dbndT, double dbpndT,
                            double dlamndT, double dlampdT) {
    
      double npt=pow(lam_n,3)/2.0*np;
      double nnt=pow(lam_p,3)/2.0*nn;
    
      // At high densities or very low densities, just use the
      // non-interacting result

      if (nnt<5.0e-6 || npt<5.0e-6) {
      
        dzndnn=nnt/nn;
        dzpdnp=npt/np;
        dzndnp=0.0;
        dzpdnn=0.0;
      
        d2zndnn2=0.0;
        d2zndnndnp=0.0;
        d2zndnp2=0.0;
        d2zpdnn2=0.0;
        d2zpdnndnp=0.0;
        d2zpdnp2=0.0;
        dzndT=1.5*lam_n*lam_n*nn*dlamndT;
        dzpdT=1.5*lam_p*lam_p*np*dlampdT;
      
        return;
      }
    
      dzndnn=-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)/
        (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
         (4*b_n*zp+2*b_pn*zn+1));
      dzndnp=b_pn*pow(lam_p,3)*zn/
        (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
         (4*b_n*zp+2*b_pn*zn+1));
      dzpdnn=b_pn*pow(lam_n,3)*zp/
        (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
         (4*b_n*zp+2*b_pn*zn+1));
      dzpdnp=-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)/
        (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
         (4*b_n*zp+2*b_pn*zn+1));
    
      double dzndnn_dzn=-b_pn*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzndnn_dzp=-2*b_n*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzndnp_dzn=b_pn*pow(lam_p,3)*zn*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2)+b_pn*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1));
      double dzndnp_dzp=b_pn*pow(lam_p,3)*zn*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnn_dzn=b_pn*pow(lam_n,3)*zp*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnn_dzp=b_pn*pow(lam_n,3)*zp*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2)+b_pn*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1));
      double dzpdnp_dzn=-2*b_n*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnp_dzp=-b_pn*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);

      d2zndnn2=dzndnn_dzn*dzndnn+dzndnn_dzp*dzpdnn;
      d2zndnndnp=dzndnn_dzn*dzndnp+dzndnn_dzp*dzpdnp;
      d2zndnp2=dzndnp_dzn*dzndnp+dzndnp_dzp*dzpdnp;
      d2zpdnn2=dzpdnn_dzn*dzndnn+dzpdnn_dzp*dzpdnn;
      d2zpdnndnp=dzpdnn_dzn*dzndnp+dzpdnn_dzp*dzpdnp;
      d2zpdnp2=dzpdnp_dzn*dzndnp+dzpdnp_dzp*dzpdnp;

      dzndT=(1.0/2.0)*(-16*b_n*dbndT*pow(zn,2)*zp-16*b_n*dbpndT*zn*pow(zp,2)+12*b_n*dlamndT*pow(lam_n,2)*nn*zp-8*b_pn*dbndT*pow(zn,3)+8*b_pn*dbndT*zn*pow(zp,2)+6*b_pn*dlamndT*pow(lam_n,2)*nn*zn-6*b_pn*dlampdT*pow(lam_p,2)*np*zn-4*dbndT*pow(zn,2)-4*dbpndT*zn*zp+3*dlamndT*pow(lam_n,2)*nn)/(16*pow(b_n,2)*zn*zp+8*b_n*b_pn*pow(zn,2)+8*b_n*b_pn*pow(zp,2)+4*b_n*zn+4*b_n*zp+2*b_pn*zn+2*b_pn*zp+1);
    
      dzpdT=(-b_pn*dzndT*zp+(3.0/4.0)*dlamndT*pow(lam_n,2)*nn-1.0/2.0*dzndT-zn*(2*b_n*dzndT+dbndT*zn+dbpndT*zp))/(b_pn*zn);
    
      return;
    }
    //@}

  };
  
}

#endif
