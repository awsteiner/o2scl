!/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2020, Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include <o2scl/nucmass_ldrop_ext.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

nucmass_ldrop_ext::nucmass_ldrop_ext() {
  dim=3.0;
  use_ame=true;
  use_moller=true;
  extra_corr=true;

  o2scl_hdf::ame_load(ame,"16",true);
  o2scl_hdf::mnmsk_load(moller);
}

int nucmass_ldrop_ext::test_derivatives() {

  test_mgr t;
  t.set_output_level(2);

  use_ame=false;
  use_moller=false;

  // First basic test
  run_test(80,124,0.01,0.01,0.02,0.01,t);

  // Test with AME
  use_ame=true;
  use_moller=false;
  extra_corr=true;
  run_test(80,124,0.01,0.01,0.02,0.01,t);
  run_test(2,2,0.01,0.01,0.02,0.01,t);
  run_test(1,0,0.01,0.01,0.02,0.01,t);
  run_test(0,1,0.01,0.01,0.02,0.01,t);
  extra_corr=false;
  run_test(80,124,0.01,0.01,0.02,0.01,t);
  run_test(2,2,0.01,0.01,0.02,0.01,t);
  run_test(1,0,0.01,0.01,0.02,0.01,t);
  run_test(0,1,0.01,0.01,0.02,0.01,t);
  extra_corr=true;
  use_ame=false;
  use_moller=false;

  // Test with Moller
  use_ame=false;
  use_moller=true;
  extra_corr=true;
  run_test(80,124,0.01,0.01,0.02,0.01,t);
  extra_corr=false;
  run_test(80,124,0.01,0.01,0.02,0.01,t);
  extra_corr=true;
  use_ame=false;
  use_moller=false;

  // Test for chi -> 0
  run_test(80,124,0.01,0.01,2.0e-4,0.01,t);

  // Other misc. tests
  run_test(20,20,0.01,0.01,0.02,0.01,t);
  run_test(124,80,0.01,0.01,0.02,0.01,t);
  run_test(80,124,0.0,0.0,0.02,0.0,t);
  run_test(80,124,0.0,0.04,0.1,0.1,t);

  if (!t.report()) {
    exit(-1);
  }
  
  return 0;
}

int nucmass_ldrop_ext::run_test(double Z, double N, double npout,
                                double nnout, double chi, double T,
                                test_mgr &t) {
  double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;

  // Derivative object
  deriv_gsl<ldrop_be_deriv> gd;
  gd.h=1.0e-4;

  double der, dere;

  cout.setf(ios::showpos);
  
  cout << "Running test with Z,N,np,nn,chi,T: " << Z << " " << N << endl;
  cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;

  double be=drip_binding_energy_full_d(Z,N,npout,nnout,chi,T);
  d1=dbulk_dchi;
  d2=dcoul_dchi;
  d3=dexc_dnn;
  d4=dexc_dnp;
  d5=dshell_dnp;
  d6=dshell_dnn;
  d7=dexc_dchi; 
  d8=dsurf_dchi;
  d9=dshell_dchi;
  d10=dcoul_dnp;

  cout << "be: " << be << " " << bulk << " " << surf << " " << coul << endl;
  cout << "    " << pair << " " << exc << " " << shell << endl;
  cout << endl;

  cout << "chi: " << endl;

  // dbulk_dchi
  ldrop_be_deriv flc1(*this,Z,N,npout,nnout,chi,T,1,4);
  gd.h/=100.0;
  gd.deriv_err(chi,flc1,der,dere);
  gd.h*=100.0;
  cout << "dbulk_dchi: " 
       << der << " " << dere << " " << d1 << " " << fabs(der-d1) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d1),0.0,dere*2.0,"dbulk_dchi");
  }
  if (fabs(der-d1)>dere*2.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dbulk_dchi",o2scl::exc_efailed);
  }
  
  // dcoul_dchi
  ldrop_be_deriv flc2(*this,Z,N,npout,nnout,chi,T,3,4);
  gd.deriv_err(chi,flc2,der,dere);
  cout << "dcoul_dchi: " 
       << der << " " << dere << " " << d2 << " " << fabs(der-d2) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d2),0.0,dere*2.0,"dcoul_dchi");
  }
  if (fabs(der-d2)>dere*2.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dcoul_dchi",o2scl::exc_efailed);
  }

  // dexc_dchi
  ldrop_be_deriv flc3(*this,Z,N,npout,nnout,chi,T,6,4);
  gd.deriv_err(chi,flc3,der,dere);
  cout << "dexc_dchi : " 
       << der << " " << dere << " " << d7 << " " << fabs(der-d7) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d7),0.0,dere*2.0,"dexc_dchi");
  }
  if (fabs(der-d7)>dere*2.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dexc_dchi ",o2scl::exc_efailed);
  }

  // dsurf_dchi
  ldrop_be_deriv flc4(*this,Z,N,npout,nnout,chi,T,2,4);
  gd.deriv_err(chi,flc4,der,dere);
  cout << "dsurf_dchi: " 
       << der << " " << dere << " " << d8 << " " << fabs(der-d8) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d8),0.0,dere*2.0,"dsurf_dchi");
  }
  if (fabs(der-d8)>dere*2.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dsurf_dchi",o2scl::exc_efailed);
  }

  // dshell_dchi
  ldrop_be_deriv flc4xx(*this,Z,N,npout,nnout,chi,T,5,4);
  gd.deriv_err(chi,flc4xx,der,dere);
  cout << "dshel_dchi: " 
       << der << " " << dere << " " << d9 << " " << fabs(der-d9) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d9),0.0,dere*2.0,"dshel_dchi");
  }
  if (fabs(der-d9)>dere*2.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dshel_dchi",o2scl::exc_efailed);
  }
  
  // Total of chi derivatives
  ldrop_be_deriv flc4b(*this,Z,N,npout,nnout,chi,T,0,4);
  gd.deriv_err(chi,flc4b,der,dere);
  cout << "dfull_dchi: " 
       << der << " " << dere << " " << d1+d2+d7+d8+d9 << " " 
       << fabs(der-d1-d2-d7-d8-d9) << endl;
  if (der!=0.0) {
    t.test_rel(fabs(der-d1-d2-d7-d8-d9),0.0,dere*4.0,"dfull_dchi");
  }
  if (fabs(der-d1-d2-d7-d8-d9)>dere*4.0) {
    cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
    cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
    O2SCL_ERR("Failed at dfull_dchi",o2scl::exc_efailed);
  }

  cout << endl;

  if (nnout>0.0) {

    cout << "nn: " << endl;

    // dexc_dnn
    ldrop_be_deriv flc5(*this,Z,N,npout,nnout,chi,T,6,3);
    gd.h*=10.0;
    gd.deriv_err(nnout,flc5,der,dere);
    gd.h/=10.0;
    cout << "dexc_dnn  : " 
	 << der << " " << dere << " " << d3 << " " << fabs(der-d3) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d3),0.0,dere*2.0,"dexc_dnn");
    }
    if (fabs(der-d3)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dexc_dnn  ",o2scl::exc_efailed);
    }

    // dshell_dnn
    ldrop_be_deriv flc6(*this,Z,N,npout,nnout,chi,T,5,3);
    gd.deriv_err(nnout,flc6,der,dere);
    cout << "dshel_dnn : " 
	 << der << " " << dere << " " << d6 << " " << fabs(der-d6) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d6),0.0,dere*2.0,"dshel_dnn");
    }
    if (fabs(der-d6)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dshel_dnn ",o2scl::exc_efailed);
    }

    // Total of nn derivatives
    ldrop_be_deriv flc7(*this,Z,N,npout,nnout,chi,T,0,3);
    gd.deriv_err(nnout,flc7,der,dere);
    cout << "dfull_dnn : " << der << " " << dere << " " << d6+d3 << " " 
	 << fabs(der-d6-d3) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d3-d6),0.0,dere*2.0,"dfull_dnn");
    }
    if (fabs(der-d3-d6)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dfull_dnn ",o2scl::exc_efailed);
    }
  
    cout << endl;

  }

  if (npout>0.0) {

    cout << "np: " << endl;

    // dexc_dnp
    ldrop_be_deriv flc8(*this,Z,N,npout,nnout,chi,T,6,2);
    gd.deriv_err(npout,flc8,der,dere);
    cout << "dexc_dnp  : " 
	 << der << " " << dere << " " << d4 << " " << fabs(der-d4) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d4),0.0,dere*20.0,"dexc_dnp");
    }
    if (fabs(der-d4)>dere*20.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dexc_dnp.",o2scl::exc_efailed);
    }

    // dshell_dnp
    ldrop_be_deriv flc9(*this,Z,N,npout,nnout,chi,T,5,2);
    gd.deriv_err(npout,flc9,der,dere);
    cout << "dshel_dnp : " 
	 << der << " " << dere << " " << d5 << " " << fabs(der-d5) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d5),0.0,dere*2.0,"dshel_dnp");
    }
    if (fabs(der-d5)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dshel_dnp ",o2scl::exc_efailed);
    }

    // dcoul_dnp
    ldrop_be_deriv flc9b(*this,Z,N,npout,nnout,chi,T,3,2);
    gd.deriv_err(npout,flc9b,der,dere);
    cout << "dcoul_dnp : " << der << " " << dere << " "
	 << d10 << " " << fabs(der-d10) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d10),0.0,dere*2.0,"dcoul_dnp");
    }
    if (fabs(der-d10)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dcoul_dnp ",o2scl::exc_efailed);
    }

    // Total of np derivatives
    ldrop_be_deriv flc10(*this,Z,N,npout,nnout,chi,T,0,2);
    gd.deriv_err(npout,flc10,der,dere);
    cout << "dfull_dnp : " 
	 << der << " " << dere << " " << d5+d4+d10 << " " 
	 << fabs(der-d5-d4-d10) << endl;
    if (der!=0.0) {
      t.test_rel(fabs(der-d4-d5-d10),0.0,dere*2.0,"dfull_dnp");
    }
    if (fabs(der-d4-d5-d10)>dere*2.0) {
      cout << "Z,N,np,nn,chi,T: " << Z << " " << N << endl;
      cout << "\t" << npout << " " << nnout << " " << chi << " " << T << endl;
      O2SCL_ERR("Failed at dfull_dnp ",o2scl::exc_efailed);
    }

    cout << endl;

  }

  cout.unsetf(ios::showpos);

  return 0;
}

double nucmass_ldrop_ext::shell_energy_new(double Z, double N, 
                                           double pfact, double nfact, double T,
                                           double dpfact, double dnfact,
                                           double dndc, double dpdc,
                                           double &dsdp, double &dsdn,
                                           double &dsdc) {
  
  double Dn=0.0, Dz=0.0, nv=0.0, zv=0.0, nvbar=0.0, zvbar=0.0;
    
  // Determine the appropriate proton shell
  if (Z<2.0) {
    Dz=2.0;
    zv=Z;
  } else {
    bool done=false;
    for(size_t i=0;i<nshells-1 && done==false;i++) {
      if (Z>=shells[i] && Z<shells[i+1]) {
	Dz=shells[i+1]-shells[i];
	zv=Z-shells[i];
	zvbar=Dz-zv;
	done=true;
      }
    }
    if (Z>=shells[nshells-1]) {
      Dz=1.0;
      zv=0.0;
      zvbar=0.0;
      done=true;
    }
    if (done==false) {
      O2SCL_ERR2("Failed to do shell model correction ",
		 "in nucmass_ldrop_ext::shell_energy_new().",o2scl::exc_esanity);
    }
  }

  // Determine the appropriate neutron shell
  if (N<2.0) {
    Dn=2.0;
    nv=N;
  } else {
    bool done=false;
    for(size_t i=0;i<nshells-1 && done==false;i++) {
      if (N>=shells[i] && N<shells[i+1]) {
	Dn=shells[i+1]-shells[i];
	nv=N-shells[i];
	nvbar=Dn-nv;
	done=true;
      }
    }
    if (N>=shells[nshells-1]) {
      Dn=1.0;
      nv=0.0;
      nvbar=0.0;
      done=true;
    }
    if (done==false) {
      O2SCL_ERR2("Failed to do shell model correction ",
		 "in nucmass_ldrop_ext::shell_energy_new().",o2scl::exc_esanity);
    }
  }

  double S2=nv*nvbar/Dn*nfact+zv*zvbar/Dz*pfact;
  double S3=nv*nvbar*(nv-nvbar)/Dn*nfact+zv*zvbar*(zv-zvbar)/Dz*pfact;
  double Snp=nv*nvbar*zv*zvbar/Dn/Dz*nfact*pfact;

  // There's an extra minus sign here b/c the shell effects
  // where designed as an extra binding
  double ret=-(s_a1*S2+s_a2*S2*S2+s_a3*S3+s_anp*Snp);

  dsdn=-(s_a1*nv*nvbar/Dn*dnfact+s_a2*2.0*S2*nv*nvbar/Dn*dnfact+
	 s_a3*nv*nvbar*(nv-nvbar)/Dn*dnfact+
	 s_anp*nv*nvbar*zv*zvbar/Dn/Dz*pfact*dnfact);
  dsdp=-(s_a1*zv*zvbar/Dz*dpfact+s_a2*2.0*S2*zv*zvbar/Dz*dpfact+
	 s_a3*zv*zvbar*(zv-zvbar)/Dz*dpfact+
	 s_anp*nv*nvbar*zv*zvbar/Dn/Dz*nfact*dpfact);
  dsdc=-(s_a1*(nv*nvbar/Dn*dndc+zv*zvbar/Dz*dpdc)+
	 s_a2*2.0*S2*(nv*nvbar/Dn*dndc+zv*zvbar/Dz*dpdc)+
	 s_a3*(nv*nvbar*(nv-nvbar)/Dn*dndc+zv*zvbar*(zv-zvbar)/Dz*dpdc)+
	 s_anp*(nv*nvbar*zv*zvbar/Dn/Dz*(nfact*dpdc+pfact*dndc)));
  
  if (!finite(ret)) {
    cout << "Shell energy not finite: " << endl;
    cout << "Z,Zv,Dz,zvbar: " << Z << " " << zv << " " << Dz << " " 
	 << zvbar << endl;
    cout << "N,Nv,Dn,nvbar: " << N << " " << nv << " " << Dn << " " 
	 << nvbar << endl;
    cout << "S2,S3,Snp: " << S2 << " " << S3 << " " << Snp << endl;
    cout << "nfact,pfact: " << nfact << " " << pfact << endl;
    exit(-1);
  }

  return ret;
}

double nucmass_ldrop_ext::drip_binding_energy_full_d
(double Z, double N, double npout, double nnout, double chi, double T) {

  if (chi<0.0) {
    O2SCL_ERR2("Variable 'chi' is less than zero in ",
	       "nucmass_ldrop_ext::drip_binding_energy_full_d().",
               o2scl::exc_einval);
    return 0.0;
  }
    
  if (!finite(nnout)) {
    O2SCL_ERR2("Variable 'nnout' not finite in ",
	       "nucmass_ldrop_ext::drip_binding_energy_full_d().",
               o2scl::exc_einval);
    return 0.0;
  }

  int err;
  double ret=0.0, A=Z+N, nL;

  // Look for nuclei in the other databases
  bool ext_mass=false;
  if (use_ame && ame.is_included(((int)(Z+1.0e-8)),((int)(N+1.0e-8)))) {
    ext_mass=true;
    bulk=ame.binding_energy(((int)(Z+1.0e-8)),((int)(N+1.0e-8)));
  } else if (use_moller && 
	     moller.is_included(((int)(Z+1.0e-8)),((int)(N+1.0e-8)))) {
    ext_mass=true;
    bulk=moller.binding_energy(((int)(Z+1.0e-8)),((int)(N+1.0e-8)));
  }
  
  if (ext_mass) {

    // Use frdm_mass to determine densities and radii
    frdm.drip_binding_energy_d(Z,N,npout,nnout,chi);
    np=frdm.np;
    nn=frdm.nn;
    nL=np+nn;
    Rp=frdm.Rp;
    Rn=frdm.Rn;

    dRn_dchi=0.0;
    dRp_dchi=0.0;
    df_dchi=0.0;

    // Excluded volume part
    if (extra_corr && exc_volume && N>0.0 && Z>0.0 && 
	(nnout>0.0 || npout>0.0)) {

      n->n=nnout;
      p->n=npout;
      n->mu=n->m;
      p->mu=p->m;
      
      double vol=4.0/3.0*pi*pow(Rn,3.0);
      
      // Ensure a consistent guess for the chemical potentials
      n->nu=n->m*1.0001;
      p->nu=p->m*1.0001;
    
      if (T<=0.0) {
	err=heos->calc_e(*n,*p,th);
	exc=-vol*(th.ed-nnout*n->m-npout*p->m)*o2scl_const::hc_mev_fm;
      } else {
	err=heos->calc_temp_e(*n,*p,T,th);
	exc=-vol*(th.ed-T*th.en-nnout*n->m-npout*p->m)*
	  o2scl_const::hc_mev_fm;
      }
      
      dexc_dchi=0.0;
      dexc_dnn=-vol*(n->mu-n->m)*o2scl_const::hc_mev_fm;
      dexc_dnp=-vol*(p->mu-p->m)*o2scl_const::hc_mev_fm;
      
      if (err!=0) {
	O2SCL_ERR2("Hadronic eos failed in ",
		   "ldrop_mass_skin::drip_binding_energy_full_d().",
		   o2scl::exc_efailed);
      }

    } else {
      exc=0.0;
      dexc_dchi=0.0;
      dexc_dnn=0.0;
      dexc_dnp=0.0;
    }

    if (extra_corr==false || Rp==0.0 || Z==0.0 || N==0.0 || Rn==0.0) {

      coul=0.0;
      dcoul_dchi=0.0;
      dcoul_dnp=0.0;

    } else {

      // Add the finite-size part of the Coulomb energy
      double chip=chi*pow(Rp/Rn,3.0);
      double fdu=0.2*chip-0.6*cbrt(chip);
      
      double dchip_dchi, dfdu_dchip;
      if (chi==0.0) {
	// These limits are probably incorrect, and in general
	// dcoul_dchi diverges when chi is zero.
	dchip_dchi=0.0;
	dfdu_dchip=1.0;
      } else {
	dchip_dchi=chip/chi;
	dfdu_dchip=0.2-0.2/cbrt(chip)/cbrt(chip);
      }
      coul=A*coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
	o2scl_const::fine_structure_f<double>()*Rp*Rp*pow(fabs(np-npout),2.0)/nL*fdu;
      dcoul_dchi=A*coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
	o2scl_const::fine_structure_f<double>()*Rp*Rp*pow(fabs(np-npout),2.0)/nL*dfdu_dchip*
	dchip_dchi;

      if (!finite(dcoul_dchi)) {
	dcoul_dchi=0.0;
	dcoul_dnp=0.0;
	return 1.0e3;
	/*
	  cout << "Point 1:" << endl;
	  cout << "extra_corr: " << extra_corr << endl;
	  cout << "Z,N: " << Z << " " << N << endl;
	  cout << "A,Rp,np,npout: " << A << " " << Rp << " " << np << " "
	  << npout << endl;
	  cout << "nL,dfdu_dchip,dchip_dchi: " << nL << " " 
	  << dfdu_dchip << " " << dchip_dchi << endl;
	  O2SCL_ERR2("Variable dcoul_dchi not finite in ",
	  "nucmass_ldrop_ext.",o2scl::exc_efailed);
	*/
      }

      dcoul_dnp=-4.0*A*coul_coeff*o2scl_const::pi*o2scl_const::hc_mev_fm*
	o2scl_const::fine_structure_f<double>()*Rp*Rp*fabs(np-npout)/nL*fdu;

    }

    // Set surface, pairing, and shell energies to zero
    surf=0.0;
    pair=0.0;
    shell=0.0;

    dbulk_dchi=0.0;
    dsurf_dchi=0.0;
    dshell_dchi=0.0;
    dshell_dnp=0.0;
    dshell_dnn=0.0;
    
    ret=bulk+coul+exc;

    return ret;
  }

  // Force inc_rest_mass to true
  n->inc_rest_mass=true;
  p->inc_rest_mass=true;
    
  double I=1.0-2.0*Z/A;
  double X=Z/A;
  double delta=I*doi;

  // Determine the inner densities
  nL=n0+n1*I*I;
  double fchi=hd_coeff*(1.0-exp(hd_exp*chi))/(1.0-exp(hd_exp));
  np=nL*(1.0-delta)/2.0-fchi;
  nn=nL*(1.0+delta)/2.0+fchi;
  // Try this instead:
  // 
  // double g2=(1.0-exp(hd_exp*(chi-1.0)))/(1.0-exp(-hd_exp));
  // np=nL*(1.0-delta)/2.0*g2;
  // nn=nL*(1.0+delta)/2.0/g2;
  // 
  // and maybe increase hd_exp to something like 8.0
  // and rewrite the derivatives. 

  df_dchi=-hd_coeff/(1.0-exp(hd_exp))*hd_exp*exp(hd_exp*chi);
  double dnn_dchi=df_dchi;
  double dnp_dchi=-df_dchi;

  if (!finite(nn) || !finite(np)) {
    cout << "Z,N,nL,I: " << Z << " " << N << " " << nL << " " << I << endl;
    cout << "delta, np, nn, chi: " << delta << " " << np << " " 
	 << nn << " " << chi << endl;
    cout << "fchi, exp, coeff: " << fchi << " " << hd_exp << " "
	 << hd_coeff << endl;
    O2SCL_ERR2("Neutron or proton density not finite in ",
	       "nucmass_ldrop_ext::drip_binding_energy_full_d().",
               o2scl::exc_efailed);
    return 0.0;
  }
      
  // Determine radii
    
  Rp=cbrt(3.0*Z/np/4.0/o2scl_const::pi);
  Rn=cbrt(3.0*N/nn/4.0/o2scl_const::pi);

  dRn_dchi=-Rn/nn/3.0*dnn_dchi;
  dRp_dchi=-Rp/np/3.0*dnp_dchi;
    
  if (new_skin_mode==true) {

    cout << "Unsupported." << endl;
    exit(-1);

  } else {

    // Bulk part of the free energy 
    n->n=nn;
    p->n=np;
    n->mu=n->m;
    p->mu=p->m;

    // Ensure a consistent guess for the chemical potentials
    n->nu=n->m*1.0001;
    p->nu=p->m*1.0001;
    
    if (T<=0.0) {
      err=heos->calc_e(*n,*p,th);
      bulk=A*(th.ed-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm;
    } else {
      err=heos->calc_temp_e(*n,*p,T,th);
      bulk=A*(th.ed-T*th.en-nn*n->m-np*p->m)/nL*o2scl_const::hc_mev_fm;
    }
    
    double dbulk_dnn=A*(n->mu-n->m)/nL*o2scl_const::hc_mev_fm-bulk/nL;
    double dbulk_dnp=A*(p->mu-p->m)/nL*o2scl_const::hc_mev_fm-bulk/nL;
    dbulk_dchi=dbulk_dnn*dnn_dchi+dbulk_dnp*dnp_dchi;

    if (err!=0) {
      O2SCL_ERR2("Hadronic eos failed in ",
		 "ldrop_mass_skin::drip_binding_energy_full_d().",
		 o2scl::exc_efailed);
    }

    if (exc_volume && (nnout>0.0 || npout>0.0)) {
      n->n=nnout;
      p->n=npout;
      n->mu=n->m;
      p->mu=p->m;

      double vol=4.0/3.0*pi*pow(Rn,3.0);

      // Ensure a consistent guess for the chemical potentials
      n->nu=n->m*1.0001;
      p->nu=p->m*1.0001;
    
      if (T<=0.0) {
	err=heos->calc_e(*n,*p,th);
	exc=-vol*(th.ed-nnout*n->m-npout*p->m)*
	  o2scl_const::hc_mev_fm;
	dexc_dchi=-4.0*pi*Rn*Rn*(th.ed-nnout*n->m-npout*p->m)*
	  o2scl_const::hc_mev_fm*dRn_dchi;
      } else {
	err=heos->calc_temp_e(*n,*p,T,th);
	exc=-vol*(th.ed-T*th.en-nnout*n->m-npout*p->m)*
	  o2scl_const::hc_mev_fm;
	dexc_dchi=-4.0*pi*Rn*Rn*(th.ed-T*th.en-nnout*n->m-npout*p->m)*
	  o2scl_const::hc_mev_fm*dRn_dchi;
      }
      
      dexc_dnn=-vol*(n->mu-n->m)*o2scl_const::hc_mev_fm;
      dexc_dnp=-vol*(p->mu-p->m)*o2scl_const::hc_mev_fm;
	
      if (err!=0) {
	O2SCL_ERR2("Hadronic eos failed in ",
		   "ldrop_mass_skin::drip_binding_energy_full_d().",
		   o2scl::exc_efailed);
      }
    } else {
      exc=0.0;
      dexc_dchi=0.0;
      dexc_dnn=0.0;
      dexc_dnp=0.0;
    }

  }

  ret+=bulk;
  ret+=exc;
 
  // Determine surface energy
    
  //double surf_fact=(pow((nn-nnout)/nn,2.0)+
  //pow((np-npout)/np,2.0))/2.0;
  double surf_fact=1.0;
      
  if (full_surface) {

    double x=np/(nn+np);
    double x3=x*x*x;
    double omx=1.0-x;
    double omx3=omx*omx*omx;
    double bcoeff;
    if (ss==0.0) bcoeff=-16.0+96.0*surften/0.5;
    else bcoeff=-16.0+96.0/ss;
    double bfun=(16.0+bcoeff)/(1.0/x3+bcoeff+1.0/omx3);

    //double y=0.5-x;
    //double y2=y*y, y4=y2*y2;
    //double a=a0+a2*y2+a4*y4;
    //double arg=1.0-3.313*y2-7.362*y4;
    //double Tc=Tchalf*sqrt(1.0-3.313*y2-7.362*y4);
    
    double dbfun_dnp=bfun*(3.0*pow(nn+np,3.0)*(pow(np,-4.0)-pow(nn,-4.0))/
			   (bcoeff+pow(nn+np,3.0)*(pow(nn,-3.0)+
						   pow(np,-3.0))));
    
    /*
      if (T<Tc) {
      bfun*=pow((1.0-T*T/Tc/Tc)/(1.0+a*T*T/Tc/Tc),pp);
      } else {
      bfun=0.0;
      dbfun_dnp=0.0;
      }
    */
      
    double Rsurf=cbrt(3.0*A/4.0/o2scl_const::pi/nL);

    // To correct if nuclei are "inside-out"
    //if (Rsurf>Rws-Rsurf) Rsurf=Rws-Rsurf;

    // Equivalent but slightly less obtuse way of writing it
    surf=surften*4.0*pi*Rsurf*Rsurf*surf_fact*bfun;
    
    dsurf_dchi=surften*4.0*pi*Rsurf*Rsurf*surf_fact*dbfun_dnp*dnp_dchi;
    
  } else {
      
    double c4=0.0;//ss-1.0+0.5;
      
    surf=A*surften*(1.0-ss*delta*delta+c4*delta*delta*delta*delta)*
      cbrt(36.0*o2scl_const::pi*nL)/nL/cbrt(A)*surf_fact;

    dsurf_dchi=0.0;
      
  }
  
  //surf*=dim/3.0;
      
  ret+=surf;
      
  // Add Coulomb energy

  double chip=chi*pow(Rp/Rn,3.0);
  
  // Before I had chip=chi if Rn<=Rp, but this causes problems with
  // the derivatives and I'm not sure it's necessary anyway. The
  // variable chip is only used for the Coulomb energy.

  //if (Rn>Rp) {
  //chip=
  //} else {
  //chip=chi;
  //}

  double dfdu_dchip, dchip_dchi, fdu;

  // This mass formula only works for d=3, but the full
  // Coulomb energy is given here for reference
  fdu=(2.0/(dim-2.0)*(1.0-0.5*dim*pow(chip,1.0-2.0/dim))+chip)/
    (dim+2.0);
  
  if (chi==0.0) {
    // These limits are probably incorrect, and in general
    // dcoul_dchi diverges when chi is zero.
    dchip_dchi=0.0;
    dfdu_dchip=1.0;
  } else {
    dchip_dchi=chip/chi+chip*3.0*(dRp_dchi/Rp-dRn_dchi/Rn);
    dfdu_dchip=(1-pow(chip,-2.0/dim))/(2.0+dim);    
  }

  coul=A*coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*Rp*pow(fabs(np-npout),2.0)/nL*fdu;

  dcoul_dchi=A*coul_coeff*4.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*pow(fabs(np-npout),2.0)/nL*fdu*dRp_dchi
    +A*coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*Rp*pow(fabs(np-npout),2.0)/nL*dfdu_dchip*
    dchip_dchi+A*coul_coeff*4.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*Rp*fabs(np-npout)/nL*fdu*dnp_dchi;

  if (!finite(dcoul_dchi)) {
    dcoul_dchi=0.0;
    dcoul_dnp=0.0;
    return 1.0e3;
    /*
      cout << "Point 2: " << endl;
      cout << "Z,N,Rn,chip: " << Z << " " << N << " " << Rn << " " 
      << chip << endl;
      cout << "Rn,fdu,dRp_dchi: " << Rn << " " << fdu << " " 
      << dRp_dchi << endl;
      cout << "Rp,np,npout,nL: " << Rp << " " << np << " "
      << npout << " " << nL << endl;
      cout << "chi,dfdu_dchip,dchip_dchi: " << chi << " " 
      << dfdu_dchip << " " << dchip_dchi << endl;
      O2SCL_ERR2("Variable dcoul_dchi not finite in ",
      "nucmass_ldrop_ext.",o2scl::exc_efailed);
    */
  }

  dcoul_dnp=-4.0*A*coul_coeff*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rp*Rp*fabs(np-npout)/nL*fdu;
    
  /*
    Preliminary code for d!=3

    if (chip>0.5) 
    chip=1.0-chip;
    double Rcoul=Rws-Rp;
    fdu=(2.0/(dim-2.0)*(1.0-0.5*dim*pow(chip,1.0-2.0/dim))+chip)/
    (dim+2.0);
    coul=A*coul_coeff*2.0*o2scl_const::pi*o2scl_const::hc_mev_fm*
    o2scl_const::fine_structure_f<double>()*Rcoul*Rcoul*pow(fabs(np-npout),2.0)/nL*fdu;
  */
    
  if (!finite(coul)) {
    coul=1.0e10;
  }

  ret+=coul;

  // Pairing corrections
    
  pair=-A*Epair*(cos(Z*o2scl_const::pi)+cos(N*o2scl_const::pi))/
    2.0/pow(A,4.0/3.0);

  ret+=pair;

  if (inc_shell==false) return ret;
    
  // Make sure N and Z are integer-like
  double N0=floor(N+1.0e-12);
  double Z0=floor(Z+1.0e-12);
  
  // In-medium corrections to the shell effects
  double pfact=pow((np-npout)/np,2.0);
  double nfact=pow((nn-nnout)/nn,2.0);
  double dpfact=-2.0/np/np*(np-npout);
  double dnfact=-2.0/nn/nn*(nn-nnout);
  double dndc=2.0*(1.0-nnout/nn)*nnout/nn/nn*dnn_dchi;
  double dpdc=2.0*(1.0-npout/np)*npout/np/np*dnp_dchi;

  if (!finite(nfact)) {
    cout << "nfact not finite." << endl;
    cout << "nn,nnout: " << nn << " " << nnout << endl;
    exit(-1);
  }
  
  shell=shell_energy_new(Z0,N0,pfact,nfact,0.0,dpfact,dnfact,dndc,dpdc,
                         dshell_dnp,dshell_dnn,dshell_dchi);
  
  ret+=shell;

  return ret;
}

double nucmass_ldrop_ext::nucleus_be(int Z, int N, double npout, double nnout, 
                                     double T, double ne, double &Rws, double &chi) {
    
  if (!finite(nnout)) {
    cout << "Z,N,npout,nnout,T,ne: ";
    cout << Z << " " << N << " " << npout << "\n\t" << nnout << " "
	 << T << " " << ne << endl;
    O2SCL_ERR2("Variable 'nnout' not ",
	       "finite in nucleus_be().",o2scl::exc_efailed);
  }
    
  // The solver is more elegant but iteration appears faster,
  // especially in full_eq() and full_eq2()
  bool iterate=true, success=false;

  // Need to set be to an initial value because of the comparison
  // in the for loop below for iteration
  double be=0.0;
    
  // Old iteration method
  if (iterate) {
      
    // Compute the binding energy with an infinite cell size
    double be_last=drip_binding_energy_full_d(Z,N,npout,nnout,0.0,T)/
      hc_mev_fm;

    for(size_t i=0;i<30 && (i==0 || fabs((be_last-be)/be)>1.0e-10);i++) {
	
      be_last=be;
      
      // Compute the WS cell size
      Rws=cbrt((3.0*Z/4.0/pi-Rp*Rp*Rp*npout)/(ne-npout));
	
      // Compute the volume taken up by the nucleus in each cell, chi
      chi=pow(Rn/Rws,3.0);
	
      // Now, recompute the binding energy (in fm^-1) with the new
      // WS cell size
      be=drip_binding_energy_full_d(Z,N,npout,nnout,chi,T)/hc_mev_fm;
	
    }

    if (fabs((be_last-be)/be)<1.0e-6) {
      success=true;
    }
  }
    
  // Solver method
  if (success==false) {
    
    if (ne<=0.0) {

      chi=0.0;

    } else {
	
      ldrop_solve_chi fsp(this,Z,N,npout,nnout,T,ne);
      
      // Get a good guess for chi
      drip_binding_energy_full_d(Z,N,npout,nnout,0.0,T);
      Rws=cbrt((3.0*Z/4.0/pi-Rp*Rp*Rp*npout)/(ne-npout));
      chi=pow(Rn/Rws,3.0);
	
      // Find a bracketing range
      double ul=chi*1.01;
      chi/=1.01;
      int it=0;
      if (ul>1.0) ul=1.0;
      while (it<100 && fsp(ul)*fsp(chi)>0.0) {
	ul*=1.5;
	chi/=1.5;
	it++;
	if (ul>1.0) ul=1.0;
      }
      if (it==100) {
	O2SCL_ERR2("Computation of chi boundary failed in ",
		   "nucleus_be().",o2scl::exc_efailed);
      }
      
      grb.solve_bkt(chi,ul,fsp);
    }
	
    if (!finite(nnout)) {
      cout << "Z,N,npout,nnout,T,ne: ";
      cout << Z << " " << N << " " << npout << "\n\t" << nnout << " "
	   << T << " " << ne << endl;
      O2SCL_ERR2("Variable 'nnout' not finite in ",
		 "nucleus_be() [point 3].",o2scl::exc_efailed);
    }
    
    // Compute final binding energy
    be=drip_binding_energy_full_d(Z,N,npout,nnout,chi,T)/hc_mev_fm;
    
    // Compute the WS cell size
    Rws=cbrt((3.0*Z/4.0/pi-Rp*Rp*Rp*npout)/(ne-npout));
    
  }

  // If we converged to an unphysical nucleus
  if (Rws<Rn || Rws<Rp) {
    return 1.0e100;
  }

  return be;
}
  
double nucmass_ldrop_ext::nucleus_be_pasta
(int Z, int N, double npout, double nnout, double T, double ne, double &Rws, 
 double &chi) {
				   
    
  double be_min=1.0e100, dim_min=0.0;
    
  for(dim=3.0;dim>=0.9999;dim-=0.5) {
    double be=nucleus_be(Z,N,npout,nnout,T,ne,Rws,chi);
    if (be<be_min) {
      be_min=be;
      dim_min=dim;
    }
  }
    
  if (dim_min<0.01) {
    O2SCL_ERR("All dimensions failed in nucleus_be_pasta().",
	      o2scl::exc_efailed);
  }
    
  dim=dim_min;
  be_min=nucleus_be(Z,N,npout,nnout,T,ne,Rws,chi);

  return be_min;
}
  
nucmass_ldrop_ext::ldrop_solve_chi::ldrop_solve_chi
(nucmass_ldrop_ext *tp, int Z, int N, double np, double nn, double T, double ne) {
  tptr=tp;
  Z_=Z;
  N_=N;
  np_=np;
  nn_=nn;
  T_=T;
  ne_=ne;
}
  
double nucmass_ldrop_ext::ldrop_solve_chi::operator()(double chi) const {
  
  // Compute the nucleus
  double x=tptr->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi,T_)/
    o2scl_const::hc_mev_fm;
  
  // Compute the WS cell size
  double Rws=cbrt((3.0*Z_/4.0/o2scl_const::pi-
		   tptr->Rp*tptr->Rp*tptr->Rp*np_)/(ne_-np_));
  
  // Solve for chi
  return (chi-pow(tptr->Rn/Rws,3.0))/chi;
}

nucmass_ldrop_ext::ldrop_be_deriv::ldrop_be_deriv
(nucmass_ldrop_ext &ld, double Z, double N, double npout, double nnout,
 double chi, double T, size_t ix, size_t jx) {
  ldp=&ld;
  Z_=Z;
  N_=N;
  np_=npout;
  nn_=nnout;
  chi_=chi;
  T_=T;
  ix_=ix;
  jx_=jx;
}

double nucmass_ldrop_ext::ldrop_be_deriv::operator()(double x) {
  if (jx_==0) {
    Z_=x;
  } else if (jx_==1) {
    N_=x;
  } else if (jx_==2) {
    np_=x;
  } else if (jx_==3) {
    nn_=x;
  } else if (jx_==4) {
    chi_=x;
  } else {
    T_=x;
  }
  if (ix_==0) {
    return ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
  } else if (ix_==1) {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->bulk;
  } else if (ix_==2) {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->surf;
  } else if (ix_==3) {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->coul;
  } else if (ix_==4) {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->pair;
  } else if (ix_==5) {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->shell;
  } else {
    ldp->drip_binding_energy_full_d(Z_,N_,np_,nn_,chi_,T_);
    return ldp->exc;
  }
}

