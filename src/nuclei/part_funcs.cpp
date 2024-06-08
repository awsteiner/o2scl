/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2022-2024, Andrew W. Steiner
  
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
#include <o2scl/part_funcs.h>
#include <o2scl/nucdist.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

part_funcs::part_funcs() : cu(o2scl_settings.get_convert_units()) {
  rtk_alpha=0.1337;
  rtk_beta=-0.06571;
  rtk_gamma=0.04884;
  spin_deg_mode=0;
}

int part_funcs::load_rt00(std::string fname, bool external) {

  if (!external) {
    std::string dir=o2scl::o2scl_settings.get_data_dir();
    fname=dir+"/pf_frdm_low.o2";
  }
  
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  hdf_input(hf,tab_rt00);
  hf.close();

  return 0;
}

int part_funcs::load_r03(std::string fname, bool external) {

  if (!external) {
    std::string dir=o2scl::o2scl_settings.get_data_dir();
    fname=dir+"/pf_frdm_high.o2";
  }
  
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  hdf_input(hf,tab_r03);
  hf.close();

  return 0;
}

int part_funcs::load_g08(std::string fname, bool external) {

  if (!external) {
    std::string dir=o2scl::o2scl_settings.get_data_dir();
    fname=dir+"/bruslib.o2";
  }
  
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  hdf_input(hf,tab_g08);
  hf.close();

  return 0;
}

int part_funcs::rt00(int Z, int N, double T_K, double &pf, double &TdpfdT) {

  int A=Z+N;
  interp_vec<vector<double>> itp(itp_linear);
  for(size_t i=0;i<tab_rt00.get_nlines();i++) {
    if (fabs(tab_rt00.get("Z",i)-Z)+fabs(tab_rt00.get("A",i)-A)<1.0e-4) {
      vector<double> x={0.01,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
        1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0};
      for(size_t j=0;j<x.size();j++) x[j]*=1.0e9;
      vector<double> y;
      for(size_t j=3;j<27;j++) {
        y.push_back(tab_rt00.get(tab_rt00.get_column_name(j),i));
      }
      itp.set(x.size(),x,y);
      pf=itp.eval(T_K);//,x.size(),x,y);
      TdpfdT=T_K*itp.deriv(T_K);//,x.size(),x,y);
      double g=tab_rt00.get("spin",i)*2.0+1.0;
      pf*=g;
      TdpfdT*=g;
      return 0;
    }
  }
  return 1;
}

int part_funcs::r03(int Z, int N, double T_K, double &pf, double &TdpfdT) {

  int A=Z+N;
  interp_vec<vector<double>> itp(itp_linear);
  for(size_t i=0;i<tab_r03.get_nlines();i++) {
    if (fabs(tab_r03.get("Z",i)-Z)+fabs(tab_r03.get("A",i)-A)<1.0e-4) {
      vector<double> x={12,14,16,18,20,22,24,26,28,30,35,40,45,50,
        55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,
        140,145,150,155,160,165,170,175,180,190,200,210,220,230,240,
        250,275};
      for(size_t j=0;j<x.size();j++) x[j]*=1.0e9;
      vector<double> y;
      for(size_t j=3;j<3+x.size();j++) {
        y.push_back(tab_r03.get(tab_r03.get_column_name(j),i));
      }
      itp.set(x.size(),x,y);
      pf=itp.eval(T_K);//,x.size(),x,y);
      TdpfdT=T_K*itp.deriv(T_K);//,x.size(),x,y);
      double g=tab_r03.get("J0",i)*2.0+1.0;
      pf*=g;
      TdpfdT*=g;
      return 0;
    }
  }
  return 1;
}

int part_funcs::ghk08(int Z, int N, double T_K, double &pf,
                      double &TdpfdT) {

  interp_vec<vector<double>> itp(itp_linear);
  for(size_t i=0;i<tab_g08.get_nlines();i++) {
    if (fabs(tab_g08.get("Z",i)-Z)+fabs(tab_g08.get("N",i)-N)<1.0e-4) {
      vector<double> x={0.001,0.005,0.010,0.050,0.100,0.150,0.200,0.250,
        0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.000,
        1.500,2.000,2.500,3.000,3.500,4.000,5.000,6.000,
        7.000,8.000,9.000,10.000};
      for(size_t j=0;j<x.size();j++) x[j]*=1.0e9;
      vector<double> y;
      for(size_t j=3;j<3+x.size();j++) {
        y.push_back(tab_g08.get(tab_g08.get_column_name(j),i));
      }
      itp.set(x.size(),x,y);
      pf=itp.eval(T_K);//,x.size(),x,y);
      TdpfdT=T_K*itp.deriv(T_K);//,x.size(),x,y);
      double g=tab_g08.get("spin",i)*2.0+1.0;
      pf*=g;
      TdpfdT*=g;
      return 0;
    }
  }
  return 1;
}

double part_funcs::delta_small_iand(double E, double T_MeV, double delta,
                                    double a) {
  
  //, int a_delta, int Z, int N,
  //double E_mic) {
  if (E<1.0e-200) return 0.0;

  //if (a_delta==1) {
  //a=(rtk_alpha*(Z+N)+rtk_beta*pow(((double)Z+N),2.0/3.0))*
  //(1.0+E_mic*(1.0-exp(-rtk_gamma*U)/U));
  //}
  
  // a is in 1/MeV, delta, E, and T_MeV are in MeV, so the
  // integrand is in 1/MeV
  double ret=sqrt(pi)/12.0*exp(2.0*sqrt(a*(E-delta)))/
    pow(a,1.0/4.0)/pow((E-delta),5.0/4.0)*exp(-E/T_MeV);
  
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    cout << exp(2.0*sqrt(a*(E-delta))) << " " << exp(-E/T_MeV) << " " 
	 << pow((E-delta),5.0/4.0) << " " << ret << endl;
    O2SCL_ERR2("Value of delta_small_iand is not finite ",
	       "in part_funcs::delta_small_iand().",o2scl::exc_efailed);
  }
  return ret;
}

double part_funcs::delta_small_iand_prime(double E, double T_MeV, double delta,
                                          double a) {
  
  if (E<1.0e-200) return 0.0;
  
  // a is in 1/MeV, delta, E, and T_MeV are in MeV, so the
  // integrand is in 1/MeV. 
  double ret=E/T_MeV*sqrt(pi)/12.0*exp(2.0*sqrt(a*(E-delta)))/
    pow(a,1.0/4.0)/pow((E-delta),5.0/4.0)*exp(-E/T_MeV);  
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    O2SCL_ERR2("Value of delta_small_iand_prime is not finite ",
	       "in part_funcs::delta_small_iand_prime().",
	       o2scl::exc_efailed);
  }
  return ret;
}
  
double part_funcs::delta_large_iand(double E, double T_MeV, double delta,
                                    double Tc, double C) {
  
  if (E<1.0e-200) return 0.0;
  // C is in 1/MeV so the integrand is in 1/MeV
  double ret=C*exp((E-delta)/Tc)*exp(-E/T_MeV); 
  if (!std::isfinite(ret)) {
    O2SCL_ERR2("Value of delta_large_iand is not finite ",
	       "in part_funcs::delta_large_iand().",o2scl::exc_efailed);
  }
  return ret;
} 

double part_funcs::delta_large_iand_prime(double E, double T_MeV, double delta,
                                          double Tc, double C) {
  if (E<1.0e-200) return 0.0;
  // C is in 1/MeV so the integrand is in 1/MeV
  double ret=E/T_MeV*C*exp((E-delta)/Tc)*exp(-E/T_MeV); 
  if (!std::isfinite(ret)) {
    O2SCL_ERR2("Value of delta_large_iand_prime is not finite ",
	       "in part_funcs::delta_large_iand_prime().",
	       o2scl::exc_efailed);
  }
  return ret;
} 

void part_funcs::compare_spin_deg() {

  vector<nucleus> dist;
  nucdist_set(dist,mnmsk);
  int comp1=0, comp2=0;
  
  for(size_t i=0;i<dist.size();i++) {
    int Z=dist[i].Z;
    int N=dist[i].N;

    double g_hfb=-1.0;
    if (hfb.is_included(Z,N)) {
      if (hfb.get_ZN(Z,N).Jexp<99) {
        g_hfb=2.0*hfb.get_ZN(Z,N).Jexp+1.0;
      } else {
        g_hfb=2.0*hfb.get_ZN(Z,N).Jth+1.0;
      }
    } else {
      if (Z%2==0 && N%2==0) {
        g_hfb=1.0;
      } else {
        g_hfb=2.0;
      }
    }

    double g_rt00=-1.0;
    for(size_t j=0;j<tab_rt00.get_nlines();j++) {
      if (fabs(Z-tab_rt00.get("Z",j))+
          fabs(Z+N-tab_rt00.get("A",j))<1.0e-3) {
        g_rt00=tab_rt00.get("spin",j)*2.0+1.0;
      }
    }

    double g_r03=-1.0;
    for(size_t j=0;j<tab_r03.get_nlines();j++) {
      if (fabs(Z-tab_r03.get("Z",j))+
          fabs(Z+N-tab_r03.get("A",j))<1.0e-3) {
        g_r03=tab_r03.get("J0",j)*2.0+1.0;
      }
    }

    bool found=false;
    if (g_rt00>=0.0 && g_hfb>=0.0 && fabs(g_rt00-g_hfb)>0.5) {
      comp1++;
      found=true;
    }
    if (g_rt00>=0.0 && g_r03>=0.0 && fabs(g_rt00-g_r03)>0.5) {
      comp2++;
      found=true;
    }
    if (false && found) {
      cout << Z << " " << N << " " << g_hfb << " "
           << g_rt00 << " " << g_r03 << endl;
    }
    
  }

  cout << dist.size() << " " << comp1 << " " << comp2 << endl;

  return;
}

double part_funcs::get_spin_deg(int Z, int N) {

  if (Z==2 && N==2) {
    return 1.0;
  } else if (Z==1 && N==1) {
    return 3.0;
  } else if (Z==1 && N==2) {
    return 2.0;
  } else if (Z==3 && N==1) {
    return 5.0;
  } else if (Z==2 && N==1) {
    return 2.0;
  }
  
  if (spin_deg_mode==0) {

    for(size_t j=0;j<tab_g08.get_nlines();j++) {
      if (fabs(Z-tab_g08.get("Z",j))+
          fabs(N-tab_g08.get("N",j))<1.0e-3) {
        return tab_g08.get("spin",j)*2.0+1.0;
      }
    }
    for(size_t j=0;j<tab_r03.get_nlines();j++) {
      if (fabs(Z-tab_r03.get("Z",j))+
          fabs(Z+N-tab_r03.get("A",j))<1.0e-3) {
        return tab_r03.get("J0",j)*2.0+1.0;
      }
    }
    
    for(size_t j=0;j<tab_rt00.get_nlines();j++) {
      if (fabs(Z-tab_rt00.get("Z",j))+
          fabs(Z+N-tab_rt00.get("A",j))<1.0e-3) {
        return tab_rt00.get("spin",j)*2.0+1.0;
      }
    }
  }

  if (hfb.is_included(Z,N)) {
    if (hfb.get_ZN(Z,N).Jexp<99) {
      return 2.0*hfb.get_ZN(Z,N).Jexp+1.0;
    } else {
      return 2.0*hfb.get_ZN(Z,N).Jth+1.0;
    }
  } else {
    if (Z%2==0 && N%2==0) {
      return 1.0;
    } else {
      return 2.0;
    }
  }
  
  O2SCL_ERR("Not found in get_spin_deg().",o2scl::exc_esanity);
  
  return 0.0;
}

int part_funcs::few78(int Z, int N, double T_K, double &pf,
                      double &TdpfdT) {
  return shen10(Z,N,T_K,pf,TdpfdT,0);
}

int part_funcs::rtk97(int Z, int N, double T_K, double &pf,
                      double &TdpfdT) {
  return shen10(Z,N,T_K,pf,TdpfdT,1);
}
  
int part_funcs::shen10(int Z, int N, double T_K, double &pf,
                       double &TdpfdT,
                       int a_delta) {

  /// Temperature in MeV
  double T_MeV=cu.convert("K","MeV",T_K);
    
  /// Nuclear level spacing in 1/MeV
  double a=0.0;
    
  /// Backshift parameter in MeV
  double delta=0.0;
    
  /// Coefficient for large \f$ \delta \f$ in 1/MeV
  double C;
    
  /// Critical temperature connecting low and high energies in MeV
  double Tc;

  double g=get_spin_deg(Z,N);

  if (Z==1 && N==1) {
    pf=3.0;
    TdpfdT=0.0;
    return 0;
  }

  // Get the neutron and proton separation energy 
  
  o2scl::nucleus nuc1, nuc2;
  // These are stored in MeV
  double Sneut=0.0, Sprot=0.0;

  if (ame.is_included(Z,N) && ame.is_included(Z-1,N) &&
      ame.is_included(Z,N-1)) {
    
    ame.get_nucleus(Z,N,nuc1);
    ame.get_nucleus(Z,N-1,nuc2);
    Sneut=-(nuc1.be-nuc2.be)*hc_mev_fm;
    ame.get_nucleus(Z-1,N,nuc2);
    Sprot=-(nuc1.be-nuc2.be)*hc_mev_fm;
    
  } else if (mnmsk.is_included(Z,N) && mnmsk.is_included(Z-1,N) &&
	     mnmsk.is_included(Z,N-1)) {
    
    mnmsk.get_nucleus(Z,N,nuc1);
    mnmsk.get_nucleus(Z,N-1,nuc2);
    Sneut=-(nuc1.be-nuc2.be)*hc_mev_fm;
    mnmsk.get_nucleus(Z-1,N,nuc2);
    Sprot=-(nuc1.be-nuc2.be)*hc_mev_fm;
    
  } else {
    Sneut=-1.0;
    Sprot=-1.0;

  }
          
  if (Sneut<0.0 || Sprot<0.0) {
    pf=g;
    TdpfdT=0.0;
    return 0;
  }

  // The nuclear mass in 1/fm
  double m=nuc1.m;//mnmsk.total_mass(Z,N);
  double n0=0.16;
  
  double res, err, ret, ret_prime, res_prime, err_prime;

  // zEd is in MeV
  double zEd=min(Sneut,Sprot)/2.0;
  // zR is in 1/fm  
  double zR=1.25*cbrt(Z+N-1.0);
  // zER is in MeV since the the nuclear mass is in 1/fm
  double zER=0.5/m/zR/zR*hc_mev_fm;
  // zEc is in MeV
  double zEc=(Z-1.0)*fine_structure_f<double>()/zR*hc_mev_fm;
  // zEt is in MeV
  double zEt=min(Sneut+zER,Sprot+zER+zEc/2.0);
  // zrA is in fm
  double zrA=cbrt(3.0*(N+Z)/4.0/pi/n0);

  // delta_p is in MeV
  double delta_p=11.0/sqrt(Z+N)*
    (1.0+pow(-1.0,Z)/2.0+pow(-1.0,N)/2.0);

  if (a_delta==0) {
    if (Z<=30) {
      // a is in 1/MeV
      a=0.052*pow(N+Z,1.2);
      // delta is in MeV
      delta=delta_p-80.0/(Z+N);
    } else {
      // a is in 1/MeV
      a=0.125*(N+Z);
      // delta is in MeV
      delta=delta_p-80.0/(Z+N)-0.5;
    }
  }
  
  if (!std::isfinite(a)) {
    cout << N << " " << Z << " " << delta_p << " " << a << endl;
    O2SCL_ERR2("Variable a not finite in ",
               " eos_nuclei::eos_fixed_dist()",
               o2scl::exc_esanity);
  }
  
  if (!std::isfinite(delta)) {
    cout << N << " " << Z << " " << delta << endl;
    cout << delta_p << endl;
    O2SCL_ERR2("Variable delta not finite",
               " eos_nuclei::eos_fixed_dist()",
               o2scl::exc_esanity);
  }
  
  if (delta>zEd) {

    // MeV
    delta=zEd;
    // MeV
    Tc=1.0/(-1.25/delta+sqrt(a)/sqrt(delta));
    // 1/MeV
    C=sqrt(pi)/12.0*pow(a,-0.25)*pow(delta,-1.25)*
      exp(1.25+sqrt(a*delta));
    
    // Set up integrands
    funct f1=std::bind(std::mem_fn<double(double,double,double,double)>
                       (&part_funcs::delta_small_iand),this,
                       std::placeholders::_1,T_MeV,delta,a);
    funct f1_prime=std::bind(std::mem_fn<double(double,double,double,double)>
                             (&part_funcs::delta_small_iand_prime),
                             this,std::placeholders::_1,T_MeV,delta,a);
    funct f2=std::bind(std::mem_fn<double(double,double,double,double,
                                          double)>
                       (&part_funcs::delta_large_iand),this,
                       std::placeholders::_1,T_MeV,delta,Tc,C);
    funct f2_prime=std::bind(std::mem_fn<double(double,double,double,
                                                double,double)>
                             (&part_funcs::delta_large_iand_prime),
                             this,std::placeholders::_1,
                             T_MeV,delta,Tc,C);
    
    if (2.0*zEd<zEt) {
      if (!std::isfinite(zEd)) {
        O2SCL_ERR2("Variable zEd not finite",
                   "eos_nuclei::eos_fixed_dist()",
                   o2scl::exc_esanity);
      }
      if (!std::isfinite(zEt)) {
        O2SCL_ERR2("Variable zEt not finite",
                   "eos_nuclei::eos_fixed_dist()",
                   o2scl::exc_esanity);
      }
      
      ret=iqg.integ_err(f1,2.0*zEd,zEt,res,err)+
        iqg.integ_err(f2,zEd,2.0*zEd,res,err);
      ret_prime=iqg.integ_err(f1_prime,2.0*zEd,zEt,res_prime,err_prime)+
        iqg.integ_err(f2_prime,zEd,2.0*zEd,res_prime,err_prime);
      
    } else {
      
      ret=iqg.integ_err(f2,zEd,zEt,res,err);
      ret_prime=iqg.integ_err(f2_prime,zEd,zEt,res_prime,err_prime);
      
    }
    
  } else {
    
    if (!std::isfinite(zEd)) {
      O2SCL_ERR2("Variable zEd not finite",
                 "eos_nuclei::eos_fixed_dist()",
                 o2scl::exc_esanity);
    }
    
    if (!std::isfinite(zEt)) {
      O2SCL_ERR2("Variable zEt not finite",
                 "eos_nuclei::eos_fixed_dist()",
                 o2scl::exc_esanity);
    }
    
    // Set up integrands
    funct f1=std::bind(std::mem_fn<double(double,double,double,double)>
                       (&part_funcs::delta_small_iand),this,
                       std::placeholders::_1,T_MeV,delta,a);
    funct f1_prime=std::bind(std::mem_fn<double(double,double,double,double)>
                             (&part_funcs::delta_small_iand_prime),
                             this,std::placeholders::_1,T_MeV,delta,a);
    
    ret=iqg.integ_err(f1,zEd,zEt,res,err);
    ret_prime=iqg.integ_err(f1_prime,zEd,zEt,res_prime,err_prime);
  }

  pf=g+res;
  // res_prime is, T \Omega^{\prime}, i.e. Eq. 26 of Shen et al. (2010),
  // which is unitless
  TdpfdT=res_prime;
  
  return 0;
}
