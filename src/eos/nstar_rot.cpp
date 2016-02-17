/*
  -------------------------------------------------------------------
  
  This file is part of O2scl. It has been adapted from RNS v1.1d
  written by N. Stergioulas and S. Morsink. The modifications made in
  this version from the original are copyright (C) 2015-2016, Andrew
  W. Steiner.
  
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
/*
  -------------------------------------------------------------------
  Relativistic models of rapidly rotating compact stars,
  using tabulated or polytropic equations of state.
  
  Author:  Nikolaos Stergioulas
  
  Version: 1.1
  
  Date:    June, 1995
  
  Changes made to code by Sharon Morsink
   
  -------------------------------------------------------------------
*/

// AWS, 11/11/15: cstdlib appears to be
// required for size_t in gsl_sf_legendre
#include <cstdlib>

#include <gsl/gsl_sf_legendre.h>

#include <o2scl/constants.h>
#include <o2scl/nstar_rot.h>

using namespace std;
using namespace o2scl;

eos_nstar_rot_interp::eos_nstar_rot_interp() {
  n_nearest=1;

  C=o2scl_cgs::speed_of_light;
  G=o2scl_cgs::gravitational_constant;
  KAPPA=1.0e-15*C*C/G;
  KSCALE=KAPPA*G/(C*C*C*C);
}

int eos_nstar_rot_interp::new_search(int n, double *x, double val) {
  int ret;
  double *xnew=x+1;
  bool inc=false;
  if (xnew[0]<xnew[n-1]) inc=true;
  if (inc) {
    if (val<xnew[0]) {
      return 0;
    } else if (val>xnew[n-1]) {
      return n;
    }
  } else {
    if (val>xnew[0]) {
      return 0;
    } else if (val<xnew[n-1]) {
      return n;
    }
  }
  sv.set_vec(n,xnew);
  return sv.find(val)+1;
}

double eos_nstar_rot_interp::ed_from_pr(double pr) {
  return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pr)));
}

double eos_nstar_rot_interp::pr_from_ed(double ed) {
  return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ed)));
}

double eos_nstar_rot_interp::nb_from_pr(double pr) {
  return pow(10.0,interp(log_p_tab,log_n0_tab,n_tab,log10(pr)));
}

double eos_nstar_rot_interp::pr_from_nb(double nb) {
  return pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(nb)));
}

double eos_nstar_rot_interp::ed_from_nb(double nb) {
  return pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(nb)));
}

double eos_nstar_rot_interp::nb_from_ed(double ed) {
  return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ed)));
}

double eos_nstar_rot_interp::enth_from_pr(double pr) {
  return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pr)));
}

double eos_nstar_rot_interp::enth_from_nb(double nb) {
  return pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(nb)));
}

double eos_nstar_rot_interp::pr_from_enth(double enth) {
  return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(enth)));
}

void eos_nstar_rot_interp::ed_nb_from_pr(double pr, double &ed, double &nb) {
  ed_from_pr(pr);
  nb_from_pr(pr);
  return;
}

double eos_nstar_rot_interp::interp(double xp[], double yp[], int np,
				    double xb) {
  // index of 1st point
  int k;        
  // degree of interpolation
  int m=4;      
 
  // intermediate value
  double y;     

  n_nearest=new_search(np,xp,xb);

  int max=n_nearest-(m-1)/2;
  if (max<1) max=1;
  k=np+1-m;
  if (max<k) k=max;

  if (xb==xp[k] || xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) {
    xb+=1.0e-12;
  }
  
  y=(xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
    ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
    ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
    ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
    ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));
  
  return (y);
}

eos_nstar_rot_C::eos_nstar_rot_C(bool rns_constants) {
  if (rns_constants) {
    C=2.9979e10;                  
    G=6.6732e-8;                  
    KAPPA=1.346790806509621e+13;  
    KSCALE=1.112668301525780e-36; 
  }
  double eosC_arr[96][4]={
    {7.800e+00,1.010e+08,1.000000000000000e+00,4.698795180722962e+24},
    {7.860e+00,1.010e+09,1.157946629165269e+08,4.734939759036205e+24},
    {7.900e+00,1.010e+10,1.269452049889617e+09,4.759036144578364e+24},
    {8.150e+00,1.010e+11,1.267708525005907e+10,4.909638554215315e+24},
    {1.160e+01,1.210e+12,1.211370595001572e+11,6.987951807098076e+24},
    {1.640e+01,1.400e+13,1.017364459510011e+12,9.879518070489597e+24},
    {4.510e+01,1.700e+14,6.076705858546280e+12,2.716867462904601e+25},
    {2.120e+02,5.820e+15,4.872391666226939e+13,1.277108403508764e+26},
    {1.150e+03,1.900e+17,3.206724388438867e+14,6.927709645088004e+26},
    {1.044e+04,9.744e+18,2.085685492452927e+15,6.289148562640985e+27},
    {2.622e+04,4.968e+19,4.300724422116231e+15,1.579513843816999e+28},
    {6.587e+04,2.431e+20,8.585327648535554e+15,3.968050678245718e+28},
    {1.654e+05,1.151e+21,1.661134940613050e+16,9.963748410271617e+28},
    {4.156e+05,5.266e+21,3.113184639693159e+16,2.503563031417219e+29},
    {1.044e+06,2.318e+22,5.637078789809274e+16,6.288917532113082e+29},
    {2.622e+06,9.755e+22,9.823793802270347e+16,1.579410809416864e+30},
    {6.588e+06,3.911e+23,1.651370193851722e+17,3.968207649843547e+30},
    {8.293e+06,5.259e+23,1.833070570850680e+17,4.995116726219748e+30},
    {1.655e+07,1.435e+24,2.575234489468157e+17,9.967984755458204e+30},
    {3.302e+07,3.833e+24,3.565410566998838e+17,1.988624478073943e+31},
    {6.589e+07,1.006e+25,4.855034973143420e+17,3.967807406359445e+31},
    {1.315e+08,2.604e+25,6.514242926503165e+17,7.917691186982454e+31},
    {2.624e+08,6.676e+25,8.653913867049318e+17,1.579648605894070e+32},
    {3.304e+08,8.738e+25,9.351655321505760e+17,1.988876577393412e+32},
    {5.237e+08,1.629e+26,1.113042991360343e+18,3.152005155076383e+32},
    {8.301e+08,3.029e+26,1.322173059425539e+18,4.995278531652059e+32},
    {1.045e+09,4.129e+26,1.440858462676231e+18,6.287859551784352e+32},
    {1.316e+09,5.036e+26,1.518045189928309e+18,7.917701445937253e+32},
    {1.657e+09,6.860e+26,1.639959584391741e+18,9.968319738044036e+32},
    {2.626e+09,1.272e+27,1.916631713149610e+18,1.579408507997411e+33},
    {4.164e+09,2.356e+27,2.239467717942762e+18,2.503766293549853e+33},
    {6.601e+09,4.362e+27,2.618527558269814e+18,3.967852390467774e+33},
    {8.312e+09,5.662e+27,2.793223911233673e+18,4.995474308724729e+33},
    {1.046e+10,7.702e+27,3.010604603009213e+18,6.285277578607203e+33},
    {1.318e+10,1.048e+28,3.246131568141483e+18,7.918132634568090e+33},
    {1.659e+10,1.425e+28,3.499915730795552e+18,9.964646988214994e+33},
    {2.090e+10,1.938e+28,3.774741006863546e+18,1.255052800774333e+34},
    {2.631e+10,2.503e+28,4.014780206768711e+18,1.579545673652798e+34},
    {3.313e+10,3.404e+28,4.317831726448436e+18,1.988488463504033e+34},
    {4.172e+10,4.628e+28,4.646291608758636e+18,2.503379640977065e+34},
    {5.254e+10,5.949e+28,4.927359777046148e+18,3.151720931652274e+34},
    {6.617e+10,8.089e+28,5.287586080904092e+18,3.968151735612910e+34},
    {8.332e+10,1.100e+29,5.677637396835295e+18,4.994995310195290e+34},
    {1.049e+11,1.495e+29,6.097968154425379e+18,6.286498800006776e+34},
    {1.322e+11,2.033e+29,6.553704076675019e+18,7.919521253825185e+34},
    {1.664e+11,2.597e+29,6.932922860174182e+18,9.964341016667146e+34},
    {1.844e+11,2.892e+29,7.100996231341637e+18,1.104024323001462e+35},
    {2.096e+11,3.290e+29,7.302912999460339e+18,1.254619611126682e+35},
    {2.640e+11,4.473e+29,7.801189133603082e+18,1.579588892045295e+35},
    {3.325e+11,5.816e+29,8.252940103235718e+18,1.988565738933728e+35},
    {4.188e+11,7.538e+29,8.710551350551819e+18,2.503561780689725e+35},
    {4.299e+11,7.805e+29,8.774012262748510e+18,2.569780082714395e+35},
    {4.460e+11,7.890e+29,8.793402718214299e+18,2.665824694449485e+35},
    {5.228e+11,8.352e+29,8.888584165828376e+18,3.123946525953616e+35},
    {6.610e+11,9.098e+29,9.015182344330039e+18,3.948222384313103e+35},
    {7.964e+11,9.831e+29,9.115886202428306e+18,4.755697604312120e+35},
    {9.728e+11,1.083e+30,9.228938242554155e+18,5.807556544067428e+35},
    {1.196e+12,1.218e+30,9.353548588340060e+18,7.138304213736713e+35},
    {1.471e+12,1.399e+30,9.489304401520411e+18,8.777653631971616e+35},
    {1.805e+12,1.683e+30,9.662916598353355e+18,1.076837272716171e+36},
    {2.202e+12,1.950e+30,9.796648174499881e+18,1.313417953138369e+36},
    {2.930e+12,2.592e+30,1.004639229994465e+19,1.747157788902558e+36},
    {3.833e+12,3.506e+30,1.031719959699455e+19,2.285004034820638e+36},
    {4.933e+12,4.771e+30,1.060626231961342e+19,2.939983642627298e+36},
    {6.248e+12,6.481e+30,1.091248816114249e+19,3.722722765704268e+36},
    {7.801e+12,8.748e+30,1.123546353605510e+19,4.646805278760175e+36},
    {9.611e+12,1.170e+31,1.157469224223553e+19,5.723413975645761e+36},
    {1.246e+13,1.695e+31,1.205104455978235e+19,7.417258934884369e+36},
    {1.496e+13,2.209e+31,1.242585565332612e+19,8.902909532230595e+36},
    {1.778e+13,2.848e+31,1.281598407175551e+19,1.057801059193907e+37},
    {2.210e+13,3.931e+31,1.335910916365704e+19,1.314278492046241e+37},
    {2.988e+13,6.178e+31,1.422481793925897e+19,1.775810743961577e+37},
    {3.767e+13,8.774e+31,1.499308970128912e+19,2.237518046976615e+37},
    {5.081e+13,1.386e+32,1.614317463895106e+19,3.015480061626022e+37},
    {6.193e+13,1.882e+32,1.702142123464848e+19,3.673108933334910e+37},
    {7.732e+13,2.662e+32,1.813939645708965e+19,4.582250451016437e+37},
    {9.826e+13,3.897e+32,1.954253029978894e+19,5.817514573447143e+37},
    {1.262e+14,5.861e+32,2.128347366138737e+19,7.462854442694524e+37},
    {1.706e+14,1.756e+33,2.885278398767732e+19,1.006639916443579e+38},
    {2.567e+14,4.565e+33,4.178979788188595e+19,1.505335697605081e+38},
    {3.458e+14,9.397e+33,5.738041973520725e+19,2.013381591984608e+38},
    {4.350e+14,1.657e+34,7.507944724919534e+19,2.512629566428361e+38},
    {5.277e+14,2.625e+34,9.423281380878824e+19,3.020920486283259e+38},
    {7.166e+14,5.550e+34,1.379858119588385e+20,4.021520391505177e+38},
    {9.163e+14,1.000e+35,1.872894652613005e+20,5.025727131203635e+38},
    {1.128e+15,1.630e+35,2.412896024049271e+20,6.030527187437616e+38},
    {1.353e+15,2.418e+35,2.951587122936742e+20,7.036003566506166e+38},
    {1.596e+15,3.385e+35,3.490539607233137e+20,8.058746778042895e+38},
    {1.847e+15,4.518e+35,4.015539689708425e+20,9.054485628019740e+38},
    {2.121e+15,5.898e+35,4.554634605761517e+20,1.007895513739328e+39},
    {3.726e+15,1.614e+36,7.092250567926366e+20,1.510986122033221e+39},
    {5.812e+15,3.289e+36,9.371241098390350e+20,2.011274032860530e+39},
    {8.468e+15,5.718e+36,1.139974260034128e+21,2.512899672966290e+39},
    {1.175e+16,8.982e+36,1.320440905946775e+21,3.014116952600337e+39},
    {2.032e+16,1.825e+37,1.626616375316661e+21,4.005670613879821e+39},
    {3.227e+16,3.204e+37,1.886033418779976e+21,5.017750644193835e+39}};

  n_tab=96;

  if (false) {
    table<> t;
    convert_units &cu=o2scl_settings.get_convert_units();
    t.line_of_names("ed pr h nb");
    for(size_t i=0;i<n_tab;i++) {
      double line[4]={cu.convert("g/cm^3","MeV/fm^3",eosC_arr[i][0]),
		      cu.convert("dyne/cm^2","MeV/fm^3",eosC_arr[i][1]),
		      eosC_arr[i][2],
		      cu.convert("1/cm^3","1/fm^3",eosC_arr[i][3])};
      t.line_of_data(4,line);
    }
    t.new_column("iand");
    t.add_constant("c",o2scl_cgs::speed_of_light);
    t.set_interp_type(itp_linear);
    t.function_column("c^2/(ed+pr)","iand");
    t.integ("pr","iand","h2");
    t.function_column("c^2*(log((ed+pr)/nb))","h3");
    double h30=t.get("h3",0);
    for(size_t i=0;i<n_tab;i++) {
      t.set("h3",i,t.get("h3",i)-h30);
    }
    for(size_t i=0;i<n_tab;i++) {
      cout << t.get("pr",i) << " "
	   << t.get("iand",i) << " "
	   << t.get("h",i) << " "
	   << t.get("h2",i) << " "
	   << t.get("h3",i) << endl;
    }
    exit(-1);
  }
  
  for(int i=1;i<=n_tab;i++) {  
    
    double rho=eosC_arr[i-1][0];
    double p=eosC_arr[i-1][1];
    double h=eosC_arr[i-1][2];
    double n0=eosC_arr[i-1][3];
    
    //cout << (rho+p/pow(o2scl_cgs::speed_of_light,2.0))/n0 << " "
    //<< 1.66e-24*exp(h/pow(o2scl_cgs::speed_of_light,2.0)) << endl;

    log_e_tab[i]=log10(rho*C*C*KSCALE);
    log_p_tab[i]=log10(p*KSCALE);
    log_h_tab[i]=log10(h/(C*C));
    log_n0_tab[i]=log10(n0);
  }
  //exit(-1);
}

eos_nstar_rot_L::eos_nstar_rot_L(bool rns_constants) {
  if (rns_constants) {
    C=2.9979e10;                  
    G=6.6732e-8;                  
    KAPPA=1.346790806509621e+13;  
    KSCALE=1.112668301525780e-36; 
  }
  double eosL_arr[64][4]={
    {7.800e+00,1.010e+08,1.000000000000000e+00,4.698795180722962e+24},
    {7.860e+00,1.010e+09,1.157946629165269e+08,4.734939759036205e+24},
    {7.900e+00,1.010e+10,1.269452049889617e+09,4.759036144578364e+24},
    {8.150e+00,1.010e+11,1.267708525005907e+10,4.909638554215315e+24},
    {1.160e+01,1.210e+12,1.211370595001572e+11,6.987951807098076e+24},
    {1.640e+01,1.400e+13,1.017364459510011e+12,9.879518070489597e+24},
    {4.510e+01,1.700e+14,6.076705858546280e+12,2.716867462904601e+25},
    {2.120e+02,5.820e+15,4.872391666226939e+13,1.277108403508764e+26},
    {1.150e+03,1.900e+17,3.206724388438867e+14,6.927709645088004e+26},
    {1.044e+04,9.744e+18,2.085685492452927e+15,6.289148562640985e+27},
    {2.622e+04,4.968e+19,4.300724422116231e+15,1.579513843816999e+28},
    {6.587e+04,2.431e+20,8.585327648535554e+15,3.968050678245718e+28},
    {1.654e+05,1.151e+21,1.661134940613050e+16,9.963748410271617e+28},
    {4.156e+05,5.266e+21,3.113184639693159e+16,2.503563031417219e+29},
    {1.044e+06,2.318e+22,5.637078789809274e+16,6.288917532113082e+29},
    {2.622e+06,9.755e+22,9.823793802270347e+16,1.579410809416864e+30},
    {6.588e+06,3.911e+23,1.651370193851722e+17,3.968207649843547e+30},
    {8.293e+06,5.259e+23,1.833070570850680e+17,4.995116726219748e+30},
    {1.655e+07,1.435e+24,2.575234489468157e+17,9.967984755458204e+30},
    {3.302e+07,3.833e+24,3.565410566998838e+17,1.988624478073943e+31},
    {6.589e+07,1.006e+25,4.855034973143420e+17,3.967807406359445e+31},
    {1.315e+08,2.604e+25,6.514242926503165e+17,7.917691186982454e+31},
    {2.624e+08,6.676e+25,8.653913867049318e+17,1.579648605894070e+32},
    {3.304e+08,8.738e+25,9.351655321505760e+17,1.988876577393412e+32},
    {5.237e+08,1.629e+26,1.113042991360343e+18,3.152005155076383e+32},
    {8.301e+08,3.029e+26,1.322173059425539e+18,4.995278531652059e+32},
    {1.045e+09,4.129e+26,1.440858462676231e+18,6.287859551784352e+32},
    {1.316e+09,5.036e+26,1.518045189928309e+18,7.917701445937253e+32},
    {1.657e+09,6.860e+26,1.639959584391741e+18,9.968319738044036e+32},
    {2.626e+09,1.272e+27,1.916631713149610e+18,1.579408507997411e+33},
    {4.164e+09,2.356e+27,2.239467717942762e+18,2.503766293549853e+33},
    {6.601e+09,4.362e+27,2.618527558269814e+18,3.967852390467774e+33},
    {8.312e+09,5.662e+27,2.793223911233673e+18,4.995474308724729e+33},
    {1.046e+10,7.702e+27,3.010604603009213e+18,6.285277578607203e+33},
    {1.318e+10,1.048e+28,3.246131568141483e+18,7.918132634568090e+33},
    {1.659e+10,1.425e+28,3.499915730795552e+18,9.964646988214994e+33},
    {2.090e+10,1.938e+28,3.774741006863546e+18,1.255052800774333e+34},
    {2.631e+10,2.503e+28,4.014780206768711e+18,1.579545673652798e+34},
    {3.313e+10,3.404e+28,4.317831726448436e+18,1.988488463504033e+34},
    {4.172e+10,4.628e+28,4.646291608758636e+18,2.503379640977065e+34},
    {5.254e+10,5.949e+28,4.927359777046148e+18,3.151720931652274e+34},
    {6.617e+10,8.089e+28,5.287586080904092e+18,3.968151735612910e+34},
    {8.332e+10,1.100e+29,5.677647984386809e+18,4.994995251352083e+34},
    {1.000e+11,1.402e+29,6.007356633978408e+18,5.993299174014897e+34},
    {2.000e+11,3.134e+29,7.206489006907356e+18,1.197281032014887e+35},
    {4.000e+11,7.157e+29,8.672720920156850e+18,2.391248847376186e+35},
    {8.000e+11,1.036e+30,9.228502681906242e+18,4.776917888094152e+35},
    {1.000e+12,1.257e+30,9.473901401216950e+18,5.969265140210954e+35},
    {2.000e+12,2.122e+30,1.008175138682551e+19,1.192786036339965e+36},
    {4.000e+12,3.780e+30,1.065412286010822e+19,2.383745936109886e+36},
    {8.000e+12,8.527e+30,1.145481617345254e+19,4.763886330347948e+36},
    {1.000e+13,1.162e+31,1.179848405023799e+19,5.953217185533753e+36},
    {2.000e+13,3.262e+31,1.321734278540367e+19,1.189384774424109e+37},
    {4.000e+13,9.407e+31,1.529633911184461e+19,2.375173071697946e+37},
    {8.000e+13,2.746e+32,1.834132290018557e+19,4.739957123069321e+37},
    {1.000e+14,3.929e+32,1.965106871694968e+19,5.919574167379903e+37},
    {2.000e+14,3.380e+33,3.747246764168948e+19,1.177348214981387e+38},
    {4.000e+14,3.255e+34,1.262888655801689e+20,2.283330744375166e+38},
    {8.000e+14,1.927e+35,3.530690068286150e+20,4.125664617034758e+38},
    {1.000e+15,2.968e+35,4.423515236827836e+20,4.898535647530366e+38},
    {2.000e+15,8.092e+35,6.966531453934924e+20,8.048311708651028e+38},
    {4.000e+15,1.863e+36,9.409246465914199e+20,1.284130443857022e+39},
    {8.000e+15,4.859e+36,1.261068058921259e+21,1.985296494438113e+39},
    {1.000e+16,6.551e+36,1.371873265801520e+21,2.263286974848659e+39}};

  n_tab=64;

  for(int i=1;i<=n_tab;i++) {  

    double rho=eosL_arr[i-1][0];
    double p=eosL_arr[i-1][1];
    double h=eosL_arr[i-1][2];
    double n0=eosL_arr[i-1][3];

    log_e_tab[i]=log10(rho*C*C*KSCALE);
    log_p_tab[i]=log10(p*KSCALE);
    log_h_tab[i]=log10(h/(C*C));
    log_n0_tab[i]=log10(n0);

  }
}

nstar_rot::nstar_rot() {
  verbose=1;
  
  SMAX=0.9999;                
  DS=(SMAX/(SDIV-1.0));
  DM=(1.0/(MDIV-1.0));          
  RMIN=1.0e-15;
  n_nearest=1;                     
  a_check=0;                       
  s_e=0.5;

  C=o2scl_cgs::speed_of_light;
  G=o2scl_cgs::gravitational_constant;
  MSUN=o2scl_cgs::solar_mass;
  PI=o2scl_const::pi;
  MB=o2scl_cgs::mass_neutron;
  KAPPA=1.0e-15*C*C/G;
  KSCALE=KAPPA*G/(C*C*C*C);

  // set program defaults
  cf=1.0;
  eq_radius_tol_rel=1.0e-5;    
  tol_abs=1.0e-4;   
 
  CL_LOW=false;
  scaled_polytrope=false;

  // Default polytropic index
  n_P=1.0;

  eos_set=false;
  
  // Create the computational mesh for variables "s" and "mu=cos theta"
  make_grid();
  
  /* Create the 2-point functions and legendre polynomials needed
     to integrate the metric potentials rho, gamma and omega 
  */
  comp_f_P();
}

void nstar_rot::constants_rns() {
  C=2.9979e10;                  
  G=6.6732e-8;                  
  KAPPA=1.346790806509621e+13;  
  KSCALE=1.112668301525780e-36; 
  MSUN=1.987e33;                
  PI=3.1415926535;              
  MB=1.66e-24;
  return;
}

void nstar_rot::constants_o2scl() {
  C=o2scl_cgs::speed_of_light;
  G=o2scl_cgs::gravitational_constant;
  MSUN=o2scl_cgs::solar_mass;
  PI=o2scl_const::pi;
  MB=o2scl_cgs::mass_neutron;
  KAPPA=1.0e-15*C*C/G;
  KSCALE=KAPPA*G/(C*C*C*C);
  return;
}

void nstar_rot::make_grid() {                                    
    
  for(int s=1;s<=SDIV;s++) {
    s_gp[s]=SMAX*(1.0*s-1.0)/(SDIV-1);
    s_1_s[s]=s_gp[s]*(1.0-s_gp[s]);
    one_s[s]=1.0-s_gp[s];
  }
  
  for(int m=1;m<=MDIV;m++) {
    mu[m]=(1.0*m-1.0)/(MDIV-1);
    one_m2[m]=1.0-mu[m]*mu[m];
    theta[m]=acos(mu[m]);
    sin_theta[m]=sqrt(one_m2[m]);
  }

  return;
}

int nstar_rot::new_search(int n, double *x, double val) {
  int ret;
  double *xnew=x+1;
  bool inc=false;
  if (xnew[0]<xnew[n-1]) inc=true;
  if (inc) {
    if (val<xnew[0]) {
      return 0;
    } else if (val>xnew[n-1]) {
      return n;
    }
  } else {
    if (val>xnew[0]) {
      return 0;
    } else if (val<xnew[n-1]) {
      return n;
    }
  }
  sv.set_vec(n,xnew);
  return sv.find(val)+1;
}

double nstar_rot::interp(double xp[], double yp[], int np, double xb) {
  // index of 1st point
  int k;        
  // degree of interpolation
  int m=4;      
 
  // intermediate value
  double y;     

  n_nearest=new_search(np,xp,xb);

  int max=n_nearest-(m-1)/2;
  if (max<1) max=1;
  k=np+1-m;
  if (max<k) k=max;

  if (xb==xp[k] || xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) {
    xb+=1.0e-12;
  }

  y=(xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
    ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
    ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
    ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))+
    (xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
    ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));
  
  return (y);
}

double nstar_rot::interp_4_k(double xp[], double yp[], int np, double xb, 
			      int k) {

  // intermediate value
  double y;     

  if (xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) {
    xb+=1.0e-14;
  }
  
  y=(xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
    ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
    ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
    ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
    ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

  return(y);
}

double nstar_rot::int_z(double f[MDIV+1], int m) {
  double x[9];

  x[1]=f[m-1];  
  x[2]=interp(mu,f,MDIV,mu[m-1]+DM/7.0);
  x[3]=interp(mu,f,MDIV,mu[m-1]+2.0*DM/7.0);
  x[4]=interp(mu,f,MDIV,mu[m-1]+3.0*DM/7.0);
  x[5]=interp(mu,f,MDIV,mu[m-1]+4.0*DM/7.0);
  x[6]=interp(mu,f,MDIV,mu[m-1]+5.0*DM/7.0);
  x[7]=interp(mu,f,MDIV,mu[m-1]+6.0*DM/7.0);
  x[8]=f[m];
  
  return((DM/17280.0)*(751.0*x[1]+3577.0*x[2]+1323.0*x[3]+2989.0*x[4]+
		       2989.0*x[5]+1323.0*x[6]+3577.0*x[7]+751.0*x[8]));
}

double nstar_rot::e_at_p(double pp) {
  if (scaled_polytrope==false) {
    if (CL_LOW==true && pp > p_match) {
      return pp+e_match+de_pt-p_match;
    } else {
      return eosp->ed_from_pr(pp);
      //return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp)));
    }
  } else {
    return pp/(Gamma_P-1.0)+pow(pp,1.0/Gamma_P); 
  }
}

double nstar_rot::p_at_e(double ee) {
  if (CL_LOW==true && ee > e_match) {
    if (ee<=e_cl) {
      return p_match; 
    } else {
      return ee-e_match-de_pt+p_match;
    }
  } else {
    return eosp->pr_from_ed(ee);
    //return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee)));
  }
} 

double nstar_rot::p_at_h(double hh) {
  if (CL_LOW==true && hh > h_match) {
    return 0.5*( (e_match+de_pt+p_match)*exp(2.0*(hh-h_match)) 
		 +p_match-e_match-de_pt);
  } else {
    return eosp->pr_from_enth(hh);
    //return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh)));
  }
}

double nstar_rot::h_at_p(double pp) {
  if (CL_LOW==true && pp > p_match) {
    return h_match+0.5*log( (2.0*pp+e_match+de_pt-p_match)/
			    (e_match+de_pt+p_match));
  } else {
    return eosp->enth_from_pr(pp);
    //return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp)));
  }
}

double nstar_rot::n0_at_e(double ee) {
  if (CL_LOW==true && ee > e_match) {
    if (ee<=e_cl) {
      return ((ee+p_match)/(MB*C*C*KSCALE))*exp(-h_match);
    } else {
      return (n0_match+(de_pt/(MB*C*C*KSCALE))*exp(-h_match))  
	*sqrt(1.0+2.0*(ee-e_match-de_pt)/(e_match+de_pt+p_match));
    }
  } else {
    return eosp->nb_from_ed(ee);
    //return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee)));
  }
}

double nstar_rot::s_deriv(double f[SDIV+1], int s) {

  double d_temp;
  
  switch(s) { 
  case 1 : 
    d_temp=(f[s+1]-f[s])/DS;
    break;
  case SDIV : 
    d_temp=(f[s]-f[s-1])/DS;
    break;
  default : 
    d_temp=(f[s+1]-f[s-1])/(2.0*DS);
    break; 
  }
 
  return d_temp;
}

double nstar_rot::m_deriv(double f[MDIV+1], int m) {

  double d_temp;

  switch(m) {  
  case 1 : 
    d_temp=(f[m+1]-f[m])/DM;
    break;
  case MDIV : 
    d_temp=(f[m]-f[m-1])/DM;
    break;
  default : 
    d_temp=(f[m+1]-f[m-1])/(2.0*DM);
    break; 
  }
 
  return d_temp;
}

double nstar_rot::deriv_s(double f[SDIV+1][MDIV+1], int s, int m) {

  double d_temp;

  switch(s) { 
  case 1 : 
    d_temp=(f[s+1][m]-f[s][m])/DS;
    break;
  case SDIV : 
    d_temp=(f[s][m]-f[s-1][m])/DS;
    break;
  default : 
    d_temp=(f[s+1][m]-f[s-1][m])/(2.0*DS);
    break; 
  } 

  return d_temp;
}

double nstar_rot::deriv_m(double f[SDIV+1][MDIV+1], int s, int m) {

  double d_temp;

  switch(m) { 
  case 1 : 
    d_temp=(f[s][m+1]-f[s][m])/DM;
    break;
  case MDIV : 
    d_temp=(f[s][m]-f[s][m-1])/DM;
    break;
  default : 
    d_temp=(f[s][m+1]-f[s][m-1])/(2.0*DM);
    break; 
  } 

  return d_temp;
}

double nstar_rot::deriv_sm(double f[SDIV+1][MDIV+1], int s, int m) {

  double d_temp;

  switch(s) {

  case 1 : 
    if (m==1) {   
      d_temp=(f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
    } else {
      if (m==MDIV) {
	d_temp=(f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
      } else {         
	d_temp=(f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/
	  (2.0*DM*DS);
      }
    }
    break;

  case SDIV : 
    if (m==1) {   
      d_temp=(f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
    } else {
      if (m==MDIV) {
	d_temp=(f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
      } else {         
	d_temp=(f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/
	  (2.0*DM*DS);
      }
    }
    break;
  
  default : 

    if (m==1) {   
      d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
    } else {
      if (m==MDIV) {
	d_temp=(f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/
	  (2.0*DM*DS);
      } else {         
	d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/
	  (4.0*DM*DS);
      }
    }
    break;

  }

  return d_temp;
}

double nstar_rot::legendre(int n, double x) {
  // counter
  int i;           

  // Legendre polynomials of order n, n-1, and n-2 
  double p;   
  double p_1=x;
  double p_2=1.0;      

  if (n>=2) { 
    for(i=2;i<=n;i++) {
      p=(x*(2.0*i-1.0)*p_1-(i-1.0)*p_2)/i;
      p_2=p_1;
      p_1=p;
    }
    return p;
  } else { 
    if (n==1) return p_1;
    else return p_2;
  }
}

void nstar_rot::comp_f_P() {
  // counter
  int n;                 
  // counter for s'
  int k;                 
  // counter for s
  int j;                 
  int s_temp;

  double sj, sk, sj1, sk1;

  if (SMAX==1.0) s_temp=SDIV-1;
  else s_temp=SDIV;

  // n=0, k>j case

  j=1;
 
  n=0; 
  for(k=2;k<=SDIV;k++) {
    sk=s_gp[k];
    f_rho[j][n][k]=1.0/(sk*(1.0-sk));
  }

  // n=1, k>j case

  n=1;
  for(k=2;k<=SDIV;k++) {
    sk=s_gp[k];
    sk1=1.0-sk;
         
    f_rho[j][n][k]=0.0;
    f_gamma[j][n][k]=1.0/(sk*sk1);
    f_omega[j][n][k]=1.0/(sk*sk1);
  }

  // n>=2, k>=j case

  for(n=2;n<=LMAX;n++) {
    for(k=1;k<=SDIV;k++) {
      f_rho[j][n][k]=0.0;
      f_gamma[j][n][k]=0.0;
      f_omega[j][n][k]=0.0;
    }
  }

  // n=0, j>=k case

  k=1;

  n=0;
  for(j=1;j<=SDIV;j++) {
    f_rho[j][n][k]=0.0;
  }

  // n>=1, j>=k case

  for(j=1;j<=SDIV;j++) {
    for(n=1;n<=LMAX;n++) {
      f_rho[j][n][k]=0.0;
      f_gamma[j][n][k]=0.0;
      f_omega[j][n][k]=0.0;
    }
  }

  // n=0, ? case

  n=0;
  for(j=2;j<=SDIV;j++) {
    for(k=2;k<=SDIV;k++) {

      if (SMAX==1.0 && (k==SDIV || j==SDIV)) {
	f_rho[j][n][k]=0.0;
      } else {      
	sk=s_gp[k];
	sj=s_gp[j];
	sk1=1.0-sk;
	sj1=1.0-sj;

	if (k<j) {
	  f_rho[j][n][k]=pow(sj1/sj,2.0*n+1.0)*
	    pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
	} else {
	  f_rho[j][n][k]=pow(sj/sj1,2.0*n)*
	    pow(sk1,2.0*n-1.0)/pow(sk,2.0*n+1.0);
	}
      }
    }
  }

  // General case

  for(j=2;j<=SDIV;j++) {
    for(n=1;n<=LMAX;n++) {
      for(k=2;k<=SDIV;k++) {

	if (SMAX==1.0 && (k==SDIV || j==SDIV)) {
	  f_rho[j][n][k]=0.0;
	  f_gamma[j][n][k]=0.0;
	  f_omega[j][n][k]=0.0;
	} else {      
	  sk=s_gp[k];
	  sj=s_gp[j];
	  sk1=1.0-sk;
	  sj1=1.0-sj;

	  if (k<j) {   
	    f_rho[j][n][k]=pow(sj1/sj,2.0*n+1.0)*
	      pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
	    f_gamma[j][n][k]=pow(sj1/sj,2.0*n)*
	      pow(sk,2.0*n-1.0)/pow(sk1,2.0*n+1.0);
	    f_omega[j][n][k]=pow(sj1/sj,2.0*n+1.0)*
	      pow(sk,2.0*n)/pow(sk1,2.0*n+2.0);
	  } else {     
	    f_rho[j][n][k]=pow(sj/sj1,2.0*n)*
	      pow(sk1,2.0*n-1.0)/pow(sk,2.0*n+1.0);
	    f_gamma[j][n][k]=pow(sj/sj1,2.0*n-2.0)*
	      pow(sk1,2.0*n-3.0)/pow(sk,2.0*n-1.0);
	    f_omega[j][n][k]=pow(sj/sj1,2.0*n-2.0)*
	      pow(sk1,2.0*n-3.0)/pow(sk,2.0*n-1.0);
	  }
	}
      }

    }
  }

  // counter for mu grid
  int i;

  // n=0 case

  n=0;
  for(i=1;i<=MDIV;i++) {
    P_2n[i][n]=legendre(2*n,mu[i]);
  }

  // n>=1 case

  for(i=1;i<=MDIV;i++) {
    for(n=1;n<=LMAX;n++) {
      P_2n[i][n]=legendre(2*n,mu[i]);
      P1_2n_1[i][n]=gsl_sf_legendre_Plm(2*n-1,1,mu[i]);
    }
  }

}

void nstar_rot::make_center(double e_center_loc) {
  double rho0_center;

  if (scaled_polytrope==false) {

    p_center=p_at_e(e_center_loc);
    h_center=h_at_p(p_center);

  } else {

    rbc.tol_abs=1.0e-16;
    polytrope_solve ps(Gamma_P,e_center_loc);
    rho0_center=0.0;
    rbc.solve_bkt(rho0_center,e_center_loc,ps);

    p_center=pow(rho0_center,Gamma_P);
    h_center=log((e_center_loc+p_center)/rho0_center);
  }

  return;
}

void nstar_rot::comp_omega() {
  int s;
  int m;

  double d_o_e[SDIV+1];
  double d_g_e[SDIV+1];
  double d_r_e[SDIV+1];
  double doe;
  double dge; 
  double dre;
  double vek;
  // rho at equator
  double rho_equator;               
  // omega at equator
  double omega_equator;             

  if (scaled_polytrope==false) {
    Omega=Omega_h*C/(r_e*sqrt(KAPPA));
  } else {
    // Omega=Omega_h*C/r_e;
    Omega=Omega_h/r_e;
  }

  // Kepler angular velocity

  for(s=1;s<=SDIV;s++) {
    rho_mu_0[s]=rho[s][1];                     
    omega_mu_0[s]=omega[s][1];                 
  }

  n_nearest=SDIV/2;
  rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

  if (r_ratio==1.0) {
    omega_equator=0.0;
  } else {
    omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);
  }

  for(s=1;s<=SDIV;s++) { 
    d_o_e[s]=deriv_s(omega,s,1);
    d_g_e[s]=deriv_s(gamma,s,1);
    d_r_e[s]=deriv_s(rho,s,1);
  }

  doe=interp(s_gp,d_o_e,SDIV,s_e);
  dge=interp(s_gp,d_g_e,SDIV,s_e);
  dre=interp(s_gp,d_r_e,SDIV,s_e);

  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator)+
    sqrt(((dge+dre)/(8.0+dge-dre))+pow((doe/(8.0+dge-dre))*
				       r_e*exp(-rho_equator),2.0));
  
  if (scaled_polytrope==false) {
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  } else {
    Omega_K=omega_equator+vek*exp(rho_equator)/r_e;
  }

  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      gamma_guess[s][m]=gamma[s][m];
      rho_guess[s][m]=rho[s][m];
      alpha_guess[s][m]=alpha[s][m];
      omega_guess[s][m]=omega[s][m];
    }
  }

  r_e_guess=r_e;
}
  
void nstar_rot::comp_M_J() {
  int s;
  int m;
  int s_temp;

  // int. quantity for M_0
  double D_m_0[SDIV+1];             
  // int. quantity for M
  double D_m[SDIV+1];               
  // int. quantity for J
  double D_J[SDIV+1];               
  // rest mass density
  double rho_0[SDIV+1][MDIV+1];     
  double rho_mu_0[SDIV+1];
  double omega_mu_0[SDIV+1];
  double rho_equator;
  double omega_equator;                      
  double d_o_e[SDIV+1];
  double d_g_e[SDIV+1];
  double d_r_e[SDIV+1];
  double doe;
  double dge;
  double dre;
  double vek;

  // Angular velocity

  if (r_ratio==1.0) {
    Omega=0.0;
  } else {
    if (scaled_polytrope==false) {
      Omega=Omega_h*C/(r_e*sqrt(KAPPA));
    } else {
      // Omega=Omega_h*C/r_e;
      Omega=Omega_h/r_e;
    }
  }   

  // Kepler angular velocity

  for(s=1;s<=SDIV;s++) {
    rho_mu_0[s]=rho[s][1];                     
    omega_mu_0[s]=omega[s][1];                 
  }

  n_nearest=SDIV/2;
  rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

  if (r_ratio==1.0) {
    omega_equator=0.0; 
  } else {
    omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);
  }
 
  for(s=1;s<=SDIV;s++) { 
    d_o_e[s]=deriv_s(omega,s,1);
    d_g_e[s]=deriv_s(gamma,s,1);
    d_r_e[s]=deriv_s(rho,s,1);
  }

  doe=interp(s_gp,d_o_e,SDIV,s_e);
  dge=interp(s_gp,d_g_e,SDIV,s_e);
  dre=interp(s_gp,d_r_e,SDIV,s_e);
  
  vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator)+
    sqrt(((dge+dre)/(8.0+dge-dre))+pow((doe/(8.0+dge-dre))*
				       r_e*exp(-rho_equator),2.0));

  if (scaled_polytrope==false) {
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  } else {
    Omega_K=omega_equator+vek*exp(rho_equator)/r_e;
  }

  // Rest mass and angular momentum

  Mass_0=0.0;
  Mass=0.0;
  J=0.0;

  if (scaled_polytrope==false) {
    for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++) {
	if (energy[s][m]>e_surface)
	  rho_0[s][m]=n0_at_e(energy[s][m])*MB*C*C*KSCALE;
	else
	  rho_0[s][m]=0.0;
      } 
  } else {
    for(s=1;s<=SDIV;s++)
      for(m=1;m<=MDIV;m++)
	rho_0[s][m]=(energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
  }


  if (SMAX==1.0) {
    s_temp=SDIV-1;
  } else {
    s_temp=SDIV;
  }

  for(s=1;s<=s_temp;s++) {
    // initialize
    D_m[s]=0.0;           
    D_m_0[s]=0.0;
    D_J[s]=0.0;

    for(m=1;m<=MDIV-2;m+=2) {

      D_m[s]+=(1.0/(3.0*(MDIV-1)))*
	( exp(2.0*alpha[s][m]+gamma[s][m])*
	  (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
	   (1.0+velocity_sq[s][m]+
	    (2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
	     (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*omega[s][m]*
	    exp(-rho[s][m]))+2.0*pressure[s][m])
	  +4.0*exp(2.0*alpha[s][m+1]+gamma[s][m+1])*
	  (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
	   (1.0+velocity_sq[s][m+1]+
	    (2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
	     (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*omega[s][m+1]*
	    exp(-rho[s][m+1]))+2.0*pressure[s][m+1]) 
	  +exp(2.0*alpha[s][m+2]+gamma[s][m+2])*
	  (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
	   (1.0+velocity_sq[s][m+2]+
	    (2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
	     (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*omega[s][m+2]*
	    exp(-rho[s][m+2]))+2.0*pressure[s][m+2]));    

      D_m_0[s]+=(1.0/(3.0*(MDIV-1)))*
	(exp(2.0*alpha[s][m]+
	     (gamma[s][m]
	      -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])
	 +4.0* exp(2.0*alpha[s][m+1]+
		   (gamma[s][m+1]
		    -rho[s][m+1])/2.0)*rho_0[s][m+1]/
	 sqrt(1.0-velocity_sq[s][m+1])
	 +exp(2.0*alpha[s][m+2]+(gamma[s][m+2]
				 -rho[s][m+2])/2.0)*rho_0[s][m+2]
	 /sqrt(1.0-velocity_sq[s][m+2])); 

      D_J[s]+=(1.0/(3.0*(MDIV-1)))*
	(sin_theta[m]*
	 exp(2.0*alpha[s][m]+gamma[s][m]-rho[s][m])*
	 (energy[s][m]+pressure[s][m])*sqrt(velocity_sq[s][m])/
	 (1.0-velocity_sq[s][m])
	 +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
	 exp(2.0*alpha[s][m+1]+gamma[s][m+1]-rho[s][m+1])
	 *(energy[s][m+1]+pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
	 (1.0-velocity_sq[s][m+1])
	 +sqrt(1.0-mu[m+2]*mu[m+2])*
	 exp(2.0*alpha[s][m+2]+gamma[s][m+2]-rho[s][m+2])
	 *(energy[s][m+2]+pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
	 (1.0-velocity_sq[s][m+2]));
    }
  }
    
  if (SMAX==1.0) {
    D_m[SDIV]=0.0;
    D_m_0[SDIV]=0.0;
    D_J[SDIV]=0.0;
  }

  for(s=1;s<=SDIV-4;s+=2) { 
    Mass+=(SMAX/(3.0*(SDIV-1)))*
      (pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
       D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
       +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);
      
    Mass_0+=(SMAX/(3.0*(SDIV-1)))*
      (pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
       D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
       +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
      
    J+=(SMAX/(3.0*(SDIV-1)))*
      ((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
       D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
       D_J[s+1]+(pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
       D_J[s+2]);
  }
    
  if (scaled_polytrope==false) {
    Mass*=4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G; 
    Mass_0*=4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
    J*=4*PI*KAPPA*C*C*C*pow(r_e,4.0)/G;
  } else {
    Mass*=4*PI*pow(r_e,3.0); 
    Mass_0*=4*PI*pow(r_e,3.0);
    J*=4*PI*pow(r_e,4.0);
  }

  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      gamma_guess[s][m]=gamma[s][m];
      rho_guess[s][m]=rho[s][m];
      alpha_guess[s][m]=alpha[s][m];
      omega_guess[s][m]=omega[s][m];
    }
  }

  r_e_guess=r_e;

  return;
}

void nstar_rot::comp() {

  // Loop indices
  int s;
  int m;
  int i;                

  // ---------------------------------------------------------------
  // Radius at pole

  r_p=r_ratio*r_e;                              

  // ---------------------------------------------------------------
  // The s-coordinate at pole and equator
  
  s_p=r_p/(r_p+r_e);                            
  s_e=0.5;

  // ---------------------------------------------------------------
  // Compute velocity over the grid

  double velocity[SDIV+1][MDIV+1];  
  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      velocity[s][m]=sqrt(velocity_sq[s][m]);
    }
  }

  // ---------------------------------------------------------------
  // Compute velocity_equator
 
  for(s=1;s<=SDIV;s++) {
    gamma_mu_1[s]=gamma[s][MDIV];                
    gamma_mu_0[s]=gamma[s][1];                   
    rho_mu_0[s]=rho[s][1];                     
    rho_mu_1[s]=rho[s][MDIV];                  
    omega_mu_0[s]=omega[s][1];                 
  }
  
  n_nearest=SDIV/2;

  // gamma at pole
  double gamma_pole=interp(s_gp,gamma_mu_1,SDIV,s_p);    
  // gamma at equator
  double gamma_equator=interp(s_gp,gamma_mu_0,SDIV,s_e); 
  // rho at pole
  double rho_pole=interp(s_gp,rho_mu_1,SDIV,s_p);      
  // rho at equator
  double rho_equator=interp(s_gp,rho_mu_0,SDIV,s_e);   

  // omega at equator
  double omega_equator;             

  if (r_ratio==1.0) {
    velocity_equator=0.0;
    omega_equator=0.0;
  } else {
    omega_equator=interp(s_gp,omega_mu_0,SDIV,s_e);
    velocity_equator=sqrt(1.0-exp(gamma_pole+rho_pole-gamma_equator
				  -rho_equator));
  }

  // ---------------------------------------------------------------
  // Circumferential radius

  if (scaled_polytrope==false) {
    R_e=sqrt(KAPPA)*r_e*exp((gamma_equator-rho_equator)/2.0);
  } else {
    R_e=r_e*exp((gamma_equator-rho_equator)/2.0);
  }

  // ---------------------------------------------------------------
  // Masses and angular momentum

  // initialize
  Mass=0.0;              
  Mass_0=0.0;
  Mass_p=0.0;
  J=0.0;

  // rest mass density
  double rho_0[SDIV+1][MDIV+1];     

  if (scaled_polytrope==false) {
    for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
	if (energy[s][m]>e_surface) {
	  rho_0[s][m]=n0_at_e(energy[s][m])*MB*C*C*KSCALE;
	} else {
	  rho_0[s][m]=0.0;
	}
      } 
    }
  } else {
    for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
	rho_0[s][m]=(energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
      }
    }
  }


  // int. quantity for M
  double D_m[SDIV+1];               
  // int. quantity for M_0
  double D_m_0[SDIV+1];             
  // int. quantity for M_p
  double D_m_p[SDIV+1];             
  // int. quantity for J
  double D_J[SDIV+1];               

  int s_temp;

  if (SMAX==1.0) s_temp=SDIV-1;
  else s_temp=SDIV;

  for(s=1;s<=s_temp;s++) {
    // initialize
    D_m[s]=0.0;           
    D_m_0[s]=0.0;
    D_m_p[s]=0.0;
    D_J[s]=0.0;

    for(m=1;m<=MDIV-2;m+=2) {
      D_m[s]+=(1.0/(3.0*(MDIV-1)))*
	(exp(2.0*alpha[s][m]+gamma[s][m])*
	 (((energy[s][m]+pressure[s][m])/(1.0-velocity_sq[s][m]))*
	  (1.0+velocity_sq[s][m]+
	   (2.0*s_gp[s]*sqrt(velocity_sq[s][m])/
	    (1.0-s_gp[s]))*sqrt(1.0-mu[m]*mu[m])*r_e*omega[s][m]*
	   exp(-rho[s][m]))+2.0*pressure[s][m])
	 +4.0*exp(2.0*alpha[s][m+1]+gamma[s][m+1])*
	 (((energy[s][m+1]+pressure[s][m+1])/(1.0-velocity_sq[s][m+1]))*
	  (1.0+velocity_sq[s][m+1]+
	   (2.0*s_gp[s]*sqrt(velocity_sq[s][m+1])/
	    (1.0-s_gp[s]))*sqrt(1.0-mu[m+1]*mu[m+1])*r_e*omega[s][m+1]*
	   exp(-rho[s][m+1]))+2.0*pressure[s][m+1]) 
	 +exp(2.0*alpha[s][m+2]+gamma[s][m+2])*
	 (((energy[s][m+2]+pressure[s][m+2])/(1.0-velocity_sq[s][m+2]))*
	  (1.0+velocity_sq[s][m+2]+
	   (2.0*s_gp[s]*sqrt(velocity_sq[s][m+2])/
	    (1.0-s_gp[s]))*sqrt(1.0-mu[m+2]*mu[m+2])*r_e*omega[s][m+2]*
	   exp(-rho[s][m+2]))+2.0*pressure[s][m+2]));    
      
      D_m_0[s]+=(1.0/(3.0*(MDIV-1)))*
	(exp(2.0*alpha[s][m]+
	     (gamma[s][m]
	      -rho[s][m])/2.0)*rho_0[s][m]/sqrt(1.0-velocity_sq[s][m])
	 +4.0* exp(2.0*alpha[s][m+1]+
		   (gamma[s][m+1]-rho[s][m+1])/2.0)*
	 rho_0[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
	 +exp(2.0*alpha[s][m+2]+
	      (gamma[s][m+2]-rho[s][m+2])/2.0)*
	 rho_0[s][m+2]/sqrt(1.0-velocity_sq[s][m+2])); 
      
      D_m_p[s]+=(1.0/(3.0*(MDIV-1)))*
	(exp(2.0*alpha[s][m]+(gamma[s][m]-rho[s][m])/2.0)*
	 energy[s][m]/sqrt(1.0-velocity_sq[s][m])
	 +4.0* exp(2.0*alpha[s][m+1]+
		   (gamma[s][m+1]-rho[s][m+1])/2.0)*
	 energy[s][m+1]/sqrt(1.0-velocity_sq[s][m+1])
	 +exp(2.0*alpha[s][m+2]+
	      (gamma[s][m+2]-rho[s][m+2])/2.0)*
	 energy[s][m+2]/sqrt(1.0-velocity_sq[s][m+2])); 
      
      D_J[s]+=(1.0/(3.0*(MDIV-1)))*
	(sqrt(1.0-mu[m]*mu[m])*
	 exp(2.0*alpha[s][m]+gamma[s][m]-rho[s][m])*
	 (energy[s][m]+pressure[s][m])*sqrt(velocity_sq[s][m])/
	 (1.0-velocity_sq[s][m])
	 +4.0*sqrt(1.0-mu[m+1]*mu[m+1])*
	 exp(2.0*alpha[s][m+1]+gamma[s][m+1]-rho[s][m+1])*
	 (energy[s][m+1]+pressure[s][m+1])*sqrt(velocity_sq[s][m+1])/
	 (1.0-velocity_sq[s][m+1])
	 +sqrt(1.0-mu[m+2]*mu[m+2])*
	 exp(2.0*alpha[s][m+2]+gamma[s][m+2]-rho[s][m+2])*
	 (energy[s][m+2]+pressure[s][m+2])*sqrt(velocity_sq[s][m+2])/
	 (1.0-velocity_sq[s][m+2]));
      
    }
  }

  if (SMAX==1.0) {
    D_m[SDIV]=0.0;
    D_m_0[SDIV]=0.0;
    D_m_p[SDIV]=0.0;
    D_J[SDIV]=0.0;
  }

  for(s=1;s<=SDIV-4;s+=2) { 

    Mass+=(SMAX/(3.0*(SDIV-1)))*
      (pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
       D_m[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m[s+1]
       +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m[s+2]);
      
    Mass_0+=(SMAX/(3.0*(SDIV-1)))*
      (pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
       D_m_0[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_0[s+1]
       +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_0[s+2]);
      
    Mass_p+=(SMAX/(3.0*(SDIV-1)))*
      (pow(sqrt(s_gp[s])/(1.0-s_gp[s]),4.0)*
       D_m_p[s]+4.0*pow(sqrt(s_gp[s+1])/(1.0-s_gp[s+1]),4.0)*D_m_p[s+1]
       +pow(sqrt(s_gp[s+2])/(1.0-s_gp[s+2]),4.0)*D_m_p[s+2]);
      
    J+=(SMAX/(3.0*(SDIV-1)))*
      ((pow(s_gp[s],3.0)/pow(1.0-s_gp[s],5.0))*
       D_J[s]+ 4.0*(pow(s_gp[s+1],3.0)/pow(1.0-s_gp[s+1],5.0))*
       D_J[s+1]+(pow(s_gp[s+2],3.0)/pow(1.0-s_gp[s+2],5.0))*
       D_J[s+2]);
  }
   
  if (scaled_polytrope==false) {
    Mass*=4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
    Mass_0*=4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
    Mass_p*=4*PI*sqrt(KAPPA)*C*C*pow(r_e,3.0)/G;
  } else {
    Mass*=4*PI*pow(r_e,3.0);
    Mass_0*=4*PI*pow(r_e,3.0);
    Mass_p*=4*PI*pow(r_e,3.0);
  }
 
  if (r_ratio==1.0) {
    J=0.0; 
    Omega=0.0;
  } else {    
    if (scaled_polytrope==false) {
      Omega=Omega_h*C/(r_e*sqrt(KAPPA));
      J*=4*PI*KAPPA*C*C*C*pow(r_e,4.0)/G;
    } else {
      Omega=Omega_h/r_e;
      J*=4*PI*pow(r_e,4.0);
    }
  }

  T=0.5*J*Omega;

  if (r_ratio==1.0) I=0.0;
  else I=J/Omega;
 
  if (scaled_polytrope==false) {
    W=Mass_p*C*C-Mass*C*C+T;
  } else {
    W=Mass_p-Mass+T;
  }

  // ---------------------------------------------------------------
  // Redshifts

  Z_p=exp(-0.5*(gamma_pole+rho_pole))-1.0;

  Z_b=sqrt((1.0+velocity_equator)/(1.0-velocity_equator))*
    (exp(-0.5*(gamma_equator+rho_equator))/
     (1.0-omega_equator*r_e*exp(-rho_equator)))-1.0;
  
  Z_f=sqrt((1.0-velocity_equator)/(1.0+velocity_equator))*
    (exp(-0.5*(gamma_equator+rho_equator))/(1.0+omega_equator*r_e
					    *exp(-rho_equator)))-1.0;

  // ---------------------------------------------------------------
  // Kepler angular velocity

  double d_o_e[SDIV+1];
  double d_g_e[SDIV+1];
  double d_r_e[SDIV+1];
  double d_v_e[SDIV+1];

  for(s=1;s<=SDIV;s++) { 
    d_o_e[s]=deriv_s(omega,s,1);
    d_g_e[s]=deriv_s(gamma,s,1);
    d_r_e[s]=deriv_s(rho,s,1);
    d_v_e[s]=deriv_s(velocity,s,1);
  }

  double doe=interp(s_gp,d_o_e,SDIV,s_e);
  double dge=interp(s_gp,d_g_e,SDIV,s_e);
  double dre=interp(s_gp,d_r_e,SDIV,s_e);
  double dve=interp(s_gp,d_v_e,SDIV,s_e);

  double vek=(doe/(8.0+dge-dre))*r_e*exp(-rho_equator)+
    sqrt(((dge+dre)/(8.0+dge-dre))+pow((doe/(8.0+dge-dre))*
				       r_e*exp(-rho_equator),2.0));
  
  if (scaled_polytrope==false) {
    Omega_K=(C/sqrt(KAPPA))*(omega_equator+vek*exp(rho_equator)/r_e);
  } else {
    Omega_K=omega_equator+vek*exp(rho_equator)/r_e;
  }

  // ---------------------------------------------------------------
  // Embedding
    
  double gamma_m[SDIV+1];
  double rho_m[SDIV+1];
  double alpha_m[SDIV+1];
  double enthalpy_m[SDIV+1];
  double s_surface[MDIV+1];
  double r_surface[MDIV+1];
  double gamma_surface[MDIV+1];
  double rho_surface[MDIV+1];
  double alpha_surface[MDIV+1]; 
  int s_surf_1=0;

  for(m=1;m<=MDIV;m++) { 

    /* Find last s grid point inside star. Construct arrays for the
       metric potentials and the enthalpy along each direction \mu.
    */
    for(s=1;s<=SDIV;s++) {
      if (energy[s][m]>0) s_surf_1=s;       
      gamma_m[s]=gamma[s][m];                 
      rho_m[s]=rho[s][m];                    
      alpha_m[s]=alpha[s][m];         
      enthalpy_m[s]=enthalpy[s][m];
    }
 
    /* If the positive enthalpy region outside the star is at a
       distance greater than two gridpoints, interpolate using two
       points on either side of the surface, else, if only one
       gridpoint has negative enthalpy, use linear interpolation, else
       use s_surf_1 as the surface.
    */
        
    if (enthalpy_m[s_surf_1+2]<0.0) {
      s_surface[m]=interp_4_k(enthalpy_m,s_gp,SDIV,0.0,s_surf_1-1);
    } else {     
      if (enthalpy_m[s_surf_1+1]<0.0) {
	s_surface[m]=s_gp[s_surf_1]-DS*enthalpy_m[s_surf_1]/
	  (enthalpy_m[s_surf_1+1]-enthalpy_m[s_surf_1]);
      } else {
	s_surface[m]=s_gp[s_surf_1];
      }
    }

    if (m==1) s_surface[m]=0.5;

    /* Construct arrays for the coordinate radius and the metric
       potentials on the surface.
    */
 
    r_surface[m]=r_e*(s_surface[m]/(1.0-s_surface[m]));

    gamma_surface[m]=interp(s_gp,gamma_m,SDIV,s_surface[m]);
    rho_surface[m]=interp(s_gp,rho_m,SDIV,s_surface[m]);
    alpha_surface[m]=interp(s_gp,alpha_m,SDIV,s_surface[m]);
  }

  double d_r_s_dm;
  double d_gamma_s_dm;
  double d_rho_s_dm;
  double f_mu[MDIV+1];
  double dpi_dm;

  for(m=1;m<=MDIV;m++) {
    d_gamma_s_dm=m_deriv(gamma_surface,m);
    d_rho_s_dm=m_deriv(rho_surface,m);
    d_r_s_dm=m_deriv(r_surface,m);

    if (m==MDIV) {
      f_mu[m]=0.0;
    } else {
      dpi_dm=exp(0.5*(gamma_surface[m]-rho_surface[m]))*
	(0.5*r_surface[m]*
	 (d_gamma_s_dm-d_rho_s_dm)*sin_theta[m]+d_r_s_dm*sin_theta[m] -
	 r_surface[m]*mu[m]/sin_theta[m]);
      f_mu[m]=sqrt(exp(2.0*alpha_surface[m])*
		   (pow(d_r_s_dm,2.0)
		    +pow(r_surface[m],2.0)/(1.0-mu[m]*mu[m]))-
		   pow(dpi_dm,2.0));
    }
  }

  double pi_bar[MDIV+1];
  double z_emb[MDIV+1];

  pi_bar[1]=R_e;
  // integrate f_mu
  z_emb[1]=0.0;                         
  for(m=2;m<=MDIV;m++) {                
    z_emb[m]=z_emb[m-1]+sqrt(KAPPA)*int_z(f_mu,m);    
    pi_bar[m]=sqrt(KAPPA)*exp((gamma_surface[m]-rho_surface[m])/2.0)
      *r_surface[m]*sin_theta[m];
    if (pi_bar[m]>pi_bar[m-1] && m>=2) {
      pi_bar[m]=pi_bar[m-1];
    }
  }  
  
  eccentricity=sqrt(1.0-pow(z_emb[MDIV]/R_e,2.0));
  
  // ---------------------------------------------------------------
  // Last stable circular orbit

  double s1;
  double s_1;
  double d_gamma_s;
  double d_rho_s;
  double d_omega_s;
  double d_gamma_m;
  double d_rho_m;
  double d_omega_m;
  double d_alpha_s;
  double d_alpha_m;
  double sqrt_v;

  for(s=1;s<=SDIV-1;s++) {
    s1=s_gp[s]*(1.0-s_gp[s]);
    s_1=1.0-s_gp[s];
    
    d_gamma_s=deriv_s(gamma,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    sqrt_v=exp(-2.0*rho[s][1])*r_e*r_e*pow(s_gp[s],4.0)*pow(d_omega_s,2.0) 
      +2*s1*(d_gamma_s+d_rho_s)+s1*s1*(d_gamma_s*d_gamma_s-d_rho_s*d_rho_s);

    if (sqrt_v>0.0) {
      sqrt_v=sqrt(sqrt_v);
    } else {
      sqrt_v=0.0;
      if (s_gp[s]>=s_e) {
	string str=((string)"Velocity imaginary at ")+
	  std::to_string(s)+" "+std::to_string(s_gp[s])+" " +
	  std::to_string(s_e)+" .";
	O2SCL_ERR(str.c_str(),o2scl::exc_efailed);
      }
    }

    v_plus[s]=(exp(-rho[s][1])*r_e*s_gp[s]*s_gp[s]*d_omega_s+sqrt_v)/
      (2.0+s1*(d_gamma_s-d_rho_s));

    v_minus[s]=(exp(-rho[s][1])*r_e*s_gp[s]*s_gp[s]*d_omega_s-sqrt_v)/
      (2.0+s1*(d_gamma_s-d_rho_s));
  }

  v_plus[SDIV]=0.0;
  v_minus[SDIV]=0.0;

  double B_st_p[SDIV+1];
  double B_st_m[SDIV+1];
  double d_v_plus_s;
  double d_v_minus_s;

  for(s=1;s<=SDIV;s++) {
    s1=s_gp[s]*(1.0-s_gp[s]);
        
    d_gamma_s=deriv_s(gamma,s,1);
    d_rho_s=deriv_s(rho,s,1);
    d_omega_s=deriv_s(omega,s,1);

    d_v_plus_s=s_deriv(v_plus,s);
    d_v_minus_s=s_deriv(v_minus,s);
 
    B_st_p[s]=v_plus[s]*(r_e*s_gp[s]*s_gp[s]*d_omega_s*exp(-rho[s][1])
			 +s1*d_v_plus_s)+0.5*s1*(d_gamma_s+d_rho_s)-
      pow(v_plus[s],4.0)*(1.0+0.5*s1*(d_gamma_s-d_rho_s));
    
    B_st_m[s]=v_minus[s]*(r_e*s_gp[s]*s_gp[s]*d_omega_s*exp(-rho[s][1])
			  +s1*d_v_minus_s)+0.5*s1*(d_gamma_s+d_rho_s)-
      pow(v_minus[s],4.0)*(1.0+0.5*s1*(d_gamma_s-d_rho_s));
  }  

  double B_st_p_surface;
  double B_st_m_surface;
  double B_st_p_out[SDIV-(SDIV/2)+2+1];
  double B_st_m_out[SDIV-(SDIV/2)+2+1];
  double s_plus;
  double s_minus;
  double r_plus=0.0;
  double gamma_plus=0.0;
  double rho_plus=0.0;
  double gamma_minus;
  double rho_minus;
  double omega_plus;
  double s_gp_out[SDIV-(SDIV/2)+2+1];

  B_st_p_surface=interp(s_gp,B_st_p,SDIV,s_surface[1]);
  B_st_m_surface=interp(s_gp,B_st_m,SDIV,s_surface[1]);

  if (B_st_p_surface>0.0) {
    h_plus=0.0;
    vel_plus=vek;
    Omega_plus=Omega_K;
  } else {
    for(i=1;i<=SDIV-(SDIV/2)+2;i++) {
      B_st_p_out[i]=B_st_p[(SDIV/2)-2+i];
      s_gp_out[i]=s_gp[(SDIV/2)-2+i];
    }

    n_nearest=SDIV/4;
    s_plus=interp(B_st_p_out, s_gp_out, SDIV-(SDIV/2)+2, 0.0);
    r_plus=r_e*s_plus/(1.0-s_plus);
    gamma_plus=interp(s_gp,gamma_mu_0,SDIV,s_plus);
    rho_plus=interp(s_gp,rho_mu_0,SDIV,s_plus);
    omega_plus=interp(s_gp,omega_mu_0,SDIV,s_plus);
    vel_plus=interp(s_gp,v_plus,SDIV,s_plus); 
 
    if (scaled_polytrope==false) {
      h_plus=sqrt(KAPPA)*(r_plus*exp((gamma_plus-rho_plus)/2.0)
			  -r_e*exp((gamma_equator-rho_equator)/2.0));
      Omega_plus=(C/sqrt(KAPPA))*(omega_plus+vel_plus*exp(rho_plus)/r_e);
    } else {
      h_plus=(r_plus*exp((gamma_plus-rho_plus)/2.0)
	      -r_e*exp((gamma_equator-rho_equator)/2.0));
      Omega_plus=(omega_plus+vel_plus*exp(rho_plus)/r_e);
    }

  }

  if (B_st_m_surface>0.0) {
    h_minus=0.0;
    vel_minus=vek;
  } else {
    for(i=1;i<=SDIV-(SDIV/2)+2;i++) {
      B_st_m_out[i]=B_st_m[(SDIV/2)-2+i];
      s_gp_out[i]=s_gp[(SDIV/2)-2+i];
    }
    
    n_nearest=SDIV/4;
    s_minus=interp(B_st_m_out,s_gp_out,SDIV-(SDIV/2)+2,0.0);
    gamma_minus=interp(s_gp,gamma_mu_0,SDIV,s_minus);
    rho_minus=interp(s_gp,rho_mu_0,SDIV,s_minus);
    vel_minus=interp(s_gp,v_plus,SDIV,s_minus); 
 
    if (scaled_polytrope==false) {
      h_minus=sqrt(KAPPA)*r_e*((s_minus/(1.0-s_minus))*
			       exp((gamma_minus-rho_minus)/2.0)-
			       exp((gamma_equator-rho_equator)/2.0));      
    } else {
      h_minus=r_e*((s_minus/(1.0-s_minus))*
		   exp((gamma_minus-rho_minus)/2.0)-
		   exp((gamma_equator-rho_equator)/2.0));
    }
  }

  double vel_p;

  if (h_plus!= 0.0) {

    u_phi=vel_plus*r_plus*exp((gamma_plus-rho_plus)/2.0)/
      sqrt(1.0-pow(vel_plus,2.0));
    vel_p=u_phi/sqrt(pow(u_phi,2.0)+pow(r_e,2.0)*
		     exp(gamma_equator-rho_equator));
    
    Omega_p=omega_equator+(vel_p/r_e)*exp(rho_equator);

    if (scaled_polytrope==false) {
      Omega_p*=C/sqrt(KAPPA);
    }

  } else {

    u_phi=vek*r_e*exp((gamma_equator-rho_equator)/2.0)/
      sqrt(1.0-pow(vek,2.0));
    Omega_p= Omega_K;

  }  

  if (scaled_polytrope==false) {
    u_phi*=C*C*sqrt(KAPPA)/(G*MSUN);
  }

  // Virial theorem

  // GRV2 spherical

  double virial1;
  double virial2;
  double S_virial1[SDIV+1][MDIV+1];
  double S_virial2[SDIV+1][MDIV+1];
  double grv2_spherical;

  virial1=0.0;
  virial2=0.0; 

  if (r_ratio==1.0) {

    if (SMAX==1.0) s_temp=SDIV-1;
    else s_temp=SDIV;

    for(s=1;s<=s_temp;s++) {
      d_gamma_s=deriv_s(gamma,s,1);
      d_rho_s=deriv_s(rho,s,1);
 
      S_virial1[s][1]=8*pow(PI*r_e,2.0)*s_gp[s]*pressure[s][1]
	*exp(2.0*alpha[s][1])/pow(1.0-s_gp[s],3.0);

      S_virial2[s][1]=PI*s_1_s[s]*pow(d_gamma_s+d_rho_s,2.0)/4.0;
    }
 
 
    if (SMAX==1.0) { 
      S_virial1[SDIV][1] =0.0;
 
      d_gamma_s=deriv_s(gamma,SDIV,1);
      d_rho_s=deriv_s(rho,SDIV,1);
      S_virial2[SDIV][1] =PI*s_1_s[s]*pow(d_gamma_s+d_rho_s,2.0)/4.0;
    }

    for(s=1;s<=SDIV-2;s+=2) {
      virial1+=(DS/3.0)*(S_virial1[s][1]+4.0*S_virial1[s+1][1]
			 +S_virial1[s+2][1]);
      virial2+=(DS/3.0)*(S_virial2[s][1]+4.0*S_virial2[s+1][1]
			 +S_virial2[s+2][1]);
    }
    
    grv2_spherical=fabs(1.0-virial1/virial2); 
  }

  // GRV2

  double virial3;
  double S_virial3[SDIV+1][MDIV+1];

  virial1=0.0;
  virial2=0.0; 
  virial3=0.0;

  double temp_x[5];

  for(i=1;i<=4;i++) {
    temp_x[i]=mu[MDIV+1-i];
  }

  double temp9=cos(0.5*(theta[MDIV]+theta[MDIV-2]));
  double temp_y1[5];
  double temp_y2[5];
  double S_ad1[SDIV+1][4];
  double S_ad2[SDIV+1][4];
  double m1;

  // Set equal to zero at center

  // AWS 4/25/14: There was some code here in the original which
  // didn't do anything so it was removed

  for(s=2;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      s1=s_1_s[s];
      m1=one_m2[m]; 
      d_gamma_s=deriv_s(gamma,s,m);
      d_rho_s=deriv_s(rho,s,m);
      d_gamma_m=deriv_m(gamma,s,m);
      d_rho_m=deriv_m(rho,s,m);
      d_omega_s=deriv_s(omega,s,m);
      d_omega_m=deriv_m(omega,s,m);
      if (m==MDIV || m==MDIV-1) { 
	S_virial1[s][m]=0.0;
	S_virial3[s][m]=0.0;
      } else {
	if (SMAX==1 && s==SDIV) {
	  S_virial1[s][m]=0.0;
	} else {
	  S_virial1[s][m]=16*PI*pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))*
	    (pressure[s][m]+(energy[s][m]+pressure[s][m])*
	     (velocity_sq[s][m]/(1.0-velocity_sq[s][m])))*
	    exp(2.0*alpha[s][m])/sqrt(m1);
	}
	S_virial3[s][m]=(s1/(2.0*sqrt(m1)))*pow(d_gamma_s+d_rho_s,2.0);
      }

      if (SMAX==1.0 && s==SDIV) {
	S_virial2[s][m]=0.0;
      } else {
	S_virial2[s][m]=2.0*
	  ( sqrt(m1)*pow(d_gamma_m+d_rho_m,2.0)/s1 
	    -(3.0*pow(r_e,2.0)
	      *s_gp[s]/(4.0*pow(1.0-s_gp[s],3.0)))*sqrt(m1)
	    *exp(-2.0*rho[s][m])
	    *(pow(s1*d_omega_s,2.0)+m1*pow(d_omega_m,2.0)));
      }
    }
      
    if (SMAX==1.0 && s==SDIV) {
      // AWS 4/25/14: This loop was modified because of array-indexing
      // confusion in the original
      for(m=1;m<4;m++) { 
	if (m!=2) {
	  S_ad1[s][m]=0.0;
	  S_ad2[s][m]=0.0;
	}
      }
      S_ad1[s][2]=0.0;
      S_ad2[s][2]=0.0;
    } else {
      for(m=1;m<=4;m++) { 
	temp_y1[m]=16*PI*pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))*
	  (pressure[s][MDIV+1-m]+
	   (energy[s][MDIV+1-m] 
	    +pressure[s][MDIV+1-m])*
	   (velocity_sq[s][MDIV+1-m]/
	    (1.0-velocity_sq[s][MDIV+1-m])))*exp(2.0*alpha[s][MDIV+1-m]);
	  
	temp_y2[m]=0.5*s1*pow(deriv_s(gamma,s,MDIV+1-m)+
			      deriv_s(rho,s,MDIV+1-m),2.0);
	
	// AWS 4/25/14: This conditional was modified from the 
	// original because of array-indexing confusion
	if (m!=2 && m!=4) {
	  S_ad1[s][m]=temp_y1[m];
	  S_ad2[s][m]=temp_y2[m];
	}            
      }  

      double temp1=interp(temp_x,temp_y1,4,temp9);
      double temp2=interp(temp_x,temp_y2,4,temp9);

      S_ad1[s][2]=temp1;
      S_ad2[s][2]=temp2;
    }

  }

  double D_virial1[SDIV+1];
  double D_virial2[SDIV+1];

  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV-4;m+=2) {
      virial1+=(DM/3.0)*(S_virial1[s][m]+4.0*S_virial1[s][m+1]
			 +S_virial1[s][m+2]);
      virial2+=(DM/3.0)*(S_virial2[s][m]+4.0*S_virial2[s][m+1]
			 +S_virial2[s][m+2]);    
      virial3+=(DM/3.0)*(S_virial3[s][m]+4.0*S_virial3[s][m+1]
			 +S_virial3[s][m+2]);    

    }

    virial1+=((theta[MDIV-2]-theta[MDIV])/6.0)*
      (S_ad1[s][1]+4.0*S_ad1[s][2]+S_ad1[s][3]);
    
    virial2+=(DM/3.0)*(S_virial2[s][MDIV-2]+4.0*S_virial2[s][MDIV-1]
		       +S_virial2[s][MDIV]);    
    virial3+=((theta[MDIV-2]-theta[MDIV])/6.0)*
      (S_ad2[s][1]+4.0*S_ad2[s][2]+S_ad2[s][3]);
    
    virial2+=virial3;

    D_virial1[s]=virial1;
    D_virial2[s]=virial2;
    virial1=0.0;
    virial2=0.0;
    virial3=0.0;
  }

  for(s=1;s<=SDIV-2;s+=2) {
    virial1+=(DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
		       +D_virial1[s+2]);
    virial2+=(DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
		       +D_virial2[s+2]);
  }

  grv2=fabs(1.0-virial1/virial2); 

  // GRV2 GAUSS-CHEBYSHEV

  double t_rho[MDIV+1];
  double t_alpha[MDIV+1];
  double t_rho_s[MDIV+1];
  double t_gamma_s[MDIV+1];
  double t_gamma_m[MDIV+1];
  double t_rho_m[MDIV+1];
  double t_omega_s[MDIV+1];
  double t_omega_m[MDIV+1];
  double muCh[MDIV+1];
  double t_pressure[MDIV+1];
  double t_energy[MDIV+1];
  double Sv1Ch[SDIV+1][MDIV+1];
  double Sv2Ch[SDIV+1][MDIV+1];
  double Dv1Ch[SDIV+1];
  double Dv2Ch[SDIV+1];
  double t_v2[MDIV+1];

  virial1=0.0;
  virial2=0.0;

  for(i=1;i<=MDIV;i++) muCh[i]=cos((1.0*i-0.5)/((2.0*MDIV-1.0))*PI);

  for(s=1;s<=s_temp;s++) {
    for(m=1;m<=MDIV;m++) {
      t_rho[m]=rho[s][m];
      t_alpha[m]=alpha[s][m];
      t_rho_s[m]=deriv_s(rho,s,m);
      t_gamma_s[m]=deriv_s(gamma,s,m);
      t_gamma_m[m]=deriv_m(gamma,s,m);
      t_rho_m[m]=deriv_m(rho,s,m);
      t_omega_s[m]=deriv_s(omega,s,m);
      t_omega_m[m]=deriv_m(omega,s,m);
      t_pressure[m]=pressure[s][m];
      t_energy[m]=energy[s][m];
      t_v2[m]=velocity_sq[s][m];
    }
    
    double rhoCh;
    double alphaCh;
    double rho_sCh;
    double gamma_sCh;
    double gamma_mCh;
    double rho_mCh;
    double omega_sCh;
    double omega_mCh;
    double pressureCh;
    double energyCh;
    double v2Ch;
    
    for(m=1;m<=MDIV;m++) {
      rhoCh=interp(mu,t_rho,MDIV,muCh[m]);
      alphaCh=interp(mu,t_alpha,MDIV,muCh[m]);
      rho_sCh=interp(mu,t_rho_s,MDIV,muCh[m]);
      gamma_sCh=interp(mu,t_gamma_s,MDIV,muCh[m]);
      gamma_mCh=interp(mu,t_gamma_m,MDIV,muCh[m]);
      rho_mCh=interp(mu,t_rho_m,MDIV,muCh[m]);
      omega_sCh=interp(mu,t_omega_s,MDIV,muCh[m]);
      omega_mCh=interp(mu,t_omega_m,MDIV,muCh[m]);
      pressureCh=interp(mu,t_pressure,MDIV,muCh[m]);
      energyCh=interp(mu,t_energy,MDIV,muCh[m]);
      v2Ch=interp(mu,t_v2,MDIV,muCh[m]);
      s1=s_1_s[s];
      m1=1.0-pow(muCh[m],2.0); 

      double temp1;
      if (s==1) temp1=0.0;
      else temp1=1.0/s1;
      
      Sv1Ch[s][m]= 8*PI*pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))*
	(pressureCh+(energyCh+pressureCh)*v2Ch/(1.0-v2Ch))*
	exp(2.0*alphaCh);
      
      Sv2Ch[s][m]=s1*pow(gamma_sCh+rho_sCh,2.0)/4.0+
	m1*pow(gamma_mCh+rho_mCh,2.0)*
	temp1-(3.0/4.0)*pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))*m1*
	exp(-2.0*rhoCh)*(pow(s1*omega_sCh,2.0)+m1*pow(omega_mCh,2.0));	 
    }


    if (SMAX==1.0) {
      for(m=1;m<=MDIV;m++) {    
	Sv1Ch[SDIV][m]=0.0;
	Sv2Ch[SDIV][m]=0.0;
      }
    }

    Dv1Ch[s]=0.0;
    Dv2Ch[s]=0.0;

    for(m=1;m<=MDIV;m++) {
      Dv1Ch[s]+=(2.0*PI/(2.0*MDIV-1.0))*Sv1Ch[s][m];
      Dv2Ch[s]+=(2.0*PI/(2.0*MDIV-1.0))*Sv2Ch[s][m];
    }
  }

  for(s=1;s<=SDIV-2;s+=2) {
    virial1+=(DS/3.0)*(Dv1Ch[s]+4.0*Dv1Ch[s+1]+Dv1Ch[s+2]);
    virial2+=(DS/3.0)*(Dv2Ch[s]+4.0*Dv2Ch[s+1]+Dv2Ch[s+2]);
  }

  grv2_new=fabs(1.0-virial1/virial2); 

  // GRV2 BY PARTS

  virial1=0.0;
  virial2=0.0; 

  /* Set equal to zero at center */
  /*
    for(m=1;m<=MDIV;m++) {
    S_virial1[1][m];
    S_virial2[1][m];
    Sv1Ch[1][m];
    Sv2Ch[1][m];
    }

    for(s=2;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
    s1=s_gp[s]*(1.0-s_gp[s]);
    m1=1.0-pow(mu[m],2.0); 
    d_gamma_s=deriv_s(gamma,s,m);
    d_rho_s=deriv_s(rho,s,m);
    d_gamma_m=deriv_m(gamma,s,m);
    d_rho_m=deriv_m(rho,s,m);
    d_omega_s=deriv_s(omega,s,m);
    d_omega_m=deriv_m(omega,s,m);
 
    S_virial1[s][m]=pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
    *( pressure[s][m]+(energy[s][m]+pressure[s][m])
    *(velocity_sq[s][m]/(1.0-velocity_sq[s][m])))
    *exp(2.0*alpha[s][m]);
 
    S_virial2[s][m]=(s1/4.0)*pow(d_gamma_s+d_rho_s,2.0)+
    m1*pow(d_gamma_m+d_rho_m,2.0)/s1-(3.0*pow(r_e,2.0)*s_gp[s]/
    (4.0*pow(1.0-s_gp[s],3.0)))*m1*exp(-2.0*rho[s][m])
    *( pow(s1*d_omega_s,2.0)+m1*pow(d_omega_m,2.0) ) ;
    }
    } 

    for(s=2;s<=SDIV;s++) 
    for(m=1;m<=MDIV;m++) {
    Sv1Ch[s][m]=deriv_m(S_virial1,s,m);
    Sv2Ch[s][m]=deriv_m(S_virial2,s,m);
    }

    D_virial1[1]=0.0;
    D_virial2[1]=0.0;

    for(s=2;s<=SDIV;s++) {
    s1=s_gp[s]*(1.0-s_gp[s]); 
    d_gamma_s=deriv_s(gamma,s,MDIV);
    d_rho_s=deriv_s(rho,s,MDIV);

    for(m=1;m<=MDIV-4;m+=2) {
    virial1+=(DM/3.0)*(asin(mu[m])*Sv1Ch[s][m]
    +4.0*asin(mu[m+1])*Sv1Ch[s][m+1]
    +asin(mu[m+2])*Sv1Ch[s][m+2]);

    virial2+=(DM/3.0)*(asin(mu[m])*Sv2Ch[s][m]
    +4.0*asin(mu[m+1])*Sv2Ch[s][m+1]
    +asin(mu[m+2])*Sv2Ch[s][m+2]);
    }

    virial1*=(-1.0);
    virial1+=(PI/2.0)*pow(r_e,2.0)*(s_gp[s]/pow(1.0-s_gp[s],3.0))
    *pressure[s][MDIV]*exp(2.0*alpha[s][MDIV]);       

    virial2*=(-1.0);
    virial2+=(PI/8.0)*s1*pow(d_gamma_s+d_rho_s,2.0);              
              
    D_virial1[s]=virial1;
    D_virial2[s]=virial2;
    virial1=0.0;
    virial2=0.0;
    }

    for(s=1;s<=SDIV-2;s+=2) {
    virial1+=(DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
    +D_virial1[s+2]);
    virial2+=(DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
    +D_virial2[s+2]);
    }

    virial1*=16.0*PI;
    virial2*=2.0;

    grv2_new=fabs(1.0-virial1/virial2); 

    printf("\n");
    printf("virial1=%6.5e, virial2=%6.5e\n",virial1,virial2);  
    printf("grv2 b.p. =%6.5e\n",grv2_new);
  */

  /* GRV3 */

  virial1=0.0;
  virial2=0.0;
 
  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      s1=s_1_s[s];
      m1=one_m2[m];
      d_gamma_s=deriv_s(gamma,s,m);
      d_gamma_m=deriv_m(gamma,s,m);
      d_rho_s=deriv_s(rho,s,m);
      d_rho_m=deriv_m(rho,s,m);
      d_omega_s=deriv_s(omega,s,m);
      d_omega_m=deriv_m(omega,s,m);
      d_alpha_s=deriv_s(alpha,s,m);
      d_alpha_m=deriv_m(alpha,s,m);


      if (SMAX==1.0 && s==SDIV) {
	S_virial1[s][m]=0.0;

	/* CORRECT THIS ! */
	
	S_virial2[s][m]=0.5*(pow(s_gp[s]*(d_gamma_s+d_rho_s),2.0)- d_alpha_s*
			     (pow(s_gp[s],2.0)*(d_gamma_s-d_rho_s)))*
	  exp(0.5*(gamma[s][m]-rho[s][m]))*r_e;
	
      } else {
	
	S_virial1[s][m]=8.0*PI*
	  (3.0*pressure[s][m]+(energy[s][m]+pressure[s][m])*
	   velocity_sq[s][m]/(1.0-velocity_sq[s][m]))*
	  exp(2.0*alpha[s][m]+0.5*(gamma[s][m]-rho[s][m]))*
	  pow(r_e,2.0)*r_e*pow(s_gp[s]/pow(1.0-s_gp[s],2.0),2.0);
	
	S_virial2[s][m]=0.5*
	  (pow(s1*(d_gamma_s+d_rho_s),2.0)+m1*pow(d_gamma_m+d_rho_m,2.0)
	   -d_alpha_s*(pow(s1,2.0)*(d_gamma_s-d_rho_s)+2.0*s1)
	   +d_alpha_m*(m1*(-d_gamma_m+d_rho_m)+2.0*mu[m]) 
	   +2.0*exp(2.0*alpha[s][m]-gamma[s][m]+rho[s][m])
	   *(s1*d_alpha_s-mu[m]*d_alpha_m)
	   +0.5*(1.0-exp(2.0*alpha[s][m]-gamma[s][m]+rho[s][m]))
	   *(s1*(d_gamma_s-d_rho_s)-mu[m]*(d_gamma_m-d_rho_m))
	   -1.5*exp(-2.0*rho[s][m])*
	   pow(r_e*s_gp[s]/(1.0-s_gp[s]),2.0)*m1*(pow(s1*d_omega_s,2.0)+
						  m1*pow(d_omega_m,2.0)))*
	  exp(0.5*(gamma[s][m]-rho[s][m]))*r_e/pow(1.0-s_gp[s],2.0);
      }
    }

  }
    
  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV-2;m+=2) {
      virial1+=(DM/3.0)*(S_virial1[s][m]+4.0*S_virial1[s][m+1]
			 +S_virial1[s][m+2]);
      virial2+=(DM/3.0)*(S_virial2[s][m]+4.0*S_virial2[s][m+1]
			 +S_virial2[s][m+2]);
    }     
    D_virial1[s]=virial1;
    D_virial2[s]=virial2;
    virial1=0.0;
    virial2=0.0;
  }

  for(s=1;s<=SDIV-2;s+=2) {
    virial1+=(DS/3.0)*(D_virial1[s]+4.0*D_virial1[s+1]
		       +D_virial1[s+2]);
    virial2+=(DS/3.0)*(D_virial2[s]+4.0*D_virial2[s+1]
		       +D_virial2[s+2]);
  }

  grv3=fabs(1.0-virial1/virial2);

  // ------------------------------------------------------
  // Code moved from print_table()

  if (r_ratio!=1.0) om_over_Om=omega[1][1]*r_e/Omega_h;
  else om_over_Om=0.0;

  // Calculate the mass quadrupole moment
  if (r_ratio!=1.0) {
    double r_infinity;
    // Value of the coordinate r at infinity
    r_infinity=r_e*s_gp[SDIV-1]/(1.0-s_gp[SDIV-1]);	
    mass_quadrupole=-0.5*pow(r_infinity,3)*D2_rho[SDIV-1][1];
  } else {
    mass_quadrupole=0.0;
  }

  // ------------------------------------------------------
  // Prepare next guess

  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      gamma_guess[s][m]=gamma[s][m];
      rho_guess[s][m]=rho[s][m];
      alpha_guess[s][m]=alpha[s][m];
      omega_guess[s][m]=omega[s][m];
    }
  }

  r_e_guess=r_e;
  
  return;
}

double nstar_rot::dm_dr_is(double r_is, double r, double m, double p) {
  double dmdr, e_d;

  if (p<p_surface) e_d=0.0;
  else e_d=e_at_p(p);
  
  if (r_is<RMIN) {
    dmdr=4.0*PI*e_center*r*r*(1.0+4.0*PI*e_center*r*r/3.0);
  } else {
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
  }
 
  return dmdr;
}
 
double nstar_rot::dp_dr_is(double r_is, double r, double m, double p) {
  double dpdr, e_d; 
    
  if (p<p_surface) e_d=0.0;
  else e_d=e_at_p(p);
  
  if (r_is<RMIN) {
    dpdr=-4.0*PI*(e_center+p)*(e_center+3.0*p)*r*
      (1.0+4.0*e_center*r*r/3.0)/3.0;
  } else {
    dpdr=-(e_d+p)*(m+4.0*PI*r*r*r*p)/(r*r_is*sqrt(1.0-2.0*m/r));
  }

  return dpdr;
}

double nstar_rot::dr_dr_is(double r_is, double r, double m) {
  double drdris;

  if (r_is<RMIN) drdris=1.0;
  else drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

  return drdris;
}

void nstar_rot::integrate(int i_check, double &r_final, double &m_final,
			   double &r_is_final) {
  int i=2;

  // radius 
  double r;                           
  // isotropic radial coordinate 
  double r_is;                        
  // mass   
  double m;                           
  // pressure 
  double p;                           
  // density 
  double e_d;                         
  // estimate on final isotr. radius
  double r_is_est;                    
  // r_is saving interval
  double dr_is_save;                  
  double r_is_check;
  double nu_s;
  double hh;
  // stepsize during integration 
  double h;                           
  // coeff. in Runge-Kutta equations 
  double a1, a2, a3, a4;
  double b1, b2, b3, b4;     
  double c1, c2, c3, c4;
  double rho0;   

  if (i_check==1) {
    if (scaled_polytrope==false) {
      r_is_est=1.5e6/sqrt(KAPPA);
    } else {
      r_is_est=2.0*sqrt(Gamma_P/(4.0*PI*(Gamma_P-1.0)))*
	pow(e_center,(Gamma_P-2.0)/2.0);
    }
    h=r_is_est/100;     
  } else {
    r_is_est=r_is_final;
    h=r_is_est/10000;     
  }

  dr_is_save=r_is_final/RDIV;
  r_is_check=dr_is_save;
  // initial isotropic radius 
  r_is=0.0;                            
  // initial radius 
  r=0.0;                               
  // initial mass
  m=0.0;                               
  // initial pressure
  p=p_center;                          

  r_is_gp[1]=0.0;
  r_gp[1]=0.0;
  m_gp[1]=0.0;
  lambda_gp[1]=0.0;
  e_d_gp[1]=e_center; 

  while (p>=p_surface) { 

    e_d=e_at_p(p);

    if ((i_check==3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      m_gp[i]=m;
      e_d_gp[i]=e_d; 
      i++;   
      r_is_check+=dr_is_save;
    }    
       
    r_is_final=r_is;
    r_final=r;
    m_final=m;
 
    a1=dr_dr_is(r_is,r,m);
    b1=dm_dr_is(r_is,r,m,p);
    c1=dp_dr_is(r_is,r,m,p);
 
    a2=dr_dr_is(r_is+h/2.0,r+h*a1/2.0,m+h*b1/2.0);
    b2=dm_dr_is(r_is+h/2.0,r+h*a1/2.0,m+h*b1/2.0,p+h*c1/2.0);
    c2=dp_dr_is(r_is+h/2.0,r+h*a1/2.0,m+h*b1/2.0,p+h*c1/2.0);

    a3=dr_dr_is(r_is+h/2.0,r+h*a2/2.0,m+h*b2/2.0);
    b3=dm_dr_is(r_is+h/2.0,r+h*a2/2.0,m+h*b2/2.0,p+h*c2/2.0);
    c3=dp_dr_is(r_is+h/2.0,r+h*a2/2.0,m+h*b2/2.0,p+h*c2/2.0);

    a4=dr_dr_is(r_is+h,r+h*a3,m+h*b3);
    b4=dm_dr_is(r_is+h,r+h*a3,m+h*b3,p+h*c3);
    c4=dp_dr_is(r_is+h,r+h*a3,m+h*b3,p+h*c3);

    r+=(h/6.0)*(a1+2*a2+2*a3+a4);
    m+=(h/6.0)*(b1+2*b2+2*b3+b4);
    p+=(h/6.0)*(c1+2*c2+2*c3+c4);

    r_is+=h;
  }

  r_is_gp[RDIV]=r_is_final;
  r_gp[RDIV]=r_final;
  m_gp[RDIV]=m_final;

  /* Rescale r_is and compute lambda */

  if (i_check==3) {
    double k_rescale=0.5*(r_final/r_is_final)*(1.0-m_final/r_final+
					       sqrt(1.0-2.0*m_final/r_final));
 
    r_is_final*=k_rescale;
 
    nu_s=log((1.0-m_final/(2.0*r_is_final))/(1.0+m_final/
					     (2.0*r_is_final)));

    for(i=1;i<=RDIV;i++) {
      r_is_gp[i]*=k_rescale;
 
      if (i==1) lambda_gp[1]=log(1.0/k_rescale);
      else lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); 

      if (e_d_gp[i]<e_surface) {
 
	hh=0.0;
 
      } else { 

	if (scaled_polytrope==false) {

	  p=p_at_e(e_d_gp[i]);
	  hh=h_at_p(p);

	} else { 
	  
	  rbc.tol_abs=1.0e-12;
	  polytrope_solve ps(Gamma_P,e_d_gp[i]);
	  rho0=0.0;
	  rbc.solve_bkt(rho0,e_d_gp[i],ps);

	  p=pow(rho0,Gamma_P);
	  hh=log((e_d_gp[i]+p)/rho0);

	}

      }
 
      nu_gp[i]=nu_s-hh;
    }
    nu_gp[RDIV]=nu_s;

  }

  return;
}

void nstar_rot::spherical_star() {
  int i;
  int s;
  int m;
  int s_temp;

  /** \brief */
  double r_final;
  /// Desc
  double m_final;
  /// Desc
  double r_is_final;

  r_final=0.0;                 
  m_final=0.0;
  r_is_final=0.0;

  double r_is_s;
  double lambda_s;
  double nu_s;
  double gamma_eq;
  double rho_eq;
  
  /* Solve the Oppenheimer-Volkov equations using Runge-Kutta. The
     function integrate() solves the equations using the r_is
     coordinate. */
  integrate(1,r_final,m_final,r_is_final);
  integrate(2,r_final,m_final,r_is_final);
  integrate(3,r_final,m_final,r_is_final);

  if (SMAX==1.0) s_temp=SDIV-1;
  else s_temp=SDIV;

  n_nearest=RDIV/2;

  for(s=1;s<=s_temp;s++) {
    r_is_s=r_is_final*(s_gp[s]/(1.0-s_gp[s]));
    // Convert the spherical solution to the 's' coordinte
    if (r_is_s<r_is_final) {
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s);
      nu_s=interp(r_is_gp,nu_gp,RDIV,r_is_s);
    } else {
      lambda_s=2.0*log(1.0+m_final/(2.0*r_is_s));
      nu_s=log((1.0-m_final/(2.0*r_is_s))/(1.0+m_final/(2*r_is_s)));
    }

    // Transform to the standard metric functions, gamma and rho
    gamma[s][1]=nu_s+lambda_s;
    rho[s][1]=nu_s-lambda_s;

    for(m=1;m<=MDIV;m++) {
      /* Since the solution is spherically symmetric, 
	 funct[m]=funct[1]=value of function on equatorial plane */
      gamma_guess[s][m]=gamma[s][1];
      rho_guess[s][m]=rho[s][1];
      alpha_guess[s][m]=(gamma[s][1]-rho[s][1])/2.0;
      omega_guess[s][m]=0.0; 
    }
 
    // gamma at \mu=0 
    gamma_mu_0[s]=gamma[s][1];                   
    // rho at \mu=0 
    rho_mu_0[s]=rho[s][1];                     
  }

  if (SMAX==1.0) {
    for(m=1;m<=MDIV;m++) {
      gamma_guess[SDIV][m]=0.0;
      rho_guess[SDIV][m]=0.0;
      alpha_guess[SDIV][m]=0.0;
      omega_guess[SDIV][m]=0.0; 
    }
 
    gamma_mu_0[SDIV]=0.0;
    rho_mu_0[SDIV]=0.0;
  }
   
  n_nearest=SDIV/2;
  // gamma at equator
  gamma_eq=interp(s_gp,gamma_mu_0,SDIV,s_e);      
  // rho at equator 
  rho_eq=interp(s_gp,rho_mu_0,SDIV,s_e);        
 
  r_e_guess=r_final*exp(0.5*(rho_eq-gamma_eq)); 

  R_e=r_final*sqrt(KAPPA);
  Mass=m_final*sqrt(KAPPA)*(C*C/G);
  Z_p=exp(-0.5*(gamma_eq+rho_eq))-1.0;

  return;
}

int nstar_rot::iterate(double r_ratio_loc) {
  int m;
  int s;
  int n;
  int k;
  int n_of_it=0;
  int s_temp;
  
  // Term in sum
  double s_term=0.0;          
  // Intermediate sum in eqn for rho
  double sum_rho=0.0;         
  // Intermediate sum in eqn for gamma
  double sum_gamma=0.0;        
  // Intermediate sum in eqn for omega
  double sum_omega=0.0;       
  // Equatorial radius in previous cycle
  double r_e_old;             
  // Difference | r_e_old -r_e |
  double r_e_diff=1.0;             
  // Derivative of gamma w.r.t. s
  double d_gamma_s;            
  // Derivative of gamma w.r.t. m
  double d_gamma_m;            
  // Derivative of rho w.r.t. s
  double d_rho_s;             
  // Derivative of rho w.r.t. m
  double d_rho_m;             
  // Derivative of omega w.r.t. s
  double d_omega_s;           
  // Derivative of omega w.r.t. m
  double d_omega_m;           
  // Second derivative of gamma w.r.t. s
  double d_gamma_ss;           
  // Second derivative of gamma w.r.t. m
  double d_gamma_mm;           
  // Derivative of gamma w.r.t. m and s
  double d_gamma_sm;           
  // Temporary term in da_dm
  double m1;                  
  double s1;
  double s2;
  double ea;
  double R_emb_eq;
  double R_emb_pole;
  double rsm;
  double gsm;
  double esm;
  double psm;
  double v2sm;
  double mum;
  double omsm;
  double sgp;
  double s_1;
  double e_gsm;
  double e_rsm; 
  double rho0sm;
 
  /// Desc
  double dgds[SDIV+1][MDIV+1];
  /// Desc
  double dgdm[SDIV+1][MDIV+1];

  if (SMAX==1.0) {
    s_temp=SDIV-1;
  } else {
    s_temp=SDIV;
  }
  
  /* The variable 'thing_guess' is the spherical solution for 'thing',
     or a previously computed star with lower angular momentum 
  */
  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      gamma[s][m]=gamma_guess[s][m];
      rho[s][m]=rho_guess[s][m];
      alpha[s][m]=alpha_guess[s][m];
      omega[s][m]=omega_guess[s][m];
    }
  }

  r_e=r_e_guess;

  while (r_e_diff>eq_radius_tol_rel || n_of_it<2) { 

    if (verbose>1) {
      cout << "r_e_diff, eq_radius_tol_rel: "
	   << r_e_diff << " " << eq_radius_tol_rel << endl;
    }
 
    /* Rescale potentials and construct arrays with the potentials
       along the equatorial and polar directions.
    */        

    for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
	rho[s][m]/=pow(r_e,2.0);
	gamma[s][m]/=pow(r_e,2.0); 
	alpha[s][m]/=pow(r_e,2.0);
	omega[s][m]*=r_e;
      }
      // The value of rho on the equatorial plane is rho_mu_0
      rho_mu_0[s]=rho[s][1];     
      gamma_mu_0[s]=gamma[s][1];   
      omega_mu_0[s]=omega[s][1]; 
      // Value of rho on the polar axis is rho_mu_1
      rho_mu_1[s]=rho[s][MDIV];  
      gamma_mu_1[s]=gamma[s][MDIV];
    }
 
    // Compute new r_e

    r_e_old=r_e;
    r_p=r_ratio_loc*r_e;                          
    s_p=r_p/(r_p+r_e);                        
  
    n_nearest=SDIV/2;
    gamma_pole_h=interp(s_gp,gamma_mu_1,SDIV,s_p); 
    gamma_equator_h=interp(s_gp,gamma_mu_0,SDIV,s_e);
    gamma_center_h=gamma[1][1];                    
  
    rho_pole_h=interp(s_gp,rho_mu_1,SDIV,s_p);   
    rho_equator_h=interp(s_gp,rho_mu_0,SDIV,s_e);
    rho_center_h=rho[1][1];                      
 
    r_e=sqrt(2*h_center/(gamma_pole_h+rho_pole_h-gamma_center_h-
			 rho_center_h));

    // Compute angular velocity Omega
 
    if (r_ratio_loc==1.0) {
      Omega_h=0.0;
      omega_equator_h=0.0;
    } else {
      // AWS: This looks like Eq. 46 from Cook '92
      omega_equator_h=interp(s_gp,omega_mu_0,SDIV,s_e);
      Omega_h=omega_equator_h+exp(pow(r_e,2.0)*rho_equator_h)*
	sqrt(1.0-exp(pow(r_e,2.0)*(gamma_pole_h+rho_pole_h-gamma_equator_h
				   -rho_equator_h)));
    }
 
    // Compute velocity, energy density and pressure
 
    for(s=1;s<=SDIV;s++) {
      sgp=s_gp[s];

      for(m=1;m<=MDIV;m++) {
	rsm=rho[s][m];
            
	if (r_ratio_loc==1.0 || s > (SDIV/2+2) ) {
	  velocity_sq[s][m]=0.0;
	} else {
	  velocity_sq[s][m]=pow((Omega_h-omega[s][m])*(sgp/(1.0-sgp))
				*sin_theta[m]*exp(-rsm*pow(r_e,2.0)),2.0);
	}
	
	if (velocity_sq[s][m]>=1.0) {
	  velocity_sq[s][m]=0.0;
	}

	enthalpy[s][m]=enthalpy_min+0.5*
	  (pow(r_e,2.0)*(gamma_pole_h+rho_pole_h
			 -gamma[s][m]-rsm)-log(1.0-velocity_sq[s][m]));
	
	if ((enthalpy[s][m]<=enthalpy_min) || (sgp>s_e)) {
	  pressure[s][m]=0.0;
	  energy[s][m]=0.0; 
	} else { 

	  if (scaled_polytrope==false) {
	    pressure[s][m]=p_at_h(enthalpy[s][m]);
	    energy[s][m]=e_at_p(pressure[s][m]);
	  } else {
	    rho0sm=pow(((Gamma_P-1.0)/Gamma_P)
		       *(exp(enthalpy[s][m])-1.0),1.0/(Gamma_P-1.0));
	    pressure[s][m]=pow(rho0sm,Gamma_P);
	    energy[s][m]=pressure[s][m]/(Gamma_P-1.0)+rho0sm;
	  }
	}  

	// Rescale back metric potentials (except omega)

	rho[s][m]*=pow(r_e,2.0);
	gamma[s][m]*=pow(r_e,2.0);
	alpha[s][m]*=pow(r_e,2.0);
      }
    }

    // Evaluation of source terms (Eqs. 30-33 of Cook, et al. (1992))

    if (SMAX==1.0) {
      s=SDIV;
      for(m=1;m<=MDIV;m++) {
	S_rho[s][m]=0.0;
	S_gamma[s][m]=0.0;
	S_omega[s][m]=0.0;
      }
    }

    for(s=1;s<=s_temp;s++) {
      for(m=1;m<=MDIV;m++) {
	rsm=rho[s][m];
	gsm=gamma[s][m];
	omsm=omega[s][m];
	esm=energy[s][m];
	psm=pressure[s][m];
	e_gsm=exp(0.5*gsm);
	e_rsm=exp(-rsm);
	v2sm=velocity_sq[s][m];
	mum=mu[m];            
	m1=1.0-pow(mum,2.0);
	sgp=s_gp[s];
	s_1=1.0-sgp;
	s1=sgp*s_1;
	s2=pow(sgp/s_1,2.0);  

	ea=16.0*PI*exp(2.0*alpha[s][m])*pow(r_e,2.0);
 
	if (s==1) {
	  d_gamma_s=0.0;
	  d_gamma_m=0.0;
	  d_rho_s=0.0;
	  d_rho_m=0.0;
	  d_omega_s=0.0;
	  d_omega_m=0.0;
	} else {
	  d_gamma_s=deriv_s(gamma,s,m);
	  d_gamma_m=deriv_m(gamma,s,m);
	  d_rho_s=deriv_s(rho,s,m);
	  d_rho_m=deriv_m(rho,s,m);
	  d_omega_s=deriv_s(omega,s,m);
	  d_omega_m=deriv_m(omega,s,m);
	}

	S_rho[s][m]=e_gsm*
	  (0.5*ea*(esm+psm)*s2*(1.0+v2sm)/(1.0-v2sm)+
	   s2*m1*pow(e_rsm,2.0)*
	   (pow(s1*d_omega_s,2.0)+m1*pow(d_omega_m,2.0))
	   +s1*d_gamma_s-mum*d_gamma_m+0.5*rsm*
	   (ea*psm*s2-s1*d_gamma_s*(0.5*s1*d_gamma_s+1.0)
	    -d_gamma_m*(0.5*m1*d_gamma_m-mum)));

	S_gamma[s][m]=e_gsm*
	  (ea*psm*s2+0.5*gsm*(ea*psm*s2-0.5*pow(s1*d_gamma_s,2.0)-
			      0.5*m1*pow(d_gamma_m,2.0)));
	
	S_omega[s][m]=e_gsm*e_rsm*
	  (-ea*(Omega_h-omsm)*(esm+psm)*s2/(1.0-v2sm)+
	   omsm*(-0.5*ea*(((1.0+v2sm)*esm+2.0*v2sm*psm)/(1.0-v2sm))*s2-
		 s1*(2*d_rho_s+0.5*d_gamma_s)
		 +mum*(2*d_rho_m+0.5*d_gamma_m)+0.25*pow(s1,2.0)*
		 (4*pow(d_rho_s,2.0)-pow(d_gamma_s,2.0))+0.25*m1*
		 (4*pow(d_rho_m,2.0)-pow(d_gamma_m,2.0))-
		 m1*pow(e_rsm,2.0)*(pow(pow(sgp,2.0)*d_omega_s,2.0)
				    +s2*m1*pow(d_omega_m,2.0))));
      }
    }
    
    // Angular integration (see Eqs. 27-29 of Cook, et al. (1992))
  
    n=0;
    for(k=1;k<=SDIV;k++) {      

      for(m=1;m<=MDIV-2;m+=2) {
	sum_rho+=(DM/3.0)*(P_2n[m][n]*S_rho[k][m]
			   +4.0*P_2n[m+1][n]*S_rho[k][m+1] 
			   +P_2n[m+2][n]*S_rho[k][m+2]);
      }

      D1_rho[n][k]=sum_rho;
      D1_gamma[n][k]=0.0;
      D1_omega[n][k]=0.0;
      sum_rho=0.0;

    }

    for(n=1;n<=LMAX;n++) {
      for(k=1;k<=SDIV;k++) {      
	for(m=1;m<=MDIV-2;m+=2) {

	  sum_rho+=(DM/3.0)*(P_2n[m][n]*S_rho[k][m]
			     +4.0*P_2n[m+1][n]*S_rho[k][m+1] 
			     +P_2n[m+2][n]*S_rho[k][m+2]);
                       
	  sum_gamma+=(DM/3.0)*(sin((2.0*n-1.0)*theta[m])*S_gamma[k][m]
			       +4.0*sin((2.0*n-1.0)*theta[m+1])*
			       S_gamma[k][m+1]
			       +sin((2.0*n-1.0)*theta[m+2])*S_gamma[k][m+2]);
  
	  sum_omega+=(DM/3.0)*(sin_theta[m]*P1_2n_1[m][n]*S_omega[k][m]
			       +4.0*sin_theta[m+1]*P1_2n_1[m+1][n]*
			       S_omega[k][m+1]
			       +sin_theta[m+2]*P1_2n_1[m+2][n]*
			       S_omega[k][m+2]);
	}
	D1_rho[n][k]=sum_rho;
	D1_gamma[n][k]=sum_gamma;
	D1_omega[n][k]=sum_omega;
	sum_rho=0.0;
	sum_gamma=0.0;
	sum_omega=0.0;
      }
    }

    // Radial integration

    n=0;
    for(s=1;s<=SDIV;s++) {
      for(k=1;k<=SDIV-2;k+=2) { 
	sum_rho+=(DS/3.0)*(f_rho[s][n][k]*D1_rho[n][k]+
			   4.0*f_rho[s][n][k+1]*D1_rho[n][k+1]+
			   f_rho[s][n][k+2]*D1_rho[n][k+2]);
      }
      D2_rho[s][n]=sum_rho;
      D2_gamma[s][n]=0.0;
      D2_omega[s][n]=0.0;
      sum_rho=0.0;
    }

    for(s=1;s<=SDIV;s++) {
      for(n=1;n<=LMAX;n++) {
	for(k=1;k<=SDIV-2;k+=2) { 
	  sum_rho+=(DS/3.0)*(f_rho[s][n][k]*D1_rho[n][k] 
			     +4.0*f_rho[s][n][k+1]*D1_rho[n][k+1]
			     +f_rho[s][n][k+2]*D1_rho[n][k+2]);
	  
	  sum_gamma+=(DS/3.0)*(f_gamma[s][n][k]*D1_gamma[n][k] 
			       + 4.0*f_gamma[s][n][k+1]*D1_gamma[n][k+1]
			       + f_gamma[s][n][k+2]*D1_gamma[n][k+2]);
	  
	  sum_omega+=(DS/3.0)*(f_omega[s][n][k]*D1_omega[n][k] 
			       +4.0*f_omega[s][n][k+1]*D1_omega[n][k+1]
			       +f_omega[s][n][k+2]*D1_omega[n][k+2]);
	}
	D2_rho[s][n]=sum_rho;
	D2_gamma[s][n]=sum_gamma;
	D2_omega[s][n]=sum_omega;
	sum_rho=0.0;
	sum_gamma=0.0;
	sum_omega=0.0;
      }
    }

    // Summation of coefficients

    for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {

	gsm=gamma[s][m];
	rsm=rho[s][m];
	omsm=omega[s][m];             
	e_gsm=exp(-0.5*gsm);
	e_rsm=exp(rsm);
	double temp1=sin_theta[m];

	sum_rho+=-e_gsm*P_2n[m][0]*D2_rho[s][0]; 

	for(n=1;n<=LMAX;n++) {

	  sum_rho+=-e_gsm*P_2n[m][n]*D2_rho[s][n]; 

	  if (m==MDIV) {             
	    sum_omega+=0.5*e_rsm*e_gsm*D2_omega[s][n]; 
	    sum_gamma+=-(2.0/PI)*e_gsm*D2_gamma[s][n];   
	  } else { 
	    sum_omega+=-e_rsm*e_gsm*(P1_2n_1[m][n]/
				     (2.0*n*(2.0*n-1.0)
				      *temp1))*D2_omega[s][n];  
	    sum_gamma+=-(2.0/PI)*e_gsm*(sin((2.0*n-1.0)*theta[m])
					/((2.0*n-1.0)*temp1))*
	      D2_gamma[s][n];   
	  }
	}
	   
	rho[s][m]=rsm+cf*(sum_rho-rsm);
	gamma[s][m]=gsm+cf*(sum_gamma-gsm);
	omega[s][m]=omsm+cf*(sum_omega-omsm);

	sum_omega=0.0;
	sum_rho=0.0;
	sum_gamma=0.0; 
      }
    }

    // Check for divergence

    if (fabs(omega[2][1])>100.0 || fabs(rho[2][1])>100.0 
	|| fabs(gamma[2][1])>300.0) {
      a_check=200; 
      break;
    }

    // Treat spherical case

    if (r_ratio_loc==1.0) {
      for(s=1;s<=SDIV;s++) {
	for(m=1;m<=MDIV;m++) {
	  rho[s][m]=rho[s][1];
	  gamma[s][m]=gamma[s][1];
	  omega[s][m]=0.0;          
	}
      }
    }
    
    // Treat infinity when SMAX=1.0
    
    if (SMAX==1.0) {
      for(m=1;m<=MDIV;m++) {
	rho[SDIV][m]=0.0;
	gamma[SDIV][m]=0.0;
	omega[SDIV][m]=0.0;
      }
    } 
      
    // Compute first order derivatives of gamma
 
    for(s=1;s<=SDIV;s++) {
      for(m=1;m<=MDIV;m++) {
	dgds[s][m]=deriv_s(gamma,s,m);
	dgdm[s][m]=deriv_m(gamma,s,m);
      }
    }

    // ALPHA (Integration of eq (39) of Cook, et al. (1992))
 
    if (r_ratio_loc==1.0) {
      for(s=1;s<=SDIV;s++) {
	for(m=1;m<=MDIV;m++) {
	  da_dm[s][m]=0.0;
	}
      }
    } else {
    
      for(s=2;s<=s_temp;s++) {
	for(m=1;m<=MDIV;m++) {

	  da_dm[1][m]=0.0; 
       
	  sgp=s_gp[s];
	  s1=sgp*(1.0-sgp);
	  mum=mu[m]; 
	  m1=1.0-pow(mum,2.0);
          
	  d_gamma_s=dgds[s][m];
	  d_gamma_m=dgdm[s][m];
	  d_rho_s=deriv_s(rho,s,m);
	  d_rho_m=deriv_m(rho,s,m);
	  d_omega_s=deriv_s(omega,s,m);
	  d_omega_m=deriv_m(omega,s,m);
	  d_gamma_ss=s1*deriv_s(dgds,s,m)+(1.0-2.0*sgp)*d_gamma_s;
	  d_gamma_mm=m1*deriv_m(dgdm,s,m)-2.0*mum*d_gamma_m;  
	  d_gamma_sm=deriv_sm(gamma,s,m);

	  double temp1=2.0*pow(sgp,2.0)*(sgp/(1.0-sgp))*m1*d_omega_s*d_omega_m
	    *(1.0+s1*d_gamma_s)-(pow(pow(sgp,2.0)*d_omega_s,2.0)-
				 pow(sgp*d_omega_m/(1.0-sgp),2.0)*m1)*
	    (-mum+m1*d_gamma_m); 
	  
	  double temp2=1.0/(m1*pow(1.0+s1*d_gamma_s,2.0)+
			    pow(-mum+m1*d_gamma_m,2.0));

	  double temp3=s1*d_gamma_ss+pow(s1*d_gamma_s,2.0);
  
	  double temp4=d_gamma_m*(-mum+m1*d_gamma_m);
	  
	  double temp5=(pow(s1*(d_rho_s+d_gamma_s),2.0)-
			m1*pow(d_rho_m+d_gamma_m,2.0))*(-mum+m1*d_gamma_m);

	  double temp6=s1*m1*(0.5*(d_rho_s+d_gamma_s)*(d_rho_m+d_gamma_m) 
			      +d_gamma_sm+d_gamma_s*d_gamma_m)*
	    (1.0+s1*d_gamma_s); 
	  
	  double temp7=s1*mum*d_gamma_s*(1.0+s1*d_gamma_s);
	  double temp8=m1*exp(-2*rho[s][m]);
	  
	  da_dm[s][m]=-0.5*(d_rho_m+d_gamma_m)-
	    temp2*(0.5*(temp3-d_gamma_mm-temp4)*(-mum+m1*d_gamma_m)+0.25*temp5-
		   temp6+temp7+0.25*temp8*temp1);	 
	}
      }
    }

    for(s=1;s<=s_temp;s++) {
      alpha[s][1]=0.0;
      for(m=1;m<=MDIV-1;m++) {
	alpha[s][m+1]=alpha[s][m]+0.5*DM*(da_dm[s][m+1]+
					  da_dm[s][m]);
      }
    } 
 
    for(s=1;s<=s_temp;s++) {
      for(m=1;m<=MDIV;m++) {     
	alpha[s][m]+=-alpha[s][MDIV]+0.5*(gamma[s][MDIV]-rho[s][MDIV]);
	
	if (alpha[s][m]>=300.0) {
	  a_check=200; 
	  break;
	}
	
	omega[s][m]/=r_e;
      } 
    }

    if (SMAX==1.0) {
      for(m=1;m<=MDIV;m++) {
	alpha[SDIV][m]=0.0;
      }
    }
    
    if (a_check==200) break;
    
    r_e_diff=fabs(r_e_old-r_e)/r_e;
    n_of_it++;

    // end of while (r_e_diff<eq_radius_tol_rel || n_of_it<2)
  }   

  for(s=1;s<=SDIV;s++) {
    for(m=1;m<=MDIV;m++) {
      gamma_guess[s][m]=gamma[s][m];
      rho_guess[s][m]=rho[s][m];
      alpha_guess[s][m]=alpha[s][m];
      omega_guess[s][m]=omega[s][m];
    }
  }
 
  r_e_guess=r_e;

  return 0;
} 

int nstar_rot::fix_cent_eden_axis_rat(double cent_eden, double axis_rat) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  
  int i;
  int j;

  double e_min;
 
  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    e_min=2.66e15*C*C*KSCALE;
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    e_min=0.34;
    r_ratio=0.58447265625;
    
  }

  e_min=cent_eden;
  r_ratio=axis_rat;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  e_center=e_min;
  make_center(e_center);
  spherical_star();
  
  if (r_ratio<0.8) {
    iterate(0.8);
  }
  
  if (r_ratio<0.6) {
    iterate(0.6);
  }
  
  iterate(r_ratio);
  comp();

  return 0;
}

int nstar_rot::fix_cent_eden_with_kepler(double cent_eden) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  e_center=e_min;

  /* First model is guess */
  
  make_center(e_center);
  spherical_star();             
  double d_Omega=1.0;           
  double sign=1.0;
  double dr=0.1;
  r_ratio=1.0;

  double diff_omega=1.0;
 
  /* Decrease r_p. Whenever Omega_K-Omega changes sign, or iteration does
     | not converge (a_check=200), cut the step dr in half and reverse 
     | its direction.
  */

  while((diff_omega>tol_abs || d_Omega<0.0) && r_ratio<=1.0) {          
    if (d_Omega*sign<0.0) { 
      sign=d_Omega;        
      dr/=(-2.0);              
    } 
    r_ratio -= dr; 
    a_check=0;
    iterate(r_ratio);
    if (a_check==200) {
      d_Omega=-1.0;
    }
    else {
      comp_omega();
      d_Omega=Omega_K-Omega;
      diff_omega=fabs(Omega-Omega_K)/Omega_K;
    }
    
    if (r_ratio>=1.0) {
      O2SCL_ERR("r_ratio>1.",o2scl::exc_efailed);
      //printf("r_ratio>=1.0 !\n");
    }
  }
  comp();
  
  return 0;
}

int nstar_rot::fix_cent_eden_non_rot(double cent_eden) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  r_ratio=1.0;
  e_center=e_min;
  make_center(e_center);
  spherical_star();
  iterate(r_ratio);
  comp();
  
  return 0;
}

int nstar_rot::fix_cent_eden_grav_mass(double cent_eden, double grav_mass) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  double M_fix=grav_mass;
  if (scaled_polytrope==false) {
    M_fix*=MSUN;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  e_center=e_min;
  
  {
    double dr;
    double diff_M;
    double d_ratio_M=0.0;

    dr=0.1;
    r_ratio=1.0-dr;

    /* Compute first rotating model */

    make_center(e_center);
    spherical_star();
    a_check=0;
    if (iterate(r_ratio)!=0) return 1;
    double sign;
    if (a_check==200) {
      diff_M=-1.0;
      sign=-1.0;
    } 
    else { 
      comp_M_J();
      diff_M=M_fix-Mass;
      sign=diff_M;
      d_ratio_M=fabs(diff_M)/M_fix;
    }

    /* If mass is greater than desired, reverse direction and cut stepsize
       | in half.
    */
 
    if (diff_M<0.0) {
      dr/=(-2.0);
    }

    while(d_ratio_M>tol_abs && r_ratio<=1.0) {
      if (diff_M*sign<0.0) {
	sign=diff_M;
	dr/=(-2.0);
      }
      r_ratio -= dr;
      a_check=0;
      if (iterate(r_ratio)!=0) return 2;
      if (a_check==200) {
	diff_M=-1.0;
      }
      else { 
	comp_M_J();      
	if (Omega>Omega_K) {
	  diff_M=-1.0;
	} else {
	  diff_M=M_fix-Mass;
	  d_ratio_M=fabs(diff_M)/M_fix;
	}     
      }
    } 
    comp();
  }

  return 0;
}

int nstar_rot::fix_cent_eden_bar_mass(double cent_eden, double bar_mass) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  double M_0const=bar_mass;
  if (scaled_polytrope==false) {
    M_0const*=MSUN;
  }
  e_center=e_min;

  double M_0=M_0const;
  {
    double dr;
    double diff_M_0;
    double d_ratio_M_0=0.0;

    dr=0.1;
    r_ratio=1.0-dr;

    /* Compute first rotating model */

    make_center(e_center);
    spherical_star();
    a_check=0;
    iterate(r_ratio);
    double sign;
    if (a_check==200) {
      diff_M_0=-1.0;
      sign=-1.0;
    } 
    else {   
      comp_M_J();
      diff_M_0=M_0-Mass_0;
      sign=diff_M_0;
      d_ratio_M_0=fabs(diff_M_0)/M_0;
    }

    /* If rest mass is greater than desired, reverse direction and cut stepsize
       | in half.
    */
 
    if (diff_M_0<0.0) {
      dr/=(-2.0);
    }

    while(d_ratio_M_0>tol_abs && r_ratio<=1.0) {
      if (diff_M_0*sign<0.0) {
	sign=diff_M_0;
	dr/=(-2.0);
      }
      r_ratio -= dr;
      a_check=0;
      iterate(r_ratio);
      if (a_check==200) {
	diff_M_0=-1.0;
      } 
      else { 
	comp_M_J();      
	if (Omega>Omega_K)
	  diff_M_0=-1.0;
	else {
	  diff_M_0=M_0-Mass_0;
	  d_ratio_M_0=fabs(diff_M_0)/M_0;
	}     
      }
    }
    comp();
  }

  return 0;
}

int nstar_rot::fix_cent_eden_ang_vel(double cent_eden, double ang_vel) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  double Omega_const=ang_vel;
  if (scaled_polytrope==false) {
    Omega_const*=1.0e4;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  e_center=e_min;
  {
    double dr;
    double diff_Omega;
    double d_ratio_Omega=0.0;

    dr=0.1;
    r_ratio=1.0-dr;

    /* Compute first rotating model */

    make_center(e_center);
    spherical_star();
    a_check=0;
    iterate(r_ratio);
    double sign;
    if (a_check==200) {
      diff_Omega=-1.0;
      sign=-1.0;
    } 
    else { 
      comp_omega();
      diff_Omega=Omega_const-Omega;
      sign=diff_Omega;
      d_ratio_Omega=fabs(diff_Omega)/Omega_const;
    }

    /* If Omega is greater than desired, reverse direction and cut stepsize
       | in half.
    */
 
    if (diff_Omega<0.0) {
      dr/=(-2.0);
    }

    while(d_ratio_Omega>tol_abs && r_ratio<=1.0) {
      if (diff_Omega*sign<0.0) {
	sign=diff_Omega;
	dr/=(-2.0);
      }
      r_ratio -= dr;
      a_check=0;
      iterate(r_ratio);
      if (a_check==200) {
	diff_Omega=-1.0;
      }
      else { 
	comp_omega();      
	if (Omega>Omega_K) {
	  diff_Omega=-1.0;
	} else {
	  diff_Omega=Omega_const-Omega;
	  d_ratio_Omega=fabs(diff_Omega)/Omega_const;
	}     
      }
    } 
    comp();
  }

  return 0;
}

int nstar_rot::fix_cent_eden_ang_mom(double cent_eden, double ang_mom) {

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "nstar_rot::fix_cent_eden...().",exc_einval);
  }
  int i;
  int j;

  if (scaled_polytrope==false) {
    e_surface=7.8*C*C*KSCALE;
    p_surface=1.01e8*KSCALE;
    enthalpy_min=1.0/(C*C);
  } else {
    e_surface=0.0;
    p_surface=0.0;
    enthalpy_min=0.0;
  }

  if (scaled_polytrope==false) {

    /* load the equation of state file */
    //load_eos(file_name);

    if (CL_LOW==true) {
      n0_match=1.0e38;
      e_cl=0.0; 
      de_pt=0.0;
    }

    /* set default values for star with tabulated eos */
    r_ratio=0.75;

  } else {

    /* set default values for polytropic star */
    r_ratio=0.58447265625;
    
  }

  double e_min=cent_eden;
  if (scaled_polytrope==false) {
    e_min*=C*C*KSCALE;
  }
  double J_const=ang_mom;
  if (scaled_polytrope==false) {
    J_const*=G*MSUN*MSUN/C;
  }
  
  Gamma_P=1.0+1.0/n_P;

  if (CL_LOW==true) {
    e_match=eosp->ed_from_nb(n0_match);
    p_match=eosp->pr_from_nb(n0_match);
    h_match=eosp->enth_from_nb(n0_match);
    //e_match=pow(10.0,interp(log_n0_tab,log_e_tab,n_tab,log10(n0_match)));
    //p_match=pow(10.0,interp(log_n0_tab,log_p_tab,n_tab,log10(n0_match)));
    //h_match=pow(10.0,interp(log_n0_tab,log_h_tab,n_tab,log10(n0_match)));
   
    if (e_cl != 0.0) de_pt=e_cl - e_match;   
  }

  e_center=e_min;

  {
    double dr;
    double diff_J;
    double d_ratio_J=0.0;

    dr=0.1;
    r_ratio=1.0-dr;

    /* Compute first rotating model */

    make_center(e_center);
    spherical_star();
    a_check=0;
    if (iterate(r_ratio)!=0) return 1;
    double sign;
    if (a_check==200) {
      diff_J=-1.0;
      sign=-1.0;
    } 
    else { 
      comp_M_J();
      diff_J=J_const-J;
      sign=diff_J;
      d_ratio_J=fabs(diff_J)/J_const;
    }

    /* If J is greater than desired, reverse direction and cut stepsize
       | in half.
    */
 
    if (diff_J<0.0) {
      dr/=(-2.0);
    }

    int new_ctr=0;

    while(d_ratio_J>tol_abs && r_ratio<=1.0) {
      if (diff_J*sign<0.0) {
	sign=diff_J;
	dr/=(-2.0);
      }
      r_ratio -= dr;
      a_check=0;
      if (iterate(r_ratio)!=0) return 2;
      if (a_check==200) {
	diff_J=-1.0;
      } else { 
	comp_M_J();      
	if (Omega>Omega_K) {
	  diff_J=-1.0;
	} else {
	  diff_J=J_const-J;
	  d_ratio_J=fabs(diff_J)/J_const;
	}     
      }
      if (new_ctr>100) return 3;
      new_ctr++;
    } 
    comp();
  }

  return 0;
}

void nstar_rot::test1(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_axis_rat(2.0e15,0.59);
  eos_set=false;

  t.test_rel(e_center,2.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,2.13324,2.0e-6,"2");
  t.test_rel(Mass_0/MSUN,2.43446,2.0e-6,"3");
  t.test_rel(R_e/1.0e5,13.9518,2.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.961702,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,1.00322,4.0e-6,"6");
  t.test_rel(T/W,0.109774,4.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,2.89653,2.0e-6,"8");
  t.test_rel(I/1.0e45,2.64698,2.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     277.074,2.0e-6,"10");
  t.test_rel(h_plus/1.0e5,0.129382,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,11.9949,4.0e-6,"12");
  t.test_rel(Z_p,0.591305,2.0e-6,"13");
  t.test_rel(Z_f,-0.296001,2.0e-6,"14");
  t.test_rel(Z_b,1.66307,2.0e-6,"15");
  t.test_rel(om_over_Om,0.722768,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,10.4129,4.0e-6,"17");
  t.test_rel(r_ratio,0.59,2.0e-6,"18");

  return;
}

void nstar_rot::test2(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_with_kepler(2.0e15);
  eos_set=false;

  t.test_rel(e_center,2.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,2.13633,2.0e-6,"2");
  t.test_rel(Mass_0/MSUN,2.43798,2.0e-6,"3");
  t.test_rel(R_e/1.0e5,14.3217,2.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.964219,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,0.964286,2.0e-6,"6");
  t.test_rel(T/W,0.110522,6.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,2.91467,2.0e-6,"8");
  t.test_rel(I/1.0e45,2.6566,2.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     279.481,2.0e-6,"10");
  t.test_abs(h_plus/1.0e5,0.0,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,11.6851,4.0e-6,"12");
  t.test_rel(Z_p,0.592893,2.0e-6,"13");
  t.test_rel(Z_f,-0.315912,2.0e-6,"14");
  t.test_rel(Z_b,1.67823,4.0e-6,"15");
  t.test_rel(om_over_Om,0.723265,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,10.7938,4.0e-6,"17");
  t.test_rel(r_ratio,0.568311,2.0e-6,"18");

  return;
}

void nstar_rot::test3(o2scl::test_mgr &t) {
  
  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_grav_mass(1.0e15,1.5);
  eos_set=false;

  t.test_rel(e_center,1.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,1.49996,4.0e-6,"2");
  t.test_rel(Mass_0/MSUN,1.64049,4.0e-6,"3");
  t.test_rel(R_e/1.0e5,13.7545,2.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.566135,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,0.87723,2.0e-6,"6");
  t.test_rel(T/W,0.0642316,2.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,1.14759,2.0e-6,"8");
  t.test_rel(I/1.0e45,1.78148,2.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     162.164,2.0e-6,"10");
  t.test_abs(h_plus/1.0e5,0.0,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,4.56495,2.0e-6,"12");
  t.test_rel(Z_p,0.272299,2.0e-6,"13");
  t.test_rel(Z_f,-0.117745,4.0e-6,"14");
  t.test_rel(Z_b,0.689716,2.0e-6,"15");
  t.test_rel(om_over_Om,0.475405,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,11.3694,4.0e-6,"17");
  t.test_rel(r_ratio,0.763672,2.0e-6,"18");

  return;
}

void nstar_rot::test4(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_bar_mass(1.0e15,1.55);
  eos_set=false;

  t.test_rel(e_center,1.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,1.41870,4.0e-6,"2");
  t.test_rel(Mass_0/MSUN,1.55003,4.0e-6,"3");
  t.test_rel(R_e/1.0e5,12.8888,4.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.439794,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,0.935245,2.0e-6,"6");
  t.test_rel(T/W,0.0364707,2.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,0.767165,2.0e-6,"8");
  t.test_rel(I/1.0e45,1.53303,4.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     86.1211,2.0e-6,"10");
  t.test_abs(h_plus/1.0e5,0.0,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,3.17552,2.0e-6,"12");
  t.test_rel(Z_p,0.247316,2.0e-6,"13");
  t.test_rel(Z_f,-0.0334667,2.0e-6,"14");
  t.test_rel(Z_b,0.542684,2.0e-6,"15");
  t.test_rel(om_over_Om,0.454775,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,10.6565,4.0e-6,"17");
  t.test_rel(r_ratio,0.863281,2.0e-6,"18");

  return;
}

void nstar_rot::test5(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_ang_vel(1.0e15,0.5);
  eos_set=false;

  t.test_rel(e_center,1.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,1.45222,4.0e-6,"2");
  t.test_rel(Mass_0/MSUN,1.58737,4.0e-6,"3");
  t.test_rel(R_e/1.0e5,13.2294,2.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.499991,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,0.912182,2.0e-6,"6");
  t.test_rel(T/W,0.0483527,2.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,0.929068,2.0e-6,"8");
  t.test_rel(I/1.0e45,1.63304,4.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     117.189,2.0e-6,"10");
  t.test_abs(h_plus/1.0e5,0.0,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,3.81883,2.0e-6,"12");
  t.test_rel(Z_p,0.257732,2.0e-6,"13");
  t.test_rel(Z_f,-0.0714265,2.0e-6,"14");
  t.test_rel(Z_b,0.606809,2.0e-6,"15");
  t.test_rel(om_over_Om,0.463524,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,10.9333,4.0e-6,"17");
  t.test_rel(r_ratio,0.820703,2.0e-6,"18");

  return;
}

void nstar_rot::test6(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_ang_mom(1.0e15,1.5);
  eos_set=false;

  t.test_rel(e_center,1.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,1.57929,2.0e-6,"2");
  t.test_rel(Mass_0/MSUN,1.72880,2.0e-6,"3");
  t.test_rel(R_e/1.0e5,14.8602,2.0e-6,"4");
  t.test_rel(Omega/1.0e4,0.644938,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,0.803770,4.0e-6,"6");
  t.test_rel(T/W,0.0882324,4.0e-6,"7");
  t.test_rel(J*C/G/MSUN/MSUN,1.50003,4.0e-6,"8");
  t.test_rel(I/1.0e45,2.04407,2.0e-6,"9");
  t.test_rel(mass_quadrupole*pow(sqrt(KAPPA),3.0)*C*C/G/1.0e42,
	     239.194,2.0e-6,"10");
  t.test_abs(h_plus/1.0e5,0.0,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,5.35325,4.0e-6,"12");
  t.test_rel(Z_p,0.295876,2.0e-6,"13");
  t.test_rel(Z_f,-0.188541,2.0e-6,"14");
  t.test_rel(Z_b,0.818747,2.0e-6,"15");
  t.test_rel(om_over_Om,0.493751,2.0e-6,"16");
  t.test_rel(r_e*sqrt(KAPPA)/1.0e5,12.3358,4.0e-6,"17");
  t.test_rel(r_ratio,0.670508,2.0e-6,"18");

  return;
}

void nstar_rot::test7(o2scl::test_mgr &t) {

  constants_rns();
  eos_nstar_rot_C p(true);
  set_eos(p);
  fix_cent_eden_non_rot(2.0e15);
  eos_set=false;

  t.test_rel(e_center,2.0,2.0e-6,"1");
  t.test_rel(Mass/MSUN,1.79249,2.0e-6,"2");
  t.test_rel(Mass_0/MSUN,2.04981,2.0e-6,"3");
  t.test_rel(R_e/1.0e5,10.7698,2.0e-6,"4");
  t.test_abs(Omega/1.0e4,0.0,2.0e-6,"5");
  t.test_rel(Omega_K/1.0e4,1.37951,2.0e-6,"6");
  t.test_abs(T/W,0.0,2.0e-6,"7");
  t.test_abs(J*C/G/MSUN/MSUN,0.0,2.0e-6,"8");
  t.test_rel(h_plus/1.0e5,5.10143,4.0e-6,"11");
  t.test_rel(h_minus/1.0e5,5.10143,2.0e-6,"12");
  t.test_rel(Z_p,0.401754,2.0e-6,"13");
  t.test_rel(Z_f,0.401754,2.0e-6,"14");
  t.test_rel(Z_b,0.401754,2.0e-6,"15");
  t.test_abs(om_over_Om,0.0,2.0e-6,"16");
  t.test_rel(r_e,0.215396,2.0e-6,"17");
  t.test_rel(r_ratio,1.0,2.0e-6,"18");
  t.test_rel(Omega_p/1.0e4,1.28671,4.0e-6,"19");
  t.test_rel(Omega_plus/1.0e4,1.27724,4.0e-6,"20");
  t.test_rel(u_phi,6.20982,2.0e-6,"21");

  return;
}

void nstar_rot::test8(o2scl::test_mgr &t) {

  polytrope_eos(1.0);
  fix_cent_eden_with_kepler(0.137);
  
  t.test_rel(e_center,0.137,2.0e-6,"1");
  t.test_rel(Mass,0.167269,2.0e-6,"2");
  t.test_rel(Mass_0,0.180264,2.0e-6,"3");
  t.test_rel(R_e,1.36735,2.0e-6,"4");
  t.test_rel(Omega,0.258135,2.0e-6,"5");
  t.test_rel(Omega_K,0.258149,2.0e-6,"6");
  t.test_rel(T/W,0.0924266,2.0e-6,"7");
  t.test_rel(J,0.0180245,2.0e-6,"8");
  t.test_rel(I,0.0698259,2.0e-6,"9");
  t.test_rel(mass_quadrupole,0.00971932,2.0e-6,"10");
  t.test_abs(h_plus,0.0,2.0e-6,"11");
  t.test_rel(h_minus,0.196584,4.0e-6,"12");
  t.test_rel(Z_p,0.250523,2.0e-6,"13");
  t.test_rel(Z_f,-0.247359,2.0e-6,"14");
  t.test_rel(Z_b,0.772557,2.0e-6,"15");
  t.test_rel(om_over_Om,0.459837,2.0e-6,"16");
  t.test_rel(r_e,1.18909,2.0e-6,"17");
  t.test_rel(r_ratio,0.575244,2.0e-6,"18");

  return;
}

