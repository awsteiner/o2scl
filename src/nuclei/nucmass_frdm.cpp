/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

#include <o2scl/nucmass_frdm.h>
// For unit conversions
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

nucmass_frdm::nucmass_frdm() {
  MH=7.289034;
  Mn=8.071431;
  e2=1.4399764;
  ael=1.433e-5;
  K=240.0;
  rp=0.8;
  r0=1.16;
  a=0.68;
  aden=0.7;
  rmac=4.8;
  h=6.6;
  W=30.0;
  L=0;
  a3=0;
  a1=16.247;
  a2=22.92;
  J=32.73;
  Q=29.21;
  a0=0.0;
  ca=0.436;
  C=60.0;
  gamma=0.831;
  amu=931.5014;
  //nfit=12;
  nfit=10;

  kg_to_invfm=o2scl_settings.get_convert_units().convert("kg","1/fm",1.0);

}

double nucmass_frdm::mass_excess_d(double Z, double N) {
  double ret;
      
  double A=Z+N;
  double cN=cbrt(N);
  double cZ=cbrt(Z);
  double cA=cbrt(A);
  double I=(N-Z)/A;
      
  // Ansatz for computing nuclear radius
  // This should imply a saturation density
  double R=r0*cA;

  double y=R/aden;
  double y0=r0*cA/aden;
  double x=R/a;
  double x0=r0*cA/a;

  double x02=x0*x0;
  B1=1.0-3.0/x02+(1.0+x0)*(2.0+3.0/x0+3.0/x02)*exp(-2.0*x0);
  B2=1.0-(1.0+2.0*x0-x02)*exp(-2.0*x0);
  double y02=y0*y0, y03=y02*y0, y04=y03*y0, y05=y04*y0;
  B3=1.0-5.0/y02*(1.0-15.0/8.0/y0+21.0/8.0/y03-0.75*
		  (1.0+4.5/y0+7.0/y02+3.5/y03)*exp(-2.0*y0));
  B4=1.0+5.0*(-3.0/y02+7.5/y03-63.0/4.0/pow(y0,5.0)
	      +0.75*(2.0/y0+12.0/y02+32.0/y03+
		     42.0/y04+21.0/y05)*exp(-2.0*y0));
      
  // Assume spherical shape
  Bs=1.0;
  Bv=1.0;
  Bw=1.0;
  Bk=1.0;
  Br=1.0;
      
  /* Pairing parameters:

     Deltan goes like N^{-1/3}
     Deltap goes like Z^{-1/3}
     deltanp goes like A^{-2/3}
  */
  Deltan=rmac*Bs/cN;
  Deltap=rmac*Bs/cZ;
  deltanp=h/Bs/cA/cA;

  // Coulomb coefficient
  c1=0.6*e2/r0;
  // Volume redistribution coefficient
  c2=1.0/336.0*(1.0/J+18.0/K)*c1*c1;
  // Coulomb exchange coefficient
  c4=1.25*pow(3.0/2.0/o2scl_const::pi,2.0/3.0)*c1;
  // Surface redistribution coefficient
  c5=c1*c1/64.0/Q;
  // Proton form-factor coefficient
  f0=-0.125*145.0/48.0*rp*rp*e2/r0/r0/r0;

  // Average bulk nuclear asymmetry, which
  // goes like \f$ (I,Z A^{-2/3} B1^{-1}) / (1,B1^{-1}) \f$.
  deltabar=(I+3.0/16.0*c1/Q*Z/cA/cA*Bv*Bs/B1)/
    (1.0+2.25*J/Q/cA*Bs*Bs/B1);

  // Average relative deviation of bulk density
  epsbar=(C*exp(-gamma*cA)-2.0*a2*B2/cA+L*deltabar*deltabar+
	  c1*Z*Z/A/cA*B4)/K;
      
  // Pairing contribution
  double tm=0.0, tm2;
  if (((int)Z)%2==1 && ((int)N)%2==1) {
    if (Z==N) tm=1.0/A;
    tm2=Deltap+Deltan-deltanp;
  } else if (((int)Z)%2==1 && ((int)N)%2==0) {
    tm2=Deltap;
  } else if (((int)N)%2==1) {
    tm2=Deltan;
  } else {
    tm2=0.0;
  }
      
  ret=MH*Z+Mn*N+(-a1+J*deltabar*deltabar-0.5*K*epsbar*epsbar)*A
    +(a2*B1+2.25*J*J/Q*deltabar*deltabar*Bs*Bs/B1)*cA*cA
    +a3*cA*Bk+c1*Z*Z/cA*B3-c2*Z*Z*cA*Br-c4*Z*cZ/cA
    -c5*Z*Z*Bw*Bs/B1+f0*Z*Z/A-ca*(N-Z)+W*(fabs(I)+tm)+tm2
    -ael*pow(Z,2.39);
      
  return ret;
}

int nucmass_frdm::fit_fun(size_t nv, const ubvector &x) {
  K=x[0]*200.0;
  r0=x[1];
  W=x[2];
  a1=x[3];
  a2=x[4];
  J=x[5];
  Q=x[6];
  ca=x[7];
  C=x[8]*60.0;
  gamma=x[9];
  //rmac=x[10];
  //h=x[11];
  return 0;

}

int nucmass_frdm::guess_fun(size_t nv, ubvector &x) {
  x[0]=K/200.0;
  x[1]=r0;
  x[2]=W;
  x[3]=a1;
  x[4]=a2;
  x[5]=J;
  x[6]=Q;
  x[7]=ca;
  x[8]=C/60.0;
  x[9]=gamma;
  //x[10]=rmac;
  //x[11]=h;
  return 0;

}

double nucmass_frdm::drip_binding_energy_d
(double Z, double N, double npout, double nnout, double chi) {
  
  double ret=(drip_mass_excess_d(Z,N,npout,nnout,chi)+
	      ((Z+N)*o2scl_mks::unified_atomic_mass-Z*o2scl_mks::mass_electron-
	       N*o2scl_mks::mass_neutron-Z*o2scl_mks::mass_proton)*
	      o2scl_const::hc_mev_fm*kg_to_invfm);
  return ret;
}

double nucmass_frdm::drip_mass_excess_d(double Z, double N,
					double np_out, double nn_out,
					double chi) {
  double ret;
      
  double dN=N;
  double dZ=Z;
  double dA=Z+N;
  double cN=cbrt(dN);
  double cZ=cbrt(dZ);
  double cA=cbrt(dA);
  double I=(dN-dZ)/dA;
      
  // Ansatz for computing nuclear radius
  // This should imply a saturation density
  double R=r0*cA;

  double y=R/aden;
  double y0=r0*cA/aden;
  double x=R/a;
  double x0=r0*cA/a;

  double x02=x0*x0;
  B1=1.0-3.0/x02+(1.0+x0)*(2.0+3.0/x0+3.0/x02)*exp(-2.0*x0);
  B2=1.0-(1.0+2.0*x0-x02)*exp(-2.0*x0);
  double y02=y0*y0, y03=y02*y0, y04=y03*y0, y05=y04*y0;
  B3=1.0-5.0/y02*(1.0-15.0/8.0/y0+21.0/8.0/y03-0.75*
		  (1.0+4.5/y0+7.0/y02+3.5/y03)*exp(-2.0*y0));
  B4=1.0+5.0*(-3.0/y02+7.5/y03-63.0/4.0/pow(y0,5.0)
	      +0.75*(2.0/y0+12.0/y02+32.0/y03+
		     42.0/y04+21.0/y05)*exp(-2.0*y0));
      
  // Assume spherical shape
  Bs=1.0;
  Bv=1.0;
  Bw=1.0;
  Bk=1.0;
  Br=1.0;
      
  /* Pairing parameters:

     Deltan goes like N^{-1/3}
     Deltap goes like Z^{-1/3}
     deltanp goes like A^{-2/3}
  */
  Deltan=rmac*Bs/cN;
  Deltap=rmac*Bs/cZ;
  deltanp=h/Bs/cA/cA;

  // Coulomb coefficient
  c1=0.6*e2/r0;
  // Volume redistribution coefficient
  c2=1.0/336.0*(1.0/J+18.0/K)*c1*c1;
  // Coulomb exchange coefficient
  c4=1.25*pow(3.0/2.0/o2scl_const::pi,2.0/3.0)*c1;
  // Surface redistribution coefficient
  c5=c1*c1/64.0/Q;
  // Proton form-factor coefficient
  f0=-0.125*145.0/48.0*rp*rp*e2/r0/r0/r0;

  // Average bulk nuclear asymmetry, which
  // goes like \f$ (I,Z A^{-2/3} B1^{-1}) / (1,B1^{-1}) \f$.
  deltabar=(I+3.0/16.0*c1/Q*dZ/cA/cA*Bv*Bs/B1)/
    (1.0+2.25*J/Q/cA*Bs*Bs/B1);

  // Average relative deviation of bulk density
  epsbar=(C*exp(-gamma*cA)-2.0*a2*B2/cA+L*deltabar*deltabar+
	  c1*dZ*dZ/dA/cA*B4)/K;
      
  // Saturation density
  double n0=0.152946;

  // Determine densities and radii
  nn=0.5*(1.0+deltabar)*(1.0-3.0*epsbar)*n0;
  np=0.5*(1.0-deltabar)*(1.0-3.0*epsbar)*n0;
      
  Rn=cbrt(3*N/2.0/o2scl_const::pi/(1.0+deltabar)/(1.0-3.0*epsbar)/n0);
  Rp=cbrt(3*Z/2.0/o2scl_const::pi/(1.0-deltabar)/(1.0-3.0*epsbar)/n0);
      
  double tm=0.0, tm2=0.0;
  if (false) {
    // Pairing contribution
    if (((int)Z)%2==1 && ((int)N)%2==1) {
      if (Z==N) tm=1.0/dA;
      tm2=Deltap+Deltan-deltanp;
    } else if (((int)Z)%2==1 && ((int)N)%2==0) {
      tm2=Deltap;
    } else if (((int)N)%2==1) {
      tm2=Deltan;
    } else {
      tm2=0.0;
    }
  }
      
  //double xp=np/(nn+np);
  double x3fact=1.0;//+(xp*xp*xp-1.0)/(1.0+exp((xp-0.25)/0.01));
      
  double surf_part=(a2*B1+2.25*J*J/Q*deltabar*deltabar*Bs*Bs/B1)*cA*cA+
    -c5*dZ*dZ*Bw*Bs/B1;
      
  ret=MH*dZ+Mn*dN+(-a1+J*deltabar*deltabar-0.5*K*epsbar*epsbar)*dA
    +a3*cA*Bk+c1*dZ*dZ/cA*B3-c2*dZ*dZ*cA*Br-c4*dZ*cZ/cA
    +f0*dZ*dZ/dA-ca*(dN-dZ)+W*(fabs(I)+tm)+tm2
    -ael*pow(Z,2.39)+surf_part*x3fact;
      
  return ret;
}

nucmass_mnmsk::nucmass_mnmsk() {
  n=0;
}

nucmass_mnmsk::~nucmass_mnmsk() {
  if (n>0) {
    delete[] mass;
  }
}

int nucmass_mnmsk::set_data(int n_mass, nucmass_mnmsk::entry *m, 
			    std::string ref) {
  n=n_mass;
  mass=m;
  reference=ref;
  last=n/2;
  return 0;
}

double nucmass_mnmsk::mass_excess(int Z, int N) {
  nucmass_mnmsk::entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.Mth;
}

bool nucmass_mnmsk::is_included(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      return false;
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    return true;
  }


  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {
    if (mass[mid].N==l_N) {
      return true;
    } else if (mass[mid].N>l_N) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==n-1) return false;
      mid++;
    }
  }
  
  return false;
}

bool nucmass_mnmsk_exp::is_included(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      return false;
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    if (fabs(mass[mid].Mexp)>1.0e-20 &&
	fabs(mass[mid].Mexp)<1.0e90) {
      return true;
    } else {
      return false;
    }
  }

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {
    if (mass[mid].N==l_N) {
      if (fabs(mass[mid].Mexp)>1.0e-20 &&
	  fabs(mass[mid].Mexp)<1.0e90) {
	return true;
      } else {
	return false;
      }
    } else if (mass[mid].N>l_N) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==n-1) return false;
      mid++;
    }
  }
  
  return false;
}

double nucmass_mnmsk_exp::mass_excess(int Z, int N) {
  nucmass_mnmsk::entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.Mexp;
}

nucmass_mnmsk::entry nucmass_mnmsk::get_ZN(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;
  
  nucmass_mnmsk::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  
  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      O2SCL_ERR((((string)"Nuclei with Z=")+itos(l_Z) 
		 +" not found in nucmass_mnmsk::get_ZN().").c_str(),
		exc_enotfound);
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    ret=mass[mid];
    last=mid;
    return ret;
  }

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {
    if (mass[mid].N==l_N) {
      ret=mass[mid];
      last=mid;
      return ret;
    } else if (mass[mid].N>l_N) {
      mid--;
    } else {
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)
	     +" not found in nucmass_mnmsk::get_ZN().").c_str(),exc_enotfound);
  return ret;
}

