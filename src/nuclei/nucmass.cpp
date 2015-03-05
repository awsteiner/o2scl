/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#include <o2scl/nucmass.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_info::nucmass_info() {
  
  std::string tlist[119]=
  // 0-10
    {"n","H","He","Li","Be","B","C","N","O","F","Ne",
     // 11-20
     "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
     // 21-30
     "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
     // 31-40
     "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
     // 41-50
     "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
     // 51-60
     "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
     // 61-70
     "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
     // 71-80
     "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
     // 81-90
     "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
     // 91-100
     "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
     // 101-110
     "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
     // 111-118
     "Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};

#ifdef O2SCL_NEVER_DEFINED
  std::string namelist[119]={
    "neutron",
    // 1-10
    "Hydrogen","Helium","Lithium","Beryllium","Boron",
    "Carbon","Nitrogen","Oxygen","Fluorine","Neon",
    // 11-20
    "Sodium","Magnesium","Aluminum","Silicon","Phosphorus",
    "Sulfur","Chlorine","Argon","Potassium","Calcium",
    // 21-30
    "Scandium","Titanium","Vanadium","Chromium","Manganese",
    "Iron","Cobalt","Nickel","Copper","Zinc",
    // 31-40
    "Gallium","Germanium","Arsenic","Selenium","Bromine",
    "Krypton","Rubidium","Strontium","Yttrium","Zirconium",
    // 41-50
    "Niobium","Molybdenum","Technetium","Ruthenium","Rhodium",
    "Palladium","Silver","Cadmium","Indium","Tin",
    // 51-60
    "Antimony","Tellurium","Iodine","Xenon","Cesium",
    "Barium","Lanthanum","Cerium","Praseodymium","Neodymium",
    // 61-70
    "Promethium","Samarium","Europium","Gadolinium","Terbium",
    "Dysprosium","Holmium","Erbium","Thulium","Ytterbium",
    // 71-80
    "Lutetium","Hafnium","Tantalum","Tungsten","Rhenium",
    "Osmium","Iridium","Platinum","Gold","Mercury",
    // 81-90
    "Thallium","Lead","Bismuth","Polonium","Astatine",
    "Radon","Francium","Radium","Actinium","Thorium",
    // 91-100
    "Protactinium","Uranium","Neptunium","Plutonium","Americium",
    "Curium","Berkelium","Californium","Einsteinium","Fermium",
    // 101-110
    "Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium",
    "Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium",
    // 111-118
    "Roentgenium","Copernicium","Ununtrium","Flerovium","Ununpentium",
    "Livermorium","Ununseptium","Ununoctium"};
#endif
    
  for(int i=0;i<nelements;i++) {
    element_list[i]=tlist[i];
    element_table.insert(make_pair(element_list[i],i));
  }

}

int nucmass_info::parse_elstring(std::string ela, int &Z, int &N, 
				      int &A) {

  bool removed=true;
  while (removed) {
    removed=false;
    for(size_t i=0;i<ela.size();i++) {
      if (!isalnum(ela[i])) {
	ela=ela.substr(0,i)+ela.substr(i+1,ela.size()-i-1);
	removed=true;
      }
    }
  }

  // If its just a neutron
  if (ela=="n") {
    Z=0;
    N=1;
    A=1;
    return 0;
  }
  // If its just a proton
  if (ela=="p") {
    Z=1;
    N=0;
    A=1;
    return 0;
  }
  // If its just a deuteron
  if (ela=="d") {
    Z=1;
    N=1;
    A=2;
    return 0;
  }
  // If its just a triton;
  if (ela=="t") {
    Z=1;
    N=2;
    A=3;
    return 0;
  }

  // If the number precedes the element, swap them
  {
    size_t inum=0, ichar=0;
    // Count number of digits
    while (inum<ela.length() && isdigit(ela[inum])) inum++;
    // Count number of letters following the digits
    while (ichar+inum<ela.length() && 
	   isalpha(ela[ichar+inum])) ichar++;
    if (inum>0 && ichar>0 && ichar<4) {
      string stemp;
      // Build a new element string with the characters first...
      for(size_t i=0;i<ichar;i++) {
	stemp+=ela[inum+i];
      }
      // and then the digits
      for(size_t i=0;i<inum;i++) {
	stemp+=ela[i];
      }
      ela=stemp;
    }
  }

  // Change the first character to upper case if necessary
  if (ela[0]>='a' && ela[0]<='z') {
    ela[0]='A'+(ela[0]-'a');
  }
  if (ela.length()<2) {
    O2SCL_ERR2("Element name too short in ",
		   "nucmass::parse_elstring().",exc_efailed);
  }
  if (ela.length()>3 && isalpha(ela[2])) {
    std::string el=ela.substr(0,3);
    Z=eltoZ(el);
    A=o2scl::stoi(ela.substr(3,ela.length()-2));
    N=A-Z;
    return 0;
  }
  if (ela.length()>2 && isalpha(ela[1])) {
    std::string el=ela.substr(0,2);
    Z=eltoZ(el);
    A=o2scl::stoi(ela.substr(2,ela.length()-2));
    N=A-Z;
    return 0;
  }
  std::string el=ela.substr(0,1);
  Z=eltoZ(el);
  A=o2scl::stoi(ela.substr(1,ela.length()-1));
  N=A-Z;
  return 0;
}

int nucmass_info::eltoZ(std::string el) {
  std::map<std::string,int,std::greater<std::string> >::iterator 
    eti=element_table.find(el);
  if (eti==element_table.end()) {
    O2SCL_ERR2("Failed to find element in ",
		   "nucmass_info::eltoZ().",-1);
  }
  return eti->second;
}

std::string nucmass_info::Ztoel(size_t Z) {
  if (((int)Z)>=nelements) {
    O2SCL_ERR2("Invalid element in ",
	       "nucmass_info::Ztoel().",exc_einval);
  }
  return element_list[Z];
}

std::string nucmass_info::tostring(size_t Z, size_t N) {
  if (((int)Z)>=nelements) {
    O2SCL_ERR2("Invalid element in ",
	       "nucmass_info::tostring().",exc_einval);
  }
  return element_list[Z]+itos(N+Z);
}

nucmass::nucmass() {
  m_neut=o2scl_mks::mass_neutron*
    o2scl_settings.get_convert_units().convert("kg","MeV",1.0);
  m_prot=o2scl_mks::mass_proton*
    o2scl_settings.get_convert_units().convert("kg","MeV",1.0);
  m_elec=o2scl_mks::mass_electron*
    o2scl_settings.get_convert_units().convert("kg","MeV",1.0);
  m_amu=o2scl_mks::unified_atomic_mass*
    o2scl_settings.get_convert_units().convert("kg","MeV",1.0);
}

int nucmass::get_nucleus(int Z, int N, nucleus &n) {
  n.Z=Z;
  n.N=N;
  n.A=Z+N;
  n.mex=mass_excess(Z,N)/o2scl_const::hc_mev_fm;
  n.m=n.mex+((Z+N)*m_amu-Z*m_elec)/o2scl_const::hc_mev_fm;
  n.ms=n.m;
  n.be=n.m-(N*m_neut+Z*m_prot)/o2scl_const::hc_mev_fm;
  if (n.A%2==0) n.g=1.0;
  else n.g=2.0;
  return 0;
}

double nucmass_table::mass_excess_d(double Z, double N) {
  int Z1=(int)Z;
  int N1=(int)N;
  double mz1n1=mass_excess(Z1,N1);
  double mz1n2=mass_excess(Z1,N1+1);
  double mz2n1=mass_excess(Z1+1,N1);
  double mz2n2=mass_excess(Z1+1,N1+1);
  double mz1=mz1n1+(N-N1)*(mz1n2-mz1n1);
  double mz2=mz2n1+(N-N1)*(mz2n2-mz2n1);
  return mz1+(Z-Z1)*(mz2-mz1);
}

nucmass_semi_empirical::nucmass_semi_empirical() {
  B=-16.0;
  Ss=18.0;
  Ec=0.7;
  Sv=23.7;
  Epair=13.0;
  nfit=5;
}

double nucmass_semi_empirical::mass_excess_d(double Z, double N) {
  double A=Z+N, cA=cbrt(A);
  double EoA=B+Ss/cA+Ec*Z*Z/cA/A+Sv*pow(1.0-2.0*Z/A,2.0);
  
  EoA+=-Epair*(cos(Z*o2scl_const::pi)+cos(N*o2scl_const::pi))/
    2.0/pow(A,1.5);
  
  double ret=EoA*A+Z*(m_prot+m_elec)+N*m_neut-A*m_amu;
  return ret;
}

int nucmass_semi_empirical::fit_fun(size_t nv, const ubvector &x) {
  B=-x[0]; Sv=x[1]; Ss=x[2]; Ec=x[3]; Epair=x[4];
  return 0;
}

int nucmass_semi_empirical::guess_fun(size_t nv, ubvector &x) {
  x[0]=-B; x[1]=Sv; x[2]=Ss; x[3]=Ec; x[4]=Epair;
  return 0;
}

nucmass_dvi::nucmass_dvi() {
  av=15.78;
  as=18.56;
  sv=31.51;
  y=2.75;
  ap=5.4;
  ac=0.71;
  nfit=10;
}

double nucmass_dvi::mass_excess_d(double Z, double N) {
  double A=Z+N, cA=cbrt(A);
  double T=fabs(N-Z)/2.0;

  double Delta=1.0;

  int iZ=((int)(Z+1.0e-10));
  int iN=((int)(N+1.0e-10));
  if (iZ%2==0 && iN%2==0) Delta=2.0;
  else if (iZ%2==1 && iN%2==1) Delta=0.0;
  
  double Lambda=(N-Z)/6.0/Z/(1.0+cA/y);//-5.0*pi2/6.0*d*d/r0/r0/cA/cA;
  double E=-av*A+as*cA*cA+sv*4.0*T*(T+1.0)/A/(1.0+y/cA)+
    ac*Z*(Z-1.0)/(1.0-Lambda)/cA-ap*Delta/cA;
  
  double ret=E+Z*(m_prot+m_elec)+N*m_neut-A*m_amu;
  ret-=shell_energy_interp(Z,N);

  return ret;
}

int nucmass_dvi::fit_fun(size_t nv, const ubvector &x) {
  av=x[0];
  as=x[1];
  sv=x[2];
  y=x[3];
  ap=x[4];
  ac=x[5];
  s_a1=x[6];
  s_a2=x[7];
  s_a3=x[8];
  s_anp=x[9];
  return 0;
}

int nucmass_dvi::guess_fun(size_t nv, ubvector &x) {
  x[0]=av;
  x[1]=as;
  x[2]=sv;
  x[3]=y;
  x[4]=ap;
  x[5]=ac;
  x[6]=s_a1;
  x[7]=s_a2;
  x[8]=s_a3;
  x[9]=s_anp;
  return 0;
}

nucmass_ibm_shell::nucmass_ibm_shell() {
  s_a1=-1.39;
  s_a2=0.02;
  s_a3=0.003;
  s_anp=0.075;
  shells[0]=2;
  shells[1]=8;
  shells[2]=14;
  shells[3]=28;
  shells[4]=50;
  shells[5]=82;
  shells[6]=126;
  shells[7]=184;
  shells[8]=228;
  shells[9]=308;
  shells[10]=406;
}

double nucmass_ibm_shell::shell_energy(int Z, int N) {

  int Dn=0, Dz=0, nv=0, zv=0;

  // Determine the appropriate proton shell
  if (Z<2) {
    Dz=2;
    zv=Z;
  } else {
    bool done=false;
    for(size_t i=0;i<nshells-1 && done==false;i++) {
      if (Z>=shells[i] && Z<shells[i+1]) {
	Dz=shells[i+1]-shells[i];
	zv=Z-shells[i];
	done=true;
      }
    }
    if (Z>=shells[nshells-1]) {
      nv=0;
      Dn=0;
      done=true;
    }
    if (done==false) {
      O2SCL_ERR2("Failed to do shell model correction ",
		     "in ldrop_shell::shell_energy().",exc_esanity);
    }
  }

  // Determine the appropriate neutron shell
  if (N<2) {
    Dn=2;
    nv=N;
  } else {
    bool done=false;
    for(size_t i=0;i<nshells-1 && done==false;i++) {
      if (N>=shells[i] && N<shells[i+1]) {
	Dn=shells[i+1]-shells[i];
	nv=N-shells[i];
	done=true;
      }
    }
    if (N>=shells[nshells-1]) {
      nv=0;
      Dn=0;
      done=true;
    }
    if (done==false) {
      O2SCL_ERR2("Failed to do shell model correction ",
		     "in ldrop_shell::shell_energy().",exc_esanity);
    }
  }
      
  double nvbar=Dn-nv;
  double zvbar=Dz-zv;
  double S2=nv*nvbar/Dn+zv*zvbar/Dz;
  double S3=nv*nvbar*(nv-nvbar)/Dn+zv*zvbar*(zv-zvbar)/Dz;
  double Snp=nv*nvbar*zv*zvbar/Dn/Dz;
  double ret=s_a1*S2+s_a2*S2*S2+s_a3*S3+s_anp*Snp;

  if (!std::isfinite(ret)) {
    cout << S2 << " " << S3 << " " << Snp << " " << ret << endl;
    O2SCL_ERR("Not finite in nucmass_ibm_shell.",exc_efailed);
  }

  return ret;
}

double nucmass_ibm_shell::shell_energy_interp(double Z, double N) {

  // Two-dimensional linear interpolation of the four
  // surrounding points

  int N0=(int)floor(N+1.0e-10);
  int Z0=(int)floor(Z+1.0e-10);
  double dN=N-N0;
  double dZ=Z-Z0;
  int Npn=N0+1;
  int Zpn=Z0;
  int N1=N0;
  int Z1=Z0+1;
  int N1pn=N0+1;
  int Z1pn=Z0+1;

  double shell0, shell1;
  {
    double shell_lo=shell_energy(Z0,N0);
    double shell_hi=shell_energy(Zpn,Npn);
    shell0=shell_lo+dN*(shell_hi-shell_lo);
  }
  {
    double shell_lo=shell_energy(Z1,N1);
    double shell_hi=shell_energy(Z1pn,N1pn);
    shell1=shell_lo+dN*(shell_hi-shell_lo);
  }

  return shell0+dZ*(shell1-shell0);
}

double nucmass_radius::density(double r, double Rfermi, double d, 
			       double rho0) {
  return rho0/(1+exp((r-Rfermi)/d));
}

double nucmass_radius::iand2_new(double r, double Rfermi, double d, 
				 double rho0) {
  return r*r*density(r,Rfermi,d,rho0);
}

double nucmass_radius::eval_N(double Rfermi, double d, double rho0) {
			      
  double N, N_err;
  eval_N_err(Rfermi,d,rho0,N,N_err);
  return N;
}

void nucmass_radius::eval_N_err(double Rfermi, double d, double rho0,
				double &N, double &N_err) {
  funct11 f=std::bind(std::mem_fn<double(double,double,double,double)>
		      (&nucmass_radius::iand2_new),
		      this,std::placeholders::_1,Rfermi,d,rho0);
  it.integ_err(f,0.0,0.0,N,N_err);
  N*=4.0*o2scl_const::pi;
  N_err*=4.0*o2scl_const::pi;
  return;
}

double nucmass_radius::iand(double r) {
  return urho0*4.0*o2scl_const::pi*pow(r,4.0)/(1+exp((r-uRfermi)/ud));
}

double nucmass_radius::iand2(double r) {
  return urho0*4.0*o2scl_const::pi*r*r/(1+exp((r-uRfermi)/ud));
}

double nucmass_radius::solve(double x) {
  uRfermi=x;
  funct11 it_fun2=std::bind(std::mem_fn<double(double)>
			    (&nucmass_radius::iand2),
			    this,std::placeholders::_1);
  return it.integ(it_fun2,0.0,0.0)-uN;
}

nucmass_radius::nucmass_radius() {
}

void nucmass_radius::eval_rms_rho(double rho0, double N, double d,
			     double &Rcd, double &Rfermi, double &Rrms) {
  
  urho0=rho0;
  ud=d;
  uN=N;

  Rcd=cbrt(3.0*N/4.0/o2scl_const::pi/rho0);
      
  uRfermi=Rcd;
  funct11 solve_fun=std::bind(std::mem_fn<double(double)>
			      (&nucmass_radius::solve),
			      this,std::placeholders::_1);
  cr.solve(uRfermi,solve_fun);
  Rfermi=uRfermi;

  funct11 it_fun=std::bind(std::mem_fn<double(double)>
			    (&nucmass_radius::iand),
			    this,std::placeholders::_1);
  Rrms=sqrt(it.integ(it_fun,0.0,0.0)/N);

  return;
}

void nucmass_radius::eval_rms_rsq(double Rfermi, double N, double d,
			     double &rho0, double &Rcd, double &Rrms) {
  
  ud=d;
  uN=N;
  
  uRfermi=Rfermi;
  urho0=1.0;
  funct11 it_fun2=std::bind(std::mem_fn<double(double)>
			    (&nucmass_radius::iand2),
			    this,std::placeholders::_1);
  rho0=N/it.integ(it_fun2,0.0,0.0);
  urho0=rho0;

  Rcd=cbrt(3.0*N/4.0/o2scl_const::pi/rho0);

  funct11 it_fun=std::bind(std::mem_fn<double(double)>
			   (&nucmass_radius::iand),
			   this,std::placeholders::_1);
  Rrms=sqrt(it.integ(it_fun,0.0,0.0)/N);

  return;
}

