/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <o2scl/nucmass_dz.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

bool nucmass_dz_table::is_included(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;
  int l_A=l_Z+l_N;

  // binary search for the correct A first
  if (data.get("A",mid)!=l_A) {
    if (data.get("A",mid)>l_A) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (data.get("A",mp)<l_A) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (data.get("A",mid)!=l_A) mid=hi;
    if (data.get("A",mid)!=l_A) {
      return false;
    }
  }

  // The cached point was the right one, so we're done
  if (data.get("Z",mid)==l_Z) {
    return true;
  }

  // Now look for the right N among all the Z's
  while (data.get("A",mid)==l_A) {

    if (data.get("Z",mid)==l_Z) {
      return true;
    } else if (data.get("Z",mid)>l_Z) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==((int)(n-1))) return false;
      mid++;
    }
  }
  
  return false;
}

double nucmass_dz_table::mass_excess(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;
  int l_A=l_Z+l_N;
  
  // binary search for the correct A first
  if (data.get("A",mid)!=l_A) {
    if (data.get("A",mid)>l_A) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (data.get("A",mp)<l_A) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (data.get("A",mid)!=l_A) mid=hi;
    if (data.get("A",mid)!=l_A) {
      O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		 " not found in nucmass_dz_table::mass_excess().").c_str(),
		exc_enotfound);
    }
  }
  
  // The cached point was the right one, so we're done
  if (data.get("Z",mid)==l_Z) {
    return data.get("ME",mid);
  }
  
  // Now look for the right N among all the Z's
  while (data.get("A",mid)==l_A) {
    
    if (data.get("Z",mid)==l_Z) {
      return data.get("ME",mid);
    } else if (data.get("Z",mid)>l_Z) {
      if (mid==0) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_dz_table::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid--;
    } else {
      if (mid==((int)(n-1))) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_dz_table::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_dz_table::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}

nucmass_dz_table::nucmass_dz_table(std::string model, bool external) {
  n=0;
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    if (model=="96") {
      fname=dir+"/nucmass/du_zu_96.o2";
    } else {
      fname=dir+"/nucmass/du_zu_95.o2";
    }
  }

  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  string name;
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input(hf,data,name);
#endif
  hf.close();
  
  n=data.get_nlines();
  last=n/2;
}

nucmass_dz_table::~nucmass_dz_table() {
}
  
nucmass_dz_fit::nucmass_dz_fit() {
  
  b.resize(10);
  
  b[0]=0.7043;
  b[1]=17.7418;
  b[2]=16.2562;
  b[3]=37.5562;
  b[4]=53.9017;
  b[5]=0.4711;
  b[6]=2.1307;
  b[7]=0.021;
  b[8]=40.5356;
  b[9]=6.0632;
  nfit=10;
  
  size_t dim[3]={9,2,2};
  onp.resize(3,dim);
  noc.resize(18,2);
  y.resize(2);
  pp.resize(2);
  oei.resize(2);
  dei.resize(2);
  qx.resize(2);
  dx.resize(2);
  op.resize(2);
  os.resize(2);
  n2.resize(2);
  dyda.resize(10);
  
  // The values 8.07138 and 7.28903 are those given in
  // the original Fortran.
  m_neut=8.07138+m_amu;
  m_prot=7.28903+m_amu-m_elec;
}

nucmass_dz_fit::~nucmass_dz_fit() {
}

bool nucmass_dz_fit::is_included(int Z, int N) {
  if (Z<=0 || N<=0 || N>240 || Z>240) return false;
  return true;
}

int nucmass_dz_fit::fit_fun(size_t nv, const ubvector &x) {
  for(size_t i=0;i<10;i++) b[i]=x[i];
  return 0;
}

int nucmass_dz_fit::guess_fun(size_t nv, ubvector &x) {
  for(size_t i=0;i<10;i++) x[i]=b[i];
  return 0;
}

double nucmass_dz_fit::binding_energy(int Z, int N) {

  // In the code below, comments with begin with "DZ:"
  // are comments from the original Fortran

  if (Z<=0 || N<=0) {
    O2SCL_ERR("Z or N <=0 in dz mass.",exc_efailed);
  }

  // Generally, arrays are shifted to be zero-indexed
  ubvector_int nn(2);
  nn[0]=N;
  nn[1]=Z;

  double a=N+Z;
  double t=fabs(((double)N-Z));
  double r=cbrt(a);
  double s=r*r;

  // Charge radius
  double rc=r*(1.0-0.25*t*t/a/a);
  double ra=rc*rc/r;
  double z2=Z*(Z-1.0);

  // Coulomb energy
  dyda[0]=(-z2+0.76*pow(z2,2.0/3.0))/rc;

  // Double check that they are all set to zero
  for(size_t ii=0;ii<y.size();ii++) y[ii]=0.0;
  for(size_t ii=0;ii<pp.size();ii++) pp[ii]=0.0;
  for(size_t ii=0;ii<oei.size();ii++) oei[ii]=0.0;
  for(size_t ii=0;ii<dei.size();ii++) dei[ii]=0.0;
  for(size_t ii=0;ii<qx.size();ii++) qx[ii]=0.0;
  for(size_t ii=0;ii<dx.size();ii++) dx[ii]=0.0;
  for(size_t ii=0;ii<op.size();ii++) op[ii]=0.0;
  for(size_t ii=0;ii<os.size();ii++) os[ii]=0.0;
  for(size_t ii=0;ii<n2.size();ii++) n2[ii]=0.0;
  for(size_t ii=0;ii<noc.size1();ii++) {
    for(size_t jj=0;jj<noc.size2();jj++) {
      noc(ii,jj)=0.0;
    }
  }
  onp.set_all(0.0);

  // 'E' is the final binding energy (here it's positive
  // for bound nuclei)
  double E;

  // 'ndef' is either 0 (spherical) or 1 (deformed), shifted by
  // one compared to DZ
  for(int ndef=0;ndef<2;ndef++) {

    // 'ju' is the same as in the DZ version
    int ju=0;
	
    y[ndef]=0.0;
    if (ndef==1) ju=4;
    for(int kk=1;kk<10;kk++) {
      dyda[kk]=0.0;
    }
    // DZ: beginning of loop over N and Z
    // 'j' and 'ell' are shifted, while 'k' is not. 
    for(int j=0;j<2;j++) {
      for(int ell=0;ell<18;ell++) {
	noc(ell,j)=0;
      }
      for(int ell=0;ell<2;ell++) {
	for(int k=0;k<=8;k++) {
	  onp.get(k,ell,j)=0.0;
	}
      }
      // DZ: For pairing calculation
      n2[j]=2*(nn[j]/2);

      // Neither 'ncum' or 'i' are shifted
      int ncum=0, i=0, id, i2;

      bool line20=true;
      while (line20) {
	line20=false;
	i++;
	i2=(i/2)*2;
	// DZ: For ssh j
	id=i+1;
	// DZ: For ssc r
	if (i2==i) id=i*(i-2)/4;
	ncum+=id;
	if (ncum<nn[j]) {
	  // DZ: nb of nucleons in each ssh
	  noc(i,j)=id;
	  line20=true;
	}
      }

      // DZ: imax = last subshell nb
      int imax=i+1;
      // DZ: HO number (p)
      if (i==0) {
	O2SCL_ERR("Bad arithmetic 1 in dz mass.",exc_efailed);
      }
      int ip=(i-1)/2;
      int ipm=i/2;
      pp[j]=ip;

      // Double check arithmetic
      if (nn[j]+id<ncum) {
	O2SCL_ERR("Bad arithmetic 2 in dz mass.",exc_efailed);
      }
      int moc=nn[j]+id-ncum;
	  
      noc(i,j)=moc-ju;
      noc(i+1,j)=ju;
      if (i2!=i) {
	oei[j]=moc+ip*(ip-1.0);
	dei[j]=ip*(ip+1.0)+2.0;
      } else {
	// DZ: ssh r
	// DZ: nb of nucleons in last EI shell
	oei[j]=moc-ju;
	// DZ: size of the EI shell
	dei[j]=(ip+1.0)*(ip+2.0)+2.0;
      }
      // DZ: n*(D-n)/D S3(j)
      qx[j]=oei[j]*(dei[j]-oei[j]-ju)/dei[j];
      // DZ: n*(D-n)*(2n-D)/D Q
      dx[j]=qx[j]*(2.0*oei[j]-dei[j]);
      // Scaling for deformed nuclei
      if (ndef==1) qx[j]/=sqrt(dei[j]);
	  
      // Variable 'i' is not shifted, and neither is ip
      // DZ: Amplitudes
      for(i=1;i<=imax;i++) {
	if (i<1) {
	  O2SCL_ERR("Bad arithmetic 3 in nucmass_dz.",exc_efailed);
	}
	ip=(i-1)/2;
	// DZ: For FM term
	double fact=sqrt((ip+1.0)*(ip+2.0));
	onp.get(ip,0,j)+=noc(i,j)/fact;
	// DZ: For spin-orbit term
	double vm=-1.0;
	if (2*(i/2)!=i) vm=0.5*ip;
	onp.get(ip,1,j)+=noc(i,j)*vm;
      }
      op[j]=0.0;
      os[j]=0.0;

      // Variable 'ip' is not shifted
      for(ip=0;ip<=ipm;ip++) {
	double pi=ip;
	double den=pow((pi+1.0)*(pi+2.0),1.5);
	op[j]+=onp.get(ip,0,j);
	os[j]+=onp.get(ip,1,j)*(1.0+onp.get(ip,0,j))*(pi*pi/den)+
	  onp.get(ip,1,j)*(1.0-onp.get(ip,0,j))*((4.0*pi-5.0)/den);
      }
      op[j]*=op[j];

      // End of loop from j=0 while j<2
    }
	
    // DZ: Master term (FM): volume
    dyda[1]=op[0]+op[1];
    // DZ: Surface
    dyda[2]=-dyda[1]/ra;
    // DZ: FM+SO
    dyda[1]+=os[0]+os[1];
    // DZ: Isospin term: volume
    dyda[3]=-t*(t+2.0)/(r*r);
    // DZ: Isospin term: surface
    dyda[4]=-dyda[3]/ra;
    if (ndef==0) {
      // Spherical case
      // DZ: S3 volume
      dyda[5]=dx[0]+dx[1];
      // DZ: S3 surface
      dyda[6]=-dyda[5]/ra;
      // DZ: QQ spherical
      double px=sqrt(pp[0])+sqrt(pp[1]);
      dyda[7]=qx[0]*qx[1]*pow(2,px);
    } else {
      // Deformed case
      // DZ: QQ deformed
      dyda[8]=qx[0]*qx[1];
    }
    // DZ: Wigner term
    dyda[4]+=t*(1.0-t)/(a*ra*ra*ra);
    // DZ: Pairing
    if ((n2[0]!=nn[0]) && (n2[1]!=nn[1])) dyda[9]=t/a;
    if (N>Z) {
      if ((n2[0]==nn[0]) && (n2[1]!=nn[1])) dyda[9]=1.0-t/a;
      if ((n2[0]!=nn[0]) && (n2[1]==nn[1])) dyda[9]=1.0;
    } else {
      if ((n2[0]==nn[0]) && (n2[1]!=nn[1])) dyda[9]=1.0;
      if ((n2[0]!=nn[0]) && (n2[1]==nn[1])) dyda[9]=1.0-t/a;
    }
    if ((n2[1]==nn[1]) && n2[0]==nn[0]) dyda[9]=2.0-t/a;
    // End of pairing
    for (int mss=1;mss<10;mss++) {
      dyda[mss]/=ra;
    }
    for (int mss=0;mss<10;mss++) {
      y[ndef]+=dyda[mss]*b[mss];
    }
    // DZ: Binding energy for deformed nuclei
    E=y[1];
    // DZ: Binding energy for spherical nuclei
    double de=y[1]-y[0];
    if (de<0.0 || Z<=50.0) E=y[0];
  }
  // Correct for O2scl sign convention
  return -E;
}

double nucmass_dz_fit::binding_energy_d(double Z, double N) {
  return binding_energy(((int)(Z+1.0e-8)),((int)(N+1.0e-8)));
}

double nucmass_dz_fit::mass_excess(int Z, int N) {
  return (binding_energy(Z,N)-((Z+N)*m_amu-Z*m_elec-N*m_neut-Z*m_prot));
}

double nucmass_dz_fit::mass_excess_d(double Z, double N) {
  return (binding_energy_d(Z,N)-((Z+N)*m_amu-Z*m_elec-N*m_neut-Z*m_prot));
}

nucmass_dz_fit_33::nucmass_dz_fit_33() {

  // The values 8.07132 and 7.28897 are those given in
  // the original Fortran.
  m_neut=8.07132+m_amu;
  m_prot=7.28897+m_amu-m_elec;

  a.resize(33);
  dyda.resize(33);
  fyda.resize(33);
  size_t dim[3]={2,3,2};
  op.resize(3,dim);
  fyd0.resize(33);
  onps.resize(2);
  n4.resize(2);
  size_t dim2[3]={9,2,2};
  onp.resize(3,dim2);
  ot.resize(3,dim2);
  oei.resize(2);
  dei.resize(2);
  nn.resize(2);
  noc.resize(18,2);
  op2.resize(2);
  jup.resize(2);
  jud.resize(2);
  ym.resize(2);
  op1.resize(2);
  n2.resize(2);
  shell.resize(2);
  sshell.resize(2);

  // F=full
  // P=partial
  // M=master
  // S=spin-orbit
  // C=cross
  // A=isoscalar
  // T=isubvector

  // FM+  
  a[0]=9.0914;
  // fm+ 
  a[1]=6.3355;
  // FS+ 
  a[2]=4.5791;
  // fs+ 
  a[3]=19.8946;
  // FS- 
  a[4]=1.7325;
  // fs- 
  a[5]=7.5247;
  // FC+ 
  a[6]=-7.1953;
  // fc+ 
  a[7]=-39.9787;
  // PM+ 
  a[8]=-0.3976;
  // pm+ 
  a[9]=0.8131;
  // PS+ 
  a[10]=-0.7435;
  // ps+ 
  a[11]=-3.7291;
  // PS- 
  a[12]=-0.1305;
  // ps- 
  a[13]=-0.6387;
  // S3  
  a[14]=0.4534;
  // s3  
  a[15]=2.0605;
  // SQ- 
  a[16]=0.3449;
  // sq- 
  a[17]=1.4727;
  // D3  
  a[18]=-1.0433;
  // d3
  a[19]=0.0000;
  // QQ+ 
  a[20]=5.2495;
  // qq+
  a[21]=0.0000;
  // D0 
  a[22]=-32.1007;
  // d0 
  a[23]=-151.1164;
  // QQ- 
  a[24]=-4.6103;
  // qq- 
  a[25]=-32.4238;
  // TT 
  a[26]=-37.3226;
  // tt 
  a[27]=-52.1673;
  // SS 
  a[28]=0.9597;
  // ss 
  a[29]=3.0024;
  // C
  a[30]=0.6977;
  // P0
  a[31]=6.0390;
  // P1
  a[32]=17.7960;

  nfit=33;
}

nucmass_dz_fit_33::~nucmass_dz_fit_33() {
}

int nucmass_dz_fit_33::fit_fun(size_t nv, const ubvector &x) {
  for(size_t i=0;i<33;i++) a[i]=x[i];
  return 0;
}
    
int nucmass_dz_fit_33::guess_fun(size_t nv, ubvector &x) {
  for(size_t i=0;i<33;i++) x[i]=a[i];
  return 0;
}

bool nucmass_dz_fit_33::is_included(int Z, int N) {
  if (Z<6 || N<6 || N>240 || Z>240) return false;
  return true;
}

double nucmass_dz_fit_33::binding_energy(int Z, int N) {

  // In the code below, comments with begin with "DZ:"
  // are comments from the original Fortran

  if (Z<=0 || N<=0) {
    O2SCL_ERR("Z or N <=0 in dz mass.",exc_efailed);
  }

  int idef=0;
  int imax=18;
  int mo=2;
  int maxp=8;
  nn[0]=N;
  nn[1]=Z;
  int nx=N;
  int nz=Z;
  double v=N+Z;
  double t=fabs(((double)N-Z));
  double r=cbrt(v);
  double s=r*r;
  double rc=r*(1.0-0.25*t*t/v/v);
  double ra=rc*rc/r;
  // ndef is shifted, 0 for spherical, 1 for deformed
  for(int ndef=0;ndef<2;ndef++) {
    ym[ndef]=0.0;
    int ju=0;
    jup[0]=0;
    jup[1]=0;
    jud[0]=0;
    jud[1]=0;
    if (ndef==1 && Z>50) ju=4;
    // kk is shifted
    for(int kk=0;kk<33;kk++) {
      fyda[kk]=0.0;
      fyd0[kk]=0.0;
      dyda[kk]=0.0;
    }
    // Initialize noc and onp to zero
    // j and m are shifted, k and i are not
    for(int j=0;j<2;j++) {
      for(int i=1;i<=imax;i++) {
	noc(i-1,j)=0;
      }
      for(int k=0;k<=maxp;k++) {
	for(int m=0;m<mo;m++) {
	  onp.get(k,m,j)=0.0;
	}
      }
    }
    // j is shifted again
    for(int j=0;j<2;j++) {
      int ncum=0;
      // None of these are shifted
      int i=0, i2, id, ip;
      // This while loop is for line 20 in the original Fortran
      bool done=false;
      while (done==false) {
	done=true;
	i++;
	i2=(i/2)*2;
	if (i2!=i) {
	  id=i+1;
	  if (ncum<nn[j]) sshell[j]=1;
	} else {
	  id=i*(i-2)/4;
	  if (ncum<nn[j]) sshell[j]=2;
	}
	ncum+=id;
	if (ncum<=nn[j]) {
	  if (i==0) {
	    std::cout << Z << " " << N << std::endl;
	    O2SCL_ERR("Indexing problem for 'i' in 'noc'.",exc_efailed);
	  }
	  noc(i-1,j)=id;
	  done=false;
	}
      }
      shell[j]=i;
      ip=(i-1)/2;
      int moc=nn[j]-ncum+id;
      if (ndef==1) {
	if (i2!=i) {
	  if (ju-moc>0) jud[j]=ju-moc;
	  else jud[j]=0;
	  jup[j]=0;
	} else {
	  // This section appears like it was written so that
	  // it didn't do anything in the original Fortran(?!),
	  // as jup[j] is just set equal to 'ju' below anyway.
	  if (ju<moc) jup[j]=ju;
	  else jup[j]=moc;
	  // 
	  jup[j]=ju;
	  jud[j]=0;
	}
      }
      noc(i-1,j)=moc-jup[j]+jud[j];
      noc(i,j)=jup[j];
      if (i<2) {
	std::cout << Z << " " << N << std::endl;
	O2SCL_ERR("Indexing problem for 'i' in 'noc'.",exc_efailed);
      }
      noc(i-2,j)-=jud[j];
      if (i2!=i) {
	oei[j]=moc+ip*(ip-1)-ju;
	dei[j]=ip*(ip+1)+2;
      } else {
	oei[j]=moc-ju;
	dei[j]=(ip+1)*(ip+2)+2;
      }
      int ipl=0;
      double vmr=0.0;
      double vmj=0.0;
      // Variable 'ii' is not shifted
      for(int ii=1;ii<=imax;ii++) {
	onps[j]=0.0;
	ip=(ii-1)/2;
	double degi=(ip+1)*(ip+2);
	if (degi<0.0) {
	  O2SCL_ERR("Arg of sqrt negative in dz 33.",exc_einval);
	}
	double fac=1.0/sqrt(degi);
	if (ip!=ipl) ipl++;
	double vm2=0.0;
	if ((2*ip+1)==ii) {
	  vm2=0.5*ip/(ip+1.0);
	  double degr=ip*(ip-1.0);
	  if (ip>2) {
	    vmr=0.5*(ip-1.0)/ip;
	    vmj=-1.0/ip;
	    if (noc(ii-1,j)<=degr) {
	      onps[j]=noc(ii-1,j)*vmr;
	    }
	    if (noc(ii-1,j)>degr) {
	      onps[j]=degr*vmr+(noc(ii-1,j)-degr)*vmj;
	    }
	  }
	}
	if ((2*ip+1)!=ii) vm2=-1.0/(ip+1.0);
	onp.get(ipl,1,j)+=noc(ii-1,j)*vm2;
	onp.get(ipl,0,j)+=noc(ii-1,j)*fac;
	fyd0[28]+=onps[j]*(onp.get(ipl,0,j)+onp.get(ipl,1,j));
      }
    }
    double facn, facz;
    if (ndef==1) {
      facn=1.0;
      facz=1.0;
    } else {
      facn=sqrt(dei[0]);
      facz=sqrt(dei[1]);
    }
    double dnnb=oei[0]*(dei[0]-oei[0])/dei[0];
    double dzzb=oei[1]*(dei[1]-oei[1])/dei[1];
    double qn=dnnb*facn/sqrt(dei[0]);
    double qz=dzzb*facz/sqrt(dei[1]);
    double di1n=dnnb*(2.0*oei[0]-dei[0])*facn*facn/dei[0];
    double di1z=dzzb*(2.0*oei[1]-dei[1])*facz*facz/dei[1];
    double s3=di1z+di1n;
    double qq0=(qn+qz)*(qn+qz);
    double qq1=(qn-qz)*(qn-qz);
    double qqp=qq0+qq1;
    double qqm=qq0-qq1;
    // Variable 'm' is shifted, and 'i' is not
    for(int m=0;m<mo;m++) {
      for(int i=0;i<=maxp;i++) {
	ot.get(i,m,0)=onp.get(i,m,0)+onp.get(i,m,1);
	ot.get(i,m,1)=onp.get(i,m,0)-onp.get(i,m,1);
      }
    }
    // Initialize all elements of 'op' to zero
    // All three of 'ell', 'm', and 'j' are shifted
    for(int ell=0;ell<2;ell++) {
      for(int m=0;m<3;m++) {
	for(int j=0;j<2;j++) {
	  op.get(ell,m,j)=0.0;
	}
      }
    }
    // Variable 'i' is not shifted, but n and m are
    for(int i=0;i<=maxp;i++) {
      double degi=(i+1.0)*(i+2.0);
      double fac=sqrt(degi);
      for(int n=0;n<2;n++) {
	for(int m=0;m<mo;m++) {
	  op.get(0,m,n)+=ot.get(i,m,n);
	  double otx=ot.get(i,m,n)*ot.get(i,m,n)*fac;
	  op.get(1,m,n)+=otx;
	}
      }
    }
    for(int n=0;n<2;n++) {
      op1[n]=0.0;
      op2[n]=0.0;
      for(int m=0;m<mo;m++) {
	double opxx=op.get(0,m,n)*op.get(0,m,n);
	op.get(0,m,n)=opxx;
      }
    }
    for(int n=0;n<2;n++) {
      for(int i=0;i<=maxp;i++) {
	double degi=(i+1.0)*(i+2.0);
	double fac=sqrt(degi);
	op1[n]+=ot.get(i,0,n)/fac;
	op2[n]+=ot.get(i,1,n)/fac;
      }
    }
    for(int n=0;n<2;n++) {
      op.get(0,2,n)=op1[n]*op2[n];
    }
    int k=-2;
    for(int ell=0;ell<2;ell++) {
      for(int m=0;m<3;m++) {
	for(int n=0;n<2;n++) {
	  k+=2;
	  fyd0[k]=op.get(ell,m,n);
	}
      }
    }
    for(int jf=0;jf<17;jf+=4) {
      fyd0[jf]+=fyd0[jf+2];
      fyd0[jf+2]=fyd0[jf]-2.0*fyd0[jf+2];
    }
    fyda[0]=fyd0[0];
    fyda[2]=fyd0[4];
    fyda[4]=fyd0[6];
    fyda[6]=fyd0[8];
    fyda[8]=fyd0[12];
    fyda[10]=fyd0[16];
    fyda[12]=fyd0[18];
    fyda[26]=t*(t+2)/s;
    fyda[28]=fyd0[28];
    if (ndef==0) {
      fyda[14]=s3;
      fyda[16]=qqm;
    } else {
      fyda[18]=s3;
      fyda[20]=qqp;
      fyda[22]=16.0-qqm;
      fyda[24]=qqm;
    }
    for(int mss=0;mss<=28;mss+=2) {
      dyda[mss]=fyda[mss]/ra;
      dyda[mss+1]=-dyda[mss]/ra;
      ym[ndef]+=dyda[mss]*a[mss]+dyda[mss+1]*a[mss+1];
    }
    double z2=Z*(Z-1);
    dyda[30]=(-z2+0.76*pow(z2,2.0/3.0))/rc;
    double rxz=1.0/ra;
    double vxz=1.0/v;
    double txz=rxz*t/v;
    double uxz=rxz-txz;
    // Variable 'ell' is shifted
    for(int ell=0;ell<2;ell++) {
      n2[ell]=2*(nn[ell]/2);
      dyda[31]+=(-rxz);
      if (sshell[ell]==2) dyda[32]+=vxz;
      if (n2[ell]==nn[ell]) dyda[31]+=uxz;
    }
    // Variable 'j' is shifted
    int j=1;
    if (nn[0]>=nn[1]) j=0;
    // Variable 'k' is shifted
    k=1-j;
    if (n2[j]==nn[j] && n2[k]!=nn[k]) dyda[31]+=(-txz);
    for(int mss=30;mss<=32;mss++) {
      ym[ndef]+=dyda[mss]*a[mss];
    }

    // End of ndef loop
  }

  double y=ym[0];
  if (Z>50 && ym[1]>ym[0]) y=ym[1];
  if (Z>50 && ym[1]>ym[0]) idef=1;
  double epair=dyda[31]*a[31]+dyda[32]*a[32];

  return -y;
}

double nucmass_dz_fit_33::binding_energy_d(double Z, double N) {
  return binding_energy(((int)(Z+1.0e-8)),((int)(N+1.0e-8)));
}

double nucmass_dz_fit_33::mass_excess(int Z, int N) {
  return (binding_energy(Z,N)-((Z+N)*m_amu-Z*m_elec-
			       N*m_neut-Z*m_prot));
}
    
double nucmass_dz_fit_33::mass_excess_d(double Z, double N) {
  return (binding_energy_d(Z,N)-((Z+N)*m_amu-Z*m_elec-N*m_neut-Z*m_prot));
}
