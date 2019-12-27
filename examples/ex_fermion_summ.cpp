/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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

/* Example: ex_fermion_summ.cpp
   -------------------------------------------------------------------
*/
#include <iostream>
#include <string>
#include <vector>

#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/contour.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

fermion_zerot fzt;

void fermion_func(double x, double y, double T,
		  fermion_deriv &f, fermion_deriv_rel &fr,
		  double &log10_Pt, double &phi, double &psi,
		  double &s, double &cv, double &cp, double &cs2,
		  double &kfom, int &lm) {

  f.m=T/pow(10.0,y);
  double nc=f.g*f.m*f.m*f.m/2.0/pi2;
  f.n=pow(10.0,x)*nc;
  fr.calc_density(f,T);
  
  log10_Pt=log10(f.pr/nc/f.m);
  phi=f.mu/T;
  psi=(f.mu-f.m)/T;
  s=f.en/f.n;
  cv=fr.heat_cap_ppart_const_vol(f,T);
  cp=fr.heat_cap_ppart_const_press(f,T);
  cs2=fr.squared_sound_speed(f,T);
  fzt.kf_from_density(f);
  kfom=f.kf/f.m;
  lm=fr.last_method;

  return;
}

void pairs_func(double x, double y, double T,
		fermion_deriv &f, fermion_deriv_rel &fr, 
		double &log10_Pt, double &phi, double &psi,
		double &s, double &cv, double &cp, double &cs2,
		double &kfom, int &lm) {
  
  f.m=T/pow(10.0,y);
  double nc=f.g*f.m*f.m*f.m/2.0/pi2;
  f.n=pow(10.0,x)*nc;
  fr.pair_density(f,T);
  
  log10_Pt=log10(f.pr/nc/f.m);
  phi=f.mu/T;
  psi=(f.mu-f.m)/T;
  s=f.en/f.n;
  cv=fr.heat_cap_ppart_const_vol(f,T);
  cp=fr.heat_cap_ppart_const_press(f,T);
  cs2=fr.squared_sound_speed(f,T);
  fzt.kf_from_density(f);
  kfom=f.kf/f.m;
  lm=fr.last_method;
  
  return;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  
  fermion_deriv f(1.0,2.0);
  f.inc_rest_mass=true;
  fermion_deriv_rel fr;
  // Choose a spin degeneracy of 2 and T=1
  double T=1.0;

  hdf_file hf;
  hf.open_or_create("ex_fermion_summ.o2");

  double refine=1.0;

  {
  
    // ------------------------------------------------------------
    // Fermions

    uniform_grid_end_width<double> lp_grid(-8.0,7.0,1.0);
    vector<double> phi_grid={-10.0,-3.0,0.0,10.0,100.0,1000.0,1.0e4};
    vector<double> psi_grid={-30.0,-20.0,-10.0,-3.0,0.0,
			     10.0,100.0,1000.0,1.0e4};
    vector<double> s_grid={1.0e-7,1.0e-5,0.001,0.01,0.1,1.0,5.0,15.0,
			   25.0,35.0};
    vector<double> cv_grid={3.0e-5,3.0e-4,3.0e-3,0.03,0.3,0.9,1.5,
			    1.8,2.1,2.7,2.9,2.99,2.999};
    vector<double> cp_grid={4.0e-5,4.0e-4,4.0e-3,0.04,0.4,1.2,2.0,
			    2.4,2.8,3.6,3.9,3.99,3.999};
    vector<double> cs2_grid={0.001,0.01,0.03,0.1,0.3,0.33,0.3333,
			     0.333333};
    vector<double> kfom_grid={0.5,1.0,2.0};
    
    cout << "Fermions: " << endl;
    
    uniform_grid_end_width<double> x_grid(-5.0,5.0,0.1/refine);
    uniform_grid_end_width<double> y_grid(-5.0,3.0,0.1/refine);
    table3d t3d;
    t3d.set_xy("x",x_grid,"y",y_grid);
    t3d.line_of_names("log10_Pt phi psi s cv cp cs2 kfom lm dndT dndmu dsdT");
    
    for(size_t i=0;i<t3d.get_nx();i++) {
      for(size_t j=0;j<t3d.get_ny();j++) {

	double log10_Pt, phi, psi, s, cv, cp, cs2, kfom;
	int lm;
	
	fermion_func(x_grid[i],y_grid[j],T,f,fr,
		     log10_Pt,phi,psi,s,cv,cp,cs2,kfom,lm);

	if (i==0) {
	  cout << i << " " << j << " " 
	       << f.dndT << " " << f.dndmu << " " << f.dsdT << " "
	       << cs2 << " ";
	  cout << f.n*f.n*f.dsdT << " "
	       << -2.0*f.n*f.en*f.dndT << " "
	       << f.en*f.en*f.dndmu << " "
	       << f.dndmu*f.dsdT << " " << f.dndT*f.dndT << " "
	       << lm << endl;
	  //if (j==t3d.get_ny()-1) exit(-1);
	}

	t3d.set(i,j,"log10_Pt",log10_Pt);
	t3d.set(i,j,"phi",phi);
	t3d.set(i,j,"psi",psi);
	t3d.set(i,j,"s",s);
	t3d.set(i,j,"cv",cv);
	t3d.set(i,j,"cp",cp);
	t3d.set(i,j,"cs2",cs2);
	t3d.set(i,j,"kfom",kfom);
	t3d.set(i,j,"lm",lm);
	t3d.set(i,j,"dndT",f.dndT);
	t3d.set(i,j,"dsdT",f.dsdT);
	t3d.set(i,j,"dndmu",f.dndmu);
      }
    }

    hdf_output(hf,(const table3d &)t3d,"fermions");
    
    for(size_t k=0;k<8;k++) {
      contour co;
      co.set_data(x_grid,y_grid,t3d.get_slice(t3d.get_slice_name(k)));
      vector<contour_line> conts;

      if (k==0) {
	vector<double> dtmp;
	lp_grid.vector(dtmp);
	co.set_levels(dtmp.size(),dtmp);
      } else if (k==1) {
	co.set_levels(phi_grid.size(),phi_grid);
      } else if (k==2) {
	co.set_levels(psi_grid.size(),psi_grid);
      } else if (k==3) {
	co.set_levels(s_grid.size(),s_grid);
      } else if (k==4) {
	co.set_levels(cv_grid.size(),cv_grid);
      } else if (k==5) {
	co.set_levels(cp_grid.size(),cp_grid);
      } else if (k==6) {
	co.set_levels(cs2_grid.size(),cs2_grid);
      } else if (k==7) {
	co.set_levels(kfom_grid.size(),kfom_grid);
      }
      
      co.calc_contours(conts);
      
      hdf_output(hf,conts,((string)"fermions_")+
		 t3d.get_slice_name(k));
      
    }

  }
  {
  
    // ------------------------------------------------------------
    // Fermions and antifermions

    uniform_grid_end_width<double> lp_grid(-8.0,8.0,1.0);
    vector<double> phi_grid={3.0,10.0,100.0,1000.0,1.0e4};
    vector<double> psi_grid={-5.0,-3.0,-1.0,0.0,3.0,10.0,100.0,1000.0,1.0e4};
    vector<double> s_grid={1.0e-7,1.0e-5,0.001,0.01,0.1,1.0,5.0,15.0,
			   25.0,35.0};
    vector<double> cv_grid={3.0e-5,3.0e-4,3.0e-3,0.03,0.3,0.9,1.5,
			    1.8,2.1,2.7,2.9,2.99,2.999};
    vector<double> cp_grid={4.0e-5,4.0e-4,4.0e-3,0.04,0.4,1.2,2.0,
			    2.4,2.8,3.6,3.9,3.99,3.999};
    vector<double> cs2_grid={0.001,0.01,0.03,0.1,0.3,0.45,0.49,
			     0.499,0.4999,0.49999};
    
    cout << "Fermions and anti-fermions: " << endl;
      
    uniform_grid_end_width<double> x_grid(-5.0,5.0,0.1/refine);
    uniform_grid_end_width<double> y_grid(-5.0,2.0,0.1/refine);
    table3d t3d;
    t3d.set_xy("x",x_grid,"y",y_grid);
    t3d.line_of_names("log10_Pt phi psi s cv cp cs2 kfom lm");
      
    for(size_t i=0;i<t3d.get_nx();i++) {
      for(size_t j=0;j<t3d.get_ny();j++) {

	double log10_Pt, phi, psi, s, cv, cp, cs2, kfom;
	int lm;

	pairs_func(x_grid[i],y_grid[j],T,f,fr,
		   log10_Pt,phi,psi,s,cv,cp,cs2,kfom,lm);
	t3d.set(i,j,"log10_Pt",log10_Pt);
	t3d.set(i,j,"phi",phi);
	t3d.set(i,j,"psi",psi);
	t3d.set(i,j,"s",s);
	t3d.set(i,j,"cv",cv);
	t3d.set(i,j,"cp",cp);
	t3d.set(i,j,"cs2",cs2);
	t3d.set(i,j,"kfom",kfom);
	t3d.set(i,j,"lm",lm);
      }
    }
      
    hdf_output(hf,(const table3d &)t3d,"pairs");

    for(size_t k=0;k<8;k++) {
      contour co;
      co.set_data(x_grid,y_grid,t3d.get_slice(t3d.get_slice_name(k)));
      vector<contour_line> conts;
      if (k==0) {
	vector<double> dtmp;
	lp_grid.vector(dtmp);
	co.set_levels(dtmp.size(),dtmp);
      } else if (k==1) {
	co.set_levels(phi_grid.size(),phi_grid);
      } else if (k==2) {
	co.set_levels(psi_grid.size(),psi_grid);
      } else if (k==3) {
	co.set_levels(s_grid.size(),s_grid);
      } else if (k==4) {
	co.set_levels(cv_grid.size(),cv_grid);
      } else if (k==5) {
	co.set_levels(cp_grid.size(),cp_grid);
      } else {
	co.set_levels(cs2_grid.size(),cs2_grid);
      }
      co.calc_contours(conts);
      
      hdf_output(hf,conts,((string)"pairs_")+
		 t3d.get_slice_name(k));
	
    }

  }
  
  hf.close();

  t.report();

  return 0;
}
