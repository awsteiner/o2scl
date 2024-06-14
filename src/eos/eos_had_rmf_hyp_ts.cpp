/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/eos_had_rmf_hyp.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(2);

  eos_had_rmf_hyp re;

  re.def_neutron.m=939.0/hc_mev_fm;
  re.def_proton.m=939.0/hc_mev_fm;
  re.mnuc=939.0/hc_mev_fm; 

  // These masses are just fiducial values, they don't impact
  // the calculation
  re.ms=500.0/hc_mev_fm; 
  re.mw=783.0/hc_mev_fm; 
  re.mr=763.0/hc_mev_fm; 

  // This is the K=240, Mstar/M=0.78 parameter set ("GM3") from GM91
  re.cs=sqrt(9.927);
  re.cw=sqrt(4.820);
  re.cr=sqrt(4.791);
  re.b=0.008659;
  re.c=-0.00241;
  re.inc_cascade=true;

  cout << "GM3: " << endl;
  re.saturation();
  cout << "  Saturation density: " << re.n0 << endl;
  t.test_rel(re.n0,0.153,1.0e-3,"sat density");
  cout << "  Effective mass: " << re.msom << " "
       << re.def_neutron.ms/re.def_neutron.m << endl;
  t.test_rel(re.msom,re.def_neutron.ms/re.def_neutron.m,1.0e-6,"msom");
  t.test_rel(re.msom,0.78,1.0e-2,"msom 2");
  cout << "  Compressibility: " << re.comp*hc_mev_fm << endl;
  t.test_rel(re.comp*hc_mev_fm,240.0,2.0,"comp");
  cout << "  Symmetry energy: " << re.esym*hc_mev_fm << endl;
  t.test_rel(re.esym*hc_mev_fm,32.5,0.1,"esym");
  cout << "  Zero pressure: " << re.def_thermo.pr << endl;
  t.test_rel(re.def_thermo.pr,0.0,1.0e-8,"zero press");
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (re.def_thermo.ed/re.n0-re.mnuc)*hc_mev_fm << endl;
  t.test_rel(re.eoa,re.def_thermo.ed/re.n0-re.mnuc,1.0e-6,"eoa");
  t.test_rel(re.eoa,-16.3/hc_mev_fm,0.2/hc_mev_fm,"eoa 2");
  cout << "  Thermomodynamic identity: " 
       << re.def_thermo.ed+re.def_thermo.pr-
    re.def_neutron.n*re.def_neutron.mu-
    re.def_proton.n*re.def_proton.mu << endl;
  t.test_rel(re.def_thermo.ed+re.def_thermo.pr-
	     re.def_neutron.n*re.def_neutron.mu-
	     re.def_proton.n*re.def_proton.mu,0.0,2.0e-9,"TI");
  cout << endl;

  // Compare with Table I in GM91.

  // This corrects an apparent typo in GM91, 568 -> 658
  vector<double> xw_dat={0.091,0.233,0.375,0.517,0.658,0.800,
			 0.942,1.08,1.23};
  
  fermion e, mu;
  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_const::mass_electron_f<double>()),2.0);
  mu.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_const::mass_muon_f<double>()),2.0);
  std::shared_ptr<table_units<> > eos_table(new table_units<>);
  nstar_cold nc;

  uniform_grid_end_width<double> nBg(0.04,1.2,0.01);
  ubvector nB_grid;
  nBg.vector(nB_grid);
  ubvector guess(5);
  fermion_rel frel;

  tov_solve ts;
  eos_tov_interp eti;
  
  table_units<> xM;
  xM.line_of_names("xs xw Mmax");
  
  size_t j=0;
  for(re.xs=0.2;re.xs<1.0001;re.xs+=0.1,j++) {
    re.xr=re.xs;
    
    re.calc_xw(-28.0/hc_mev_fm);
    
    cout << j << " " << re.xs << " " << re.xw << " "
	 << fabs(re.xw-xw_dat[j])/xw_dat[j] << " ";
    
    if (j<7) {
      t.test_rel(re.xw,xw_dat[j],0.002,"xw");
    } else {
      t.test_rel(re.xw,xw_dat[j],0.005,"xw");
    }

    //cout.precision(10);
    eos_leptons elep;
    elep.include_muons=true;
    cout << "H1." << endl;
    re.beta_eq_T0(nB_grid,guess,elep,eos_table);
    cout << "H2." << endl;
    /*
    for(size_t j=0;j<eos_table->get_nlines();j+=10) {
      cout << j << " ";
      cout << eos_table->get("nb",j) << " ";
      cout << eos_table->get("ed",j) << " ";
      cout << eos_table->get("pr",j) << endl;
    }
    */

    eti.read_table(*eos_table,"ed","pr","nb");
    ts.set_eos(eti);
    ts.max();
    shared_ptr<table_units<> > tov_max=ts.get_results();
    double m_max=tov_max->max("gm");
    double nb_max=tov_max->get("nb",0);
					    
    cout << m_max << " " << nb_max << endl;

    xM.line_of_data(vector<double>({re.xs,re.xw,m_max}));

    if (false && j==4) {
      cout << endl;
      for(size_t i=0;i<eos_table->get_nlines();i++) {
	cout << eos_table->get("nb",i) << " ";
	cout << eos_table->get("ne",i) << " ";
	cout << eos_table->get("nmu",i) << " ";
	cout << eos_table->get("nn",i) << " ";
	cout << eos_table->get("np",i) << " ";
	cout << eos_table->get("nlam",i) << " ";
	cout << eos_table->get("nsigm",i) << " ";
	cout << eos_table->get("nsigz",i) << " ";
	cout << eos_table->get("nsigp",i) << endl;
      }
      cout << endl;
    }      
    
    
  }
  cout << endl;

  // GM1, to compare with Fig. 3 in GM91

  re.cs=sqrt(11.79);
  re.cw=sqrt(7.149);
  re.cr=sqrt(4.411);
  re.b=0.002947;
  re.c=-0.001070;
  re.xs=0.6;
  re.xr=0.6;
  re.calc_xw(-28.0/hc_mev_fm);

  re.saturation();
  cout << "  Saturation density: " << re.n0 << endl;
  t.test_rel(re.n0,0.153,2.0e-3,"sat density");
  cout << "  Effective mass: " << re.msom << " "
       << re.def_neutron.ms/re.def_neutron.m << endl;
  t.test_rel(re.msom,re.def_neutron.ms/re.def_neutron.m,1.0e-6,"msom");
  t.test_rel(re.msom,0.7,1.0e-2,"msom 2");
  cout << "  Compressibility: " << re.comp*hc_mev_fm << endl;
  t.test_rel(re.comp*hc_mev_fm,300.0,2.0,"comp");
  cout << "  Symmetry energy: " << re.esym*hc_mev_fm << endl;
  t.test_rel(re.esym*hc_mev_fm,32.5,0.1,"esym");
  cout << "  Zero pressure: " << re.def_thermo.pr << endl;
  t.test_rel(re.def_thermo.pr,0.0,1.0e-8,"zero press");
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (re.def_thermo.ed/re.n0-re.mnuc)*hc_mev_fm << endl;
  t.test_rel(re.eoa,re.def_thermo.ed/re.n0-re.mnuc,1.0e-6,"eoa");
  t.test_rel(re.eoa,-16.3/hc_mev_fm,0.8/hc_mev_fm,"eoa 2");
  cout << endl;

  eos_leptons elep;
  elep.include_muons=true;
  re.beta_eq_T0(nB_grid,guess,elep,eos_table);

  if (false) {
    re.verbose=2;
    re.calc_hyp_e_nobeta(0.48,0.4,0.1,
                         re.def_neutron,re.def_proton,re.def_lambda,
                         re.def_sigma_p,re.def_sigma_z,re.def_sigma_m,
                         re.def_cascade_z,re.def_cascade_m,
                         re.def_thermo);
    exit(-1);
  }
  
  eti.read_table(*eos_table,"ed","pr","nb");
  ts.set_eos(eti);
  ts.max();
  shared_ptr<table_units<> > tov_max=ts.get_results();
  double m_max=tov_max->max("gm");
  double nb_max=tov_max->get("nb",0);
    
  cout.precision(4);
  for(size_t i=0;i<eos_table->get_nlines();i++) {
    cout << eos_table->get("nb",i)/re.n0 << " ";
    cout << eos_table->get("nn",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("np",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("ne",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("nmu",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("nlam",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("nsigm",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("nsigz",i)/eos_table->get("nb",i) << " ";
    cout << eos_table->get("nsigp",i)/eos_table->get("nb",i) << endl;
  }
  cout << endl;
  
  cout << m_max << " " << nb_max << endl;

  tov_max->add_col_from_table(*eos_table,"nb","ne");
  tov_max->add_col_from_table(*eos_table,"nb","nmu");
  tov_max->add_col_from_table(*eos_table,"nb","nn");
  tov_max->add_col_from_table(*eos_table,"nb","np");
  tov_max->add_col_from_table(*eos_table,"nb","nlam");
  tov_max->add_col_from_table(*eos_table,"nb","nsigp");
  tov_max->add_col_from_table(*eos_table,"nb","nsigz");
  tov_max->add_col_from_table(*eos_table,"nb","nsigm");
  
  /*
    o2graph -read eos_had_rmf_hyp_ts.o2 fig3 \
    -set xlo 0.0 -set xhi 12.0 -set ylo 1.0e-3 \
    -set yhi 1 -set logy 1 \
    -function nlam/nb xlam -plot r xlam \
    -function nn/nb xn -plot r xn \
    -function np/nb xp -plot r xp \
    -function ne/nb xe -plot r xe \
    -function nmu/nb xmu -plot r xmu \
    -function nsigp/nb xsigp -plot r xsigp \
    -function nsigz/nb xsigz -plot r xsigz \
    -function nsigm/nb xsigm -plot r xsigm \
    -show
  */

  hdf_file hf;
  hf.open_or_create("eos_had_rmf_hyp_ts.o2");
  hdf_output(hf,*eos_table,"fig3_eos");
  hdf_output(hf,*tov_max,"fig3");
  hdf_output(hf,xM,"xM");
  hf.close();
  
  t.report();

  return 0;
}


