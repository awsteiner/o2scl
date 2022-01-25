/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/** \file ex_eos_had_skyrme.cpp
    \brief File defining \ref ex_eos_had_skyrme class
*/
/* Example: ex_eos_had_skyrme.cpp
   -------------------------------------------------------------------
*/
#include <o2scl/test_mgr.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/format_float.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/cli.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

/** \brief Output data for a Skyrme EOS [Example class]
 */
class ex_skyrme_data {
  
public:

  /// Name of model
  string name;
  /// Saturation density
  double n0;
  /// Binding energy
  double B;
  /// Compressibility
  double K;
  /// Symmetry energy
  double S;
  /// Symmetry energy slope parameter
  double L;
  /// Maximum mass
  double m_max;
  /// Radius of maximum mass star
  double r_max;
  /// Central baryon density of maximum mass star
  double nb_max;
  /// Radius of a 1.4 solar mass star
  double r_14;
  /// Central baryon density of a 1.4 solar mass neutron star
  double nb_14;
  /// True if the EOS is acausal
  double acausal_nb;
  /// True if the pressure is flat
  double pressure_dec_nb;
  /// Quality of neutron matter
  double neut_qual;
  /// True if the maximum mass is large enough
  bool max_mass_ok;
  /// True if the pressure is nondecreasing
  bool inc_pressure;
  /// True if neutron matter is always positive
  bool pos_neut;
  /// True if saturation is good
  bool good_sat;
  /** \brief The baryon density at which the beta-equilibrium EOS becomes 
      pure neutron matter
  */
  double pure_neut;
  /// Desc
  int other;
  /// If true, the EOS calculation was successful
  bool success;
  /// Alternate description of symmetry energy
  double alt_S;
  /// Energy of neutron matter at saturation
  double E_neut_n0;
  /// Pressure of neutron matter at saturation
  double P_neut_n0;
};

/** \brief Class to analyze Skyrme EOSs and output the results
    [Example class]
*/
class ex_eos_had_skyrme {

protected:

  /// To compute cold neutron stars
  nstar_cold nst;

  /// To compute hot neutron stars
  nstar_hot nh;
  
  /// Neutron
  fermion n;

  /// Proton
  fermion p;

  /// Neutron
  fermion_deriv nd;

  /// Proton
  fermion_deriv pd;

  /// Thermodynamics
  thermo th;

  /// File for I/O
  hdf_file hf;

  /// Results
  ex_skyrme_data res;

  /// Formatting output
  format_float fd;

  /// Neutron matter
  table_units<> tneut;

  /// Neutron matter
  table_units<> tneut2;

public:

  /// Base EOS model
  eos_had_skyrme sk;

  /// Model name
  string name;

  /// Parameter t0 in MeV
  double t0hc;

  /// Parameter t1 in MeV
  double t1hc;

  /// Parameter t2 in MeV
  double t2hc;

  /// Parameter t3 in MeV
  double t3hc;

  /// Parameter W0 in MeV
  double W0hc;

  /// Parameter b4
  double b4hc;

  /// Parameter b4p
  double b4phc;

  /// Verbose parameter
  int verbose;

  /// If true, create output files for individual EOSs
  bool output_files;

  /// Prefix for output files
  string file_prefix;

  /// If true, output JSON
  bool json_mode;

  ex_eos_had_skyrme() {

    n.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
    p.init(o2scl_settings.get_convert_units().convert
	   ("kg","1/fm",o2scl_mks::mass_proton),2.0);
    n.non_interacting=false;
    p.non_interacting=false;
    
    nst.set_eos(sk);
    nst.def_eos_tov.verbose=0;
    nst.def_tov.verbose=0;
    nst.include_muons=true;
    verbose=1;

    nh.set_eos_T(sk);
    nh.def_eos_tov.verbose=0;
    nh.def_tov.verbose=0;
    nh.include_muons=true;
    
    fd.html_mode();

    output_files=false;
    file_prefix="skyrme_";

    json_mode=false;
  }

  /** \brief Generate a table comparing neutron matter and neutron 
      star matter
   */
  void compare_neut_nstar() {
    
    tneut.clear_table();
    tneut2.clear_table();

    tneut.new_column("nb");
    tneut.new_column("pr_neut");
    tneut.new_column("pr_nstar");

    tneut2.new_column("ed");
    tneut2.new_column("pr_neut");
    tneut2.new_column("pr_nstar");

    tneut.set_unit("nb","fm^-3");
    tneut.set_unit("pr_neut","fm^-4");
    tneut.set_unit("pr_nstar","fm^-4");

    tneut2.set_unit("ed","fm^-4");
    tneut2.set_unit("pr_neut","fm^-4");
    tneut2.set_unit("pr_nstar","fm^-4");

    for(double nbt=0.06;nbt<0.9001;nbt+=0.02) {
      n.n=nbt;
      p.n=0;
      int ret=sk.calc_e(n,p,th);
      vector<double> line={nbt,th.pr,0.0};
      tneut.line_of_data(line.size(),line);
      line={th.ed,th.pr,0.0};
      tneut2.line_of_data(line.size(),line2);
    }
    
    return;
  }

  /// Check if the pressure of neutron matter is positive
  int check_pressure() {
    bool failed_once=false;

    res.pos_neut=true;
    res.inc_pressure=true;

    double prlast=0.0;
    for(double nbt=0.005;nbt<0.1601 && failed_once==false;nbt+=0.01) {
      n.n=nbt;
      p.n=0;
      int ret=sk.calc_e(n,p,th);
      
      if (th.pr<prlast) {
	res.inc_pressure=false;
	res.success=false;
	prlast=th.pr;
      }
      if (th.ed-n.n*n.m<0.0) {
	res.pos_neut=false;
	res.success=false;
      }
    }
    if (verbose>0) {
      if (failed_once==false) {
	cout << "Pressure increases and energy of neutrons is positive.\n" 
	     << endl;
      } else {
	cout << "Pressure decreased or energy of neutrons is negative.\n" 
	     << endl;
      }
    }
    return 0;
  }
  
  /** \brief Check low-density neutron matter

      I'm not sure where the data originally came from, but it was
      based on something similar to the left panel of Fig. 1 in
      Gandolfi et al. (2015). 
      
      \todo AWS 12/21/2020: I need to document here what "quality"
      values are reasonable. I should probably make this optional as
      well, as many Skyrme interactions are not intended to be used at
      this density.
   */
  int low_neutron_mat() {
    
    double g_kf[4]={0.06,0.14,0.26,0.54};
    double g_rat[4]={0.80,0.67,0.57,0.52};

    table_units<> ln;
    ln.line_of_names("nb epb kf fermi rat");
    if (verbose>0) cout << "Low density neutron matter: " << endl;
    for(double nbt=0.000003;nbt<0.0051;nbt*=1.2) {
      n.n=nbt;
      p.n=0;
      int ret=sk.calc_e(n,p,th);

      // The Skyrme energy per baryon
      double epb=(th.ed-n.n*n.m)/nbt*hc_mev_fm;

      // The Fermi gas energy per baryon
      double fermi=pow(n.kf,5.0)/10.0/o2scl_const::pi2/n.ms/nbt*hc_mev_fm;

      if (verbose>1) {
        cout << nbt << " " << epb << " "
             << n.kf << " " << fermi << " " << epb/fermi << endl;
      }
      vector<double> line={nbt,epb,n.kf,fermi,epb/fermi};
      ln.line_of_data(line.size(),line);
    }
    if (verbose>1) cout << endl;

    // At the saturation density, evaluate the properties of 
    // neutron matter
    n.n=sk.n0;
    p.n=0.0;
    sk.calc_e(n,p,th);
    res.E_neut_n0=(th.ed/sk.n0-n.m)*hc_mev_fm;
    res.P_neut_n0=th.pr*hc_mev_fm;

    res.neut_qual=0.0;
    for(size_t i=0;i<4;i++) {
      if (verbose>0) {
        cout << g_kf[i] << " " << g_rat[i] << " "
	     << ln.interp("kf",g_kf[i],"rat") << " "
	     << ln.interp("kf",g_kf[i],"epb") << endl;
      }
      res.neut_qual+=g_kf[i]*fabs(g_rat[i]-ln.interp("kf",g_kf[i],"rat"));
    }
    ln.add_constant("neut_qual",res.neut_qual);
    if (verbose>0) {
      cout << "Quality: " << res.neut_qual << endl;
      cout << endl;
    }

    return 0;
  }

  /// Test saturation density
  int saturation_prop() {

    n.n=0.08;
    p.n=0.08;
    sk.saturation();

    // Collect saturation results
    res.n0=sk.n0;
    res.B=sk.eoa*hc_mev_fm;
    res.K=sk.comp*hc_mev_fm;
    res.S=sk.esym*hc_mev_fm;
    res.alt_S=sk.fesym_diff(sk.n0)*hc_mev_fm;
    res.L=sk.fesym_slope(sk.n0)*hc_mev_fm;

    if (json_mode) {
      cout << "{\"n0\":" << sk.n0 << "," << endl;
      cout << "\"EoA\":" << sk.eoa*hc_mev_fm << "," << endl;
      cout << "\"K\":" << sk.comp*hc_mev_fm << "," << endl;
      cout << "\"msom\":" << sk.msom << "," << endl;
      cout << "\"S\":" << sk.esym*hc_mev_fm << "," << endl;
      cout << "\"S2\":" << res.alt_S << "," << endl;
      cout << "\"L\":" << res.L << "}" << endl;
    } else if (verbose>0) {
      cout << "Saturation: " << endl;
      cout << "n_0=" << sk.n0 << " fm^-3" << endl;
      cout << "E_B=" << sk.eoa*hc_mev_fm << " MeV" << endl;
      cout << "K=" << sk.comp*hc_mev_fm << " MeV" << endl;
      cout << "M^*/M=" << sk.msom << endl;
      cout << "S=" << sk.esym*hc_mev_fm << " MeV" << endl;
      cout << "S2=" << res.alt_S << " MeV" << endl;
      cout << "L=" << res.L << " MeV" << endl;
      cout << endl;
    }

    // Check saturation properties

    res.good_sat=true;

    if (fabs(sk.n0-0.16)>0.013) {
      res.good_sat=false;
    }
    if (fabs(sk.eoa*hc_mev_fm+16.0)>1.2) {
      res.good_sat=false;
    }
    if (fabs(sk.comp*hc_mev_fm-220.0)>20.0) {
      res.good_sat=false;
    }
    if (fabs(sk.esym*hc_mev_fm-32.0)>4.0) {
      res.good_sat=false;
    }
    if (res.good_sat==false) {
      cout << "Bad saturation." << endl;
      res.success=false;
    }

    return 0;
  }
  
  /// Compute the M vs. R curve
  int mvsr() {
    
    res.other=0;
    
    cout << "EOS:" << endl;
    {
      double nb_last=10.0;
      for(double nb=0.16;nb<=2.0001;nb+=0.001) {
	p.n=0.0;
	n.n=nb-p.n;
	sk.calc_e(n,p,th);
	double me=o2scl_settings.get_convert_units().convert
	  ("kg","1/fm",o2scl_mks::mass_electron);
	if (n.mu-p.mu-me<0.0) {
	  nb_last=nb-0.001;
	  nb=2.1;
	}
      }
      if (nb_last<10.0) {
	cout << "Pure neutron matter after nb=" << nb_last << endl;
	nst.nb_end=nb_last;
      } else {
	nst.nb_end=2.0;
      }
    }
    nst.verbose=2;
    nst.calc_eos();
    cout << endl;

    cout << "Neutron stars:" << endl;
    nst.calc_nstar();
    std::shared_ptr<table_units<> > te=nst.get_eos_results();
    std::shared_ptr<table_units<> > tr=nst.get_tov_results();

    // Add new columns
    te->line_of_names("msn msp nun nup dnndnun dnpdnup ");
    te->line_of_names("dmundnn dmudn_mixed dmupdnp");
    te->set_unit("msn","1/fm");
    te->set_unit("msp","1/fm");
    te->set_unit("nun","1/fm");
    te->set_unit("nup","1/fm");
    te->set_unit("dnndnun","1/fm^2");
    te->set_unit("dnpdnup","1/fm^2");
    te->set_unit("dmundnn","fm^2");
    te->set_unit("dmudn_mixed","fm^2");
    te->set_unit("dmupdnp","fm^2");
    for(size_t i=0;i<te->get_nlines();i++) {
      n.n=te->get("nn",i);
      p.n=te->get("np",i);
      int ret=sk.calc_e(n,p,th);
      nd=n;
      pd=p;
      thermo_np_deriv_helm thd;
      int ret2=sk.calc_deriv_e(nd,pd,th,thd);
      te->set("msn",i,n.ms);
      te->set("msp",i,p.ms);
      te->set("nun",i,n.nu);
      te->set("nup",i,p.nu);
      te->set("dnndnun",i,nd.dndmu);
      te->set("dnpdnup",i,pd.dndmu);
      te->set("dmundnn",i,thd.dmundnn);
      te->set("dmudn_mixed",i,thd.dmudn_mixed);
      te->set("dmupdnp",i,thd.dmupdnp);
    }
    
    if (output_files) {
      // Output EOS and M vs. R curve to file
      string fn=file_prefix+res.name+"_eos.o2";
      hf.open_or_create(fn.c_str());
      hdf_output(hf,*te,"eos");
      hf.close();
      fn=file_prefix+res.name+"_mvsr.o2";
      hf.open_or_create(fn.c_str());
      hdf_output(hf,*tr,"mvsr");
      hf.close();
    }
      
    if (verbose>0) {
      cout << "M_{max} = " << tr->max("gm") << " R_{max} = "
	   << tr->get("r",tr->lookup("gm",tr->max("gm"))) 
	   << " cent. density = "
	   << tr->get("nb",tr->lookup("gm",tr->max("gm"))) << endl;
      cout << "R_{1.4} = "
	   << tr->get("r",tr->lookup("gm",1.4)) << " cent. density = "
	   << tr->get("nb",tr->lookup("gm",1.4)) << endl;
    }

    res.m_max=tr->max("gm");
    res.r_max=tr->get("r",tr->lookup("gm",tr->max("gm")));
    res.nb_max=tr->get("nb",tr->lookup("gm",tr->max("gm")));
    res.r_14=tr->get("r",tr->lookup("gm",1.4));
    res.nb_14=tr->get("nb",tr->lookup("gm",1.4));

    // Fix the maximum density to check if necessary
    double nbtop=res.nb_max;
    if (nbtop>2.0) nbtop=2.0;
    if (nbtop<0.4) nbtop=0.7;

    // Check for pure neutron matter, now that we have the central
    // density of the maximum mass star. If we find pure neutron
    // matter, set the neutron star parameters to zero since they're
    // likely invalid.
    res.pure_neut=true;
    for(double nb=0.1;nb<=nbtop;nb+=0.01) {
      if (te->get("np",te->lookup("nb",nb))<1.0e-5) {
	res.success=false;
	res.m_max=0.0;
	res.r_max=0.0;
	res.nb_max=0.0;
	res.r_14=0.0;
	res.nb_14=0.0;
	res.pure_neut=false;
      }
    }

    // Check that the maximum mass is larger than 1.6
    res.max_mass_ok=true;
    if (tr->max("gm")<1.6) {
      res.max_mass_ok=false;
      res.success=false;
    }

    // Find where the EOS becomes acausal
    res.acausal_nb=nst.acausal_nb;
    res.pressure_dec_nb=nst.pressure_dec_nb;
    if (nst.pressure_dec_nb>0.0 && 
	nst.pressure_dec_nb<tr->get("nb",tr->lookup("gm",tr->max("gm")))) {
      cout << "Pressure decreases in maximum mass star" << endl;
      cout << "pressure_dec_nb: " << nst.pressure_dec_nb << endl;
      res.success=false;
    }
    if (nst.acausal_nb>0.0 &&
	nst.acausal_nb<tr->get("nb",tr->lookup("gm",tr->max("gm")))) {
      cout << "Acausal before central density of maximum mass star." << endl;
      cout << "acausal_nb: " << nst.acausal_nb << endl;
      res.success=false;
    }

    if (verbose>0) {
      cout << endl;
    }

    return 0;
  }

  /// Summarize the results of one model
  int json(vector<string> &sv, bool itive_com) {
    if (sv.size()<1) {
      cout << "No model to summarize." << endl;
      return exc_einval;
    }

    cout << "Model: " << sv[1] << endl;
    cout << "-----------------------------------------------------" 
	 << "-----------------------" << endl;
    cout << endl;

    skyrme_load(sk,sv[1]);
    name=sv[1];
    res.name=sv[1];
    
    W0hc=sk.W0*hc_mev_fm;
    t0hc=sk.t0*hc_mev_fm;
    t1hc=sk.t1*hc_mev_fm;
    t2hc=sk.t2*hc_mev_fm;
    t3hc=sk.t3*hc_mev_fm;

    nst.verbose=0;
    nst.nb_end=1.18;

    json_mode=true;
    saturation_prop();
    json_mode=false;
    
    return 0;
  }

  /// Summarize the results of one model
  int summary(vector<string> &sv, bool itive_com) {
    if (sv.size()<1) {
      cout << "No model to summarize." << endl;
      return exc_einval;
    }

    cout << "Model: " << sv[1] << endl;
    cout << "-----------------------------------------------------" 
	 << "-----------------------" << endl;
    cout << endl;

    skyrme_load(sk,sv[1]);
    name=sv[1];
    res.name=sv[1];
    
    W0hc=sk.W0*hc_mev_fm;
    t0hc=sk.t0*hc_mev_fm;
    t1hc=sk.t1*hc_mev_fm;
    t2hc=sk.t2*hc_mev_fm;
    t3hc=sk.t3*hc_mev_fm;

    nst.verbose=0;
    nst.nb_end=1.18;

    saturation_prop();
    check_pressure();
    low_neutron_mat();

    compare_neut_nstar();

    mvsr();

    cout << "-----------------------------------------------------" 
	 << "-----------------------" << endl;

    return 0;
  }

  /// Test the code
  int test(vector<string> &sv, bool itive_com) {
    test_mgr t;
    t.set_output_level(1);

    bool of_old=output_files;
    output_files=true;
    
    string file_prefix_old=file_prefix;
    file_prefix="ex_eos_had_skyrme_";

    // Just summarize SLy4
    vector<string> args;
    args.push_back("summary");
    args.push_back("SLy4");
    summary(args,0);

    // And compare to expected results
    t.test_rel(res.n0,0.1595468,1.0e-5,"n0");
    t.test_rel(res.m_max,2.050391,4.0e-4,"m_max");
    t.test_rel(res.r_14,11.72476,4.0e-3,"R_1.4");
    
    t.report();

    output_files=of_old;
    file_prefix=file_prefix_old;

    return 0;
  }

  /// Write to a file
  int store(vector<string> &sv, bool itive_com) {

    if (sv.size()<2) {
      cout << "No filename specified in 'store'." << endl;
      return exc_efailed;
    }

    hf.open_or_create(sv[1]);
    skyrme_write(hf,sk,name);
    hf.close();
    
    if (verbose>0) {
      cout << "Wrote model '" << name << "' to file named '" 
	   << sv[1] << "'." << endl;
    }

    return 0;
  }

  /// Create data files for the UNEDF forces
  int unedf(vector<string> &sv, bool itive_com) {

    test_mgr t;
    t.set_output_level(1);

    /*
      These couplings are in the supplemental material
      for Kortelainen et al. PRC 89 (2014)
     */
    double coups[13][3]={
      {-706.382928878428856,-779.3730087208653,-650.796319465688839},
      {240.049520427681131,287.722131583286796,291.664014339185485},
      {868.871771539645351,891.47789044234969,768.32770588203357},
      {-69.0518957481631617,-200.587774317884879,-283.187292227492492},
      {-12.9172408208016112,-0.989915057807676746,9.78520558824309639},
      {-45.1894169426759476,-33.6320970701835549,-23.3573612977035339},
      {0.321955989588264435,0.2700180115027076,0.351455132555483607},
      {-55.2606,-45.135131022237303,-46.8314091470605973},
      {-55.6226,-145.382167908057,-113.163790795259004},
      {-79.5308,-74.0263331764599,-64.3088624157838069},
      {45.6302,-35.6582611147917,-38.6501946851355029},
      {0.0,0.0,-54.4333635973721002},
      {0.0,0.0,-65.9030310445938028}
    };

    sk.def_neutron.m=hc_mev_fm/20.73553/2.0;
    sk.def_proton.m=hc_mev_fm/20.73553/2.0;
    
    for(size_t i=0;i<3;i++) {
      cout << "UNEDF" << i << ":" << endl;

      // Convert to O2scl units
      for(size_t j=0;j<13;j++) {
	if (j!=6) {
	  coups[j][i]/=hc_mev_fm;
	}
      }

      sk.alt_params_set(coups[0][i],coups[1][i],coups[2][i],coups[3][i],
			coups[4][i],coups[5][i],coups[7][i],coups[8][i],
			coups[9][i],coups[10][i],coups[6][i]);

      // Test coefficients from Kortelainen et al. PRC 85 (2012) 024304
      if (i==0) {
	t.test_rel(sk.t0*hc_mev_fm,-1883.68781034,1.0e-4,"");
	t.test_rel(sk.t1*hc_mev_fm,277.50021224,1.0e-4,"");
	t.test_rel(sk.t2*hc_mev_fm,608.43090559,1.0e-4,"");
	t.test_rel(sk.t3*hc_mev_fm,13901.94834463,1.0e-4,"");
	t.test_rel(sk.x0,0.00974375,1.0e-4,"");
	t.test_rel(sk.x1,-1.77784395,1.0e-4,"");
	t.test_rel(sk.x2,-1.67699035,1.0e-4,"");
	t.test_rel(sk.x3,-0.38079041,1.0e-4,"");
	t.test_rel(sk.b4*hc_mev_fm,125.16100000,1.0e-4,"");
	t.test_rel(sk.b4p*hc_mev_fm,-91.2604000,1.0e-4,"");
      } else if (i==1) {
	t.test_rel(sk.t0*hc_mev_fm,-2078.32802326,1.0e-4,"");
	t.test_rel(sk.t1*hc_mev_fm,239.40081204,1.0e-4,"");
	t.test_rel(sk.t2*hc_mev_fm,1575.11954190,1.0e-4,"");
	t.test_rel(sk.t3*hc_mev_fm,14263.64624708,1.0e-4,"");
	t.test_rel(sk.x0,0.05375692,1.0e-4,"");
	t.test_rel(sk.x1,-5.07723238,1.0e-4,"");
	t.test_rel(sk.x2,-1.36650561,1.0e-4,"");
	t.test_rel(sk.x3,-0.16249117,1.0e-4,"");
	t.test_rel(sk.b4*hc_mev_fm,38.36807206,1.0e-4,"");
	t.test_rel(sk.b4p*hc_mev_fm,71.31652223,1.0e-4,"");
      }

      if (i==0) {
	sk.alt_params_saturation
	  (0.160526,-16.0559/hc_mev_fm,230.0/hc_mev_fm,1.0/0.9,
	   30.5429/hc_mev_fm,45.0804/hc_mev_fm,1.0/1.249838,coups[7][i],
	   coups[8][i],coups[9][i],coups[10][i]);
	
	sk.saturation();
	cout << "n0: " << sk.n0 << endl;
	cout << "E/A: " << sk.eoa*hc_mev_fm << endl;
	cout << "K: " << sk.comp*hc_mev_fm << endl;
	cout << "Esym: " << sk.esym*hc_mev_fm << endl;
	cout << "M_n^*: " << sk.f_effm_neut(sk.n0) << endl;
	cout << "M_p^*: " << sk.f_effm_prot(sk.n0) << endl;
	cout << "1/M_s^*: " << 1.0/sk.f_effm_scalar(sk.n0) << endl;
	cout << "1/M_v^*: " << 1.0/sk.f_effm_vector(sk.n0) << endl;
	cout << "L: " << sk.fesym_slope(sk.n0)*hc_mev_fm << endl;
	cout << "alpha: " << sk.alpha << endl;
      }

      sk.saturation();

      // These numbers from Table IV in Kortelainen et al. (2014)
      if (i==0) {
	t.test_rel(sk.n0,0.16053,1.0e-4,"n0 0");
      } else if (i==1) {
	t.test_rel(sk.n0,0.15871,1.0e-4,"n0 1");
      } else {
	t.test_rel(sk.n0,0.15631,1.0e-4,"n0 2");
      }
      if (i==0) {
	t.test_rel(sk.eoa*hc_mev_fm,-16.056,1.0e-4,"eoa 0");
      } else if (i==1) {
	t.test_rel(sk.eoa*hc_mev_fm,-15.8,1.0e-4,"eoa 1");
      } else {
	t.test_rel(sk.eoa*hc_mev_fm,-15.8,1.0e-4,"eoa 2");
      }
      if (i==0) {
	t.test_rel(sk.fesym_slope(sk.n0)*hc_mev_fm,45.080,1.0e-4,"L 0");
      } else if (i==1) {
	t.test_rel(sk.fesym_slope(sk.n0)*hc_mev_fm,40.005,1.0e-4,"L 1");
      } else {
	t.test_rel(sk.fesym_slope(sk.n0)*hc_mev_fm,40.0,1.0e-4,"L 2");
      }
      
      if (i==0) {
	sk.reference=((string)"M. Kortelainen, T. Lesinski, ")+
	  "J. More, W. Nazarewicz, J. Sarich, N. Schunck, "+
	  "M. V. Stoitsov, and S. Wild, Phys. Rev. C 82, "+
	  "024313 (2010).";
      } else if (i==1) {
	sk.reference=((string)"M. Kortelainen, J. McDonnell, ")+
	  "W. Nazarewicz, P.-G. Reinhard, J. Sarich, N. Schunck, "+
	  "M. V. Stoitsov, and S. M. Wild, Phys. Rev. C 85, "+
	  "024304 (2012).";
      } else {
	sk.reference=((string)"M. Kortelainen, J. McDonnell, ")+
	  "W. Nazarewicz, E. Olsen, P.-G. Reinhard, J. Sarich, "+
	  "N. Schunck, S. M. Wild, D. Davesne, J. Erler, "+
	  "and A. Pastore, Phys. Rev. C 89, 054314 (2014)";
      }

      hf.open_or_create(((string)"UNEDF")+szttos(i)+".o2");
      skyrme_write(hf,sk,((string)"UNEDF")+szttos(i));
      hf.close();
      cout << endl;
    }

    t.report();

    return 0;
  }
  
  /// Load internally stored model
  int load(vector<string> &sv, bool itive_com) {
    
    if (sv.size()<2) {
      cout << "No model specified in 'load'." << endl;
      return exc_efailed;
    }
    
    name=sv[1];
    skyrme_load(sk,name);
    
    if (verbose>0) {
      cout << "Loaded model '" << name << "'." << endl;
    }

    return 0;
  }

  /// Run all the models
  int run_all(vector<string> &sv, bool itive_com) {

    size_t nmods=0;
    std::string mlist[200], stemp;

    string fname=o2scl_settings.get_data_dir()+"/skdata/model_list";
    ifstream fin(fname.c_str());
    fin >> nmods;
    for(size_t i=0;i<nmods;i++) {
      fin >> mlist[i];
    }
    fin.close();

    ofstream fouu("table.csv");
    fouu << "Name, n0, B, K, ";
    fouu << "S, L, Mmax, ";
    fouu << "Rmax, nB_cent_max, R1.4, ";
    fouu << "nB_cent_14, acausal_nb, pressure_dec_nb, ";
    fouu << "neut_qual, max_mass_ok, ";
    fouu << "inc_pressure, pos_neut, good_sat, ";
    fouu << "pure_neut, other, success" << endl;

    ofstream fout("table.html");
    fout << "<html><body>" << endl;
    fout << "<table border=0 cellspacing=0><tr bgcolor=\"#bbbbbb\">" << endl;
    fout << "<td>Name&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>n<sub>0</sub>&nbsp;(fm<sup>-3</sup>)"
	 << "&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>B&nbsp;(MeV)&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>K&nbsp;(MeV)&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>S&nbsp;(MeV)&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>L&nbsp;(MeV)&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>M<sub>max</sub>&nbsp;(M<sub>sun</sub>)"
	 << "&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>R<sub>max</sub>&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>n<sub>B,cent,max</sub>&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>R<sub>1.4</sub>&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>n<sub>B,cent,1.4</sub>&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>acausal_nb&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>pressure_dec_nb&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>neut_qual&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>max_mass_ok&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>inc_pressure&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>pos_neut&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>good_sat&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>pure_neut&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>other&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "<td>success&nbsp;&nbsp;&nbsp;&nbsp;</td>" << endl;
    fout << "</tr>" << endl;
    
    // Currently only models which don't give pure neutron matter work.
    // This list rules out those which have pure neutron matter and
    // the PeHF-type models which don't work.
    static const size_t N=66;
    int list[N]={0,100,101,102,103,109,110,113,114,115,121,122,123,124,
		 125,126,127,128,129,130,131,132,133,134,135,145,147,17,1,
		 25,26,27,28,29,3,40,41,42,43,44,4,51,53,54,5,
		 62,63,64,64,66,67,68,69,6,71,73,75,76,77,7,81,82,84,
		 97,98,99};

    for(size_t j=0;j<N;j++) {

      size_t i=list[j];

      i=2;
      
      cout << "Running model: " << i << endl;
      vector<string> tmp;
      tmp.push_back("x");
      tmp.push_back(mlist[i]);

      // These models are the ones that I used to compute pressure
      // corrections for neutron -> neutron star matter
      // 
      //if (mlist[i]=="SkT3" || mlist[i]=="SLy3" || mlist[i]=="SLy4" ||
      //mlist[i]=="SLy7" || mlist[i]=="SLy230a" || mlist[i]=="SV-mas07" ||
      //mlist[i]=="SV-sym34" || mlist[i]=="BSk9" || mlist[i]=="KDE0v" ||
      //mlist[i]=="KDE0v1" || mlist[i]=="Ly5" || mlist[i]=="mst0.81" ||
      //mlist[i]=="NRAPR" || mlist[i]=="SkMP" || mlist[i]=="BSk14" ||
      //mlist[i]=="SkO" || mlist[i]=="SkOp" || mlist[i]=="SkT1") {

      {

	res.success=true;
	summary(tmp,true);
	exit(-1);

	ofstream fx("jim.dat",ios::app);
	fx.setf(ios::scientific);
	fx.setf(ios::showpos);
	fx.width(10);
	fx << mlist[i] << " ";
	fx << res.S << " " << res.alt_S << " "
	   << res.E_neut_n0 << " " << res.P_neut_n0 << " "
	   << res.L << " ";
	fx << sk.t0 << " " << sk.t1 << " " << sk.t2 << " " << sk.t3 << " ";
	fx << sk.x0 << " " << sk.x1 << " " << sk.x2 << " " << sk.x3 << " ";
	fx << sk.alpha << " " << sk.n0 << " " << sk.comp*hc_mev_fm << " ";
	
	double c2;
	{
	  double h=1.0e-4;
	  n.n=sk.n0-h;
	  p.n=0.0;
	  sk.calc_e(n,p,th);
	  double pr1=th.pr;
	  n.n=sk.n0+h;
	  p.n=0.0;
	  sk.calc_e(n,p,th);
	  double pr2=th.pr;
	  c2=(pr2-pr1)/2.0/h*9.0*hc_mev_fm;
	}
	fx << c2;

	fx << endl; 
	fx.close();

	cout << sk.t0 << " " << sk.t1 << " " << sk.t2 << " " << sk.t3 << endl;
	cout << sk.x0 << " " << sk.x1 << " " << sk.x2 << " " << sk.x3 << endl;
	cout << res.good_sat << endl;
	if (res.success==true) {
	  fout << "<tr bgcolor=\"#dddddd\">";
	} else {
	  fout << "<tr>";
	}

	// Output HTML row
	fout << "<td><a href=\"http://o2scl.svn.sourceforge.net/viewvc/"
	     << "o2scl/trunk/data/o2scl/skdata/" << res.name 
	     << "\">" << res.name << "</a></td>";
	fout << "<td>" << fd.convert(res.n0)  << "</td>";
	fout << "<td>" << fd.convert(res.B)  << "</td>";
	fout << "<td>" << fd.convert(res.K)  << "</td>";
	fout << "<td>" << fd.convert(res.S)  << "</td>";
	fout << "<td>" << fd.convert(res.L)  << "</td>";
	fout << "<td>" << fd.convert(res.m_max)  << "</td>";
	fout << "<td>" << fd.convert(res.r_max)  << "</td>";
	fout << "<td>" << fd.convert(res.nb_max)  << "</td>";
	fout << "<td>" << fd.convert(res.r_14)  << "</td>";
	fout << "<td>" << fd.convert(res.nb_14)  << "</td>";
	fout << "<td>" << fd.convert(res.acausal_nb)  << "</td>";
	fout << "<td>" << fd.convert(res.pressure_dec_nb)  << "</td>";
	fout << "<td>" << fd.convert(res.neut_qual)  << "</td>";
	if (res.max_mass_ok) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	if (res.inc_pressure) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	if (res.pos_neut) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	if (res.good_sat) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	if (res.pure_neut) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	fout << "<td>" << res.other << "</td>";
	if (res.success) fout << "<td>True</td>";
	else fout << "<td>False</td>";
	fout << "</tr>" << endl;

	// Output CSV row
	fouu << res.name << ", ";
	fouu << res.n0 << ", ";
	fouu << res.B << ", ";
	fouu << res.K << ", ";
	fouu << res.S << ", ";
	fouu << res.L << ", ";
	fouu << res.m_max << ", ";
	fouu << res.r_max << ", ";
	fouu << res.nb_max << ", ";
	fouu << res.r_14 << ", ";
	fouu << res.nb_14 << ", ";
	fouu << res.acausal_nb << ", ";
	fouu << res.pressure_dec_nb << ", ";
	fouu << res.neut_qual << ", ";
	if (res.max_mass_ok) fouu << "True, ";
	else fouu << "False, ";
	if (res.inc_pressure) fouu << "True, ";
	else fouu << "False, ";
	if (res.pos_neut) fouu << "True, ";
	else fouu << "False, ";
	if (res.good_sat) fouu << "True, ";
	else fouu << "False, ";
	if (res.pure_neut) fouu << "True, ";
	else fouu << "False, ";
	fouu << res.other << ", ";
	if (res.success) fouu << "True ";
	else fouu << "False ";
	fouu << endl;

      }
    }

    fout << "</table></body></html>" << endl;
    fout.close();

    fouu.close();
    

    return 0;
  }

};


int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  
  ex_eos_had_skyrme se;

  // ---------------------------------------
  // Specify command-line option object

  cli cl;
  cl.prompt="ex_eos_had_skyrme> ";

  int comm_option_cl_param=1;
  int comm_option_both=2;

  static const int narr=7;
  comm_option_s options_arr[narr]={
    {0,"run-all","Run all internally stored Skyrme models.",0,0,"","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::run_all),
     comm_option_both},
    {'s',"store","Store current model.",1,1,"","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::store),
     comm_option_both},
    {'l',"load","Load internally stored model.",1,1,"","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::load),
     comm_option_both},
    {0,"unedf","Desc.",0,0,"","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::unedf),
     comm_option_both},
    {'t',"test","Test ex_eos_had_skyrme.",0,0,"",
     ((string)"This command temporarily sets output_files to true, ")+
     "file_prefix to \"ex_eos_had_skyrme_\", and then runs the "+
     "'summary' command for SLy4.",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::test),
     comm_option_both},
    {'u',"summary","Summarize the properties of a Skyrme model.",
     1,1,"<model>","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::summary),
     comm_option_both},
    {'j',"json","Summarize the properties of a Skyrme model.",
     1,1,"<model>","",
     new comm_option_mfptr<ex_eos_had_skyrme>(&se,&ex_eos_had_skyrme::json),
     comm_option_both}
  };

  cl.set_comm_option_vec(narr,options_arr);
  cl.cmd_name="ex_eos_had_skyrme";

  // ---------------------------------------
  // Set the parameters

  cli::parameter_int p_verbose;
  p_verbose.i=&se.verbose;
  p_verbose.help="Verbose (default 1).";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  cli::parameter_bool p_output_files;
  p_output_files.b=&se.output_files;
  p_output_files.help="Output files (default 0).";
  cl.par_list.insert(make_pair("output_files",&p_output_files));

  cli::parameter_string p_file_prefix;
  p_file_prefix.str=&se.file_prefix;
  p_file_prefix.help="File prefix (default \"\").";
  cl.par_list.insert(make_pair("file_prefix",&p_file_prefix));

  cli::parameter_string p_name;
  p_name.str=&se.name;
  p_name.help="Model name (default \"\").";
  cl.par_list.insert(make_pair("name",&p_name));

  cli::parameter_string p_reference;
  p_reference.str=&se.sk.reference;
  p_reference.help="Model reference (default \"\").";
  cl.par_list.insert(make_pair("reference",&p_reference));

  cli::parameter_double p_t0hc;
  p_t0hc.d=&se.t0hc;
  p_t0hc.help="Model parameter t0 in MeV.";
  cl.par_list.insert(make_pair("t0hc",&p_t0hc));

  cli::parameter_double p_t1hc;
  p_t1hc.d=&se.t1hc;
  p_t1hc.help="Model parameter t1 in MeV.";
  cl.par_list.insert(make_pair("t1hc",&p_t1hc));

  cli::parameter_double p_t2hc;
  p_t2hc.d=&se.t2hc;
  p_t2hc.help="Model parameter t2 in MeV.";
  cl.par_list.insert(make_pair("t2hc",&p_t2hc));

  cli::parameter_double p_t3hc;
  p_t3hc.d=&se.t3hc;
  p_t3hc.help="Model parameter t3 in MeV.";
  cl.par_list.insert(make_pair("t3hc",&p_t3hc));

  cli::parameter_double p_x0;
  p_x0.d=&se.sk.x0;
  p_x0.help="Model parameter x0.";
  cl.par_list.insert(make_pair("x0",&p_x0));

  cli::parameter_double p_x1;
  p_x1.d=&se.sk.x1;
  p_x1.help="Model parameter x1.";
  cl.par_list.insert(make_pair("x1",&p_x1));

  cli::parameter_double p_x2;
  p_x2.d=&se.sk.x2;
  p_x2.help="Model parameter x2.";
  cl.par_list.insert(make_pair("x2",&p_x2));

  cli::parameter_double p_x3;
  p_x3.d=&se.sk.x3;
  p_x3.help="Model parameter x3.";
  cl.par_list.insert(make_pair("x3",&p_x3));

  cli::parameter_double p_a;
  p_a.d=&se.sk.a;
  p_a.help="Model parameter a.";
  cl.par_list.insert(make_pair("a",&p_a));

  cli::parameter_double p_b;
  p_b.d=&se.sk.b;
  p_b.help="Model parameter b.";
  cl.par_list.insert(make_pair("b",&p_b));

  cli::parameter_double p_W0hc;
  p_W0hc.d=&se.W0hc;
  p_W0hc.help="Model parameter W0hc.";
  cl.par_list.insert(make_pair("W0hc",&p_W0hc));

  cli::parameter_double p_alpha;
  p_alpha.d=&se.sk.alpha;
  p_alpha.help="Model parameter alpha.";
  cl.par_list.insert(make_pair("alpha",&p_alpha));

  cli::parameter_double p_b4;
  p_b4.d=&se.sk.b4;
  p_b4.help="Model parameter b4.";
  cl.par_list.insert(make_pair("b4",&p_b4));

  cli::parameter_double p_b4p;
  p_b4p.d=&se.sk.b4p;
  p_b4p.help="Model parameter b4p.";
  cl.par_list.insert(make_pair("b4p",&p_b4p));

  // ---------------------------------------
  // Process command-line options

  cl.run_auto(argc,argv);
  
  return 0;
}

