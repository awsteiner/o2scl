/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

#include <o2scl/nucleus_rmf.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucleus_rmf::nucleus_rmf() : 
  profiles(new table_units<>),chden_table(new table_units<>)  {

  verbose=1;
  err_nonconv=true;
  generic_ode=false;

  fields.resize(grid_size,4);
  field0.resize(grid_size,3);
  xrho.resize(grid_size,4);
  xrhosp.resize(grid_size);
  gin.resize(grid_size,4);
  gout.resize(grid_size,4);
  xrhos.resize(grid_size);
  xrhov.resize(grid_size);
  xrhor.resize(grid_size);
  chden1.resize(grid_size);
  chdenc.resize(grid_size);

  ode_y.resize(2);
  ode_dydx.resize(2);
  ode_yerr.resize(2);

  energy.resize(grid_size+6);
  arho.resize(grid_size+6);

  shell pshells[n_internal_levels]=
    {
      {2,-1,-60.0,0.5,"^{1}s_{1/2}",4.0,0,0,0},
      // Z=2
      {4,-2,1.0,0.5,"^{1}p_{3/2}",4.0,0,0,0},
      {2,1,1.0,0.5,"^{1}p_{1/2}",4.0,0,0,0},
      // Z=8
      {6,-3,1.0,0.5,"^{1}d_{5/2}",5.0,0,0,0},
      {4,2,1.0,0.5,"^{1}d_{3/2}",5.0,0,0,0},
      {2,-1,1.0,0.5,"^{2}s_{1/2}",5.0,0,0,0},
      // Z=20
      {8,-4,1.0,0.5,"^{1}f_{7/2}",6.0,0,0,0},
      // Z=28
      {6,3,1.0,0.5,"^{1}f_{5/2}",6.0,0,0,0},
      {4,-2,1.0,0.5,"^{2}p_{3/2}",6.0,0,0,0},
      {2,1,1.0,0.5,"^{2}p_{1/2}",6.0,0,0,0},
      // Z=40
      {10,-5,1.0,0.5,"^{1}g_{9/2}",6.0,0,0,0},
      // Z=50
      {8,4,1.0,0.5,"^{1}g_{7/2}",6.0,0,0,0},
      {6,-3,1.0,0.5,"^{2}d_{5/2}",7.0,0,0,0},
      {12,-6,1.0,0.5,"^{1}h_{11/2}",6.0,0,0,0}, 
      {4,2,1.0,0.5,"^{2}d_{3/2}",7.0,0,0,0},
      {2,-1,1.0,0.5,"^{3}s_{1/2}",7.0,0,0,0},
      // Z=82
      {10,5,-4.0,0.5,"^{1}h_{9/2}",6.0,0,0,0}, 
      {8,-4,-3.0,0.5,"^{2}f_{7/2}",7.0,0,0,0}, 
      {14,-7,-2.0,0.5,"^{1}i_{13/2}",7.0,0,0,0},
      {6,3,-0.5,0.5,"^{2}f_{5/2}",7.0,0,0,0},
      {4,-2,-1.0,0.5,"^{3}p_{3/2}",7.0,0,0,0},
      {2,1,1.0,-0.5,"^{3}p_{1/2}",7.0,0,0,0},
      // Z=126
      {10,-5,-4.0,-0.5,"^{2}g_{9/2}",7.0,0,0,0},
      {12,6,-3.0,-0.5,"^{1}i_{11/2}",7.0,0,0,0},
      {16,-8,-2.5,-0.5,"^{1}j_{15/2}",7.0,0,0,0},
      {6,-3,-2.5,-0.5,"^{3}d_{5/2}",7.0,0,0,0},
      {2,-1,-2.0,-0.5,"^{4}s_{1/2}",7.0,0,0,0},
      {8,4,-1.5,-0.5,"^{2}g_{7/2}",7.0,0,0,0},
      {4,2,-1.5,-0.5,"^{3}d_{3/2}",7.0,0,0,0}
      // Z=184
    };
  shell nshells[n_internal_levels]=
    {
      {2,-1,-60.0,-0.5,"^{1}s_{1/2}",3.0,0,0,0},
      // N=2
      {4,-2,1.0,-0.5,"^{1}p_{3/2}",4.0,0,0,0},
      {2,1,1.0,-0.5,"^{1}p_{1/2}",4.0,0,0,0},
      // N=8
      {6,-3,1.0,-0.5,"^{1}d_{5/2}",5.0,0,0,0},
      {4,2,1.0,-0.5,"^{1}d_{3/2}",5.0,0,0,0},
      {2,-1,1.0,-0.5,"^{2}s_{1/2}",5.0,0,0,0},
      // N=20
      {8,-4,1.0,-0.5,"^{1}f_{7/2}",6.0,0,0,0},
      // N=28
      {6,3,1.0,-0.5,"^{1}f_{5/2}",6.0,0,0,0},
      {4,-2,1.0,-0.5,"^{2}p_{3/2}",6.0,0,0,0},
      {2,1,1.0,-0.5,"^{2}p_{1/2}",6.0,0,0,0},
      // N=40
      {10,-5,1.0,-0.5,"^{1}g_{9/2}",6.0,0,0,0},
      // N=50
      {8,4,1.0,-0.5,"^{1}g_{7/2}",6.0,0,0,0},
      {6,-3,1.0,-0.5,"^{2}d_{5/2}",6.0,0,0,0},
      {12,-6,1.0,-0.5,"^{1}h_{11/2}",6.0,0,0,0}, 
      {4,2,1.0,-0.5,"^{2}d_{3/2}",6.0,0,0,0},
      {2,-1,1.0,-0.5,"^{3}s_{1/2}",6.0,0,0,0},
      // N=82
      {10,5,1.0,0.5,"^{1}h_{9/2}",6.0,0,0,0}, 
      {8,-4,1.0,-0.5,"^{2}f_{7/2}",7.0,0,0,0}, 
      {14,-7,1.0,-0.5,"^{1}i_{13/2}",7.0,0,0,0},
      {6,3,1.0,-0.5,"^{2}f_{5/2}",7.0,0,0,0},
      {4,-2,1.0,-0.5,"^{3}p_{3/2}",7.0,0,0,0},
      {2,1,1.0,-0.5,"^{3}p_{1/2}",7.0,0,0,0},
      // N=126
      {10,-5,-4.0,-0.5,"^{2}g_{9/2}",7.0,0,0,0},
      {12,6,-3.0,-0.5,"^{1}i_{11/2}",7.0,0,0,0},
      {16,-8,-2.5,-0.5,"^{1}j_{15/2}",7.0,0,0,0},
      {6,-3,-2.5,-0.5,"^{3}d_{5/2}",7.0,0,0,0},
      {2,-1,-2.0,-0.5,"^{4}s_{1/2}",7.0,0,0,0},
      {8,4,-1.5,-0.5,"^{2}g_{7/2}",7.0,0,0,0},
      {4,2,-1.5,-0.5,"^{3}d_{3/2}",7.0,0,0,0}
      // N=184
    };

  for(int k=0;k<n_internal_levels;k++) neutron_shells[k]=nshells[k];
  for(int k=0;k<n_internal_levels;k++) proton_shells[k]=pshells[k];

  // Allocate memory for the vector<shells> objects.
  // Allocate too much memory here to be safe.
  for(size_t i=0;i<2*n_internal_levels;i++) {
    levels.push_back(pshells[0]);
    unocc_levels.push_back(pshells[0]);
  }

  ostep=&def_step;

  rmf=&def_rmf;

  itmax=70;
  meson_itmax=10000;
  dirac_itmax=100;
  meson_tol=1.0e-6;

  step_size=0.04;

  dirac_tol=5.0e-3;
  dirac_tol2=5.0e-4;

  profiles->line_of_names(((string)"r rhop rhon sig ome rho ")+
			 "coul chden rhosp rhosn ");
  profiles->set_unit("r","fm");
  profiles->set_unit("rhop","1/fm^3");
  profiles->set_unit("rhon","1/fm^3");
  profiles->set_unit("sig","1/fm");
  profiles->set_unit("ome","1/fm");
  profiles->set_unit("rho","1/fm");
  profiles->set_unit("rhosp","1/fm^3");
  profiles->set_unit("rhosn","1/fm^3");
  chden_table->line_of_names("x chden1 chdenc");
  chden_table->set_unit("chden1","1/fm^3");
  chden_table->set_unit("chdenc","1/fm^3");
  
  init_called=false;
  
  initial_guess ig_tmp={310/hc_mev_fm,240/hc_mev_fm,
			-6/hc_mev_fm,25.9/hc_mev_fm,6.85,0.6};
  ig=ig_tmp;

  // From gsl-1.15/integration/gslfixed.c
  x12[0]=0.1252334085114689154724414;
  x12[1]=0.3678314989981801937526915;
  x12[2]=0.5873179542866174472967024;
  x12[3]=0.7699026741943046870368938;
  x12[4]=0.9041172563704748566784659;
  x12[5]=0.9815606342467192506905491;
  w12[0]=0.2491470458134027850005624;
  w12[1]=0.2334925365383548087608499;
  w12[2]=0.2031674267230659217490645;
  w12[3]=0.1600783285433462263346525;
  w12[4]=0.1069393259953184309602547;
  w12[5]=0.0471753363865118271946160;

  x100[0]=0.0156289844215430828722167;
  x100[1]=0.0468716824215916316149239; 
  x100[2]=0.0780685828134366366948174;
  x100[3]=0.1091892035800611150034260; 
  x100[4]=0.1402031372361139732075146;
  x100[5]=0.1710800805386032748875324; 
  x100[6]=0.2018898640957359972360489; 
  x100[7]=0.2323024818449739696495100; 
  x100[8]=0.2625881203715034791689293;
  x100[9]=0.2926171880384719647375559; 
  x100[10]=0.3223603439005291517224766;
  x100[11]=0.3517885263724217209723438; 
  x100[12]=0.3808729816246299567633625;
  x100[13]=0.4095852916783015425288684; 
  x100[14]=0.4378974021720315131089780;
  x100[15]=0.4657816497733580422492166; 
  x100[16]=0.4932107892081909335693088;
  x100[17]=0.5201580198817630566468157; 
  x100[18]=0.5465970120650941674679943;
  x100[19]=0.5725019326213811913168704; 
  x100[20]=0.5978474702471787212648065; 
  x100[21]=0.6226088602037077716041908; 
  x100[22]=0.6467619085141292798326303;
  x100[23]=0.6702830156031410158025870; 
  x100[24]=0.6931491993558019659486479;
  x100[25]=0.7153381175730564464599671; 
  x100[26]=0.7368280898020207055124277;
  x100[27]=0.7575981185197071760356680; 
  x100[28]=0.7776279096494954756275514;
  x100[29]=0.7968978923903144763895729; 
  x100[30]=0.8153892383391762543939888;
  x100[31]=0.8330838798884008235429158; 
  x100[32]=0.8499645278795912842933626;
  x100[33]=0.8660146884971646234107400; 
  x100[34]=0.8812186793850184155733168;
  x100[35]=0.8955616449707269866985210; 
  x100[36]=0.9090295709825296904671263;
  x100[37]=0.9216092981453339526669513; 
  x100[38]=0.9332885350430795459243337;
  x100[39]=0.9440558701362559779627747; 
  x100[40]=0.9539007829254917428493369;
  x100[41]=0.9628136542558155272936593; 
  x100[42]=0.9707857757637063319308979; 
  x100[43]=0.9778093584869182885537811; 
  x100[44]=0.9838775407060570154961002;
  x100[45]=0.9889843952429917480044187; 
  x100[46]=0.9931249370374434596520099;
  x100[47]=0.9962951347331251491861317; 
  x100[48]=0.9984919506395958184001634;
  x100[49]=0.9997137267734412336782285;
  
  w100[0]=0.0312554234538633569476425; 
  w100[1]=0.0312248842548493577323765;
  w100[2]=0.0311638356962099067838183; 
  w100[3]=0.0310723374275665165878102;
  w100[4]=0.0309504788504909882340635; 
  w100[5]=0.0307983790311525904277139;
  w100[6]=0.0306161865839804484964594; 
  w100[7]=0.0304040795264548200165079;
  w100[8]=0.0301622651051691449190687; 
  w100[9]=0.0298909795933328309168368; 
  w100[10]=0.0295904880599126425117545; 
  w100[11]=0.0292610841106382766201190;
  w100[12]=0.0289030896011252031348762; 
  w100[13]=0.0285168543223950979909368;
  w100[14]=0.0281027556591011733176483; 
  w100[15]=0.0276611982207923882942042;
  w100[16]=0.0271926134465768801364916; 
  w100[17]=0.0266974591835709626603847; 
  w100[18]=0.0261762192395456763423087; 
  w100[19]=0.0256294029102081160756420;
  w100[20]=0.0250575444815795897037642; 
  w100[21]=0.0244612027079570527199750;
  w100[22]=0.0238409602659682059625604; 
  w100[23]=0.0231974231852541216224889;
  w100[24]=0.0225312202563362727017970; 
  w100[25]=0.0218430024162473863139537;
  w100[26]=0.0211334421125276415426723; 
  w100[27]=0.0204032326462094327668389;
  w100[28]=0.0196530874944353058653815; 
  w100[29]=0.0188837396133749045529412;
  w100[30]=0.0180959407221281166643908; 
  w100[31]=0.0172904605683235824393442; 
  w100[32]=0.0164680861761452126431050; 
  w100[33]=0.0156296210775460027239369;
  w100[34]=0.0147758845274413017688800; 
  w100[35]=0.0139077107037187726879541;
  w100[36]=0.0130259478929715422855586; 
  w100[37]=0.0121314576629794974077448;
  w100[38]=0.0112251140231859771172216; 
  w100[39]=0.0103078025748689695857821;
  w100[40]=0.0093804196536944579514182; 
  w100[41]=0.0084438714696689714026208;
  w100[42]=0.0074990732554647115788287; 
  w100[43]=0.0065469484508453227641521;
  w100[44]=0.0055884280038655151572119; 
  w100[45]=0.0046244500634221193510958;
  w100[46]=0.0036559612013263751823425; 
  w100[47]=0.0026839253715534824194396;
  w100[48]=0.0017093926535181052395294; 
  w100[49]=0.0007346344905056717304063;

  load_nl3(def_rmf);

  a_proton=sqrt(0.71)*1.0e3/o2scl_const::hc_mev_fm;
}

int nucleus_rmf::load_nl3(eos_had_rmf &r) {

  r.ms=508.194;
  r.mw=782.501;
  r.mr=763.0;
  r.mnuc=939.0;
  r.ms/=hc_mev_fm; 
  r.mw/=hc_mev_fm; 
  r.mr/=hc_mev_fm; 
  r.mnuc/=hc_mev_fm;
    
  double gs, gw, gr;
  gs=10.217;
  gw=12.868;
  gr=4.474;
  r.b=-10.431;
  r.c=-28.885;
  r.b/=-r.mnuc*pow(fabs(gs),3.0);
  r.c/=pow(gs,4.0);
  gr*=2.0;
  r.cs=gs/r.ms;
  r.cw=gw/r.mw;
  r.cr=gr/r.mr;
    
  r.xi=0.0; 
  r.zeta=0.0;
  r.a1=0.0;
  r.a2=0.0;
  r.a3=0.0;
  r.a4=0.0;
  r.a5=0.0;
  r.a6=0.0;
  r.b1=0.0;
  r.b2=0.0;
  r.b3=0.0;
    
  return 0;
}

nucleus_rmf::~nucleus_rmf() {
  fields.clear();
  xrho.clear();
  xrhosp.clear();
  gin.clear();
  gout.clear();
  xrhos.clear();
  xrhov.clear();
  xrhor.clear();
  chden1.clear();
  chdenc.clear();
  field0.clear();
  energy.clear();
  arho.clear();

  levels.clear();
  unocc_levels.clear();
}

int nucleus_rmf::run_nucleus(int nucleus_Z, int nucleus_N,
			      int unocc_Z, int unocc_N) {
  
  init_run(nucleus_Z,nucleus_N,unocc_Z,unocc_N);

  int iteration=1, iconverged=0;

  /*
    We allow the Dirac equations and meson field equations to 
    not converge temporarily, but require that they successfully
    converged by the time iconverged is nonzero.
  */
  int dirac_converged=0, meson_converged=0;
  
  while (iconverged==0) {
    
    if (iteration>itmax) {
      O2SCL_CONV_RET((((string)"Failed to converge after ")+
		      itos(itmax)+" iterations.").c_str(),exc_emaxiter,
		     err_nonconv);
    }
    
    int iret=iterate(nucleus_Z,nucleus_N,unocc_Z,unocc_N,iconverged,
		     dirac_converged,meson_converged);

    if (verbose>0) {
      cout << "Iteration: " << iteration << " ret: " << iret
	   << " dirac: " << dirac_converged
	   << " meson: " << meson_converged << endl;
    }
    
    if (iret!=0) {
      O2SCL_CONV_RET("Function iterate() failed.",exc_efailed,
		     err_nonconv);
    }
    
    iteration++;
  }

  if (dirac_converged!=0) {
    O2SCL_CONV2_RET("Dirac equations did not converge during final ",
		    "iteration in nucleus_rmf::run_nucleus().",
		    exc_efailed,err_nonconv);
  }
  if (meson_converged!=0) {
    O2SCL_CONV2_RET("Meson field equations did not converge during final ",
		    "iteration in nucleus_rmf::run_nucleus().",
		    exc_efailed,err_nonconv);
  }
  
  post_converge(nucleus_Z,nucleus_N,unocc_Z,unocc_N);
  
  return success;
}

void nucleus_rmf::init_run(int nucleus_Z, int nucleus_N,
			   int unocc_Z, int unocc_N) {

  int i, iteration;

  // The number of (occupied) neutron, proton, and total and
  // unoccupied states
  int nistate, pistate;
  nuolevels=unocc_Z+unocc_N;
  
  n.non_interacting=false;
  p.non_interacting=false;

  profiles->clear_data();
  chden_table->clear_data();

  //--------------------------------------------------------------
  // Assign the appropriate neutron and proton 
  // levels

  nlevels=0, i=0;
  int nprotons=nucleus_Z, iprot=0;
  int nneutrons=nucleus_N, ineut=0;
  while(iprot<nprotons && i<n_internal_levels) {
    levels[nlevels]=proton_shells[i];
    iprot+=levels[nlevels].twojp1;
    levels[nlevels].isospin=0.5;
    nlevels++;
    i++;
  }
  pistate=nlevels-1;
  i=0;
  while(ineut<nneutrons && i<n_internal_levels) {
    levels[nlevels]=neutron_shells[i];
    ineut+=levels[nlevels].twojp1;
    levels[nlevels].isospin=-0.5;
    nlevels++;
    i++;
  }
  nistate=nlevels-pistate;

  if (ineut!=nneutrons || iprot!=nprotons) {
    O2SCL_ERR2("Not a closed-shell nucleus or too many nucleons in ",
	       "nucleus_rmf::init_run().",exc_einval);
  }
  
  //--------------------------------------------------------------
  // Assign unoccupied levels
  
  for(int ile=0;ile<unocc_Z;ile++) {
    if (ile+pistate>=n_internal_levels) {
      O2SCL_ERR2("Requested too many unoccupied protons in ",
		 "nucleus_rmf::init_run().",exc_einval);
    }
    unocc_levels[ile]=proton_shells[ile+pistate];
  }
  for(int ile=unocc_Z;ile<unocc_Z+unocc_N;ile++) {
    if (ile+nistate>=n_internal_levels) {
      O2SCL_ERR2("Requested too many unoccupied neutrons in ",
		 "nucleus_rmf::init_run().",exc_einval);
    }
    unocc_levels[ile]=neutron_shells[ile-unocc_Z+nistate];
  }
  if (nuolevels==0) {
    if (verbose>1) {
      cout << "No unoccupied levels." << endl;
    }
  }
  
  //--------------------------------------------------------------
  // Initialize some masses, coupling constants
  // and convergence parameters
  
  n.init(rmf->mnuc,2.0);
  p.init(rmf->mnuc,2.0);

  mnuc=rmf->mnuc;
  
  //--------------------------------------------------------------
  // Initial guess for fields

  for (i=0;i<grid_size;i++) {
    double ex=exp((step_size*((double)(i+1))-ig.fermi_radius)/ig.fermi_width);
    fields(i,0)=ig.sigma0/(1.0+ex);
    fields(i,1)=ig.omega0/(1.0+ex);
    fields(i,2)=ig.rho0/(1.0+ex);
    fields(i,3)=ig.A0/(1.0+ex);
    field0(i,0)=fields(i,0);
    field0(i,1)=fields(i,1);
    field0(i,2)=fields(i,2);
  }
  surf_index=((int)(ig.fermi_radius/step_size+1.0e-6));
  
  //--------------------------------------------------------------

  // Initialize meson greens functions
  meson_init();

  init_called=true;

  return;
}

int nucleus_rmf::iterate(int nucleus_Z, int nucleus_N,
			 int unocc_Z, int unocc_N, int &iconverged,
			 int &dirac_converged, int &meson_converged) {
  iconverged=0;
  
  if (init_called==false) {
    O2SCL_ERR2("Function init_run() has not been called in ",
	       "nucleus_rmf::iterate().",exc_efailed);
  }
  
  //--------------------------------------------------------------
  // Calculate new fields, densities, and RHSs

  init_meson_density();
  
  levp=&levels;

  //--------------------------------------------------------------
  // Solve Dirac equations
    
  if (verbose>1) cout << "Solving Dirac equations. " << endl;
  for (int ilevel=0;ilevel<nlevels;ilevel++) {
    int dret=dirac(ilevel);
    if (dret!=0) dirac_converged=dret;
    (*levp)[ilevel].eigenc=(*levp)[ilevel].eigen-(*levp)[ilevel].energy;
  }
    
  //--------------------------------------------------------------
  // Test the Dirac equations for convergence
	
  iconverged=1;
  for (int i=0;i<nlevels;i++) {
    if (fabs((*levp)[i].eigenc)>dirac_tol) iconverged=0;
  }

  //--------------------------------------------------------------
  // Solve the meson field equations

  meson_converged=meson_solve();

  //--------------------------------------------------------------
  // Modify densities according to solution of
  // meson field equations

  for (int i=0;i<grid_size;i++) {
    double xr=pow(((double)(i+1))*step_size,2.0);
    xrho(i,0)=xrho(i,0)-xrhos[i]*xr;
    xrho(i,1)=xrho(i,1)-xrhov[i]*xr;
    xrho(i,2)=xrho(i,2)-xrhor[i]*xr;
  }
    
  //--------------------------------------------------------------
  // Record iteration information

  profiles->clear_data();

  for (int i=0;i<grid_size;i++) {

    double x=((double)(i+1))*step_size;

    // If we're done, fold proton form factor
    if (iconverged==1) {
      pfold(x,chden1[i]);
    }

    // baryon densities
    double xp=xrho(i,3)/x/x;
    double xn=(xrho(i,1)-xrho(i,3))/x/x;

    // scalar densities 
    double xnns=(xrho(i,0)-xrhosp[i])/x/x;
    double xpns=xrhosp[i]/x/x;

    if (!std::isfinite(xn) ||
	!std::isfinite(xp)) {
      return 1;
    }
    
    double line[10]={x,xp,xn,fields(i,0),fields(i,1),
		     fields(i,2),fields(i,3),chden1[i],
		     xpns,xnns};
    profiles->line_of_data(10,line);
    
  }

  //--------------------------------------------------------------
  // Compute the energy, first by summing the single 
  // particle energies, and then by adding the mean-field
  // contribution

  double sum_sp_energies=0.0;
  for (int i=0;i<nlevels;i++) {
    sum_sp_energies+=(*levp)[i].twojp1*(*levp)[i].eigen;
  }

  //--------------------------------------------------------------
  // Compute total energy and radii
  
  int eret=energy_radii(nucleus_Z,nucleus_N,sum_sp_energies);
  if (eret!=0) return eret;
      
  //--------------------------------------------------------------
  // Set eigenvalues for next iteration
    
  for(int i=0;i<nlevels;i++) {
    (*levp)[i].eigenc=0.0;
    (*levp)[i].energy=(*levp)[i].eigen;
  }

  return 0;
}

int nucleus_rmf::post_converge(int nucleus_Z, int nucleus_N, int unocc_Z, 
			       int unocc_N) {

  //--------------------------------------------------------------
  // Correct for center of mass motion 
  
  center_mass_corr(nucleus_Z+nucleus_N);
  
  //--------------------------------------------------------------
  // Compute unoccupied levels
  
  if (nuolevels>0) {

    levp=&unocc_levels;
    
    //--------------------------------------------------------------
    // Calculate new fields, densities, and RHSs
    
    init_meson_density();
    
    //--------------------------------------------------------------
    // Solve Dirac equations
    
    if (verbose>1) cout << "Solving Dirac equations for unoccupied levels. " 
			<< endl;
    
    for (int ilevel=0;ilevel<nuolevels;ilevel++) {
      int dret=dirac(ilevel);
      if (dret!=0) return dret;
      (*levp)[ilevel].eigenc=(*levp)[ilevel].eigen-(*levp)[ilevel].energy;
    }

  }
      
  //--------------------------------------------------------------
  // Calculate surface tension
  
  // Old integration constants
  ubvector fac;
  
  fac.resize(4);

  fac[0]=3.0/8.0;
  fac[1]=7.0/6.0-3.0/8.0;
  fac[2]=23.0/24.0-7.0/6.0;
  fac[3]=1.0/24.0;

  for(int i=4;i<295;i++) {
    for(int j=0;j<4;j++) {
      stens+=step_size*(energy[i+j-1]-energy[2]*arho[i+j-1]/arho[2])*fac[j];
    }
  }
  if (!std::isfinite(stens)) {
    stens=0.0;
  }
  
  return 0;
}

int nucleus_rmf::meson_solve() {
  int i;

  meson_iter(4);

  double t1=xrho(0,0);
  double t2=xrho(surf_index-1,0);
  double t3=xrho(0,1);
  double t4=xrho(surf_index-1,1);
  double t5=xrho(0,2);
  double xitno=0.0;
  
  if (verbose>1) cout << "Solving meson field equations. " << endl;
	
  bool mesondone=false;
  int mesonfieldcount=0;
  while(mesondone==false) {

    for (i=0;i<grid_size;i++) {
      double xr=pow(((double)(i+1))*step_size,2.0);
      xrho(i,0)=xrho(i,0)-xrhos[i]*xr;
      xrho(i,1)=xrho(i,1)-xrhov[i]*xr;
      xrho(i,2)=xrho(i,2)-xrhor[i]*xr;
      xrhos[i]=(sigma_rhs(fields(i,0),fields(i,1),fields(i,2))
		+xitno*xrhos[i])/(xitno+1.0);
      xrhov[i]=(omega_rhs(fields(i,0),fields(i,1),fields(i,2))
		+xitno*xrhov[i])/(xitno+1.0);
      xrhor[i]=(rho_rhs(fields(i,0),fields(i,1),fields(i,2))
		+xitno*xrhor[i])/(xitno+1.0);
      xrho(i,0)=xrho(i,0)+xrhos[i]*xr;
      xrho(i,1)=xrho(i,1)+xrhov[i]*xr;
      xrho(i,2)=xrho(i,2)+xrhor[i]*xr;
    }
    
    mesondone=true;
    if (fabs(t1/xrho(0,0)-1.0)>meson_tol) mesondone=false;
    if (fabs(t2/xrho(surf_index-1,0)-1.0)>meson_tol) mesondone=false;
    if (fabs(t3/xrho(0,1)-1.0)>meson_tol) mesondone=false;
    if (fabs(t4/xrho(surf_index-1,1)-1.0)>meson_tol) mesondone=false;
    if (fabs(t5/xrho(0,2)-1.0)>meson_tol) mesondone=false;

    if (mesondone==false) {
      xitno=1.0;
      t1=xrho(0,0);
      t2=xrho(surf_index-1,0);
      t3=xrho(0,1);
      t4=xrho(surf_index-1,1);
      t5=xrho(0,2);

      meson_iter(3);
    }

    mesonfieldcount++;
    if (mesonfieldcount>meson_itmax) {
      O2SCL_CONV2_RET("Failed to solve meson field equations in ",
		      "nucleus_rmf::meson_solve().",exc_efailed,err_nonconv);
    }

  }

  return 0;
}

int nucleus_rmf::energy_radii(double xpro, double xnu, double e) {

  rprms=0.0;
  rnrms=0.0;

  double etemp=0.0;

  for (int i=0;i<grid_size;i++) {
    double x=((double)(i+1))*step_size;
    etemp+=fields(i,0)*xrho(i,0)-fields(i,1)*xrho(i,1)-
      fields(i,2)*xrho(i,2)-fields(i,3)*xrho(i,3);
    double sig=fields(i,0)/rmf->ms/rmf->cs;
    double ome=fields(i,1)/rmf->mw/rmf->cw;
    double pefac=rmf->a1*sig+2.0*rmf->a2*sig*sig+
      3.0*rmf->a3*pow(sig,3.0)+4.0*rmf->a4*pow(sig,4.0)+
      5.0*rmf->a5*pow(sig,5.0)+6.0*rmf->a6*pow(sig,6.0)+
      2.0*rmf->b1*pow(ome,2.0)+4.0*rmf->b2*pow(ome,4.0)+
      6.0*rmf->b3*pow(ome,6.0);
    etemp+=(rmf->zeta/12.0*pow(fields(i,1),4.0)+
	    rmf->xi/12.0*pow(fields(i,2),4.0)+
	    fields(i,2)*fields(i,2)*pefac-
	    rmf->b*rmf->mnuc/3.0*pow(fields(i,0),3.0)-
	    rmf->c/2.0*pow(fields(i,0),4.0))*x*x;
    if (std::isnan(etemp)) {
      // This quantity can become NaN if the computation fails
      // to converge, so we don't designate this as a fatal
      // error and only throw if err_nonconv is true.
      O2SCL_CONV2_RET("Energy contribution is not finite in ",
		      "nucleus_rmf::energy_radii().",exc_efailed,err_nonconv);
    }
    rprms=rprms+x*x*xrho(i,3);
    rnrms=rnrms+x*x*(xrho(i,1)-xrho(i,3));
  }

  etot=e/(xpro+xnu)+etemp*2.0*pi*step_size/(xpro+xnu)*hc_mev_fm;
  if (xnu<1.0e-3) xnu=1.0e12;
  if (xpro<1.0e-3) xpro=1.0e12;
  rprms=sqrt(rprms/xpro*4.0*pi*step_size);
  rnrms=sqrt(rnrms/xnu*4.0*pi*step_size);
  rnrp=rnrms-rprms;

  if (verbose>0) {
    cout << "p rms radius=" << rprms << " fm n rms=" 
	 << rnrms << " E/A=" << etot << " MeV " << endl;
  }
    
  return 0;
}

void nucleus_rmf::center_mass_corr(double atot) {
  int nn,i,j;
  double xqmax, hw, b, factor, x;
  double charge=0.0;
  double chcm=0.0;
  
  ubvector fq(100), rhoq(100), xq(100), w(100);

  // Does not work with nn!=100
  nn=100;
  xqmax=8.0;
  // AWS: The harmonic-oscillator energy in inverse fm,
  // see e.g. 2.30 in Negele (1970), but the coefficients
  // here may have been updated by Paul Ellis
  hw=(3.923+23.265/cbrt(atot))/hc_mev_fm;
  b=hc_mev_fm/sqrt(mnuc*hw);
  factor=b*b/(4.0*atot);
  
  // readjust the abscissas and weights
  for(int iw=0;iw<nn;iw++) {
    if (iw<nn/2) {
      w[iw]=w100[49-iw]*xqmax/2.0;
      xq[iw]=(-x100[49-iw]+1.0)/2.0*xqmax;
    } else { 
      w[iw]=w100[iw-50]*xqmax/2.0;
      xq[iw]=(x100[iw-50]+1.0)/2.0*xqmax;
    }
  }
  
  // Correction e^{b^2 q^2/4/A}, see Tassie and Barker (1958)
  for (i=0;i<nn;i++) {
    fq[i]=exp(xq[i]*xq[i]*factor);
  }

  // Transform to momentum space
  for (i=0;i<nn;i++) {
    rhoq[i]=0.0;
    for (j=0;j<grid_size;j++) {
      x=((double)(j+1))*step_size;
      rhoq[i]=rhoq[i]+x*sin(xq[i]*x)/xq[i]*chden1[j];
    }
    rhoq[i]=step_size*rhoq[i]*4.0*pi;
  }

  // Multiply by correction e^{b^2 q^2/4/A} and transform back
  // to position space
  for (i=0;i<grid_size;i++) {
    chdenc[i]=0.0;
    x=((double)(i+1))*step_size;
    for (j=0;j<nn;j++) {
      chdenc[i]=chdenc[i]+w[j]*xq[j]*sin(xq[j]*x)/x*rhoq[j]*fq[j];
    }
    chdenc[i]=chdenc[i]/(2.0*pow(pi,2.0));
  }

  r_charge=0.0;
  r_charge_cm=0.0;

  for (i=0;i<grid_size;i++) {

    x=((double)(i+1))*step_size;
    r_charge=r_charge+pow(x,4.0)*chden1[i];
    r_charge_cm=r_charge_cm+pow(x,4.0)*chdenc[i];
    chcm=chcm+chdenc[i]*x*x;
    charge=charge+chden1[i]*x*x;
    
    if (!std::isfinite(chdenc[i])) chdenc[i]=0.0;
    double line[3]={x,chden1[i],chdenc[i]};
    chden_table->line_of_data(3,line);
    
  }

  r_charge=sqrt(r_charge/charge);
  r_charge_cm=sqrt(r_charge_cm/chcm);

  return;
}

void nucleus_rmf::init_meson_density() {

  for (int i=0;i<grid_size;i++) {
    
    fields(i,0)=0.75*field0(i,0)+0.25*fields(i,0);
    field0(i,0)=fields(i,0);
    fields(i,1)=0.75*field0(i,1)+0.25*fields(i,1);
    field0(i,1)=fields(i,1);
    fields(i,2)=0.75*field0(i,2)+0.25*fields(i,2);
    field0(i,2)=fields(i,2);
    
    xrhosp[i]=0.0;
    
    xrhos[i]=sigma_rhs(fields(i,0),fields(i,1),fields(i,2));
    xrhov[i]=omega_rhs(fields(i,0),fields(i,1),fields(i,2));
    xrhor[i]=rho_rhs(fields(i,0),fields(i,1),fields(i,2));
    
    xrho(i,0)=xrhos[i]*pow(((double)(i+1))*step_size,2.0);
    xrho(i,1)=xrhov[i]*pow(((double)(i+1))*step_size,2.0);
    xrho(i,2)=xrhor[i]*pow(((double)(i+1))*step_size,2.0);
    xrho(i,3)=0.0;
  }

  return;
}

double nucleus_rmf::sigma_rhs(double sig, double ome, double rho) {
  double xm, ret, sig2, gs, dfdphi;

  //--------------------------------------------------------------
  // The fields are in fm^{-1}

  gs=rmf->cs*rmf->ms;

  // sigma in fm^-1
  sig2=sig/gs;
  dfdphi=rmf->a1+2.0*rmf->a2*sig2+3.0*rmf->a3*sig2*sig2+
    4.0*rmf->a4*pow(sig2,3.0)+5.0*rmf->a5*pow(sig2,4.0)+
    6.0*rmf->a6*pow(sig2,5.0);

  // ret should be in MeV^3 here (kappa, sig, and rho are
  // in MeV):
  ret=-rmf->b*rmf->mnuc*sig*sig-rmf->c*pow(sig,3.0)+
    rho*rho/gs*dfdphi;
    
  return ret;
}

double nucleus_rmf::omega_rhs(double sig, double ome, double rho) {
  double ret, gw, omet, dfdome;

  //--------------------------------------------------------------
  // The fields are in fm^{-1}

  gw=rmf->cw*rmf->mw;
  omet=ome/gw;

  dfdome=2.0*rmf->b1*omet+4.0*rmf->b2*pow(omet,3.0)+
    6.0*rmf->b3*pow(omet,5.0);

  // ret should be in MeV^3 here (rho is in MeV):
  ret=-rho*rho/gw*dfdome-pow(ome,3.0)*rmf->zeta/6.0;

  return ret;
}

double nucleus_rmf::rho_rhs(double sig, double ome, double rho) {
  double ret, f, sigt, omet, gs, gw;

  //--------------------------------------------------------------
  // The fields are in fm^{-1}

  gs=rmf->cs*rmf->ms;
  gw=rmf->cw*rmf->mw;
  sigt=sig/gs;
  omet=ome/gw;
  f=rmf->a1*sigt+rmf->a2*sigt*sigt+rmf->a3*pow(sigt,3.0)+
    rmf->a4*pow(sigt,4.0)+rmf->a5*pow(sigt,5.0)+rmf->a6*pow(sigt,6.0)+
    rmf->b1*pow(omet,2.0)+rmf->b2*pow(omet,4.0)+rmf->b3*pow(omet,6.0);

  ret=-2.0*rho*f-pow(rho,3.0)*rmf->xi/6.0;
    
  return ret;
}

int nucleus_rmf::dirac(int ilevel) {
  int iturn, i, j, jmatch, jtop, no=0;
  double deltae, x, yfs, ygs, alpha, xk, v0, e;
  double scale, xnorm=0.0, factor, x1, x2, yf1, yf2, yg1, yg2;
  
  ubvector g(grid_size), v(grid_size), f(grid_size), xz(6);
  
  if ((*levp)[ilevel].energy>0) {
    (*levp)[ilevel].energy=(*levp)[ilevel-1].eigen;
  }      
  (*levp)[ilevel].eigen=(*levp)[ilevel].energy;

  //--------------------------------------------------------------
  // combine vector and rho fields plus photons
    
  for (i=0;i<grid_size;i++) {
    v[i]=(fields(i,1)+(*levp)[ilevel].isospin*
	  (fields(i,2)+fields(i,3))+fields(i,3)/2.0)*hc_mev_fm;
  }

  iturn=1;
  deltae=dirac_tol2*3.0;
  while (iturn<=dirac_itmax && fabs(deltae)>dirac_tol2) {

    //--------------------------------------------------------------
    // Small r solutions

    g[0]=10.0*pow(step_size,-(*levp)[ilevel].kappa);
    f[0]=(v[0]-fields(0,0)*hc_mev_fm-(*levp)[ilevel].eigen)*g[0]*
      step_size/hc_mev_fm;
    f[0]=f[0]/(1.0-2.0*(*levp)[ilevel].kappa);

    if ((*levp)[ilevel].kappa>=0) {
      g[0]=10.0*pow(0.04,1.0+(*levp)[ilevel].kappa);
      f[0]=hc_mev_fm/0.04*g[0]*(2.0*(*levp)[ilevel].kappa+1.0);
      f[0]=f[0]/((*levp)[ilevel].eigen-v[0]-fields(0,0)*hc_mev_fm+
		 2.0*mnuc*hc_mev_fm);
    }

    jmatch=(int)(25.0*(*levp)[ilevel].match_point+0.5+1.0e-6);
    ode_y[0]=f[0];
    ode_y[1]=g[0];
    x=step_size;

    double dtmptmp=(*levp)[ilevel].kappa;
    
    for (i=1;i<jmatch;i++) {
      dirac_step(x,step_size,(*levp)[ilevel].eigen,(*levp)[ilevel].kappa,v);
      f[i]=ode_y[0];
      g[i]=ode_y[1];
    }
    
    //--------------------------------------------------------------
    // store end values for latter matching
  
    yfs=ode_y[0];
    ygs=ode_y[1];
  
    //--------------------------------------------------------------
    // Large r solutions
    // here the code does not use coulomb wavefunctions
    // but simply uses assym. expansions in 1/r which seem
    // to be just fine for bound state coulomb wavefunctions
  
    alpha=sqrt(-(*levp)[ilevel].eigen*((*levp)[ilevel].eigen+
				       2.0*mnuc*hc_mev_fm))/hc_mev_fm;
    g[grid_size-1]=exp(-12.0*alpha);
    xz[1]=-sqrt(-(*levp)[ilevel].eigen/((*levp)[ilevel].eigen+
					2.0*mnuc*hc_mev_fm));
    xk=(*levp)[ilevel].kappa;
    x=12.0;
    v0=12.0*v[grid_size-1]/hc_mev_fm;
    e=hc_mev_fm/2.0/((*levp)[ilevel].eigen+2.0*mnuc*hc_mev_fm);
    
    f[grid_size-1]=xz[1]*g[grid_size-1];
    ode_y[0]=f[grid_size-1];
    ode_y[1]=g[grid_size-1];
    jtop=grid_size-jmatch;

    for (j=1;j<=jtop;j++) {
      dirac_step(x,-step_size,(*levp)[ilevel].eigen,
		 (*levp)[ilevel].kappa,v);
      i=grid_size-j;
      f[i-1]=ode_y[0];
      g[i-1]=ode_y[1];
    }
  
    //--------------------------------------------------------------
    // Match solutions
  
    scale=ode_y[1]/ygs;
  
    jtop=jmatch-1;
    for (i=0;i<jtop;i++) {
      f[i]=f[i]*scale;
      g[i]=g[i]*scale;
    }

    //--------------------------------------------------------------
    // Compute normalization integral and count nodes

    no=0;
    xnorm=pow(f[0],2.0)+pow(g[0],2.0);
    for (i=1;i<grid_size;i++) {
      if (g[i]*g[i-1]<0.0) no++;
      xnorm=xnorm+pow(f[i-1],2.0)+pow(g[i-1],2.0);
    }
    xnorm=xnorm*step_size;

    //--------------------------------------------------------------
    // Adjust eigenvalue

    deltae=-g[jmatch-1]*(f[jmatch-1]-yfs*scale)*hc_mev_fm/xnorm;

    if (iturn==1) deltae=deltae/2.0;

    (*levp)[ilevel].eigen=(*levp)[ilevel].eigen+deltae;

    //--------------------------------------------------------------
    // If the level is unbound, then set to a
    // default value.

    if ((*levp)[ilevel].eigen>0.0) {
      (*levp)[ilevel].eigen=-4.0/((double)(iturn));
    }
    
    iturn++;
  }

  if (iturn>dirac_itmax) {
    if (verbose>0) {
      cout << "No Dirac convergence (" << dirac_itmax << ") tries. " << endl;
      cout << "\tLevel: " << (*levp)[ilevel].state << ", deltae: " 
	   << deltae << endl;
    }
    O2SCL_CONV2_RET("Dirac failed to converge in ",
		    "nucleus_rmf::dirac().",exc_efailed,err_nonconv);
  }

  //--------------------------------------------------------------
  // sum up the densities
  // rho = (2j+1)/norm(f*f+g*g)/(4pi*x*x)
  
  (*levp)[ilevel].nodes=no;

  factor=(*levp)[ilevel].twojp1/xnorm/4.0/pi;
  for (i=0;i<grid_size;i++) {
    
    xrhosp[i]=xrhosp[i]+factor*((*levp)[ilevel].isospin+0.5)*
      (g[i]*g[i]-f[i]*f[i]);
    xrho(i,0)=xrho(i,0)+factor*(g[i]*g[i]-f[i]*f[i]);
    xrho(i,1)=xrho(i,1)+factor*(g[i]*g[i]+f[i]*f[i]);
    xrho(i,2)=xrho(i,2)+factor*(*levp)[ilevel].isospin*
      (g[i]*g[i]+f[i]*f[i]);
    xrho(i,3)=xrho(i,3)+factor*((*levp)[ilevel].isospin+0.5)*
      (g[i]*g[i]+f[i]*f[i]);

  }

  return 0;
}

double nucleus_rmf::dirac_rk4(double x, double g1, double f1, double &funt, 
			 double eigent, double kappat, ubvector &varr) {
  double ret,v,s;
  field(x,s,v,varr);
  ret=f1*(eigent+mnuc*hc_mev_fm*2.0-s-v)/hc_mev_fm-kappat*g1/x;
  funt=g1*(v-s-eigent)/hc_mev_fm+kappat*f1/x;
  return ret;
}

void nucleus_rmf::dirac_step(double &x, double ht, 
			     double eigent, double kappat, ubvector &varr) {

  if (generic_ode) {
    
    odparms op={eigent,kappat,&fields,&varr};
    odefun(x,2,ode_y,ode_dydx,op);

    ode_funct ofm=std::bind
      (std::mem_fn<int(double,size_t,const ubvector &,ubvector &,odparms &)>
       (&nucleus_rmf::odefun),
       this,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,std::placeholders::_4,std::ref(op));

    //ode_funct_mfptr_param<nucleus_rmf,odparms,ubvector> 
    //ofm(this,&nucleus_rmf::odefun,op);
    ostep->step(x,ht,2,ode_y,ode_dydx,ode_y,ode_yerr,ode_dydx,ofm);

    x+=ht;

  } else {
    
    double g1,f1,gg2,f2,g3,f3,g4,f4,funt;
    
    g1=ht*dirac_rk4(x,ode_y[1],ode_y[0],funt,eigent,kappat,varr);
    f1=ht*funt;
    
    gg2=ht*dirac_rk4(x+ht/2.0,ode_y[1]+g1/2.0,ode_y[0]+f1/2.0,funt,
		eigent,kappat,varr);
    f2=ht*funt;
    
    g3=ht*dirac_rk4(x+ht/2.0,ode_y[1]+gg2/2.0,ode_y[0]+f2/2.0,funt,
	       eigent,kappat,varr);
    f3=ht*funt;
    
    g4=ht*dirac_rk4(x+ht,ode_y[1]+g3,ode_y[0]+f3,funt,eigent,kappat,varr);
    f4=ht*funt;
    
    ode_y[0]+=(f1+2.0*(f2+f3)+f4)/6.0;
    ode_y[1]+=(g1+2.0*(gg2+g3)+g4)/6.0;
    x=x+ht;
    
  }

  return;

}

int nucleus_rmf::odefun(double x, size_t nv, const ubvector &y,
			ubvector &dydx, odparms &op) {
  double s, v;
  
  field(x,s,v,*op.varr);
  
  dydx[0]=y[1]*(v-s-op.eigen)/hc_mev_fm+op.kappa*y[0]/x;
  dydx[1]=y[0]*(op.eigen+mnuc*hc_mev_fm*2.0-s-v)/hc_mev_fm-op.kappa*y[1]/x;
  
  return 0;
}

void nucleus_rmf::field(double x, double &s, double &v1, ubvector &v) {

  int i=(int)(x*25.0+0.1+1.0e-6);
  s=fields(i-1,0)*hc_mev_fm;
  v1=v[i-1];
  if (x*25.0-((double)(i))<0.3) return;
  s=(s+fields(i,0)*hc_mev_fm)*0.5;
  v1=(v1+v[i])*0.5;
  return;
}

void nucleus_rmf::meson_init() {

  for (size_t i=0;i<3;i++) {
    double mass;
    if (i==0) mass=rmf->ms;
    else if (i==1) mass=rmf->mw;
    else mass=rmf->mr;
    for (size_t j=0;((int)j)<grid_size;j++) {
      double x=((double)(j+1))*step_size;
      double ex=exp(mass*x);
      gin(j,i)=0.5/mass/x*ex;
      gout(j,i)=1.0/x/ex;
      gin(j,3)=1.0;
      gout(j,3)=1.0/x;
    }
  }
  
  return;
}

void nucleus_rmf::meson_iter(double ic) {
  double xi20, xx, xmin, xmax, x;

  ubvector xi1(grid_size), xi2(grid_size), xf(4);

  xf[0]=3.0/8.0;
  xf[1]=7.0/6.0-3.0/8.0;
  xf[2]=23.0/24.0-7.0/6.0;
  xf[3]=1.0/24.0;
  
  for (int l=0;l<ic;l++) {

    xi1[0]=xf[1]*gin(0,l)*xrho(0,l)+
      xf[2]*gin(1,l)*xrho(1,l)+
      xf[3]*gin(2,l)*xrho(2,l);
    xi2[grid_size-1]=xf[1]*gout(grid_size-1,l)*xrho(grid_size-1,l)+
      xf[2]*gout(grid_size-2,l)*xrho(grid_size-2,l)+
      xf[3]*gout(grid_size-3,l)*xrho(grid_size-3,l);
    
    //--------------------------------------------------------------
    // start doing the inside and outside edges

    for (int i=1;i<3;i++) {
      xi1[i]=xi1[i-1];
      int j=grid_size-i;
      xi2[j-1]=xi2[j];
      for (int k=0;k<4;k++) {
	xi1[i]=xi1[i]+xf[k]*gin(i-1+k,l)*xrho(i-1+k,l);
	xi2[j-1]=xi2[j-1]+xf[k]*gout(j-k,l)*xrho(j-k,l);
      }
    }

    for (int i=3;i<grid_size;i++) {
      int j=grid_size-i;
      xi1[i]=xi1[i-1];
      xi2[j-1]=xi2[j];
      for (int k=0;k<4;k++) {
	xi1[i]=xi1[i]+xf[k]*gin(i-k,l)*xrho(i-k,l);
	xi2[j-1]=xi2[j-1]+xf[k]*gout(j+k-1,l)*xrho(j+k-1,l);
      }
    }
    
    xi20=xi2[0];
    for (int k=1;k<4;k++) {
      xi20=xi20+xf[k]*gout(k-1,l)*xrho(k-1,l);
    }

    double mass;
    if (l==0) mass=rmf->ms;
    else if (l==1) mass=rmf->mw;
    else mass=rmf->mr;

    if (l!=3) xx=-0.5/mass;
    if (l==3) xx=0.0;
    xi20=xi20*xx;
    
    for (int i=0;i<grid_size;i++) {
      fields(i,l)=(gout(i,l)*(xi1[i]+xi20)+gin(i,l)*xi2[i]);
      if (l==0) {
	fields(i,0)=fields(i,0)*step_size*rmf->cs*rmf->ms*rmf->cs*
	  rmf->ms;
      } else if (l==1) {
	fields(i,1)=fields(i,1)*step_size*rmf->cw*rmf->mw*rmf->cw*
	  rmf->mw;
      } else if (l==2) {
	fields(i,2)=fields(i,2)*step_size*rmf->cr*rmf->mr*rmf->cr*
	  rmf->mr;
      } else {
	fields(i,3)=fields(i,3)*step_size*4.0*pi*o2scl_const::fine_structure;
      }
    }

    /*
      ofstream fout("temp.txt");
      fout.precision(4);
      fout << "r rho gin gout xi1 xi2 f" << endl;
      for(int kk=0;kk<grid_size-1;kk++) {
      fout << ((double)(kk+1))*step_size << " " 
      << xrho[kk][l] << " " << gin[kk][l] << " " << gout[kk][l] << " "
      << xi1[kk] << " " << xi2[kk] << " " << fields[kk][l] << endl;
      }
      fout.precision(6);
      exit(-1);
    */

  }

  return;
}

void nucleus_rmf::pfold(double x, double &xrhof) {

  double xmin, xmax, xi, xi1;

  xmin=x-2.7;
  xmax=x;
  if (xmin<0.0) xmin=0.0001;
  gauss(xmin,xmax,x,xi);
  xmax=x+2.7;
  xmin=x;
  gauss(xmin,xmax,x,xi1);
  xrhof=(xi+xi1)*a_proton/4.0/x;

  return;
}

double nucleus_rmf::xpform(double x, double xp, double a) {
  double xra, ret, xr;

  // Need to document 65.0

  /*
    From Paul's notes:
    
    ((1+mu*|x-xp|*e^{-mu*|x-xp|}-(1+mu*(x+xp)*e^{-mu*(x+xp)}))/xp
  */

  xra=x-xp;
  xra=fabs(xra)*a;
  if (fabs(xra)>65.0) xra=65.0*xra/fabs(xra);
  xr=(x+xp)*a;
  if (fabs(xr)>65.0) xr=65.0*xr/fabs(xr);
  ret=((1.0+xra)/exp(xra)-(1.0+xr)/exp(xr))/xp;

  return ret;
}

void nucleus_rmf::gauss(double xmin, double xmax, double x, double &xi) {
  double xdelta, x1, x2;
  
  ubvector xpnew(12), wnew(12);
  for(size_t ii=0;ii<12;ii++) {
    if (ii<6) {
      xpnew[ii]=-x12[5-ii];
      wnew[ii]=w12[5-ii];
    } else {
      xpnew[ii]=x12[ii-6];
      wnew[ii]=w12[ii-6];
    }
   
  }

  ubvector charge(grid_size), ax(grid_size);
  for (size_t i=0;((int)i)<grid_size;i++) {
    ax[i]=((double)(i+1))*step_size;
    charge[i]=xrho(i,3)/(ax[i]*ax[i]);
  }
  gi=new interp_vec<ubvector>(grid_size,ax,charge);
  
  xi=0.0;
  xdelta=xmax-xmin;
  for (size_t i=6;i<12;i++) {
    x1=xmin+xdelta*(0.5+xpnew[i]/2.0);
    x2=xmin+xdelta*(0.5-xpnew[i]/2.0);
    xi=xi+xpform(x,x1,a_proton)*xrhop(x1)*wnew[i]*xdelta/2.0;
    xi=xi+xpform(x,x2,a_proton)*xrhop(x2)*wnew[i]*xdelta/2.0;
  }
  
  delete gi;

  return;
}

double nucleus_rmf::xrhop(double x1) {
  
  double xrhopret;
  
  xrhopret=gi->eval(x1);
  if (x1>12.0) xrhopret=0.0;
  xrhopret=xrhopret*x1*x1;

  return xrhopret;
}

