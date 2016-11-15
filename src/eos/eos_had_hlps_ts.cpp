/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/eos_had_hlps.h>
#include <o2scl/interp.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  eos_had_hlps he;
  fermion n(939.0/hc_mev_fm,2.0), pr(939.0/hc_mev_fm,2.0);
  thermo hb;

  {
    typedef boost::numeric::ublas::vector<double> ubvector;

    double K=220.0/hc_mev_fm;
    
    ubvector vcomp, valpha, vgamma, veta;
    o2scl::interp_vec<> iv_alpha;
    o2scl::interp_vec<> iv_gamma;
    o2scl::interp_vec<> iv_eta;
    
    // Compressibility, alpha, gamma, eta
    double arr[21][4]={
      {1.114901e+00,7.018317e+00,1.256284e+00,4.949759e+00},
      {1.125036e+00,6.842056e+00,1.265751e+00,4.773499e+00},
      {1.135172e+00,6.677136e+00,1.275261e+00,4.608578e+00},
      {1.145307e+00,6.523232e+00,1.284771e+00,4.454674e+00},
      {1.155443e+00,6.379277e+00,1.294281e+00,4.310720e+00},
      {1.165578e+00,6.244337e+00,1.303790e+00,4.175780e+00},
      {1.175714e+00,6.117591e+00,1.313300e+00,4.049033e+00},
      {1.185849e+00,5.998313e+00,1.322809e+00,3.929755e+00},
      {1.195985e+00,5.885851e+00,1.332320e+00,3.817294e+00},
      {1.206120e+00,5.779659e+00,1.341829e+00,3.711102e+00},
      {1.216255e+00,5.679206e+00,1.351339e+00,3.610649e+00},
      {1.226391e+00,5.584049e+00,1.360849e+00,3.515491e+00},
      {1.236526e+00,5.493788e+00,1.370358e+00,3.425231e+00},
      {1.246662e+00,5.408039e+00,1.379868e+00,3.339481e+00},
      {1.256797e+00,5.326478e+00,1.389378e+00,3.257920e+00},
      {1.266933e+00,5.248807e+00,1.398887e+00,3.180249e+00},
      {1.277068e+00,5.174753e+00,1.408397e+00,3.106195e+00},
      {1.287204e+00,5.104070e+00,1.417907e+00,3.035512e+00},
      {1.297339e+00,5.036533e+00,1.427416e+00,2.967975e+00},
      {1.307475e+00,4.971936e+00,1.436926e+00,2.903378e+00},
      {1.317610e+00,4.910091e+00,1.446435e+00,2.841533e+00}};
    vcomp.resize(21);
    valpha.resize(21);
    vgamma.resize(21);
    veta.resize(21);
    for(size_t i=0;i<21;i++) {
      vcomp[i]=arr[i][0];
      valpha[i]=arr[i][1];
      vgamma[i]=arr[i][2];
      veta[i]=arr[i][3];
    }
    iv_alpha.set(21,vcomp,valpha,itp_cspline);
    iv_gamma.set(21,vcomp,vgamma,itp_cspline);
    iv_eta.set(21,vcomp,veta,itp_cspline);

    // Set up couplings from parameters
    he.alpha=iv_alpha.eval(K);
    he.gamma=iv_gamma.eval(K);
    he.eta=iv_eta.eval(K);
  }

  double S=32.0/hc_mev_fm;
  double L=50.0/hc_mev_fm;
  
  double C=pow(1.5*0.16*pi2,2.0/3.0)/4.0/939.0*hc_mev_fm;
  double den=18.0*C*(he.gamma-1.0);
  he.alphaL=(C*(-4.0+9.0*he.alpha*(he.gamma-1.0)+6.0*he.gamma)+
             3.0*(L-3.0*S*he.gamma))/den;
  he.etaL=(3*(L-3.0*S)+C*(2.0+9.0*(he.gamma-1.0)*he.eta))/den;

  t.test_rel(he.fesym_slope(0.16)*hc_mev_fm,50.0,1.0e-4,"L");
  cout << he.alphaL << " " << he.etaL << endl;
  he.fix_SL(939.0/hc_mev_fm,0.16,S,L);
  cout << he.alphaL << " " << he.etaL << endl;
	    
  t.report();
  return 0;
}

