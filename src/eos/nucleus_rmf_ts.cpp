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
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/nucleus_rmf.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  double eigen1_pb[39], eigenf_pb[39], rprms_pb[40], eoa_pb[40];
  
  // These numbers were obtained from an earlier and more pristine
  // version of the code
  eigen1_pb[1]=-4.090492e+01;
  eigen1_pb[2]=-3.505993e+01;
  eigen1_pb[3]=-3.458594e+01;
  eigen1_pb[4]=-2.826713e+01;
  eigen1_pb[5]=-2.715770e+01;
  eigen1_pb[6]=-2.503351e+01;
  eigen1_pb[7]=-2.078490e+01;
  eigen1_pb[8]=-1.883148e+01;
  eigen1_pb[9]=-1.608946e+01;
  eigen1_pb[10]=-1.538950e+01;
  eigen1_pb[11]=-1.284224e+01;
  eigen1_pb[12]=-9.924422e+00;
  eigen1_pb[13]=-7.327572e+00;
  eigen1_pb[14]=-4.678966e+00;
  eigen1_pb[15]=-6.162912e+00;
  eigen1_pb[16]=-5.519752e+00;
  eigen1_pb[17]=-6.024013e+01;
  eigen1_pb[18]=-5.380996e+01;
  eigen1_pb[19]=-5.342634e+01;
  eigen1_pb[20]=-4.626541e+01;
  eigen1_pb[21]=-4.535558e+01;
  eigen1_pb[22]=-4.285220e+01;
  eigen1_pb[23]=-3.785370e+01;
  eigen1_pb[24]=-3.622119e+01;
  eigen1_pb[25]=-3.273068e+01;
  eigen1_pb[26]=-3.210403e+01;
  eigen1_pb[27]=-2.878172e+01;
  eigen1_pb[28]=-2.626986e+01;
  eigen1_pb[29]=-2.234329e+01;
  eigen1_pb[30]=-1.923667e+01;
  eigen1_pb[31]=-2.119351e+01;
  eigen1_pb[32]=-1.987945e+01;
  eigen1_pb[33]=-1.575232e+01;
  eigen1_pb[34]=-1.202816e+01;
  eigen1_pb[35]=-9.403301e+00;
  eigen1_pb[36]=-1.038858e+01;
  eigen1_pb[37]=-9.091618e+00;
  eigen1_pb[38]=-8.497554e+00;

  eigenf_pb[1]=-4.800858e+01;
  eigenf_pb[2]=-4.242098e+01;
  eigenf_pb[3]=-4.174174e+01;
  eigenf_pb[4]=-3.542398e+01;
  eigenf_pb[5]=-3.386970e+01;
  eigenf_pb[6]=-3.013840e+01;
  eigenf_pb[7]=-2.745251e+01;
  eigenf_pb[8]=-2.469045e+01;
  eigenf_pb[9]=-2.020373e+01;
  eigenf_pb[10]=-1.914679e+01;
  eigenf_pb[11]=-1.884223e+01;
  eigenf_pb[12]=-1.465960e+01;
  eigenf_pb[13]=-1.047102e+01;
  eigenf_pb[14]=-9.852898e+00;
  eigenf_pb[15]=-8.847358e+00;
  eigenf_pb[16]=-7.723940e+00;
  eigenf_pb[17]=-5.910526e+01;
  eigenf_pb[18]=-5.291661e+01;
  eigenf_pb[19]=-5.232676e+01;
  eigenf_pb[20]=-4.546109e+01;
  eigenf_pb[21]=-4.405166e+01;
  eigenf_pb[22]=-4.095795e+01;
  eigenf_pb[23]=-3.711248e+01;
  eigenf_pb[24]=-3.452583e+01;
  eigenf_pb[25]=-3.061302e+01;
  eigenf_pb[26]=-2.955342e+01;
  eigenf_pb[27]=-2.820276e+01;
  eigenf_pb[28]=-2.419604e+01;
  eigenf_pb[29]=-2.062937e+01;
  eigenf_pb[30]=-1.900494e+01;
  eigenf_pb[31]=-1.898367e+01;
  eigenf_pb[32]=-1.823106e+01;
  eigenf_pb[33]=-1.354131e+01;
  eigenf_pb[34]=-1.116584e+01;
  eigenf_pb[35]=-9.745716e+00;
  eigenf_pb[36]=-9.153783e+00;
  eigenf_pb[37]=-8.433499e+00;
  eigenf_pb[38]=-7.657131e+00;

  rprms_pb[1]=5.565969e+00;
  rprms_pb[2]=5.211760e+00;
  rprms_pb[3]=5.218356e+00;
  rprms_pb[4]=5.216396e+00;
  rprms_pb[5]=5.225783e+00;
  rprms_pb[6]=5.241865e+00;
  rprms_pb[7]=5.262153e+00;
  rprms_pb[8]=5.284449e+00;
  rprms_pb[9]=5.307085e+00;
  rprms_pb[10]=5.328938e+00;
  rprms_pb[11]=5.349285e+00;
  rprms_pb[12]=5.367741e+00;
  rprms_pb[13]=5.384144e+00;
  rprms_pb[14]=5.398489e+00;
  rprms_pb[15]=5.410871e+00;
  rprms_pb[16]=5.421439e+00;
  rprms_pb[17]=5.430373e+00;
  rprms_pb[18]=5.437859e+00;
  rprms_pb[19]=5.444082e+00;
  rprms_pb[20]=5.449217e+00;
  rprms_pb[21]=5.453424e+00;
  rprms_pb[22]=5.456850e+00;
  rprms_pb[23]=5.459616e+00;
  rprms_pb[24]=5.461844e+00;
  rprms_pb[25]=5.463618e+00;
  rprms_pb[26]=5.465024e+00;
  rprms_pb[27]=5.466130e+00;
  rprms_pb[28]=5.466996e+00;
  rprms_pb[29]=5.467666e+00;
  rprms_pb[30]=5.468180e+00;
  rprms_pb[31]=5.468574e+00;
  rprms_pb[32]=5.468874e+00;
  rprms_pb[33]=5.469092e+00;
  rprms_pb[34]=5.469253e+00;
  rprms_pb[35]=5.469371e+00;
  rprms_pb[36]=5.469452e+00;
  rprms_pb[37]=5.469512e+00;
  rprms_pb[38]=5.469547e+00;
  rprms_pb[39]=5.469573e+00;

  eoa_pb[1]=-3.116991e+00;
  eoa_pb[2]=-4.150183e+00;
  eoa_pb[3]=-5.829608e+00;
  eoa_pb[4]=-7.134314e+00;
  eoa_pb[5]=-7.984302e+00;
  eoa_pb[6]=-8.475642e+00;
  eoa_pb[7]=-8.701679e+00;
  eoa_pb[8]=-8.750137e+00;
  eoa_pb[9]=-8.691963e+00;
  eoa_pb[10]=-8.579273e+00;
  eoa_pb[11]=-8.446599e+00;
  eoa_pb[12]=-8.314761e+00;
  eoa_pb[13]=-8.194935e+00;
  eoa_pb[14]=-8.092035e+00;
  eoa_pb[15]=-8.007210e+00;
  eoa_pb[16]=-7.939548e+00;
  eoa_pb[17]=-7.887119e+00;
  eoa_pb[18]=-7.847654e+00;
  eoa_pb[19]=-7.818827e+00;
  eoa_pb[20]=-7.798603e+00;
  eoa_pb[21]=-7.785132e+00;
  eoa_pb[22]=-7.776786e+00;
  eoa_pb[23]=-7.772341e+00;
  eoa_pb[24]=-7.770798e+00;
  eoa_pb[25]=-7.771259e+00;
  eoa_pb[26]=-7.773136e+00;
  eoa_pb[27]=-7.775910e+00;
  eoa_pb[28]=-7.779221e+00;
  eoa_pb[29]=-7.782725e+00;
  eoa_pb[30]=-7.786308e+00;
  eoa_pb[31]=-7.789828e+00;
  eoa_pb[32]=-7.793166e+00;
  eoa_pb[33]=-7.796194e+00;
  eoa_pb[34]=-7.798981e+00;
  eoa_pb[35]=-7.801453e+00;
  eoa_pb[36]=-7.803684e+00;
  eoa_pb[37]=-7.805622e+00;
  eoa_pb[38]=-7.807310e+00;
  eoa_pb[39]=-7.808767e+00;

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  nucleus_rmf el;

  eos_had_rmf rmf;
  o2scl_hdf::rmf_load(rmf,"../../data/o2scl/rmfdata/NL4.o2",true);
  el.set_eos(rmf);

  el.init_run(82,126,2,2);

  int iteration=1, iconverged=0;
  int dirac_converged, meson_converged;
  while(iconverged==0) {
    
    el.iterate(82,126,2,2,iconverged,dirac_converged,meson_converged);

    if (iteration==1) {
      for(int i=0;i<el.nlevels;i++) {
	t.test_rel(el.levels[i].eigen,eigen1_pb[i+1],1.0e-4,
		   "eigenvalues at i=1");
      }
    }
    t.test_rel(el.rprms,rprms_pb[iteration],1.0e-4,"proton radii");
    t.test_rel(el.etot,eoa_pb[iteration],1.0e-4,"total energy");

    iteration++;
  }
  cout << "Iterations: " << iteration << endl;
  t.test_gen(iteration==40,"number of iterations");

  for(int i=0;i<el.nlevels;i++) {
    t.test_rel(el.levels[i].eigen,eigenf_pb[i+1],1.0e-4,
	       "final eigenvalues");
  }
  
  el.post_converge(82,126,2,2);

  t.report();

  return 0;
}


