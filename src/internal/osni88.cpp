#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <o2scl/table3d.h>
#include <o2scl/constants.h>
#include <o2scl/contour.h>
#include <o2scl/gsl_mroot_hybrids.h>
#include <o2scl/vec_stats.h>
#include <o2scl/poly.h>
#include <o2scl/format_float.h>
#include <o2scl/cern_minimize.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl_ext/table3d_fp.h>
#include <o2scl/convert_units_gnu.h>
#include <o2scl/gsl_inte_qag.h>
#include <o2scl/latswe_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_ext;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  double tab[9][6]={
    {0.0003,50,38,200,146,0.000940},
    {0.001,50,39,460,385,0.00179},
    {0.005,50,39,1140,831,0.00813},
    {0.01,40,38,1215,1115,0.0185},
    {0.02,40,35,1485,1302,0.0448},
    {0.03,40,33,1590,1303,0.0784},
    {0.04,40,31,1610,1261,0.121},
    {0.05,20,30,800,1171,0.175},
    {0.06,20,29,780,1105,0.243}
  };

  table_units t;
  t.line_of_names("rho Z_shell Z A_shell A P_shell");
  t.set_unit("rho","1/fm^3");
  t.set_unit("P_shell","MeV/fm^3");
  
  for(size_t i=0;i<9;i++) {
    double line[6];
    for(size_t j=0;j<6;j++) line[j]=tab[i][j];
    t.line_of_data(6,line);
  }

  hdf_file hf;

  string s="Inner crust from Table 3 of Osni et al. 1988.";

  hf.open("odcgcp88.o2");
  hf.sets("comment",s);
  hdf_output(hf,t,"Table3");
  hf.close();


  return 0;
}

