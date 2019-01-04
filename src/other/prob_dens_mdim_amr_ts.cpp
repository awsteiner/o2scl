/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2019, Andrew W. Steiner
  
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
#include <iostream>

#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/prob_dens_mdim_amr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  std::vector<double> low={1.0,1.0};
  std::vector<double> high={2.0,2.0};

  prob_dens_mdim_amr<std::vector<double>,
		     matrix_view_table<std::vector<double> > > amr(low,high);

  table<> t;
  t.line_of_names("x y z");
  {
    double line[3]={1.8,1.8,1.0};
    t.line_of_data(3,line);
  }
  {
    double line[3]={1.2,1.3,2.0};
    t.line_of_data(3,line);
  }
  {
    double line[3]={1.2,1.31,3.0};
    t.line_of_data(3,line);
  }
  matrix_view_table<std::vector<double> > mvt(t,{"x","y","z"});
  amr.initial_parse(mvt);
  tm.test_rel(amr.total_volume(),1.0,1.0e-8,"total volume 1");
  cout << amr.total_volume() << endl;
  vector<double> v1={1.21,1.31};
  cout << amr.pdf(v1) << endl;

  table<> t2;
  t2.line_of_names("x y z");
  vector<double> v2(3);
  for(size_t i=0;i<10000;i++) {
    amr(v2);
    v2[2]=amr.pdf(v2);
    t2.line_of_data(3,v2);
  }

  if (false) {
    hdf_file hf;
    hf.open_or_create("amr.o2");
    hdf_output(hf,t2,"amr");
    hf.close();
  }

  rng_gsl r;

  static const size_t N=20;
  
  table<> t3;
  t3.line_of_names("x y z");
  for(size_t i=0;i<N;i++) {
    double line[3]={r.random(),r.random(),r.random()};
    t3.line_of_data(3,line);
  }
  
  std::vector<double> low2={0.0,0.0};
  std::vector<double> high2={1.0,1.0};
  
  matrix_view_table<std::vector<double> > mvt2(t3,{"x","y","z"});

  prob_dens_mdim_amr<std::vector<double>,
		     matrix_view_table<std::vector<double> > >
    amr2(low2,high2);
  amr2.verbose=2;
  amr2.initial_parse(mvt2);
  tm.test_rel(amr2.total_volume(),1.0,1.0e-8,"total volume 2");
  cout << amr2.total_volume() << endl;

  if (false) {
    ofstream fout;
    fout.open("temp.scr");
    fout << "o2graph -set xlo 0 -set xhi 1 -set ylo 0 -set yhi 1 \\" << endl;
    for(size_t i=0;i<N;i++) {
      fout << "-point " << t3.get("x",i) << " "
	   << t3.get("y",i) << " marker=x \\" << endl;
    }
    for(size_t i=0;i<N;i++) {
      fout << "-rect " << amr2.mesh[i].low[0] << " "
	   << amr2.mesh[i].low[1] << " "
	   << amr2.mesh[i].high[0] << " "
	   << amr2.mesh[i].high[1] << " \\" << endl;
    }
    fout << "-show" << endl;
    fout.close();
  }
  
  tm.report();
  
  return 0;
}
