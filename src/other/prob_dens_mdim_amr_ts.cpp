/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018, Andrew W. Steiner
  
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
#include <o2scl/prob_dens_mdim_amr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  std::vector<double> low={1.0,1.0};
  std::vector<double> high={2.0,2.0};

  cout << "H1." << endl;
  prob_dens_mdim_amr<std::vector<double>,
		     matrix_view_table<std::vector<double> > > amr(low,high);

  cout << "H2." << endl;
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
  cout << "H3." << endl;
  matrix_view_table<std::vector<double> > mvt(t,{"x","y","z"});
  cout << "H4." << endl;
  amr.initial_parse(mvt);
  cout << "H5." << endl;
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

  hdf_file hf;
  hf.open_or_create("amr.o2");
  hdf_output(hf,t2,"amr");
  hf.close();

  return 0;
}
