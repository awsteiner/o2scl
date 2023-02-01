/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/uniform_grid.h>

typedef boost::numeric::ublas::vector<double> ubvector;

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  uniform_grid_end<> g1(0,1,5);
  ubvector v1(g1.get_npoints());
  g1.vector(v1);

  //cout << v1 << endl;
  cout << "g1: ";
  for(size_t i=0;i<g1.get_npoints();i++) {
    cout << g1[i] << " ";
  }
  cout << endl;

  uniform_grid<> g1b=g1;
  cout << "g1b: ";
  for(size_t i=0;i<g1.get_npoints();i++) {
    cout << g1b[i] << " ";
  }
  cout << endl;

  // Show that this works with floats as well as doubles
  typedef boost::numeric::ublas::vector<float> ubvector_float;
  uniform_grid_end<float> g2(0,1,5);
  ubvector_float v2(g2.get_npoints());
  g2.vector(v2);

  cout << "g2: ";
  for(size_t i=0;i<g2.get_npoints();i++) {
    cout << g2[i] << " ";
  }
  cout << endl;

  // Logarithmic grid
  uniform_grid_log_end<> g3(0.1,1,5);
  ubvector v3(g3.get_npoints());
  g3.vector(v3);

  cout << "g3: ";
  for(size_t i=0;i<g3.get_npoints();i++) {
    cout << g3[i] << " ";
  }
  cout << endl;

  double wid=g3[1]/g3[0];

  // Logarithmic grid specifying end and width
  uniform_grid_log_end_width<> g4(0.1,1,wid);
  ubvector v4(g4.get_npoints());
  g4.vector(v4);

  cout << "g4: ";
  for(size_t i=0;i<g4.get_npoints();i++) {
    cout << g4[i] << " ";
  }
  cout << endl;

  // Linear decreasing grid
  uniform_grid_end<> g5(1,0,5);

  cout << "g5: ";
  for(size_t i=0;i<g5.get_npoints();i++) {
    cout << g5[i] << " ";
  }
  cout << endl;
  
  // Logarthmic decreasing grid
  uniform_grid_log_end<> g6(1,0.1,5);

  cout << "g6: ";
  for(size_t i=0;i<g6.get_npoints();i++) {
    cout << g6[i] << " ";
  }
  cout << endl;

  // Logarthmic grid with negative values
  uniform_grid_log_end<> g7(-0.1,-1.0,5);
  
  cout << "g7: ";
  for(size_t i=0;i<g7.get_npoints() && i<10;i++) {
    cout << g7[i] << " ";
  }
  cout << endl;
  
  // Logarthmic decreasing grid with negative values
  uniform_grid_log_end<> g8(-1.0,-0.1,5);

  cout << "g8: ";
  for(size_t i=0;i<g8.get_npoints() && i<10;i++) {
    cout << g8[i] << " ";
  }
  cout << endl;

  // Show how end_width handles inexact width parameters
  
  uniform_grid_end_width<> g9(0,1,0.2);
  ubvector v9(g9.get_npoints());
  g9.vector(v9);

  //cout << v9 << endl;
  cout << "g9: ";
  for(size_t i=0;i<g9.get_npoints();i++) {
    cout << g9[i] << " ";
  }
  cout << endl;

  uniform_grid_end_width<> ga(0,1,0.2-1.0e-12);
  ubvector va(ga.get_npoints());
  ga.vector(va);

  //cout << va << endl;
  cout << "ga: ";
  for(size_t i=0;i<ga.get_npoints();i++) {
    cout << ga[i] << " ";
  }
  cout << endl;

  uniform_grid_end_width<> gb(0,1,0.2+1.0e-12);
  ubvector vb(gb.get_npoints());
  gb.vector(vb);

  //cout << vb << endl;
  cout << "gb: ";
  for(size_t i=0;i<gb.get_npoints();i++) {
    cout << gb[i] << " ";
  }
  cout << endl;

  t.report();
  return 0;
}
