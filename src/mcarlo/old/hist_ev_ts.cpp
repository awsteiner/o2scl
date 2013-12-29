/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2013, Andrew W. Steiner
  
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
#include <o2scl/rng_gsl.h>
#include <o2scl/hist_ev.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

#ifdef O2SCL_NEVER_DEFINED

  rng_gsl gr;

  ubvector avg, sd, avge, reps;
  ubvector_int m_block, m_per_block;
  hist h;
  hist_ev he(uniform_grid_end<>(0.0,1.0,10),5,5);
  for(size_t i=0;i<10000;i++) {
    double x=gr.random();
    he.add(x,x);
    he.current_avg_stats(reps,avg,sd,avge,m_block,m_per_block);
    if (false) {
      cout << "Adding: " << x << endl;
      cout.precision(1);
      cout << "               ";
      for(size_t j=0;j<reps.size();j++) {
	cout.width(2);
	cout << m_block[j] << " ";
      }
      cout << endl;
      cout << "               ";
      for(size_t j=0;j<reps.size();j++) {
	cout.width(2);
	cout << m_per_block[j] << " ";
      }
      cout << endl;
      cout << "                ";
      for(size_t j=0;j<reps.size();j++) {
	cout << avg[j] << " ";
      }
      cout << endl;
      cout << "                ";
      for(size_t j=0;j<reps.size();j++) {
	cout << sd[j] << " ";
      }
      cout << endl;
      cout.precision(6);
      char ch;
      cin >> ch;
    }
  }

  size_t tot=0;
  for(size_t i=0;i<10;i++) {
    cout << i << " " << avg[i] << " " << sd[i] << " "
	 << avge[i] << " " << fabs(reps[i]-avg[i])/avge[i] << " "
	 << m_block[i] << " " << m_per_block[i] << endl;
    t.test_gen(fabs(reps[i]-avg[i])/avge[i]<5.0,"avg");
    tot+=m_block[i]*m_per_block[i];
  }
  cout << tot << endl;
  
  //hist_2d_ev he2(h.grid_end_size(0.0,1.0,0.1),
  //h.grid_end_size(0.0,1.0,0.1),5);
  
#endif

  t.report();
  return 0;
}
