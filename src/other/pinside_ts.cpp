/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#include <o2scl/pinside.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  pinside pi;
  pi.test(t);

  ubvector x(12), y(12);
  x[0]=1; y[0]=0;
  x[1]=2; y[1]=0;
  x[2]=2; y[2]=1;
  x[3]=3; y[3]=1;
  x[4]=3; y[4]=2;
  x[5]=2; y[5]=2;
  x[6]=2; y[6]=3;
  x[7]=1; y[7]=3;
  x[8]=1; y[8]=2;
  x[9]=0; y[9]=2;
  x[10]=0; y[10]=1;
  x[11]=1; y[11]=1;

  // This generates a file which should contain only the points
  // inside the region defined above (a "plus sign").

  ofstream fout("pinside.dat");
  for(double r=0.1;r<=1.501;r+=0.1) {
    for(double th=0.0;th<=2.0*acos(-1.0);th+=2.0*acos(-1.0)/300.0) {
      if (pi.inside(1.5+r*cos(th),1.5+r*sin(th),x,y)) {
	fout << 1.5+r*cos(th) << " " << 1.5+r*sin(th) << endl;
      }
    }
  }
  fout.close();

  t.report();
  return 0;
}
