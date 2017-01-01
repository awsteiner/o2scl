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

#include <o2scl/test_mgr.h>
#include <o2scl/bin_size.h>

using namespace std;
using namespace o2scl;

int main (void) {
  test_mgr t;
  t.set_output_level(1);

  double a1, a2, bl, bh, bwid;
  int naa, nb;
  bin_size bs;

  a1=0.0143;
  a2=0.9931;
  naa=10;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.0144458;
  naa=10;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  naa=100;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0001;
  a2=0.9931;
  naa=100;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  naa=2;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  naa=3;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  bwid=0.2;
  naa=1;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  bwid=0.5;
  naa=1;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  a1=0.0143;
  a2=0.9931;
  bwid=1.0;
  naa=1;
  bs.calc_bin(a1,a2,naa,bl,bh,nb,bwid);
  cout << bl << " " << bh << " " << nb << " " << bwid << endl;

  t.report();
  return 0;
}
