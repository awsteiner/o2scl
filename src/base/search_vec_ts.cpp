/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gsl/gsl_spline.h>
#include <o2scl/search_vec.h>
#include <o2scl/test_mgr.h>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  size_t i;

  ubvector x(10);
  double xa[10];

  for(i=0;i<10;i++) {
    x[i]=((double)(i));
  }
  search_vec<ubvector> ib(10,x);

  t.test_gen(ib.ordered_lookup(-5.5)==0,"ordered_lookup (inc)");
  t.test_gen(ib.ordered_lookup(15.5)==9,"ordered_lookup (inc)");
  t.test_gen(ib.ordered_lookup(5.0)==5,"ordered_lookup (inc)");
  t.test_gen(ib.ordered_lookup(5.49)==5,"ordered_lookup (inc)");
  t.test_gen(ib.ordered_lookup(5.51)==6,"ordered_lookup (inc)");

  for(i=0;i<10;i++) {
    x[i]=9.0-((double)(i));
  }

  t.test_gen(ib.ordered_lookup(-5.5)==9,"ordered_lookup (dec)");
  t.test_gen(ib.ordered_lookup(15.5)==0,"ordered_lookup (dec)");
  t.test_gen(ib.ordered_lookup(5.0)==4,"ordered_lookup (dec)");
  t.test_gen(ib.ordered_lookup(5.5)==3,"ordered_lookup (dec)");
  t.test_gen(ib.ordered_lookup(5.51)==3,"ordered_lookup (dec)");
  cout << endl;
  
  cout.setf(ios::scientific);
  cout << "Index Data " << endl;
  for(i=0;i<10;i++) {
    x[i]=((double)(i));
    cout << i << "     " << x[i] << endl;
    xa[i]=x[i];
  }
  cout << endl;

  double idata[11]={-1.0,0.0,0.9999,1.0,1.5,5.0,5.5,8.0,8.8,9.0,9.5};

  cout << "find(),ordered_lookup(),"
       << "bsearch_inc()" << endl;
  cout.setf(ios::showpos);
  for(size_t k=0;k<11;k++) {
    cout << idata[k] << " " << ib.find(idata[k]) << " " 
	 << ib.ordered_lookup(idata[k]) << " " 
	 << vector_bsearch_inc(idata[k],x,0,9) << endl;
  }
  cout.unsetf(ios::showpos);
  cout << endl;

  cout << "Index Data " << endl;
  for(i=0;i<10;i++) {
    x[i]=9.0-((double)(i));
    cout << i << "     " << x[i] << endl;
    xa[i]=x[i];
  }
  cout << endl;

  double idata2[11]={-1.0,0.0,0.9999,1.0,1.5,5.0,5.5,8.0,8.5,9.0,9.5};
  
  cout << "find(), ordered_lookup(), "
       << "bsearch_dec()." 
       << endl;
  cout.setf(ios::showpos);
  for(size_t k=0;k<11;k++) {
    cout << idata2[k] << " " << ib.find_dec(idata2[k]) << " " 
	 << ib.ordered_lookup(idata2[k]) << " " 
	 << vector_bsearch_dec(idata2[k],x,0,9) << endl; 
  }
  cout.unsetf(ios::showpos);
  cout << endl;

  t.report();
  return 0;
}

