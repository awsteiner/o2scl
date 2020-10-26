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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix<std::string> ubmatrix_string;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  size_t nr=10;
  std::vector<std::string> table[5], ctable(nr);
  std::vector<int> align;

  columnify c;

  align.push_back(1);
  align.push_back(2);
  align.push_back(3);
  align.push_back(4);
  align.push_back(6);
  table[0].push_back("left");
  table[1].push_back("right");
  table[2].push_back("lmid");
  table[3].push_back("rmid");
  table[4].push_back("lnum");
  for(size_t i=0;i<nr;i++) {
    for(size_t j=0;j<5;j++) {
      table[j].push_back(dtos(sin((i+j)*100.0)));
    }
  }
  
  c.align(table,5,nr,ctable,align);
  
  for(size_t i=0;i<nr;i++) {
    cout << ctable[i] << endl;
    t.test_gen(ctable[i].size()==69,"row size");
  }
  cout << endl;

  align[4]=5;
  table[4][0]="dp";

  for(size_t i=1;i<nr+1;i++) {
    table[4][i]=dtos(sqrt(2.0)*pow(10.0,((double)(((int)i)-5)*2)),5,true);
  }

  c.table_lines=1;
  c.align(table,5,nr,ctable,align);
  c.table_lines=0;

  size_t test_size=0;
  for(size_t i=0;i<nr;i++) {
    cout << ctable[i] << endl;
    if (i==0) test_size=ctable[i].size();
    else t.test_gen(ctable[i].size()==test_size,"row size");
  }

  ubmatrix om(3,3);
  om(0,0)=-3.0;
  om(0,1)=1.0;
  om(0,2)=4.0;
  om(1,0)=1.0;
  om(1,1)=5.0;
  om(1,2)=9.0;
  om(2,0)=2.0;
  om(2,1)=6.0;
  om(2,2)=53.0;
  matrix_out<ubmatrix>(cout,om);

  ubmatrix ox(3,4);
  ox(0,0)=-3.0;
  ox(0,1)=1.0;
  ox(0,2)=4.0;
  ox(0,3)=1.0;
  ox(1,0)=5.0;
  ox(1,1)=9.0;
  ox(1,2)=2.0;
  ox(1,3)=6.0;
  ox(2,0)=5.0;
  ox(2,1)=3.0;
  ox(2,2)=5.0;
  ox(2,3)=89.0;
  matrix_out<ubmatrix>(cout,ox);
  cout << endl;

  ubmatrix oy(4,3);
  oy(0,0)=-3.0;
  oy(0,1)=1.0;
  oy(0,2)=4.0;
  oy(1,0)=1.0;
  oy(1,1)=5.0;
  oy(1,2)=9.0;
  oy(2,0)=2.0;
  oy(2,1)=6.0;
  oy(2,2)=5.0;
  oy(3,0)=3.0;
  oy(3,1)=5.0;
  oy(3,2)=89.0;
  cout << "Normal: " << endl;
  matrix_out<ubmatrix>(cout,oy);
  cout << endl;

  cout << "Transpose: " << endl;
  matrix_trans_out<ubmatrix>(cout,oy);
  cout << endl;

  t.report();
  return 0;
}
