/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
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
  terminal ter;

  columnify c;

  if (true) {
    std::vector<std::string> table[6], ctable;
    std::vector<int> align;

    // Test data
    table[0].push_back("left");
    table[1].push_back("right");
    table[2].push_back("lmid");
    table[3].push_back("rmid");
    table[4].push_back("dec. point");
    table[5].push_back("lnum");
    for(size_t i=0;i<nr-1;i++) {
      for(size_t j=0;j<6;j++) {
        if (j==4 && i==0) {
          table[j].push_back("dp");
        } else {
          table[j].push_back(dtos(sin((i+j)*100.0)));
        }
      }
    }
    // Demonstrate all of the various alignments
    align.push_back(columnify::align_left);
    align.push_back(columnify::align_right);
    align.push_back(columnify::align_lmid);
    align.push_back(columnify::align_rmid);
    align.push_back(columnify::align_dp);
    align.push_back(columnify::align_lnum);
    
    // No headers, no lines
    c.align(table,6,nr,ctable,align);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
    
    // No headers, simple lines
    c.table_lines=1;
    
    c.align(table,6,nr,ctable,align);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
    
    // No headers, terminal lines
    c.table_lines=2;
    
    c.align(table,6,nr,ctable,align);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
    
    for(size_t i=0;i<6;i++) {
      table[i].clear();
    }
  }
  
  if (true) {
    
    std::vector<std::string> table[1], ctable;
    std::vector<int> align;

    // Test data
    table[0].push_back("dec. point");
    for(size_t j=0;j<1;j++) {
      for(size_t i=0;i<nr-1;i++) {
        table[j].push_back(dtos(sin((i+j)*100.0)));
      }
    }
    align.push_back(columnify::align_dp);
    
    // Two lines of headers, terminal lines
    c.table_lines=2;
    c.align(table,1,nr,ctable,align,2);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << " " << ter.str_len(ctable[i]) << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
    
  }
  
  if (true) {

    std::vector<std::string> table[2], ctable;
    std::vector<int> align;
    
    // Test data
    table[0].push_back("dec. point");
    table[1].push_back("lnum");
    for(size_t j=0;j<2;j++) {
      for(size_t i=0;i<nr-1;i++) {
        table[j].push_back(dtos(sin((i+j)*100.0)));
      }
    }
    align.push_back(columnify::align_dp);
    align.push_back(columnify::align_lnum);

    // Two lines of headers, terminal lines
    c.table_lines=2;
    c.align(table,2,nr,ctable,align,2);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << " " << ter.str_len(ctable[i]) << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
  }
  
  if (true) {

    c.verbose=1;
    std::vector<std::vector<std::string>> table;
    std::vector<std::string> ctable;
    std::vector<int> align;
    table.push_back({"-1","+2","-100","+100","1234532"});
    table.push_back({"-1","10.0","-1.3e4","1e5","-1.234"});
    table.push_back({"-1.0e4","1.345e2","-1.22e2","1.2345e2","-1.44e3"});
    table.push_back({"this","is","a","test","blah"});
    
    // Test data
    align.push_back(columnify::align_auto);
    align.push_back(columnify::align_auto);
    align.push_back(columnify::align_auto);
    align.push_back(columnify::align_auto);

    // Two lines of headers, terminal lines
    c.table_lines=2;
    c.align(table,table.size(),5,ctable,align);
    
    for(size_t i=0;i<ctable.size();i++) {
      cout << ctable[i] << " " << ter.str_len(ctable[i]) << endl;
      t.test_gen(ter.str_len(ctable[i])==ter.str_len(ctable[0]),
                 "row size");
    }
    cout << endl;
    exit(-1);
  }
  
  if (true) {
    
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
    
  }

  t.report();
  return 0;
}
