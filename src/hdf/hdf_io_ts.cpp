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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/hdf_io.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // Test of table HDF I/O
  {
    table<> tab, tab2;
    tab.add_constant("pi",acos(-1.0));
    tab.line_of_names("a b c");
    for(size_t ii=0;ii<10;ii++) {
      double d=((double)ii);
      double line[3]={d,sin(d),cos(d)};
      tab.line_of_data(3,line);
    }

    hdf_file hf;
    hf.open_or_create("table.o2");
    hdf_output(hf,tab,"table_test");
    hf.close();

    hf.open("table.o2");
    hdf_input(hf,tab2,"table_test");
    hf.close();

    t.test_gen(tab.get_nlines()==tab2.get_nlines(),"lines");
    t.test_gen(tab.get_ncolumns()==tab2.get_ncolumns(),"cols");
    t.test_gen(tab.get_nconsts()==tab2.get_nconsts(),"cols");
    t.test_rel(tab.get("b",4),tab2.get("b",4),1.0e-8,"data");
  }

  // Test of table_units I/O
  {
    table_units<> tab, tab2;
    tab.add_constant("pi",acos(-1.0));
    tab.line_of_names("a b c");
    tab.set_unit("a","m");
    tab.set_unit("b","cm");
    tab.set_unit("c","km");
    for(size_t i2=0;i2<10;i2++) {
      double d=((double)i2);
      double line[3]={d,sin(d),cos(d)};
      tab.line_of_data(3,line);
    }

    hdf_file hf;
    hf.open_or_create("table_units.o2");
    hdf_output(hf,tab,"table_test");
    hf.close();

    hf.open("table_units.o2");
    hdf_input(hf,tab2,"table_test");
    hf.close();

    t.test_gen(tab.get_nlines()==tab2.get_nlines(),"lines");
    t.test_gen(tab.get_ncolumns()==tab2.get_ncolumns(),"cols");
    t.test_gen(tab.get_nconsts()==tab2.get_nconsts(),"cols");
    t.test_gen(tab.get_unit("a")==tab2.get_unit("a"),"unit");
  }

  // Tests for vector_spec()
  std::vector<double> v=vector_spec("list:1,2,3,4");
  t.test_gen(v.size()==4,"vector_spec().");

  std::vector<double> v2=vector_spec("text:hdf_io_ts_table.txt:0");
  cout << "v2 ";
  vector_out(cout,v2,true);

  std::vector<double> v3=vector_spec("grid:1,5,1.5,log");
  cout << "v3 ";
  vector_out(cout,v3,true);

  std::vector<double> v3b=vector_spec("grid:1,5,0.99");
  cout << "v3b ";
  vector_out(cout,v3b,true);

  std::vector<double> v4=vector_spec("hdf5:table_units.o2:table_test:c");
  cout << "v4 ";
  vector_out(cout,v4,true);

  std::vector<double> v5=vector_spec("func:5:exp(i)");
  cout << "v5 ";
  vector_out(cout,v5,true);

  // Tests for value_spec()
  double d1, d2, d3, d4, d5;
  value_spec("sin(2.0)",d1,2);
  cout << d1 << endl;
  value_spec("shell:ls -l | awk '{print $5}' | tail -n 1",d2,2);
  cout << d2 << endl;

  table<> tx;
  tx.line_of_names("col");
  for(size_t j=0;j<5;j++) {
    double line[]={((double)j)*2.1};
    tx.line_of_data(1,line);
  }
  
  hdf_file hf;
  hf.open_or_create("hdf_io_value_ts.o2");
  hf.seti("testi",2.0);
  uniform_grid_end<double> ug(1,10,5);
  hdf_output(hf,ug,"ug");
  hdf_output(hf,tx,"tx");
  hf.close();

  value_spec("hdf5:hdf_io_value_ts.o2:testi",d3,2);
  cout << d3 << endl;
  value_spec("hdf5:hdf_io_value_ts.o2:ug:2",d4,2);
  cout << d4 << endl;
  value_spec("hdf5:hdf_io_value_ts.o2:tx:col,3",d5,2);
  cout << d5 << endl;

  // Tests for mult_vector_spec()
  vector<std::vector<double>> vv1;
  mult_vector_spec("func:3:10+i:sin(i+1)*cos(j)",vv1);
  for(size_t i=0;i<vv1.size();i++) {
    cout << "vv1 " << i << " ";
    vector_out(cout,vv1[i],true);
  }

  vector<std::vector<double>> vv2;
  mult_vector_spec("text:hdf_io_ts_table.txt:0,2-3",vv2,3);
  for(size_t i=0;i<vv2.size();i++) {
    cout << "vv2 " << i << " ";
    vector_out(cout,vv2[i],true);
  }
  
  vector<std::vector<double>> vv3;
  mult_vector_spec("hdf5:table_units.o2:table_test:?",vv3);
  for(size_t i=0;i<vv3.size();i++) {
    cout << "vv3 " << i << " ";
    vector_out(cout,vv3[i],true);
  }
  
  t.report();

  return 0;
}
