/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2019, Andrew W. Steiner

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
  cout.precision(10);

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

  t.report();

  return 0;
}
