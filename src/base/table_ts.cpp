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
#include <o2scl/table.h>
#include <o2scl/test_mgr.h>

typedef boost::numeric::ublas::vector<double> ubvector;

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  {
    table<> at(20);
    ofstream fout;
    int i;

    // -------------------------------------------------------------
    // Test some simple get and set functions

    at.new_column("col1");
    at.new_column("col2");

    t.test_str(at.get_column_name(0),"col1","get_column_name() (1)");
    t.test_str(at.get_column_name(1),"col2","get_column_name() (2)");

    double line1[2]={1.0,2.5};
    double line2[2]={3.0,4.5};
    double line3[2]={5.0,3.5};
    at.line_of_data(2,line1);
    at.line_of_data(2,line2);
    at.line_of_data(2,line3);
    t.test_rel(at.get("col1",0),1.0,1.0e-14,"line_of_data(new) and get() (1)");
    t.test_rel(at.get("col2",at.get_nlines()-1),3.5,1.0e-14,
	       "line_of_data(new) and get() (2)");
    //t.test_abs(at.get("col1",3),0.0,1.0e-14,"get() (3)");

    at.functions_columns("m1=col1+col2");
    at.functions_columns("m2=m1+col1");
    for(size_t ii=0;ii<at.get_nlines();ii++) {
      cout << at.get("col1",ii) << " ";
      cout << at.get("col2",ii) << " ";
      cout << at.get("m1",ii) << " ";
      cout << at.get("m2",ii) << " ";
      cout << endl;
      t.test_rel(at.get("m1",ii),at.get("col1",ii)+
		 at.get("col2",ii),1.0e-12,"fc1");
      t.test_rel(at.get("m2",ii),at.get("col1",ii)+
		 at.get("m1",ii),1.0e-12,"fc2");
    }

    // -------------------------------------------------------------
    // Test constants

    at.add_constant("hc",197.33);
    at.add_constant("pi",3.14);
    string tnam;
    double tval=0.0;
    for(size_t ii=0;ii<at.get_nconsts();ii++) {
      at.get_constant(ii,tnam,tval);
      cout << ii << " " << tnam << " " << tval << endl;
    }
  
    // -------------------------------------------------------------
    // Test copy constructors

    table<> at4(at);
    table<> at5;
    at5=at;
    t.test_gen(at4.get_nlines()==at.get_nlines(),"cc 1");
    t.test_gen(at4.get_ncolumns()==at.get_ncolumns(),"cc 2");
    t.test_gen(at5.get_nlines()==at.get_nlines(),"cc 3");
    t.test_gen(at5.get_ncolumns()==at.get_ncolumns(),"cc 4");

    // -------------------------------------------------------------
    // Test swap column data and rename

    string cn=at5.get_column_name(0);
    vector<double> v(at5.get_maxlines());
    for(size_t i=0;i<v.size();i++) {
      v[i]=((double)i)+1.0;
    }
    at5.swap_column_data("col1",v);
    t.test_rel(v[0],1.0,1.0e-12,"scd 1");
    t.test_rel(v[1],3.0,1.0e-12,"scd 2");
    t.test_rel(v[2],5.0,1.0e-12,"scd 3");

    at5.rename_column("col1","rcol1");
    t.test_rel(at5.get("rcol1",0),1.0,1.0e-12,"ren 1");
    t.test_rel(at5.get("rcol1",1),2.0,1.0e-12,"ren 2");
    t.test_rel(at5.get("rcol1",2),3.0,1.0e-12,"ren 3");

  }

  {
    table<boost::numeric::ublas::vector<double> > at(20);
    ofstream fout;
    int i;

    // -------------------------------------------------------------
    // Test some simple get and set functions

    at.new_column("col1");
    at.new_column("col2");

    t.test_str(at.get_column_name(0),"col1","get_column_name() (1)");
    t.test_str(at.get_column_name(1),"col2","get_column_name() (2)");

    double line1[2]={1.0,2.5};
    double line2[2]={3.0,4.5};
    double line3[2]={5.0,3.5};
    at.line_of_data(2,line1);
    at.line_of_data(2,line2);
    at.line_of_data(2,line3);
    t.test_rel(at.get("col1",0),1.0,1.0e-14,"line_of_data(new) and get() (1)");
    t.test_rel(at.get("col2",at.get_nlines()-1),3.5,1.0e-14,
	       "line_of_data(new) and get() (2)");
    //t.test_abs(at.get("col1",3),0.0,1.0e-14,"get() (3)");

    at.functions_columns("m1=col1+col2");
    at.functions_columns("m2=m1+col1");
    for(size_t ii=0;ii<at.get_nlines();ii++) {
      cout << at.get("col1",ii) << " ";
      cout << at.get("col2",ii) << " ";
      cout << at.get("m1",ii) << " ";
      cout << at.get("m2",ii) << " ";
      cout << endl;
      t.test_rel(at.get("m1",ii),at.get("col1",ii)+
		 at.get("col2",ii),1.0e-12,"fc1");
      t.test_rel(at.get("m2",ii),at.get("col1",ii)+
		 at.get("m1",ii),1.0e-12,"fc2");
    }

    // -------------------------------------------------------------
    // Test constants

    at.add_constant("hc",197.33);
    at.add_constant("pi",3.14);
    string tnam;
    double tval=0.0;
    for(size_t ii=0;ii<at.get_nconsts();ii++) {
      at.get_constant(ii,tnam,tval);
      cout << ii << " " << tnam << " " << tval << endl;
    }
  
    // -------------------------------------------------------------
    // Test copy constructors

    table<boost::numeric::ublas::vector<double> > at4(at);
    table<boost::numeric::ublas::vector<double> > at5;
    at5=at;
    t.test_gen(at4.get_nlines()==at.get_nlines(),"cc 1");
    t.test_gen(at4.get_ncolumns()==at.get_ncolumns(),"cc 2");
    t.test_gen(at5.get_nlines()==at.get_nlines(),"cc 3");
    t.test_gen(at5.get_ncolumns()==at.get_ncolumns(),"cc 4");

    // -------------------------------------------------------------
    // Test swap column data and rename

    string cn=at5.get_column_name(0);
    boost::numeric::ublas::vector<double> v(at5.get_maxlines());
    for(size_t i=0;i<v.size();i++) {
      v[i]=((double)i)+1.0;
    }
    at5.swap_column_data("col1",v);
    t.test_rel(v[0],1.0,1.0e-12,"scd 1");
    t.test_rel(v[1],3.0,1.0e-12,"scd 2");
    t.test_rel(v[2],5.0,1.0e-12,"scd 3");

    at5.rename_column("col1","rcol1");
    t.test_rel(at5.get("rcol1",0),1.0,1.0e-12,"ren 1");
    t.test_rel(at5.get("rcol1",1),2.0,1.0e-12,"ren 2");
    t.test_rel(at5.get("rcol1",2),3.0,1.0e-12,"ren 3");

  }

  t.report();

  return 0;
}


