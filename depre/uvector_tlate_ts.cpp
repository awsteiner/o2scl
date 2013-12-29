/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/uvector_tlate.h>
#include <o2scl/test_mgr.h>
#include <o2scl/cx_arith.h>
#include <o2scl/vec_arith.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  uvector u(3), u2(3); 
  u[0]=3.0;
  u[1]=1.0;
  u[2]=4.0;

  u2[0]=3.0;
  u2[1]=1.0;
  u2[2]=4.0;

  // Test uvector_view functions
  { 
    // Copy-constructors
    uvector_view av(u);
    uvector_view aw=u;
    
    // Indexing operators
    const double d1=av[0];
    double d2=av[1];
    const double d3=av(2);
    cout << d1 << " " << d2 << " " << d3 << endl;
    t.test_rel(d1,3.0,1.0e-10,"t1a");
    t.test_rel(d2,1.0,1.0e-10,"t2a");
    t.test_rel(d3,4.0,1.0e-10,"t3a");

    // Get and set
    av.set(1,2.0);
    cout << av[0] << " " << av.get(1)-1 << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av.get(1)-1,1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");
    av.set(1,1.0);
    cout << av[0] << " " << av.get(1) << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av.get(1),1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");

    // get size and stride
    cout << av.size() << " " << 1.0 << " " << av[2] << endl;
    t.test_rel(av.size(),3.0,1.0e-10,"t1a");
    t.test_rel(1.0,1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");

    // lookup, max and min
    cout << av.get(av.lookup(3.1)) << " " << av.min() << " " 
	 << av.max() << endl;
    t.test_rel(av.get(av.lookup(3.1)),3.0,1.0e-10,"t1a");
    t.test_rel(av.min(),1.0,1.0e-10,"t2a");
    t.test_rel(av.max(),4.0,1.0e-10,"t3a");

    // ostream operator
    cout << av << endl;
    cout << aw << endl;

    // operator *=
    av*=2;
    cout << av[0]/2 << " " << av[1]/2 << " " << av[2]/2 << endl;
    t.test_rel(av[0]/2,3.0,1.0e-10,"t1a");
    t.test_rel(av[1]/2,1.0,1.0e-10,"t2a");
    t.test_rel(av[2]/2,4.0,1.0e-10,"t3a");

    // operator +=, -=
    av+=u2;
    av-=u2;
    av-=u2;
    cout << av[0] << " " << av[1] << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av[1],1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");
  }

  // Test uvector functions
  {
    // Copy-constructors
    uvector av(u);
    uvector aw=u;

    // Indexing operators
    const double d1=av[0];
    double d2=av[1];
    const double d3=av(2);
    cout << d1 << " " << d2 << " " << d3 << endl;
    t.test_rel(d1,3.0,1.0e-10,"t1a");
    t.test_rel(d2,1.0,1.0e-10,"t2a");
    t.test_rel(d3,4.0,1.0e-10,"t3a");

    // Get and set
    av.set(1,2.0);
    cout << av[0] << " " << av.get(1)-1 << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av.get(1)-1,1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");
    av.set(1,1.0);
    cout << av[0] << " " << av.get(1) << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av.get(1),1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");

    // get size and stride
    cout << av.size() << " " << 1.0 << " " << av[2] << endl;
    t.test_rel(av.size(),3.0,1.0e-10,"t1a");
    t.test_rel(1.0,1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");

    // lookup, max and min
    cout << av.get(av.lookup(3.1)) << " " << av.min() << " " 
	 << av.max() << endl;
    t.test_rel(av.get(av.lookup(3.1)),3.0,1.0e-10,"t1a");
    t.test_rel(av.min(),1.0,1.0e-10,"t2a");
    t.test_rel(av.max(),4.0,1.0e-10,"t3a");

    // ostream operator
    cout << av << endl;
    cout << aw << endl;

    // operator *=
    av*=2;
    cout << av[0]/2 << " " << av[1]/2 << " " << av[2]/2 << endl;
    t.test_rel(av[0]/2,3.0,1.0e-10,"t1a");
    t.test_rel(av[1]/2,1.0,1.0e-10,"t2a");
    t.test_rel(av[2]/2,4.0,1.0e-10,"t3a");

    // operator +=, -=
    av+=u2;
    av-=u2;
    av-=u2;
    cout << av[0] << " " << av[1] << " " << av[2] << endl;
    t.test_rel(av[0],3.0,1.0e-10,"t1a");
    t.test_rel(av[1],1.0,1.0e-10,"t2a");
    t.test_rel(av[2],4.0,1.0e-10,"t3a");
    
    // operator +, -, *
    uvector ax=av-aw;
    cout << ax[0]+av[0] << " " << ax[1]+av[1] << " " 
	 << ax[2]+av[2] << endl;
    t.test_rel(ax[0]+av[0],3.0,1.0e-10,"t1a");
    t.test_rel(ax[1]+av[1],1.0,1.0e-10,"t2a");
    t.test_rel(ax[2]+av[2],4.0,1.0e-10,"t3a");
    uvector ay=ax+av;
    cout << ay[0] << " " << ay[1] << " " << ay[2] << endl;
    t.test_rel(ay[0],3.0,1.0e-10,"t1a");
    t.test_rel(ay[1],1.0,1.0e-10,"t2a");
    t.test_rel(ay[2],4.0,1.0e-10,"t3a");

  }

  // Compare versions of operator=:
  uvector o1(3), o2(3), o3;

  o1[0]=3.0;
  o1[1]=1.0;
  o1[2]=4.0;
  o2[0]=1.0;
  o2[1]=5.0;
  o2[2]=9.0;

  uvector_view ov1(o1), ov2(o2), ov3(o1);

  // vector=vector
  o3=o1;
  cout << o3[0] << " " << o3[1] << " " << o3[2] << endl;
  t.test_rel(o3[0],3.0,1.0e-10,"operator=1");
  t.test_rel(o3[1],1.0,1.0e-10,"operator=2");
  t.test_rel(o3[2],4.0,1.0e-10,"operator=3");
  // vector_view=vector
  ov2=o1;
  cout << ov2[0] << " " << ov2[1] << " " << ov2[2] << endl;
  t.test_rel(ov2[0],3.0,1.0e-10,"operator=4");
  t.test_rel(ov2[1],1.0,1.0e-10,"operator=5");
  t.test_rel(ov2[2],4.0,1.0e-10,"operator=6");
  //vector_view=vector_view
  ov1=ov3;
  cout << ov1[0] << " " << ov1[1] << " " << ov1[2] << endl;
  t.test_rel(ov1[0],3.0,1.0e-10,"operator=7");
  t.test_rel(ov1[1],1.0,1.0e-10,"operator=8");
  t.test_rel(ov1[2],4.0,1.0e-10,"operator=9");

  // Test sort_unique()
  uvector xo(10);
  for(size_t i=0;i<10;i++) xo[i]=i%4;
  xo.sort_unique();
  t.test_gen(xo.size()==4,"sort unique 0");
  t.test_rel(xo[0],0.0,1.0e-12,"sort unique 1");
  t.test_rel(xo[3],3.0,1.0e-12,"sort unique 1");

  // Test uvector_array and operator==
  {
    double arr[3]={3.0,1.0,4.0};
    uvector_array uarr(3,arr);
    uvector uv(3);
    uv[0]=3.0;
    uv[1]=1.0;
    uv[2]=4.0;
    t.test_gen(uv==uarr,"uvector_array and operator==");
  }

  // Test uvector_subvector and operator==
  {
    uvector uv6(6);
    uv6[0]=3.0;
    uv6[1]=1.0;
    uv6[2]=4.0;
    uv6[3]=1.0;
    uv6[4]=5.0;
    uv6[5]=9.0;
    uvector_subvector usub(uv6,2,4);
    uvector uv(4);
    uv[0]=4.0;
    uv[1]=1.0;
    uv[2]=5.0;
    uv[3]=9.0;
    t.test_gen(uv==usub,"uvector_subvector and operator==");
  }

  // Test uvector_const_array and operator==
  {
    double arr[3]={3.0,1.0,4.0};
    const double *carr=arr;
    uvector_const_array uarr(3,carr);
    uvector uv(3);
    uv[0]=3.0;
    uv[1]=1.0;
    uv[2]=4.0;
    t.test_gen(uv==uarr,"uvector_array and operator==");
  }

  // Test uvector_const_subvector and operator==
  {
    uvector uv6(6);
    uv6[0]=3.0;
    uv6[1]=1.0;
    uv6[2]=4.0;
    uv6[3]=1.0;
    uv6[4]=5.0;
    uv6[5]=9.0;
    const uvector &cu=uv6;
    uvector_const_subvector usub(cu,2,4);
    uvector uv(4);
    uv[0]=4.0;
    uv[1]=1.0;
    uv[2]=5.0;
    uv[3]=9.0;
    t.test_gen(uv==usub,"uvector_const_subvector and operator==");
  }

  t.report();
  return 0;
}
