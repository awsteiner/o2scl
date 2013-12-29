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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/vec_arith.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/ovector_rev_tlate.h>

using namespace std;
using namespace o2scl;

template<class vec_t, class alloc_t> class vtest {
public:
  int fun(int t=0) {
    vec_t v;
    alloc_t a;
    a.allocate(v,3);
    v[0]=3.0;
    v[1]=1.0;
    v[2]=4.0;
    cout << v[0] << " " << v[1] << " " << v[2] << endl;
    a.free(v);
    return 0;
  }
    
};

int show_data(std::string s, ovector_const_view &v) {
  gsl_vector *g=&v;
  size_t i1=g->size;
  size_t i2=0, i3=0;
  size_t i4=g->owner;
  if (g->block!=0) {
    i2=1;
    i3=g->block->size;
  }
  cout.width(30);
  cout << s << " ";
  cout.width(4);
  cout << i1 << " ";
  cout.width(5);
  cout << i4 << " ";
  cout.width(5);
  cout << i2 << " ";
  cout.width(10);
  cout << i3 << " ";
  cout << v << endl;
  return 0;
}

int main(void) {
  
  // ---------------------------------------------------------
  // Demonstrate how memory allocation works

  vtest<double[3],array_alloc<double[3]> > v1;
  v1.fun();
  cout << "HereA." << endl;
  vtest<double *,pointer_alloc<double> > v2;
  cout << "HereB." << endl;
  v2.fun();
  cout << "HereC." << endl;
  vtest<ovector,ovector_alloc> v3;
  cout << "HereD." << endl;
  v3.fun();

  // ---------------------------------------------------------

  test_mgr t;
  t.set_output_level(1);

  {
    ovector a(3), a2(3); 

    a[0]=3.0;
    a[1]=1.0;
    a[2]=4.0;

    a2[0]=3.0;
    a2[1]=1.0;
    a2[2]=4.0;

    // Test ovector_view functions
    { 
      // Copy-constructors
      ovector_base av(a);
      ovector_base aw=a;
    
      // Indexing operators
      cout << "Here1." << endl;
      const double d1=av[0];
      double d2=av[1];
      const double d3=av(2);
      cout << "Here2." << endl;
      cout << d1 << " " << d2 << " " << d3 << endl;
      t.test_rel(d1,3.0,1.0e-10,"t1a");
      t.test_rel(d2,1.0,1.0e-10,"t2a");
      t.test_rel(d3,4.0,1.0e-10,"t3a");

      // Get and set
      av.set(1,2.0);
      cout << "Here3." << endl;
      cout << av[0] << " " << av.get(1)-1 << " " << av[2] << endl;
      t.test_rel(av[0],3.0,1.0e-10,"t1a");
      t.test_rel(av.get(1)-1,1.0,1.0e-10,"t2a");
      t.test_rel(av[2],4.0,1.0e-10,"t3a");
      cout << "Here4." << endl;
      av.set(1,1.0);
      cout << av[0] << " " << av.get(1) << " " << av[2] << endl;
      t.test_rel(av[0],3.0,1.0e-10,"t1a");
      t.test_rel(av.get(1),1.0,1.0e-10,"t2a");
      t.test_rel(av[2],4.0,1.0e-10,"t3a");
      cout << "Here." << endl;

      // get size and stride
      cout << av.size() << " " << av.stride() << " " << av[2] << endl;
      t.test_rel(av.size(),3.0,1.0e-10,"t1a");
      t.test_rel(av.stride(),1.0,1.0e-10,"t2a");
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
      av+=a2;
      av-=a2;
      av-=a2;
      cout << av[0] << " " << av[1] << " " << av[2] << endl;
      t.test_rel(av[0],3.0,1.0e-10,"t1a");
      t.test_rel(av[1],1.0,1.0e-10,"t2a");
      t.test_rel(av[2],4.0,1.0e-10,"t3a");
    }

    // Test ovector functions
    {
      // Copy-constructors
      ovector av(a);
      ovector aw=a;

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
      cout << av.size() << " " << av.stride() << " " << av[2] << endl;
      t.test_rel(av.size(),3.0,1.0e-10,"t1a");
      t.test_rel(av.stride(),1.0,1.0e-10,"t2a");
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
      av+=a2;
      av-=a2;
      av-=a2;
      cout << av[0] << " " << av[1] << " " << av[2] << endl;
      t.test_rel(av[0],3.0,1.0e-10,"t1a");
      t.test_rel(av[1],1.0,1.0e-10,"t2a");
      t.test_rel(av[2],4.0,1.0e-10,"t3a");
    
      // operator +, -, * (from vec_arith.h)
      ovector ax=av-aw;
      cout << ax[0]+av[0] << " " << ax[1]+av[1] << " " 
	   << ax[2]+av[2] << endl;
      t.test_rel(ax[0]+av[0],3.0,1.0e-10,"t1a");
      t.test_rel(ax[1]+av[1],1.0,1.0e-10,"t2a");
      t.test_rel(ax[2]+av[2],4.0,1.0e-10,"t3a");
      ovector ay=ax+av;
      cout << ay[0] << " " << ay[1] << " " << ay[2] << endl;
      t.test_rel(ay[0],3.0,1.0e-10,"t1a");
      t.test_rel(ay[1],1.0,1.0e-10,"t2a");
      t.test_rel(ay[2],4.0,1.0e-10,"t3a");
      ovector az=ay*2.0;
      ovector aa=az-ay;
      cout << aa[0] << " " << aa[1] << " " << aa[2] << endl;
      t.test_rel(aa[0],3.0,1.0e-10,"t1a");
      t.test_rel(aa[1],1.0,1.0e-10,"t2a");
      t.test_rel(aa[2],4.0,1.0e-10,"t3a");

    }

    // test using gsl_vector_get derived from a 'ovector':
    gsl_vector *gv=(gsl_vector *)(&a);
  
    t.test_rel(gsl_vector_get(gv,0),3.0,1.0e-12,"get2a");
    t.test_rel(gsl_vector_get(gv,1),1.0,1.0e-12,"get2b");
    t.test_rel(gsl_vector_get(gv,2),4.0,1.0e-12,"get2c");

    gsl_vector_set(gv,2,4.1);

    // Show that changes are reflected in the original
    t.test_rel(a[0],3.0,1.0e-12,"get3a");
    t.test_rel(a[1],1.0,1.0e-12,"get3b");
    t.test_rel(a[2],4.1,1.0e-12,"get3c");

    // Test using a gsl_vector as a 'ovector'
    gsl_vector *gv2=gsl_vector_alloc(3);
    gsl_vector_set(gv2,0,3.0);
    gsl_vector_set(gv2,1,1.0);
    gsl_vector_set(gv2,2,4.0);

    ovector *ap=(ovector *)gv2;
    t.test_rel((*ap)[0],3.0,1.0e-12,"get4a");
    t.test_rel((*ap)[1],1.0,1.0e-12,"get4b");
    t.test_rel((*ap)[2],4.0,1.0e-12,"get4c");
  
    gsl_vector_free(gv2);

    /// Test ovector_array and ovector_array_stride
    double arr[3]={3.0,1.0,4.0};
    double arrs[6]={3.0,1.0,4.0,1.0,5.0,9.0};

    ovector_array ava(3,arr);
    t.test_rel(ava[0],3.0,1.0e-12,"get5a");
    t.test_rel(ava[1],1.0,1.0e-12,"get5b");
    t.test_rel(ava[2],4.0,1.0e-12,"get5c");

    ovector_array_stride avas(3,arrs,2);
    t.test_rel(avas[0],3.0,1.0e-12,"get6a");
    t.test_rel(avas[1],4.0,1.0e-12,"get6b");
    t.test_rel(avas[2],5.0,1.0e-12,"get6c");

    // Test lookup
    ovector x(10);

    for(size_t i=0;i<10;i++) {
      x[i]=((double)(i));
    }
    t.test_gen(x.lookup(-5.5)==0,"x.lookup (inc)");
    t.test_gen(x.lookup(15.5)==9,"x.lookup (inc)");
    t.test_gen(x.lookup(5.0)==5,"x.lookup (inc)");
    t.test_gen(x.lookup(5.5)==5,"x.lookup (inc)");
    t.test_gen(x.lookup(5.51)==6,"x.lookup (inc)");

    // Test reverse
    ovector_reverse rx(x);
    for(size_t i=0;i<10;i++) {
      t.test_rel(rx[i],((double)(9-i)),1.0e-10,"reverse");
    }
  }

  {
    // Compare versions of operator=:
    ovector o1(3), o2(3), o3;

    o1[0]=3.0;
    o1[1]=1.0;
    o1[2]=4.0;
    o2[0]=1.0;
    o2[1]=5.0;
    o2[2]=9.0;
  
    ovector_base ov1(o1), ov2(o2), ov3(o1);
  
    // vector=vector (Deep)
    o3=o1;
    cout << o3[0] << " " << o3[1] << " " << o3[2] << endl;
    t.test_rel(o3[0],3.0,1.0e-10,"operator=1");
    t.test_rel(o3[1],1.0,1.0e-10,"operator=2");
    t.test_rel(o3[2],4.0,1.0e-10,"operator=3");
    // vector_view=vector (Shallow)
    ov2=o1;
    cout << ov2[0] << " " << ov2[1] << " " << ov2[2] << endl;
    t.test_rel(ov2[0],3.0,1.0e-10,"operator=4");
    t.test_rel(ov2[1],1.0,1.0e-10,"operator=5");
    t.test_rel(ov2[2],4.0,1.0e-10,"operator=6");
    //vector_view=vector_view (Shallow)
    ov1=ov3;
    cout << ov1[0] << " " << ov1[1] << " " << ov1[2] << endl;
    t.test_rel(ov1[0],3.0,1.0e-10,"operator=7");
    t.test_rel(ov1[1],1.0,1.0e-10,"operator=8");
    t.test_rel(ov1[2],4.0,1.0e-10,"operator=9");
    //vector=vector_view (Deep)
    ovector onew;
    onew=ov1;
    t.test_rel(onew[0],3.0,1.0e-10,"operator=10");
    t.test_rel(onew[1],1.0,1.0e-10,"operator=11");
    t.test_rel(onew[2],4.0,1.0e-10,"operator=12");
    cout << onew[0] << " " << onew[1] << " " << onew[2] << endl;
  }

  {
    /// Test push_back and pop_back()
    ovector ove(3);
    ove[0]=3.0;
    ove[1]=1.0;
    ove[2]=4.0;

    ove.push_back(1.0);
    ove.push_back(5.0);
    ove.push_back(9.0);
    ove.push_back(2.0);
    t.test_gen(ove.size()==7,"push size");
    t.test_rel(ove.pop_back(),2.0,1.0e-10,"pop_back");
    t.test_gen(ove.size()==6,"pop size");
    cout << ove << endl;
    ove.free();

    /// Test push_back and pop_back() 
    ofvector<3> ove2;
    ove2[0]=3.0;
    ove2[1]=1.0;
    ove2[2]=4.0;

    ove2.push_back(1.0);
    ove2.push_back(5.0);
    ove2.push_back(9.0);
    ove2.push_back(2.0);
    t.test_gen(ove2.size()==7,"push size");
    t.test_rel(ove2.pop_back(),2.0,1.0e-10,"pop");
    t.test_gen(ove2.size()==6,"pop size");
    cout << ove2 << endl;
    ove2.free();

    /// Try reversed subvectors
    ovector ors(3);
    ors[0]=3.0;
    ors[1]=1.0;
    ors[2]=4.0;
  
    ovector_subvector_reverse osr(ors,1,2);
    t.test_rel(osr[0],4.0,1.0e-5,"reverse subvector 1");
    t.test_rel(osr[1],1.0,1.0e-5,"reverse subvector 2");
  
    const ovector_const_subvector ocs(ors,1,2);
    cout << ocs[0] << endl;

    // Try reserve()

    ovector qx;
    cout << qx.size() << " " << qx.capacity() << endl;
    qx.reserve(10);
    cout << qx.size() << " " << qx.capacity() << endl;

    // Try const
    double ca[3]={2,1,0};
    const ovector_const_array oca(3,ca);
    const double cax=oca[2];
    cout << cax << endl;

    // Test sort_unique()
    ovector xo(10);
    for(size_t i=0;i<10;i++) xo[i]=i%4;
    xo.sort_unique();
    t.test_gen(xo.size()==4,"sort unique 0");
    t.test_rel(xo[0],0.0,1.0e-12,"sort unique 1");
  }

  {
    // Different way of testing push(), pop_back(), reserve(), erase(), etc.
    cout.width(30);
    cout << " " << " size owner block block_size contents" << endl;

    ovector dvec;

    show_data("Empty vector:",dvec);

    dvec.push_back(2.0);
    show_data("Push an element:",dvec);
    dvec.push_back(3.0);
    show_data("Push an element:",dvec);
    dvec.pop_back();
    show_data("Pop:",dvec);
    dvec.pop_back();
    show_data("Pop:",dvec);
    dvec.free();
    show_data("Free:",dvec);

    dvec.push_back(2.0);
    show_data("Push an element:",dvec);
    dvec.reserve(2);
    show_data("Reserve another spot:",dvec);
    dvec.push_back(3.0);
    show_data("Push an element:",dvec);
    dvec.erase(0);
    show_data("Erase the first element:",dvec);
    dvec.erase(0);
    show_data("Erase the first element:",dvec);
    dvec.free();
    show_data("Free:",dvec);

    dvec.reserve(3);
    show_data("Reserve three spots:",dvec);
    dvec.push_back(2.0);
    show_data("Push an element:",dvec);
    dvec.reserve(1);
    show_data("Reserve one spot:",dvec);
    dvec.reserve(4);
    show_data("Reserve four spots:",dvec);
    dvec.push_back(3.0);
    show_data("Push an element:",dvec);
    dvec.pop_back();
    show_data("Pop:",dvec);
    dvec.erase(0);
    show_data("Erase the first element:",dvec);
    dvec.free();
    show_data("Free:",dvec);
    dvec.allocate(2);
    dvec[0]=1.0;
    dvec[1]=3.0;
    show_data("Allocate two spots, and set:",dvec);
    dvec.reserve(3);
    show_data("Reserve three spots:",dvec);
    dvec.push_back(4.0);
    show_data("Push an element:",dvec);
    dvec.erase(2);
    show_data("Erase the last element:",dvec);
    dvec.pop_back();
    show_data("Pop:",dvec);
    dvec.free();
    show_data("Free:",dvec);
  }

  // Test const_iterator
  {
    ovector a(3);

    a[0]=3.0;
    a[1]=1.0;
    a[2]=4.0;

    for(ovector::const_iterator it=a.begin();it!=a.end();it++) {
      cout << *it << " ";
    }
    cout << endl;

    // Show that norm works
    cout << a.norm() << endl;
    
    // This doesn't yet work
    //std::sort(a.begin(),a.end());
    //cout << a << endl;
  }

  // Show that dynamic_cast to gsl_vector works
  {
    ovector a(3);

    a[0]=3.0;
    a[1]=1.0;
    a[2]=4.0;

    gsl_vector *gv=dynamic_cast<gsl_vector *>(&a);
    t.test_gen(gv!=0,"dynamic cast 1");
    cout << gsl_vector_get(gv,0) << endl;
    t.test_rel(gsl_vector_get(gv,0),3.0,1.0e-10,"dynamic cast 2");

    gsl_vector *gv2=a.get_gsl_vector();
    t.test_rel(gsl_vector_get(gv2,0),3.0,1.0e-10,"dynamic cast 3");
    const gsl_vector *gv3=a.get_gsl_vector_const();
    t.test_rel(gsl_vector_get(gv3,0),3.0,1.0e-10,"dynamic cast 4");

  }
  t.report();
  return 0;
}
