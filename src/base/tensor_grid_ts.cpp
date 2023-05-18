/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/tensor.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace o2scl_hdf;

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  {
  // Test a rank three tensor
  tensor_grid<> m3;
  size_t i3[3], j3[3], k3[3];
  i3[0]=4;
  i3[1]=3;
  i3[2]=3;
  m3.resize(3,i3);
  for(size_t i=0;i<i3[0];i++) {
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	j3[0]=i;
	j3[1]=j;
	j3[2]=k;
	m3.set(j3,((double)i)+j+k);
	size_t jp=m3.pack_indices(j3);
	m3.unpack_index(jp,k3);
	t.test_rel(m3.get(j3),((double)i+j+k),1.0e-12,"R2 element.");
	t.test_gen(j3[0]==k3[0],"R3 pack/unpack.");
	t.test_gen(j3[1]==k3[1],"R3 pack/unpack.");
	t.test_gen(j3[2]==k3[2],"R3 pack/unpack.");
      }
    }
  }
  t.test_gen(m3.total_size()==36,"R3 size.");
  
  // Create a grid, and set data for later interpolation
  std::vector<double> grid;
  size_t j4[3];
  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  grid.push_back(4.0);
  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  grid.push_back(1.0);
  grid.push_back(2.0);
  grid.push_back(3.0);
  m3.set_grid_packed(grid);
  for(size_t i=0;i<i3[0];i++) {
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	double x=m3.get_grid(0,i);
	double y=m3.get_grid(1,j);
	double z=m3.get_grid(2,k);
	j3[0]=i;
	j3[1]=j;
	j3[2]=k;
	m3.set(j3,2.0*x*x-y-3.0*z*z);
	size_t ix;
	ix=m3.pack_indices(j3);
	m3.unpack_index(ix,j4);
	t.test_gen(j4[0]==i,"pack/unpack 1");
	t.test_gen(j4[1]==j,"pack/unpack 2");
	t.test_gen(j4[2]==k,"pack/unpack 3");
      }
    }
  }
  }

  // -------------------------------------------------------
  // Test interpolation with a tensor object
  // built upon std:: vectors

  {

    // Create a sample tensor_grid object with a grid
    tensor_grid<> m3;
    size_t i3[3], j3[3], k3[3];
    i3[0]=4;
    i3[1]=3;
    i3[2]=3;
    m3.resize(3,i3);
    std::vector<double> grid;
    size_t j4[3];
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    grid.push_back(4.0);
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    m3.set_grid_packed(grid);
    for(size_t i=0;i<i3[0];i++) {
      for(size_t j=0;j<i3[1];j++) {
	for(size_t k=0;k<i3[2];k++) {
	  double x=m3.get_grid(0,i);
	  double y=m3.get_grid(1,j);
	  double z=m3.get_grid(2,k);
	  j3[0]=i;
	  j3[1]=j;
	  j3[2]=k;
	  m3.set(j3,2.0*x*x-y-3.0*z*z);
	}
      }
    }

    std::vector<double> v(3), res, res3;
    double res2;

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3.interp_linear_vec0(v,res);
      m3.interp_linear_vec(v,0,res3);
      v[0]=1.0;
      res2=m3.interp_linear(v);
      t.test_rel(res[0],res2,1.0e-12,"interp_linear_vec0 1");
      t.test_rel(res3[0],res2,1.0e-12,"interp_linear_vec 1");
      v[0]=2.0;
      res2=m3.interp_linear(v);
      t.test_rel(res[1],res2,1.0e-12,"interp_linear_vec0 2");
      t.test_rel(res3[1],res2,1.0e-12,"interp_linear_vec 2");
      v[0]=3.0;
      res2=m3.interp_linear(v);
      t.test_rel(res[2],res2,1.0e-12,"interp_linear_vec0 3");
      t.test_rel(res3[2],res2,1.0e-12,"interp_linear_vec 3");
      v[0]=4.0;
      res2=m3.interp_linear(v);
      t.test_rel(res[3],res2,1.0e-12,"interp_linear_vec0 4");
      t.test_rel(res3[3],res2,1.0e-12,"interp_linear_vec 4");
    }

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3.interp_linear_vec(v,1,res3);
      v[1]=1.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[0],res2,1.0e-12,"interp_linear_vec 5");
      v[1]=2.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[1],res2,1.0e-12,"interp_linear_vec 6");
      v[1]=3.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[2],res2,1.0e-12,"interp_linear_vec 7");
    }

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3.interp_linear_vec(v,2,res3);
      v[2]=1.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[0],res2,1.0e-12,"interp_linear_vec 8");
      v[2]=2.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[1],res2,1.0e-12,"interp_linear_vec 9");
      v[2]=3.0;
      res2=m3.interp_linear(v);
      t.test_rel(res3[2],res2,1.0e-12,"interp_linear_vec 10");
    }

  }

  // -------------------------------------------------------
  // Test slicing and interpolation with a tensor object
  // built upon ublas vectors

  {

    // Create a sample tensor_grid object with a grid
    tensor_grid<ubvector,ubvector_size_t> m3u;
    size_t i3[3], j3[3], k3[3];
    i3[0]=4;
    i3[1]=3;
    i3[2]=3;
    m3u.resize(3,i3);
    std::vector<double> grid;
    size_t j4[3];
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    grid.push_back(4.0);
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    m3u.set_grid_packed(grid);
    for(size_t i=0;i<i3[0];i++) {
      for(size_t j=0;j<i3[1];j++) {
	for(size_t k=0;k<i3[2];k++) {
	  double x=m3u.get_grid(0,i);
	  double y=m3u.get_grid(1,j);
	  double z=m3u.get_grid(2,k);
	  j3[0]=i;
	  j3[1]=j;
	  j3[2]=k;
	  m3u.set(j3,2.0*x*x-y-3.0*z*z);
	}
      }
    }

    ubvector v(3), res, res3;
    double res2;

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3u.interp_linear_vec0(v,res);
      m3u.interp_linear_vec(v,0,res3);
      v[0]=1.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res[0],res2,1.0e-12,"u interp_linear_vec0 1");
      t.test_rel(res3[0],res2,1.0e-12,"u interp_linear_vec 1");
      v[0]=2.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res[1],res2,1.0e-12,"u interp_linear_vec0 2");
      t.test_rel(res3[1],res2,1.0e-12,"u interp_linear_vec 2");
      v[0]=3.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res[2],res2,1.0e-12,"u interp_linear_vec0 3");
      t.test_rel(res3[2],res2,1.0e-12,"u interp_linear_vec 3");
      v[0]=4.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res[3],res2,1.0e-12,"u interp_linear_vec0 4");
      t.test_rel(res3[3],res2,1.0e-12,"u interp_linear_vec 4");
    }

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3u.interp_linear_vec(v,1,res3);
      v[1]=1.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[0],res2,1.0e-12,"u interp_linear_vec 5");
      v[1]=2.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[1],res2,1.0e-12,"u interp_linear_vec 6");
      v[1]=3.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[2],res2,1.0e-12,"u interp_linear_vec 7");
    }

    v[0]=3.1;
    v[1]=2.2;
    v[2]=1.3;

    if (true) {
      m3u.interp_linear_vec(v,2,res3);
      v[2]=1.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[0],res2,1.0e-12,"u interp_linear_vec 8");
      v[2]=2.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[1],res2,1.0e-12,"u interp_linear_vec 9");
      v[2]=3.0;
      res2=m3u.interp_linear(v);
      t.test_rel(res3[2],res2,1.0e-12,"u interp_linear_vec 10");
    }
    
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;

    // Test vector slicing
    for(size_t i=0;i<i3[0];i++) {
      for(size_t j=0;j<i3[1];j++) {
	j3[0]=i;
	j3[1]=j;
	ubvector_slice vx=m3u.vector_slice(2,j3);
	for(size_t k=0;k<i3[2];k++) {
	  j3[2]=k;
	  t.test_rel(m3u.get(j3),vx[k],1.0e-12,"vector slice 1.");
	}
      }
    }
    for(size_t i=0;i<i3[0];i++) {
      for(size_t k=0;k<i3[2];k++) {
	j3[0]=i;
	j3[2]=k;
	ubvector_slice vx=m3u.vector_slice(1,j3);
	for(size_t j=0;j<i3[1];j++) {
	  j3[1]=j;
	  t.test_rel(m3u.get(j3),vx[j],1.0e-12,"vector slice 2.");
	}
      }
    }
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	j3[1]=j;
	j3[2]=k;
	ubvector_slice vx=m3u.vector_slice(0,j3);
	for(size_t i=0;i<i3[0];i++) {
	  j3[0]=i;
	  t.test_rel(m3u.get(j3),vx[i],1.0e-12,"vector slice 3.");
	}
      }
    }

    // -------------------------------------------------------
    // Test successive and linear interpolation
    
    t.test_gen(m3u.get_size(0)==4,"size 1");
    t.test_gen(m3u.get_size(1)==3,"size 2");
    t.test_gen(m3u.get_size(2)==3,"size 3");
    
    double vals[3]={2.5,2.5,1.5};
    t.test_rel(m3u.interpolate(vals),m3u.interp_linear(vals),
	       1.0e-12,"interp 1");
    double vals2[3]={1.8,2.2,1.0};
    t.test_rel(m3u.interpolate(vals),m3u.interp_linear(vals),
	       1.0e-12,"interp 2");
    double vals3[3]={2.0,2.0,1.1};
    t.test_rel(m3u.interpolate(vals),m3u.interp_linear(vals),
	       1.0e-12,"interp 3");
  }

  // -------------------------------------------------------
  // Test tensor_grid3 HDF5 I/O

  {

    // Construct a small rank-3 tensor
    tensor_grid3<> tg(3,2,1);
    double grid2[6]={4,5,6,7,8,9};
    tg.set_grid_packed(grid2);
    for(size_t j=0;j<3;j++) {
      for(size_t k=0;k<2;k++) {
	for(size_t ell=0;ell<1;ell++) {
	  tg.set(j,k,ell,((double)(j+k+ell)));
	}
      }
    }

    // Save it to a file
    hdf_file hf;
    hf.open_or_create("tens_grid.o2");
    hdf_output(hf,tg,"tens_grid_test");
    hf.close();

    // Read that file into a new tensor object
    tensor_grid3<> tg2;
    hf.open("tens_grid.o2");
    string name_temp="tens_grid_test";
    hdf_input(hf,tg2,name_temp);
    hf.close();

    // Check that the second tensor is the same as the first
    t.test_gen(tg2.get_rank()==3,"rank");
    t.test_rel(tg2.get_grid(0,1),5.0,1.0e-12,"grid 1");
    t.test_rel(tg2.get_grid(2,0),9.0,1.0e-12,"grid 2");
    for(size_t j=0;j<3;j++) {
      for(size_t k=0;k<2;k++) {
	for(size_t ell=0;ell<1;ell++) {
	  double x=tg.get(j,k,ell);
	  t.test_rel(((double)(j+k+ell)),x,1.0e-12,"element");
	}
      }
    }

  }

  if (true) {
    
    // Test rearrange_and_copy()
    tensor_grid<> tx, tx2, tx3, tx2b, tx3b;
    size_t sz[5]={3,3,3,3,3};
    vector<size_t> ix_old, ix_new;

    // Create the test data
    tx.resize(5,sz);
    vector<double> grid;
    for(size_t i=0;i<12;i++) {
      grid.push_back(((double)i));
    }
    grid.push_back(0.1);
    grid.push_back(0.01);
    grid.push_back(0.001);
    tx.set_grid_packed(grid);
    vector<double> data(243);
    for(size_t i=0;i<243;i++) {
      data[i]=((double)i);
    }
    tx.swap_data(data);

    hdf_file hf;
    hf.open_or_create("tensor_ts.o2");
    hf.setd_ten("rk5",tx);
    hf.close();

    // First test. ix_interp(1,3.0) has the same effect as
    // ix_fixed(1,0) and ix_grid(4,...) has the same effect
    // as ix_index(4)
    tx2=grid_rearrange_and_copy<tensor_grid<>,double>
      (tx,{ix_interp(1,3.0),ix_gridw(4,10.0,15.0,1.0),
        ix_fixed(3,2),ix_sum(0),ix_sum(2)},3);
    cout << tx2.get_size(0) << " " << tx2.get_rank() << endl;
    for(size_t k=0;k<tx2.get_size(0);k++) {
      cout << tx2.get_grid(0,k) << " ";
    }
    cout << endl;

    if (false) {
      size_t sz2b[2]={3,3};
      tx2b.resize(2,sz2b);
      for(size_t i1=0;i1<3;i1++) {
	for(size_t i2=0;i2<3;i2++) {
	  double val=0.0;
	  for(size_t i3=0;i3<3;i3++) {
	    for(size_t i4=0;i4<3;i4++) {
	      ix_old={i3,i1,i4,2,2-i2};
	      val+=tx.get(ix_old);
	    }
	  }
	  ix_new={i1,i2};
	  tx2b.set(ix_new,val);
	}
      }
      t.test_gen(tx2==tx2b,"rearrange 1");
      
      // Second test
      tx3=grid_rearrange_and_copy<tensor_grid<>,double>
        (tx,{ix_index(1),ix_range(4,1,0),
          ix_interp(3,9.5),ix_trace(0,2)},2);
      size_t sz3b[2]={3,2};
      tx3b.resize(2,sz3b);
      for(size_t i1=0;i1<3;i1++) {
	for(size_t i2=0;i2<2;i2++) {
	  double val=0.0;
	  for(size_t i3=0;i3<3;i3++) {
	    ix_old={i3,i1,i3,2,1-i2};
	    val+=tx.get(ix_old);
	  }
	  ix_new={i1,i2};
	  tx3b.set(ix_new,val);
	}
      }
      t.test_gen(tx3==tx3b,"rearrange 2");
    }
  }

  /*
  {
    index_spec i1=ix_index2(4);
    ix_index i2(4);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_fixed2(4,1);
    ix_fixed i2(4,1);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_trace2(4,1);
    ix_trace i2(4,1);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_sum2(4);
    ix_sum i2(4);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_reverse2(4);
    ix_reverse i2(4);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_range2(4,1,2);
    ix_range i2(4,1,2);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_interp2(4,1.0);
    ix_interp i2(4,1.0);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_grid2(4,10.0,15.0,1.1);
    ix_grid i2(4,10.0,15.0,1.1);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_grid2(4,10.0,15.0,1.1,1);
    ix_grid i2(4,10.0,15.0,1.1,1);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_gridw2(4,10.0,15.0,1.1);
    ix_gridw i2(4,10.0,15.0,1.1);
    cout << i1.equal(i2) << endl;
  }
  {
    index_spec i1=ix_gridw2(4,10.0,15.0,1.1,1);
    ix_gridw i2(4,10.0,15.0,1.1,1);
    cout << i1.equal(i2) << endl;
  }
  */
  
  t.report();

  return 0;
}
