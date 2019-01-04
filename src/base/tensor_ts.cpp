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
#include <o2scl/tensor.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/test_mgr.h>
#if O2SCL_HDF
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
using namespace o2scl_hdf;
#endif

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  // Test a rank one tensor
  tensor<> arr;
  size_t i1[1], j1[1];
  i1[0]=3;
  arr.resize(1,i1);
  for(size_t i=0;i<i1[0];i++) {
    j1[0]=i;
    arr.set(j1,((double)i));
    t.test_rel(arr.get(j1),((double)i),1.0e-12,"R1 element.");
  }
  t.test_gen(arr.total_size()==3,"R1 size.");
  arr.resize(0,i1);

  // Test a rank two tensor
  tensor<> mat;
  size_t i2[2], j2[2];
  i2[0]=4;
  i2[1]=3;
  mat.resize(2,i2);
  for(size_t i=0;i<i2[0];i++) {
    for(size_t j=0;j<i2[1];j++) {
      j2[0]=i;
      j2[1]=j;
      mat.set(j2,((double)i)+j);
      t.test_rel(mat.get(j2),((double)i+j),1.0e-12,"R2 element.");
    }
  }
  t.test_gen(mat.total_size()==12,"R2 size.");
  mat.resize(0,i2);
    
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
	ubvector_slice v=m3u.vector_slice(2,j3);
	for(size_t k=0;k<i3[2];k++) {
	  j3[2]=k;
	  t.test_rel(m3u.get(j3),v[k],1.0e-12,"vector slice 1.");
	}
      }
    }
    for(size_t i=0;i<i3[0];i++) {
      for(size_t k=0;k<i3[2];k++) {
	j3[0]=i;
	j3[2]=k;
	ubvector_slice v=m3u.vector_slice(1,j3);
	for(size_t j=0;j<i3[1];j++) {
	  j3[1]=j;
	  t.test_rel(m3u.get(j3),v[j],1.0e-12,"vector slice 2.");
	}
      }
    }
    for(size_t j=0;j<i3[1];j++) {
      for(size_t k=0;k<i3[2];k++) {
	j3[1]=j;
	j3[2]=k;
	ubvector_slice v=m3u.vector_slice(0,j3);
	for(size_t i=0;i<i3[0];i++) {
	  j3[0]=i;
	  t.test_rel(m3u.get(j3),v[i],1.0e-12,"vector slice 3.");
	}
      }
    }

#ifdef O2SCL_NEVER_DEFINED
    
    for(size_t j=0;j<i3[1];j++) {
      j3[1]=j;
      ubmatrix_view mt=m3.matrix_slice(0,j3);
      for(size_t m=0;m<i3[0];m++) {
	for(size_t n=0;n<i3[2];n++) {
	  j3[0]=m;
	  j3[2]=n;
	  t.test_rel(m3.get(j3),mt[m][n],1.0e-12,"matrix slice 0,2.");
	}
      }
    }
    
    for(size_t j=0;j<i3[0];j++) {
      j3[0]=j;
      ubmatrix_view mt=m3.matrix_slice(1,j3);
      for(size_t m=0;m<i3[1];m++) {
	for(size_t n=0;n<i3[2];n++) {
	  j3[1]=m;
	  j3[2]=n;
	  t.test_rel(m3.get(j3),mt[m][n],1.0e-12,"matrix slice 1,2.");
	}
      }
    }
    
#endif
    
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

#if O2SCL_HDF

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
    hdf_input(hf,tg2,"tens_grid_test");
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

#endif
    
  // -------------------------------------------------------

  {
    // Construct a sample rank-3 tensor
    tensor3<> temp;
    size_t sz[3]={2,2,2};
    temp.resize(3,sz);
    temp.set(0,0,0,1.0);
    temp.set(0,0,1,4.0);
    temp.set(0,1,0,9.0);
    temp.set(0,1,1,16.0);
    temp.set(1,0,0,25.0);
    temp.set(1,0,1,36.0);
    temp.set(1,1,0,49.0);
    temp.set(1,1,1,64.0);

    // First fix first index to one and output the corresponding matrix
    typedef std::function<double &(size_t,size_t)> data_t;
    data_t temp2=std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
			   (&tensor3<>::get),&temp,
			   1,std::placeholders::_1,std::placeholders::_2);
    t.test_rel(temp2(0,0),25.0,1.0e-12,"mat part 1");
    t.test_rel(temp2(0,1),36.0,1.0e-12,"mat part 2");
    t.test_rel(temp2(1,0),49.0,1.0e-12,"mat part 3");
    t.test_rel(temp2(1,1),64.0,1.0e-12,"mat part 4");

    // Now show how to use matrix_column_gen to select the second
    // column of the matrix created above. We fix last index to 1
    // (i.e. the second column), giving 36 and 64.
    matrix_column_gen<data_t> column=
      o2scl::matrix_column<data_t,matrix_column_gen<data_t> >(temp2,1);
    t.test_rel(column[0],36.0,1.0e-12,"mat column 1");
    t.test_rel(column[1],64.0,1.0e-12,"mat column 2");

  }
  
  if (true) {
    
    // Test rearrange_and_copy()
    tensor<> tx, tx2, tx3, tx2b, tx3b;
    size_t sz[5]={3,3,3,3,3};
    vector<size_t> ix_old, ix_new;

    // Create the test data
    tx.resize(5,sz);
    vector<double> data(243);
    for(size_t i=0;i<243;i++) {
      data[i]=((double)i);
    }
    tx.swap_data(data);

    hdf_file hf;
    hf.open_or_create("tensor_ts.o2");
    hf.setd_ten("rk5",tx);
    hf.close();

    // First test
    tx2=tx.rearrange_and_copy({ix_index(1),ix_reverse(4),
	  ix_fixed(3,2),ix_sum(0),ix_sum(2)},1);
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
    tx3=tx.rearrange_and_copy({ix_index(1),ix_range(4,1,0),
	  ix_fixed(3,2),ix_trace(0,2)},2);
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
  
  t.report();

  return 0;
}
