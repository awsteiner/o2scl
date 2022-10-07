/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

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
    tx2=rearrange_and_copy2<double,tensor<>>(tx,{ix_index(1),ix_reverse(4),
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
    tx3=rearrange_and_copy2<double,tensor<>>(tx,{ix_index(1),ix_range(4,1,0),
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
