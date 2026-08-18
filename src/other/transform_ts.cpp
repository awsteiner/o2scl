/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2026, Andrew W. Steiner
  
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
#include <o2scl/transform.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr tm;
  tm.set_output_level(2);

  transform_base<tensor2<>> *tb;
  
  transform_quantile tq("normal",1000);
  transform_quantile tq2("uniform",1000);
  transform_standard ts;
  transform_minmax tmm;

  // Instantiate with a different type to make sure it works
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  transform_standard<ubmatrix> ts2;
  
  // ───────────────────────────────────────────────────────────────────
  // Create a couple different data sets to test with (n_samples=200,
  // n_features=3)
  // 
  // 0: exponential distribution with lambda=1
  // 1: uniform distribution between 5 and 15
  // 2: normal distribution (mu=100, sigma=20)                  
  
  const size_t N=200, F=3;
  std::vector<size_t> dims={N,F};
  o2scl::tensor2<> t_orig(N,F);
  
  std::mt19937 rng(42);
  std::exponential_distribution<> exp_d(1.0);
  std::uniform_real_distribution<> uni_d(5.0,15.0);
  std::normal_distribution<> nor_d(100.0,20.0);
  
  std::vector<size_t> idx(2);
  for (size_t s=0;s<N;s++) {
    t_orig(s,0)=exp_d(rng);
    t_orig(s,1)=uni_d(rng);
    t_orig(s,2)=nor_d(rng);
  }

  if (true) {
    o2scl::tensor2<> train, test;
    train_test_split(t_orig,0.62,train,test);
    tm.test_gen(train.size1()==76,"train size");
    tm.test_gen(test.size1()==124,"test size");
  }
  
  for(size_t i=0;i<4;i++) {

    if (i==0) {
      tb=&ts;
    } else if (i==1) {
      tb=&tmm;
    } else if (i==2) {
      tb=&tq;
    } else {
      tb=&tq2;
    }

    // Perform the forward and inverse transformation and
    // make sure the result is equal to the original.
    
    o2scl::tensor2<> t=t_orig;
    vector<double> orig_data=t.get_data();

    o2scl::tensor2<> t_norm=tb->fit_transform(t);
    tb->inverse_transform(t_norm);
    const vector<double> &t_norm_data=t_norm.get_data();

    if (i!=2) {
      tm.test_abs_vec(t_norm_data.size(),
                      t_norm_data,orig_data,1.0e-13,"max_err");
    } else {
      tm.test_abs_vec(t_norm_data.size(),
                      t_norm_data,orig_data,1.0e-3,"max_err");
    }

    // Perform the forward and inverse transformation using the
    // alternative interface and make sure the result is equal to the
    // original.
    
    o2scl::tensor2<> t_copy=t_orig;
    
    tb->fit(t_copy);
    tb->transform(t_copy);
    tb->inverse_transform(t_copy);
    const vector<double> &t_copy_data=t_copy.get_data();

    if (i!=2) {
      tm.test_abs_vec(t_copy_data.size(),
                      t_copy_data,orig_data,1.0e-13,"max_err");
    } else {
      tm.test_abs_vec(t_copy_data.size(),
                      t_copy_data,orig_data,1.0e-3,"max_err");
    }
    
  }
  
  tm.report();
  
  return 0;
}

