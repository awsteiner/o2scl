/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2025, Andrew W. Steiner

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
#include <iostream>
#include <ctime>

#include <o2scl/invert_fast.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // We choose a nearly diagonal positive symmetric matrix which
  // is easy to invert
  for(size_t n=10;n<10000;n*=2) {

    size_t mult=10;
    
    matrix_invert_cholesky_fast micf;
    
    cout << ((double)n) << " " << std::flush;
    
    for(size_t ell=0;ell<3;ell++) {
      
      if (ell==0) {
        micf.mode=micf.force_o2;
      } else if (ell==1) {
        micf.mode=micf.force_arma;
      } else {
        micf.mode=micf.force_cuda;
      }
      
      if (ell!=0 || n<=640) {
        tensor2<> t1[mult], t2[mult];
        vector<size_t> size={n,n};
        
        struct timespec ts0;
        timespec_get(&ts0,TIME_UTC);
        
        for(size_t k=0;k<mult;k++) {
          t1[k].resize(2,size);
          t2[k].resize(2,size);
          for(size_t i=0;i<n;i++) {
            for(size_t j=0;j<n;j++) {
              if (i==j) t1[k](i,j)=((double)(i+2));
              else t1[k](i,j)=1.0e-2*exp(-2.0*pow(((double)i)+((double)j),2.0));
            }
          }
        }
        
        struct timespec ts1;
        timespec_get(&ts1,TIME_UTC);
        
        double diff0=(double)(ts1.tv_sec-ts0.tv_sec);
        double ndiff0=(double)(ts1.tv_nsec-ts0.tv_nsec);
        //cout << (diff0+ndiff0*1.0e-9)/((double)mult) << " " << std::flush;
        
        for(size_t k=0;k<mult;k++) {
          micf.invert(n,t1[k],t2[k]);
        }
        
        struct timespec ts2;
        timespec_get(&ts2,TIME_UTC);
        
        double diff=(double)(ts2.tv_sec-ts1.tv_sec);
        double ndiff=(double)(ts2.tv_nsec-ts1.tv_nsec);
        cout << (diff+ndiff*1.0e-9)/((double)mult) << " " << std::flush;
      } else {
        cout << 0.0 << " " << std::flush;
      }
    }
    cout << endl;

  }

  t.report();
  return 0;
}

