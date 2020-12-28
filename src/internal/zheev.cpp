/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#include <iostream>

using namespace std;
using namespace o2scl;

/**
   The lapack zheev function to compute the eigenvalues
   of a complex Hermitian matrix
*/
extern "C" {
  void zheev_(char *c1, char *c2, int *n, double *a, int *n2, double *w, 
	      double *work, int *lwork, double *rwork, int *info);
}

int main(void) {

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);

  // Create a Hermitian matrix
  double *a=new double[16*2];

  a[(0+4*0)*2]=1.0; a[(0+4*0)*2+1]=0.0;
  a[(0+4*1)*2]=2.0; a[(0+4*1)*2+1]=-1.0;
  a[(0+4*2)*2]=3.0; a[(0+4*2)*2+1]=-1.0;
  a[(0+4*3)*2]=4.0; a[(0+4*3)*2+1]=-1.0;
  a[(1+4*1)*2]=2.0; a[(1+4*1)*2+1]=0.0;
  a[(1+4*2)*2]=3.0; a[(1+4*2)*2+1]=-2.0;
  a[(1+4*3)*2]=4.0; a[(1+4*3)*2+1]=-2.0;
  a[(2+4*2)*2]=3.0; a[(2+4*2)*2+1]=0;
  a[(2+4*3)*2]=4.0; a[(2+4*3)*2+1]=-3.0;
  a[(3+4*3)*2]=4.0; a[(3+4*3)*2+1]=0;

  char vec='V';
  char uplo='U';
  int N=4;
  int info=0;
  int nb=64;
  int lwork=(nb+1)*4;
  
  double *w=new double[4];
  double *work=new double[lwork*2];
  double *rwork=new double[10];

  zheev_(&vec,&uplo,&N,a,&N,w,work,&lwork,rwork,&info);
  
  // Print out the results

  cout << "Info: " << info << endl;
  
  cout << "Eigenvalues: " << endl;
  cout << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;
  cout << endl;

  cout << "Eigenvectors: " << endl;
  cout << "1: " << endl;
  cout << a[(0+4*0)*2] << " " 
       << a[(0+4*1)*2] << " " 
       << a[(0+4*2)*2] << " " 
       << a[(0+4*3)*2] << endl;
  cout << a[(0+4*0)*2+1] << " " 
       << a[(0+4*1)*2+1] << " " 
       << a[(0+4*2)*2+1] << " " 
       << a[(0+4*3)*2+1] << endl;
  cout << endl;
  cout << "2: " << endl;
  cout << a[(1+4*0)*2] << " " 
       << a[(1+4*1)*2] << " " 
       << a[(1+4*2)*2] << " " 
       << a[(1+4*3)*2] << endl;
  cout << a[(1+4*0)*2+1] << " " 
       << a[(1+4*1)*2+1] << " " 
       << a[(1+4*2)*2+1] << " " 
       << a[(1+4*3)*2+1] << endl;
  cout << endl;
  cout << "3: " << endl;
  cout << a[(2+4*0)*2] << " " 
       << a[(2+4*1)*2] << " " 
       << a[(2+4*2)*2] << " " 
       << a[(2+4*3)*2] << endl;
  cout << a[(2+4*0)*2+1] << " " 
       << a[(2+4*1)*2+1] << " " 
       << a[(2+4*2)*2+1] << " " 
       << a[(2+4*3)*2+1] << endl;
  cout << endl;
  cout << "4: " << endl;
  cout << a[(3+4*0)*2] << " " 
       << a[(3+4*1)*2] << " " 
       << a[(3+4*2)*2] << " " 
       << a[(3+4*3)*2] << endl;
  cout << a[(3+4*0)*2+1] << " " 
       << a[(3+4*1)*2+1] << " " 
       << a[(3+4*2)*2+1] << " " 
       << a[(3+4*3)*2+1] << endl;
  cout << endl;
  
  delete[] a;
  delete[] w;
  delete[] work;
  delete[] rwork;

  ubmatrix_cx u(4,4);

  u.set(0,0,1,0);
  u.set(1,0,2,-1);
  u.set(2,0,3,-1);
  u.set(3,0,4,-1);

  u.set(0,1,0,0);
  u.set(1,1,2,0);
  u.set(2,1,3,-2);
  u.set(3,1,4,-2);

  u.set(0,2,0,0);
  u.set(1,2,0,0);
  u.set(2,2,3,0);
  u.set(3,2,4,-3);

  u.set(0,3,0,0);
  u.set(1,3,0,0);
  u.set(2,3,0,0);
  u.set(3,3,4,0);

  double *a2=&(u.get_ptr(0,0)->dat[0]);
  
  double *w2=new double[4];
  double *work2=new double[lwork*2];
  double *rwork2=new double[10];
  
  zheev_(&vec,&uplo,&N,a2,&N,w2,work2,&lwork,rwork2,&info);
  
  cout << "Info: " << info << endl;
  
  cout << "Eigenvalues: " << endl;
  cout << w2[0] << " " << w2[1] << " " << w2[2] << " " << w2[3] << endl;
  cout << endl;

  /*
    ubvector_cx evec[4];
    for(size_t i=0;i<4;i++) {
    evec[i]=ubmatrix_cx_col(a2,i);
    }
  */

  cout << "Eigenvectors: " << endl;
  cout << "1: " << endl;
  cout << u.get(0,0).dat[0] << " "
       << u.get(1,0).dat[0] << " "
       << u.get(2,0).dat[0] << " "
       << u.get(3,0).dat[0] << endl;
  cout << u.get(0,0).dat[1] << " "
       << u.get(1,0).dat[1] << " "
       << u.get(2,0).dat[1] << " "
       << u.get(3,0).dat[1] << endl;
  cout << endl;
  /*
  cout << "1: " << endl;
  cout << evec[0].get(0).dat[0] << " "
       << evec[0].get(1).dat[0] << " "
       << evec[0].get(2).dat[0] << " "
       << evec[0].get(3).dat[0] << endl;
  cout << evec[0].get(0).dat[1] << " "
       << evec[0].get(1).dat[1] << " "
       << evec[0].get(2).dat[1] << " "
       << evec[0].get(3).dat[1] << endl;
  cout << endl;
  */
  cout << "2: " << endl;
  cout << u.get(0,1).dat[0] << " "
       << u.get(1,1).dat[0] << " "
       << u.get(2,1).dat[0] << " "
       << u.get(3,1).dat[0] << endl;
  cout << u.get(0,1).dat[1] << " "
       << u.get(1,1).dat[1] << " "
       << u.get(2,1).dat[1] << " "
       << u.get(3,1).dat[1] << endl;
  cout << endl;
  cout << "3: " << endl;
  cout << u.get(0,2).dat[0] << " "
       << u.get(1,2).dat[0] << " "
       << u.get(2,2).dat[0] << " "
       << u.get(3,2).dat[0] << endl;
  cout << u.get(0,2).dat[1] << " "
       << u.get(1,2).dat[1] << " "
       << u.get(2,2).dat[1] << " "
       << u.get(3,2).dat[1] << endl;
  cout << endl;
  cout << "4: " << endl;
  cout << u.get(0,3).dat[0] << " "
       << u.get(1,3).dat[0] << " "
       << u.get(2,3).dat[0] << " "
       << u.get(3,3).dat[0] << endl;
  cout << u.get(0,3).dat[1] << " "
       << u.get(1,3).dat[1] << " "
       << u.get(2,3).dat[1] << " "
       << u.get(3,3).dat[1] << endl;
  cout << endl;

  return 0;
}
