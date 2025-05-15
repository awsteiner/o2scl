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
#include "solve_cuda.h"

#include <iostream>

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

using namespace std;
using namespace o2scl_linalg;

int linear_solver_LU_cuda_base::solve_base
(int n, const std::vector<double> &A, const std::vector<double> &b,
 std::vector<double> &x) {
  
  double *d_A;
  double *d_b;
  int *d_pivots;
  int *d_info;
  double *d_work;
  int work_size=0;

  x.resize(n);
  
  cusolverDnHandle_t cusolverH;
  if (cusolverDnCreate(&cusolverH)!=CUSOLVER_STATUS_SUCCESS) {
    return 1;
  }

  if (cudaMalloc(&d_A,n*n*sizeof(double))!=cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    return 5;
  }
  if (cudaMalloc(&d_b,n*sizeof(double))!=cudaSuccess) {
    cudaFree(d_A);
    cusolverDnDestroy(cusolverH);
    return 6;
  }
  if (cudaMalloc(&d_pivots,n*sizeof(int))!=cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cusolverDnDestroy(cusolverH);
    return 7;
  }
  if (cudaMalloc(&d_info,sizeof(int))!=cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cusolverDnDestroy(cusolverH);
    return 8;
  }

  if (cudaMemcpy(d_A,&(A[0]),n*n*sizeof(double),
                 cudaMemcpyHostToDevice)!=cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cusolverDnDestroy(cusolverH);
    return 9;
  }
  if (cudaMemcpy(d_b,&(b[0]),n*sizeof(double),cudaMemcpyHostToDevice)!=
      cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cusolverDnDestroy(cusolverH);
    return 10;
  }

  // Query working space
  if (cusolverDnDgetrf_bufferSize(cusolverH,n,n,d_A,n,&work_size)!=
      CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cusolverDnDestroy(cusolverH);
    return 2;
  }
  if (cudaMalloc(&d_work,work_size*sizeof(double))!=cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cusolverDnDestroy(cusolverH);
    return 11;
  }

  // LU decomposition (with partial pivoting)
  if (cusolverDnDgetrf(cusolverH,n,n,d_A,n,d_work,d_pivots,d_info)!=
      CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    return 3;
  }

  // Solve using LU
  if (cusolverDnDgetrs(cusolverH,CUBLAS_OP_N,n,1,d_A,n,d_pivots,
                       d_b,n,d_info)!=CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    return 4;
  }

  // Copy result
  if (cudaMemcpy(&(x[0]),d_b,n*sizeof(double),
                 cudaMemcpyDeviceToHost)!=cudaSuccess) {
    cudaFree(d_A);
    cudaFree(d_b);
    cudaFree(d_pivots);
    cudaFree(d_info);
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    return 12;
  }

  // Cleanup
  cudaFree(d_A);
  cudaFree(d_b);
  cudaFree(d_pivots);
  cudaFree(d_info);
  cudaFree(d_work);
  cusolverDnDestroy(cusolverH);

  return 0;
}

int main(void) {

  vector<double> A={1,0,0,0,2,0,0,0,3};
  vector<double> b={4,5,6}, x(3);

  linear_solver_LU_cuda lslc;
  int ret=lslc.solve(3,A,b,x);
  
  cout << ret << " " << x[0] << " " << x[1] << " " << x[2] << endl;

  return 0;
}
