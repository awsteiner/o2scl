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
#include "invert_cuda.h"

#include <cuda_runtime.h>
#include <cusolverDn.h>

using namespace o2scl_linalg;

int o2scl_linalg::cholesky_decomp_cuda(const size_t n, std::vector<double> &A) {

  // Note that the function cusolverDnDpotrf presumes that the matrix
  // is stored in column-major, rather than row-major order. This
  // isn't a problem for this function, however, because the matrix
  // is assumed to be symmetric. 
  
  // Allocate device memory
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) {
    return 1;
  }
  
  // Copy data to device
  cudaStat=cudaMemcpy(d_A,A.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 2;
  }
  
  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 3;
  }
  
  // Get buffer size for potrf (Cholesky)
  int work_size=0;
  cusolver_status=cusolverDnDpotrf_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,&work_size);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 4;
  }
  
  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 5;
  }
  
  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 6;
  }
  
  // Cholesky decomposition (A=L*L^T)
  cusolver_status=cusolverDnDpotrf
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 7;
  }
  
  // Copy result back to host
  cudaStat=cudaMemcpy(A.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);
  if (cudaStat != cudaSuccess) {
    return 8;
  }
  
  return 0;
}

int matrix_invert_det_cholesky_cuda::invert
(size_t n, const std::vector<double> &A,
 std::vector<double> &A_inv) {

  // Make a copy of the original matrix since the original
  // will be destroyed?
  std::vector<double> Acopy(n*n);
  for(size_t i=0;i<n*n;i++) {
    Acopy[i]=A[i];
  }
  
  // Allocate device memory
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) {
    return 1;
  }
  
  // Copy data to device
  cudaStat=cudaMemcpy(d_A,Acopy.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 2;
  }
  
  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 3;
  }
  
  // Get buffer size for potrf (Cholesky)
  int work_size=0;
  cusolver_status=cusolverDnDpotrf_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,&work_size);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 4;
  }
  
  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 5;
  }
  
  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 6;
  }
  
  // Cholesky decomposition (A=L*L^T)
  cusolver_status=cusolverDnDpotrf
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 7;
  }
  
  // Invert using Cholesky result
  cusolver_status=cusolverDnDpotri
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 8;
  }
  
  // Copy result back to host
  A_inv.resize(n*n);
  cudaStat=cudaMemcpy(A_inv.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);
  if (cudaStat != cudaSuccess) {
    return 9;
  }
  
  // Symmetrize the result (only lower triangle is filled)
  for (int i=0;i<n;++i) {
    for (int j=i+1;j<n;++j) {
      A_inv[i*n+j]=A_inv[j*n+i];
    }
  }
  
  // Clean up
  cudaFree(d_A);
  cudaFree(d_work);
  cudaFree(devInfo);
  cusolverDnDestroy(cusolverH);
  
  return 0;
}

int matrix_invert_det_cholesky_cuda::invert_det
(size_t n, const std::vector<double> &A,
 std::vector<double> &A_inv, double &A_det) {
  
  // Allocate device memory
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) {
    return 1;
  }
  
  // Copy data to device
  cudaStat=cudaMemcpy(d_A,A.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 2;
  }
  
  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 3;
  }
  
  // Get buffer size for potrf (Cholesky)
  int work_size=0;
  cusolver_status=cusolverDnDpotrf_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,&work_size);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 4;
  }
  
  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 5;
  }
  
  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 6;
  }
  
  // Cholesky decomposition (A=L*L^T)
  cusolver_status=cusolverDnDpotrf
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 7;
  }

  // Copy Cholesky decomposition back to host to compute
  // determinant
  std::vector<double> chol(n*n);
  cudaStat=cudaMemcpy(chol.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);

  // Compute determinant
  double sqrt_det=1.0;
  for(size_t i=0;i<n;i++) sqrt_det*=chol[i*n+i];
  A_det=sqrt_det*sqrt_det;
  
  // Invert using Cholesky result
  cusolver_status=cusolverDnDpotri
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 8;
  }
  
  // Copy result back to host
  A_inv.resize(n*n);
  cudaStat=cudaMemcpy(A_inv.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);
  if (cudaStat != cudaSuccess) {
    return 9;
  }
  
  // Symmetrize the result (only lower triangle is filled)
  for (int i=0;i<n;++i) {
    for (int j=i+1;j<n;++j) {
      A_inv[i*n+j]=A_inv[j*n+i];
    }
  }
  
  // Clean up
  cudaFree(d_A);
  cudaFree(d_work);
  cudaFree(devInfo);
  cusolverDnDestroy(cusolverH);
  
  return 0;
}
  
double matrix_invert_det_cholesky_cuda::det
(size_t n, const std::vector<double> &A) {

  // Allocate device memory
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) {
    return 1;
  }
  
  // Copy data to device
  cudaStat=cudaMemcpy(d_A,A.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 2;
  }
  
  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 3;
  }
  
  // Get buffer size for potrf (Cholesky)
  int work_size=0;
  cusolver_status=cusolverDnDpotrf_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,&work_size);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 4;
  }
  
  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 5;
  }
  
  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 6;
  }
  
  // Cholesky decomposition (A=L*L^T)
  cusolver_status=cusolverDnDpotrf
    (cusolverH,CUBLAS_FILL_MODE_LOWER,n,d_A,n,d_work,
     work_size,devInfo);
  
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 7;
  }

  // Copy Cholesky decomposition back to host to compute
  // determinant
  std::vector<double> chol(n*n);
  cudaStat=cudaMemcpy(chol.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);

  // Compute determinant
  double sqrt_det=1.0;
  for(size_t i=0;i<n;i++) sqrt_det*=chol[i*n+i];
  double A_det=sqrt_det*sqrt_det;
  
  cudaFree(d_work);
  cudaFree(devInfo);
  cusolverDnDestroy(cusolverH);
  cudaFree(d_A);

  return A_det;
}
  
int matrix_invert_det_cholesky_cuda::invert_inplace
(size_t n, std::vector<double> &A) {
  return invert(n,A,A);
}

