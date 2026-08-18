/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2025-2026, Andrew W. Steiner

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

int o2scl_linalg::cholesky_decomp_cuda_base(const size_t n,
                                            std::vector<double> &A) {
  
  // Note that the function cusolverDnDpotrf presumes that the matrix
  // is lower triangular and stored in column-major order, but the
  // argument A is in row-major order, so we use CUBLAS_FILL_MODE_UPPER.
  
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
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,&work_size);
  
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
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,d_work,
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
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 8;
  }

  cudaFree(d_work);
  cudaFree(devInfo);
  cusolverDnDestroy(cusolverH);
  cudaFree(d_A);
  
  return 0;
}

int matrix_invert_det_cholesky_cuda_base::invert
(size_t n, const std::vector<double> &A,
 std::vector<double> &A_inv) {

  // Make a copy of the original matrix so that we
  // can modify it
  std::vector<double> Acopy(n*n);

  for (size_t j=0;j<n;j++) {
    for (size_t i=0;i<n;i++) {
      if (i<j) {
        Acopy[i*n+j]=A[j*n+i];
      } else {
        Acopy[i*n+j]=A[i*n+j];
      }
    }
  }

  // Perform the Cholesky decomposition
  int ret=o2scl_linalg::cholesky_decomp_cuda_base(n,Acopy);
  if (ret != 0) return ret;
  
  
  // Re-upload the Cholesky factor to device for potri
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) return 10;

  cudaStat=cudaMemcpy(d_A,Acopy.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 11;
  }

  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 12;
  }

  // Get buffer size for potri
  int work_size=0;
  cusolver_status=cusolverDnDpotri_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,&work_size);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 13;
  }

  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 14;
  }

  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 15;
  }

  // Invert using Cholesky result
  cusolver_status=cusolverDnDpotri
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,d_work,
     work_size,devInfo);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 16;
  }

  // Copy result back to host
  A_inv.resize(n*n);
  cudaStat=cudaMemcpy(A_inv.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 17;
  }

  // Symmetrize the result (potri only fills the upper triangle in
  // column major order)
  for (size_t j=0;j<n;j++) {
    for (size_t i=0;i<n;i++) {
      if (i<j) {
        A_inv[i*n+j]=A_inv[j*n+i];
      }
    }
  }

  // Clean up
  cudaFree(d_A);
  cudaFree(d_work);
  cudaFree(devInfo);
  cusolverDnDestroy(cusolverH);
  
  return 0;
}

int matrix_invert_det_cholesky_cuda_base::invert_det
(size_t n, const std::vector<double> &A,
 std::vector<double> &A_inv, double &A_det) {
  
  // Make a copy of the original matrix since the original
  // will be destroyed
  std::vector<double> chol(n*n);
  for(size_t i=0;i<n*n;i++) {
    chol[i]=A[i];
  }

  // Perform the Cholesky decomposition
  int ret=o2scl_linalg::cholesky_decomp_cuda_base(n,chol);
  if (ret != 0) return ret;

  // Compute determinant from the diagonal of the Cholesky factor
  double sqrt_det=1.0;
  for (size_t i=0;i<n;i++) sqrt_det*=chol[i*n+i];
  A_det=sqrt_det*sqrt_det;

  // Re-upload the Cholesky factor to device for potri
  double *d_A=0;
  cudaError_t cudaStat=cudaMalloc((void**)&d_A,n*n*sizeof(double));
  if (cudaStat != cudaSuccess) return 10;

  cudaStat=cudaMemcpy(d_A,chol.data(),n*n*sizeof(double),
                      cudaMemcpyHostToDevice);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_A);
    return 11;
  }

  // Create cuSolver handle
  cusolverDnHandle_t cusolverH=0;
  cusolverStatus_t cusolver_status=cusolverDnCreate(&cusolverH);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_A);
    return 12;
  }

  // Get buffer size for potri
  int work_size=0;
  cusolver_status=cusolverDnDpotri_bufferSize
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,&work_size);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 13;
  }

  // Allocate workspace and devInfo
  double *d_work=0;
  int *devInfo=0;
  cudaStat=cudaMalloc((void**)&d_work,work_size*sizeof(double));
  if (cudaStat != cudaSuccess) {
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 14;
  }

  cudaStat=cudaMalloc((void**)&devInfo,sizeof(int));
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 15;
  }

  // Invert using Cholesky result
  cusolver_status=cusolverDnDpotri
    (cusolverH,CUBLAS_FILL_MODE_UPPER,n,d_A,n,d_work,
     work_size,devInfo);
  if (cusolver_status != CUSOLVER_STATUS_SUCCESS) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 16;
  }

  // Copy result back to host
  A_inv.resize(n*n);
  cudaStat=cudaMemcpy(A_inv.data(),d_A,n*n*sizeof(double),
                      cudaMemcpyDeviceToHost);
  if (cudaStat != cudaSuccess) {
    cudaFree(d_work);
    cudaFree(devInfo);
    cusolverDnDestroy(cusolverH);
    cudaFree(d_A);
    return 17;
  }

  // Symmetrize the result (potri only fills the upper triangle)
  for (int i=0;i<(int)n;++i) {
    for (int j=i+1;j<(int)n;++j) {
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
  
double matrix_invert_det_cholesky_cuda_base::det
(size_t n, const std::vector<double> &A) {

  // Make a copy of the original matrix since the original
  // will be destroyed
  std::vector<double> chol(n*n);
  for(size_t i=0;i<n*n;i++) {
    chol[i]=A[i];
  }

  // Perform the Cholesky decomposition
  int ret=o2scl_linalg::cholesky_decomp_cuda_base(n,chol);
  if (ret != 0) return 0.0;

  // Compute determinant from the diagonal of the Cholesky factor
  double sqrt_det=1.0;
  for (size_t i=0;i<n;i++) sqrt_det*=chol[i*n+i];

  return sqrt_det*sqrt_det;
}
  
int matrix_invert_det_cholesky_cuda_base::invert_inplace
(size_t n, std::vector<double> &A) {
  return invert(n,A,A);
}

