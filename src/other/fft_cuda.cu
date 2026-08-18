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
#include "fft_cuda.h"
#include <cufft.h>
#include <iostream>

using namespace std;
using namespace o2scl;

int o2scl::vector_forward_fft_cuda_nothrow
(size_t N, const std::vector<double> &data,
 std::vector<std::complex<double>> &fft) {
 
  // Allocate host memory
  fft.resize(N/2+1);
  
  // Allocate device memory
  double *d_input;
  cufftDoubleComplex *d_output;
  cudaMalloc((void**)&d_input,sizeof(double)*N);
  cudaMalloc((void**)&d_output,sizeof(cufftDoubleComplex)*(N/2+1));
  
  // Copy data to device
  cudaMemcpy(d_input,(&(data[0])),
             sizeof(double)*N,cudaMemcpyHostToDevice);
  
  // Create FFT plan
  cufftHandle plan;
  if (cufftPlan1d(&plan,N,CUFFT_D2Z,1) != CUFFT_SUCCESS) {
    //std::cerr << "CUFFT error: Plan creation failed" << std::endl;
    return -1;
  }
  
  // Execute FFT
  if (cufftExecD2Z(plan,d_input,d_output) != CUFFT_SUCCESS) {
    //std::cerr << "CUFFT error: ExecD2Z failed" << std::endl;
    return -2;
  }
  
  // Copy result back to host
  void *fftp=&(fft[0]);
  cudaMemcpy(fftp,d_output,sizeof(cufftDoubleComplex)*(N/2+1),
             cudaMemcpyDeviceToHost);

  // Cleanup
  cufftDestroy(plan);
  cudaFree(d_input);
  cudaFree(d_output);
  
  return 0;
}

int o2scl::vector_backward_fft_cuda_nothrow
(size_t N, const std::vector<std::complex<double>> &fft,
 std::vector<double> &data) {
 
  // Allocate host memory
  int Ndata=(N-1)*2;
  data.resize(Ndata);
  
  // Allocate device memory
  cufftDoubleComplex *d_input;
  
  double *d_output;
  cudaMalloc((void**)&d_input,sizeof(cufftDoubleComplex)*N);
  cudaMalloc((void**)&d_output,sizeof(double)*Ndata);
  
  // Copy data to device
  void *fftp=(void *)&(fft[0]);
  cudaMemcpy(d_input,fftp,
             sizeof(cufftDoubleComplex)*N,cudaMemcpyHostToDevice);
  
  // Create FFT plan
  cufftHandle plan;
  if (cufftPlan1d(&plan,Ndata,CUFFT_Z2D,1) != CUFFT_SUCCESS) {
    //std::cerr << "CUFFT error: Plan creation failed" << std::endl;
    return -1;
  }
  
  // Execute FFT
  if (cufftExecZ2D(plan,d_input,d_output) != CUFFT_SUCCESS) {
    //std::cerr << "CUFFT error: ExecD2Z failed" << std::endl;
    return -2;
  }
  
  // Copy result back to host
  cudaMemcpy(&data[0],d_output,sizeof(double)*Ndata,
             cudaMemcpyDeviceToHost);

  // Cleanup
  cufftDestroy(plan);
  cudaFree(d_input);
  cudaFree(d_output);
  
  return 0;
}

#ifdef O2SCL_NEVER_DEFINED

int main(void) {

  vector<double> data={0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0};
  vector<complex<double>> fft2;
  
  int ret=vector_forward_fft_cuda_nothrow(8,data,fft2);
  std::cout << "ret: " << ret << std::endl;
  for(size_t i=0;i<fft2.size();i++) {
    cout << fft2[i].real() << " " << fft2[i].imag() << endl;
  }
  cout << endl;

  for(size_t i=0;i<data.size();i++) {
    data[i]=2.0;
  }

  ret=vector_backward_fft_cuda_nothrow(fft2.size(),fft2,data);
  std::cout << "ret: " << ret << std::endl;
  for(size_t i=0;i<data.size();i++) {
    cout << data[i] << endl;
  }
  cout << endl;
  
  return 0;
}

/*
  FFT output:
  Freq 0: 1       + 0j
  Freq 1: 0.7071  - 0.7071j
  Freq 2: 0       - 1j
  Freq 3: -0.7071 - 0.7071j
  Freq 4: -1      + 0j
*/

#endif
