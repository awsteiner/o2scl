#include <iostream>
#include <cusolverDn.h>
#include <cuda_runtime.h>

void checkCuda(cudaError_t result) {
    if (result != cudaSuccess) {
        std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
        exit(EXIT_FAILURE);
    }
}

void checkCusolver(cusolverStatus_t result) {
    if (result != CUSOLVER_STATUS_SUCCESS) {
        std::cerr << "cuSOLVER Error: " << result << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main() {
    cusolverDnHandle_t cusolverH;
    checkCusolver(cusolverDnCreate(&cusolverH));

    const int N = 3;  // Matrix size
    double A[N * N] = {3.0, 2.0, 4.0, 
                       2.0, 0.0, 2.0, 
                       4.0, 2.0, 3.0};  // Symmetric matrix

    double *d_A, *d_W;
    int *d_info, Lwork;
    double *d_work;

    checkCuda(cudaMalloc((void**)&d_A, N * N * sizeof(double)));
    checkCuda(cudaMalloc((void**)&d_W, N * sizeof(double)));
    checkCuda(cudaMalloc((void**)&d_info, sizeof(int)));

    checkCuda(cudaMemcpy(d_A, A, N * N * sizeof(double), cudaMemcpyHostToDevice));

    // Query workspace size
    checkCusolver(cusolverDnDsyevd_bufferSize(cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, N, d_A, N, d_W, &Lwork));

    checkCuda(cudaMalloc((void**)&d_work, Lwork * sizeof(double)));

    // Compute eigenvalues and eigenvectors
    checkCusolver(cusolverDnDsyevd(cusolverH, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, N, d_A, N, d_W, d_work, Lwork, d_info));

    // Copy results back to host
    double W[N];
    checkCuda(cudaMemcpy(W, d_W, N * sizeof(double), cudaMemcpyDeviceToHost));

    // Print eigenvalues
    std::cout << "Eigenvalues: ";
    for (int i = 0; i < N; i++) {
        std::cout << W[i] << " ";
    }
    std::cout << std::endl;

    // Clean up
    checkCuda(cudaFree(d_A));
    checkCuda(cudaFree(d_W));
    checkCuda(cudaFree(d_info));
    checkCuda(cudaFree(d_work));
    checkCusolver(cusolverDnDestroy(cusolverH));

    return 0;
}
