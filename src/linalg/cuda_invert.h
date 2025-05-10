#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <iostream>
#include <vector>

class matrix_invert_det_cholesky_cuda {
public:
  
  int invert(size_t n, const std::vector<double> &A,
             std::vector<double> &A_inv);
  
  int invert_det(size_t n, const std::vector<double> &A,
                 std::vector<double> &A_inv, double &A_det);
  
  double det(size_t n, const std::vector<double> &A);
  
  int invert_inplace(size_t n, std::vector<double> &A);
  
};

