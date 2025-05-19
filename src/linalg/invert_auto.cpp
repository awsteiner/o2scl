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

#include <o2scl/invert_auto.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;

int o2scl_linalg::cholesky_decomp_cuda(const size_t n,
                                       std::vector<double> &A) {
#ifdef O2SCL_SET_CUDA
  int ret=cholesky_decomp_cuda_base(n,A);
  if (ret!=0) {
    std::string err=((std::string)"Error number ")+o2scl::itos(ret)+
      " in o2scl_linalg::cholesky_decomp_cuda().";
    O2SCL_ERR(err.c_str(),o2scl::exc_einval);
  }
#else
  O2SCL_ERR("Cuda support not included in this O2scl installation.",
            o2scl::exc_eunimpl);
#endif
  return ret;
}

