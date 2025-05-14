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
/** \file base_cuda.h
    \brief File for CUDA solver
*/
#ifndef O2SCL_BASE_CUDA_H
#define O2SCL_BASE_CUDA_H

#include <vector>

namespace o2scl {
  
  /** \brief Given a device ID, obtain the compute mode, major and
      minor version
   */
  int cuda_get_mode(int dev_id, int &mode, int &major, int &minor,
                    int verbose=0);
  
  /** \brief Obtain the number of CUDA cores given the major and
      minor version numbers
   */
  int cuda_cores_per_sm(int major, int minor);

  /** \brief Loop through \c dev_count devices, determining which
      has the maximum number of GFLOPs
   */
  int cuda_get_dev_max_gflops(int dev_count, int verbose=0);

  /** \brief Find the best GPU device and return the ID, compute
      mode, major, and minor version
   */
  int cuda_find_device_nothrow(int &dev_id, int &mode, int &major,
                               int &minor, int verbose=0);
  
}

#endif

