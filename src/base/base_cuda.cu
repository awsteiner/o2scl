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
#include "base_cuda.h"

#include <cstdio>
#include <cstdint>
#include <cuda_runtime.h>

using namespace o2scl;

// AWS, 5/13/25: We avoid cout and use printf because I've had 
// warnings regarding ABI compatibility.

int o2scl::cuda_cores_per_sm(int major, int minor) {
    
  // Defines for GPU Architecture types (using the SM version to
  // determine the # of cores per SM
  typedef struct {
    int SM;  // 0xMm (hexidecimal notation), M = SM Major version,
    // and m = SM minor version
    int Cores;
  } sSMtoCores;
    
  sSMtoCores nGpuArchCoresPerSM[] = {
    {0x30, 192},
    {0x32, 192},
    {0x35, 192},
    {0x37, 192},
    {0x50, 128},
    {0x52, 128},
    {0x53, 128},
    {0x60,  64},
    {0x61, 128},
    {0x62, 128},
    {0x70,  64},
    {0x72,  64},
    {0x75,  64},
    {0x80,  64},
    {0x86, 128},
    {0x87, 128},
    {0x89, 128},
    {0x90, 128},
    {-1, -1}};

  int index = 0;
    
  while (nGpuArchCoresPerSM[index].SM != -1) {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchCoresPerSM[index].Cores;
    }
      
    index++;
  }
    
  // If we don't find the values, we default use the previous one
  // to run properly
  printf("MapSMtoCores for SM %d.%d is undefined."
         "  Default to use %d Cores/SM\n",
         major, minor, nGpuArchCoresPerSM[index - 1].Cores);
  return nGpuArchCoresPerSM[index - 1].Cores;
}

int o2scl::cuda_get_mode(int dev_id, int &mode, int &major, int &minor,
                         int verbose) {
  
  mode=-1;
  major=0;
  minor=0;
  
  int ret=cudaDeviceGetAttribute(&mode,
                                 cudaDevAttrComputeMode,dev_id);
  if (ret!=0) {
    if (verbose>0) {
      printf("cuda_get_mode(): ");
      printf("Failed to obtain device compute mode.\n");
    }
    return -4;
  }
  ret=cudaDeviceGetAttribute(&major,cudaDevAttrComputeCapabilityMajor,
                             dev_id);
  if (ret!=0) {
    if (verbose>0) {
      printf("cuda_get_mode(): ");
      printf("Failed to get device major mode.\n");
    }
    return -5;
  }
  ret=cudaDeviceGetAttribute(&minor,cudaDevAttrComputeCapabilityMinor,
                             dev_id);
  if (ret!=0) {
    if (verbose>0) {
      printf("cuda_get_mode(): ");
      printf("Failed to get device minor mode.\n");
    }
    return -6;
  }

  return 0;
}

int o2scl::cuda_get_dev_max_gflops(int dev_count, int verbose) {

  int curr_device=0;
  int sm_per_multiproc=0;
  int max_perf_device=0;
  int devices_prohibited=0;

  uint64_t max_compute_perf=0;

  int ret;
  
  while (curr_device<dev_count) {
    
    // Get the compute mode
    int mode=-1, major=0, minor=0;
    ret=cuda_get_mode(curr_device,mode,major,minor,verbose);
    if (ret!=0) return ret;
    
    // If this GPU is not running on Compute Mode prohibited,
    // then we can add it to the list
    
    if (mode==cudaComputeModeProhibited) {
      
      devices_prohibited++;
      if (verbose>1) {
        printf("cuda_get_dev_max_gflops(): ");
        printf("Device %d compute mode prohibited.\n",curr_device);
      }
      
    } else {
      
      if (major==9999 && minor==9999) {
        sm_per_multiproc=1;
      } else {
        sm_per_multiproc=cuda_cores_per_sm(major,minor);
      }
      int multi_proc_count=0;
      int clock_rate=0;

      ret=cudaDeviceGetAttribute(&multi_proc_count,
                                 cudaDevAttrMultiProcessorCount,
                                 curr_device);
      if (ret!=0) {
        if (verbose>0) {
          printf("Failed to get multiprocessor count.");
        }
        return -1;
      }
                       
      cudaError_t result=cudaDeviceGetAttribute
        (&clock_rate,cudaDevAttrClockRate,curr_device);
      
      if (result!=cudaSuccess) {
        if (result==cudaErrorInvalidValue) {
          // If cudaDevAttrClockRate attribute is not supported we set
          // clock_rate as 1.
          clock_rate=1;
        } else {
          if (verbose>0) {
            printf("Failed to get clock rate.");
          }
          return -12;
        }
      }
      
      uint64_t compute_perf=(uint64_t)multi_proc_count*
        sm_per_multiproc*clock_rate;

      if (verbose>1) {
        printf("cuda_get_dev_max_gflops(): ");
        printf("Device %d: compute performance %lu,\n",curr_device,
               compute_perf);
        printf("  multi_proc_count: %d, sm_per_multiproc: %d\n"
               ,multi_proc_count,sm_per_multiproc);
        printf("  clock_rate: %d\n",clock_rate);
      }
      
      if (compute_perf>max_compute_perf) {
        max_compute_perf=compute_perf;
        max_perf_device=curr_device;
      }

    }

    ++curr_device;
  }

  if (devices_prohibited == dev_count) {
    if (verbose>0) {
      printf("All devices prohibited.");
    }
    return -13;
  }

  return max_perf_device;
}

int o2scl::cuda_find_device_nothrow(int &dev_id, int &mode, int &major,
                                    int &minor, int verbose) {
  int ret;

  // Make sure that we have a GPU somewhere
  int dev_count;
  
  ret=cudaGetDeviceCount(&dev_count);
  if (ret!=0) {
    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Attempt to obtain device count failed.\n");
    }
    return -2;
  }
  
  if (dev_count==0) {
    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Device count is zero.\n");
    }
    return -10;
  }

  // If a device ID was requested, then select it
  if (dev_id>=0) {

    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Attempting to set device %d.\n",dev_id);
    }

    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Device count is: %d.\n",dev_count);
    }

    if (dev_id>dev_count-1) {
      if (verbose>0) {
        printf("cuda_find_device_nothrow(): ");
        printf("Device count is %d but user requested device %d.\n",
               dev_count,dev_id);
      }
      return -3;
    }

  } else {

    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Looking for fastest device.\n");
    }
    dev_id=cuda_get_dev_max_gflops(dev_count,verbose);
    
  }

  // Get the compute mode
  ret=cuda_get_mode(dev_id,mode,major,minor,verbose);
  if (ret!=0) return ret;
  
  // Check that the requested device can run CUDA
  if (mode==cudaComputeModeProhibited) {
    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Device reported compute mode was prohibited.\n");
    }
    return -7;
  }
  
  if (major<1) {
    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Device reported major mode less than 1.\n" );
    }
    return -8;
  }
  
  // Now finally set the CUDA device
  
  ret=cudaSetDevice(dev_id);
  if (ret!=0) {
    if (verbose>0) {
      printf("cuda_find_device_nothrow(): ");
      printf("Failed to set device %d.\n",dev_id);
    }
    return -9;
  }
  
  if (verbose>0) {
    printf("cuda_find_device_nothrow(): Success!\n");
  }
  return 0;
}

