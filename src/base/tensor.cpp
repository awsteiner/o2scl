/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/tensor.h>

using namespace std;
using namespace o2scl;

index_spec o2scl::ix_index2(size_t ix) {
  return index_spec(index_spec::index,ix,0,0,0.0,0.0,0.0);
}

index_spec o2scl::ix_fixed2(size_t ix, size_t fix) {
  return index_spec(index_spec::fixed,ix,fix,0,0.0,0.0,0.0);
}

index_spec o2scl::ix_sum2(size_t ix) {
  return index_spec(index_spec::sum,ix,0,0,0.0,0.0,0.0);
}

index_spec o2scl::ix_trace2(size_t ix, size_t ix2) {
  return index_spec(index_spec::trace,ix,ix2,0,0.0,0.0,0.0);
}
  
index_spec o2scl::ix_reverse2(size_t ix) {
  return index_spec(index_spec::reverse,ix,0,0,0.0,0.0,0.0);
}
  
index_spec o2scl::ix_range2(size_t ix, size_t start, size_t end) {
  return index_spec(index_spec::range,ix,start,end,0.0,0.0,0.0);
}

index_spec o2scl::ix_interp2(size_t ix, double v) {
  return index_spec(index_spec::interp,ix,0,0,v,0.0,0.0);
}
  
index_spec o2scl::ix_grid2(size_t ix, double begin, double end, 
			  size_t n_bins, bool log) {
  if (log) {
    return index_spec(index_spec::grid,ix,n_bins,1,begin,end,0.0);
  }
  return index_spec(index_spec::grid,ix,n_bins,0,begin,end,0.0);
}

index_spec o2scl::ix_gridw2(size_t ix, double begin, double end, 
			   double width, bool log) {
  if (log) {
    return index_spec(index_spec::gridw,ix,0,1,begin,end,width);
  }
  return index_spec(index_spec::gridw,ix,0,0,begin,end,width);
}
