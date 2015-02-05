/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_BIN_SIZE_H
#define O2SCL_BIN_SIZE_H

#include <iostream>
#include <cmath>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Determine bin size (CERNLIB)

      This is adapted from the KERNLIB routine \c binsiz.f written by
      F. James. 
      
      This class computes an appropriate set of histogram bins given the
      upper and lower limits of the data and the maximum number of bins.
      The bin width is always an integral power of ten times 1, 2, 2.5 
      or 5. The bin width may also be specified by the user, in which
      case the class only computes the appropriate limits.
      
      \note This class is not working yet.

      \future Finish this.
  */
  class bin_size {

  public:

    bin_size() {
      cern_mode=true;
    }

    /// (default true) 
    bool cern_mode;

    /** \brief Compute bin size
	
	- \c al - Lower limit of data
	- \c ah - Upper limit of data
	- \c na - Maximum number of bins desired. 
	- \c bl - Lower limit (BL<=AL)
	- \c bh - Upper limit (BH>=AH)
	- \c nb - Number of bins determined by BINSIZ (NA/2<NB<=NA)
	- \c bwid - Bin width (BH-BL)/NB
	
	If \c na=0 or \c na=-1, this function always makes exactly 
	one bin.
	
	If \c na=1, this function takes \c bwid as input and determines
	only \c bl, \c hb, and \c nb. This is especially useful when it is
	desired to have the same bin width for several histograms (or
	for the two axes of a scatter-plot).
	
	If \c al > \c ah, this function takes \c al to be the upper
	limit and \c ah to be the lower limit, so that in fact \c al
	and \c ah may appear in any order. They are not changed by
	calc_bin(). If \c al = \c ah, the lower limit is taken to be \c al,
	and the upper limit is set to \c al+1.
	
	If \ref cern_mode is true (which is the default) the starting
	guess for the number of bins is \c na-1. Otherwise, the 
	starting guess for the number of bins is \c na. 
    */
    int calc_bin(double al, double ah, int na, double &bl, double &bh, 
		 int &nb, double &bwid);
      
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
