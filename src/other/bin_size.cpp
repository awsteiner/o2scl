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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/bin_size.h>

using namespace std;
using namespace o2scl;

int bin_size::calc_bin(double a1, double a2, int naa, double &bl, double &bh, 
		       int &nb, double &bwid) {
  
  // AWS - 7/27/07 - I've thrown in some comments here and there to
  // try to better document what's going on, but I don't fully
  // understand it yet.

  double al, ah;
  if (a1<a2) {
    al=a1;
    ah=a2;
  } else {
    al=a2;
    ah=a1;
  }
  if (al==ah) ah=al+1;

  bool goto20=false;

  // Set these three variables to zero to avoid uninitialized 
  // variable warnings
  int na=0, lg=0;
  double sigrnd=0.0;

  do {
    
    if (goto20 || naa!=-1 || bwid<=0.0) {

      if (!goto20) {
	//if (cern_mode) 
	na=naa-1;
	//else na=naa;
	if (na<1) na=1;
      }
      goto20=false;

      // Compute the width of the interval and
      // the value of sigrnd

      double awid=(ah-al)/na;
      lg=((int)log10(awid));
      if (awid<1.0) lg--;
      double sigfig=awid*pow(10.0,-lg);
      if (sigfig<2) {
	sigrnd=2.0;
      } else if (sigfig<2.5) {
	sigrnd=2.5;
      } else if (sigfig<5.0) {
	sigrnd=5.0;
      } else {
	lg++;
	sigrnd=1.0;
      }

    }

    // Using sigrnd, compute the bin size and the corresponding
    // limits bl, and bh. Compute the new number of bins

    bwid=sigrnd*pow(10.0,lg);
    double alb=al/bwid;
    int lwid=((int)alb);
    if (alb<0.0) lwid--;
    bl=bwid*lwid;
    alb=ah/bwid+1.0;
    int kwid=((int)alb);
    if (alb<0.0) kwid--;
    bh=bwid*kwid;
    nb=kwid-lwid;

    // The user requested only one bin, so we're done
    if (naa==-1) return 0;

    if (naa>5) {

      // Check to make sure we didn't get exactly twice as many bins
      // as the user requested. (Why?) If not we're done, otherwise, try
      // again

      if (2*nb!=naa) return 0;
      na++;
      goto20=true;

    }

  } while (goto20==true);

  // If the user requested more than one bin, or if one bin was
  // requested and one bin was computed, then we're done

  if (naa>1 || nb==1) return 0;

  // Otherwise, make the bin size larger and return
  // only 1 bin (Why?)

  bwid*=2.0;
  nb=1;

  return 0;
}

