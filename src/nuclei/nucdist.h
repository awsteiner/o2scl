/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCDIST_H
#define O2SCL_NUCDIST_H

/** \file nucdist.h
    \brief File defining \ref o2scl::nucdist_set().
*/

#include <iostream>
#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Set a distribution of nuclei from a mass formula
      and a function string

      Given a nuclear mass formula \c nm, this adds nuclei to a
      <tt>std::vector</tt> object. The function \c expr, a function of
      \c Z and \z N, determines which nuclei will be added to the
      distribution.
  */
  void nucdist_set(std::vector<nucleus> &dist, nucmass &nm, 
		   std::string expr="1", int maxA=400,
		   bool include_neutron=false);

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
