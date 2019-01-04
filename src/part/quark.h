/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_QUARK_H
#define O2SCL_QUARK_H

/** \file quark.h
    \brief File defining \ref o2scl::quark
*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>

#include <o2scl/fermion.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Quark class

  */
  class quark : public fermion {

  public:

    /// Contribution to the bag constant
    double B;

    /// Quark condensate
    double qq;

    /// Create a boson with mass \c m and degeneracy \c g 
    quark(double m=0.0, double g=0.0);
    
    /// Return string denoting type ("quark")
    virtual const char *type() { return "quark"; }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
