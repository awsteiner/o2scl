/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#ifndef O2SCL_QUARK_H
#define O2SCL_QUARK_H

/** \file quark.h
    \brief File defining \ref o2scl::quark_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>

#include <o2scl/fermion.h>

namespace o2scl {

  /** \brief Quark class
   */
  template<class fp_t=double> class quark_tl : public fermion_tl<fp_t> {
    
  public:
    
    /// Contribution to the bag constant
    fp_t B;
    
    /// Quark condensate
    fp_t qq;
    
    /// Create a boson with mass \c m and degeneracy \c g 
    quark_tl(fp_t mass=0, fp_t dof=0) : fermion_tl<double>(mass,dof) {
      qq=0;
      B=0;      
    }
    
    /// Return string denoting type ("quark")
    virtual const char *type() { return "quark"; }

  };

  typedef quark_tl<double> quark;

}

#endif
