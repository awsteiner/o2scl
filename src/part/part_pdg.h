/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_PART_PDG_H
#define O2SCL_PART_PDG_H

#include <string>
#include <iostream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/part.h>

/** \file part_pdg.h
    \brief File defining \ref o2scl::thermo_tl and \ref o2scl::part_pdg_tl 
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A particle object compatible with the PDG
   */
  class part_pdg_db {

  public:

    part_pdg_db();
    
    typedef struct pdg_entry_s {
      int id;
      double mass;
      double mass_errp;
      double mass_errm;
      double width;
      double width_errp;
      double width_errm;
      std::string name;
      int charge;
    } pdg_entry;

    void output_text();

  protected:
    
    std::vector<pdg_entry> db;
    
  };
  
  /** \brief Particle from the PDG
   */
  template<class fp_t=double> class part_pdg_tl :
    public o2scl::part_tl<fp_t> {
    
  public:

  /// Copy constructor
  part_pdg_tl(const part_pdg_tl &p) {
    this->g=p.g;
    this->m=p.m;
    this->ms=p.ms;
    this->n=p.n;
    this->ed=p.ed;
    this->pr=p.pr;
    this->mu=p.mu;
    this->en=p.en;
    this->nu=p.nu;
    this->inc_rest_mass=p.inc_rest_mass;
    this->non_interacting=p.non_interacting;
  }

  /// Copy construction with operator=()
  part_pdg_tl &operator=(const part_pdg_tl &p) {
    if (this!=&p) {
      this->g=p.g;
      this->m=p.m;
      this->ms=p.ms;
      this->n=p.n;
      this->ed=p.ed;
      this->pr=p.pr;
      this->mu=p.mu;
      this->en=p.en;
      this->nu=p.nu;
      this->inc_rest_mass=p.inc_rest_mass;
      this->non_interacting=p.non_interacting;
    }
    return *this;
  }
    
  /// Make a particle of mass \c mass and degeneracy \c dof.
  part_pdg_tl(fp_t mass=0.0, fp_t dof=0.0) {
    this->m=mass; 
    this->ms=mass; 
    this->g=dof;
    
    this->non_interacting=true;
    this->inc_rest_mass=true;
  }    
  
  virtual ~part_pdg_tl() {
  }
  
  /// Set the mass \c mass and degeneracy \c dof.
  virtual void init(fp_t mass, fp_t dof) {
    this->m=mass; 
    this->ms=mass; 
    this->g=dof;
    return;
  }

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
