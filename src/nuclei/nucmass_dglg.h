/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCMASS_DGLG_H
#define O2SCL_NUCMASS_DGLG_H

/** \file nucmass_dglg.h
    \brief File defining \ref o2scl::nucmass_dglg
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucmass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nuclear properties from Delaroche et al. 

      See \ref Delaroche10 .
  */
  class nucmass_dglg : public nucmass_table {
    
  public:
    
    nucmass_dglg(std::string model="", bool external=false);

    ~nucmass_dglg();
    
    /** \brief Entry structure
	
    */
    class entry {

    public:

      /// Proton number
      int Z;
      /// Neutron number
      int N;
      /// HFB energy minimum (MeV) 
      double EHFB;
      /// β deformation at HFB energy minimum  
      double BMIN;
      /// γ deformation (degree) at HFB energy minimum
      double GMIN;
      /// Charge radius (fm) at HFB energy minimum  
      double RCHFB;
      /// Point proton radius (fm) at HFB energy minimum
      double RPHFB;
      /// Neutron radius (fm) at HFB energy minimum
      double RNHFB;
      /// CHFB+5DCH ground state (g.s.) energy (MeV) 
      double EABS;
      /// Correlation energy (MeV) 
      double ECORR;
      /// Mean g.s. beta deformation  
      double BET01;
      /// Mean g.s. gamma deformation (degree) 
      double GAM01;
      /// Variance of g.s. beta deformation  
      double DELB01;
      /// Variance of g.s. gamma deformation (degree) 
      double DELG01;
      /// Yrast 2^+ state energy (MeV) 
      double E21;
      /// Yrast 4^+ state energy (MeV) 
      double E41;
      /// Yrast  6^+ state energy (MeV) 
      double E61;
      /// First 0^+ excited state energy (MeV) 
      double E02;
      /// Second 2^+ excited state energy (MeV) 
      double E22;
      /// Third 2^+ excited state energy (MeV) 
      double E23;
      /// P(K=0) for yrast 2^+ state (%) 
      double PK0_2_1;
      /// P(K=2) for second 2^+ excited state (%) 
      double PK2_2_2;
      /// P(K=2) for third 2^+ excited state (%) 
      double PK2_2_3;
      /// B(E2; 2^+_1 --> 0^+_1) (e**2 b**2) 
      double BE2_2_1_0_1;
      /// B(E2; 2^+_3 --> 0^+_1) (e**2 b**2) 
      double BE2_2_3_0_1;
      /// B(E2; 2^+_1 --> 0^+_2) (e**2 b**2) 
      double BE2_2_1_0_2;
      /// B(E2; 4^+_1 --> 2^+_1) (e**2 b**2) 
      double BE2_4_1_2_1;
      /// B(E2; 2^+_3 --> 2^+_1) (e**2 b**2) 
      double BE2_2_3_2_1;
      /// B(E2; 2^+_3 --> 0^+_2) (e**2 b**2) 
      double BE2_2_3_0_2;
      /// CHFB+5DCH charge radius (fm) 
      double RC5DCH;
      /// CHFB+5DCH point proton radius (fm)
      double RP5DCH;
      /// CHFB+5DCH neutron radius (fm)
      double RN5DCH;
      /// CHFB+5DCH squared E0 matrix element  
      double ROE0TH;
      /// For each Z : Neutron number at 5DCH proton drip-line 
      int NMIN;
      /// For each Z : Neutron number at 5DCH neutron drip-line 
      int NMAX;

    };
  
    /// Return the type, \c "nucmass_dglg".
    virtual const char *type() { return "nucmass_dglg"; }

    /// Returns true if data has been loaded
    bool is_loaded() { return (n>0); }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Return number of entries
    int get_nentries() { return n; }
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The number of entries (about 3000).
    int n;
    
    /// The reference for the original data
    std::string reference;
    
    /// The array containing the mass data of length n
    entry *mass;

    /// The last table index for caching
    int last;
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
