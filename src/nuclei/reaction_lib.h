/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2015, Andrew W. Steiner
  
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
#ifndef REACTION_LIB_H
#define REACTION_LIB_H

/** \file reaction_lib.h
    \brief File defining \ref o2scl::nuclear_reaction
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>
#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A simple nuclear reaction specification

      This class is very experimental.
  */
  class nuclear_reaction {

  public:

    /// Chapter
    size_t chap;
    /// Names of the participating nuclei
    std::string name[6];
    /// Reference
    std::string ref;
    /// Type of rate (resonant/non-resonant/weak)
    char type;
    /// Forward or reverse
    char rev;
    /// Q value
    double Q;
    /// Coefficients
    double a[7];
    /// Proton number of participating nuclei
    size_t Z[6];
    /// Mass number of participating nuclei
    size_t A[6];
    /// Isomer designation of participating nuclei
    size_t isomer[6];

    nuclear_reaction() {
    }

    /** \brief Convert the reaction to a string for screen output
     */
    std::string to_string() {
      std::ostringstream outs;
      outs.width(2);
      outs << chap << ": ";
      nucmass_info nmi;
      for(size_t i=0;i<6;i++) {
	if (name[i].length()>0 && A[i]>0) {
	  if (Z[i]==0 && A[i]==1) outs << "n ";
	  else if (Z[i]==1 && A[i]==1) outs << "p ";
	  else if (Z[i]==1 && A[i]==2) outs << "d ";
	  else if (Z[i]==1 && A[i]==3) outs << "t ";
	  else outs << nmi.Ztoel(Z[i]) << A[i] << " ";
	  if (i==0 && (chap==1 || chap==2 || chap==3 || chap==11)) {
	    outs << "-> ";
	  } else if (i==1 && (chap==4 || chap==5 || chap==6 || chap==7)) {
	    outs << "-> ";
	  } else if (i==2 && (chap==8 || chap==9)) {
	    outs << "-> ";
	  } else if (i==3 && chap==10) {
	    outs << "-> ";
	  } else if (i<5 && name[i+1].length()>0 && A[i+1]>0) {
	    outs << "+ ";
	  }
	} else {
	  i=6;
	}
      }
      return outs.str();
    }

    /// Clear the rate
    int clear() {
      chap=0;
      ref="";
      type=0;
      rev=0;
      Q=0.0;
      for(size_t i=0;i<6;i++) {
	name[i]="";
	a[i]=0.0;
	Z[i]=0;
	A[i]=0;
	isomer[i]=0;
      }
      a[6]=0.0;
      
      return 0;
    }

    /// Copy constructor
    nuclear_reaction(const nuclear_reaction &nr) {
      chap=nr.chap;
      ref=nr.ref;
      type=nr.type;
      rev=nr.rev;
      Q=nr.Q;
      for(size_t i=0;i<6;i++) {
	name[i]=nr.name[i];
	a[i]=nr.a[i];
	Z[i]=nr.Z[i];
	A[i]=nr.A[i];
	isomer[i]=nr.isomer[i];
      }
      a[6]=nr.a[6];
    }
    
    /// Copy constructor
    nuclear_reaction &operator=(const nuclear_reaction &nr) {
      
      // Check for self-assignment
      if (this==&nr) return *this;
      
      chap=nr.chap;
      ref=nr.ref;
      type=nr.type;
      rev=nr.rev;
      Q=nr.Q;
      for(size_t i=0;i<6;i++) {
	name[i]=nr.name[i];
	a[i]=nr.a[i];
	Z[i]=nr.Z[i];
	A[i]=nr.A[i];
	isomer[i]=nr.isomer[i];
      }
      a[6]=nr.a[6];

      return *this;
    }

    /** \brief Compute the reaction rate from the temperature in units of 
	\f$ 10^9 K \f$ 
    */
    double rate(double T9) {
      double ret;
      double T913=cbrt(T9);
      ret=exp(a[0]+a[1]/T9+a[2]/T913+a[3]*T913+a[4]*T9+
	      a[5]*T9*T913*T913)*pow(T9,a[6]);
      return ret;
    }
    
  };

  /** \brief Simple reaction library 

      This class is very experimental.
      
      Units:
      - Chapters 1,2,3, and 11: 1/s
      - Chapters 4,5,6, and 7: cm^3/g/s
      - Chapter 8 and 9: cm^6/g^2/s
      - Chapter 10: cm^9/g^3/s

      Chapters:
      - 1: nuc1 -> nuc2
      - 2: nuc1 -> nuc2 + nuc3
      - 3: nuc1 -> nuc2 + nuc3 + nuc4
      - 4: nuc1 + nuc2 -> nuc3
      - 5: nuc1 + nuc2 -> nuc3 + nuc4
      - 6: nuc1 + nuc2 -> nuc3 + nuc4 + nuc5
      - 7: nuc1 + nuc2 -> nuc3 + nuc4 + nuc5 + nuc6
      - 8: nuc1 + nuc2 + nuc3 -> nuc4
      - 9: nuc1 + nuc2 + nuc3 -> nuc4 + nuc5
      - 10: nuc1 + nuc2 + nuc3 + nuc4 -> nuc5 + num6
      - 11: nuc1 -> nuc2 + nuc3 + nuc4 + nuc5

      Original FORTRAN format:
      \verbatim
      FORMAT(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5) 
      FORMAT(4e13.6) 
      FORMAT(3e13.6) 
      \endverbatim
  */
  class reaction_lib {

  public:

    /// The library
    std::vector<nuclear_reaction> lib;
    
    /** \brief Read from a file in the REACLIB2 format

	\note This function does not check that the chapter numbers
	are correct for the subsequent reaction.
     */
    int read_file_reaclib2(std::string fname);
    
    /** \brief Find a set of nuclear reactions in a specified chapter

	\comment
	This function doesn't make any assumptions about the ordering of 
	the rates and the chapters.
	\endcomment
     */
    int find_in_chap(std::vector<nuclear_reaction> &nrl,
		     size_t chap, std::string nuc1, std::string nuc2="", 
		     std::string nuc3="", std::string nuc4="", 
		     std::string nuc5="", std::string nuc6="");
    
  protected:
    
    /// \name Storage for the find function
    //@{
    int fN[6], fZ[6], fA[6];
    size_t fi;
    //@}
    
    /// Test if entry \c ul in the arrays matches the library reaction
    bool matches(size_t ul, size_t ri);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
