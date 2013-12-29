/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_NUCLEAR_DIST_H
#define O2SCL_NUCLEAR_DIST_H

#include <iostream>
#include <o2scl/nucleus.h>
#include <o2scl/nuclear_mass.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A distribution of nuclei [abstract base]
      
      The virtual base class for a collection of objects of type \ref
      nucleus . See \ref full_dist and \ref arb_dist
      for implementations of this base class.

      Generally, children need only specify an implementation of
      \ref begin(), \ref end(), \ref size() and \ref next() to be
      functional.

      \future The iterator class has only limited functionality,
      and could be extended to a full STL iterator. Alternatively,
      maybe this class should be replaced with something like
      either vector<nucleus> or map<nucleus>?
   */
  class nuclear_dist {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The function the iterator uses to compute the next nucleus
    virtual nucleus *next(nucleus *np)=0;

#endif
    
  public:

    virtual ~nuclear_dist() {}
    
    /** \brief An iterator for the nuclear distribution
	
	\comment
	The standard usage of this iterator is something of the form:
	\code
	mnmsk_mass mth;
	simple_dist sd(5,6,10,12,&mth);
	for(nuclear_dist::iterator ndi=sd.begin();ndi!=sd.end();ndi++) {
	// do something here for each nucleus
	}
	\endcode
	which would create a list consisting of three isotopes (A=10,
	11, and 12) of boron and three isotopes carbon for a total
	of six nuclei.
	\endcomment
    */
    class iterator {
      
#ifndef DOXYGEN_INTERNAL

    protected:

      /// A pointer to the current nucleus 
      nucleus *np;

      /// A pointer to the distribution
      nuclear_dist *ndp;

#endif

    public:

      /** \brief Create an iterator from the given distribution using the
	  nucleus specified in \c npp.
      */
      iterator(nuclear_dist *ndpp, nucleus *npp) {
	ndp=ndpp;
	np=npp;
      }

      /// Proceed to the next nucleus
      iterator operator++() {
	np=ndp->next(np);
	return iterator(ndp,np);
      }

      /// Proceed to the next nucleus
      iterator operator++(int unused) {
	nucleus *tmp=np;
	np=ndp->next(np);
	return iterator(ndp,tmp);
      }

      /// Pointing at operator
      nucleus *operator->() const {
	return np;
      }

      /// Dereference the iterator
      nucleus &operator*() const {
	return *np;
      }

      /// Give access to the == operator
      friend bool operator==(const nuclear_dist::iterator &i1,
			     const nuclear_dist::iterator &i2);
      
      /// Give access to the != operator
      friend bool operator!=(const nuclear_dist::iterator &i1,
			    const nuclear_dist::iterator &i2);
    };

    /// The beginning of the distribution
    virtual iterator begin()=0;

    /// The end of the distribution
    virtual iterator end()=0;

    /// The number of nuclei in the distribution
    virtual size_t size()=0;
  };

  /// Compare two nuclei
  bool operator==(const nuclear_dist::iterator &i1,
		  const nuclear_dist::iterator &i2);

  /// Compare two nuclei
  bool operator!=(const nuclear_dist::iterator &i1,
		  const nuclear_dist::iterator &i2);
  
  /** \brief Full distribution including all nuclei from a
      discrete mass formula

      For example, to create a collection of all nuclei from the most
      recent (2012) Atomic Mass Evaluation, and then output all the
      nuclei in the collection
      \code
      ame_mass ame;
      full_dist fd(&ame);
      for(nuclear_dist::iterator ndi=fd.begin();ndi!=fd.end();ndi++) {
        cout << ndi->Z << " " << ndi->A << " " << ndi->m << endl;
      }
      \endcode
      
  */
  class full_dist : public nuclear_dist {

  public:
    
    full_dist() {
      list_size=0;
    }

    virtual ~full_dist() {
      if (list_size>0) delete[] list;
    }
    
    /** \brief Create a distribution including all nuclei with atomic
	numbers less than \c maxA from the mass formula \c nm
    */
    full_dist(nuclear_mass &nm, int maxA=400, bool include_neutron=false);
    
    /** \brief Set the distribution to all nuclei with atomic
	numbers less than \c maxA from the mass formula \c nm

	The information for the previous distribution is cleared 
	before a new distribution is set. 
    */
    int set_dist(nuclear_mass &nm, int maxA=400, bool include_neutron=false);
  
    /// The beginning of the distribution
    virtual iterator begin() {
      return iterator(this,&list[0]);
    };

    /// The end of the distribution
    virtual iterator end() {
      return iterator(this,&list[list_size]);
    };

    /// The number of nuclei in the distribution
    virtual size_t size() {
      return list_size;
    };

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Specify how to move to the next nucleus
    virtual nucleus *next(nucleus *np) {
      return np+1;
    }

    /// The distribution of nuclei as an array
    nucleus *list;

    /// The size of \ref list array
    size_t list_size;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
