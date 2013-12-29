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
#ifndef O2SCL_ARB_DIST_H
#define O2SCL_ARB_DIST_H

#include <iostream>
#include <o2scl/err_hnd.h>
#include <o2scl/nucleus.h>
#include <o2scl/nuclear_mass.h>
#include <o2scl/nuclear_dist.h>
#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Arbitrary distribution of nuclei chosen either with a 
      function or a list
  */
  class arb_dist : public nuclear_dist {

  public:

    /// Create an empty distribution
    arb_dist() {
      list_size=0;
    }

    virtual ~arb_dist() {
      if (list_size>0) delete[] list;
    }
    
    /** \brief Create a distribution directly from a list of values
	for Z and N
    */
    template<class int_vec_t, class int_vec2_t>
      void set_dist(nuclear_mass &nm, size_t sz, 
		    int_vec_t &Zvec, int_vec2_t &Nvec) {
      
      // Delete old list if present
      if (list_size>0) delete[] list;

      // Allocate memory
      list_size=sz;
      list=new nucleus[list_size];

      // Copy selected nuclei over
      for(size_t i=0;i<list_size;i++) {
	nm.get_nucleus(Zvec[i],Nvec[i],list[i]);
      }

      return;
    }
    
    /** \brief Add the nuclei in a distribution to the current one

	\note This is relatively slow because it has to count the new
	distribution first, then create a new pointer and copy all of
	the old nuclei over.
     */
    void add_dist(nuclear_dist &nd) {
      
      // First, count nuclei in the new distribution which
      // are not already included in the old one
      size_t new_nuclei=0;
      for(nuclear_dist::iterator ndi=nd.begin();ndi!=nd.end();ndi++) {
	bool found=false;
	for(nuclear_dist::iterator ndi2=this->begin();ndi2!=this->end();
	    ndi2++) {
	  if (ndi->Z==ndi2->Z && ndi->N == ndi2->N) {
	    found=true;
	  }
	}
	if (!found) new_nuclei++;
      }
      if (new_nuclei==0) return;

      // Make space for new distribution
      nucleus *new_list=new nucleus[list_size+new_nuclei];

      // Copy old nuclei
      for(size_t i=0;i<list_size;i++) {
	new_list[i]=list[i];
      }
      // Copy new nuclei not in the old distribution
      size_t j=list_size;
      for(nuclear_dist::iterator ndi=nd.begin();ndi!=nd.end();ndi++) {
	bool found=false;
	for(nuclear_dist::iterator ndi2=this->begin();ndi2!=this->end();
	    ndi2++) {
	  if (ndi->Z==ndi2->Z && ndi->N == ndi2->N) {
	    found=true;
	  }
	}
	if (!found) {
	  new_list[j]=*ndi;
	  j++;
	}
      }

      // Now switch distributions
      nucleus *old_list=list;
      list=new_list;
      list_size+=new_nuclei;
      delete[] old_list;
      
      return;
    }
    
    /** \brief Set the distribution to include nuclei from masses 
	in \c nm constrained by function \c expr

	The string \c expr can contain any function of \c Z and \c N.
	All nuclei for which the function evaluates (using \ref
	FunctionParser) to a number greater than 0.5 are included in
	the distribution. An example is <tt>"Z>=10 & Z<=20 & N>=10 & N
	<=20"</tt> for a square distribution, or "(Z%2=0) & (N%2=0)"
	for all even-even nuclei.
    */
    void set_dist(nuclear_mass &nm, std::string expr, int maxA=400,
		  bool include_neutron=false) {
  
      if (list_size>0) delete[] list;
      
      nucleus n;
      full_dist fd(nm,maxA,include_neutron);
      double vals[2];

      // Parse the formula
      int ret=fp.Parse(expr,"Z,N");
      if (ret!=-1) {
	O2SCL_ERR("Failed to parse in arb_dist::arb_dist().",exc_einval);
      }

      // Count up nuclei
      list_size=0;
      for(nuclear_dist::iterator ndi=fd.begin();ndi!=fd.end();ndi++) {
	vals[0]=ndi->Z;
	vals[1]=ndi->N;
	if (fp.Eval(vals)>0.5) {
	  list_size++;
	}
      }

      // Allocate memory
      list=new nucleus[list_size];

      // Copy selected nuclei over
      int ix=0;
      for(nuclear_dist::iterator ndi=fd.begin();ndi!=fd.end();ndi++) {
	vals[0]=ndi->Z;
	vals[1]=ndi->N;
	if (fp.Eval(vals)>0.5) {
	  nm.get_nucleus(ndi->Z,ndi->N,list[ix]);
	  ix++;
	}
      }

      return;
    }
  
    /// The beginning of the distribution
    virtual iterator begin() {
      return iterator(this,&list[0]);
    };

    /// The end of the distribution
    virtual iterator end() {
      return iterator(this,&list[list_size]);
    };

    /// A location in the middle of the distribution
    virtual iterator index(size_t i) {
      return iterator(this,&list[i]);
    }

    /// The number of nuclei in the distribution
    virtual size_t size() {
      return list_size;
    };

#ifndef DOXYGEN_NO_O2NS

  protected:

    /// The function parser
    FunctionParser fp;

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
