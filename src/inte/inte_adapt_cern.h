/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
/** \file inte_adapt_cern.h
    \brief File defining \ref o2scl::inte_adapt_cern
*/
#ifndef O2SCL_CERN_ADAPT_H
#define O2SCL_CERN_ADAPT_H

#include <o2scl/inte.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/string_conv.h>
 
#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Adaptive integration (CERNLIB)
    
      Uses a base integration object (default is \ref
      inte_gauss56_cern) to perform adaptive integration by
      automatically subdividing the integration interval. At each
      step, the interval with the largest absolute uncertainty is
      divided in half. The routine succeeds if the absolute tolerance
      is less than \ref tol_abs or if the relative tolerance is less
      than \ref tol_rel, i.e.
      \f[
      \mathrm{err}\leq\mathrm{tol\_abs}~\mathrm{or}~
      \mathrm{err}\leq\mathrm{tol\_rel}\cdot|I|
      \f]
      where \f$ I \f$ is the current estimate for the integral and \f$
      \mathrm{err} \f$ is the current estimate for the uncertainty. If
      the number of subdivisions exceeds the template parameter \c
      nsub, the error handler is called, since the integration may not
      have been successful. The number of subdivisions used in the
      last integration can be obtained from get_nsubdivisions().

      The template parameter \c nsub, is the maximum number of
      subdivisions. It is automatically set to 100 in the original
      CERNLIB routine, and defaults to 100 here. The default base
      integration object is of type \ref inte_gauss56_cern. This is the
      CERNLIB default, but can be modified by calling set_inte().
      
      This class is based on the CERNLIB routines RADAPT and
      DADAPT which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d102/top.html
      
      \future 
      - Allow user to set the initial subdivisions?
      - It might be interesting to directly compare the performance
      of this class to \ref o2scl::inte_qag_gsl .
      - There is a fixme entry in the code which could be resolved.
      - Output the point where most subdividing was required?
  */
  template<class func_t=funct, size_t nsub=100> 
    class inte_adapt_cern : public inte<func_t> {

#ifndef DOXYGEN_INTERNAL

      protected:

      /// Lower end of subdivision
      double xlo[nsub];

      /// High end of subdivision
      double xhi[nsub];

      /// Value of integral for subdivision
      double tval[nsub];

      /// Squared error for subdivision
      double ters[nsub];
      
      /// Previous number of subdivisions
      int prev_subdiv;

      /// The base integration object
      inte<func_t> *it;
      
#endif
      
      public:
  
      inte_adapt_cern() {
	nsubdiv=1;
	prev_subdiv=0;

	it=&def_inte;
      }

      /// \name Basic usage
      //@{
      /** \brief Integrate function \c func from \c a to \c b
	  giving result \c res and error \c err
      */
      virtual int integ_err(func_t &func, double a, double b,
			    double &res, double &err) {
	
	double tvals=0.0, terss, xlob, xhib, yhib=0.0, te, root=0.0;
	int i, nsubdivd;

	if (nsubdiv==0) {
	  if (prev_subdiv==0) {
	    // If the previous binning was requested, but
	    // there is no previous binning stored, then
	    // just shift to automatic binning
	    nsubdivd=1;
	  } else {
	    tvals=0.0;
	    terss=0.0;
	    for(i=0;i<prev_subdiv;i++) {
	      it->integ_err(func,xlo[i],xhi[i],tval[i],te);
	      ters[i]=te*te;
	      tvals+=tval[i];
	      terss+=ters[i];
	    }
	    err=sqrt(2.0*terss);
	    res=tvals;
	    return 0;
	  }
	}
  
	// Ensure we're not asking for too many divisions
	if (nsub<nsubdiv) {
	  nsubdivd=nsub;
	} else {
	  nsubdivd=nsubdiv;
	}

	// Compute the initial set of intervals and integral values
	xhib=a;
	double bin=(b-a)/((double)nsubdivd);
	for(i=0;i<nsubdivd;i++) {
	  xlo[i]=xhib;
	  xlob=xlo[i];
	  xhi[i]=xhib+bin;
	  if (i==nsubdivd-1) xhi[i]=b;
	  xhib=xhi[i];
	  it->integ_err(func,xlob,xhib,tval[i],te);
	  ters[i]=te*te;
	}
	prev_subdiv=nsubdivd;

	for(size_t iter=1;iter<=nsub;iter++) {

	  // Compute the total value of the integrand
	  // and the squared uncertainty
	  tvals=tval[0];
	  terss=ters[0];
	  for(i=1;i<prev_subdiv;i++) {
	    tvals+=tval[i];
	    terss+=ters[i];
	  }
	  
	  // Output iteration information
	  if (this->verbose>0) {
	    std::cout << "inte_adapt_cern Iter: " << iter;
	    std::cout.setf(std::ios::showpos);
	    std::cout << " Res: " << tvals;
	    std::cout.unsetf(std::ios::showpos);
	    std::cout << " Err: " << sqrt(2.0*terss);
	    if (this->tol_abs>this->tol_rel*fabs(tvals)) {
	      std::cout << " Tol: " << this->tol_abs << std::endl;
	    } else {
	      std::cout << " Tol: " << this->tol_rel*fabs(tvals) << std::endl;
	    }
	    if (this->verbose>1) {
	      char ch;
	      std::cout << "Press a key and type enter to continue. " ;
	      std::cin >> ch;
	    }
	  }

	  // See if we're finished
	  root=sqrt(2.0*terss);
	  if (root<=this->tol_abs || root<=this->tol_rel*fabs(tvals)) {
	    res=tvals;
	    err=root;
	    this->last_iter=iter;
	    return 0;
	  }

	  // Test if we've run out of intervals
	  if (prev_subdiv==nsub) {
	    res=tvals;
	    err=root;
	    this->last_iter=iter;
	    std::string s="Reached maximum number ("+itos(nsub)+
	      ") of subdivisions in inte_adapt_cern::integ_err().";
	    O2SCL_CONV_RET(s.c_str(),exc_etable,this->err_nonconv);
	  }

	  // Find the subdivision with the largest error
	  double bige=ters[0];
	  int ibig=0;
	  for(i=1;i<prev_subdiv;i++) {
	    if (ters[i]>bige) {
	      bige=ters[i];
	      ibig=i;
	    }
	  }

	  // Subdivide that subdivision further
	  xhi[prev_subdiv]=xhi[ibig];
	  double xnew=(xlo[ibig]+xhi[ibig])/2.0;
	  xhi[ibig]=xnew;
	  xlo[prev_subdiv]=xnew;
	  it->integ_err(func,xlo[ibig],xhi[ibig],tval[ibig],te);
	  ters[ibig]=te*te;
	  it->integ_err(func,xlo[prev_subdiv],
			xhi[prev_subdiv],tval[prev_subdiv],te);
	  ters[prev_subdiv]=te*te;
	  prev_subdiv++;

	}

	// FIXME: Should we set an error here, or does this
	// only happen if we happen to need exactly nsub
	// intervals?
	res=tvals;
	err=root;
	return 0;
      }
      //@}

      /// \name Integration object
      //@{
      /// Set the base integration object to use
      int set_inte(inte<func_t> &i) {
	it=&i;
	return 0;
      }
      
      /// Default integration object
      inte_gauss56_cern<func_t> def_inte;
      //@}
      
      /// \name Subdivisions
      //@{
      /** \brief Number of subdivisions 

	  The options are
	  - 0: Use previous binning and do not subdivide further \n
	  - 1: Automatic - adapt until tolerance is attained (default) \n
	  - n: (n>1) split first in n equal subdivisions, then adapt
	  until tolerance is obtained.
      */
      size_t nsubdiv;

      /// Return the number of subdivisions used in the last integration
      size_t get_nsubdivisions() {
	return prev_subdiv;
      }

      /// Return the ith subdivision
      int get_ith_subdivision(size_t i, double &xlow, double &xhigh, 
			      double &value, double &errsq) {
	if (i<prev_subdiv) {
	  xlow=xlo[i];
	  xhigh=xhi[i];
	  value=tval[i];
	  errsq=ters[i];
	}
	return 0;
      }
      
      /// Return all of the subdivisions
      template<class vec_t> 
      int get_subdivisions(vec_t &xlow, vec_t &xhigh, vec_t &value, 
			   vec_t &errsq) {

	for(int i=0;i<prev_subdiv;i++) {
	  xlow[i]=xlo[i];
	  xhigh[i]=xhi[i];
	  value[i]=tval[i];
	  errsq[i]=ters[i];
	}
	return 0;
      }
      //@}

    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
