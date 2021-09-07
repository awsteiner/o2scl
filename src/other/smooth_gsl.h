 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2009-2021, Marco Cammarata and Andrew W. Steiner
  
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
#ifndef SMOOTH_GSL_H
#define SMOOTH_GSL_H

/** \file smooth_gsl.h
    \brief File defining \ref o2scl::smooth_gsl 
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Smooth a GSL vector using GSL bsplines 
      
      \verbatim embed:rst
      
      .. todo::

         In function pair_density():

         - Needs a bit more error checking and more documentation.

         - Future: Generalize to generic vector types. (Does this require
         reworking the GSL linear fitting routines? Doesn't matter now,
         the GSL linear fitting routines are now reworked.)

         - Future: Possibly create a new gsl_bspline class which replaces
         the GSL bspline workspace 

         - Future: Allow user to probe chi squared and the covariance?

      \endverbatim
  */
  class smooth_gsl {

  public:

    smooth_gsl() {
      init_pointers_and_defs();
    }
    
    /// Begin using x-values from vector \c ix 
    smooth_gsl(const gsl_vector *ix) {
      init_pointers_and_defs();
      x=ix;
      init();
    }
    
    ~smooth_gsl() {
      free();
    }
    
    /// Set the number of coefficients
    void set_ncoeff(int incoeffs) { 
      ncoeffs=incoeffs; 
    }

    /// Set order
    void set_order(int order) { 
      norder=order;  
    }

    /// Set parameters
    void set_pars(int incoeffs, int order) { 
      ncoeffs=incoeffs;
      norder=order;
    }

    /** \brief Set the x-values

	\comment 
	This is just a reminder to me that this function can't
	be virtual since it is called in a constructor.
	\endcomment
    */
    void set_x(const gsl_vector *ix) {
      x=ix;
      init();
    }
    
    /** \brief Smooth data in \c y with errors \c e returning result \c ys
     */
    int smooth_data(const gsl_vector *y, const gsl_vector *e, gsl_vector *ys);
    
    /** \brief Smooth data in \c y returning result \c ys
     */
    int smooth_data(const gsl_vector *y, gsl_vector *ys);
    
  protected:

    /// Number of free coefficients for spline
    size_t ncoeffs;

    /// Order of spline to be used (4=cubic)
    size_t norder;

    /// internally calculated, number of "segment" to split the data into
    size_t nbreak;

    /// True of the x values have been set
    bool x_set;

    /// Spline workspace
    gsl_bspline_workspace *bw;

    /// Spline temporary vector
    gsl_vector *B;

    /// Parameters of linear fit, y=X*c 
    gsl_vector *c;

    /// Linear fit workspace
    gsl_multifit_linear_workspace *mw;

    /// Workspace for spline fitting
    gsl_matrix *X;

    /// Values of the independent variable
    const gsl_vector *x;

    /// Covariance matrix
    gsl_matrix *cov;

    /// Construct un-weighted fit
    int fit(const gsl_vector *y);
    
    /// Construct weighted fit
    int fit_errors(const gsl_vector *y, const gsl_vector *e);

    /// calculate smoothed curve value for a certain xi
    double calc_for_x(double xi);
    
    /** \brief Allocate memory and initialize splines
	
	\comment 
	This is just a reminder to me that this function can't
	be virtual since it is called in a constructor.
	\endcomment
    */
    int init();
    
    /// Set default values and zero pointers
    void init_pointers_and_defs();

    /// Free memory
    int free();

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
