/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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
/** \file frac.h
    \brief Desc
*/
#ifndef O2SCL_FRACT_H
#define O2SCL_FRACT_H

#include <iostream>
#include <string>

#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/table3d.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Jacobian function (not necessarily square) in src/root/jacobian.h
  typedef std::function<
    int(size_t,const boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::matrix<double> &) > nrf_funct;

  /** \brief Desc
   */
  class fract {
    
  public:

    int verbose;

    fract() {
      verbose=1;
    }
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /** \brief Desc

	If the iteration converges to a root, then the <tt>"root"</tt>
	slice in the \ref o2scl::table3d object is set to the index of
	the root (beginning with 1 not 0) and the <tt>"it"</tt> slice
	in the \ref o2scl::table3d object is set to the iteration
	number at which convergence was achieved. If \f$
	x^2+y^2>r_{\mathrm{max}}^2 \f$ then the iteration is presumed
	to have diverged, then <tt>"root"</tt> is set to -1 and
	<tt>"it"</tt> is set to \c kmax. If, after \c kmax iterations,
	the iteration has not converged, then <tt>"root"</tt> is set
	to 0 and <tt>"it"</tt> is set to \c kmax.
     */
    template<class func_t=nrf_funct, class mat_t=ubmatrix,
      class vec_t=ubvector, class vec2_t=std::vector<double>,
      class vec2_size_t=std::vector<size_t>, class fp_t=double>
      int nrf(func_t &f, uniform_grid<fp_t> &gx,
	      uniform_grid<fp_t> &gy, size_t kmax, fp_t rmax,
	      o2scl::table3d &t3d, vec2_t &roots_x,
	      vec2_t &roots_y, vec2_size_t &min_count,
	      vec2_size_t &max_count) {
      
      t3d.clear();
      roots_x.clear();
      roots_y.clear();
      min_count.clear();
      max_count.clear();

      t3d.set_xy("x",gx,"y",gy);
      t3d.new_slice("it");
      t3d.new_slice("root");
      t3d.set_slice_all("it",0.0);
      t3d.set_slice_all("root",0.0);

      vec_t x0(2), x1(2), fx(2);
      mat_t J(2,2);

      //#ifdef O2SCL_OPENMP
      //#pragma omp parallel for
      //#endif
      for(size_t i=0;i<t3d.get_nx();i++) {

	if (verbose>0 && (i+1)%10==0) {
	  std::cout << "nrf progress: " << i+1 << "/" << t3d.get_nx()
		    << std::endl;
	}
	
	for(size_t j=0;j<t3d.get_ny();j++) {
	  
	  x0[0]=gx[i];
	  x0[1]=gy[j];
	  
	  bool found=false;
	  
	  for(size_t k=0;k<kmax && found==false;k++) {

	    // x_{n+1} = x_n - J(x_n)^{-1} f(x_n)
	    f(2,x0,fx,J);
	    // The inverse of the Jacobian
	    fp_t det=J(0,0)*J(1,1)-J(0,1)*J(1,0);
	    fp_t ai=J(1,1)/det;
	    fp_t bi=-J(0,1)/det;
	    fp_t ci=-J(1,0)/det;
	    fp_t di=J(0,0)/det;
	    // The NR step
	    x1[0]=x0[0]-ai*fx[0]-bi*fx[1];
	    x1[1]=x0[1]-ci*fx[0]-di*fx[1];

	    // Test for divergence or convergence
	    fp_t dx=x1[0]-x0[0];
	    fp_t dy=x1[1]-x0[1];
	    if (dx*dx+dy*dy>rmax*rmax) {
	      found=true;
	      t3d.set(i,j,"it",kmax);
	      t3d.set(i,j,"root",-1);
	    } else if (fabs(dx)+fabs(dy)<1.0e-14) {
	      found=true;
	      bool root_found=false;
	      for(size_t m=0;m<roots_x.size();m++) {
		if (fabs(x1[0]-roots_x[m])<1.0e-13 &&
		    fabs(x1[1]-roots_y[m])<1.0e-13) {
		  root_found=true;
		  t3d.set(i,j,"it",k+1);
		  t3d.set(i,j,"root",m+1);
		  if (k+1>max_count[m]) {
		    max_count[m]=k+1;
		  }
		  if (k+1<min_count[m]) {
		    min_count[m]=k+1;
		  }
		}
	      }
	      if (root_found==false) {
		size_t m=roots_x.size();
		roots_x.push_back(x1[0]);
		roots_y.push_back(x1[1]);
		max_count.push_back(k+1);
		min_count.push_back(k+1);
		t3d.set(i,j,"it",k+1);
		t3d.set(i,j,"root",m+1);
	      }
	    }

	    // Proceed to next iteration
	    x0[0]=x1[0];
	    x0[1]=x1[1];
	  }

	  // In the case of neither convergence nor divergence
	  if (found==false) {
	    t3d.set(i,j,"it",kmax);
	  }

	  // Proceed to next point
	}
      }

      return 0;
    }

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
