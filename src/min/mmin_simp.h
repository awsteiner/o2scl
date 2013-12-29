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
/* multimin/simplex.c
 * 
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_GSL_MMIN_SIMP_H
#define O2SCL_GSL_MMIN_SIMP_H

#include <gsl/gsl_multimin.h>
#include <o2scl/mmin.h>
#include <o2scl/cblas.h>
#include <o2scl/vec_stats.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimization by the Simplex method (GSL)

      This class mins a function using Nelder and Mead's Simplex
      algorithm. A simplex in a N-dimensional space is defined as a
      set of N+1 points which describe an N-dimensional volume
      surrounding the minimum. The algorithm proceeds by shifting the
      simplex points until the simplex is sufficiently small and thus
      the minimum is known with sufficient accuracy.

      For a slightly improved method, see \ref mmin_simp2 .

      This class has a high-level interface using mmin(),
      mmin_twovec() or mmin_simplex() which automatically performs the
      memory allocation and minimization, or a GSL-like interface
      using allocate(), free(), interate() and set() or set_simplex().

      The simplex can be completely specified by the user (see
      mmin_simplex() and set_simplex()). Alternatively, the simplex is
      automatically specified given initial guess \f$ x_j \f$ and a
      step size vector \f$ s_k \f$ for \f$ 0\leq k<n_s \f$. The
      simplex \f$ p_{ij} \f$ with \f$ 0\leq i\leq n \f$ and \f$ 0\leq
      j<n \f$ is chosen with \f$ p_{0j} = x_j \f$ and
      \f{eqnarray*}
      p_{i+1,j} &=& x_j\quad\mathrm{for}\quad i\neq j \\
      p_{i+1,j} &=& x_j+s_{j~\mathrm{mod}~n_s}\quad\mathrm{for}\quad i=j
      \f}
      for \f$ 0<i<n \f$. The step size vector \f$ s \f$ is set by the
      set_step() member function. The prsence of \f$ \mathrm{mod} \f$
      in the recipe above just indicates that elements of the step
      size vector are automatically re-used if there are less step
      sizes than dimensions in the minimization problem.
      
      \note It is important that the initial simplex contains
      sufficient variation in every direction of the parameter space
      over which one is minimizing. For example, if all three points
      in a simplex for minimizing over a two-dimensional space contain
      nearly the same value for the second parameter, then the
      minimizer may only min the function with respect to the
      first parameter. 

      Default template arguments
      - \c func_t - \ref multi_funct\<\>
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>

      Based on \ref Nelder65 .
  */
#ifndef O2SCL_CPP11
  template<class func_t=multi_funct<>,
    class vec_t=boost::numeric::ublas::vector<double> > class mmin_simp :
    public mmin_base<func_t,func_t,vec_t>
#else
    template<class func_t=multi_funct11,
    class vec_t=boost::numeric::ublas::vector<double> > class mmin_simp :
    public mmin_base<func_t,grad_funct11,vec_t>
#endif
    {
    
#ifndef DOXYGEN_INTERNAL
    
    public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
    protected:
  
    /// An array of n+1 vectors containing the simplex
    vec_t *x1;
    /** \brief The n+1 function values at the simplex points
      
	\comment
	This can't be the template type because it has to be 
	a different size (n+1) rather than n.
	\endcomment
    */
    ubvector y1;
    /// Workspace vector 1
    vec_t ws1;
    /// Workspace vector 2
    vec_t ws2;
    /// Workspace vector 3
    vec_t ws3;
    
    /// Compute the center of the simplex and store in \c mp
    int nmsimplex_calc_center(vec_t &mp) {
      size_t i, j;
      double val;
      
      for (j = 0; j < dim; j++) {
	val = 0.0;
	for (i = 0; i < dim+1; i++) {
	  val += x1[i][j];
	}
	val /= dim+1;
	mp[j]=val;
      }
      
      return 0;
    }
    
    /** \brief Compute the size of the simplex
	
	Calculates simplex size as average sum of length of vectors
	from simplex center to corner points:
	
	\f$ (1/n) \sum || y - y_{\mathrm{middlepoint}} || \f$
    */
    double nmsimplex_size() {
      
      double ss=0.0;
      nmsimplex_calc_center(ws2);
      
      for(size_t i=0;i<dim+1;i++) {
	for(size_t j=0;j<dim;j++) ws1[j]=x1[i][j];
	o2scl_cblas::daxpy(-1.0,dim,ws2,ws1);
	ss += o2scl_cblas::dnrm2(dim,ws1);
      }
      
      return ss/(double)(dim+1);
    }

    /** \brief Move a corner of a simplex

	Moves a simplex corner scaled by coeff (negative value
	represents mirroring by the middle point of the "other" corner
	points) and gives new corner in xc and function value at xc 
	in \c newval.
    */
    virtual int move_corner_err(const double coeff, size_t corner, 
				vec_t &xc, func_t &f, 
				size_t nvar, double &newval) {
      
      size_t i,j;
      double mp;
      
      for (j=0;j<dim;j++) {
	mp=0.0;
	for (i=0;i<dim+1;i++) {
	  if (i!=corner) {
	    mp+=x1[i][j];
	  }
	}
	mp/=(double)dim;
	xc[j]=mp-coeff*(mp-x1[corner][j]);
      }
	
      newval=f(nvar,xc);

      return 0;
    }
      
    /** \brief Contract the simplex towards the best point
	  
	Function contracts the simplex in respect to best valued
	corner. All corners besides the best corner are moved.
	
	The vector, \c xc, is used as work space.
    */
    virtual int contract_by_best(size_t best, vec_t &xc,
				 func_t &f, size_t nvar) {
      
      size_t i,j, it;
      double newval;
      bool failed;

      for (i=0;i<dim+1;i++) {
	
	if (i!=best) {
	  
	  for (j=0;j<dim;j++) {
	    x1[i][j]=0.5*(x1[i][j]+x1[best][j]);
	  }
	  y1[i]=f(nvar,x1[i]);
	  if (!o2scl::is_finite(y1[i])) {
	    O2SCL_ERR2_RET("Function not finite in ",
			   "mmin_simp::contract_by_best().",
			   exc_ebadfunc);
	  }

	  if (avoid_nonzero==true) {
	    std::cout << "Found problem in contract." << std::endl;

	    // copy the old point to ws3
	    for (j=0;j<dim;j++) {
	      ws3[j]=2.0*x1[i][j]-x1[best][j];
	    }
	    
	    // Try (21*best+20*xold)/41, (22*best+19*xold)/41, ...,
	    // (40*best+xold)/41 to see if we can find a contraction
	    // that works
	    for(it=0;it<20;it++) {
	      for (j=0;j<dim;j++) {
		x1[i][j]=((20-it)*ws3[j]+(it+21)*x1[best][j])/41.0;
	      }
	      y1[i]=f(nvar,x1[i]);
	      //std::cout << "it, x: " << it << " "
	      //		<< x << std::endl;
	    }
	    /*
	      if (ret!=0) {
	      O2SCL_CONV2_RET("Failed to find suitable contraction ",
	      "in mmin_simp::contract_by_best().",
	      exc_efailed,this->err_nonconv);
	      }
	    */
	  }
	}
      }
	
      return success;
    }

    /// Number of variables to be mind over
    size_t dim;

    /// Function
    func_t *func;

    /// True if set() has been called
    bool set_called;

    /// Vector of step sizes
    ubvector step_vec;

    /** \brief If true, try to automatically avoid regions where the 
	function returns a non-zero value (default false)
	  
	\note This option doesn't work yet, so I've made the variable
	protected to prevent the user from changing it.
    */
    bool avoid_nonzero;

#endif

    public:

    mmin_simp() {
      dim=0;
      print_simplex=0;
      step_vec.resize(1);
      step_vec[0]=1.0;
      avoid_nonzero=false;
    }
    
    virtual ~mmin_simp() {
      free();
      step_vec.clear();
    }

    /// Set the step sizes for each independent variable
    template<class vec2_t> int set_step(size_t nv, vec2_t &step) {
      if (nv>0) {
	step_vec.clear();
	step_vec.resize(nv);
	for(size_t i=0;i<nv;i++) step_vec[i]=step[i];
      }
      return 0;
    }

    /// Size of current simplex computed by iterate()
    double size;

    /// Present minimum vector computed by iterate()
    vec_t x;

    /// Function value at minimum computed by iterate()
    double fval;

    /** \brief Print simplex information in print_iter() (default 0)
	
	If this is 1 and \ref verbose is greater than 0, 
	then print_iter() will print the function values at all 
	the simplex points. If this is 2 and \ref verbose is
	greater than 0, then print_iter() will print the 
	simplex coordinates in addition to the function values.
    */
    int print_simplex;

    /** \brief Calculate the minimum \c min of \c func w.r.t the
	array \c x of size \c nvar.
    */
    virtual int mmin(size_t nn, vec_t &xx, double &fmin, 
		     func_t &ufunc) {
	
      if (nn==0) {
	O2SCL_ERR2_RET("Tried to min over zero variables ",
		       " in mmin_simp::mmin().",exc_einval);
      }

      int ret=0,status,iter=0;
      
      allocate(nn);

      vec_t ss(nn);
      for (size_t is=0;is<nn;is++) ss[is]=step_vec[is % step_vec.size()];
      ret=set(ufunc,nn,xx,ss);

      if(ret!=0) {
	free();
	return ret;
      }
  
      do {
	iter++;
	  
	status=iterate();
	if(status) break;
	
	if(this->verbose>0) {
	  print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
		     "mmin_simp");
	}
    
	status=gsl_multimin_test_size(size,this->tol_abs);
	  
      } while(status == GSL_CONTINUE && iter<this->ntrial);
	
      for (size_t i=0;i<nn;i++) xx[i]=x[i];
      fmin=fval;
  
      free();
      this->last_ntrial=iter;
      
      if(iter>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	  itos(this->ntrial)+") in mmin_simp2::mmin().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }
    
    /** \brief Calculate the minimum \c min of \c func w.r.t the
	array \c x of size \c nvar, using \c xx and \c xx2 to specify
	the simplex
    */
    virtual int mmin_twovec(size_t nn, vec_t &xx, vec_t &xx2, double &fmin, 
			    func_t &ufunc) {
      
      int ret=0,i,status,iter=0;
      
      allocate(nn);

      vec_t ss(nn);

      for (size_t is=0;is<nn;is++) ss[is]=xx2[is]-xx[is];
      ret=set(ufunc,nn,xx,ss);
      
      if(ret!=0) {
	free();
	return ret;
      }
  
      do {
	iter++;
	  
	status=iterate();
	if(status) break;

	if(this->verbose>0) {
	  print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
		     "mmin_simp");
	}
	
	status=gsl_multimin_test_size(size,this->tol_abs);
	  
      } while(status == GSL_CONTINUE && iter<this->ntrial);
	
      for (i=0;i<((int)nn);i++) xx[i]=x[i];
      fmin=fval;
  
      free();
      this->last_ntrial=iter;

      if(iter>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	  itos(this->ntrial)+") in mmin_simp2::mmin_twovec().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }
    
    /** \brief Calculate the minimum \c min of \c func w.r.t the
	array \c x of size \c nvar, given an initial simplex
    */
    template<class mat_t> 
    int mmin_simplex(size_t nn, mat_t &sx, double &fmin, 
		     func_t &ufunc) {
      
      int ret=0,i,status,iter=0;
      
      allocate(nn);
	
      ret=set_simplex(ufunc,sx);
      if(ret!=0) {
	free();
	return ret;
      }
  
      do {
	iter++;
	  
	status=iterate();
	if(status) break;

	if(this->verbose>0) {
	  print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
		     "mmin_simp");
	}
    
	status=gsl_multimin_test_size(size,this->tol_abs);
	  
      } while(status == GSL_CONTINUE && iter<this->ntrial);
	
      for (i=0;i<((int)nn);i++) sx(0,i)=x[i];
      fmin=fval;
  
      free();
      this->last_ntrial=iter;

      if(iter>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	  itos(this->ntrial)+") in mmin_simp::mmin_simplex().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }
    
    /// Allocate the memory
    virtual int allocate(size_t n) {
      int status;
      if(dim!=0) free();
      set_called=false;
      dim=n;

      x.resize(n);
      x1=new vec_t[n+1];
      for(size_t i=0;i<n+1;i++) x1[i].resize(n);
      y1.resize(n+1);
      ws1.resize(n);
      ws2.resize(n);
      ws3.resize(n);

      return success;
    }
    
    /// Free the allocated memory
    virtual int free() {

      if (dim>0) {
	for(size_t i=0;i<dim+1;i++) x1[i].clear();
	delete[] x1;
	y1.clear();
	ws1.clear();
	ws2.clear();
	ws3.clear();
	x.clear();
      }

      dim=0;

      return 0;
    }

    /// Set the function and initial guess
    virtual int set(func_t &ufunc, size_t n, vec_t &ax, 
		    vec_t &step_size) {
      size_t i;
      
      if(dim!=n) allocate(n);
      func=&ufunc;

      // Copy initial guess to x
      for (i=0;i<dim;i++) x[i]=ax[i];
      
      // first point is the original x0 
      
      y1[0]=ufunc(dim,ax);
      if (!o2scl::is_finite(y1[0])) {
	O2SCL_ERR2_RET("Function not finite in ",
		       "mmin_simp::set().",exc_ebadfunc);
      }
      for(i=0;i<dim;i++) x1[0][i]=ax[i];
  
      /* following points are initialized to x0+step_size */
      
      for (i=1;i<dim+1;i++) {
	for(size_t j=0;j<dim;j++) x1[i][j]=x[j];
	x1[i][i-1]=x1[i][i-1]+step_size[i-1];
	y1[i]=ufunc(dim,x1[i]);
      }
 
      /* Initialize simplex size */
      
      size=nmsimplex_size();
	
      set_called=true;

      return success;
    }

    /// Set the function and initial simplex
    template<class mat_t> 
    int set_simplex(func_t &ufunc, mat_t &sx) {

      if(dim==0) {
	O2SCL_ERR_RET("Memory not allocated in mmin_simp::set().",
		      exc_ebadlen);
      }
      
      func=&ufunc;

      for(size_t i=0;i<dim+1;i++) {
	for(size_t j=0;j<dim;j++) {
	  x1[i][j]=sx(i,j);
	}
	y1[i]=ufunc(dim,x1[i]);
	if (!o2scl::is_finite(y1[i])) {
	  O2SCL_ERR2_RET("Function not finite in ",
			 "mmin_simp::set_simplex().",exc_ebadfunc);
	}
      }
	
      /* Initialize simplex size */
	
      size=nmsimplex_size();
	
      set_called=true;
 
      return success;
    }

    /// Perform an iteration
    virtual int iterate() {

      size_t n=dim+1;
      size_t i;
      size_t lo, hi, s_hi;
      double dhi,ds_hi,dlo;
      int status;
      double val,val2;
 
      /* get index of highest,second highest and lowest point */

      lo=hi=0;
      dlo=dhi=y1[0];
      s_hi=1;
      ds_hi=y1[1];
      
      for (i=1;i<n;i++) {
	val=y1[i];
	if(val<dlo) {
	  dlo=val;
	  lo=i;
	} else if (val > dhi) {
	  ds_hi=dhi;
	  s_hi=hi;
	  dhi=val;
	  hi=i;
	} else if(val > ds_hi) {
	  ds_hi=val;
	  s_hi=i;
	}
      }

      //std::cout << y1 << std::endl;
 
      /* reflect the highest value */

      int ret1=move_corner_err(-1.0,hi,ws1,*func,dim,val);

      if (avoid_nonzero && ret1!=0) {
	std::cout << "Found problem move1: " << std::endl;
	for (size_t it=0;it<20 && ret1!=0;it++) {
	  ret1=move_corner_err(-1.0+((double)it)/10.0,hi,ws1,
			       *func,dim,val);
	  std::cout << "it,ret: " << it << " " << ret1 << std::endl;
	}
	if (ret1!=0) {
	  O2SCL_ERR2_RET("Failed to move corner (1) in ",
			 "mmin_simp::iterate().",exc_efailed);
	}
      }
      
      if (o2scl::is_finite(val) && val<y1[lo]) {

	//std::cout << "1" << std::endl;

	/* reflected point becomes lowest point,try expansion */

	int ret2=move_corner_err(-2.0,hi,ws2,*func,dim,val2);
	if (avoid_nonzero && ret2!=0) {
	  std::cout << "Found problem move2: " << std::endl;
	  for (size_t it=0;it<20 && ret2!=0;it++) {
	    ret2=move_corner_err(-2.0+((double)it)/6.7,hi,ws2,
				 *func,dim,val2);
	    std::cout << "it,ret: " << it << " " << ret2 << std::endl;
	  }
	  if (ret2!=0) {
	    O2SCL_ERR2_RET("Failed to move corner (2) in ",
			   "mmin_simp::iterate().",exc_efailed);
	  }
	}

	if (o2scl::is_finite(val2) && val2<y1[lo]) {
	  for(i=0;i<dim;i++) x1[hi][i]=ws2[i];
	  y1[hi]=val2;
	} else {
	  for(i=0;i<dim;i++) x1[hi][i]=ws1[i];
	  y1[hi]=val;
	}
	  
      } else if (!o2scl::is_finite(val) || val > y1[s_hi]) {
	
	//std::cout << "2: " << hi << " " << s_hi << std::endl;

	/* reflection does not improve things enough */
	
	if (o2scl::is_finite(val) && val <= y1[hi]) {
	    
	  /* if trial point is better than highest point,replace
	     highest point */
	  
	  for(i=0;i<dim;i++) x1[hi][i]=ws1[i];
	  y1[hi]=val;
	}
      
	/* try one dimensional contraction */
	
	int ret3=move_corner_err(0.5,hi,ws2,*func,dim,val2);
	if (avoid_nonzero && ret3!=0) {
	  std::cout << "Found problem move3: " << std::endl;
	  for (size_t it=0;it<20 && ret3!=0;it++) {
	    ret3=move_corner_err(0.025*((double)(19-it)),hi,ws2,
				 *func,dim,val2);
	    std::cout << "it,ret: " << it << " " << ret3 << std::endl;
	  }
	  if (ret3!=0) {
	    O2SCL_ERR2_RET("Failed to move corner (2) in ",
			   "mmin_simp::iterate().",exc_efailed);
	  }
	}
	  
	if (o2scl::is_finite(val2) && val2 <= y1[hi]) {
	  //std::cout << "3" << std::endl;
	    
	  for(i=0;i<dim;i++) x1[hi][i]=ws2[i];
	  y1[hi]=val2;
	    
	} else {
	  //std::cout << "4" << std::endl;
	    
	  /* contract the whole simplex in respect to the best point */
	  status=contract_by_best(lo,ws1,*func,dim);
	  if(status != 0) {
	    O2SCL_ERR("Function contract_by_best() failed in iterate().",
		      exc_efailed);
	  }
	  
	}

      } else {

	/* trial point is better than second highest point.
	   Replace highest point by it */
      
	for(i=0;i<dim;i++) x1[hi][i]=ws1[i];
	y1[hi]=val;
      }
  
      //std::cout << y1 << std::endl;

      /* return lowest point of simplex as x */
      
      vector_min(dim+1,y1,lo,val);
      for(i=0;i<dim;i++) x[i]=x1[lo][i];
      fval=y1[lo];
      
      /* Update simplex size */
  
      size=nmsimplex_size();

      return success;
    }
      
    /** \brief Print out iteration information.
	  
	Depending on the value of the variable verbose, this prints out
	the iteration information. If verbose=0, then no information is
	printed, while if verbose>1, then after each iteration, the
	present values of x and y are output to std::cout along with the
	iteration number. If verbose>=2 then each iteration waits for a
	character.
    */
    virtual int print_iter(size_t nv, vec_t &xx, vec_t *simp,
			   double y, int iter, double value, 
			   double limit, std::string comment) {

      if (this->verbose<=0) return 0;
      
      size_t i;
      char ch;
      
      (*this->outs) << comment << " Iteration: " << iter << std::endl;
      (*this->outs) << "x: ";
      for(i=0;i<nv;i++) (*this->outs) << xx[i] << " ";
      (*this->outs) << std::endl;
      if (print_simplex>0) {
	if (print_simplex==2) {
	  (*this->outs) << "Simplex Values:" << std::endl;
	  for(i=0;i<nv+1;i++) {
	    (*this->outs) << i << ": ";
	    for(size_t j=0;j<nv;j++) {
	      (*this->outs) << simp[i][j] << " ";
	    }
	    (*this->outs) << ": " << y1[i] << std::endl;
	  }
	} else {
	  (*this->outs) << "Simplex Values:" << std::endl;
	  for(i=0;i<nv+1;i++) {
	    (*this->outs) << y1[i] << " ";
	  }
	  (*this->outs) << std::endl;
	}
      }
      (*this->outs) << "y: " << y << " Val: " << value << " Lim: " 
      << limit << std::endl;
      if (this->verbose>1) {
	(*this->outs) << "Press a key and type enter to continue. ";
	(*this->ins) >> ch;
      }
	
      return 0;
    }

    /// Return string denoting type("mmin_simp")
    virtual const char *type() { return "mmin_simp";}

    };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
