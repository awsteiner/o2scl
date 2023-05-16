/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2023, Andrew W. Steiner

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
/* multimin/simplex2.c
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
#ifndef O2SCL_MMIN_SIMP2_H
#define O2SCL_MMIN_SIMP2_H

/** \file mmin_simp2.h
    \brief File defining \ref o2scl::mmin_simp2
*/
#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include <gsl/gsl_multimin.h>

#include <o2scl/mmin.h>
#include <o2scl/cblas.h>
#include <o2scl/vec_stats.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimization by the Simplex method (v2) (GSL)

      This class mins a function using Nelder and Mead's Simplex
      algorithm. A simplex in a N-dimensional space is defined as a
      set of N+1 points which describe an N-dimensional volume
      surrounding the minimum. The algorithm proceeds by shifting the
      simplex points until the simplex is sufficiently small and thus
      the minimum is known with sufficient accuracy.

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
      minimizer may only minimize the function with respect to the
      first parameter.

      \note The algorithm used to estimate the simplex size does not
      work well any of the parameters in the minimization problem has
      a scale which is not close to 1.

      \verbatim embed:rst
      See an example for the usage of this class in 
      :ref:`Multidimensional minimizer example`.
      \endverbatim

      Default template arguments
      - \c param_t - no default
      - \c func_t - \ref multi_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>

      \verbatim embed:rst
      Based on [Nelder65]_.
      \endverbatim

      A variable <tt>count</tt> originally defined in the GSL simplex
      state is not present here, because it was unused.

      \future Double check that the updates in gsl-1.13 are included
      here, and also add support for the nmsimplex2rand algorithm
      in GSL.
  */
  template<class func_t=multi_funct,
           class vec_t=boost::numeric::ublas::vector<double> > class mmin_simp2 :
    public mmin_base<func_t,func_t,vec_t> {
      
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

#ifndef DOXYGEN_INTERNAL
      
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
    /// Center of simplex
    vec_t center;
    /// Desc
    vec_t delta;
    /// Distance of vector from center
    vec_t xmc;
    /// Squared simplex size
    double S2;
    
    /// Compute the center of the simplex 
    int compute_center() {
      size_t i, j;
      double val;
      
      for (j = 0; j < dim; j++) {
        val = 0.0;
        for (i = 0; i < dim+1; i++) {
          val += x1[i][j];
        }
        val /= ((double)dim)+1;
        center[j]=val;
      }
      
      return 0;
    }
    
    /** \brief Compute the size of the simplex
        
        Calculates simplex size as average sum of length of vectors
        from simplex center to corner points:
          
        \f$ (1/n) \sum || y - y_{\mathrm{center}} || \f$
    */
    double compute_size() {
        
      size_t P=dim+1;
      double ss=0.0;
        
      for(size_t i=0;i<P;i++) {
        for(size_t j=0;j<dim;j++) ws1[j]=x1[i][j];
        o2scl_cblas::daxpy(-1.0,dim,center,ws1);
        double t=o2scl_cblas::dnrm2(dim,ws1);
        ss += t*t;
      }
        
      /* Store squared size in the state */
      S2=ss/((double)P);
      return sqrt(ss/(double)(P));
    }

    /** \brief Move a corner of a simplex

        Moves a simplex corner scaled by coeff (negative value
        represents mirroring by the middle point of the "other" corner
        points) and gives new corner in xc and function value at xc 
        in \c newval.
    */
    virtual int try_corner_move(const double coeff, size_t corner, 
                                vec_t &xc, func_t &f, size_t nvar, 
                                double &newval) {
        
      size_t P=dim+1;
        
      /* 
         GSL: xc = (1-coeff)*(P/(P-1)) * center(all) + 
         ((P*coeff-1)/(P-1))*x_corner 
      */
      double alpha=(1.0-coeff)*((double)P)/((double)dim);
      double beta=(((double)P)*coeff-1.0)/((double)dim);
        
      for(size_t j=0;j<dim;j++) {
        xc[j]=center[j]*alpha;
        xc[j]+=x1[corner][j]*beta;
      }
        
      newval=f(nvar,xc);

      return 0;
    }

    /// Update point \c i in the simplex with values \c xx
    virtual int update_point(size_t i, vec_t &xx, double val) {

      const size_t P=dim+1;

      /* GSL: Compute delta = x - x_orig */
      for(size_t j=0;j<dim;j++) {
        delta[j]=xx[j];
        delta[j]-=x1[i][j];
      }

      /* GSL: Compute xmc = x_orig - c */
      for(size_t j=0;j<dim;j++) {
        xmc[j]=x1[i][j];
        xmc[j]-=center[j];
      }
        
      /* GSL: Update size: S2' = S2 + (2/P) * 
         (x_orig - c).delta + (P-1)*(delta/P)^2 
      */
      double d=o2scl_cblas::dnrm2(dim,delta);
      double xmcd=o2scl_cblas::ddot(dim,xmc,delta);
      S2 += (2.0 / ((double)P)) * xmcd + 
        ((((double)P) - 1.0) / ((double)P)) * (d * d / ((double)P));

      /* GSL: Update center:  c' = c + (x - x_orig) / P */
      double alpha=1.0/((double)P);
      for(size_t j=0;j<dim;j++) {
        center[j]-=alpha*x1[i][j];
        center[j]+=alpha*xx[j];
      }
        
      for(size_t j=0;j<dim;j++) {
        x1[i][j]=xx[j];
      }
      y1[i]=val;

      return 0;
    }
      
    /** \brief Contract the simplex towards the best point
          
        Function contracts the simplex in respect to best valued
        corner. All corners besides the best corner are moved.
        
        \comment 
        The GSL version requires a vector for workspace
        named 'xc' which is not required here.
        \endcomment
    */
    virtual int contract_by_best(size_t best, func_t &f, 
                                 size_t nvar) {
        
      /* GSL: Function contracts the simplex in respect to best valued
         corner. That is, all corners besides the best corner are
         moved. (This function is rarely called in practice, since
         it is the last choice, hence not optimised - BJG) 
      */
      size_t i,j, it;
      double newval;
      bool failed;

      int status=success;
        
      for (i=0;i<dim+1;i++) {
        
        if (i!=best) {
          
          for (j=0;j<dim;j++) {
            x1[i][j]=0.5*(x1[i][j]+x1[best][j]);
          }
          y1[i]=f(nvar,x1[i]);
          if (!std::isfinite(y1[i])) {
            std::string err=((std::string)"Function not finite (returned ")+
              dtos(y1[i])+" in mmin_simp2::contract_by_best().";
            O2SCL_ERR(err.c_str(),exc_ebadfunc);
          }

          /*
            if (avoid_nonzero==true && ret!=0) {
            std::cout << "Found problem in contract." << std::endl;

            // copy the old point to ws3
            for (j=0;j<dim;j++) {
            ws3[j]=2.0*x1[i][j]-x1[best][j];
            }
            
            // Try (21*best+20*xold)/41, (22*best+19*xold)/41, ...,
            // (40*best+xold)/41 to see if we can find a contraction
            // that works
            for(it=0;ret!=0 && it<20;it++) {
            for (j=0;j<dim;j++) {
            x1[i][j]=((20-it)*ws3[j]+(it+21)*x1[best][j])/41.0;
            }
            y1[i]=f(nvar,x1[i]);
            std::cout << "it, ret, x: " << it << " " << ret << " "
            << x << std::endl;
            }
            if (ret!=0) {
            O2SCL_CONV2_RET("Failed to find suitable contraction ",
            "in mmin_simp2::contract_by_best().",
            exc_efailed,this->err_nonconv);
            }
            }
          */
        }
      }
        
      compute_center();
      size=compute_size();
        
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
        protected to prevent the casual user from changing it.
    */
    bool avoid_nonzero;

#endif

  public:

    mmin_simp2() {
      dim=0;
      print_simplex=0;
      step_vec.resize(1);
      step_vec[0]=1.0;
      avoid_nonzero=false;
    }
    
    virtual ~mmin_simp2() {
    }

    /// Set the step sizes for each independent variable
    template<class vec2_t> int set_step(size_t nv, vec2_t &step) {
      if (nv>0) {
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
        O2SCL_ERR2("Tried to min over zero variables ",
                   " in mmin_simp2::mmin().",exc_einval);
      }

      int ret=0,status,iter=0;
      
      allocate(nn);

      vec_t ss(nn);
      for (size_t is=0;is<nn;is++) ss[is]=step_vec[is % step_vec.size()];
      ret=set(ufunc,nn,xx,ss);

      if(ret!=0) {
        return ret;
      }
  
      do {
        iter++;
          
        status=iterate();
        if(status) break;
          
        if(this->verbose>0) {
          print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
                     "mmin_simp2");
        }
    
        status=gsl_multimin_test_size(size,this->tol_abs);
          
      } while(status == GSL_CONTINUE && iter<this->ntrial);
        
      for (size_t i=0;i<nn;i++) xx[i]=x[i];
      fmin=fval;
  
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
        return ret;
      }
  
      do {
        iter++;
          
        status=iterate();
        if(status) break;

        if(this->verbose>0) {
          print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
                     "mmin_simp2");
        }
        
        status=gsl_multimin_test_size(size,this->tol_abs);
          
      } while(status == GSL_CONTINUE && iter<this->ntrial);
        
      for (i=0;i<((int)nn);i++) xx[i]=x[i];
      fmin=fval;
  
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

        The matrix \c sx should have the \c dim+1 initial points
        arranged in rows so that there are \c dim+1 rows and each row
        has \c dim columns.
    */
    template<class mat_t> 
    int mmin_simplex(size_t nn, mat_t &sx, double &fmin, 
                     func_t &ufunc) {
      
      int ret=0,i,status,iter=0;
      
      allocate(nn);
        
      ret=set_simplex(ufunc,sx);
      if(ret!=0) {
        return ret;
      }
  
      do {
        iter++;
          
        status=iterate();
        if(status) break;

        if(this->verbose>0) {
          print_iter(nn,x,x1,fval,iter,size,this->tol_abs,
                     "mmin_simp2");
        }
    
        status=gsl_multimin_test_size(size,this->tol_abs);
          
      } while(status == GSL_CONTINUE && iter<this->ntrial);
        
      for (i=0;i<((int)nn);i++) sx(0,i)=x[i];
      fmin=fval;
  
      this->last_ntrial=iter;

      if (iter>=this->ntrial) {
        std::string str="Exceeded maximum number of iterations ("+
          itos(this->ntrial)+") in mmin_simp2::mmin_simplex().";
        O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }
    
    /// Allocate the memory
    virtual int allocate(size_t n) {
        
      set_called=false;
      dim=n;

      x.resize(n);
      x1=new vec_t[n+1];
      for(size_t i=0;i<n+1;i++) x1[i].resize(n);
      y1.resize(n+1);
      ws1.resize(n);
      ws2.resize(n);
      ws3.resize(n);
      center.resize(n);
      delta.resize(n);
      xmc.resize(n);

      return success;
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
      if (!std::isfinite(y1[0])) {
        std::string err=((std::string)"Function not finite (returned ")+
          dtos(y1[0])+" in mmin_simp2::set().";
        O2SCL_ERR(err.c_str(),exc_ebadfunc);
      }
      for(i=0;i<dim;i++) x1[0][i]=ax[i];
  
      /* following points are initialized to x0+step_size */
      
      for (i=1;i<dim+1;i++) {
        for(size_t j=0;j<dim;j++) x1[i][j]=x[j];
        x1[i][i-1]=x1[i][i-1]+step_size[i-1];
        y1[i]=ufunc(dim,x1[i]);
      }
 
      /* Initialize simplex size */
      
      compute_center();
      size=compute_size();
        
      set_called=true;

      return success;
    }

    /** \brief Set the function and initial simplex

        The matrix \c sx should have the \c dim+1 initial points
        arranged in rows so that there are \c dim+1 rows and each row
        has \c dim columns.
     */
    template<class mat_t> 
    int set_simplex(func_t &ufunc, mat_t &sx) {

      if(dim==0) {
        O2SCL_ERR2("Memory not allocated in ",
                   "mmin_simp2::set_simplex().",exc_ebadlen);
      }
        
      func=&ufunc;

      for(size_t i=0;i<dim+1;i++) {
        for(size_t j=0;j<dim;j++) {
          x1[i][j]=sx(i,j);
        }
        y1[i]=ufunc(dim,x1[i]);
        if (!std::isfinite(y1[i])) {
          std::string err=((std::string)"Function not finite (returned ")+
            dtos(y1[i])+" in mmin_simp2::set_simplex().";
          O2SCL_ERR(err.c_str(),exc_ebadfunc);
        }
      }
        
      /* Initialize simplex size */
        
      compute_center();
      size=compute_size();
        
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
 
      /* get index of highest, second highest and lowest point */

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

      /* reflect the highest value */
        
      int ret1=try_corner_move(-1.0,hi,ws1,*func,dim,val);

      if (avoid_nonzero && ret1!=0) {
        std::cout << "Found problem move1: " << std::endl;
        for (size_t it=0;it<20 && ret1!=0;it++) {
          ret1=try_corner_move(-1.0+((double)it)/10.0,hi,ws1,
                               *func,dim,val);
          std::cout << "it,ret: " << it << " " << ret1 << std::endl;
        }
        if (ret1!=0) {
          O2SCL_ERR2("Failed to move corner (1) in ",
                     "mmin_simp2::iterate().",exc_efailed);
        }
      }
      
      if (std::isfinite(val) && val<y1[lo]) {

        /* reflected point becomes lowest point, try expansion */
        
        int ret2=try_corner_move(-2.0,hi,ws2,*func,dim,val2);

        if (avoid_nonzero && ret2!=0) {
          std::cout << "Found problem move2: " << std::endl;
          for (size_t it=0;it<20 && ret2!=0;it++) {
            ret2=try_corner_move(-2.0+((double)it)/6.7,hi,ws2,
                                 *func,dim,val2);
            std::cout << "it,ret: " << it << " " << ret2 << std::endl;
          }
          if (ret2!=0) {
            O2SCL_ERR2("Failed to move corner (2) in ",
                       "mmin_simp2::iterate().",exc_efailed);
          }
        }

        if (std::isfinite(val2) && val2<y1[lo]) {
          update_point(hi,ws2,val2);
        } else {
          update_point(hi,ws1,val);
        }
          
      } else if (!std::isfinite(val) || val > y1[s_hi]) {
        
        /* reflection does not improve things enough */
        
        if (std::isfinite(val) && val <= y1[hi]) {
            
          /* if trial point is better than highest point, replace
             highest point */
            
          update_point(hi,ws1,val);
        }
      
        /* try one-dimensional contraction */
        
        int ret3=try_corner_move(0.5,hi,ws2,*func,dim,val2);

        if (avoid_nonzero && ret3!=0) {
          std::cout << "Found problem move3: " << std::endl;
          for (size_t it=0;it<20 && ret3!=0;it++) {
            ret3=try_corner_move(0.025*((double)(19-it)),hi,ws2,
                                 *func,dim,val2);
            std::cout << "it,ret: " << it << " " << ret3 << std::endl;
          }
          if (ret3!=0) {
            O2SCL_ERR2("Failed to move corner (2) in ",
                       "mmin_simp2::iterate().",exc_efailed);
          }
        }
          
        if (std::isfinite(val2) && val2 <= y1[hi]) {

          update_point(hi,ws2,val2);

        } else {

          /* contract the whole simplex in respect to the best point */
          status=contract_by_best(lo,*func,dim);
          if(status != 0) {
            O2SCL_ERR("Function contract_by_best() failed in iterate().",
                      exc_efailed);
          }
            
        }

      } else {

        /* trial point is better than second highest point.
           Replace highest point by it */
          
        update_point(hi,ws1,val);
      }
  
      /* return lowest point of simplex as x */
      
      vector_min(dim+1,y1,lo,val);
      for(i=0;i<dim;i++) x[i]=x1[lo][i];
      fval=y1[lo];
      
      /* Update simplex size */
        
      if (S2 > 0) {
        size=sqrt(S2);
      } else {
        /* recompute if accumulated error has made size invalid */
        size=compute_size();
      }

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
                           double y, int iter,
                           double value, double limit,
                           std::string comment) {
        
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

    /// Return string denoting type("mmin_simp2")
    virtual const char *type() { return "mmin_simp2";}

#ifndef DOXYGEN_INTERNAL

  private:
  
    mmin_simp2<func_t,vec_t>
    (const mmin_simp2<func_t,vec_t> &);
    mmin_simp2<func_t,vec_t>& operator=
    (const mmin_simp2<func_t,vec_t>&);

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
