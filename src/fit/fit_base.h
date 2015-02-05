/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_FIT_BASE_H
#define O2SCL_FIT_BASE_H

/** \file fit_base.h
    \brief File defining \ref o2scl::fit_base and fitting functions
*/

#include <string>

#include <o2scl/jacobian.h>
#include <o2scl/fparser.h>
#include <o2scl/mm_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Array of multi-dimensional functions typedef (C++11 version)
  typedef std::function<
    double(size_t,const boost::numeric::ublas::vector<double> &, 
	   double)> fit_funct11;

  /** \brief String fitting function
      
      Default template arguments
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class fit_funct11_strings {

    public:

    /** \brief Specify a fitting function through a string
     */
    fit_funct11_strings(std::string formula, std::string parms, 
			 std::string var, int nauxp=0, 
			 std::string auxp="") {
      if(nauxp<1) {
	std::string both=parms+","+var;
	fpw.Parse(formula,both);
	naux=0;
	st_auxp="";
      } else {
	std::string all=parms+","+var+","+auxp;
	fpw.Parse(formula,all);
	naux=nauxp;
	aux=new double[naux];
	st_auxp=auxp;
      }
      st_form=formula;
      st_parms=parms;
      st_var=var;
    }

    virtual ~fit_funct11_strings() {
      if (naux>0) {
	delete[] aux;
      }
    }
  
    /** \brief Set the values of the auxilliary parameters that were
	specified in \c auxp in the constructor
    */
    int set_aux_parms(vec_t &ap) {
      for(int i=0;i<naux;i++) {
	aux[i]=ap[i];
      }
      return 0;
    }

    /** \brief Using parameters in \c p, predict \c y given \c x
     */
    virtual double operator()(size_t np, const vec_t &p, double x) {

      size_t i;
      double y;
      if(naux<1) {
	double *both=new double[np+1];
	for(i=0;i<np;i++) both[i]=p[i];
	both[np]=x;
	y=fpw.Eval(both);
	delete[] both;
      } else {
	double *all=new double[np+1+naux];
	for(i=0;i<np;i++) all[i]=p[i];
	all[np]=x;
	for(i=np+1;i<np+1+naux;i++) all[i]=aux[i-np-1];
	y=fpw.Eval(all);
	delete[] all;
      }
      return y;
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// The function evaluation object
    FunctionParser fpw;

    /// \name The auxillary parameters
    //@{
    int naux;
    double *aux;
    std::string st_auxp;
    //@}

    /// The formula
    std::string st_form; 
    /// The parameters
    std::string st_parms; 
    /// The variables
    std::string st_var; 

    fit_funct11_strings() {};

    /// Specify the strings which define the fitting function
    int set_function(std::string formula, std::string parms, 
		     std::string var, int nauxp=0, std::string auxp="");

    private:
    
    fit_funct11_strings(const fit_funct11_strings &);
    fit_funct11_strings& operator=(const fit_funct11_strings&);

#endif

  };

  /** \brief Generalized fitting function [abstract base]

      Default template arguments
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c mat_t - \ref boost::numeric::ublas::matrix \< double \>
  */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class gen_fit_funct {

    public:  
    
    gen_fit_funct() {}

    virtual ~gen_fit_funct() {}
    
    /** \brief Using parameters in \c p, compute the 
	relative deviations in \c f
    */
    virtual void operator()(size_t np, const vec_t &p, size_t nd,
			    vec_t &f)=0;

    /** \brief Using parameters in \c p, compute the Jacobian
	in \c J
    */
    virtual void jac(size_t np, vec_t &p, size_t nd, vec_t &f,
		     mat_t &J)=0;

    /// Return the number of data points
    virtual size_t get_ndata()=0;
    
#ifndef DOXYGEN_INTERNAL
    
    private:
    
    gen_fit_funct(const gen_fit_funct &);
    gen_fit_funct& operator=(const gen_fit_funct&);
    
#endif

  };

  /** \brief Standard fitting function based on one-dimensional
      data with a numerical Jacobian

      This class specifies the deviations (in <tt>operator()</tt> ) and
      Jacobian (in \ref jac()) for a fitting class like \ref fit_nonlin.
      It assumes a one-dimensional data set with no uncertainty in the
      abcissae and a fitting function specified in a form similar to
      \ref fit_funct11.
      \comment
      For some reason the reference to operator() above doesn't work
      in doxygen.
      \endcomment

      The default method for numerically computing the Jacobian is
      from \ref jacobian_gsl. This default is identical to the GSL
      approach, except that the default value of \ref
      jacobian_gsl::epsmin is non-zero. See \ref jacobian_gsl for more
      details.

      Default template arguments
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c mat_t - \ref boost::numeric::ublas::matrix \< double \>
      - \c func_t - \ref fit_funct11

      \future Allow a user-specified Jacobian or make that into
      a separate class?
      \future Default constructor?
  */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double>, 
    class fit_func_t=fit_funct11> class chi_fit_funct : 
    public gen_fit_funct<vec_t,mat_t> {
    
  public:
  
  /** \brief Create an object with specified data and specified 
      fitting function
  */
  chi_fit_funct(size_t ndat, const vec_t &xdat, const vec_t &ydat, 
		const vec_t &yerr, fit_func_t &fun) {
    ndat_=ndat;
    xdat_=&xdat;
    ydat_=&ydat;
    yerr_=&yerr;
    fun_=&fun;
    
    mfm=std::bind(std::mem_fn<int(size_t,const vec_t &,vec_t &)>
		  (&chi_fit_funct::jac_mm_funct),this,
		  std::placeholders::_1,std::placeholders::_2,
		  std::placeholders::_3);

    auto_jac.set_function(mfm);
    //double sqrt_dbl_eps=sqrt(std::numeric_limits<double>::epsilon());
    //auto_jac.set_epsrel(sqrt_dbl_eps);
  }

  /** \brief Set the data to be fit 
   */
  void set_data(size_t ndat, const vec_t &xdat, const vec_t &ydat, 
		const vec_t &yerr) {
    ndat_=ndat;
    xdat_=&xdat;
    ydat_=&ydat;
    yerr_=&yerr;
    return;
  }

  /** \brief Set the fitting function
   */
  void set_func(fit_func_t &fun) {
    fun_=&fun;
    return;
  }

  /** \brief Return \f$ \chi^2 \f$
   */
  virtual double chi2(size_t np, const vec_t &p) {
    double ret=0.0;
    for(size_t i=0;i<ndat_;i++) {
      double yi=((*fun_)(np,p,(*xdat_)[i])-(*ydat_)[i])/((*yerr_)[i]);
      ret+=yi*yi;
    }
    return ret;
  }
  
  /** \brief Using parameters in \c p, compute the 
      relative deviations in \c f
  */
  virtual void operator()(size_t np, const vec_t &p, size_t nd, vec_t &f) {

    for(size_t i=0;i<nd;i++) {
      double yi=(*fun_)(np,p,(*xdat_)[i]);
      f[i]=(yi-(*ydat_)[i])/((*yerr_)[i]);
    }
    return;
  }
  
  /** \brief Using parameters in \c p, compute the Jacobian
      in \c J
  */
  virtual void jac(size_t np, vec_t &p, size_t nd, vec_t &f,
		   mat_t &J) {
    
    auto_jac(np,p,nd,f,J);
    
    return;
  }

  /// Return the number of data points
  virtual size_t get_ndata() {
    return ndat_;
  }

  /// Automatic Jacobian object
  jacobian_gsl<std::function<int(size_t,const vec_t &,vec_t &)>,
  vec_t,mat_t> auto_jac;

#ifndef DOXYGEN_INTERNAL
  
  protected:

  /// Reformulate <tt>operator()</tt> into a \ref mm_funct11 object
  int jac_mm_funct(size_t np, const vec_t &p, vec_t &f) {
    operator()(np,p,ndat_,f);
    return 0;
  }
  
  /// Function object for Jacobian object
  std::function<int(size_t,const vec_t &,vec_t &)> mfm;
  
  /// \name Data and uncertainties
  //@{
  size_t ndat_;
  const vec_t *xdat_;
  const vec_t *ydat_;
  const vec_t *yerr_;
  //@}

  /// Fitting function
  fit_func_t *fun_;
  
  private:
  
  chi_fit_funct(const chi_fit_funct &);
  chi_fit_funct& operator=(const chi_fit_funct&);
  
#endif
  
  };

  // ----------------------------------------------------------------
  // Fitter base
  // ----------------------------------------------------------------

  /** \brief Non-linear least-squares fitting [abstract base]
   */
  template<class func_t=fit_funct11,
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class fit_base {

    public:

    fit_base() {
      ntrial=500;
      tol_abs=1.0e-4;
      tol_rel=1.0e-4;
      verbose=0;
    }

    virtual ~fit_base() {}
    
    /// Maximum number of iterations (default 500)
    size_t ntrial;
    
    /// Absolute tolerance (default 1.0e-4)
    double tol_abs;
    
    /// (default 1.0e-4)
    double tol_rel;

    /** \brief Print out iteration information.
	
	Depending on the value of the variable verbose, this prints out
	the iteration information. If verbose=0, then no information is
	printed, while if verbose>1, then after each iteration, the
	present values of x and y are output to std::cout along with the
	iteration number. If verbose>=2 then each iteration waits for a
	character.  
    */
    virtual int print_iter(size_t nv, vec_t &x, double y, int iter,
			   double value=0.0, double limit=0.0) {
      if (verbose<=0) return 0;
  
      size_t i;
      char ch;

      std::cout << "Iteration: " << iter << std::endl;
      std::cout << "x: ";
      for(i=0;i<nv;i++) std::cout << x[i] << " ";
      std::cout << std::endl;
      std::cout << "y: " << y << " Val: " << value << " Lim: " << limit 
      << std::endl;
      if (verbose>1) {
	std::cout << "Press a key and type enter to continue. ";
	std::cin >> ch;
      }

      return 0;
    }

    /** \brief Fit function \c fitfun using parameters in \c parms as
	initial guesses
	
	The covariance matrix for the parameters is returned in \c covar
	and the value of \f$ \chi^2 \f$ is returned in \c chi2.
	
    */
    virtual int fit(size_t npar, vec_t &parms, mat_t &covar, 
		    double &chi2, func_t &fitfun)=0;
    
    /** \brief An integer describing the verbosity of the output 
     */
    int verbose;

    /// Return string denoting type ("fit_base")
    virtual const char *type() { return "fit_base"; }

    /// The number of data points
    size_t n_dat;

    /// The number of parameters
    size_t n_par;


  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
