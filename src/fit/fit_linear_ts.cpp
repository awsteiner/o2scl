/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2013-2023, Andrew W. Steiner
  
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
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
// For gsl_blas_dnrm2
#include <gsl/gsl_blas.h>

#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>
#include <o2scl/fit_nonlin.h>
#include <o2scl/fit_linear.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

struct data {
  size_t n;
  double *y;
  double *sigma;
};

int expb_f(const gsl_vector *x, void *params, gsl_vector *f) {
  
  size_t n=((struct data *)params)->n;
  double *y=((struct data *)params)->y;
  double *sigma=((struct data *) params)->sigma;
  
  double a=gsl_vector_get(x,0);
  double b=gsl_vector_get(x,1);

  /* Model Yi=A*exp(-lambda*i)+b */
  for (size_t i=0;i<n;i++) {
    double Yi=a*exp(((double)i)/10.0)+b*sqrt(((double)i)/10.0);
    gsl_vector_set(f,i,(Yi-y[i])/sigma[i]);
  }
  
  return success;
}

int expb_df(const gsl_vector *x, void *params, gsl_matrix *J) {

  size_t n=((struct data *)params)->n;
  double *y=((struct data *)params)->y;
  double *sigma=((struct data *) params)->sigma;

  double a=gsl_vector_get(x,0);
  double b=gsl_vector_get(x,1);

  for (size_t i=0;i<n;i++) {
    gsl_matrix_set(J,i,0,exp(((double)i)/10.0));
    gsl_matrix_set(J,i,1,sqrt(((double)i)/10.0));
  }

  return success;
}

int expb_fdf(const gsl_vector *x, void *params,
	     gsl_vector *f, gsl_matrix *J) {
  expb_f(x,params,f);
  expb_df(x,params,J);
  return success;
}

#ifdef O2SCL_ARMA
#include <armadillo>

template<class vec_t, class mat_t>
class fit_linear_arma : public fit_linear<vec_t,mat_t> {
protected:
  virtual void fit_svd(size_t ndat, size_t npar) {
    svd(this->A,this->S,this->Q,this->A);
    return;
  }
};

#endif

#ifdef O2SCL_SET_EIGEN
#include <eigen3/Eigen/Dense>

template<class vec_t, class mat_t>
class fit_linear_eigen : public fit_linear<vec_t,mat_t> {

protected:

  virtual void fit_svd(size_t ndat, size_t npar) {

    Eigen::JacobiSVD<Eigen::MatrixXd> 
      svd(this->A,Eigen::ComputeThinU | Eigen::ComputeThinV);
    this->S=svd.singularValues();
    this->A=svd.matrixU();
    this->Q=svd.matrixV();

    return;
  }
};

#endif
  
int main(void) {
  test_mgr tm;
  tm.set_output_level(2);

  cout.setf(ios::scientific);

  static const size_t ndat=20;
  static const size_t npar=2;

  // Set up data
  ubvector xdat(ndat), ydat(ndat), yerr(ndat), parms(npar), parms_bench(npar);
  ubmatrix xpred(ndat,npar);
  for(size_t i=0;i<ndat;i++) {
    xdat[i]=((double)i)/10.0;
    ydat[i]=3.0*exp(xdat[i])-sqrt(xdat[i])+3.0e-1*sin(1000.0*xdat[i]);
    yerr[i]=1.0;
    xpred(i,0)=exp(xdat[i]);
    xpred(i,1)=sqrt(xdat[i]);
  }

  // Data for GSL version
  double y_gsl[ndat], sigma_gsl[ndat];
  gsl_vector *y_gsl2=gsl_vector_alloc(ndat), *p_gsl2=gsl_vector_alloc(npar);
  gsl_matrix *xpred_gsl2=gsl_matrix_alloc(ndat,npar);
  struct data d={ndat,y_gsl,sigma_gsl};
  for(size_t i=0;i<ndat;i++) {
    double tmp=((double)i)/10.0;
    y_gsl[i]=3.0*exp(tmp)-sqrt(tmp)+3.0e-1*sin(1000.0*tmp);
    sigma_gsl[i]=1.0;
    gsl_vector_set(y_gsl2,i,y_gsl[i]);
    gsl_matrix_set(xpred_gsl2,i,0,exp(tmp));
    gsl_matrix_set(xpred_gsl2,i,1,sqrt(tmp));
  }
  
  // Chi-squared and covariance matrices
  double chi2, chi2_bench;
  ubmatrix covar(npar,npar), covar_bench(npar,npar);
  gsl_matrix *covar_gsl=gsl_matrix_alloc(npar,npar);
    
  // O2scl linear fit
  {
    cout << "O2scl linear fit: " << endl;
    parms_bench[0]=2.0;
    parms_bench[1]=-0.5;
    fit_linear<> gfl;
    gfl.fit(npar,ndat,ydat,xpred,parms_bench,covar_bench,chi2_bench);

    // Output results
    cout << "Parameters: " << parms_bench[0] << " " << parms_bench[1] << endl;
    cout << "Covariance matrix: " << endl;
    matrix_out(cout,npar,npar,covar_bench);
    cout << "Chi-squared: " << chi2_bench << " rank: " << gfl.rank << endl;

    // Test
    tm.test_rel(parms_bench[0],3.0,3.0e-2,"o2scl linear p0");
    tm.test_rel(parms_bench[1],-1.0,3.0e-2,"o2scl linear p1");
    tm.test_gen(gfl.rank==2,"o2scl linear rank");
  
    // Final values
    if (false) {
      cout << "Final values: " << endl;
      for(size_t i=0;i<ndat;i++) {
	cout << xdat[i] << " " << ydat[i] << " " 
	     << parms_bench[0]*exp(xdat[i])-
	  parms_bench[1]*sqrt(xdat[i]) << endl;
      }
    }
    cout << endl;
  }

#ifdef O2SCL_NEVER_DEFINED
#ifdef O2SCL_ARMA

  {
    arma::rowvec aparms(npar), aydat(ndat);
    arma::mat axpred(ndat,npar), acovar(npar,npar);

    aparms[0]=2.0;
    aparms[1]=-0.5;

    cout << "Armadillo linear fit: " << endl;
    fit_linear_arma<arma::rowvec,arma::mat> afl;
    
    for(size_t i=0;i<ndat;i++) {
      eydat[i]=ydat[i];
      expred(i,0)=xpred(i,0);
      expred(i,1)=xpred(i,1);
    }

    cout << "Armadillo linear fit: " << endl;
    afl.fit(npar,ndat,aydat,axpred,aparms,acovar,chi2);

    cout << "Parameters: " << aparms[0] << " " << aparms[1] << endl;
    cout << "Covariance matrix: " << endl;
    cout << " " << acovar(0,0);
    cout << " " << acovar(0,1) << endl;
    cout << acovar(1,0) << "  ";
    cout << acovar(1,1) << endl;
    cout << "Chi-squared: " << chi2 << endl;

    tm.test_rel_vec(2,parms_bench,aparms,4.0e-11,
		    "Armadillo linear parms vs. O2scl linear parms");
    tm.test_rel(chi2_bench,chi2,4.0e-11,
		"Armadillo linear chi2 vs. O2scl linear chi2");
    tm.test_rel_mat(2,2,covar_bench,acovar,4.0e-11,
		    "Armadillo linear covar vs. O2scl linear covar");

    cout << endl;
  }

#endif
#endif

#ifdef O2SCL_SET_EIGEN

  {
    Eigen::VectorXd eparms(npar), eydat(ndat);
    Eigen::MatrixXd expred(ndat,npar), ecovar(npar,npar);

    eparms[0]=2.0;
    eparms[1]=-0.5;

    for(size_t i=0;i<ndat;i++) {
      eydat[i]=ydat[i];
      expred(i,0)=xpred(i,0);
      expred(i,1)=xpred(i,1);
    }

    cout << "Eigen linear fit: " << endl;
    fit_linear_eigen<Eigen::VectorXd,Eigen::MatrixXd> efl;
    efl.fit(npar,ndat,eydat,expred,eparms,ecovar,chi2);

    cout << "Parameters: " << eparms[0] << " " << eparms[1] << endl;
    cout << "Covariance matrix: " << endl;
    cout << " " << ecovar(0,0);
    cout << " " << ecovar(0,1) << endl;
    cout << ecovar(1,0) << "  ";
    cout << ecovar(1,1) << endl;
    cout << "Chi-squared: " << chi2 << endl;

    tm.test_rel_vec(2,parms_bench,eparms,4.0e-11,
		    "Eigen linear parms vs. O2scl linear parms");
    tm.test_rel(chi2_bench,chi2,4.0e-11,
		"Eigen linear chi2 vs. O2scl linear chi2");
    tm.test_rel_mat(2,2,covar_bench,ecovar,4.0e-11,
		    "Eigen linear covar vs. O2scl linear covar");

    cout << endl;
  }

#endif

  // GSL linear fit
  {
    cout << "GSL linear fit: " << endl;
    gsl_vector_set(p_gsl2,0,2.0);
    gsl_vector_set(p_gsl2,1,-0.5);

    gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(ndat,npar);
    gsl_multifit_linear(xpred_gsl2,y_gsl2,p_gsl2,covar_gsl,&chi2,work);

    cout << "Parameters: " << gsl_vector_get(p_gsl2,0) << " "
	 << gsl_vector_get(p_gsl2,1) << endl;
    cout << "Covariance matrix: " << endl;
    cout << " " << gsl_matrix_get(covar_gsl,0,0);
    cout << " " << gsl_matrix_get(covar_gsl,0,1) << endl;
    cout << gsl_matrix_get(covar_gsl,1,0) << "  ";
    cout << gsl_matrix_get(covar_gsl,1,1) << endl;
    cout << "Chi-squared: " << chi2 << endl;

    tm.test_rel_vec(2,parms_bench,gsl_vector_wrap(p_gsl2),1.0e-12,
		    "gsl linear parms vs. O2scl linear parms");
    tm.test_rel(chi2,chi2_bench,1.0e-12,
		"gsl linear chi2 vs. O2scl linear chi2");
    tm.test_rel_mat(2,2,covar_bench,gsl_matrix_wrap(covar_gsl),1.0e-12,
		    "gsl linear covar vs. O2scl linear covar");

    cout << endl;

    gsl_multifit_linear_free(work);
  }
  
  // O2scl nonlinear version
  {
    cout << "O2scl nonlinear fit:" << endl;
    fit_nonlin<> gf;
    vector<string> vars={"a","b"};
    fit_funct_strings ffs("a*exp(x)+b*sqrt(x)",vars,"x");
    chi_fit_funct<ubvector,ubmatrix,fit_funct_strings> 
      cff(ndat,xdat,ydat,yerr,ffs);
    cff.auto_jac.set_epsrel(1.0e-4);
    
    parms[0]=2.0;
    parms[1]=-0.5;
  
    gf.fit(npar,parms,covar,chi2,cff);

    double variance=0.0;
    for(size_t i=0;i<ndat;i++) {
      variance+=pow(ydat[i]-parms[0]*exp(xdat[i])-parms[1]*sqrt(xdat[i]),2.0);
    }
    variance/=(ndat-npar);
    covar*=variance;

    cout << "Parameters: " << parms[0] << " " << parms[1] << endl;
    cout << "Covariance matrix: " << endl;
    matrix_out(cout,npar,npar,covar);
    cout << "Chi-squared: " << chi2 << endl;
    
    tm.test_rel_vec(2,parms_bench,parms,4.0e-11,
		    "O2scl nonlinear parms vs. O2scl linear parms");
    tm.test_rel(chi2,chi2_bench,4.0e-11,
		"O2scl nonlinear chi2 vs. O2scl linear chi2");
    tm.test_rel_mat(2,2,covar_bench,covar,4.0e-11,
		    "O2scl nonlinear covar vs. O2scl nonlinear covar");

    cout << endl;
  }
  
  // GSL nonlinear version 
  {
    cout << "GSL nonlinear fit: " << endl;
    const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;
    gsl_matrix *J=gsl_matrix_alloc(ndat,2);

    int status;
    size_t iter=0;

    const size_t p=2;

    gsl_multifit_function_fdf f;

    double x_init[2]={2.0,-0.5};

    gsl_vector_view x=gsl_vector_view_array(x_init,p);
    
    f.f=&expb_f;
    f.df=&expb_df;
    f.fdf=&expb_fdf;
    f.n=ndat;
    f.p=p;
    f.params=&d;

    T=gsl_multifit_fdfsolver_lmsder;
    s=gsl_multifit_fdfsolver_alloc(T,ndat,p);
    gsl_multifit_fdfsolver_set(s,&f,&x.vector);

    do {
      iter++;
      status=gsl_multifit_fdfsolver_iterate(s);
      if (status)
	break;
      status=gsl_multifit_test_delta(s->dx,s->x,1e-4,1e-4);
				     
    } while (status == gsl_continue && iter < 500);
    
    gsl_multifit_fdfsolver_jac(s,J);
    gsl_multifit_covar(J,0.0,covar_gsl);

    double variance=0.0;
    for(size_t i=0;i<ndat;i++) {
      variance+=pow(ydat[i]-gsl_vector_get(s->x,0)*exp(xdat[i])-
		    gsl_vector_get(s->x,1)*sqrt(xdat[i]),2.0);
    }
    variance/=(ndat-npar);
    gsl_matrix_set(covar_gsl,0,0,gsl_matrix_get(covar_gsl,0,0)*variance);
    gsl_matrix_set(covar_gsl,0,1,gsl_matrix_get(covar_gsl,0,1)*variance);
    gsl_matrix_set(covar_gsl,1,0,gsl_matrix_get(covar_gsl,1,0)*variance);
    gsl_matrix_set(covar_gsl,1,1,gsl_matrix_get(covar_gsl,1,1)*variance);
    
    chi2=pow(gsl_blas_dnrm2(s->f),2.0);
    
    cout << "Parameters: " << gsl_vector_get(s->x,0) << " "
	 << gsl_vector_get(s->x,1) << endl;
    cout << "Covariance matrix: " << endl;
    cout << " " << gsl_matrix_get(covar_gsl,0,0);
    cout << " " << gsl_matrix_get(covar_gsl,0,1) << endl;
    cout << gsl_matrix_get(covar_gsl,1,0) << "  ";
    cout << gsl_matrix_get(covar_gsl,1,1) << endl;
    cout << "Chi-squared: " << chi2 << endl;

    tm.test_rel_vec(2,parms_bench,gsl_vector_wrap(p_gsl2),1.0e-12,
		    "gsl nonlinear parms vs. O2scl linear parms");
    tm.test_rel(chi2,chi2_bench,1.0e-12,
		"gsl nonlinear chi2 vs. O2scl linear chi2");
    tm.test_rel_mat(2,2,covar_bench,gsl_matrix_wrap(covar_gsl),1.0e-12,
		    "gsl nonlinear covar vs. O2scl linear covar");

    cout << endl;

    gsl_multifit_fdfsolver_free(s);

  } 

  // O2scl linear fit without balancing
  if (true) {
    cout << "O2scl unbalanced linear fit: " << endl;
    parms_bench[0]=2.0;
    parms_bench[1]=-0.5;
    fit_linear<> gfl;
    gfl.tol=0.99;
    gfl.column_scaling=false;
    gfl.fit(npar,ndat,ydat,xpred,parms_bench,covar_bench,chi2_bench);

    // Output results
    cout << "Parameters: " << parms_bench[0] << " " << parms_bench[1] << endl;
    cout << "covariance matrix: " << endl;
    matrix_out(cout,npar,npar,covar_bench);
    cout << "Chi-squared: " << chi2_bench << " rank: " << gfl.rank << endl;

    // Test
    tm.test_gen(gfl.rank==1,"o2scl unbalanced rank");
  
    // Final values
    if (false) {
      cout << "Final values: " << endl;
      for(size_t i=0;i<ndat;i++) {
	cout << xdat[i] << " " << ydat[i] << " " 
	     << parms_bench[0]*exp(xdat[i])-
	  parms_bench[1]*sqrt(xdat[i]) << endl;
      }
    }
    cout << endl;
  }
  
  // FIXME: Commented out for gsl-2.0 upgrade. 
  /*
  // GSL linear fit without scaling
  {
    cout << "GSL unbalanced linear fit: " << endl;
    gsl_vector_set(p_gsl2,0,2.0);
    gsl_vector_set(p_gsl2,1,-0.5);

    gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(ndat,npar);
    size_t rank;
    // Use a large tolerance to force it to ignore a fit component
    gsl_multifit_linear_usvd(xpred_gsl2,y_gsl2,0.99,&rank,
        p_gsl2,covar_gsl,&chi2,work);

    cout << "Parameters: " << gsl_vector_get(p_gsl2,0) << " "
	 << gsl_vector_get(p_gsl2,1) << endl;
    cout << "Covariance matrix: " << endl;
    cout << gsl_matrix_get(covar_gsl,0,0) << " ";
    cout << gsl_matrix_get(covar_gsl,0,1) << endl;
    cout << gsl_matrix_get(covar_gsl,1,0) << " ";
    cout << gsl_matrix_get(covar_gsl,1,1) << endl;
    cout << "Chi-squared: " << chi2 << " rank: " << rank << endl;

    // FIXME: temporarily commented out for gsl-2.0 upgrade
    tm.test_rel_vec(2,parms_bench,gsl_vector_wrap(p_gsl2),1.0e-12,
    "gsl unbalanced linear parms vs. O2scl unbalanced");
    tm.test_rel(chi2,chi2_bench,1.0e-12,
    "gsl unbalanced linear chi2 vs. O2scl unbalanced");
    tm.test_rel_mat(2,2,covar_bench,gsl_matrix_wrap(covar_gsl),1.0e-12,
    "gsl unbalanced linear covar vs. O2scl unbalanced");

    cout << endl;

    gsl_multifit_linear_free(work);
  }
  */

  tm.report();
  return 0;
}


