 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2021, Julien Garaud and Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_BV_MULTISHOOT_H
#define O2SCL_ODE_BV_MULTISHOOT_H

/** \file ode_bv_multishoot.h
    \brief File defining \ref o2scl::ode_bv_multishoot 
*/

#include <string>
#include <o2scl/astep.h>
#include <o2scl/astep_gsl.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/gsl_mroot_hybrids.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Solve a ODE boundary value problem by multishooting 

      This class is experimental.

      \todo Improve documentation a little and create
      testing code
  */
  template <class func_t=ode_funct<>,
    class vec_t=ubvector, class alloc_vec_t=ubvector, 
    class alloc_t=ubvector_alloc,class vec_int_t=ubvector_int_base, 
    class mat_t=ubmatrix> class ode_bv_multishoot {
    
    public: 
  
    ode_bv_multishoot(){
      oisp=&def_ois;
      mrootp=&def_mroot;
    }
  
    virtual ~ode_bv_multishoot() {
    }
 
    /* Solve the boundary value problem on a mesh */
    virtual int solve(vec_t &mesh, int &n_func, vec_t &y_start, 
		      func_t &left_b, func_t &right_b, func_t &extra_b, 
		      func_t &derivs, vec_t &x_save, mat_t &y_save) {
    
      /* Make copies of the input for later access */
      this->l_mesh=&mesh;
      this->l_n_func=&n_func;
      this->l_y_start=&y_start;
      this->l_left_b=&left_b;
      this->l_right_b=&right_b;
      this->l_extra_b=&extra_b;
      this->l_derivs=&derivs;
      this->l_x_save=&x_save;
      this->l_y_save=&y_save;

      /* vector of variables */
      int nsolve=y_start.size();
      ubvector sx(nsolve),sy(nsolve);
      sx=y_start;
    
      /* Equation solver */
      mm_funct_mfptr<ode_bv_multishoot<func_t,vec_t,
      alloc_vec_t,alloc_t,vec_int_t> >
      mfm(this,&ode_bv_multishoot<func_t,vec_t,
	  alloc_vec_t,alloc_t,vec_int_t>::solve_fun);

      /* Run multishooting and save at the last step */
      this->save=false;
      int ret=this->mrootp->msolve(nsolve,sx,mfm);
      this->save=true;
      solve_fun(nsolve,sx,sy);

      return ret;
    }

    /* Set initial value solver */
    int set_iv(ode_iv_solve<func_t,vec_t,alloc_vec_t,alloc_t> &ois) {
      oisp=&ois;
      return 0;
    }

    /* Set equation solver */
    int set_mroot(mroot<mm_funct<> > &root) {
      mrootp=&root;
      return 0;
    }  

    /* Default initial value solver */
    ode_iv_solve<func_t,vec_t,alloc_vec_t,alloc_t> def_ois;
  
    /* Default equation solver */
    gsl_mroot_hybrids<mm_funct<> > def_mroot;

#ifndef DOXYGEN_INTERNAL

    protected: 

    /// The initial value solver
    ode_iv_solve<func_t,vec_t,alloc_vec_t,alloc_t> *oisp;

    /// The equation solver
    gsl_mroot_hybrids<mm_funct<> > *mrootp;

    /// Desc
    vec_t *l_mesh;
    /// Desc
    vec_t *l_y_start;
    /// Desc
    func_t *l_left_b;
    /// Desc
    func_t *l_right_b;
    /// Desc
    func_t *l_extra_b;
    /// Desc
    func_t *l_derivs;
    /// Desc
    int *l_n_func;
    /// Desc
    vec_t *l_x_save;
    /// Desc
    mat_t *l_y_save;
    /// Desc
    bool save;

    /// Function to solve
    int solve_fun(size_t nv, const vec_t &sx, vec_t &sy) {
    
      double xa,xb=0.0,h;
      ubvector y((*this->l_n_func)),y2((*this->l_n_func));
    
      /* We update y_start in order that derivs knows 
	 all the values of parameters */
      for(size_t i=0;i<(*this->l_y_start).size();i++) {
	(*this->l_y_start)[i]=sx[i];
      }

      /* A loop on each subinterval */
      for(size_t k=0;k<(*this->l_mesh).size()-1;k++) {
      
	xa=(*this->l_mesh)[k];
	xb=(*this->l_mesh)[k+1];
	h=(xb-xa)/100.0;
      
	/* We load function's value at the left point of the sub-interval */
	if (k==0) {
	  (*this->l_left_b)(xa,(*this->l_n_func),sx,y);
	} else {
	  for(int i=0;i<(*this->l_n_func);i++)
	    y[i]=sx[i+(*this->l_n_func)*(1+k)];
	}
      
	if (this->save) {

	  /* iv_solver if we save */
	  int ngrid=((*this->l_x_save).size()-1)/((*this->l_mesh).size()-1)+1;
	  ubvector xxsave(ngrid);
	  ubmatrix yysave(ngrid,(*this->l_n_func));

	  if (k!=((*this->l_mesh).size()-2)) {
	    xb=(*this->l_mesh)[k+1]-((*this->l_mesh)[k+1]-
				     (*this->l_mesh)[k])/ngrid;
	  }
	
	  this->oisp->solve_grid(xa,xb,h,(*this->l_n_func),y,ngrid,
				 xxsave,yysave,(*this->l_derivs));    

	  for(int i=0;i<ngrid;i++) {
	    (*this->l_x_save)[i+k*(ngrid)]=xxsave[i];
	    for(int j=0;j<(*this->l_n_func);j++) {
	      (*this->l_y_save)[i+k*(ngrid)][j]=yysave[i][j];
	    }
	  }	  
	
	} else {
	
	  /* iv_solver if we don't save */
	  this->oisp->solve_final_value
	    (xa,xb,h,(*this->l_n_func),y,y2,(*this->l_derivs));

	}
      
	/* Then we load values at the end of sub-interval */
	if (k==(*this->l_mesh).size()-2) {
	  (*this->l_right_b)(xb,(*this->l_n_func),sx,y);
	} else {
	  for(int i=0;i<(*this->l_n_func);i++) {
	    y[i]=sx[i+(*this->l_n_func)*(2+k)];
	  }
	}

	/* Now we take the difference */
	for(int i=0;i<(*this->l_n_func);i++) {
	  //sy[i+k*(*this->l_n_func)]=(y2[i]-y[i])/y[i];
	  sy[i+k*(*this->l_n_func)]= y2[i]-y[i];
	}

      }

      /* Then load Extra boundary condition */
      (*this->l_extra_b)(xb,(*this->l_n_func),sx,y);
      for(int i=0;i<(*this->l_n_func);i++) {
	sy[i+(int((*this->l_mesh).size()-1))*(*this->l_n_func)]=y[i];
      }
  
      return 0;
    }

#endif
  
  };
 
#ifndef DOXYGEN_NO_O2NS 
}
#endif

#endif
