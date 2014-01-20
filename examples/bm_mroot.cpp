/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/uvector_tlate.h>
#include <o2scl/gsl_mroot_hybrids.h>

/*
  This compares the speed of gsl_multiroot_fdfsolver_hybridsj to 
  several different implementations of o2scl::gsl_mroot_hybrids.
  
  Should examine cost of virtual functions and average over a few
  timings. Also, should try with derivatives and arrays in addition to
  using both separately.
*/

int gsl_f(const gsl_vector *x, void *params,  gsl_vector *f) {
  gsl_vector_set(f,0,sin(gsl_vector_get(x,1)-0.2));
  gsl_vector_set(f,1,sin(gsl_vector_get(x,0)-0.25));
  return 0;
}

int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  gsl_matrix_set(J,0,0,0.0);
  gsl_matrix_set(J,0,1,cos(gsl_vector_get(x,1)-0.2));
  gsl_matrix_set(J,1,0,cos(gsl_vector_get(x,0)-0.25));
  gsl_matrix_set(J,1,1,0.0);
  return 0;
}

int gsl_fdf(const gsl_vector *x, void *params,
	    gsl_vector *f, gsl_matrix *J) {
  gsl_vector_set(f,0,sin(gsl_vector_get(x,1)-0.2));
  gsl_vector_set(f,1,sin(gsl_vector_get(x,0)-0.25));
  gsl_matrix_set(J,0,0,0.0);
  gsl_matrix_set(J,0,1,cos(gsl_vector_get(x,1)-0.2));
  gsl_matrix_set(J,1,0,cos(gsl_vector_get(x,0)-0.25));
  gsl_matrix_set(J,1,1,0.0);
  return 0;
}
  
int gfn(size_t nv, const o2scl::ovector_base &x, 
	o2scl::ovector_base &y, int &pa) {
  y[0]=sin(x[1]-0.2);
  y[1]=sin(x[0]-0.25);
  return 0;
}

class cl {
public:
  
  int operator()(size_t nv, const double x[2], double y[2], int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

  int mfn(size_t nv, const o2scl::ovector_base &x, o2scl::ovector_base &y, 
	  int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }
  int mfnu(size_t nv, const o2scl::uvector_base &x, o2scl::uvector_base &y, 
	   int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }
  int mfnd(size_t nv, o2scl::ovector_base &x, o2scl::ovector_base &y,
	   o2scl::omatrix_base &j, int &pa) {
    j[0][0]=0.0;
    j[0][1]=cos(x[1]-0.2);
    j[1][0]=cos(x[0]-0.25);
    j[1][1]=0.0;
    return 0;
  }

  int mfna(size_t nv, const double x[2], double y[2], int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }

  int mfnap(size_t nv, const double *x, double *y, int &pa) {
    y[0]=sin(x[1]-0.2);
    y[1]=sin(x[0]-0.25);
    return 0;
  }
};

using namespace std;
using namespace o2scl;

int main(void) {
  cl acl;
  ovector x(2);
  uvector xu(2);
  double xa[2];
  int i;
  int vp=0;
  size_t tmp;
  int N=1000;
  int t1=0, t2=0, t3=0, t4=0, t1b=0;
  test_mgr t;
  
  cout.setf(ios::scientific);

  t.set_output_level(1);

  for(int kk=0;kk<2;kk++) {

    if (true) {
      gsl_multiroot_fdfsolver *s;
      // for tolf
      gsl_mroot_hybrids<int> cr1;
     
      int status;
      size_t iter;
     
      gsl_multiroot_function_fdf f = {&gsl_f,&gsl_df,&gsl_fdf,2,0};
      gsl_vector *gx = gsl_vector_alloc(2);
      s = gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj,2);
      
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {

	  iter=0;

	  gsl_vector_set(gx,0,0.5);
	  gsl_vector_set(gx,1,0.5);
	  gsl_multiroot_fdfsolver_set(s,&f,gx);
	  do {
	    
	    iter++;
	    status = gsl_multiroot_fdfsolver_iterate(s);
	    if(status) break;
	    status = gsl_multiroot_test_residual(s->f,cr1.tolf);
	    
	  } while (status == GSL_CONTINUE && iter < 1000);
	}
      }
      t1+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "GSL original: "
	   << (clock()-tmp)/10000 << " ";
      cout << gsl_vector_get(s->x,0) << " ";
      cout << gsl_vector_get(s->x,1) << " " << iter << endl;
      
      gsl_multiroot_fdfsolver_free(s);
      gsl_vector_free(gx);
    }

    {                                                         

      // 1. Standard approach
      mm_funct_mfptr<cl,int,ovector_base> fmf(&acl,&cl::mfn);
      gsl_mroot_hybrids<int> cr1;
    
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  x[0]=0.5;
	  x[1]=0.5;
	  cr1.msolve(2,x,vp,fmf);
	}
      }
      t1+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "Standard with ovector (gsl): "
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(x[0],0.25,1.0e-6,"1a");
      t.test_rel(x[1],0.2,1.0e-6,"1b");

      // 1. Standard approach
      mm_funct_mfptr<cl,int,uvector_base> fmfu(&acl,&cl::mfnu);
      gsl_mroot_hybrids<int,mm_funct_mfptr<cl,int,uvector_base>,
	uvector_base,uvector,uvector_alloc> cr1u;
    
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  xu[0]=0.5;
	  xu[1]=0.5;
	  cr1u.msolve(2,xu,vp,fmfu);
	}
      }
      t1+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "Standard with uvector (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(x[0],0.25,1.0e-6,"2a");
      t.test_rel(x[1],0.2,1.0e-6,"2b");

      // 1b
      jac_funct_mfptr<cl,int,ovector_base,omatrix_base> 
	fmfd(&acl,&cl::mfnd);
    
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  x[0]=0.5;
	  x[1]=0.5;
	  cr1.msolve_de(2,x,vp,fmf,fmfd);
	}
      }
      t1b+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "With analytical derivatives (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << " "
	   << cr1.last_ntrial << endl;
      t.test_rel(x[0],0.25,1.0e-6,"3a");
      t.test_rel(x[1],0.2,1.0e-6,"3b");

      // 2
      mm_vfunct_mfptr<cl,int,2> fmf2(&acl,&cl::mfna);
      gsl_mroot_hybrids<int,mm_vfunct_mfptr<cl,int,2>,double[2],
	double[2],array_alloc<double[2]> > cr2;
    
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  xa[0]=0.5;
	  xa[1]=0.5;
	  cr2.msolve(2,xa,vp,fmf2);
	}
      }
      t2+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "With arrays (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(xa[0],0.25,1.0e-6,"4a");
      t.test_rel(xa[1],0.2,1.0e-6,"4b");

      // 2b
      gsl_mroot_hybrids<int,cl,double[2],
	double[2],array_alloc<double[2]> > cr2b;
    
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  xa[0]=0.5;
	  xa[1]=0.5;
	  cr2b.msolve(2,xa,vp,acl);
	}
      }
      t2+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "Arrays w/user func. (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(xa[0],0.25,1.0e-6,"5a");
      t.test_rel(xa[1],0.2,1.0e-6,"5b");
    
      // 3 
      mm_vfunct_mfptr<cl,int,2> fmf3(&acl,&cl::mfnap);
      gsl_mroot_hybrids<int,mm_vfunct_mfptr<cl,int,2>,double *,
	double *,pointer_alloc<double> > cr3;
    
      double *xp=new double[2];
      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  xp[0]=0.5;
	  xp[1]=0.5;
	  cr3.msolve(2,xp,vp,fmf3);
	}
      }
      t2+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "With pointers (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(xa[0],0.25,1.0e-6,"6a");
      t.test_rel(xa[1],0.2,1.0e-6,"6b");
    
      // 4
      typedef int (*gfnt)(size_t, const ovector_base &, ovector_base &, 
			  int &);
      gsl_mroot_hybrids<int,gfnt,ovector_base> cr4;
      gfnt gfnv=&gfn;

      tmp=clock();
      for(int j=0;j<N;j++) {
	for(int k=0;k<N;k++) {
	  x[0]=0.5;
	  x[1]=0.5;
	  cr4.msolve(2,x,vp,gfnv);
	}
      }
      t4+=(clock()-tmp)/10000;
      cout.width(35);
      cout << "Global function pointer (gsl): " 
	   << (clock()-tmp)/10000 << " " << x[0] << " " << x[1] << endl;
      t.test_rel(x[0],0.25,1.0e-6,"7a");
      t.test_rel(x[1],0.2,1.0e-6,"7b");
      cout << endl;
    }

  }

  t.report();
  return 0;
}
