/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Gerard
 * Jungman, Brian Gough
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
 * 
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/series_acc.h>

using namespace std;
using namespace o2scl;

series_acc::series_acc(size_t usize) {
  if (usize>0) {
    w=gsl_sum_levin_u_alloc(usize);
    wt=gsl_sum_levin_utrunc_alloc(usize);
  }
  size=usize;
}

void series_acc::set_size(size_t new_size) {
  if (size>0) {
    gsl_sum_levin_u_free(w);
    gsl_sum_levin_utrunc_free(wt);
  }
  size=new_size;
  w=gsl_sum_levin_u_alloc(size);
  wt=gsl_sum_levin_utrunc_alloc(size);
  return;
}

series_acc::~series_acc() {
  if (size>0) {
    gsl_sum_levin_u_free(w);
    gsl_sum_levin_utrunc_free(wt);
  }
}

int series_acc::levin_u_step(const double term, const size_t n, 
			     const size_t nmax, double &sum_accel) {
  
  if (n==0) {

    sum_accel=term;
    w->sum_plain=term;

    w->q_den[0]=1.0/term;
    w->q_num[0]=1.0;

    w->dq_den[series_index(0,0,nmax)]=-1.0/(term*term);
    w->dq_num[series_index(0,0,nmax)]=0.0;

    w->dsum[0]=1.0;

  } else {

    double result;
    double factor=1.0;
    double ratio=(double) n/(n+1.0);
    unsigned int i;
    int j;

    w->sum_plain += term;

    w->q_den[n]=1.0/(term*(n+1.0)*(n+1.0));
    w->q_num[n]=w->sum_plain*w->q_den[n];

    for (i=0; i < n; i++) {
      w->dq_den[series_index(i,n,nmax)]=0;
      w->dq_num[series_index(i,n,nmax)]=w->q_den[n];
    }

    w->dq_den[series_index(n,n,nmax)]=-w->q_den[n]/term;
    w->dq_num[series_index(n,n,nmax)]=
      w->q_den[n]+w->sum_plain*(w->dq_den[series_index(n,n,nmax)]);

    for (j=n-1; j >= 0; j--) {
      double c=factor*(j+1)/(n+1);
      factor*=ratio;
      w->q_den[j]=w->q_den[j+1]-c*w->q_den[j];
      w->q_num[j]=w->q_num[j+1]-c*w->q_num[j];

      for (i=0; i < n; i++) {
	w->dq_den[series_index(i,j,nmax)]=
	  w->dq_den[series_index(i,j+1,nmax)]-
	  c*w->dq_den[series_index(i,j,nmax)];
	w->dq_num[series_index(i,j,nmax)]=
	  w->dq_num[series_index(i,j+1,nmax)]-
	  c*w->dq_num[series_index(i,j,nmax)];
      }

      w->dq_den[series_index(n,j,nmax)]=
	w->dq_den[series_index(n,j+1,nmax)];
      w->dq_num[series_index(n,j,nmax)]=
	w->dq_num[series_index(n,j+1,nmax)];
    }

    result=w->q_num[0]/w->q_den[0];

    sum_accel=result;

    for (i=0; i <= n; i++) {
      w->dsum[i]=(w->dq_num[series_index(i,0,nmax)]-
		  result*w->dq_den[series_index(i,0,nmax)])/w->q_den[0];
    }
	
  }

  return success;
}
    
int series_acc::levin_utrunc_step(const double term, const size_t n, 
				  double &sum_accel) {
  
  if (term == 0.0) {

    /* This is actually harmless when treated in this way. A term
       which is exactly zero is simply ignored; the state is not
       changed. We return GSL_EZERODIV as an indicator that this
       occured. */
    return exc_ezerodiv;

  } else if (n == 0) {

    sum_accel=term;
    wt->sum_plain=term;
    wt->q_den[0]=1.0/term;
    wt->q_num[0]=1.0;

  } else {

    double factor=1.0;
    double ratio=(double) n/(n+1.0);
    int j;
	
    wt->sum_plain += term;
    wt->q_den[n]=1.0/(term*(n+1.0)*(n+1.0));
    wt->q_num[n]=wt->sum_plain*wt->q_den[n];
	
    for (j=n-1; j >= 0; j--) {
      double c=factor*(j+1)/(n+1);
      factor *= ratio;
      wt->q_den[j]=wt->q_den[j+1]-c*wt->q_den[j];
      wt->q_num[j]=wt->q_num[j+1]-c*wt->q_num[j];
    }
	
    sum_accel=wt->q_num[0]/wt->q_den[0];
  }

  return success;
}
