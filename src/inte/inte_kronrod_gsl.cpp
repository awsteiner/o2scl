/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner and Jerry Gagelman
  
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
/* 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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

#include <iostream>

#include <o2scl/misc.h>
#include <o2scl/inte_kronrod_gsl.h>

using namespace std;
using namespace o2scl;

inte_workspace_gsl::inte_workspace_gsl() {
  limit=0;
}

inte_workspace_gsl::~inte_workspace_gsl() {
  if (limit>0) free();
}

void inte_workspace_gsl::make_table(table_units<> &t) {
  t.clear();
  t.line_of_names("a b r e order lev");
  t.add_constant("size",this->size);
  t.add_constant("nrmax",this->nrmax);
  t.add_constant("i",this->i);
  t.add_constant("limit",this->limit);
  t.add_constant("maximum_level",this->maximum_level);
  for(size_t j=0;j<this->size;j++) {
    vector<double> line={this->alist[j],
                         this->blist[j],
                         this->rlist[j],
                         this->elist[j],
                         ((double)this->order[j]),
                         ((double)this->level[j])};
    t.line_of_data(line.size(),line);
  }
  return;
}

int inte_workspace_gsl::allocate(size_t sz) {

  if (sz==0) {
    O2SCL_ERR2("Tried to allocate with zero size in ",
	       "inte_workspace_gsl::allocate().",exc_einval);
  }
  if (limit>0) free();

  this->alist=new double[sz];
  this->blist=new double[sz];
  this->rlist=new double[sz];
  this->elist=new double[sz];

  this->order=new size_t[sz];
  this->level=new size_t[sz];

  size=0;
  limit=sz;
  maximum_level=0;
  i=0;
  return 0;
}

int inte_workspace_gsl::free() {
  if (limit>0) {
    delete[] alist;
    delete[] blist;
    delete[] rlist;
    delete[] elist;
    delete[] order;
    delete[] level;
  }
  limit=0;
  return 0;
}

int inte_workspace_gsl::initialise(double a, double b) {
  
  size=0;
  nrmax=0;
  i=0;
  alist[0]=a;
  blist[0]=b;
  rlist[0]=0.0;
  elist[0]=0.0;
  order[0]=0;
  level[0]=0;
  
  maximum_level=0;

  return 0;
}

int inte_workspace_gsl::set_initial_result(double result, double error) {
  size=1;
  rlist[0]=result;
  elist[0]=error;
  
  return 0;
}

int inte_workspace_gsl::retrieve(double *a, double *b, double *r, 
				 double *e) const {
  *a=alist[i];
  *b=blist[i];
  *r=rlist[i];
  *e=elist[i];

  return 0;
}

int inte_workspace_gsl::qpsrt() {

  const size_t last= size - 1;
      
  double errmax;
  double errmin;
  int ii, k, top;
      
  size_t i_nrmax= nrmax;
  size_t i_maxerr=order[i_nrmax];
      
  /* Check whether the list contains more than two error estimates */
      
  if (last < 2) {
    order[0]=0;
    order[1]=1;
    ii=i_maxerr;
    return 0;
  }
      
  errmax=elist[i_maxerr];
      
  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. 
  */
  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) {
    order[i_nrmax]=order[i_nrmax - 1];
    i_nrmax--;
  }
  
  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. 
  */
  if (last < (limit/2 + 2)) {
    top=last;
  } else {
    top=limit - last + 1;
  }
  
  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). 
  */
  ii=i_nrmax + 1;
      
  /* The order of the tests in the following line is important to
     prevent a segmentation fault 
  */
  while (ii < top && errmax < elist[order[ii]]) {
    order[ii-1]=order[ii];
    ii++;
  }
  
  order[ii-1]=i_maxerr;
  
  /* Insert errmin by traversing the list bottom-up */
      
  errmin=elist[last];
  
  k=top - 1;
  
  while (k > ii - 2 && errmin >= elist[order[k]]) {
    order[k+1]=order[k];
    k--;
  }
      
  order[k+1]=last;
      
  /* Set i_max and e_max */
      
  i_maxerr=order[i_nrmax];
      
  i=i_maxerr;
  nrmax=i_nrmax;
  
  return 0;
}

double inte_workspace_gsl::sum_results() {

  const size_t n=size;
      
  size_t k;
  double result_sum=0;
      
  for (k=0; k < n; k++) {
    result_sum += rlist[k];
  }
      
  return result_sum;
}
    
int inte_workspace_gsl::subinterval_too_small(double a1, double a2, 
					      double b2) {

  const double e=std::numeric_limits<double>::epsilon();
  const double u=std::numeric_limits<double>::min();
  
  double tmp=(1 + 100 *e) *(fabs (a2) + 1000 *u);
	
  int status=fabs (a1) <= tmp && fabs (b2) <= tmp;
	
  return status;
}
    
void inte_workspace_gsl::append_interval(double a1, double b1, double area1, 
					 double error1) {
  
  const size_t i_new=size;
  
  alist[i_new]=a1;
  blist[i_new]=b1;
  rlist[i_new]=area1;
  elist[i_new]=error1;
  order[i_new]=i_new;
  level[i_new]=0;
  
  size++;
  
  return;
}

int inte_workspace_gsl::update(double a1, double b1, double area1, 
			       double error1, double a2, double b2, 
			       double area2, double error2) {
  
  const size_t i_max=i;
  const size_t i_new=size;
      
  const size_t new_level= level[i_max] + 1;
      
  /* append the newly-created intervals to the list */
      
  if (error2 > error1) {
    alist[i_max]=a2;        
    /* blist[maxerr] is already == b2 */
    rlist[i_max]=area2;
    elist[i_max]=error2;
    level[i_max]=new_level;
	  
    alist[i_new]=a1;
    blist[i_new]=b1;
    rlist[i_new]=area1;
    elist[i_new]=error1;
    level[i_new]=new_level;
  } else {
    /* alist[maxerr] is already == a1 */
    blist[i_max]=b1;        
    rlist[i_max]=area1;
    elist[i_max]=error1;
    level[i_max]=new_level;
    
    alist[i_new]=a2;
    blist[i_new]=b2;
    rlist[i_new]=area2;
    elist[i_new]=error2;
    level[i_new]=new_level;
  }
  
  size++;
      
  if (new_level >  maximum_level) {
    maximum_level=new_level;
  }
  qpsrt();

  return 0;
}

