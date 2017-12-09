/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <random>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

#include <o2scl/prob_dens_func.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mcarlo_plain.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

prob_dens_mdim_gaussian<> pdmg;

double test_fun(size_t nv, const ubvector &x) {
  //cout << nv << " " << x[0] << " " << x[1] << " "
  //<< pdmg.pdf(x) << endl;
  //char ch;
  //cin >> ch;
  return pdmg.pdf(x);
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  {
    const gsl_rng_type *T=gsl_rng_default;
    gsl_rng *r=gsl_rng_alloc(T);
    gsl_rng_env_setup();
    
    gsl_histogram *h=gsl_histogram_alloc(10);
    gsl_histogram_set_ranges_uniform(h,0.0,1.0);
    gsl_histogram *h2=gsl_histogram_alloc(10);
    gsl_histogram_set_ranges_uniform(h2,0.0,1.0);

    for (size_t i=0;i<10000;i++) {
      gsl_histogram_accumulate(h,sin((double)i),1);
    }

    gsl_histogram_pdf *p=gsl_histogram_pdf_alloc(h->n);
    gsl_histogram_pdf_init(p,h);

    for(size_t i=0;i<11;i++) {
      cout.width(2);
      cout << i << " " << p->range[i] << " " << p->sum[i] << endl;
    }
    cout << endl;

    for (size_t i=0;i<10000;i++) {
      double u=gsl_rng_uniform(r);
      double x=gsl_histogram_pdf_sample(p,u);
      gsl_histogram_accumulate(h2,x,1);
    }
    
    gsl_histogram_pdf_free(p);
    gsl_histogram_free(h2);
    gsl_histogram_free(h);
    gsl_rng_free(r);
  }

  hist h2;
  {
    hist h;

    uniform_grid_end<> ug(0.0,1.0,10);
    h.set_bin_edges(ug);
    h2.set_bin_edges(ug);
    
    for(size_t i=0;i<10000;i++) {
      if (sin((double)i)>=0.0) {
	h.update(sin((double)i));
      }
    }
    
    prob_dens_hist p;
    p.init(h);
    
    for(size_t i=0;i<10000;i++) {
      double x=p();
      h2.update(x);
    }
    
    for(size_t i=0;i<10;i++) {
      cout << i << " " << h2[i] << endl;
    }
    cout << endl;
  }

  hist h4;
  {
    hist h3;

    std::random_device rd;
    std::normal_distribution<double> normdist(1.0,2.0);
    
    uniform_grid_end<> ug(0.0,1.0,10);
    h3.set_bin_edges(ug);
    h4.set_bin_edges(ug);

    for(size_t i=0;i<10000;i++) {
      if (sin((double)i)>=0.0) {
	h3.update(sin((double)i));
      }
    }

    prob_dens_hist p2;
    p2.init(h3);
    
    for(size_t i=0;i<10000;i++) {
      double x=p2();
      h4.update(x);
    }
    
    for(size_t i=0;i<10;i++) {
      cout << i << " " << h4[i] << endl;
    }
    cout << endl;

  }

  for(size_t i=0;i<10;i++) {
    t.test_rel(h2[i],h4[i],1.0e-10,"o2scl vs std::normal");
  }

  ubvector cent(2);
  ubmatrix covar(2,2);
  cent[0]=2.0;
  cent[1]=3.0;
  covar(0,0)=1.0;
  covar(1,1)=4.0;
  covar(0,1)=-1.0;
  covar(1,0)=-1.0;
  pdmg.set(2,cent,covar);

  // Test the gaussian PDF normalization
  {
    mcarlo_plain<> gm;
    
    ubvector a(2), b(2);
    a[0]=-15.0;
    a[1]=-15.0;
    b[0]=15.0;
    b[1]=15.0;
    
    multi_funct tf=test_fun;
    
    gm.n_points=100000;
    double res, err;
    gm.minteg_err(tf,2,a,b,res,err);
    
    cout << "O2scl res,err,rel: " 
	 << res << " " << err << endl;
    t.test_rel(res,1.0,err,"normalization");
  }

  {
    ubvector xx(2);
    prob_dens_mdim_biv_gaussian<ubvector> biv;
    biv.set(2.0,3.0,1.0,0.5,-0.5);
    xx[0]=2.0;
    xx[1]=3.0;
    cout << biv.pdf(xx) << endl;
    biv.contour(0.3,0.0,xx);
    cout << xx[0] << " " << xx[1] << endl;
    cout << biv.pdf(xx) << endl;
    t.test_rel(biv.pdf(xx),0.3,1.0e-6,"biv contour");
  }
  
  t.report();

  return 0;
}

