/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2014, Andrew W. Steiner
  
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

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

#include <o2scl/prob_dens_func.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

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

  prob_dens_gaussian(1.0,2.0);
  prob_dens_lognormal(1.0,2.0);
  prob_dens_uniform(1.0,2.0);
  
  {
    hist h, h2;
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
      double x=p.sample();
      h2.update(x);
    }
    
    for(size_t i=0;i<10;i++) {
      cout << i << " " << h2[i] << endl;
    }
    cout << endl;

  }

  t.report();

  return 0;
}

