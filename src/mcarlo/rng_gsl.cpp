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
/* rng/random.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// For time(0)
#include <ctime> 

#include <o2scl/rng_gsl.h>
#include <o2scl/err_hnd.h>

using namespace std;
using namespace o2scl;

rng_gsl::rng_gsl(const gsl_rng_type *gtype) {
  rng=gtype;
  gr=gsl_rng_alloc(gtype);
  gsl_rng_set(gr,0);
  seed=time(0);
}

rng_gsl::rng_gsl(unsigned long int lseed, const gsl_rng_type *gtype) {
  rng=gtype;
  gr=gsl_rng_alloc(gtype);
  seed=lseed;
  gsl_rng_set(gr,seed);
}

rng_gsl::~rng_gsl() {
  gsl_rng_free(gr);
}

unsigned long int rng_gsl::random_int(unsigned long int n) {
  unsigned long int offset = gr->type->min;
  unsigned long int range = gr->type->max - offset;
  unsigned long int scale = range / n;
  unsigned long int k;
  
  if (n > range) {
    O2SCL_ERR("n exceeds maximum value of generator",GSL_EINVAL);
    return 0;
  }
  
  do {
    k = (((gr->type->get) (gr->state)) - offset) / scale;
  } while (k >= n);
  
  return k;
}
