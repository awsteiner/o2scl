/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2024, Andrew W. Steiner

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

#include <o2scl/givens.h>

using namespace std;
using namespace o2scl;

void o2scl_linalg::create_givens(const double a, const double b, 
				 double &c, double &s) {
  if (b==0) {
    c=1;
    s=0;
  } else if (fabs(b)>fabs(a)) {
    double t=-a/b;
    double s1=1.0/sqrt(1+t*t);
    s=s1;
    c=s1*t;
  } else {
    double t=-b/a;
    double c1=1.0/sqrt(1+t*t);
    c=c1;
    s=c1*t;
  }

  return;
}

