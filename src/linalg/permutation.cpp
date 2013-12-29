/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2013, Andrew W. Steiner

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

#include <o2scl/permutation.h>

using namespace std;
using namespace o2scl;

/// Output operator for permutations
std::ostream &o2scl::operator<<
  (std::ostream &os, const permutation &p) {
  if (p.size()>0) {
    for (size_t i=0;i<p.size()-1;i++) {
      os << p.get(i) << " ";
    }
    os << p.get(p.size()-1);
  }
  return os;
}
