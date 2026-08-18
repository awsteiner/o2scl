"""
  -------------------------------------------------------------------
  
  Copyright (C) 2022-2024, Andrew W. Steiner
  
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
"""

# This python code is currently used by src/base/funct_ts to test the
# C++ interface to Python functions.

import numpy

def fun(x):
    return x*numpy.pi

class c:

    def fun2(self,x):
        return x*x*numpy.pi
