/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2017, Andrew W. Steiner
  
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
#include <o2scl/mcmc_mpi.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_OLDER_COMPILER

int main(void) {
  test_mgr t;
  t.report();
  return 0;
}

#else

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(1);

  // ----------------------------------------------------------------

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  typedef std::function<int(size_t,const ubvector &,double &,
			    std::array<double,1> &)> point_funct;
  
  typedef std::function<int(const ubvector &,double,size_t,bool,
			    std::array<double,1> &)> measure_funct;
  
  typedef std::function<int(const ubvector &,double,std::vector<double> &,
			    std::array<double,1> &)> fill_funct;
  
  mcmc_cli<point_funct,fill_funct,std::array<double,1>,ubvector> mcc;
  
  mcc.run(argc,argv);

  // ----------------------------------------------------------------

  tm.report();
  
  return 0;
}
#endif
