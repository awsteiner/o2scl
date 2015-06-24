/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  hdf_file hf;
  hf.open(o2scl_settings.get_data_dir()+"/rmfdata/FSUGold.o2");
  hf.setd("gs",sqrt(99.4266));
  hf.setd("gw",sqrt(169.8349));
  hf.setd("gr",sqrt(184.6877)/2.0);
  hf.seti("oakstyle",1);
  hf.sets("reference",((string)"F.J. Fattoyev, C.J. Horowitz,")+
	  "J. Piekarewicz, and G. Shen, Phys. Rev. C 82 (2010) 055803.");
  hf.setd("zeta",0.03);
  hf.setd("g2",0.0);
  hf.setd("g3",0.0);
  hf.setd("b1",0.0);

  t.report();

  return 0;
}


