 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2023-2024, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/nflows_python.h>
#include <o2scl/set_python.h>
#include <o2scl/rng.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

#ifdef O2SCL_SET_PYTHON
  
  // Construct the data
  static const size_t N=1000;

  rng<> r;
  
  table<> tab;
  tab.line_of_names("x y");
  for(size_t i=0;i<N;i++) {
    if (i%2==0) {
      double x=0.7+0.1*r.random();
      double y=0.5+0.1*r.random();
      vector<double> line={x,y};
      tab.line_of_data(line.size(),line);
    } else {
      double x=-0.7+0.1*r.random();
      double y=0.2+0.1*r.random();
      vector<double> line={x,y};
      tab.line_of_data(line.size(),line);
    }
  }
  
  hdf_file hf;
  hf.open_or_create("nflows_python_data.o2");
  hdf_output(hf,tab,"tab");
  hf.close();
  
  if (true) {

    cout << "-------------------------------------------------"
         << endl;
    
    tensor<> tin;
    vector<size_t> in_size={N,2};
    tin.resize(2,in_size);
    for(size_t j=0;j<N;j++) {
      vector<size_t> ix;
      ix={j,0};
      tin.get(ix)=tab.get(0,j);
      ix={j,1};
      tin.get(ix)=tab.get(1,j);
    }

    nflows_python<> nfp("o2sclpy",tin,"max_iter=200,verbose=0",
                        "nflows_nsf",0);
    cout << endl;

    // Try sampling
    vector<double> x(2);
    for(size_t j=0;j<10;j++) {
      nfp(x);
      cout << j << " ";
      cout.setf(ios::showpos);
      cout << x[0] << " " << x[1] << endl;
      cout.unsetf(ios::showpos);
    }
    cout << endl;
    
    // Output the log pdf
    x[0]=0.7;
    x[1]=0.5;
    cout << nfp.log_pdf(x) << endl;
    x[0]=0.0;
    x[1]=0.0;
    cout << nfp.log_pdf(x) << endl;
    x[0]=-0.7;
    x[1]=0.2;
    cout << nfp.log_pdf(x) << endl;
    
  }
    
#endif
    
  t.report();
  return 0;
}

