/*
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
*/
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/emulator.h>
#include <o2scl/hdf_io.h>

#ifdef O2SCL_ARMA
#include <armadillo>
#endif
#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y;
}

void generate_table(table<> &tab, size_t N=100) {

  tab.clear();
  
  tab.line_of_names("x y z d");
  
  for(size_t i=0;i<N;i++) {
    double x=3.0*sin(i*i);
    double y=5.0*cos(pow(i,4.0));
    vector<double> line={x,y,ft(x,y),ft(x,y)*2.0};
    tab.line_of_data(4,line);
  }
  
  return;
}

int main(void) {

  test_mgr t;
  t.set_output_level(2);
  
  cout.setf(ios::scientific);

#ifdef O2SCL_PYTHON
  
    if (true) {
    
      vector<string> col_list={"x","y","z","d"};
      
      o2scl_settings.py_init();
      o2scl_settings.add_python_path(o2scl_settings.get_data_dir()+"python/");
  
      table<> tab;
      generate_table(tab);

      hdf_file hf;
      hf.open_or_create("emu_data.o2");
      hdf_output(hf,tab,"tab");
      hf.close();

      {
        cout << "Emulator from emu_sklearn.py" << endl;
        emulator_python<ubvector,ubvector> em1;
        em1.verbose=3;
        em1.set("emu_sklearn","emu_py","train","point",2,
                "emu_data.o2",0,col_list,false);
    
        for(size_t j=0;j<10;j++) {
      
          ubvector p(2);
          p[0]=1.0+sin(j*4);
          p[1]=2.0+sin(j*5);
          double z;
          double dz;
          ubvector dat(2), datu(2);
          em1.eval_unc(2,p,z,dz,dat,datu);
          cout << ft(p[0],p[1]) << " " << z << endl;
      
        }
        cout << endl;

      }
      
      if (false) {
        cout << "Emulator from emu_sklearn2.py" << endl;
        emulator_python<ubvector,ubvector> em2;
        em2.verbose=3;
        em2.set("emu_sklearn2","emu_gpr","read_data","predict",2,
                "emu_data.o2",0,col_list,false);
    
        for(size_t j=0;j<10;j++) {
      
          ubvector p(2);
          p[0]=1.0+sin(j*4);
          p[1]=2.0+sin(j*5);
          double z;
          double dz;
          ubvector dat(2), datu(2);
          em2.eval_unc(2,p,z,dz,dat,datu);
          cout << ft(p[0],p[1]) << " " << z << endl;
      
        }
      }
      
      if (true) {
        cout << "Emulator from emu_tf.py" << endl;
        emulator_python<ubvector,ubvector> em3;
        em3.verbose=3;
        em3.set("emu_tf","emu_dnn","read_data","predict",2,
                "emu_data.o2",0,col_list,false);
    
        for(size_t j=0;j<10;j++) {
      
          ubvector p(2);
          p[0]=1.0+sin(j*4);
          p[1]=2.0+sin(j*5);
          double z;
          double dz;
          ubvector dat(2), datu(2);
          em3.eval_unc(2,p,z,dz,dat,datu);
          cout << ft(p[0],p[1]) << " " << z << endl;
      
        }

      }
      
      o2scl_settings.py_final();
    }

#endif
  
  t.report();

  return 0;
}

