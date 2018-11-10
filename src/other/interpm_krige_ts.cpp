/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#include <o2scl/interpm_krige.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef o2scl::matrix_view_vec_vec<ubvector> mat_t;
typedef vector<function<double(const matrix_row_gen<mat_t> &,
			       const matrix_row_gen<mat_t> &) > > f1_t;
typedef vector<function<double(const matrix_row_gen<mat_t> &,
			       const ubvector &) > > f2_t;

template<class vec_t, class vec2_t>
double covar(const vec_t &x, const vec2_t &y, double len) {
  double ret=exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0);
  cout << len << " " << x[0] << " " << y[0] << " "
       << x[1] << " " << y[1] << " " << ret << endl;
  return ret;
}

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  {
    // Construct the data
    vector<ubvector> x;
    ubvector tmp(2);
    tmp[0]=1.04; tmp[1]=0.02;
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.01; 
    x.push_back(tmp);
    tmp[0]=0.81; tmp[1]=0.23; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.83; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.99; 
    x.push_back(tmp);
    tmp[0]=0.82; tmp[1]=0.84; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.24; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.02; 
    x.push_back(tmp);
    mat_t x2(x);

    vector<ubvector> y;
    tmp.resize(8);
    for(size_t i=0;i<8;i++) {
      tmp[i]=ft(x[i][0],x[i][1]);
    }
    y.push_back(tmp);
    mat_t y2(y);

    interpm_krige<ubvector,mat_t,matrix_column_gen<mat_t>,
		  matrix_row_gen<mat_t> > ik;
    ik.verbose=2;
    f1_t fa1={std::bind(&covar<matrix_row_gen<mat_t>,matrix_row_gen<mat_t> >,
			std::placeholders::_1,std::placeholders::_2,0.70)};
    f2_t fa2={std::bind(&covar<matrix_row_gen<mat_t>,ubvector>,
			std::placeholders::_1,std::placeholders::_2,0.70)};
    ik.set_data<f1_t>(2,1,8,x2,y2,fa1);
  
    cout << "Interpolate at a point and compare the three methods:" << endl;
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
  }

  {
    // Construct the data
    vector<ubvector> x;
    ubvector tmp(2);
    tmp[0]=1.04; tmp[1]=0.02;
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.01; 
    x.push_back(tmp);
    tmp[0]=0.81; tmp[1]=0.23; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.83; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.99; 
    x.push_back(tmp);
    tmp[0]=0.82; tmp[1]=0.84; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.24; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.02; 
    x.push_back(tmp);
    mat_t x2(x);

    vector<ubvector> y;
    tmp.resize(8);
    for(size_t i=0;i<8;i++) {
      tmp[i]=ft(x[i][0],x[i][1]);
    }
    y.push_back(tmp);
    mat_t y2(y);

    interpm_krige<ubvector,mat_t,matrix_column_gen<mat_t>,
		  matrix_row_gen<mat_t> > ik;
    ik.verbose=2;
    f1_t fa1={std::bind(&covar<matrix_row_gen<mat_t>,matrix_row_gen<mat_t> >,
			std::placeholders::_1,std::placeholders::_2,0.70)};
    f2_t fa2={std::bind(&covar<matrix_row_gen<mat_t>,ubvector>,
			std::placeholders::_1,std::placeholders::_2,0.70)};
    ik.set_data<f1_t>(2,1,8,x2,y2,fa1,true);
  
    cout << "Interpolate at a point and compare the three methods:" << endl;
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
  }

  {
    // Construct the data
    vector<ubvector> x;
    ubvector tmp(2);
    tmp[0]=1.04; tmp[1]=0.02;
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.01; 
    x.push_back(tmp);
    tmp[0]=0.81; tmp[1]=0.23; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.83; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.99; 
    x.push_back(tmp);
    tmp[0]=0.82; tmp[1]=0.84; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.24; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.02; 
    x.push_back(tmp);
    mat_t x2(x);

    vector<ubvector> y;
    tmp.resize(8);
    for(size_t i=0;i<8;i++) {
      tmp[i]=ft(x[i][0],x[i][1]);
    }
    y.push_back(tmp);
    mat_t y2(y);

    interpm_krige_optim<ubvector,mat_t,matrix_column_gen<mat_t>,
			matrix_row_gen<mat_t> > iko;
    iko.verbose=2;
    
    iko.set_data(2,1,8,x2,y2,true);
    
    cout << "Interpolate at a point and compare the three methods:" << endl;
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
  }

  /*
  if (false) {
    // Construct the data
    vector<ubvector> x;
    ubvector tmp(2);
    tmp[0]=1.04; tmp[1]=0.02;
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.01; 
    x.push_back(tmp);
    tmp[0]=0.81; tmp[1]=0.23; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.83; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.99; 
    x.push_back(tmp);
    tmp[0]=0.82; tmp[1]=0.84; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.24; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.02; 
    x.push_back(tmp);
    mat_t x2(x);

    vector<ubvector> y;
    tmp.resize(8);
    for(size_t i=0;i<8;i++) {
      tmp[i]=ft(x[i][0],x[i][1]);
    }
    y.push_back(tmp);
    mat_t y2(y);

    interpm_krige_nn<ubvector,mat_t,matrix_column_gen<mat_t>,
		     matrix_row_gen<mat_t> > ik;
    vector<function<double(const ubvector &,const ubvector &)> > fa={covar};
    ik.set_data(2,1,8,x2,y2,fa,4);
  
    cout << "Interpolate at a point and compare the three methods:" << endl;
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;

  }
  */

  t.report();
  return 0;
}

