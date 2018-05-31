/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2018, Andrew W. Steiner
  
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

/* Example: ex_interp2_nogrid.cpp
   -------------------------------------------------------------------
   A simple example for two-dimensional interpolation using
   data not specified on a grid
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/rng_gsl.h>
#include <o2scl/interp2_planar.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/test_mgr.h>
#include <o2scl/prob_dens_mdim_amr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

// A function for filling the data and comparing results
double f(double x, double y) {
  return sin(x)*sin(1.0/(y+0.1));
}

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  rng_gsl rg;
  rg.set_seed(10);

  cout << "        "
       << "i2n        "
       << "i2p        "
       << "imw1       "
       << "imw3       "
       << "imw5       "
       << "amr1       "
       << "amr2       " << endl;
  
  for(size_t k=0;k<2;k++) {
    
    // Create the sample data
    size_t N;
    if (k==0) N=100;
    else N=1000;

    // Sample data
    table<> tdata;
    tdata.line_of_names("x y z");
    for(size_t i=0;i<N;i++) {
      double dx=rg.random()*4.0;
      double dy=rg.random()*4.0;
      double line[3]={dx,dy,f(dx,dy)};
      tdata.line_of_data(3,line);
    }

    // Full output table
    table3d t3d;
    t3d.set_xy("x",uniform_grid_end<double>(0.0,4.0,99),
	       "y",uniform_grid_end<double>(0.0,4.0,99));
    t3d.line_of_names((std::string)("exact interp2_neigh interp2_planar ")+
		      "interpm_idw1 interpm_idw3 "+
		      "interpm_idw5 prob_dens_mdim_amr1 prob_dens_mdim_amr2");
  
    // Various rearrangments of the data

    const double *dx=&(tdata["x"][0]);
    const double *dy=&(tdata["y"][0]);
    const double *dz=&(tdata["z"][0]);
    vector<const double *> darr(3);
    std::vector<double> low={0.0,0.0};
    std::vector<double> high={4.0,4.0};
    matrix_view_table<std::vector<double> > mvt(tdata,{"x","y","z"});

    // Initialize interpolation objects
    
    interp2_neigh<const double *> i2n;
    i2n.set_data(N,dx,dy,dz);

    interp2_planar<const double *> i2p;
    i2p.set_data(N,dx,dy,dz);

    interpm_idw<const double *> imw1, imw3, imw5;
    darr={&(tdata["x"][0]),&(tdata["y"][0]),&(tdata["z"][0])};
    imw1.set_data(2,1,N,darr);
    darr={&(tdata["x"][0]),&(tdata["y"][0]),&(tdata["z"][0])};
    imw3.set_data(2,1,N,darr);
    darr={&(tdata["x"][0]),&(tdata["y"][0]),&(tdata["z"][0])};
    imw5.set_data(2,1,N,darr);
    imw1.set_points(1);
    imw3.set_points(3);
    imw5.set_points(5);

    interpm_krige<> imk;

    prob_dens_mdim_amr<> pdma1;
    pdma1.set(low,high);
    pdma1.initial_parse_new(mvt);

    prob_dens_mdim_amr<> pdma2;
    pdma2.dim_choice=prob_dens_mdim_amr<>::random;
    pdma2.set(low,high);
    pdma2.initial_parse(mvt);

    static const size_t n_methods=7;
    double qual[n_methods];
    for(size_t i=0;i<n_methods;i++) qual[i]=0.0;

    for(size_t i=0;i<t3d.get_nx();i++) {
      for(size_t j=0;j<t3d.get_ny();j++) {
	double x=t3d.get_grid_x(i), y=t3d.get_grid_y(j);
	std::vector<double> arr(2);
	arr[0]=x;
	arr[1]=y;
	double v_exact=f(x,y);
	double v_i2n=i2n.eval(x,y);
	double v_i2p=i2p.eval(x,y);
	double v_imw1=imw1.eval(arr);
	double v_imw3=imw3.eval(arr);
	double v_imw5=imw5.eval(arr);
	double v_pdma1=pdma1.pdf(arr);
	double v_pdma2=pdma2.pdf(arr);
	qual[0]+=fabs(v_i2n-v_exact);
	qual[1]+=fabs(v_i2p-v_exact);
	qual[2]+=fabs(v_imw1-v_exact);
	qual[3]+=fabs(v_imw3-v_exact);
	qual[4]+=fabs(v_imw5-v_exact);
	qual[5]+=fabs(v_pdma1-v_exact);
	qual[6]+=fabs(v_pdma2-v_exact);
	t3d.set(i,j,"exact",v_exact);
	t3d.set(i,j,"interp2_neigh",v_i2n);
	t3d.set(i,j,"interp2_planar",v_i2p);
	t3d.set(i,j,"interpm_idw1",v_imw1);
	t3d.set(i,j,"interpm_idw3",v_imw3);
	t3d.set(i,j,"interpm_idw5",v_imw5);
	t3d.set(i,j,"prob_dens_mdim_amr1",v_pdma1);
	t3d.set(i,j,"prob_dens_mdim_amr2",v_pdma2);
      }
    }

    if (k==0) cout << "N=100:  ";
    else cout << "N=1000: ";
    cout.precision(4);
    vector_out(cout,n_methods,qual,true);
    cout.precision(6);
    
    hdf_file hf;
    if (k==0) {
      hf.open_or_create("ex_interp2_nogrid_100.o2");
    } else {
      hf.open_or_create("ex_interp2_nogrid_1000.o2");
    }
    hdf_output(hf,tdata,"tdata");
    hdf_output(hf,t3d,"t3d");
    hf.close();
    
  }
  
  t.report();
  return 0;
}
// End of example
