/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <o2scl/contour.h>
#include <o2scl/test_mgr.h>
#include <o2scl/rng.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double fun(double x, double y) {
  return 15.0*exp(-pow(x-20.0,2.0)/400.0-pow(y-5.0,2.0)/25.0)
    +40.0*exp(-pow(x-70.0,2.0)/4900.0-pow(y-2.0,2.0)/4.0);
}

double fun2(double x, double y) {
  return 2.0*pow(x-4.0,2.0)+6.0*pow(y-2.0,2.0);
}

template<class vec_t, class mat_t> 
int print_data_xhoriz(int nx, int ny, vec_t &x, vec_t &y, mat_t &data) {
  int j, k;

  cout.setf(ios::showpos);
  cout.precision(3);

  for(k=ny-1;k>=0;k--) {
    string stx="y";
    stx+=itos(k);
    cout.width(3);
    cout << stx << " ";
    cout << y[k] << " ";
    for(j=0;j<nx;j++) {
      cout << data(j,k) << " ";
    }
    cout << endl;
  }

  // X grid
  cout << "               ";
  for(j=0;j<nx;j++) {
    cout << x[j] << " ";
  }
  cout << endl;

  // X labels
  cout << "              ";
  for(j=0;j<nx;j++) {
    string stx="x";
    stx+=itos(j);
    cout.width(11);
    cout << stx;
  }
  cout << endl;

  cout.unsetf(ios::showpos);
  cout.precision(6);
  return 0;
}

template<class vec_t, class mat_t> 
int print_data_yhoriz(int nx, int ny, vec_t &x, vec_t &y, mat_t &data) {
  int j, k;

  cout.setf(ios::showpos);
  cout.precision(3);

  for(j=nx-1;j>=0;j--) {
    string stx="x";
    stx+=itos(j);
    cout.width(3);
    cout << stx << " ";
    cout << x[j] << " ";
    for(k=0;k<ny;k++) {
      cout << data(j,k) << " ";
    }
    cout << endl;
  }

  // Y grid
  cout << "               ";
  for(k=0;k<ny;k++) {
    cout << y[k] << " ";
  }
  cout << endl;

  // Y labels
  cout << "              ";
  for(k=0;k<ny;k++) {
    string stx="y";
    stx+=itos(k);
    cout.width(11);
    cout << stx;
  }
  cout << endl;

  cout.unsetf(ios::showpos);
  cout.precision(6);
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // This file is used by ../internal/contour_gr.cpp to generate
  // the plots for the documentation
  ofstream fout("contour.dat");
  
  // Initialization

  ubvector x(10), y(10), xs(5), ys(5);
  ubmatrix data(10,10);
  ubmatrix datas(5,5);
  int i, j, nsx, nsy, k;

  contour co;

  // Create initial array

  for(i=0;i<5;i++) {
    ys[i]=((double)i*2);
    xs[i]=ys[i]*ys[i];
  }

  for(j=0;j<5;j++) {
    for(i=0;i<5;i++) {
      datas(j,i)=15.0*exp(-pow(xs[i]-20.0,2.0)/400.0-
			  pow(ys[j]-5.0,2.0)/25.0);
      datas(j,i)+=40.0*exp(-pow(xs[i]-70.0,2.0)/4900.0-
			   pow(ys[j]-2.0,2.0)/4.0);
    }
  }

  co.set_data(5,5,xs,ys,datas);
  print_data_xhoriz(5,5,xs,ys,datas);
  
  // Regrid
  cout << "In regrid: " << endl;
  co.regrid_data(2,2);
  cout << "Done with regrid." << endl;
  
  size_t ngx=0, ngy=0;
  ubvector *xg=0, *yg=0;
  ubmatrix *datag=0;
  co.get_data(ngx,ngy,xg,yg,datag);
  
  print_data_xhoriz(ngx,ngy,*xg,*yg,*datag);

  // Now test contour lines

  for(i=0;i<10;i++) {
    y[i]=((double)i);
    x[i]=y[i]*y[i];
  }
  
  for(j=0;j<10;j++) {
    for(k=0;k<10;k++) {
      data(j,k)=fun(x[j],y[k]);
    }
  }
  
  cout << "Square 1:" << endl;
  print_data_xhoriz(10,10,x,y,data);

  ubvector levels(7);
  levels[0]=5.0;
  levels[1]=10.0;
  levels[2]=15.0;
  levels[3]=20.0;
  levels[4]=25.0;
  levels[5]=30.0;
  levels[6]=35.0;
  
  int nc;

  co.set_data(10,10,x,y,data);
  co.set_levels(7,levels);
  vector<contour_line> conts;
  co.calc_contours(conts);
  nc=conts.size();
  
  fout << nc << endl;
  for(i=0;i<nc;i++) {
    cout << "Contour " << i << " at level " 
	 << conts[i].level << ":" << endl;
    int cs=conts[i].x.size();
    fout << cs << " ";
    for(j=0;j<cs;j++) {
      cout << "(" << conts[i].x[j] << ", " 
	   << conts[i].y[j] << ")" << endl;
      fout << conts[i].x[j] << " " << conts[i].y[j] << " ";
      t.test_rel(fun(conts[i].x[j],conts[i].y[j]),conts[i].level,
		 1.5e-1,"curve");
    }
    fout << endl;
  }
  cout << endl;
  
  // ------------------------------------------------------------

  cout << "Non-square 1:" << endl;

  // Test non-square data
  ubvector sqx(10), sqy(8), sqlev(4);
  ubmatrix sqd(10,8);

  for(i=0;i<10;i++) sqx[i]=i*(8.0/9.0);
  for(i=0;i<8;i++) sqy[i]=i*(3.0/7.0);
  for(i=0;i<8;i++) {
    for(j=0;j<10;j++) {
      sqd(j,i)=2.0*pow(sqx[j]-4.0,2.0)+6.0*pow(sqy[i]-2.0,2.0);
    }
  }
  sqlev[0]=4.0;
  sqlev[1]=10.0;
  sqlev[2]=20.0;
  sqlev[3]=40.0;

  print_data_xhoriz(10,8,sqx,sqy,sqd);

  co.set_data(10,8,sqx,sqy,sqd);
  co.set_levels(4,sqlev);
  vector<contour_line> conts2;
  co.calc_contours(conts2);
  nc=conts2.size();

  fout << nc << endl;
  for(i=0;i<nc;i++) {
    cout << "Contour " << i << " at level " 
	 << conts2[i].level << ":" << endl;
    int cs=conts2[i].x.size();
    fout << cs << " ";
    for(j=0;j<cs;j++) {
      cout << "(" << conts2[i].x[j] << ", " 
	   << conts2[i].y[j] << ")" << endl;
      fout << conts2[i].x[j] << " " << conts2[i].y[j] << " ";
      t.test_rel(fun2(conts2[i].x[j],conts2[i].y[j]),conts2[i].level,
		 1.5e-1,"curve");
    }
    fout << endl;
  }
  cout << endl;
  
  // ------------------------------------------------------------

  cout << "Non-square 2:" << endl;
    
  // Test non-square data
  ubvector srx(8), sry(10), srlev(4);
  ubmatrix srd(8,10);
    
  for(i=0;i<8;i++) srx[i]=i*(8.0/7.0);
  for(i=0;i<10;i++) sry[i]=i*(3.0/9.0);
  for(i=0;i<10;i++) {
    for(j=0;j<8;j++) {
      srd(j,i)=2.0*pow(srx[j]-4.0,2.0)+6.0*pow(sry[i]-2.0,2.0);
    }
  }
  srlev[0]=4.0;
  srlev[1]=10.0;
  srlev[2]=20.0;
  srlev[3]=40.0;

  co.set_data(8,10,srx,sry,srd);
  co.set_levels(4,srlev);
  vector<contour_line> conts3;
  co.calc_contours(conts3);
  nc=conts3.size();
  
  vector<edge_crossings> xed, yed;
  co.get_edges(xed,yed);
  
  print_data_xhoriz(8,10,srx,sry,srd);

  cout << "Print edges: " << endl;
  for(size_t ir=0;ir<yed.size();ir++) {
    cout << "Level: " << srlev[ir] << endl;
    co.print_edges_xhoriz(xed[ir],yed[ir]);
  }
  cout << endl;

  print_data_yhoriz(8,10,srx,sry,srd);

  cout << "Print edges: " << endl;
  for(size_t ir=0;ir<yed.size();ir++) {
    cout << "Level: " << srlev[ir] << endl;
    co.print_edges_yhoriz(xed[ir],yed[ir]);
  }
  cout << endl;

  fout << nc << endl;
  for(i=0;i<nc;i++) {
    cout << "Contour " << i << " at level " 
	 << conts3[i].level << ":" << endl;
    int cs=conts3[i].x.size();
    fout << cs << " ";
    for(j=0;j<cs;j++) {
      cout << "(" << conts3[i].x[j] << ", " 
	   << conts3[i].y[j] << ")" << endl;
      fout << conts3[i].x[j] << " " << conts3[i].y[j] << " ";
      t.test_rel(fun2(conts3[i].x[j],conts3[i].y[j]),conts3[i].level,
		 1.7e-1,"curve");
    }
    fout << endl;
  }
  cout << endl;
  
  // ------------------------------------------------------------

  if (false) {
  
    cout << "Stress test:" << endl;

    rng<> ran;
  
    // Test non-square data
    ubvector tsx(3), tsy(3), tslev(1);
    ubmatrix tsd(3,3);
    
    for(i=0;i<3;i++) tsx[i]=((double)i);
    for(i=0;i<3;i++) tsy[i]=((double)i);
    tslev[0]=0.5;

    for(size_t count=0;count<100;count++) {
      cout << "count: " << count << endl;

      for(i=0;i<3;i++) {
	for(j=0;j<3;j++) {
	  tsd(j,i)=ran.random();
	}
      }
    
      co.set_data(3,3,tsx,tsy,tsd);
      co.set_levels(1,tslev);
      vector<contour_line> conts4;
      if (count==50) {
	co.debug_next_point=true;
      }
      co.calc_contours(conts4);
      nc=conts4.size();
      size_t cs=conts4[0].x.size();
    
      if (nc>0 && (conts4[0].x[0]!=conts4[0].x[cs-1] ||
		   conts4[0].y[0]!=conts4[0].y[cs-1])) {
	bool xedge1=false, yedge1=false;
	bool xedge2=false, yedge2=false;
	if (conts4[0].x[0]==0.0 || conts4[0].x[0]==2.0) {
	  xedge1=true;
	}
	if (conts4[0].x[cs-1]==0.0 || conts4[0].x[cs-1]==2.0) {
	  xedge2=true;
	}
	if (conts4[0].y[0]==0.0 || conts4[0].y[0]==2.0) {
	  yedge1=true;
	}
	if (conts4[0].y[cs-1]==0.0 || conts4[0].y[cs-1]==2.0) {
	  yedge2=true;
	}
	cout << xedge1 << " " << yedge1 << " "
	     << xedge2 << " " << yedge2 << endl;
	if (!((xedge1 || yedge1) && (xedge2 || yedge2))) {

	  for(size_t i=0;i<cs;i++) {
	    cout << conts4[0].x[i] << " ";
	    cout << conts4[0].y[i] << endl;
	  }

	  for(size_t i=0;i<cs;i++) {
	    for(size_t j=0;j<cs;j++) {
	      cout << i << " " << j << " "
		   << conts4[0].x[i] << " " << conts4[0].y[i] << " "
		   << conts4[0].x[j] << " " << conts4[0].y[j] << " "
		   << sqrt(pow(conts4[0].x[i]-conts4[0].x[j],2.0)+
			   pow(conts4[0].y[i]-conts4[0].y[j],2.0)) << endl;
	    }
	  }
	
	  xed.clear();
	  yed.clear();
	  co.get_edges(xed,yed);
	
	  print_data_xhoriz(3,3,tsx,tsy,tsd);
	
	  cout << "Print edges: " << endl;
	  for(size_t ir=0;ir<yed.size();ir++) {
	    cout << "Level: " << tslev[ir] << endl;
	    co.print_edges_xhoriz(xed[ir],yed[ir]);
	  }
	  cout << endl;
	
	  O2SCL_ERR("Contour failed.",exc_esanity);
	}
      }
    
    }

  }
  
  // ------------------------------------------------------------
  
  fout.close();
  
  t.report();
  return 0;
}
