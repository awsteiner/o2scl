/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
int print_data(int nx, int ny, vec_t &x, vec_t &y, mat_t &data) {
  int j, k;
  cout.setf(ios::showpos);
  cout.precision(3);
  cout << "              ";
  for(k=0;k<ny;k++) {
    string stx="y";
    stx+=itos(k);
    cout.width(11);
    cout << stx;
  }
  cout << endl;
  cout << "               ";
  for(k=0;k<ny;k++) {
    cout << y[k] << " ";
  }
  cout << endl;
  for(j=0;j<nx;j++) {
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
  print_data(5,5,xs,ys,datas);
  
  // Regrid
  cout << "In regrid: " << endl;
  co.regrid_data(2,2);
  cout << "Done with regrid." << endl;
  
  size_t ngx=0, ngy=0;
  ubvector *xg=0, *yg=0;
  ubmatrix *datag=0;
  co.get_data(ngx,ngy,xg,yg,datag);
  
  print_data(ngx,ngy,*xg,*yg,*datag);

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
  print_data(10,10,x,y,data);

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
      if (i==5 && j==5) {
	cout << i << " " << j << " " << sqd(j,i) << endl;
      }
    }
  }
  sqlev[0]=4.0;
  sqlev[1]=10.0;
  sqlev[2]=20.0;
  sqlev[3]=40.0;

  print_data(10,8,sqx,sqy,sqd);

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
  
  print_data(8,10,srx,sry,srd);

  co.set_data(8,10,srx,sry,srd);
  co.set_levels(4,srlev);
  vector<contour_line> conts3;
  co.calc_contours(conts3);
  nc=conts3.size();
  
  cout << "Print edges: " << endl;
  vector<edge_crossings> red, bed;
  co.get_edges(red,bed);
  for(size_t ir=0;ir<red.size();ir++) {
    co.print_edges(red[ir],bed[ir]);
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

  /*
  cout << "Contour with corners:" << endl;
  double cordat[7][7]={{1,2,3,4,5,6,7},
		       {0,1,2,3,4,5,6},
		       {1,0,1,2,3,4,5},
		       {2,1,0,1,2,3,4},
		       {3,2,1,0,1,2,3},
		       {4,3,2,1,0,1,2},
		       {5,4,3,2,1,0,1}};
  double corx[7]={0,1,2,3,4,5,6};
  double cory[7]={0,1,2,3,4,5,6};
  double corlev[4]={1,2,3,4};

  print_data(7,7,corx,cory,cordat);

  co.set_data(7,7,corx,cory,cordat);
  co.set_levels(4,corlev);
  vector<contour_line> conts4;
  co.calc_contours(conts4);
  nc=conts4.size();

  fout << nc << endl;
  fout.precision(10);
  for(i=0;i<nc;i++) {
    cout << "Contour " << i << " at level " 
	 << conts4[i].level << ":" << endl;
    int cs=conts4[i].x.size();
    fout << cs << " ";
    for(j=0;j<cs;j++) {
      fout << conts4[i].x[j] << " " << conts4[i].y[j] << " ";
      cout << "(" << conts4[i].x[j] << ", " 
	   << conts4[i].y[j] << ")" << endl;
    }
    fout << endl;
  }
  cout << endl;
  */

  fout.close();
  
  t.report();
  return 0;
}
