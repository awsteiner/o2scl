/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011, Andrew W. Steiner
  
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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/planar_intp.h>
#include <o2scl/contour.h>
#include <o2scl/gsl_rnga.h>
#include <o2scl/graph.h>

using namespace std;
using namespace o2scl;

double fun(double x, double y) {
  return exp(-x*x-y*y);
}

int main(void) {
  gsl_rnga gr;
  ovector ox, oy, *od=new ovector[1];

  for(size_t i=0;i<10;i++) {
    double x=gr.random()*2.0-1.0;
    double y=gr.random()*2.0-1.0;
    double f=fun(x,y);
    ox.push_back(x);
    oy.push_back(y);
    od[0].push_back(f);
  }

  planar_intp<ovector,ovector *> pi;
  pi.set_data(ox.size(),ox,oy,1,od);

  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","p1",-1.01,-1.01,1.01,1.01,
			 0,0,700,700);

  p1->SetRightMargin(0.07);
  p1->SetTopMargin(0.07);
  
  vector<string> colors;

  ovector res(1);
  double dd=0.003;
  for(double x=-1.0;x<=1.0+dd/10.0;x+=dd) {
    for(double y=-1.0;y<=1.0+dd/10.0;y+=dd) {
      size_t i1,i2,i3;
      double x1,x2,x3;
      double y1,y2,y3;
      pi.interp(x,y,res,i1,x1,y1,i2,x2,y2,i3,x3,y3);
      size_t tmp;
      if (i2<i1) {
	tmp=i1;
	i1=i2;
	i2=tmp;
      }
      if (i3<i1) {
	tmp=i1;
	i1=i3;
	i3=tmp;
      }
      if (i3<i2) {
	tmp=i2;
	i2=i3;
	i3=tmp;
      }
      string scol=itos(i1)+itos(i2)+itos(i3);
      bool found=false;
      for(size_t j=0;j<colors.size();j++) {
	if (colors[j]==scol) found=true;
      }
      if (found==false) {
	colors.push_back(scol);
      }
      size_t tcol=0;
      for(size_t j=0;j<colors.size();j++) {
	if (colors[j]==scol) {
	  tcol=j+12; 
	  j=colors.size()+1;
	}
      }

      
      TBox *b1=new TBox(x-dd/2.0,y-dd/2.0,x+dd/2.0,y+dd/2.0);
      //b1->SetFillStyle(4100);
      //TMarker *m1=new TMarker(x,y,o2scl_graph::m_fill_square);
      b1->SetFillColor(tcol);
      b1->Draw();
    }
  }

  for(size_t i=0;i<ox.size();i++) {
    double x=ox[i];
    double y=oy[i];
    TBox *b1=new TBox(x-dd*3.0,y-dd*3.0,x+dd*3.0,y+dd*3.0);
    b1->SetFillColor(1);
    b1->Draw();
  }


  c1->Print("ex_planar_plot.png");
  system("convert ex_planar_plot.eps ex_planar_plot.png");

  return 0;
}
