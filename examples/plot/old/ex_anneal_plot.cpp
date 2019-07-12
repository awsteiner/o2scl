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
#include <gsl/gsl_sf_bessel.h>
#include <o2scl/graph.h>
#include <o2scl/test_mgr.h>
#include <o2scl/contour.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_graph;

double function(size_t nvar, const ovector_base &x) {
  double a, b;
  a=(x[0]-2.0);
  b=(x[1]+3.0);
  return -gsl_sf_bessel_J0(a)*gsl_sf_bessel_J0(b);
}

int main(void) {

  contour co;

  // Initialize the data
  size_t ngrid=100;
  ovector x(ngrid), y(ngrid), v(2);
  omatrix data(ngrid,ngrid);

  for(size_t i=0;i<ngrid;i++) {
    x[i]=((double)i)/ngrid*14.0-4.0;
    y[i]=((double)i)/ngrid*14.0-4.0;
  }
  
  for(size_t j=0;j<ngrid;j++) {
    for(size_t k=0;k<ngrid;k++) {
      v[0]=x[j];
      v[1]=y[k];
      data[k][j]=function(2,v);
    }
  }
  co.set_data(ngrid,ngrid,x,y,data);
  
  // Set the contour levels

  size_t nlev=9;
  ovector levels(nlev);
  levels[0]=-0.9;
  levels[1]=-0.7;
  levels[2]=-0.5;
  levels[3]=-0.3;
  levels[4]=-0.2;
  levels[5]=-0.1;
  levels[6]=-0.05;
  levels[7]=0.05;
  levels[8]=0.1;
  
  co.set_levels(nlev,levels);

  // Compute the contours

  vector<contour_line> conts;
  co.calc_contours(conts);

  // Initialize the plot

  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","p1",-4.0,-4.0,10.0,10.0,
                         0,0,700,700);
  p1->SetTopMargin(0.05);
  p1->SetBottomMargin(0.07);
  p1->SetLeftMargin(0.07);
  p1->SetRightMargin(0.05);

  // Plot colors

  double lcolr=0.2;
  double lcolg=0.2;
  double lcolb=1.0;
  
  double hcolr=1.0;
  double hcolg=0.2;
  double hcolb=0.2;

  for(size_t i=0;i<40;i++) {
    TColor *colort=(TColor *)(gROOT->GetListOfColors()->At(11+i));
    double cr=lcolr+((double)i)*(hcolr-lcolr)/39.0;
    double cg=lcolg+((double)i)*(hcolg-lcolg)/39.0;
    double cb=lcolb+((double)i)*(hcolb-lcolb)/39.0;
    colort->SetRGB(cr,cg,cb);
  }
  double bsize=0.04;
  ovector ov(2);
  for(double x=-4.0+bsize/2.0;x<10.0;x+=bsize) {
    for(double y=-4.0+bsize/2.0;y<10.0;y+=bsize) {
      TBox *b1=new TBox(x-bsize/2.0,y-bsize/2.0,
			x+bsize/2.0,y+bsize/2.0);
      ov[0]=x;
      ov[1]=y;
      b1->SetFillColor(function(2,ov)*20.0+31);
      b1->Draw();
    }
  }

  // Plot the contours

  size_t nc=conts.size();
  for(size_t i=0;i<nc;i++) {
    size_t cs=conts[i].x.size();
    TGraph *g1=new TGraph(cs);
    for(size_t j=0;j<cs;j++) {
      g1->SetPoint(j,conts[i].x[j],conts[i].y[j]);
    }
    g1->SetName(dtos(conts[i].level).c_str());
    g1->Draw();
  }

  // Hack to get the last value from ex_anneal.scr
  string s="grep \"x:\" ../ex_anneal.scr | tail -n 1 > temp";
  system(s.c_str());
  ifstream fin("temp");
  double x0, x1;
  fin >> s >> x0 >> x1;
  fin.close();

  // Plot the markers and text

  TMarker *m1=new TMarker(6,7,m_fill_square);
  m1->Draw();

  TMarker *m2=new TMarker(x0,x1,m_fill_circle);
  m2->Draw();

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);
  
  tt.DrawLatex(3.04663,-4.782197,"x");
  tt.DrawLatex(-4.74791,3.053977,"y");
  tt.DrawLatex(5.835293,6.368371,"Initial guess");
  tt.DrawLatex(1.995167,-2.391098,"Computed minimum");
  
  // 
  
  theApp.Run(kTRUE);
  
  c1->Print("ex_anneal_plot.eps");
  system("convert ex_anneal_plot.eps ex_anneal_plot.png");
  
  // The final .eps file was actually constructed by hand from the
  // .png file created by this code

  return 0;
}
