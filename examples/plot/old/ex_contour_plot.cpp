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
#include <o2scl/graph.h>
#include <o2scl/test_mgr.h>
#include <o2scl/contour.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

double fun(double x, double y) {
  return 15.0*exp(-pow(x-20.0,2.0)/400.0-pow(y-5.0,2.0)/25.0)
    +40.0*exp(-pow(x-70.0,2.0)/500.0-pow(y-2.0,2.0)/4.0);
}

int main(void) {

  cout.setf(ios::scientific);

  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  // ------------------------------------------------------------
  // Read contour information from file

  ovector x1, x2, y1, y2;
  omatrix data1, data2;
  vector<contour_line> conts1, conts2;
  vector<edge_crossings> recs1, recs2;
  vector<edge_crossings> becs1, becs2;

  hdf_file hf;
  hf.open("../ex_contour.o2");
  hf.getd_vec("c1_x",x1);
  hf.getd_vec("c2_x",x2);
  hf.getd_vec("c1_y",y1);
  hf.getd_vec("c2_y",y2);
  hf.getd_mat("c1_data",data1);
  hf.getd_mat("c2_data",data2);
  hdf_input(hf,conts1,"c1_cl");
  hdf_input(hf,recs1,"c1_re");
  hdf_input(hf,becs1,"c1_be");
  hdf_input(hf,conts2,"c2_cl");
  hdf_input(hf,recs2,"c2_re");
  hdf_input(hf,becs2,"c2_be");
  hf.close();

  // ------------------------------------------------------------
  
  double xlow=x1[0]-8.0;
  double ylow=y1[0]-1.0;
  double xhigh=x1[x1.size()-1]+8.0;
  double yhigh=y1[y1.size()-1]+1.0;

  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","p1",xlow,ylow,xhigh,yhigh,
			 0,0,650,650);
  p1->SetTopMargin(0.05);
  p1->SetBottomMargin(0.07);
  p1->SetLeftMargin(0.07);
  p1->SetRightMargin(0.05);
  
  for(size_t i=0;i<conts1.size();i++) {
    TGraph *g1=new TGraph(conts1[i].x.size());
    for(size_t j=0;j<conts1[i].x.size();j++) {
      g1->SetPoint(j,conts1[i].x[j],conts1[i].y[j]);
    }
    g1->Draw();
  }
  for(size_t i=0;i<x1.size();i++) {
    for(size_t j=0;j<y1.size();j++) {
      TMarker *m1=new TMarker(x1[i],y1[j],o2scl_graph::m_plus);
      m1->Draw();
    }
  }

  c1->Print("ex_contour_plot1.eps");
  system("convert ex_contour_plot1.eps ex_contour_plot1.png");

  // ------------------------------------------------------------

  o2scl_graph::new_graph(c2,p2,th2,"c2","cc2","p2",xlow,ylow,xhigh,yhigh,
			 0,0,650,650);
  p2->SetTopMargin(0.05);
  p2->SetBottomMargin(0.07);
  p2->SetLeftMargin(0.07);
  p2->SetRightMargin(0.05);
  
  for(size_t i=0;i<conts2.size();i++) {
    TGraph *g1=new TGraph(conts2[i].x.size());
    for(size_t j=0;j<conts2[i].x.size();j++) {
      g1->SetPoint(j,conts2[i].x[j],conts2[i].y[j]);
    }
    g1->Draw();
  }
  for(size_t i=0;i<x1.size();i++) {
    for(size_t j=0;j<y1.size();j++) {
      TMarker *m1=new TMarker(x1[i],y1[j],o2scl_graph::m_plus);
      m1->Draw();
    }
  }

  c2->Print("ex_contour_plot2.eps");
  system("convert ex_contour_plot2.eps ex_contour_plot2.png");

  // ------------------------------------------------------------

  o2scl_graph::new_graph(c3,p3,th3,"c3","cc3","p3",xlow,ylow,xhigh,yhigh,
			 0,0,650,650);
  p3->SetTopMargin(0.05);
  p3->SetBottomMargin(0.07);
  p3->SetLeftMargin(0.07);
  p3->SetRightMargin(0.05);

  for(size_t i=0;i<x1.size();i++) {
    for(size_t j=0;j<y1.size();j++) {
      TMarker *m1=new TMarker(x1[i],y1[j],o2scl_graph::m_plus);
      m1->Draw();
    }
  }

  for(size_t i=0;i<recs1[1].status.rows();i++) {
    for(size_t j=0;j<recs1[1].status.cols();j++) {
      if (recs1[1].status[i][j]==contour::endpoint) {
	TLine *a1=new TLine(x1[j],y1[i],
			    x1[j],recs1[1].values[i][j]);
	a1->SetLineColor(kRed);
	a1->Draw();
	TMarker *m1=new TMarker(x1[j],recs1[1].values[i][j],
				o2scl_graph::m_fill_square);
	m1->SetMarkerColor(kRed);
	m1->Draw();
      } else if (recs1[1].status[i][j]==contour::contourp) {
	TLine *a1=new TLine(x1[j],y1[i],
			      x1[j],recs1[1].values[i][j]);
	a1->SetLineColor(kBlue);
	a1->Draw();
	TMarker *m1=new TMarker(x1[j],recs1[1].values[i][j],
				o2scl_graph::m_fill_square);
	m1->SetMarkerColor(kBlue);
	m1->Draw();
      }
    }
  }

  for(size_t i=0;i<becs1[1].status.rows();i++) {
    for(size_t j=0;j<becs1[1].status.cols();j++) {
      if (becs1[1].status[i][j]==contour::endpoint) {
	TLine *a1=new TLine(x1[j],y1[i],
			    becs1[1].values[i][j],y1[i]);
	a1->SetLineColor(kRed);
	a1->Draw();
	TMarker *m1=new TMarker(becs1[1].values[i][j],y1[i],
				o2scl_graph::m_fill_circle);
	m1->SetMarkerColor(kRed);
	m1->Draw();
      } else if (becs1[1].status[i][j]==contour::contourp) {
	TLine *a1=new TLine(x1[j],y1[i],
			    becs1[1].values[i][j],y1[i]);
	a1->SetLineColor(kBlue);
	a1->Draw();
	TMarker *m1=new TMarker(becs1[1].values[i][j],y1[i],
				o2scl_graph::m_fill_circle);
	m1->SetMarkerColor(kBlue);
	m1->Draw();
      }
    }
  }

  theApp.Run(kTRUE);
  c3->Print("ex_contour_plot3.eps");
  system("convert ex_contour_plot3.eps ex_contour_plot3.png");

  return 0;
}
