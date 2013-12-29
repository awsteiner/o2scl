/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "graph.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_graph;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

void o2scl_graph::axis_labels
(double left, double bottom, double right, double top,
 int talign, string xtitle, string ytitle, bool logx, bool logy,
 double xfactor, double yfactor) {

  TLatex tt;
  tt.SetTextFont(132);
  tt.SetTextAlign(talign);

  if (logx) {
    if (logy) {
      tt.SetTextAngle(0);
      if (xtitle.length()>0) {
	tt.DrawLatex(sqrt(left*right),
		     bottom/pow(top/bottom,1.0/xfactor),xtitle.c_str());
      }
      tt.SetTextAngle(90);
      if (ytitle.length()>0) {
	tt.DrawLatex(left/pow(right/left,1.0/yfactor),
		     sqrt(top*bottom),ytitle.c_str());
      }
    } else {
      tt.SetTextAngle(0);
      if (xtitle.length()>0) {
	tt.DrawLatex(sqrt(left*right),
		     bottom-(top-bottom)/xfactor,xtitle.c_str());
      }
      tt.SetTextAngle(90);
      if (ytitle.length()>0) {
	tt.DrawLatex(left/pow(right/left,1.0/yfactor),
		     (top+bottom)/2.0,ytitle.c_str());
      }
    }
  } else {
    if (logy) {
      tt.SetTextAngle(0);
      if (xtitle.length()>0) {
	tt.DrawLatex((left+right)/2.0,
		     bottom/pow(top/bottom,1.0/xfactor),xtitle.c_str());
      }
      tt.SetTextAngle(90);
      if (ytitle.length()>0) {
	tt.DrawLatex(left-(right-left)/yfactor,
		     sqrt(top*bottom),ytitle.c_str());
      }
    } else {
      tt.SetTextAngle(0);
      if (xtitle.length()>0) {
	tt.DrawLatex((left+right)/2.0,
		     bottom-(top-bottom)/xfactor,xtitle.c_str());
      }
      tt.SetTextAngle(90);
      if (ytitle.length()>0) {
	tt.DrawLatex(left-(right-left)/yfactor,
		     (top+bottom)/2.0,ytitle.c_str());
      }
    }
  }
  tt.SetTextAngle(0);

  return;
}

void o2scl_graph::arrow(double x1, double y1, double x2, double y2,
			TLine *&line, TPolyLine *&poly, 
			double size, double size2, double alpha1) {
  
  line=new TLine(x1,y1,x2,y2);
  line->Draw();

  double theta1=atan2(y1-y2,x1-x2);
  double l1=size*sqrt(pow(x1-x2,2.0)+pow(y1-y2,2.0));
  double l2=l1/cos(alpha1);

  double x3=l2*cos(theta1+alpha1)+x2;
  double y3=l2*sin(theta1+alpha1)+y2;
  //  m1=new TMarker(x3,y3,o2scl_graph::m_fill_circle);
  //  m1->Draw();

  double x4=(x1-x2)*size2*size+x2;
  double y4=(y1-y2)*size2*size+y2;
  //  m1=new TMarker(x4,y4,o2scl_graph::m_fill_circle);
  //  m1->Draw();

  double x5=l2*cos(theta1-alpha1)+x2;
  double y5=l2*sin(theta1-alpha1)+y2;
  //  m1=new TMarker(x5,y5,o2scl_graph::m_fill_circle);
  //  m1->Draw();

  double xx[5]={x2,x3,x4,x5,x2};
  double yy[5]={y2,y3,y4,y5,y2};
  
  poly=new TPolyLine(5,xx,yy);
  poly->SetFillStyle(1001);
  poly->SetFillColor(1);
  poly->SetLineColor(1);
  poly->Draw("f");
  
  return;
}

TGraph *o2scl_graph::table_graph
(const table_units<> &at, string scolx, string scoly, int style, int color) {
  TGraph *g1=new TGraph(at.get_nlines());
  for(size_t i=0;i<at.get_nlines();i++) {
    g1->SetPoint(i,at.get(scolx,i),at.get(scoly,i));
  }
  g1->SetLineStyle(style);
  g1->SetLineColor(color);
  g1->SetName(scoly.c_str());
  g1->Draw();
  return g1;
}

TGraphErrors *o2scl_graph::table_graph_errors
(const table_units<> &at, string scolx, string scoly, string xerr, 
 string yerr, int style, int color) {
  TGraphErrors *g1=new TGraphErrors(at.get_nlines());
  for(size_t i=0;i<at.get_nlines();i++) {
    g1->SetPoint(i,at.get(scolx,i),at.get(scoly,i));
    if (xerr.length()>0) {
      if (yerr.length()>0) {
	g1->SetPointError(i,at.get(xerr,i),at.get(yerr,i));
      } else {
	g1->SetPointError(i,at.get(xerr,i),0.0);
      }
    } else {
      if (yerr.length()>0) {
	g1->SetPointError(i,0.0,at.get(yerr,i));
      } else {
	g1->SetPointError(i,0.0,0.0);
      }
    }
  }
  g1->SetLineStyle(style);
  g1->SetLineColor(color);
  g1->SetName(scoly.c_str());
  g1->Draw();
  return g1;
}

void o2scl_graph::new_graph
(TCanvas *&c1, TPad *&p1, TH1 *&th1, std::string canvas_name, 
 std::string window_name, std::string pad_name, 
 double lleft, double lbottom, double lright, double ltop, 
 int left, int top, int right, int bottom, bool logx, bool logy) {

  c1=new TCanvas(canvas_name.c_str(),window_name.c_str(),
		 left,top,right,bottom);
  c1->SetFillColor(10);
  p1=new TPad(pad_name.c_str(),"",0.0,0.0,1.0,1.0);
  p1->SetTopMargin(0.05);
  p1->SetRightMargin(0.05);
  p1->SetLeftMargin(0.11);
  p1->SetBottomMargin(0.1);
  p1->SetFillColor(10);
  p1->Draw();
  p1->cd();
  if (logx) p1->SetLogx();
  if (logy) p1->SetLogy();
  th1=p1->DrawFrame(lleft,lbottom,lright,ltop);
  th1->GetXaxis()->SetLabelFont(132);
  th1->GetYaxis()->SetLabelFont(132);
  th1->GetXaxis()->CenterTitle(kTRUE);
  th1->GetYaxis()->CenterTitle(kTRUE);
  th1->GetXaxis()->SetTitleFont(132);
  th1->GetYaxis()->SetTitleFont(132);
  // We explicitly set ndivisions here, and we want to make
  // sure that it's consistent with new_graph_ticks().
  th1->GetXaxis()->SetNdivisions(510);
  th1->GetYaxis()->SetNdivisions(510);

  return;
}

void o2scl_graph::new_graph_ticks
(TCanvas *&c1, TPad *&p1, TH1 *&th1, std::string canvas_name, 
 std::string window_name, std::string pad_name, 
 double lleft, double lbottom, double lright, double ltop, 
 TGaxis *&ax, TGaxis *&ay, int left, int top, int right, int bottom, 
 bool logx, bool logy) {

  new_graph(c1,p1,th1,canvas_name,window_name,pad_name,lleft,lbottom,
	    lright,ltop,left,top,right,bottom,logx,logy);

  if (logx) {
    ax=new TGaxis(lleft,ltop,lright,ltop,lleft,lright,510,"-G");
  } else {
    ax=new TGaxis(lleft,ltop,lright,ltop,lleft,lright,510,"-");
  }
  ax->SetLabelFont(132);
  ax->SetLabelSize(0.0);
  ax->CenterTitle(kTRUE);
  ax->Draw();

  if (logy) {
    ay=new TGaxis(lright,lbottom,lright,ltop,lbottom,ltop,510,"+G");
  } else {
    ay=new TGaxis(lright,lbottom,lright,ltop,lbottom,ltop,510,"+");
  }
  ay->SetLabelFont(132);
  ay->SetLabelSize(0.0);
  ay->CenterTitle(kTRUE);
  ay->Draw();

  return;
}

void o2scl_graph::two_up_graphy
(TCanvas *&c1, TPad *&p1, TPad *&p2, TH1 *&th1, TH1 *&th2,
 std::string canvas_name, std::string window_name, 
 std::string Pad1_name, std::string Pad2_name,
 double lowx, double highx, double lowy1, double highy1, 
 double lowy2, double highy2, int left, int top, int right, int bottom, 
 bool logx, bool logy1, bool logy2, double alpha, double margin) {

  c1=new TCanvas(canvas_name.c_str(),window_name.c_str(),0,0,700,700);
  c1->SetFillColor(10);
  
  double y1top, y2bottom;
  
  // These calculate the proper value of the pad limits 
  // assuming that both pads have margins of width 'margin' except
  // for the bottom margin of pad 'p1' which has a margin of
  // width alpha*margin. This ensures that the graphs are connected
  // and that they are the same width. These relations are the
  // solution of the two equations:
  //
  // size of graph 1 = size of graph 2:
  // (1-m1l-m1r)(x1r-x1l) = (1-m2l-m2r)(x2r-x2l)
  //
  // right edge of graph 1 = left edge of graph 2:
  // x1l+(x1r-x1l)*(1-m1r) = x2l+m2l*(x2r-x2l)
  //
  // In this particular case, we have taken x1l=0 and x2r=1.
  // Allowing alpha > 1 ensures that there is enough space on
  // the left hand side of graph 1 for the axis label.
  
  y1top=(1-2*margin)/(margin-1)/(margin*(3+alpha)-2);
  y2bottom=(1-4*margin+margin*margin*(3+alpha))/(margin-1)/
    (margin*(3+alpha)-2);
  
  p1=new TPad(Pad1_name.c_str(),"",0.0,0.0,1.0,y1top);
  p1->SetFillColor(10);
  p1->SetLeftMargin(margin*alpha);
  p1->SetRightMargin(margin);
  p1->SetBottomMargin(margin*alpha);
  p1->SetTopMargin(margin);
  p1->SetFillStyle(4000);
  p1->Draw();

  if (logx) p1->SetLogx();
  if (logy1) p1->SetLogy();
  
  p1->cd();
  th1=p1->DrawFrame(lowx,lowy1,highx,highy1);
  th1->GetXaxis()->SetLabelFont(132);
  th1->GetYaxis()->SetLabelFont(132);
  th1->GetXaxis()->SetLabelSize(th1->GetXaxis()->GetLabelSize()*2.0);
  th1->GetYaxis()->SetLabelSize(th1->GetYaxis()->GetLabelSize()*2.0);
  
  c1->cd();
  p2=new TPad(Pad2_name.c_str(),"",0.0,y2bottom,1.0,1.0);
  p2->SetFillColor(10);
  // Set the fill style to 4000 so that the pads are transparent
  p2->SetFillStyle(4000);
  p2->SetLeftMargin(margin*alpha);
  p2->SetRightMargin(margin);
  p2->SetBottomMargin(margin);
  p2->SetTopMargin(margin);
  p2->Draw();
  
  if (logx) p2->SetLogx();
  if (logy2) p2->SetLogy();

  p2->cd();
  th2=p2->DrawFrame(lowx,lowy2,highx,highy2);
  th2->GetYaxis()->SetLabelFont(132);
  th2->GetXaxis()->SetLabelFont(132);
  th2->GetXaxis()->SetLabelColor(10);
  th2->GetXaxis()->SetLabelSize(0);
  th2->GetYaxis()->SetLabelSize(th2->GetYaxis()->GetLabelSize()*2.0);
  th2->GetXaxis()->SetTicks("+-");

  return;
}

void o2scl_graph::two_up_graph
(TCanvas *&c1, TPad *&p1, TPad *&p2, TH1 *&th1, TH1 *&th2,
 string canvas_name, string window_name, string Pad1_name, string Pad2_name,
 double lowx1, double highx1, double lowx2, double highx2,
 double lowy, double highy, int left, int top, int right, int bottom,
 bool logx1, bool logx2, bool logy, double alpha, double margin) {
  
  c1=new TCanvas(canvas_name.c_str(),window_name.c_str(),0,0,1000,700);
  c1->SetFillColor(10);
  
  double x1right, x2left;
  
  // These calculate the proper value of the pad limits 
  // assuming that both pads have margins of width 'margin' except
  // for the left margin of pad 'p1' which has a margin of
  // width alpha*margin. This ensures that the graphs are connected
  // and that they are the same width. These relations are the
  // solution of the two equations:
  //
  // size of graph 1 = size of graph 2:
  // (1-m1l-m1r)(x1r-x1l) = (1-m2l-m2r)(x2r-x2l)
  //
  // right edge of graph 1 = left edge of graph 2:
  // x1l+(x1r-x1l)*(1-m1r) = x2l+m2l*(x2r-x2l)
  //
  // In this particular case, we have taken x1l=0 and x2r=1.
  // Allowing alpha > 1 ensures that there is enough space on
  // the left hand side of graph 1 for the axis label.
  x1right=(1-2*margin)/(margin-1)/(margin*(3+alpha)-2);
  x2left=(1-4*margin+margin*margin*(3+alpha))/(margin-1)/
    (margin*(3+alpha)-2);

  p1=new TPad(Pad1_name.c_str(),"",0.0,0.04,x1right,1.0);
  p1->SetFillColor(10);
  p1->SetLeftMargin(margin*alpha);
  p1->SetRightMargin(margin);
  p1->SetBottomMargin(margin*alpha);
  p1->SetTopMargin(margin);
  p1->SetFillStyle(4000);
  p1->Draw();

  if (logx1) p1->SetLogx();
  if (logy) p1->SetLogy();
  
  p1->cd();
  th1=p1->DrawFrame(lowx1,lowy,highx1,highy);
  th1->GetXaxis()->SetLabelFont(132);
  th1->GetYaxis()->SetLabelFont(132);
  
  c1->cd();
  p2=new TPad(Pad2_name.c_str(),"",x2left,0.04,1.0,1.0);
  p2->SetFillColor(10);
  // Set the fill style to 4000 so that the pads are transparent
  p2->SetFillStyle(4000);
  p2->SetLeftMargin(margin);
  p2->SetRightMargin(margin);
  p2->SetBottomMargin(margin*alpha);
  p2->SetTopMargin(margin);
  p2->Draw();
  
  if (logx2) p2->SetLogx();
  if (logy) p2->SetLogy();

  p2->cd();
  th2=p2->DrawFrame(lowx2,lowy,highx2,highy);
  th2->GetXaxis()->SetLabelFont(132);
  th2->GetYaxis()->SetLabelFont(132);
  th2->GetYaxis()->SetLabelColor(10);
  th2->GetYaxis()->SetLabelSize(0);
  th2->GetYaxis()->SetTicks("+-");

  return;
}

root_color_manager::root_color_manager() {
  n_colors=0;
  min_allocate=500;
  n_allocated=0;
}
  
void root_color_manager::allocate(size_t na) {

  if (n_allocated<na) {

    if (na<min_allocate) n_allocated=min_allocate;
    else n_allocated=na;
      
    // Hack to force ROOT to allocate a larger palette
    Double_t r[]={0.0,1.0};
    Double_t g[]={0.0,1.0};
    Double_t b[]={0.0,1.0};
    Double_t stop[]={0.0,1.0};
    col_start=TColor::CreateGradientColorTable
      (2,stop,r,g,b,n_allocated);
  }
    
  return;
}

void root_color_manager::set_min_allocate(size_t n) {
  min_allocate=n;
  return;
}

void root_color_manager::colors_rgb
(size_t n, std::string red_func, std::string green_func, 
 std::string blue_func) {

  allocate(n);

  std::vector<double> color_r, color_g, color_b;
    
  FunctionParser fp[3];
  fp[0].Parse(red_func,"x");
  fp[1].Parse(green_func,"x");
  fp[2].Parse(blue_func,"x");
  for(double x=0.0;x<1.01;x+=1.0/((double)(n-1))) {
    double rt=fp[0].Eval(&x);
    if (rt<0.0) rt=0.0;
    if (rt>1.0) rt=1.0;
    double gt=fp[1].Eval(&x);
    if (gt<0.0) gt=0.0;
    if (gt>1.0) gt=1.0;
    double bt=fp[2].Eval(&x);
    if (bt<0.0) bt=0.0;
    if (bt>1.0) bt=1.0;
    color_r.push_back(rt);
    color_g.push_back(gt);
    color_b.push_back(bt);
  }
    
  // Redefine that palette with the user-specified colors
  for (size_t i=0;i<n;i++) {
    TColor *c=gROOT->GetColor(col_start+i);
    c->SetRGB(color_r[i],color_g[i],color_b[i]);
  }

  n_colors=n;
    
  return;
}

void root_color_manager::colors_hsv
(size_t n, std::string hue_func, std::string sat_func, 
 std::string val_func) {

  allocate(n);

  std::vector<double> color_r, color_g, color_b;
    
  FunctionParser fp[3];
  fp[0].Parse(hue_func,"x");
  fp[1].Parse(sat_func,"x");
  fp[2].Parse(val_func,"x");
  for(double x=0.0;x<1.01;x+=1.0/((double)(n-1))) {
    double h=fp[0].Eval(&x);
    double s=fp[1].Eval(&x);
    double v=fp[2].Eval(&x);
    double r, g, b;
    HSVtoRGB(h,s,v,r,g,b);
    color_r.push_back(r);
    color_g.push_back(g);
    color_b.push_back(b);
  }
    
  // Redefine that palette with the user-specified colors
  for (size_t i=0;i<n;i++) {
    TColor *c=gROOT->GetColor(col_start+i);
    c->SetRGB(color_r[i],color_g[i],color_b[i]);
  }

  n_colors=n;
    
  return;
}

void root_color_manager::colors_rainbow(size_t n) {
  
  allocate(n);

  std::vector<double> color_r, color_g, color_b;

  double h, s=1.0, v=1.0, r, g, b;
  for(size_t j=0;j<n;j++) {
    h=300.0*((double)(j))/((double)(n-1));
    HSVtoRGB(h,s,v,r,g,b);
    color_r.push_back(r);
    color_g.push_back(g);
    color_b.push_back(b);
  }

  // Redefine that palette with the user-specified colors
  for (size_t i=0;i<n;i++) {
    TColor *c=gROOT->GetColor(col_start+i);
    c->SetRGB(color_r[i],color_g[i],color_b[i]);
  }
    
  n_colors=n;
    
  return;
}
  
html_colors::html_colors() {

  struct color_s colors_array[147]={
    {"aliceblue",0,0xF0,0xF8,0xFF},{"antiquewhite",0,0xFA,0xEB,0xD7},
    {"aqua",0,0x00,0xFF,0xFF},{"aquamarine",0,0x7F,0xFF,0xD4},
    {"azure",0,0xF0,0xFF,0xFF},{"beige",0,0xF5,0xF5,0xDC},
    {"bisque",0,0xFF,0xE4,0xC4},{"black",0,0x00,0x00,0x00},
    {"blanchedalmond",0,0xFF,0xEB,0xCD},{"blue",0,0x00,0x00,0xFF},
    {"blueviolet",0,0x8A,0x2B,0xE2},{"brown",0,0xA5,0x2A,0x2A},
    {"burlywood",0,0xDE,0xB8,0x87},{"cadetblue",0,0x5F,0x9E,0xA0},
    {"chartreuse",0,0x7F,0xFF,0x00},{"chocolate",0,0xD2,0x69,0x1E},
    {"coral",0,0xFF,0x7F,0x50},{"cornflowerblue",0,0x64,0x95,0xED},
    {"cornsilk",0,0xFF,0xF8,0xDC},{"crimson",0,0xDC,0x14,0x3C},
    {"cyan",0,0x00,0xFF,0xFF},{"darkblue",0,0x00,0x00,0x8B},
    {"darkcyan",0,0x00,0x8B,0x8B},{"darkgoldenrod",0,0xB8,0x86,0x0B},
    {"darkgray",0,0xA9,0xA9,0xA9},{"darkgrey",0,0xA9,0xA9,0xA9},
    {"darkgreen",0,0x00,0x64,0x00},{"darkkhaki",0,0xBD,0xB7,0x6B},
    {"darkmagenta",0,0x8B,0x00,0x8B},{"darkolivegreen",0,0x55,0x6B,0x2F},
    {"darkorange",0,0xFF,0x8C,0x00},{"darkorchid",0,0x99,0x32,0xCC},
    {"darkred",0,0x8B,0x00,0x00},{"darksalmon",0,0xE9,0x96,0x7A},
    {"darkseagreen",0,0x8F,0xBC,0x8F},{"darkslateblue",0,0x48,0x3D,0x8B},
    {"darkslategray",0,0x2F,0x4F,0x4F},{"darkslategrey",0,0x2F,0x4F,0x4F},
    {"darkturquoise",0,0x00,0xCE,0xD1},{"darkviolet",0,0x94,0x00,0xD3},
    {"deeppink",0,0xFF,0x14,0x93},{"deepskyblue",0,0x00,0xBF,0xFF},
    {"dimgray",0,0x69,0x69,0x69},{"dimgrey",0,0x69,0x69,0x69},
    {"dodgerblue",0,0x1E,0x90,0xFF},{"firebrick",0,0xB2,0x22,0x22},
    {"floralwhite",0,0xFF,0xFA,0xF0},{"forestgreen",0,0x22,0x8B,0x22},
    {"fuchsia",0,0xFF,0x00,0xFF},{"gainsboro",0,0xDC,0xDC,0xDC},
    {"ghostwhite",0,0xF8,0xF8,0xFF},{"gold",0,0xFF,0xD7,0x00},
    {"goldenrod",0,0xDA,0xA5,0x20},{"gray",0,0x80,0x80,0x80},
    {"grey",0,0x80,0x80,0x80},{"green",0,0x00,0x80,0x00},
    {"greenyellow",0,0xAD,0xFF,0x2F},{"honeydew",0,0xF0,0xFF,0xF0},
    {"hotpink",0,0xFF,0x69,0xB4},{"indianred",0,0xCD,0x5C,0x5C},
    {"indigo",0,0x4B,0x00,0x82},{"ivory",0,0xFF,0xFF,0xF0},
    {"khaki",0,0xF0,0xE6,0x8C},{"lavender",0,0xE6,0xE6,0xFA},
    {"lavenderblush",0,0xFF,0xF0,0xF5},{"lawngreen",0,0x7C,0xFC,0x00},
    {"lemonchiffon",0,0xFF,0xFA,0xCD},{"lightblue",0,0xAD,0xD8,0xE6},
    {"lightcoral",0,0xF0,0x80,0x80},{"lightcyan",0,0xE0,0xFF,0xFF},
    {"lightgoldenrodyellow",0,0xFA,0xFA,0xD2},
    {"lightgray",0,0xD3,0xD3,0xD3},
    {"lightgrey",0,0xD3,0xD3,0xD3},{"lightgreen",0,0x90,0xEE,0x90},
    {"lightpink",0,0xFF,0xB6,0xC1},{"lightsalmon",0,0xFF,0xA0,0x7A},
    {"lightseagreen",0,0x20,0xB2,0xAA},{"lightskyblue",0,0x87,0xCE,0xFA},
    {"lightslategray",0,0x77,0x88,0x99},{"lightslategrey",0,0x77,0x88,0x99},
    {"lightsteelblue",0,0xB0,0xC4,0xDE},{"lightyellow",0,0xFF,0xFF,0xE0},
    {"lime",0,0x00,0xFF,0x00},{"limegreen",0,0x32,0xCD,0x32},
    {"linen",0,0xFA,0xF0,0xE6},{"magenta",0,0xFF,0x00,0xFF},
    {"maroon",0,0x80,0x00,0x00},{"mediumaquamarine",0,0x66,0xCD,0xAA},
    {"mediumblue",0,0x00,0x00,0xCD},{"mediumorchid",0,0xBA,0x55,0xD3},
    {"mediumpurple",0,0x93,0x70,0xDB},{"mediumseagreen",0,0x3C,0xB3,0x71},
    {"mediumslateblue",0,0x7B,0x68,0xEE},
    {"mediumspringgreen",0,0x00,0xFA,0x9A},
    {"mediumturquoise",0,0x48,0xD1,0xCC},
    {"mediumvioletred",0,0xC7,0x15,0x85},
    {"midnightblue",0,0x19,0x19,0x70},{"mintcream",0,0xF5,0xFF,0xFA},
    {"mistyrose",0,0xFF,0xE4,0xE1},{"moccasin",0,0xFF,0xE4,0xB5},
    {"navajowhite",0,0xFF,0xDE,0xAD},{"navy",0,0x00,0x00,0x80},
    {"oldlace",0,0xFD,0xF5,0xE6},{"olive",0,0x80,0x80,0x00},
    {"olivedrab",0,0x6B,0x8E,0x23},{"orange",0,0xFF,0xA5,0x00},
    {"orangered",0,0xFF,0x45,0x00},{"orchid",0,0xDA,0x70,0xD6},
    {"palegoldenrod",0,0xEE,0xE8,0xAA},{"palegreen",0,0x98,0xFB,0x98},
    {"paleturquoise",0,0xAF,0xEE,0xEE},
    {"palevioletred",0,0xDB,0x70,0x93},
    {"papayawhip",0,0xFF,0xEF,0xD5},{"peachpuff",0,0xFF,0xDA,0xB9},
    {"peru",0,0xCD,0x85,0x3F},{"pink",0,0xFF,0xC0,0xCB},
    {"plum",0,0xDD,0xA0,0xDD},{"powderblue",0,0xB0,0xE0,0xE6},
    {"purple",0,0x80,0x00,0x80},{"red",0,0xFF,0x00,0x00},
    {"rosybrown",0,0xBC,0x8F,0x8F},{"royalblue",0,0x41,0x69,0xE1},
    {"saddlebrown",0,0x8B,0x45,0x13},{"salmon",0,0xFA,0x80,0x72},
    {"sandybrown",0,0xF4,0xA4,0x60},{"seagreen",0,0x2E,0x8B,0x57},
    {"seashell",0,0xFF,0xF5,0xEE},{"sienna",0,0xA0,0x52,0x2D},
    {"silver",0,0xC0,0xC0,0xC0},{"skyblue",0,0x87,0xCE,0xEB},
    {"slateblue",0,0x6A,0x5A,0xCD},{"slategray",0,0x70,0x80,0x90},
    {"slategrey",0,0x70,0x80,0x90},{"snow",0,0xFF,0xFA,0xFA},
    {"springgreen",0,0x00,0xFF,0x7F},{"steelblue",0,0x46,0x82,0xB4},
    {"tan",0,0xD2,0xB4,0x8C},{"teal",0,0x00,0x80,0x80},
    {"thistle",0,0xD8,0xBF,0xD8},{"tomato",0,0xFF,0x63,0x47},
    {"turquoise",0,0x40,0xE0,0xD0},{"violet",0,0xEE,0x82,0xEE},
    {"wheat",0,0xF5,0xDE,0xB3},{"white",0,0xFF,0xFF,0xFF},
    {"whitesmoke",0,0xF5,0xF5,0xF5},{"yellow",0,0xFF,0xFF,0x00},
    {"yellowgreen",0,0x9A,0xCD,0x32}};

  for(size_t i=0;i<147;i++) {
    colors_array[i].index=i;
    cmap.insert(std::make_pair<>(colors_array[i].name,colors_array[i]));
  }

  col_start=0;

  return;
}

int html_colors::get_color_index(std::string s) const {
  std::map<std::string,struct color_s,
    o2scl::string_comp>::const_iterator citer=cmap.find(s);
  if (citer==cmap.end()) {
    O2SCL_ERR("Color not found in html_colors:;get_color_index().",
	      exc_einval);
  }
  return citer->second.index;
}

void html_colors::get_color_rgb(std::string s, double &r, double &g,
				double &b) const {

  std::map<std::string,struct color_s,
	   o2scl::string_comp>::const_iterator citer=cmap.find(s);
  if (citer==cmap.end()) {
    O2SCL_ERR("Color not found in html_colors:;get_color_index().",
	      exc_einval);
  }
  r=((double)citer->second.r)/255.0;
  g=((double)citer->second.g)/255.0;
  b=((double)citer->second.b)/255.0;
  cout << "Found color '" << s << "' (" << r << "," << g << ","
       << b << ")." << endl;
  return;
}

int html_colors::operator[](std::string s) const {
  std::map<std::string,struct color_s,
    o2scl::string_comp>::const_iterator citer=cmap.find(s);
  return citer->second.index+col_start;
}

void html_colors::add_colors(root_color_manager &rcm) {

  std::vector<double> red(147), green(147), blue(147);

  std::map<std::string,struct color_s,
    o2scl::string_comp>::const_iterator citer;
      
  for(citer=cmap.begin();citer!=cmap.end();citer++) {
    size_t i=citer->second.index;
    red[i]=((double)citer->second.r)/255.0;
    green[i]=((double)citer->second.g)/255.0;
    red[i]=((double)citer->second.b)/255.0;
  }
      
  rcm.set_colors(147,red,green,blue);
  col_start=rcm.get_col_start();
      
  return;
}

table3d_density_plot::table3d_density_plot() {
  logx=false;
  logy=false;
  logz=false;
  wleft=0;
  wtop=0;
  wright=680;
  wbottom=640;
  canvas_name="Canvas 1";
  window_name="Window 1";
  pad_name="Pad 1";
  xset=false;
  xleft=0.0;
  xright=0.0;
  yset=false;
  ybottom=0.0;
  ytop=0.0;
  zset=false;
  zbottom=0.0;
  ztop=0.0;
  xtitle="";
  ytitle="";
  prmar=0.19;
}

void table3d_density_plot::plot_canvas(TCanvas *&c1, TPad *&p1, TH1 *&th1,
				       o2scl::table3d &t) {
  
  // -----------------------------------------------------------------
  // Create the canvas

  c1=new TCanvas(canvas_name.c_str(),window_name.c_str(),
		 wleft,wtop,wright,wbottom);
  c1->SetFillColor(10);
  p1=new TPad(pad_name.c_str(),"",0.0,0.0,1.0,1.0);
  p1->SetFillColor(10);
  p1->SetTopMargin(0.05);
  p1->SetRightMargin(prmar);
  p1->SetLeftMargin(0.13);
  p1->SetBottomMargin(0.12);
  if (logx) p1->SetLogx();
  if (logy) p1->SetLogy();
  p1->Draw();
  p1->cd();

  // -----------------------------------------------------------------
  // Set limits

  size_t nx, ny;
  t.get_size(nx,ny);

  double left, right;
  if (xset) {
    left=xleft;
    right=xright;
  } else {
    left=t.get_grid_x(0);
    right=t.get_grid_x(nx-1);
  }
  double bottom, top;
  if (yset) {
    bottom=ybottom;
    top=ytop;
  } else {
    bottom=t.get_grid_y(0);
    top=t.get_grid_y(ny-1);
  }

  // -----------------------------------------------------------------
  // Draw axes
  
  th1=p1->DrawFrame(left,bottom,right,top);
  th1->GetXaxis()->SetLabelFont(132);
  th1->GetYaxis()->SetLabelFont(132);
  th1->GetXaxis()->CenterTitle(kTRUE);
  th1->GetYaxis()->CenterTitle(kTRUE);
  th1->GetXaxis()->SetTitleFont(132);
  th1->GetYaxis()->SetTitleFont(132);

  o2scl_graph::axis_labels(left,bottom,right,top,22,
			   xtitle,ytitle,logx,logy);

  // -----------------------------------------------------------------
  // Update the canvas

  c1->Update();

  return;
}

void table3d_density_plot::plot(TPad *pad, o2scl::table3d &t, 
				std::string slice, root_color_manager &rcm) {

  size_t nx, ny;
  t.get_size(nx,ny);

  double pad_left=pad->GetUxmin();
  double pad_bottom=pad->GetUymin();
  double pad_right=pad->GetUxmax();
  double pad_top=pad->GetUymax();

  if (logy) {
    pad_top=pow(10.0,pad_top);
    pad_bottom=pow(10.0,pad_bottom);
  }

  // -----------------------------------------------------------------
  // Set colors if necessary

  if (rcm.get_n_colors()==0) {
    rcm.colors_rainbow(rcm.get_min_allocate());
  }
  int col_start=rcm.get_col_start();
  size_t ncol=rcm.get_n_colors();

  // -----------------------------------------------------------------
  // Compute z range

  double zmin, zmax;
  if (zset) {
    zmin=zbottom;
    zmax=ztop;
  } else {
    const ubmatrix &m=t.get_slice(slice);
    matrix_minmax(m.size1(),m.size2(),m,zmin,zmax);
  }

  // -----------------------------------------------------------------
  // Construct box boundaries

  ubvector xvleft(nx), xvright(nx);
  for(size_t i=0;i<nx;i++) {
    double dx;
    if (i<nx-1) dx=(t.get_grid_x(i+1)-t.get_grid_x(i))/2.0;
    else dx=(t.get_grid_x(i)-t.get_grid_x(i-1))/2.0;
    xvleft[i]=t.get_grid_x(i)-dx;
    xvright[i]=t.get_grid_x(i)+dx;
    if (i==0) {
      xvleft[i]=t.get_grid_x(0);
    }
    if (i==nx-1) {
      xvright[i]=t.get_grid_x(nx-1);
    }
  }

  ubvector yvleft(ny), yvright(ny);
  for(size_t i=0;i<ny;i++) {
    double dy;
    if (i<ny-1) dy=(t.get_grid_y(i+1)-t.get_grid_y(i))/2.0;
    else dy=(t.get_grid_y(i)-t.get_grid_y(i-1))/2.0;
    yvleft[i]=t.get_grid_y(i)-dy;
    yvright[i]=t.get_grid_y(i)+dy;
    if (i==0) {
      yvleft[i]=t.get_grid_y(0);
    }
    if (i==ny-1) {
      yvright[i]=t.get_grid_y(ny-1);
    }
  }

  // -----------------------------------------------------------------
  // Plot data

  for(size_t i=0;i<nx;i++) {
    for(size_t j=0;j<ny;j++) {
      
      // Compute color
      double cd;
      if (logz) {
	cd=(log10(t.get(i,j,slice))-log10(zmin))/
	  (log10(zmax)-log10(zmin))*(ncol-1);
      } else {
	cd=(t.get(i,j,slice)-zmin)/(zmax-zmin)*(ncol-1);
      }
      int color=((int)(cd+col_start));
      
      // Draw box
      if (xvleft[i]>=pad_left && xvright[i]<=pad_right && 
	  yvleft[j]>=pad_bottom && yvright[j]<=pad_top) {
	box1=new TBox(xvleft[i],yvleft[j],xvright[i],yvright[j]);
	box1->SetFillColor(color);
	box1->Draw();
      }
      
    }
  }

  // -----------------------------------------------------------------
  // Replot axes (we have to replot them because the boxes tend to
  // draw over the original axes given by DrawFrame().

  if (logy) {
    aright=new TGaxis(pad_right,pad_bottom,pad_right,pad_top,
		      pad_bottom,pad_top,510,"+G");
    aleft=new TGaxis(pad_left,pad_bottom,pad_left,pad_top,
		     pad_bottom,pad_top,510,"-G");
  } else {
    aright=new TGaxis(pad_right,pad_bottom,pad_right,pad_top,
		      pad_bottom,pad_top,510,"+");
    aleft=new TGaxis(pad_left,pad_bottom,pad_left,pad_top,
		     pad_bottom,pad_top,510,"-");
  }
  if (logx) {
    atop=new TGaxis(pad_left,pad_top,pad_right,pad_top,
		    pad_left,pad_right,510,"-G");
    abottom=new TGaxis(pad_left,pad_bottom,pad_right,pad_bottom,
		       pad_left,pad_right,510,"+G");
  } else {
    atop=new TGaxis(pad_left,pad_top,pad_right,pad_top,
		    pad_left,pad_right,510,"-");
    abottom=new TGaxis(pad_left,pad_bottom,pad_right,pad_bottom,
		       pad_left,pad_right,510,"+");
  }
  
  aleft->SetLabelSize(0.0);
  aright->SetLabelSize(0.0);
  atop->SetLabelSize(0.0);
  abottom->SetLabelSize(0.0);
  aleft->Draw();
  aright->Draw();
  atop->Draw();
  abottom->Draw();

  // -----------------------------------------------------------------
  // Plot z-axis legend

  // compute location
  double zleg_left=pad_right+(pad_right-pad_left)*0.02;
  double zleg_right=pad_right+(pad_right-pad_left)*0.09;

  for(size_t i=0;i<ncol;i++) {

    // Bottom and top margins of this box
    double bottom3, top3;
    if (logy) {
      bottom3=pad_bottom*pow(pad_top/pad_bottom,((double)i)/ncol);
      top3=pad_bottom*pow(pad_top/pad_bottom,((double)(i+1))/ncol);
    } else {
      bottom3=((double)i)/ncol*(pad_top-pad_bottom)+pad_bottom;
      top3=((double)(i+1))/ncol*(pad_top-pad_bottom)+pad_bottom;
    }
    
    // Compute color
    int color;
    if (logz) {
      double value=
	pow(10.0,((double)i+1)/ncol*(log10(zmax)-log10(zmin))+log10(zmin));
      double cd=(log10(value)-log10(zmin))/(log10(zmax)-log10(zmin))*ncol-1;
      color=((int)(cd+col_start));
    } else {
      color=((int)((i)*ncol/((double)ncol)+((double)col_start)));
    }

    box1=new TBox(zleg_left,bottom3,zleg_right,top3);
    box1->SetFillColor(color);
    box1->Draw();
  }

  if (logz) {
    ascale=new TGaxis(zleg_right,pad_bottom,zleg_right,
		      pad_top,zmin,zmax,510,"+LG");
  } else {
    ascale=new TGaxis(zleg_right,pad_bottom,zleg_right,
		      pad_top,zmin,zmax,510,"+L");
  }
  ascale->SetLabelFont(132);
  ascale->SetLabelSize(0.04);
  ascale->SetLabelOffset(0.01);
  ascale->SetTickSize(0.13);
  ascale->CenterTitle(kTRUE);
  ascale->Draw();
  
  l1=new TLine(zleg_left,pad_bottom,zleg_left,pad_top);
  l1->Draw();
  l2=new TLine(zleg_left,pad_bottom,zleg_right,pad_bottom);
  l2->Draw();
  l3=new TLine(zleg_left,pad_top,zleg_right,pad_top);
  l3->Draw();

  // -----------------------------------------------------------------
  // Plot contours if present
  
  for(size_t k=0;k<conts.size();k++) {
    gt.resize(conts.size());
    gt[k]=new TGraph(conts[k].x.size());
    for(size_t i=0;i<conts[k].x.size();i++) {
      gt[k]->SetPoint(i,conts[k].x[i],conts[k].y[i]);
    }
    gt[k]->Draw();
  }
  
  // -----------------------------------------------------------------

  return;
}

table3d_multi_density::table3d_multi_density() {
  n_bins=30;
  combine=add_colors;
}

void table3d_multi_density::multi_plot(root_color_manager &rcm) {

  if (tables.size()<1) {
    O2SCL_ERR2("No tables to plot in table3d_multi_density::",
	       "multi_plot().",exc_einval);
  }
  
  // -----------------------------------------------------------------
  // Construct grid for boxes later
  
  size_t nx, ny;
  tables[0].get_size(nx,ny);

  double left, right;
  if (xset) {
    left=xleft;
    right=xright;
  } else {
    left=tables[0].get_grid_x(0);
    right=tables[0].get_grid_x(nx-1);
  }
  double bottom, top;
  if (yset) {
    bottom=ybottom;
    top=ytop;
  } else {
    bottom=tables[0].get_grid_y(0);
    top=tables[0].get_grid_y(ny-1);
  }

  uniform_grid_end<double> xgrid(left,right,n_bins);
  uniform_grid_end<double> ygrid(bottom,top,n_bins);
  
  // -----------------------------------------------------------------
  // Text object

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);

  // -----------------------------------------------------------------
  // Compute slice maxima and minima

  ubvector max(tables.size()), min(tables.size());
  for(size_t i=0;i<tables.size();i++) {
    const ubmatrix &sl=tables[i].get_slice(slices[i]);
    min[i]=sl(0,0);
    max[i]=sl(0,0);
    for(size_t j=0;j<sl.size1();j++) {
      for(size_t k=0;k<sl.size2();k++) {
	if (sl(j,k)<min[i]) min[i]=sl(j,k);
	if (sl(j,k)>max[i]) max[i]=sl(j,k);
      }
    }
  }

  // -----------------------------------------------------------------
  // Allocate colors

  int n_colors;
  if (combine==add_colors) {
    n_colors=n_bins*n_bins+n_bins*tables.size();
  } else {
    n_colors=n_bins*tables.size();
  }
  rcm.allocate(n_colors);
  int col_index_start=rcm.get_col_start();
  
  double r_scale=0.0, g_scale=0.0, b_scale=0.0;
  double r_shift=0.0, g_shift=0.0, b_shift=0.0;
  
  double col_start[3]={1.0,1.0,1.0};
    
  if (combine==add_colors) {

    // -----------------------------------------------------------------
    // First pass to compute color scaling
    
    double r_max=0.0, g_max=0.0, b_max=0.0;
    double r_min=0.0, g_min=0.0, b_min=0.0;
    
    // Process colors in main plot area

    double progress=0.1;
    cout << "First pass: " << endl;
    for(size_t i=0;i<xgrid.get_nbins();i++) {
      for(size_t j=0;j<ygrid.get_nbins();j++) {

	double ratio=((double)(i*j))/
	  ((double)(xgrid.get_nbins()*ygrid.get_nbins()));
	if (ratio>progress) {
	  cout << 100.0*progress << " percent complete." << endl;
	  progress+=0.1;
	}
	
	double xc=(xgrid[i]+xgrid[i+1])/2.0;
	double yc=(ygrid[j]+ygrid[j+1])/2.0;
	
	double tr=col_start[0], tg=col_start[1], tb=col_start[2];
	for(size_t ik=0;ik<tables.size();ik++) {
	  double weight=(tables[ik].interp(xc,yc,slices[ik])-min[ik])/
	    (max[ik]-min[ik]);
	  tr+=weight*(colors[ik*3]-col_start[0]);
	  tg+=weight*(colors[ik*3+1]-col_start[1]);
	  tb+=weight*(colors[ik*3+2]-col_start[2]);
	}
	if (i==0 && j==0) {
	  r_max=tr;
	  g_max=tg;
	  b_max=tb;
	  r_min=tr;
	  g_min=tg;
	  b_min=tb;
	} else {
	  if (tr>r_max) r_max=tr;
	  if (tg>g_max) g_max=tg;
	  if (tb>b_max) b_max=tb;
	  if (tr<r_min) r_min=tr;
	  if (tg<g_min) g_min=tg;
	  if (tb<b_min) b_min=tb;
	}
      }
    }
    cout << 100.0 << " percent complete." << endl;
    cout << endl;

    // Process colors in z-axis legend
    
    for(size_t ik=0;ik<tables.size();ik++) {
      
      for(size_t i=0;i<n_bins;i++) {
	
	double val=((double)i)/((double)(n_bins-1));
	double tr=col_start[0]+val*(colors[ik*3]-col_start[0]);
	double tg=col_start[1]+val*(colors[ik*3+1]-col_start[1]);
	double tb=col_start[2]+val*(colors[ik*3+2]-col_start[2]);
	
	if (tr>r_max) r_max=tr;
	if (tg>g_max) g_max=tg;
	if (tb>b_max) b_max=tb;
	if (tr<r_min) r_min=tr;
	if (tg<g_min) g_min=tg;
	if (tb<b_min) b_min=tb;
      }
    }

    // Compute scales and shifts
    
    r_scale=fabs(r_max-r_min);
    if (r_max>r_min) r_shift=-r_min/r_scale;
    else r_shift=-r_max/r_scale;
    g_scale=fabs(g_max-g_min);
    if (g_max>g_min) g_shift=-g_min/g_scale;
    else g_shift=-g_max/g_scale;
    b_scale=fabs(b_max-b_min);
    if (b_max>b_min) b_shift=-b_min/b_scale;
    else b_shift=-b_max/b_scale;

    if (false) {
      cout << "r_max, r_min, r_scale, r_shift:\n\t" << r_max << " " 
	   << r_min << " " << r_scale << " " << r_shift << endl;
      cout << "g_max, g_min, g_scale, g_shift:\n\t" << g_max << " " 
	   << g_min << " " << g_scale << " " << g_shift << endl;
      cout << "b_max, b_min, b_scale, b_shift:\n\t" << b_max << " " 
	   << b_min << " " << b_scale << " " << b_shift << endl;
    }

  } else {
    
    // -----------------------------------------------------------------
    // Set colors
    
    for(size_t ik=0;ik<tables.size();ik++) {
      for(size_t i=0;i<n_bins;i++) {
	
	size_t col_index=ik*n_bins+i;
	
	// Compute and set the color
	TColor *c=gROOT->GetColor(col_index_start+col_index);
	double val=((double)i)/((double)(n_bins-1));
	double tr=col_start[0]+val*(colors[ik*3]-col_start[0]);
	double tg=col_start[1]+val*(colors[ik*3+1]-col_start[1]);
	double tb=col_start[2]+val*(colors[ik*3+2]-col_start[2]);
	c->SetRGB(tr,tg,tb);
	
      }
    }
    
  }
  
  // -----------------------------------------------------------------
  // Plot the boxes
  
  if (combine==add_colors) {

    double progress=0.1;
    cout << "Second pass: " << endl;
    for(size_t i=0;i<xgrid.get_nbins();i++) {
      for(size_t j=0;j<ygrid.get_nbins();j++) {
	
	double ratio=((double)(i*j))/
	  ((double)(xgrid.get_nbins()*ygrid.get_nbins()));
	if (ratio>progress) {
	  cout << 100.0*progress << " percent complete." << endl;
	  progress+=0.1;
	}

	double xc=(xgrid[i]+xgrid[i+1])/2.0;
	double yc=(ygrid[j]+ygrid[j+1])/2.0;
	
	// Compute the color index
	size_t col_index=i*xgrid.get_nbins()+j;
	
	// Compute the color
	TColor *c=gROOT->GetColor(col_index_start+col_index);
	double tr=col_start[0], tg=col_start[1], tb=col_start[2];

	for(size_t ik=0;ik<tables.size();ik++) {
	  double weight=(tables[ik].interp(xc,yc,slices[ik])-min[ik])/
	    (max[ik]-min[ik]);
	  tr+=weight*(colors[ik*3]-col_start[0]);
	  tg+=weight*(colors[ik*3+1]-col_start[1]);
	  tb+=weight*(colors[ik*3+2]-col_start[2]);
	}

	tr=tr/r_scale+r_shift;
	tg=tg/g_scale+g_shift;
	tb=tb/b_scale+b_shift;
	c->SetRGB(tr,tg,tb);

	box1=new TBox(xgrid[i],ygrid[j],xgrid[i+1],ygrid[j+1]);
	box1->SetFillColor(col_index_start+col_index);
	box1->Draw();
      }
    }
    cout << 100.0 << " percent complete." << endl;
    cout << endl;
    
  } else {
    
    for(size_t i=0;i<xgrid.get_nbins();i++) {
      for(size_t j=0;j<ygrid.get_nbins();j++) {
	
	double xc=(xgrid[i]+xgrid[i+1])/2.0;
	double yc=(ygrid[j]+ygrid[j+1])/2.0;
	
	for(size_t ik=0;ik<tables.size();ik++) {
	  
	  // Compute the color index
	  double val=tables[ik].interp(xc,yc,slices[ik]);
	  if (val>0.1) {
	    size_t col_index=ik*n_bins+((size_t)(val*((double)(n_bins-1))));
	    
	    box1=new TBox(xgrid[i],ygrid[j],xgrid[i+1],ygrid[j+1]);
	    box1->SetFillColor(col_index_start+col_index);
	    box1->SetFillStyle(patterns[ik]);
	    box1->Draw();
	  }
	}
      }
    }
    
  }

  // -----------------------------------------------------------------
  // Replot axes (we have to replot them because the boxes tend to
  // draw over the original axes given by DrawFrame().

  aright=new TGaxis(right,bottom,right,
		    top,bottom,top,510,"+");
  aleft=new TGaxis(left,bottom,left,top,
		   bottom,top,510,"-");
  atop=new TGaxis(left,top,right,top,left,
		  right,510,"-");
  abottom=new TGaxis(left,bottom,right,bottom,
		     left,right,510,"+");

  aleft->SetLabelSize(0.0);
  aright->SetLabelSize(0.0);
  atop->SetLabelSize(0.0);
  abottom->SetLabelSize(0.0);
  aleft->Draw();
  aright->Draw();
  atop->Draw();
  abottom->Draw();

  // -----------------------------------------------------------------
  // Plot z-axis legend
  
  // Horizontal location of the legend
  // FIXME: These need to be recomputed if logx is true!
  
  double scale=(right-left)*0.36/((double)tables.size());
  double margin=(right-left)*0.02;

  for(size_t ik=0;ik<tables.size();ik++) {

    double left3=right+margin+((double)ik)*scale;
    double right3=left3+scale;

    for(size_t i=0;i<n_bins;i++) {
          
      // Bottom and top margins of this box
      double bottom3, top3;
      if (logy) {
	bottom3=bottom*pow(top/bottom,((double)i)/n_bins);
	top3=bottom*pow(top/bottom,((double)(i+1))/n_bins);
      } else {
	bottom3=((double)i)/n_bins*(top-bottom)+bottom;
	top3=((double)(i+1))/n_bins*(top-bottom)+bottom;
      }
          
      // Compute color
      size_t col_index;

      if (combine==add_colors) {
      
	col_index=n_bins*n_bins+ik*n_bins+i;
	TColor *c=gROOT->GetColor(col_index_start+col_index);
      
	double val=((double)i)/((double)(n_bins-1));
	double tr=col_start[0]+val*(colors[ik*3]-col_start[0]);
	double tg=col_start[1]+val*(colors[ik*3+1]-col_start[1]);
	double tb=col_start[2]+val*(colors[ik*3+2]-col_start[2]);

	tr=tr/r_scale+r_shift;
	tg=tg/g_scale+g_shift;
	tb=tb/b_scale+b_shift;
	c->SetRGB(tr,tg,tb);

	// Draw the box
	box1=new TBox(left3,bottom3,right3,top3);
	box1->SetFillColor(col_index_start+col_index);
	box1->Draw();

      } else {

	// Compute color index
	col_index=ik*n_bins+i;
          
	// Draw the box
	box1=new TBox(left3,bottom3,right3,top3);
	box1->SetFillColor(col_index_start+col_index);
	box1->SetFillStyle(patterns[ik]);
	box1->Draw();
      }
          
    }

    tt.SetTextAngle(90);
    tt.DrawLatex((left3+right3)/2.0,(17.0*bottom+3.0*top)/20.0,
		 labels[ik].c_str());
    tt.SetTextAngle(0);
  }
  
  double left4=right+margin;
  double right4=left4+((double)tables.size())*scale;

  // Plot the RHS axis for the legend
  ascale=new TGaxis(right4,bottom,right4,
		    top,0.0,1.0,510,"+L");
  ascale->SetLabelFont(132);
  ascale->SetLabelSize(0.04);
  ascale->SetLabelOffset(0.01);
  ascale->SetTickSize(0.13);
  ascale->CenterTitle(kTRUE);
  ascale->Draw();

  l1=new TLine(left4,bottom,left4,top);
  l1->Draw();
  l2=new TLine(left4,bottom,right4,bottom);
  l2->Draw();
  l3=new TLine(left4,top,right4,top);
  l3->Draw();
      
  return;
}

