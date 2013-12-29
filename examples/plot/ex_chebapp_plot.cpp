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
#include <o2scl/table.h>
#include <o2scl/graph.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_ext;
using namespace o2scl_graph;

int main(void) {
  cout.setf(ios::scientific);

  // Fill a table with the data
  table tab;
  tab.line_of_names("x exact apx100 apx50 apx25");
  ifstream fin("../ex_chebapp.out");
  for(size_t i=0;i<101;i++) {
    double line[5];
    for(size_t j=0;j<5;j++) fin >> line[j];
    tab.line_of_data(5,line);
  }

  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2;
  TPad *p1, *p2;
  TH1 *th1, *th2;
  TGaxis *ax1, *ay1, *ax2, *ay2;

  new_graph_ticks(c1,p1,th1,"c1","d1","e1",0,-1.1,1.0,1.1,ax1,ay1,
		  0,0,700,700,false,false);
  th1->GetXaxis()->SetLabelSize(th1->GetXaxis()->GetLabelSize()*1.3);
  th1->GetYaxis()->SetLabelSize(th1->GetYaxis()->GetLabelSize()*1.3);
  p1->SetTopMargin(0.04);
  p1->SetBottomMargin(0.12);
  p1->SetLeftMargin(0.12);
  p1->SetRightMargin(0.04);
  
  table_fp *tfp=(table_fp *)&tab;
  
  TGraph *gr1=table_graph(tfp,"x","exact",1,kBlack);
  gr1->SetName("exact");
  gr1->Draw();
  TGraph *gr2=table_graph(tfp,"x","apx50",2,kRed);
  gr2->SetName("apx100");
  gr2->Draw();
  //TGraph *gr3=table_graph(tfp,"x","apx50",3,3);
  //gr3->Draw();
  TGraph *gr4=table_graph(tfp,"x","apx25",4,kBlue);
  gr4->SetName("apx25");
  gr4->Draw();

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);
  //tt.SetTextSize(tt.GetTextSize()*1.2);
  tt.DrawLatex(0.5,-1.27,"x");
  tt.SetTextAlign(32);
  tt.DrawLatex(0.77,-0.5,"Exact");
  tt.SetTextColor(kRed);
  tt.DrawLatex(0.77,-0.65,"Approx. (n=50)");
  tt.SetTextColor(kBlue);
  tt.DrawLatex(0.77,-0.80,"Approx. (n=25)");
  tt.SetTextColor(1);

  TLine *l1=new TLine(0.8,-0.5,0.95,-0.5);
  l1->Draw();
  TLine *l2=new TLine(0.8,-0.65,0.95,-0.65);
  l2->SetLineColor(kRed);
  l2->SetLineStyle(2);
  l2->Draw();
  TLine *l3=new TLine(0.8,-0.8,0.95,-0.8);
  l3->SetLineColor(kBlue);
  l3->SetLineStyle(4);
  l3->Draw();
  
  c1->Update();

  // Output to file or screen
  theApp.Run(kTRUE);

  c1->Print("ex_chebapp_plot.eps");
  
  system("convert ex_chebapp_plot.eps ex_chebapp_plot.png");

  return 0;
}
