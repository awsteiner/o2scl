/*
  Plot the results of ex_nucmass.cpp
*/
#include <iostream>

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/graph.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);
  cout.precision(5);

  table_units tu;
  hdf_file hf;
  hf.open("../ex_nucmass_table.o2");
  hdf_input(hf,tu);
  hf.close();

  // The ROOT objects
  TApplication theApp("App",0,NULL);

  for(size_t ik=0;ik<6;ik++) {

    TCanvas *c1, *c2, *c3, *c4;
    TPad *p1, *p2, *p3, *p4;
    TH1 *th1, *th2, *th3, *th4;

    TLatex tt;
    tt.SetTextAlign(22);
    tt.SetTextFont(132);

    o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",0,0,120,180,
			   0,0,700,700,false,false);
    p1->SetLeftMargin(0.11);
    p1->SetRightMargin(0.19);
    p1->SetTopMargin(0.04);
    p1->SetBottomMargin(0.11);
    p1->cd();
    c1->Update();
  
    int min_col=11;
    int max_col=90;
    int n_col=max_col-min_col+1;
    for(int i=1;i<=n_col;i++) {
      TColor *c=(TColor *)(gROOT->GetListOfColors()->At(i+min_col-1));
      double cx=(((double)(i-1))/((double)(n_col-1)));
      c->SetRGB(1.0-cx,0.0,cx);
    }

    size_t magic[6]={8,20,28,50,82,126};
    for(size_t i=0;i<6;i++) {
      if (i<5) {
	TLine *l1=new TLine(magic[i],0,magic[i],180);
	l1->SetLineStyle(2);
	l1->Draw();
      }
      TLine *l2=new TLine(0,magic[i],120,magic[i]);
      l2->SetLineStyle(2);
      l2->Draw();
    }

    TBox *b1=new TBox(10,135,60,165);
    b1->SetFillStyle(1001);
    b1->SetFillColor(10);
    b1->Draw();

    if (false) {
      TBox *b1b=new TBox(10,140,60,160);
      b1b->SetFillStyle(0);
      b1b->SetLineColor(1);
      b1b->Draw();
    }
    
    string model="mn";
    if (ik==1) model="sm";
    if (ik==2) model="hfb";
    if (ik==3) model="ame13";
    if (ik==4) model="dz";
    if (ik==5) model="ktuy";

    // First make a pass to find the range
    double min_diff=1.0e100;
    double max_diff=0.0;
    for(size_t i=0;i<tu.get_nlines();i++) {
      if (fabs(tu.get(model,i))>1.0e-20) {
	double Z=tu.get("Z",i);
	double N=tu.get("N",i);
	double diff=tu.get(model,i)-tu.get("ame",i);
	if (diff>max_diff) max_diff=diff;
	if (diff<min_diff) min_diff=diff;
      }
    }
    //cout << min_diff << " " << max_diff << endl;

    for(int i=0;i<n_col;i++) {
      TBox *bx=new TBox(125,180.0/n_col*i,135,180.0/n_col*(i+1));
      bx->SetFillStyle(1001);
      bx->SetFillColor(i+min_col);
      bx->Draw();
    }

    TBox *b2=new TBox(125,0,135,180);
    b2->SetFillStyle(0);
    b2->SetLineColor(1);
    b2->Draw();

    TGaxis *ax=new TGaxis(135,0,135,180,min_diff,max_diff,510,"-");
    ax->SetLabelFont(132);
    ax->SetLabelOffset(-0.07);
    ax->CenterTitle(kTRUE);
    ax->Draw();

    // Total range on a log scale
    double range=max_diff-min_diff;

    for(size_t i=0;i<tu.get_nlines();i++) {
      if (fabs(tu.get(model,i))>1.0e-20 && tu.get("N",i)>1) {
	double Z=tu.get("Z",i);
	double N=tu.get("N",i);
	double diff=tu.get(model,i)-tu.get("ame",i);
	if (diff<min_diff) diff=min_diff;

	double current;
	current=diff-min_diff;
	int col=((int)(current/range*(max_col-min_col)+min_col));
      
	if (col<min_col || col>max_col) {
	  cout << "Color problem: " << endl;
	  cout << "min_col,max_col: " << min_col << " " << max_col << endl;
	  cout << "min,diff,max: " 
	       << min_diff << " " << diff << " " << max_diff << endl;
	  cout << "col: " << col << endl;
	  exit(-1);
	}

	if (true) {
	  TBox *tb=new TBox(Z-0.55,N-0.55,Z+0.55,N+0.55);
	  tb->SetFillStyle(1001);
	  tb->SetFillColor(col);
	  tb->SetLineColor(col);
	  tb->Draw();
	}
      }
    }

    tt.DrawLatex(70,-12,"Z");
    tt.DrawLatex(-12,110,"N");
    tt.SetTextAlign(12);
    tt.SetTextSize(tt.GetTextSize()/1.5);
    if (ik==0) {
      tt.DrawLatex(10,158,"Deviation in mass excess with");
      tt.DrawLatex(10,150,"finite-range droplet model");
      tt.DrawLatex(10,142,"(Moller et al. 1995)");
    } else if (ik==1) {
      tt.DrawLatex(10,155,"Deviation in mass excess with");
      tt.DrawLatex(10,145,"simple semi-empirical formula");
    } else if (ik==2) {
      tt.DrawLatex(10,155,"Deviation in mass excess with");
      tt.DrawLatex(10,145,"HFB-14");
    } else if (ik==3) {
      tt.DrawLatex(10,155,"Comparison to AME13");
    } else if (ik==4) {
      tt.DrawLatex(10,155,"Deviation in mass excess with");
      tt.DrawLatex(10,145,"Duflo and Zuker (1995)");
    } else {
      tt.DrawLatex(10,155,"Deviation in mass excess with");
      tt.DrawLatex(10,145,"Koura et al. (2005)");
    }
  
    tt.SetTextAlign(22);

    c1->Update();
	       
    //theApp.Run(kTRUE);
	       
    string fname=((string)"ex_nucmass_")+model+".eps";
    c1->Print(fname.c_str());
    fname=((string)"ex_nucmass_")+model+".png";
    c1->Print(fname.c_str());
    //system("convert ex_nucmass.eps ex_nucmass.png");

    delete c1;

  }

  return 0;
}
