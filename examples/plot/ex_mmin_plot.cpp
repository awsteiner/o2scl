/*
  Plot the results of ex_mmin.cpp
*/
#include <iostream>

#include <o2scl/table.h>
#include <o2scl/graph.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  cout.precision(5);

  // The ROOT objects
  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);

  if (true) {

    o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",
			   -2.0,-1.0,3.0,15.0,0,0,700,700,false,false);
    p1->SetLeftMargin(0.10);
    p1->SetRightMargin(0.05);
    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.10);

    static const size_t nfiles=6;

    string filenames[nfiles]={"../ex_mmin1.dat","../ex_mmin2.dat",
			      "../ex_mmin2g.dat","../ex_mmin3.dat",
			      "../ex_mmin3g.dat","../ex_mmin4.dat"};
    string labels[nfiles]={"Simplex","Fletcher-Reeves","FR with gradient",
			   "Polak-Ribere","PR with gradient","x","x"};
  
    tt.SetTextSize(tt.GetTextSize()/1.4);

    for(size_t k=0;k<nfiles;k++) {

      ovector ox, oy, oz, oret;
      double dtemp;
      ifstream fin;
      fin.open(filenames[k].c_str());
      while (fin >> dtemp) {
	ox.push_back(dtemp);
	fin >> dtemp;
	oy.push_back(dtemp);
	fin >> dtemp;
	oz.push_back(dtemp);
	fin >> dtemp;
	oret.push_back(dtemp);
      }
      fin.close();

      size_t n=ox.size();
      cout << n << endl;

      TGraph *g1=new TGraph(n);
      for(size_t i=0;i<n;i++) {
	g1->SetPoint(i,ox[i],oz[i]);
      }
      if (k==4) g1->SetLineColor(k+3);
      else if (k==6) g1->SetLineColor(k+2);
      else g1->SetLineColor(k+1);
      g1->SetLineStyle(5-k);
      g1->Draw();

      TLine *l1=new TLine(0.8,1.64+k*0.7,1.3,1.64+k*0.7);
      if (k==4) l1->SetLineColor(k+3);
      else l1->SetLineColor(k+1);
      l1->SetLineStyle(5-k);
      l1->Draw();
      tt.SetTextAlign(12);
      tt.DrawLatex(1.67,1.64+k*0.7,labels[k].c_str());
      tt.SetTextAlign(22);


    }

    tt.SetTextSize(tt.GetTextSize()*1.4);

    tt.DrawLatex(0.23,-2.14,"x");
    tt.DrawLatex(-2.35,6.88,"z");
    tt.DrawLatex(0.573,15.96,"Minimizer trajectories in x-z plane");

    //tt.DrawLatex(0.5,1.0e-15,"x");

    c1->Update();

    theApp.Run(kTRUE);

    c1->Print("ex_mmin.eps");
    system("convert ex_mmin.eps ex_mmin.png");

  } else {

    o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",
			   -2.0,-2.0,3.0,3.0,0,0,700,700,false,false);
    p1->SetLeftMargin(0.12);
    p1->SetRightMargin(0.05);
    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.12);

    static const size_t nfiles=5;

    string filenames[nfiles]={"../ex_mmin1.dat","../ex_mmin2.dat",
			      "../ex_mmin2g.dat","../ex_mmin3.dat",
			      "../ex_mmin3g.dat"};
    string labels[nfiles]={"Simplex","Fletcher-Reeves","FR with gradient",
			   "Polak-Ribere","PR with gradient"};
  
    tt.SetTextSize(tt.GetTextSize()/1.4);

    for(size_t k=0;k<nfiles;k++) {

      ovector ox, oy, oz, oret;
      double dtemp;
      ifstream fin;
      fin.open(filenames[k].c_str());
      while (fin >> dtemp) {
	ox.push_back(dtemp);
	fin >> dtemp;
	oy.push_back(dtemp);
	fin >> dtemp;
	oz.push_back(dtemp);
	fin >> dtemp;
	oret.push_back(dtemp);
      }
      fin.close();

      size_t n=ox.size();

      TGraph *g1=new TGraph(n);
      for(size_t i=0;i<n;i++) {
	g1->SetPoint(i,ox[i],oy[i]);
      }
      if (k==4) g1->SetLineColor(k+3);
      else g1->SetLineColor(k+1);
      g1->SetLineStyle(5-k);
      g1->Draw();

      TLine *l1=new TLine(0.6,1.7+k*0.25,1.3,1.7+k*0.25);
      if (k==4) l1->SetLineColor(k+3);
      else l1->SetLineColor(k+1);
      l1->SetLineStyle(5-k);
      l1->Draw();
      tt.SetTextAlign(12);
      tt.DrawLatex(1.42,1.7+k*0.25,labels[k].c_str());
      tt.SetTextAlign(22);

    }

    tt.SetTextSize(tt.GetTextSize()*1.4);

    tt.DrawLatex(0.45,-2.45,"x");
    tt.DrawLatex(-2.51,0.131,"y");

    c1->Update();

    theApp.Run(kTRUE);

    c1->Print("ex_mmin2.eps");
    system("convert ex_mmin2.eps ex_mmin2.png");

  }

  return 0;
}
