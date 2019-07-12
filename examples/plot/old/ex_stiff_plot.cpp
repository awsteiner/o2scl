/*
  Plot the results of ex_stiff.cpp
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

  table tab[2];
  hdf_file hf;

  // Load the files
  hf.open("../ex_stiff.o2");
  for(size_t i=0;i<2;i++) {
    hdf_input(hf,tab[i],((string)"table_")+itos(i));
  }
  hf.close();

  // The ROOT objects
  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);
  //tt.SetTextSize(tt.GetTextSize()*1.2);

  // Compare the Cash-Karp and Prince-Dormand results for
  // the Bessel function

  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",
			 0.0,1.0e-20,4.0,1.0e17,0,0,700,700,
			 false,true);
  p1->SetLeftMargin(0.1);
  p1->SetRightMargin(0.05);
  p1->SetTopMargin(0.05);
  p1->SetBottomMargin(0.1);

  // The Cash-Karp results
  size_t n0=tab[0].get_nlines();

  TGraph *g1=new TGraph(n0);
  TGraph *g2=new TGraph(n0);
  for(size_t i=0;i<n0;i++) {
    g1->SetPoint(i,tab[0]["x"][i],fabs(tab[0]["rel_err"][i]));
    g2->SetPoint(i,tab[0]["x"][i],fabs(tab[0]["rel_diff"][i]));
    cout << tab[0]["x"][i] << " ";
    cout << fabs(tab[0]["rel_err"][i]) << " ";
    cout << fabs(tab[0]["rel_diff"][i]) << endl;
  }
  g1->SetName("Bsimp_est");
  g2->SetName("Bsimp_act");
  g1->SetLineStyle(5);
  g1->SetLineColor(1);
  g2->SetLineStyle(2);
  g2->SetLineColor(2);
  g1->Draw();
  g2->Draw();
  cout << endl;

  // The Prince-Dormand results

  size_t n1=tab[1].get_nlines();

  TGraph *g3=new TGraph(n1);
  TGraph *g4=new TGraph(n1);
  for(size_t i=0;i<n1;i++) {
    g3->SetPoint(i,tab[1]["x"][i],fabs(tab[1]["rel_err"][i]));
    g4->SetPoint(i,tab[1]["x"][i],fabs(tab[1]["rel_diff"][i]));
    cout << tab[1]["x"][i] << " ";
    cout << fabs(tab[1]["rel_err"][i]) << " ";
    cout << fabs(tab[1]["rel_diff"][i]) << endl;
  }
  g3->SetName("CK_est");
  g4->SetName("CK_act");
  g3->SetLineStyle(3);
  g3->SetLineColor(kMagenta);
  g4->SetLineStyle(4);
  g4->SetLineColor(4);
  g3->Draw();
  g4->Draw();
  cout << endl;

  tt.DrawLatex(0.5,1.0e-15,"x");
  tt.DrawLatex(0.5,1.0e-10,"Bader-Deuflhard Est. Rel. Error");
  tt.DrawLatex(0.5,1.0e-5,"Bader-Deuflhard Act. Rel. Error");
  tt.DrawLatex(0.5,1.0e5,"Simple Adaptive Est. Rel. Error");
  tt.DrawLatex(0.5,1.0e10,"Simple Adaptive Act. Rel. Error");

  c1->Update();

  theApp.Run(kTRUE);

  //c1->Print("ex_stiff.png");
  c1->Print("ex_stiff.eps");
  system("convert ex_stiff.eps ex_stiff.png");

  return 0;
}
