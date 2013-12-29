/*
  Plot the results of ex_ode.cpp
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

  table tab[8];
  hdf_file hf;

  // Load the files
  hf.open("../ex_ode.o2");
  for(size_t i=0;i<8;i++) {
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

  if (false) {

    // Compare the Cash-Karp and Prince-Dormand results for
    // the Bessel function

    o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",
			   0.0,1.0e-15,1.2,2.0e-4,0,0,700,700,
			   false,true);
    p1->SetLeftMargin(0.1);
    p1->SetRightMargin(0.05);
    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.1);

    // The Cash-Karp results
    size_t n0=tab[0].get_nlines();

    TGraph *g1=new TGraph(n0);
    TGraph *g2=new TGraph(n0);
    for(size_t i=0;i<n0;i++) {
      g1->SetPoint(i,tab[0]["x"][i],tab[0]["err"][i]);
      g2->SetPoint(i,tab[0]["x"][i],
		   fabs(tab[0]["calc"][i]-tab[0]["exact"][i]));
      cout << tab[0]["x"][i] << " "
	   << tab[0]["calc"][i] << " "
	   << tab[0]["exact"][i] << " "
	   << tab[0]["calc"][i]-tab[0]["exact"][i] << " " 
	   << tab[0]["err"][i] << endl;
    }
    g1->SetName("CK_est");
    g2->SetName("CK_act");
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
      g3->SetPoint(i,tab[1]["x"][i],tab[1]["err"][i]);
      g4->SetPoint(i,tab[1]["x"][i],
		   fabs(tab[1]["calc"][i]-tab[1]["exact"][i]));
      cout << tab[1]["x"][i] << " "
	   << tab[1]["calc"][i] << " "
	   << tab[1]["exact"][i] << " "
	   << tab[1]["calc"][i]-tab[1]["exact"][i] << " " 
	   << tab[1]["err"][i] << endl;
    }
    g3->SetName("PD_est");
    g4->SetName("PD_act");
    g3->SetLineStyle(3);
    g3->SetLineColor(kMagenta);
    g4->SetLineStyle(4);
    g4->SetLineColor(4);
    g3->Draw();
    g4->Draw();
    cout << endl;

    tt.DrawLatex(0.5,1.0e-15,"x");
    tt.DrawLatex(0.5,2.0e-10,"Bessel function with non-adaptive steppers");
    tt.DrawLatex(0.5,1.0e-10,"Cash-Karp Est. Error");
    tt.DrawLatex(0.5,1.0e-10,"Cash-Karp Act. Error");
    tt.DrawLatex(0.5,1.0e-10,"Prince-Dormand Est. Error");
    tt.DrawLatex(0.5,1.0e-10,"Prince-Dormand Act. Error");

    c1->Update();

    theApp.Run(kTRUE);

    //c1->Print("ex_ode_bessel.png");
    c1->Print("ex_ode_bessel.eps");
    system("convert ex_ode_bessel.eps ex_ode_bessel.png");
  }

  if (false) {
    
    o2scl_graph::new_graph(c2,p2,th2,"c2","cc2","d2",0,1.0e-18,1.2,3.0e-9,
			   0,0,700,700,false,true);
    p2->SetLeftMargin(0.1);
    p2->SetRightMargin(0.05);
    p2->SetTopMargin(0.1);
    p2->SetBottomMargin(0.1);
    
    size_t n2=tab[2].get_nlines();
    size_t n3=tab[3].get_nlines();
    
    TGraph *g1=new TGraph(n2);
    TGraph *g2=new TGraph(n2);
    for(size_t i=0;i<n2;i++) {
      g1->SetPoint(i,tab[2]["x"][i],fabs(tab[2]["err"][i]));
      g2->SetPoint(i,tab[2]["x"][i],
		   fabs(tab[2]["calc"][i]-tab[2]["exact"][i]));
    }
    g1->SetName("CK_est");
    g2->SetName("CK_act");
    g1->SetLineStyle(5);
    g1->SetLineColor(1);
    g2->SetLineStyle(2);
    g2->SetLineColor(2);
    g1->Draw();
    g2->Draw();

    TGraph *g3=new TGraph(n3);
    TGraph *g4=new TGraph(n3);
    for(size_t i=0;i<n3;i++) {
      g3->SetPoint(i,tab[3]["x"][i],fabs(tab[3]["err"][i]));
      g4->SetPoint(i,tab[3]["x"][i],
		   fabs(tab[3]["calc"][i]-tab[3]["exact"][i]));
    }
    g3->SetName("PD_est");
    g4->SetName("PD_act");
    g3->SetLineStyle(3);
    g3->SetLineColor(kMagenta);
    g4->SetLineStyle(4);
    g4->SetLineColor(4);
    g3->Draw();
    g4->Draw();

    tt.DrawLatex(0.5,2.0e-19,"x");
    tt.DrawLatex(0.5,2.0e-10,"Airy function with non-adaptive steppers");
    tt.DrawLatex(0.5,1.0e-10,"Cash-Karp Est. Error");
    tt.DrawLatex(0.5,1.0e-10,"Cash-Karp Act. Error");
    tt.DrawLatex(0.5,1.0e-10,"Prince-Dormand Est. Error");
    tt.DrawLatex(0.5,1.0e-10,"Prince-Dormand Act. Error");

    theApp.Run(kTRUE);

    c2->Print("ex_ode_airy.eps");
    system("convert ex_ode_airy.eps ex_ode_airy.png");
  }

  if (false) {

    o2scl_graph::new_graph(c3,p3,th3,"c3","cc3","d3",0,1.0e-10,9.0,3.0e-6,
			   0,0,700,700,false,true);
    p3->SetLeftMargin(0.1);
    p3->SetRightMargin(0.05);
    p3->SetTopMargin(0.1);
    p3->SetBottomMargin(0.1);
    
    size_t n4=tab[4].get_nlines();
    size_t n5=tab[5].get_nlines();
    size_t n6=tab[6].get_nlines();
    
    TGraph *g1=new TGraph(n4);
    //TGraph *g1b=new TGraph(n4);
    TGraph *g2=new TGraph(n4);
    for(size_t i=0;i<n4;i++) {
      g1->SetPoint(i,tab[4]["x"][i],fabs(tab[4]["err0"][i]));
      //g1b->SetPoint(i,tab[4]["x"][i],fabs(tab[4]["err1"][i]));
      g2->SetPoint(i,tab[4]["x"][i],
		   fabs(tab[4]["calc"][i]-tab[4]["exact"][i]));
    }
    g1->SetName("CK err0");
    //g1b->SetName("CK err1");
    g2->SetName("CK act");
    //g1b->SetLineStyle(2);
    g2->SetLineStyle(3);
    g1->Draw();
    //g1b->Draw();
    g2->Draw();

    TGraph *g3=new TGraph(n5);
    //TGraph *g3b=new TGraph(n5);
    TGraph *g4=new TGraph(n5);
    for(size_t i=0;i<n5;i++) {
      g3->SetPoint(i,tab[5]["x"][i],fabs(tab[5]["err0"][i]));
      //g3b->SetPoint(i,tab[5]["x"][i],fabs(tab[5]["err1"][i]));
      g4->SetPoint(i,tab[5]["x"][i],
		   fabs(tab[5]["calc"][i]-tab[5]["exact"][i]));
    }
    g3->SetName("CK2 err0");
    //g3b->SetName("CK2 err1");
    g4->SetName("CK2 act");
    g3->SetLineColor(2);
    //g3b->SetLineColor(2);
    g4->SetLineColor(2);
    //g3b->SetLineStyle(2);
    g4->SetLineStyle(3);
    g3->Draw();
    //g3b->Draw();
    g4->Draw();

    TGraph *g5=new TGraph(n6);
    //TGraph *g5b=new TGraph(n6);
    TGraph *g6=new TGraph(n6);
    for(size_t i=0;i<n6;i++) {
      g5->SetPoint(i,tab[6]["x"][i],fabs(tab[6]["err0"][i]));
      //g5b->SetPoint(i,tab[6]["x"][i],fabs(tab[6]["err1"][i]));
      g6->SetPoint(i,tab[6]["x"][i],
		   fabs(tab[6]["calc"][i]-tab[6]["exact"][i]));
    }
    g5->SetName("PD err0");
    //g5b->SetName("PD err1");
    g6->SetName("PD act");
    g5->SetLineColor(4);
    //g5b->SetLineColor(4);
    g6->SetLineColor(4);
    //g5b->SetLineStyle(2);
    g6->SetLineStyle(3);
    g5->Draw();
    //g5b->Draw();
    g6->Draw();

    tt.DrawLatex(0.5,2.0e-10,"x");
    tt.DrawLatex(0.5,2.0e-10,"Bessel function with adaptive steppers");
    tt.SetTextSize(tt.GetTextSize()/1.4);
    tt.DrawLatex(0.5,1.0e-10,"Cash-Karp Est. Error");
    tt.DrawLatex(0.5,1.0e-9,"Cash-Karp Act. Error");
    tt.SetTextColor(2);
    tt.DrawLatex(0.5,1.0e-8,"Cash-Karp(2) Est. Error");
    tt.DrawLatex(0.5,1.0e-7,"Cash-Karp(2) Act. Error");
    tt.SetTextColor(4);
    tt.DrawLatex(0.5,1.0e-6,"Prince-Dormand Est. Error");
    tt.DrawLatex(4.5,1.0e-7,"Prince-Dormand Act. Error");

    theApp.Run(kTRUE);

    c3->Print("ex_ode_bessel2.eps");
    system("convert ex_ode_bessel2.eps ex_ode_bessel2.png");
  }

  if (true) {
    
    o2scl_graph::new_graph(c4,p4,th4,"c4","cc4","d4",0,1.0e-11,9.0,3.0e-7,
			   0,0,700,700,false,true);
    p4->SetLeftMargin(0.1);
    p4->SetRightMargin(0.05);
    p4->SetTopMargin(0.1);
    p4->SetBottomMargin(0.1);
    
    size_t n7=tab[7].get_nlines();
    
    TGraph *g1=new TGraph(n7);
    TGraph *g2=new TGraph(n7);
    for(size_t i=0;i<n7;i++) {
      g1->SetPoint(i,tab[7]["x"][i],fabs(tab[7]["err0"][i]));
      g2->SetPoint(i,tab[7]["x"][i],
		   fabs(tab[7]["calc"][i]-tab[7]["exact"][i]));
    }
    g1->SetName("CK_est");
    g2->SetName("CK_act");
    g1->SetLineStyle(5);
    g1->SetLineColor(1);
    g2->SetLineStyle(2);
    g2->SetLineColor(2);
    g1->Draw();
    g2->Draw();

    tt.DrawLatex(0.5,4.0e-12,"x");
    tt.DrawLatex(0.5,2.0e-10,"Bessel function, high-level adaptive stepper");
    tt.DrawLatex(0.5,1.0e-10,"Est. Error");
    tt.DrawLatex(0.5,1.0e-10,"Act. Error");

    theApp.Run(kTRUE);

    c4->Print("ex_ode_bessel3.eps");
    system("convert ex_ode_bessel3.eps ex_ode_bessel3.png");
  }

  return 0;
}
