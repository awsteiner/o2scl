/*
  Ugly hack.
*/
#include <iostream>
#include <o2scl/base_ioc.h>
#include <o2scl/table.h>
#include <o2scl/constants.h>
#include <o2scl/user_io.h>
#include <o2scl_ext/table_fp.h>
#include <o2scl/graph.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_ext;
using namespace o2scl_graph;

int main(void) {
  cout.setf(ios::scientific);

  table *tab1, *tab2;

  collection co;
  base_ioc bio;
  text_in_file *tif;
  void *vp=0;
  string name;
  int sz, sz2;

  system("echo c1 c2 c3 c4 > exftemp1");
  system(((std::string)("cd ..; cat ex_fptr.scr | grep")+
	  " -v terat | grep -v value >> plot/exftemp1").c_str());
  system("acol -generic exftemp1 -internal exftemp1.o2");
  system(((std::string)("cd ..; acol -generic ex_fptr.out ")+
	  "-internal plot/exftemp2.o2").c_str());
  
  o2scl_input_text("exftemp1.o2",tab1);
  o2scl_input_text("exftemp2.o2",tab2);
  cout << "Table 1: " << tab1->get_nlines() << endl;
  cout << "Table 2: " << tab2->get_nlines() << endl;

  TApplication theApp("App",0,NULL);
  TCanvas *c1, *c2;
  TPad *p1, *p2;
  TH1 *th1, *th2;
  TGaxis *ax1, *ay1, *ax2, *ay2;
  TLatex tt;
  
  tt.SetTextAlign(22);
  tt.SetTextFont(132);
  //tt.SetTextSize(tt.GetTextSize()*1.2);
  
  new_graph_ticks(c1,p1,th1,"c1","d1","e1",-1,-3,2,3,ax1,ay1,
		  0,0,700,700,false,false);
  th1->GetXaxis()->SetLabelSize(th1->GetXaxis()->GetLabelSize()*1.3);
  th1->GetYaxis()->SetLabelSize(th1->GetYaxis()->GetLabelSize()*1.3);
  p1->SetTopMargin(0.04);
  p1->SetBottomMargin(0.12);
  p1->SetLeftMargin(0.12);
  p1->SetRightMargin(0.04);
  
  table_fp tb1=*tab1;
  table_fp tb2=*tab2;

  TGraph *gr1=table_graph(&tb2,"x","y",1,2);
  gr1->Draw();

  for(size_t i=0;i<tab1->get_nlines();i++) {
    TMarker *m1=new TMarker(tab1->get("c1",i),tab1->get("c2",i),
			    m_fill_square);
    m1->SetMarkerColor(4);
    m1->Draw();
    tt.DrawLatex(tab1->get("c1",i),tab1->get("c2",i)+0.2,
		 itos(i+1).c_str());
  }

  c1->Update();
  c1->Print("ex_fptr_plot.eps");
  
  system("rm -f exftemp1");
  system("rm -f exftemp1.o2");
  system("rm -f exftemp2.o2");
  system("convert ex_fptr_plot.eps ex_fptr_plot.png");

  return 0;
}
