#include <iostream>
#include <fstream>

#include <o2scl/table3d.h>
#include <o2scl/graph.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(void) {

  // Output to file
  table3d t;
  hdf_file hf;
  hf.open("../ex_mandel.o2");
  hdf_input(hf,t);
  hf.close();

  double delta=t.get_constant("delta");
  double minx=t.get_constant("minx");
  double maxx=t.get_constant("maxx");
  double miny=t.get_constant("miny");
  double maxy=t.get_constant("maxy");
  size_t nx=t.get_nx();
  size_t ny=t.get_ny();
  size_t maxcol=((size_t)(t.get_constant("maxtime")+1.0e-4));

  TApplication theApp("App",0,NULL);
  TCanvas *c1;
  TPad *p1;
  TH1 *th1;
  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","p1",
			 minx-0.0*delta/2,miny-delta/2,
			 maxx+0.0*delta/2,maxy+delta/2,
                         0,0,700,700);
  p1->SetLeftMargin(0);
  p1->SetRightMargin(0);
  p1->SetTopMargin(0);
  p1->SetBottomMargin(0);
  
  for(size_t i=1;i<=maxcol;i++) {
    TColor *c=(TColor *)(gROOT->GetListOfColors()->At(i));
    double cx=sqrt(((double)(i-1))/((double)(maxcol-1)));
    c->SetRGB(1.0-cx,1.0-cx,1.0);
  }
  
  for(size_t i=0;i<nx;i++) {
    for(size_t j=0;j<ny;j++) {
      double x=t.get_grid_x(i);
      double y=t.get_grid_y(j);
      TBox *b=new TBox(x-delta/2,y-delta/2,x+delta/2,y+delta/2);
      b->SetFillColor(((int)(t.get(i,j,"time")+1.0e-4)));
      b->SetFillStyle(1001);
      b->Draw();
    }
  }

  c1->Update();
  c1->Print("ex_mandel_plot.png");

  return 0;
}
