/*
  Plot the results of ex_rmf_nuc.cpp

  ----------------------------------------------------------------------  

  # occupied
  # proton - Tl 207
  elevel[] p_occ 1.0 41 5
  8.0     4.0      1.0     0.5    1G7/2       6.0 8.0120	3.474
  6.0     -3.0     1.0     0.5    2D5/2       7.0 8.0120	1.6827
  12.0    -6.0     1.0     0.5    1H5.5       6.0 8.0120	1.3481
  4.0     2.0      1.0     0.5    2D3/2       7.0 8.0120	0.351059
  2.0     -1.0     1.0     0.5    3S1/2       7.0 8.0120  0.0
  # neutron - Pb 207
  elevel[] n_occ 1.0 49 6
  10.0    5.0      1.0     -0.5   1H9/2       6.0	7.3674	3.477
  8.0     -4.0     1.0     -0.5   2F7/2       7.0	7.3674	2.339948
  14.0    -7.0     1.0     -0.5   1I6.5       7.0	7.3674	1.633368
  4.0     -2.0     1.0     -0.5   3P3/2       7.0	7.3674	0.89780
  6.0     3.0      1.0     -0.5   2F5/2       7.0	7.3674	0.569703
  2.0     1.0      1.0     -0.5   3P1/2       7.0	7.3674	0.0
  #
  # unoccupied
  # proton - Bi 209
  elevel[] p_unocc 1.0 41 5
  10.0    5.0      -4.0    0.5    1H9/2       6.0	3.7989	0.0
  8.0     -4.0     -3.0    0.5    2F7/2       7.0 3.7989	-0.896
  14.0    -7.0     -2.0    0.5    1I6.5       7.0	3.7989	-1.60858
  6.0     3.0      -0.5    0.5    2F5/2       7.0	3.7989	-2.82619
  4.0     -2.0     -1.0    0.5    3P3/2       7.0	3.7989	-3.11954
  # neutron - Pb 209
  elevel[] n_unocc 1.0 57 7
  10.0    -5.0     -4.0    -0.5   2G9/2       7.0	3.9372	0.0
  12.0    6.0      -3.0    -0.5   1I5.5       7.0	3.9372	-0.7788
  16.0    -8.0     -2.5    -0.5   1J7.5       7.0	3.9372	-1.423
  6.0     -3.0     -2.5    -0.5   3D5/2       7.0	3.9372	-1.56709
  2.0     -1.0     -2.0    -0.5   4S1/2       7.0	3.9372	-2.03222
  8.0     4.0      -1.5    -0.5   2G7/2       7.0	3.9372	-2.14943
  4.0     2.0      -1.5    -0.5   3D3/2       7.0	3.9372	-2.538

  ----------------------------------------------------------------------  

  double crcm 5.52573 
  double diff 0.0396098 
  double e_n1D3/2 -44.0517 
  double e_n1D5/2 -45.4611 
  double e_n1F5/2 -34.5258 
  double e_n1F7/2 -37.1125 
  double e_n1G7/2 -24.196 
  double e_n1G9/2 -28.2028 
  double e_n1H5.5 -19.0049 
  double e_n1H9/2 -13.5413 
  double e_n1I6.5 -9.74572 
  double e_n1P1/2 -52.3268 
  double e_n1P3/2 -52.9166 
  double e_n1S1/2 -59.1053 
  double e_n2D3/2 -18.9837 
  double e_n2D5/2 -20.6294 
  double e_n2F5/2 -9.15378 
  double e_n2F7/2 -11.1658 
  double e_n2P1/2 -29.5534 
  double e_n2P3/2 -30.613 
  double e_n2S1/2 -40.958 
  double e_n3P1/2 -7.65713 
  double e_n3P3/2 -8.4335 
  double e_n3S1/2 -18.2311 
  double e_p1D3/2 -33.8697 
  double e_p1D5/2 -35.424 
  double e_p1F5/2 -24.6905 
  double e_p1F7/2 -27.4525 
  double e_p1G7/2 -14.6596 
  double e_p1G9/2 -18.8422 
  double e_p1H5.5 -9.8529 
  double e_p1P1/2 -41.7417 
  double e_p1P3/2 -42.421 
  double e_p1S1/2 -48.0086 
  double e_p2D3/2 -8.84736 
  double e_p2D5/2 -10.471 
  double e_p2P1/2 -19.1468 
  double e_p2P3/2 -20.2037 
  double e_p2S1/2 -30.1384 
  double e_p3S1/2 -7.72394 
  double eb -7.80877 
  double rch 5.52919 
  double rnrms 5.74303 
  double rnrp 0.273457 
  double rprms 5.46957 
  double stens 0.00956719 
  double xnu 126 

  ----------------------------------------------------------------------  

  #include "../lib/params.h"
  #include "../lib/constants.h"
  #include "../lib/graph.h"
  #include "TGaxis.h"

  using namespace std;
  using namespace aws_const;
  using namespace root_graph;

  class elevel : public aio {
  public:
  double d1, d2, d3, isospin;
  string level;
  double d4, ebase, eshift;
  
  elevel() : aio(1) {
  if (added==false) {
  addtype();
  added=true;
  }
  };
  static bool added;
  
  virtual void *create(int sz) {
  elevel *ns=new elevel[sz];
  return (void *)ns;
  };
  virtual int in(std::istream *in, void *vp, int sz);
  virtual int out(std::ostream *out, void *vp, int sz);
  virtual int remove(void *vp) {
  elevel *el=(elevel *)vp;
  delete[] el;
  return 0;
  };
  virtual std::string type() { return "elevel"; };
  virtual std::string version() { return "1.0"; };
  virtual int nw(void *vp, int sz, int sz2) { return 8*sz*sz2; };

  };

  bool elevel::added=false;

  int elevel::in(istream *ins, void *vp, int sz) {
  elevel *el=(elevel *)vp+sz;
  doublein(ins,el->d1);
  doublein(ins,el->d2);
  doublein(ins,el->d3);
  doublein(ins,el->isospin);
  wordin(ins,el->level);
  doublein(ins,el->d4);
  doublein(ins,el->ebase);
  doublein(ins,el->eshift);
  return 0;
  }

  int elevel::out(ostream *outs, void *vp, int sz) {
  elevel *el=(elevel *)vp+sz;
  doubleout(outs,el->d1);
  doubleout(outs,el->d2);
  doubleout(outs,el->d3);
  doubleout(outs,el->isospin);
  wordout(outs,el->level);
  doubleout(outs,el->d4);
  doubleout(outs,el->ebase);
  doubleout(outs,el->eshift);
  return 0;
  }

  int main(void) {
  elevel tst, *n_occ, *p_occ, *n_unocc, *p_unocc;
  int i, nno, npo, nnu, npu;
  void *vp;

  double e_p1S12, e_p1P32, e_p1P12, e_p1D52, e_p2S12, e_p1D32;
  double e_p1F72, e_p2P32, e_p1F52, e_p2P12, e_p1G92, e_p1G72;
  double e_p2D52, e_p2D32, e_p1H55, e_p3S12, e_n1S12, e_n1P32;
  double e_n1P12, e_n1D52, e_n2S12, e_n1D32, e_n1F72, e_n2P32;
  double e_n1F52, e_n2P12, e_n1G92, e_n1G72, e_n1H55, e_n2D52;
  double e_n2D32, e_n1H92, e_n3S12, e_n2F72, e_n3P32, e_n2F52;
  double e_n3P12, e_n1I65;

  double eu_p1H92, eu_p2F72, eu_p1I65, eu_p2F52, eu_p3P32;
  double eu_n2G92, eu_n1I55, eu_n1J75, eu_n3D52, eu_n4S12;
  double eu_n2G72, eu_n3D32;

  params *pa=params_ns::new_params();
  pa->fin("pblevel.dat");
  pa->get("n_occ",vp,nno);
  n_occ=(elevel *)vp;
  pa->get("p_occ",vp,npo);
  p_occ=(elevel *)vp;
  pa->get("n_unocc",vp,nnu);
  n_unocc=(elevel *)vp;
  pa->get("p_unocc",vp,npu);
  p_unocc=(elevel *)vp;
  
  TApplication theApp("App",0,NULL);

  TCanvas *c1;
  TPad *p1;
  TH1 *th1;
  c1=new TCanvas("c1","Energy levels",0,0,700,700);
  c1->SetFillColor(10);
  p1=new TPad("p1","",0.02,0.02,1.0,1.0);
  p1->SetFillColor(10);
  p1->SetFrameLineColor(10);
  p1->Draw();
  p1->cd();
  th1=p1->DrawFrame(0.0,-12.0,11.0,0.0);
  th1->GetXaxis()->SetLabelFont(132);
  th1->GetYaxis()->SetLabelFont(132);
  th1->GetXaxis()->CenterTitle(kTRUE);
  th1->GetYaxis()->CenterTitle(kTRUE);
  th1->GetXaxis()->SetTitleFont(132);
  th1->GetYaxis()->SetTitleFont(132);

  th1->GetXaxis()->SetTitleColor(10);
  th1->GetXaxis()->SetAxisColor(10);
  th1->GetXaxis()->SetLabelColor(10);

  TLatex tt;
  TLine *l1;
  tt.SetTextAlign(12);
  tt.SetTextSize(tt.GetTextSize()/2.0);
  tt.SetTextFont(132);

  params *pa2=params_ns::new_params();
  pa2->fin("../rnrp/hpnew/elsum.out");

  // Occupied neutron levels

  for(i=0;i<nno;i++) {
  l1=new TLine(3.0,-n_occ[i].ebase-n_occ[i].eshift,4.0,
  -n_occ[i].ebase-n_occ[i].eshift);
  l1->Draw();
  tt.DrawLatex(4.0,-n_occ[i].ebase-n_occ[i].eshift,
  n_occ[i].level.c_str());
  }

  pa2->get("e_n1S1/2",e_n1S12);
  pa2->get("e_n1P3/2",e_n1P32);
  pa2->get("e_n1P1/2",e_n1P12);
  pa2->get("e_n1D5/2",e_n1D52);
  pa2->get("e_n2S1/2",e_n2S12);
  pa2->get("e_n1D3/2",e_n1D32);
  pa2->get("e_n1F7/2",e_n1F72);
  pa2->get("e_n2P3/2",e_n2P32);
  pa2->get("e_n1F5/2",e_n1F52);
  pa2->get("e_n2P1/2",e_n2P12);
  pa2->get("e_n1G9/2",e_n1G92);
  pa2->get("e_n1G7/2",e_n1G72);
  pa2->get("e_n1H5.5",e_n1H55);
  pa2->get("e_n2D5/2",e_n2D52);
  pa2->get("e_n2D3/2",e_n2D32);
  pa2->get("e_n1H9/2",e_n1H92);
  pa2->get("e_n3S1/2",e_n3S12);
  pa2->get("e_n2F7/2",e_n2F72);
  pa2->get("e_n3P3/2",e_n3P32);
  pa2->get("e_n2F5/2",e_n2F52);
  pa2->get("e_n3P1/2",e_n3P12);
  pa2->get("e_n1I6.5",e_n1I65);

  l1=new TLine(1.0,e_n1S12,2.0,e_n1S12);
  l1->Draw();
  l1=new TLine(1.0,e_n1P32,2.0,e_n1P32);
  l1->Draw();
  l1=new TLine(1.0,e_n1P12,2.0,e_n1P12);
  l1->Draw();
  l1=new TLine(1.0,e_n1D52,2.0,e_n1D52);
  l1->Draw();
  l1=new TLine(1.0,e_n2S12,2.0,e_n2S12);
  l1->Draw();
  l1=new TLine(1.0,e_n1D32,2.0,e_n1D32);
  l1->Draw();
  l1=new TLine(1.0,e_n1F72,2.0,e_n1F72);
  l1->Draw();
  l1=new TLine(1.0,e_n2P32,2.0,e_n2P32);
  l1->Draw();
  l1=new TLine(1.0,e_n1F52,2.0,e_n1F52);
  l1->Draw();
  l1=new TLine(1.0,e_n2P12,2.0,e_n2P12);
  l1->Draw();
  l1=new TLine(1.0,e_n1G92,2.0,e_n1G92);
  l1->Draw();
  l1=new TLine(1.0,e_n1G72,2.0,e_n1G72);
  l1->Draw();
  l1=new TLine(1.0,e_n1H55,2.0,e_n1H55);
  l1->Draw();
  l1=new TLine(1.0,e_n2D52,2.0,e_n2D52);
  l1->Draw();
  l1=new TLine(1.0,e_n2D32,2.0,e_n2D32);
  l1->Draw();
  l1=new TLine(1.0,e_n1H92,2.0,e_n1H92);
  l1->Draw();
  l1=new TLine(1.0,e_n3S12,2.0,e_n3S12);
  l1->Draw();
  l1=new TLine(1.0,e_n2F72,2.0,e_n2F72);
  l1->Draw();
  l1=new TLine(1.0,e_n3P32,2.0,e_n3P32);
  l1->Draw();
  l1=new TLine(1.0,e_n2F52,2.0,e_n2F52);
  l1->Draw();
  l1=new TLine(1.0,e_n3P12,2.0,e_n3P12);
  l1->Draw();
  l1=new TLine(1.0,e_n1I65,2.0,e_n1I65);
  l1->Draw();

  l1=new TLine(2.0,e_n1H92,3.0,
  -n_occ[0].ebase-n_occ[0].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,e_n2F72,3.0,
  -n_occ[1].ebase-n_occ[1].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,e_n1I65,3.0,
  -n_occ[2].ebase-n_occ[2].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,e_n3P32,3.0,
  -n_occ[3].ebase-n_occ[3].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,e_n2F52,3.0,
  -n_occ[4].ebase-n_occ[4].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,e_n3P12,3.0,
  -n_occ[5].ebase-n_occ[5].eshift);
  l1->SetLineStyle(4);
  l1->Draw();

  // Occupied neutron levels

  for(i=0;i<npo;i++) {
  l1=new TLine(9.0,-p_occ[i].ebase-p_occ[i].eshift,10.0,
  -p_occ[i].ebase-p_occ[i].eshift);
  l1->Draw();
  tt.DrawLatex(10.0,-p_occ[i].ebase-p_occ[i].eshift,
  p_occ[i].level.c_str());
  }

  pa2->get("e_p1S1/2",e_p1S12);
  pa2->get("e_p1P3/2",e_p1P32);
  pa2->get("e_p1P1/2",e_p1P12);
  pa2->get("e_p1D5/2",e_p1D52);
  pa2->get("e_p2S1/2",e_p2S12);
  pa2->get("e_p1D3/2",e_p1D32);
  pa2->get("e_p1F7/2",e_p1F72);
  pa2->get("e_p2P3/2",e_p2P32);
  pa2->get("e_p1F5/2",e_p1F52);
  pa2->get("e_p2P1/2",e_p2P12);
  pa2->get("e_p1G9/2",e_p1G92);
  pa2->get("e_p1G7/2",e_p1G72);
  pa2->get("e_p2D5/2",e_p2D52);
  pa2->get("e_p2D3/2",e_p2D32);
  pa2->get("e_p1H5.5",e_p1H55);
  pa2->get("e_p3S1/2",e_p3S12);

  l1=new TLine(7.0,e_p1S12,8.0,e_p1S12);
  l1->Draw();
  l1=new TLine(7.0,e_p1P32,8.0,e_p1P32);
  l1->Draw();
  l1=new TLine(7.0,e_p1P12,8.0,e_p1P12);
  l1->Draw();
  l1=new TLine(7.0,e_p1D52,8.0,e_p1D52);
  l1->Draw();
  l1=new TLine(7.0,e_p2S12,8.0,e_p2S12);
  l1->Draw();
  l1=new TLine(7.0,e_p1D32,8.0,e_p1D32);
  l1->Draw();
  l1=new TLine(7.0,e_p1F72,8.0,e_p1F72);
  l1->Draw();
  l1=new TLine(7.0,e_p2P32,8.0,e_p2P32);
  l1->Draw();
  l1=new TLine(7.0,e_p1F52,8.0,e_p1F52);
  l1->Draw();
  l1=new TLine(7.0,e_p2P12,8.0,e_p2P12);
  l1->Draw();
  l1=new TLine(7.0,e_p1G92,8.0,e_p1G92);
  l1->Draw();
  l1=new TLine(7.0,e_p1G72,8.0,e_p1G72);
  l1->Draw();
  l1=new TLine(7.0,e_p2D52,8.0,e_p2D52);
  l1->Draw();
  l1=new TLine(7.0,e_p2D32,8.0,e_p2D32);
  l1->Draw();
  l1=new TLine(7.0,e_p1H55,8.0,e_p1H55);
  l1->Draw();
  l1=new TLine(7.0,e_p3S12,8.0,e_p3S12);
  l1->Draw();

  l1=new TLine(8.0,e_p1G72,9.0,
  -p_occ[0].ebase-p_occ[0].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,e_p2D52,9.0,
  -p_occ[1].ebase-p_occ[1].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,e_p1H55,9.0,
  -p_occ[2].ebase-p_occ[2].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,e_p2D32,9.0,
  -p_occ[3].ebase-p_occ[3].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,e_p3S12,9.0,
  -p_occ[4].ebase-p_occ[4].eshift);
  l1->SetLineStyle(4);
  l1->Draw();

  // Unoccupied neutron levels

  for(i=0;i<nnu;i++) {
  if (i==5 || i==3) {
  l1=new TLine(3.0,-n_unocc[i].ebase-n_unocc[i].eshift,4.8,
  -n_unocc[i].ebase-n_unocc[i].eshift);
  l1->Draw();
  tt.DrawLatex(4.85,-n_unocc[i].ebase-n_unocc[i].eshift,
  n_unocc[i].level.c_str());
  } else if (i==4 || i==2) {
  l1=new TLine(3.0,-n_unocc[i].ebase-n_unocc[i].eshift,4.0,
  -n_unocc[i].ebase-n_unocc[i].eshift);
  l1->Draw();
  tt.DrawLatex(4.0,-n_unocc[i].ebase-n_unocc[i].eshift-0.05,
  n_unocc[i].level.c_str());
  } else {
  l1=new TLine(3.0,-n_unocc[i].ebase-n_unocc[i].eshift,4.0,
  -n_unocc[i].ebase-n_unocc[i].eshift);
  l1->Draw();
  tt.DrawLatex(4.0,-n_unocc[i].ebase-n_unocc[i].eshift,
  n_unocc[i].level.c_str());
  }
  }

  pa2->get("eu_n2G9/2",eu_n2G92);
  pa2->get("eu_n1I5.5",eu_n1I55);
  pa2->get("eu_n1J7.5",eu_n1J75);
  pa2->get("eu_n3D5/2",eu_n3D52);
  pa2->get("eu_n4S1/2",eu_n4S12);
  pa2->get("eu_n2G7/2",eu_n2G72);
  pa2->get("eu_n3D3/2",eu_n3D32);

  l1=new TLine(1.0,eu_n2G92,2.0,eu_n2G92);
  l1->Draw();
  l1=new TLine(1.0,eu_n1I55,2.0,eu_n1I55);
  l1->Draw();
  l1=new TLine(1.0,eu_n1J75,2.0,eu_n1J75);
  l1->Draw();
  l1=new TLine(1.0,eu_n3D52,2.0,eu_n3D52);
  l1->Draw();
  l1=new TLine(1.0,eu_n4S12,2.0,eu_n4S12);
  l1->Draw();
  l1=new TLine(1.0,eu_n2G72,2.0,eu_n2G72);
  l1->Draw();
  l1=new TLine(1.0,eu_n3D32,2.0,eu_n3D32);
  l1->Draw();

  l1=new TLine(2.0,eu_n2G92,3.0,
  -n_unocc[0].ebase-n_unocc[0].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n1I55,3.0,
  -n_unocc[1].ebase-n_unocc[1].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n1J75,3.0,
  -n_unocc[2].ebase-n_unocc[2].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n3D52,3.0,
  -n_unocc[3].ebase-n_unocc[3].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n4S12,3.0,
  -n_unocc[4].ebase-n_unocc[4].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n2G72,3.0,
  -n_unocc[5].ebase-n_unocc[5].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(2.0,eu_n3D32,3.0,
  -n_unocc[6].ebase-n_unocc[6].eshift);
  l1->SetLineStyle(4);
  l1->Draw();

  // Unoccupied proton levels

  for(i=0;i<npu;i++) {
  l1=new TLine(9.0,-p_unocc[i].ebase-p_unocc[i].eshift,10.0,
  -p_unocc[i].ebase-p_unocc[i].eshift);
  l1->Draw();
  tt.DrawLatex(10.0,-p_unocc[i].ebase-p_unocc[i].eshift,
  p_unocc[i].level.c_str());
  }

  pa2->get("eu_p1H9/2",eu_p1H92);
  pa2->get("eu_p2F7/2",eu_p2F72);
  pa2->get("eu_p1I6.5",eu_p1I65);
  pa2->get("eu_p2F5/2",eu_p2F52);
  pa2->get("eu_p3P3/2",eu_p3P32);

  l1=new TLine(7.0,eu_p1H92,8.0,eu_p1H92);
  l1->Draw();
  l1=new TLine(7.0,eu_p2F72,8.0,eu_p2F72);
  l1->Draw();
  l1=new TLine(7.0,eu_p1I65,8.0,eu_p1I65);
  l1->Draw();
  l1=new TLine(7.0,eu_p2F52,8.0,eu_p2F52);
  l1->Draw();
  l1=new TLine(7.0,eu_p3P32,8.0,eu_p3P32);
  l1->Draw();

  l1=new TLine(8.0,eu_p1H92,9.0,
  -p_unocc[0].ebase-p_unocc[0].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,eu_p2F72,9.0,
  -p_unocc[1].ebase-p_unocc[1].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,eu_p1I65,9.0,
  -p_unocc[2].ebase-p_unocc[2].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,eu_p2F52,9.0,
  -p_unocc[3].ebase-p_unocc[3].eshift);
  l1->SetLineStyle(4);
  l1->Draw();
  l1=new TLine(8.0,eu_p3P32,9.0,
  -p_unocc[4].ebase-p_unocc[4].eshift);
  l1->SetLineStyle(4);
  l1->Draw();

  // Final stuff

  tt.SetTextAlign(22);
  tt.SetTextSize(tt.GetTextSize()*1.2);
  tt.DrawLatex(5.5,0.5,"Pb Energy Levels");
  tt.DrawLatex(1.5,-12.5,"Calculated");
  tt.DrawLatex(3.5,-12.5,"Experimental");
  tt.DrawLatex(7.5,-12.5,"Calculated");
  tt.DrawLatex(9.5,-12.5,"Experimental");

  theApp.Run(kTRUE);
  c1->Print("level.eps");
  c1->Print("level.C");
  //  system("scp level.eps stein@physics.umn.edu:html");
	    
  params_ns::free_params(pa);
  return 0;
  }

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
  hf.open("../ex_rmf_nuc_prof.o2");
  hdf_input(hf,tu);
  hf.close();

  // The ROOT objects
  TApplication theApp("App",0,NULL);

  TCanvas *c1, *c2, *c3, *c4;
  TPad *p1, *p2, *p3, *p4;
  TH1 *th1, *th2, *th3, *th4;

  TLatex tt;
  tt.SetTextAlign(22);
  tt.SetTextFont(132);

  // The experimental charge density for Lead 208 from de Vries, et
  // al. At. Data Nucl. Data Tables 36 (1987) 495 using the
  // Sum-of-Gaussians method
  double rr[12]={0.1,0.7,1.6,2.1,2.7,3.5,4.2,
                 5.1,6.0,6.6,7.6,8.7};
  double qq[12]={0.003845,0.009724,0.033093,0.000120,
                 0.083107,0.080869,0.139957,0.260892,
                 0.336013,0.0033637,0.018729,0.000020};
  double g=1.7/sqrt(1.5);
  double a[12];
  for(size_t i=0;i<12;i++) {
    a[i]=82.0*qq[i]/2.0/pow(o2scl_const::pi,1.5)/
      pow(g,3.0)/(1.0+2.0*rr[i]*rr[i]/g/g);
  }
  
  // Add experimental profile to table
  tu.new_column("data");
  for(size_t i=0;i<tu.get_nlines();i++) {
    double val=0.0;
    for(size_t j=0;j<12;j++) {
      val+=a[j]*(exp(-pow((tu.get("r",i)-rr[j])/g,2.0))+
		 exp(-pow((tu.get("r",i)+rr[j])/g,2.0)));
    }
    tu.set("data",i,val);
  }
  
  o2scl_graph::new_graph(c1,p1,th1,"c1","cc1","d1",0,0,12.0,0.1,
			 0,0,700,700,false,false);
  p1->SetLeftMargin(0.13);
  p1->SetRightMargin(0.02);
  p1->SetTopMargin(0.04);
  p1->SetBottomMargin(0.09);
  p1->cd();
  c1->Update();

  o2scl_graph::table_graph(tu,"r","data",1,1);
  o2scl_graph::table_graph(tu,"r","rhop",2,2);
  o2scl_graph::table_graph(tu,"r","rhon",3,4);

  tt.DrawLatex(6.0,-0.007,"fm");
  tt.SetTextAngle(90);
  tt.DrawLatex(-1.4,0.05,"Number density (fm^{-3})");
  tt.SetTextAngle(0);
  tt.DrawLatex(4.38,0.090914,"Neutrons");
  tt.DrawLatex(3.458824,0.05687603,"Protons");
  tt.DrawLatex(3.174848,0.06645457,"Protons (expt.)");
  
  c1->Update();
	       
  //theApp.Run(kTRUE);
  
  c1->Print("ex_nuc_prof.eps");
  c1->Print("ex_nuc_prof.png");
    
  delete c1;

  return 0;
}
