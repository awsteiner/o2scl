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
#include "agraphm.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_acol;
using namespace o2scl_agraph;

agraph_manager::agraph_manager() {

  xlo=0.0;
  xhi=-1.0;
  ylo=0.0;
  yhi=-1.0;
  zlo=0.0;
  zhi=-1.0;

  xtitle="";
  ytitle="";
  wtitle="";
  wtitle2="";
  ptitle="";
  logx=false;
  logy=false;
  logz=false;
  prmar=0.05;
  ptmar=0.05;
  plmar=0.11;
  pbmar=0.1;

  env_var_name="AGRAPH_DEFAULTS";

  int_pad=0;
  int_canvas=0;
  int_th1=0;
  top_axis=0;
  right_axis=0;
  theApp=0;

  line_ix=0;
  marker_ix=2;
  color_ix=0;

  int colors[9]={1,kRed,kGreen,kBlue,kMagenta+2,kCyan+2,
		 kOrange,kSpring-7,kViolet+5};
  std_colors.resize(9);
  vector_copy(9,colors,std_colors);

  lwidth=1;
  multi_bins=30;
  talign=22;
  tt.SetTextFont(132);
  tt.SetTextAlign(talign);

  surf_hist=0;

  xset=false;
  yset=false;
  zset=false;
}

agraph_manager::~agraph_manager() {
  if (int_th1!=0) {
    cout << "Deleting axis." << endl;
    delete int_th1;
    int_th1=0;
  }
  if (top_axis!=0) {
    cout << "Deleting top axis." << endl;
    delete top_axis;
    top_axis=0;
  }
  if (right_axis!=0) {
    cout << "Deleting right axis." << endl;
    delete right_axis;
    right_axis=0;
  }
  if (int_pad!=0) {
    cout << "Deleting pad." << endl;
    delete int_pad;
    int_pad=0;
  }
  if (int_canvas!=0) {
    cout << "Deleting canvas." << endl;
    delete int_canvas;
    int_canvas=0;
  }
}

int agraph_manager::setup_cli() {
  
  //---------------------------------------
  // Get HOME directory and command history
  
  char *hd=getenv("HOME");
  std::string histfile;
  histfile=hd;
  histfile+="/.agraph_hist";
  
  //---------------------------------------
  // Specify command-line option object
  
  //#ifdef O2SCL_READLINE
  cl=new cli_readline(histfile);
  //#else
  //cl=new cli;
  //#endif
  
  return 0;
}

int agraph_manager::setup_options() {

  acol_manager::setup_options();

  int cl_param=cli::comm_option_cl_param;
  int both=cli::comm_option_both;

  static const size_t narr=15;
  comm_option_s options_arr[narr]={
    {0,"text","Output text to a specified location.",1,4,
     "<string> [x y] [size]",
     ((string)"Add a text item to the current list of text items. ")+
     "The first parameter is the text string, the second and third "+
     "parameters are the x- and y-values of the center of the text "+
     "item, and the fourth parameter is the size of the text item "+
     "relative to the default. If the locations are not specified, "+
     "the text items are put in the center of the plot. To include "+
     "spaces in the text, enclose the string in double quotes.",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_text),both},
    {0,"density-colors","Set density colors for 'den-plot'.",
     2,-1,"<type> <n_colors> [...]",
     ((string)"Set color palette to use in ")+
     "den-plot. The <type> parameter can be either 'rgb', 'hsv', "+
     "'rainbow', or 'gradients'. For 'rgb', the 3rd, 4th, and 5th "+
     "arguments are functions which specify the red, green and blue "+
     "components as a function of x=[0,1]. The maximum value of each "+
     "color component is 1. For 'hsv', the 3rd, 4th, and 5th arguments "+
     "are functions which specify the hue, saturation, and value as "+
     "a function of x=[0,1]. The maximum values of HSV are 360, 1, and "+
     "1. Type 'rainbow' takes no additional arguments. Type 'gradients' "+
     "takes 2 or more sets of 3 additional arguments, specifying the "+
     "red, green, and blue components of the index colors to interpolate "+
     "between.",new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_density_colors),both},
    {0,"canvas","Create a new canvas.",0,0,"","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_canvas),both},
    {0,"plot","Plot two columns of data (2d only).",2,4,
     "<x> <y> [style] [color]","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_plot),both},
    {0,"hist-plot","Plot two columns of data (2d only).",2,5,
     "<n_bins> <x> [weights] [style] [color]","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_hist_plot),both},
    {0,"surf-plot","Surface plot of a slice (3d only).",2,2,"<slice> <option>",
     ((string)"Options are: scat (scatter plot), box (boxes with ")+
     "different sizes), arr (arrows pointing to higher values), text "+
     "(numerical weight at each point), contz (colored and filled "+
     "contours), cont1 (colored contour lines), cont2 (contour lines "+
     "with different styles), cont3 (black solid contour lines), lego "+
     "(3d bar plot), lego1 (3d bar plot w/ one edge filled), surf1 "+
     "(colored surface plot with grid, surf1z for color axis), surf2 "+
     "(colored surface plot, surf2z for color axis), surf3 (colored and "+
     "filled contours in 2D with grid surface plot), surf4 (monochromatic "+
     "surface with simple ray tracing), surf1pol (like 'surf1' but in "+
     "polar coords.), surf1cyl (like surf1 but in cylindrical coords.).",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_surf_plot),both},
    {0,"line","Draw a line.",4,6,"<x1> <y1> <x2> <y2> [style] [color]","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_line),both},
    {0,"arrow","Draw an arrow.",4,6,
     "<tail x> <tail y> <head x> <head y> [style] [color]","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_arrow),both},
    {0,"den-plot","Density plot (3d only).",1,1,"<slice>","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_den_plot),both},
    {0,"mden-plot","Multiple-density plot (3d only).",0,0,"","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_mden_plot),both},
    {0,"add-density","Add density to the list (3d only).",3,3,
     "<slice> <label> <color>","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_add_density),both},
    {0,"points","Plot two columns of data with points (2d only).",2,4,
     "<x> <y> [style] [color]","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_points),both},
    {0,"plotrun","Run the interactive ROOT interface.",0,0,"","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_plotrun),both},
    {0,"axis","Reformat the current plot axes.",4,4,
     "<left> <bottom> <right> <top>","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_axis),both},
    {0,"print","Write current plot to a file.",1,1,"<filename>","",
     new comm_option_mfptr<agraph_manager>
     (this,&agraph_manager::comm_print),both}
  };
  
  cl->set_comm_option_vec(narr,options_arr);
  
  return 0;
}

int agraph_manager::setup_help() {

  cl->cmd_name="agraph";
  
  cl->desc=(std::string)("agraph: A data table viewing and ")+
    "processing program for O2scl.\n";
  
  string dsc="\nNotes:\n\n";
  dsc+="1. Help for individual commands may be obtained with 'help ";
  dsc+="command'.\n   Required arguments are surrounded by <>'s and ";
  dsc+="optional arguments are\n   surrounded by []'s.\n\n";
  dsc+="2. Options may also be specified in the environment variable ";
  dsc+="AGRAPH_DEFAULTS.\n\n";
  dsc+="3. Long options may be preceeded by two dashes.\n\n";
  dsc+="4. In order to avoid argument confusion use \"(-x*2)\", ";
  dsc+="not \"-x*2\"\n\n";
  dsc+="Known operators:\n() - ^ * / % + - = < > & |\n\n";
  dsc+="Known functions:\n";
  dsc+="abs(x) acos(x) acosh(x) asin(x) asinh(x) atan(x) atan2(x,y)\n";
  dsc+="atanh(x) ceil(x) cos(x) cosh(x) cot(x) csc(x) eval(...) exp(x)\n";
  dsc+="floor(x) if(x,y,z) int(x) log(x) log10(x) max(x,y) min(x,y)\n";
  dsc+="sec(x) sin(x) sinh(x) sqrt(x) tan(x) tanh(x)\n\n";
  
  dsc+=((string)"Compiled at ")+((string)__TIME__)+" on "+
    ((string)__DATE__)+" for "+((string)PACKAGE)+", version "+
    ((string)VERSION)+".\n";
  
  cl->addl_help_cmd=dsc;
  cl->addl_help_cli=dsc;

  return 0;
}

int agraph_manager::setup_parameters() {

  acol_manager::setup_parameters();

  p_xlo.d=&xlo;
  p_xhi.d=&xhi;
  p_ylo.d=&ylo;
  p_yhi.d=&yhi;
  p_zlo.d=&zlo;
  p_zhi.d=&zhi;

  p_xlo.help="The lower limit for the x-axis";
  p_xhi.help="The upper limit for the x-axis";
  p_ylo.help="The lower limit for the y-axis";
  p_yhi.help="The upper limit for the y-axis";
  p_zlo.help="The lower limit for the z-axis";
  p_zhi.help="The upper limit for the z-axis";

  cl->par_list.insert(make_pair("xlo",&p_xlo));
  cl->par_list.insert(make_pair("xhi",&p_xhi));
  cl->par_list.insert(make_pair("ylo",&p_ylo));
  cl->par_list.insert(make_pair("yhi",&p_yhi));
  cl->par_list.insert(make_pair("zlo",&p_zlo));
  cl->par_list.insert(make_pair("zhi",&p_zhi));

  p_lmar.d=&plmar;
  p_tmar.d=&ptmar;
  p_rmar.d=&prmar;
  p_bmar.d=&pbmar;

  cl->par_list.insert(make_pair("lmar",&p_lmar));
  cl->par_list.insert(make_pair("tmar",&p_tmar));
  cl->par_list.insert(make_pair("rmar",&p_rmar));
  cl->par_list.insert(make_pair("bmar",&p_bmar));

  p_lmar.help="Left margin for plot (default 0.11)";
  p_tmar.help="Top margin for plot (default 0.05)";
  p_rmar.help="Right margin for plot (default 0.05)";
  p_bmar.help="Bottom margin for plot (default 0.1)";

  p_xtitle.str=&xtitle;
  p_ytitle.str=&ytitle;
  p_wtitle.str=&wtitle;
  p_wtitle2.str=&wtitle2;
  p_ptitle.str=&ptitle;

  p_xtitle.help="The x-axis title";
  p_ytitle.help="The y-axis title";
  p_wtitle.help="The window title";
  p_wtitle2.help="The other window title";
  p_ptitle.help="The pad title";

  cl->par_list.insert(make_pair("xtitle",&p_xtitle));
  cl->par_list.insert(make_pair("ytitle",&p_ytitle));
  cl->par_list.insert(make_pair("wtitle",&p_wtitle));
  cl->par_list.insert(make_pair("wtitle2",&p_wtitle2));
  cl->par_list.insert(make_pair("ptitle",&p_ptitle));

  p_logx.b=&logx;
  p_logy.b=&logy;
  p_logz.b=&logz;

  p_logx.help="If true, use a log scale for the x-axis";
  p_logy.help="If true, use a log scale for the y-axis";
  p_logz.help="If true, use a log scale for the z-axis";

  cl->par_list.insert(make_pair("logx",&p_logx));
  cl->par_list.insert(make_pair("logy",&p_logy));
  cl->par_list.insert(make_pair("logz",&p_logz));

  p_xset.b=&xset;
  p_yset.b=&yset;
  p_zset.b=&zset;

  p_xset.help="If true, xlo and xhi are x-axis limits";
  p_yset.help="If true, ylo and yhi are y-axis limits";
  p_zset.help="If true, zlo and zhi are z-axis limits";

  cl->par_list.insert(make_pair("xset",&p_xset));
  cl->par_list.insert(make_pair("yset",&p_yset));
  cl->par_list.insert(make_pair("zset",&p_zset));
  
  p_lwidth.i=&lwidth;
  p_lwidth.help="Line width (default 1).";
  cl->par_list.insert(make_pair("lwidth",&p_lwidth));

  p_multi_bins.i=&multi_bins;
  p_multi_bins.help="Number of bins for mden-plot (default 30).";
  cl->par_list.insert(make_pair("multi-bins",&p_multi_bins));

  p_talign.i=&talign;
  p_talign.help="Line width (default 22).";
  cl->par_list.insert(make_pair("talign",&p_talign));

  return 0;
}

int agraph_manager::comm_set(std::vector<std::string> &sv, bool itive_com) {
  
  if (sv.size()>2 && (sv[1]==((string)"xlo") || sv[1]==((string)"xhi"))) {
    xset=true;
  }
  if (sv.size()>2 && (sv[1]==((string)"ylo") || sv[1]==((string)"yhi"))) {
    yset=true;
  }
  if (sv.size()>2 && (sv[1]==((string)"zlo") || sv[1]==((string)"zhi"))) {
    zset=true;
  }

  return 0;
}

int agraph_manager::comm_density_colors(std::vector<std::string> &sv, 
					bool itive_com) {
  
  // This should be taken care of by the 'cli' class, but we
  // double-check just in case.
  if (sv.size()<3) {
    cout << "Not enough arguments in 'density-colors'." << endl;
    return exc_efailed;
  }

  size_t nc=o2scl::stoi(sv[2]), c_max=1000;
  if (nc>c_max) {
    if (verbose>0) {
      cout << "More than " << c_max 
	   << " colors specified. Setting to " << c_max << " colors." 
	   << endl;
    }
    nc=c_max;
  }

  if (nc<5) {
    cerr << "Number of colors below useful range (less than 5) "
	 << "in 'density-colors'." << endl;
    return exc_efailed;
  }

  if (sv[1]=="rgb") {

    if (sv.size()<6) {
      cout << "Not enough arguments in 'density-colors rgb'." << endl;
      return exc_efailed;
    }
    rcm.colors_rgb(nc,sv[3],sv[4],sv[5]);

  } else if (sv[1]=="hsv") {

    if (sv.size()<6) {
      cout << "Not enough arguments in 'density-colors hsv'." << endl;
      return exc_efailed;
    }
    rcm.colors_hsv(nc,sv[3],sv[4],sv[5]);

  } else if (sv[1]=="rainbow") {

    rcm.colors_rainbow(nc);

  } else if (sv[1]=="gradients") {

    bool debug=false;
    if (debug) cout << "size: " << sv.size() << endl;
    size_t nbase=sv.size()/3-1;
    if (nbase<2) {
      cout << "Less than two base gradient colors." << endl;
      return exc_efailed;
    }
    if ((nbase+1)*3!=sv.size()) {
      cout << "Didn't specify colors properly for gradients." << endl;
      return exc_efailed;
    }
    if (debug) cout << nbase << " base colors." << endl;

    // Create color vectors from the arguments
    std::vector<double> r, g, b;
    std::vector<int> index;
    for(size_t i=0;i<nbase;i++) {
      if (o2scl::stod(sv[3*i+3])>1.0 || o2scl::stod(sv[3*i+3])<0.0 ||
	  o2scl::stod(sv[3*i+4])>1.0 || o2scl::stod(sv[3*i+4])<0.0 ||
	  o2scl::stod(sv[3*i+5])>1.0 || o2scl::stod(sv[3*i+5])<0.0) {
	cout << "Colors out of range for gradients." << endl;
	return exc_efailed;
      }
      r.push_back(o2scl::stod(sv[3*i+3]));
      g.push_back(o2scl::stod(sv[3*i+4]));
      b.push_back(o2scl::stod(sv[3*i+5]));
      if (i==0) {
	index.push_back(0);
      } else if (i==nbase-1) {
	index.push_back(nc);
      } else {
	index.push_back(((int)(((double)i)/((double)(nbase-1))*
			       ((double)nc))));
      }
    }
    if (debug) {
      for(size_t i=0;i<nbase;i++) {
	cout << index[i] << " " << r[i] << " " << g[i] << " " << b[i] << endl;
      }
      cout << endl;
      cout << "Full list: " << endl;
    }
    rcm.colors_gradients(nc,nbase,r,g,b);

  }

  return 0;
}

int agraph_manager::comm_canvas(std::vector<std::string> &sv, bool itive_com) {

  if (theApp==0) {
    cout << "Creating new TApplication object." << endl;
    theApp=new TApplication("App",0,NULL);
    if (verbose>0 && itive_com) {
      cout << "Plot window will not close until you quit agraph." << endl;
    }
  }

  if (int_th1!=0) {
    cout << "Deleting old axis." << endl;
    delete int_th1;
    int_th1=0;
  }
  if (top_axis!=0) {
    cout << "Deleting top axis." << endl;
    delete top_axis;
    top_axis=0;
  }
  if (right_axis!=0) {
    cout << "Deleting right axis." << endl;
    delete right_axis;
    right_axis=0;
  }
  if (int_pad!=0) {
    cout << "Deleting pad." << endl;
    delete int_pad;
    int_pad=0;
  }
  if (int_canvas!=0) {
    cout << "Deleting canvas." << endl;
    delete int_canvas;
    int_canvas=0;
  }

  cout << "Creating canvas." << endl;
  string title_temp, title2_temp;
  if (wtitle.length()==0) {
    title_temp="window_title";
  } else {
    title_temp=wtitle;
  }
  if (wtitle2.length()==0) {
    title2_temp="window_title2";
  } else {
    title2_temp=wtitle2;
  }
  int_canvas=new TCanvas(title_temp.c_str(),title2_temp.c_str(),0,0,640,640);
  int_canvas->SetFillColor(10);

  cout << "Creating pad." << endl;
  if (ptitle.length()==0) {
    int_pad=new TPad("pad_title","",0.0,0.0,1.0,1.0);
  } else {
    int_pad=new TPad(ptitle.c_str(),"",0.0,0.0,1.0,1.0);
  }
  int_pad->SetFillColor(10);
  int_pad->SetTopMargin(ptmar);
  int_pad->SetRightMargin(prmar);
  int_pad->SetLeftMargin(plmar);
  int_pad->SetBottomMargin(pbmar);
  int_pad->Draw();
  int_pad->cd();

  if (logx) int_pad->SetLogx();
  if (logy) int_pad->SetLogy();

  line_ix=0;
  color_ix=0;

  return 0;
}

int agraph_manager::comm_axis(std::vector<std::string> &sv, 
			      bool itive_com) {
  if (sv.size()<5) {
    cerr << "Not enough arguments to 'axis' command." << endl;
    return exc_efailed;
  }

  internal_axis(o2scl::stod(sv[1]),o2scl::stod(sv[2]),
		o2scl::stod(sv[3]),o2scl::stod(sv[4]));
  return 0;
}

int agraph_manager::internal_axis(double left, double bottom, double right,
				  double top) {
  
  // If there's no canvas, create one
  if (int_canvas==0) {
    std::vector<std::string> sv2;
    bool itive_com=true;
    comm_canvas(sv2,itive_com);
  }

  if (int_th1!=0) {
    cout << "Deleting axis." << endl;
    delete int_th1;
  }
  if (top_axis!=0) {
    cout << "Deleting top axis." << endl;
    delete top_axis;
    top_axis=0;
  }
  if (right_axis!=0) {
    cout << "Deleting right axis." << endl;
    delete right_axis;
    right_axis=0;
  }

  a_left=left;
  a_bottom=bottom;
  a_right=right;
  a_top=top;

  // -----------------------------------------------------------------
  // Draw axes
  
  int_th1=int_pad->DrawFrame(a_left,a_bottom,a_right,a_top);
  int_th1->GetXaxis()->SetLabelFont(132);
  int_th1->GetYaxis()->SetLabelFont(132);
  int_th1->GetXaxis()->CenterTitle(kTRUE);
  int_th1->GetYaxis()->CenterTitle(kTRUE);
  int_th1->GetXaxis()->SetTitleFont(132);
  int_th1->GetYaxis()->SetTitleFont(132);
  int_th1->GetXaxis()->SetNdivisions(510);
  int_th1->GetYaxis()->SetNdivisions(510);

  // Axis labels
  o2scl_graph::axis_labels(left,bottom,right,top,talign,
			   xtitle,ytitle,logx,logy);
  
  if (true) {
    if (logx) {
      top_axis=new TGaxis(a_left,a_top,a_right,a_top,
			  a_left,a_right,510,"-G");
    } else {
      top_axis=new TGaxis(a_left,a_top,a_right,a_top,
			  a_left,a_right,510,"-");
    }
    top_axis->SetLabelFont(132);
    top_axis->SetLabelSize(0.0);
    top_axis->CenterTitle(kTRUE);
    top_axis->Draw();
    
    if (logy) {
      right_axis=new TGaxis(a_right,a_bottom,a_right,a_top,
			    a_bottom,a_top,510,"+G");
    } else {
      right_axis=new TGaxis(a_right,a_bottom,a_right,a_top,
			    a_bottom,a_top,510,"+");
    }
    right_axis->SetLabelFont(132);
    right_axis->SetLabelSize(0.0);
    right_axis->CenterTitle(kTRUE);
    right_axis->Draw();
  }

  int_canvas->Update();
  
  return 0;
}

int agraph_manager::comm_plot(std::vector<std::string> &sv, 
			      bool itive_com) {

  if (tabp==0) {
    cerr << "No table with data to plot." << endl;
    return exc_efailed;
  }

  if (sv.size()<3) {
    cerr << "Not enough arguments spectified in plot." << endl;
    return exc_efailed;
  }

  // Double check that the columns are there
  if (!tabp->is_column(sv[1])) {
    cerr << "Column '" << sv[1] << "' is not in table." << endl;
    return exc_einval;
  }
  if (!tabp->is_column(sv[2])) {
    cerr << "Column '" << sv[2] << "' is not in table." << endl;
    return exc_einval;
  }

  // sv[1] is the first column name
  // sv[2] is the second column name
  // sv[3] (optional) is the line style
  // sv[4] (optional) is the line color

  if (int_canvas==0) {
    std::vector<std::string> sv2;
    comm_canvas(sv2,itive_com);
  }

  if (int_th1==0) {

    // Initialize to avoid warnings about uninit'ed vars
    double left=0.0, right=0.0, top=0.0, bottom=0.0;
    
    if (xset) {
      
      left=xlo;
      right=xhi;

      if (yset) {

	/// Both X and Y ranges specified, so use them
	if (verbose>1) {
	  cout << "Both x and y ranges were specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {

	/// Only X-range is specified, automatically compute y-range
	if (verbose>1) {
	  cout << "The x range was specified. Computing y range." << endl;
	}

	bottom=tabp->min(sv[2]);
	top=tabp->max(sv[2]);

      }


    } else {

      // X-range was not specified, compute it.
      if (verbose>1) {
	cout << "The x range was not specified. Computing x range." << endl;
      }
    
      left=tabp->min(sv[1]);
      right=tabp->max(sv[1]);
      
      if (yset) {
      
	// Y-range was specified
	if (verbose>1) {
	  cout << "The y range was specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {
      
	// Y-range was not specified
	if (verbose>1) {
	  cout << "Computing y range." << endl;
	}

	bottom=tabp->min(sv[2]);
	top=tabp->max(sv[2]);

      }
    
    }

    if (verbose>1) {
      cout << "X range: " << left << " " << right << endl;
      cout << "Y range: " << bottom << " " << top << endl;
    }

    // -----------------------------------------------------------------
    // Draw axes

    internal_axis(left,bottom,right,top);
  }

  // -----------------------------------------------------------------
  // Plot data

  ubvector x1, x2;
  tabp->column_to_vector(sv[1],x1);
  tabp->column_to_vector(sv[2],x2);

  TGraph *gr=new TGraph(tabp->get_nlines());
  for(size_t j=0;j<tabp->get_nlines();j++) {
    gr->SetPoint(j,x1[j],x2[j]);
  }
  gr->SetName(sv[2].c_str());

  // Line style
  if (sv.size()>=4) {
    gr->SetLineStyle(o2scl::stoi(sv[3]));
  } else {
    gr->SetLineStyle(line_ix);
    line_ix++;
    if (line_ix==6) line_ix=0;
  }

  // Line color
  if (sv.size()>=5) {
    gr->SetLineColor(o2scl::stoi(sv[4]));
  } else {
    gr->SetLineColor(std_colors[color_ix]);
    color_ix++;
    if (color_ix==9) color_ix=0;
  }

  // Line width
  gr->SetLineWidth(lwidth);

  gr->Draw();

  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_hist_plot(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (tabp==0) {
    cerr << "No table with data to hist-plot." << endl;
    return exc_efailed;
  }

  if (sv.size()<3) {
    cerr << "Not enough arguments spectified in hist-plot." << endl;
    return exc_efailed;
  }

  // sv[1] is the number of histogram points
  // sv[2] is the column name
  // sv[3] (optional) is the column for the weights
  // sv[4] (optional) is the line style
  // sv[5] (optional) is the line color

  if (int_canvas==0) {
    std::vector<std::string> sv2;
    comm_canvas(sv2,itive_com);
  }

  int nbins=o2scl::stoi(sv[1]);
  if (nbins<=0) {
    cerr << "Number of bins less than 1 in hist-plot." << endl;
    return exc_efailed;
  }

  // Fill histogram
  hist h;
  uniform_grid<double> ug;
  if (logx) {
    ug=uniform_grid_log_end<double>(tabp->min(sv[2]),tabp->max(sv[2]),nbins);
  } else {
    ug=uniform_grid_end<double>(tabp->min(sv[2]),tabp->max(sv[2]),nbins);
  }
  h.set_bin_edges(ug);
  if (sv.size()==4) {
    for(size_t i=0;i<tabp->get_nlines();i++) {
      h.update(tabp->get(sv[2],i),tabp->get(sv[3],i));
    }
  } else {
    for(size_t i=0;i<tabp->get_nlines();i++) {
      h.update(tabp->get(sv[2],i));
    }
  }

  if (int_th1==0) {

    // Initialize to avoid warnings about uninit'ed vars
    double left=0.0, right=0.0, top=0.0, bottom=0.0;
    
    if (xset) {
      
      left=xlo;
      right=xhi;

      if (yset) {

	/// Both X and Y ranges specified, so use them
	if (verbose>1) {
	  cout << "Both x and y ranges were specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {

	/// Only X-range is specified, automatically compute y-range
	if (verbose>1) {
	  cout << "The x range was specified. Computing y range." << endl;
	}

	bottom=h.get_min_wgt();
	top=h.get_max_wgt();

      }


    } else {

      // X-range was not specified, compute it.
      if (verbose>1) {
	cout << "The x range was not specified. Computing x range." << endl;
      }
    
      left=tabp->min(sv[2]);
      right=tabp->max(sv[2]);
      
      if (yset) {
      
	// Y-range was specified
	if (verbose>1) {
	  cout << "The y range was specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {
      
	// Y-range was not specified
	if (verbose>1) {
	  cout << "Computing y range." << endl;
	}

	bottom=h.get_min_wgt();
	top=h.get_max_wgt();

      }
    
    }

    if (verbose>1) {
      cout << "X range: " << left << " " << right << endl;
      cout << "Y range: " << bottom << " " << top << endl;
    }

    // -----------------------------------------------------------------
    // Draw axes

    internal_axis(left,bottom,right,top);
  }

  // -----------------------------------------------------------------
  // Plot data

  TGraph *gr=new TGraph(h.size());
  for(size_t j=0;j<h.size();j++) {
    gr->SetPoint(j,h.get_rep_i(j),h.get_wgt_i(j));
  }
  gr->SetName(sv[2].c_str());

  // Line style
  if (sv.size()>=4) {
    gr->SetLineStyle(o2scl::stoi(sv[4]));
  } else {
    gr->SetLineStyle(line_ix);
    line_ix++;
    if (line_ix==6) line_ix=0;
  }

  // Line color
  if (sv.size()>=5) {
    gr->SetLineColor(o2scl::stoi(sv[5]));
  } else {
    gr->SetLineColor(std_colors[color_ix]);
    color_ix++;
    if (color_ix==9) color_ix=0;
  }

  // Line width
  gr->SetLineWidth(lwidth);

  gr->Draw();

  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_line(std::vector<std::string> &sv, 
			      bool itive_com) {

  // sv[1],sv[2] is the first coordinate
  // sv[3],sv[4] is the second coordinate
  // sv[5] (optional) is the line style
  // sv[6] (optional) is the line color

  if (int_canvas==0 || int_th1==0) {
    cerr << "Cannot draw line without a current plot and axes." << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------------
  // Plot data

  double x1[2]={o2scl::stod(sv[1]),o2scl::stod(sv[3])};
  double x2[2]={o2scl::stod(sv[2]),o2scl::stod(sv[4])};

  TGraph *gr=new TGraph(2,x1,x2);

  if (sv.size()>5) {
    gr->SetLineStyle(o2scl::stoi(sv[5]));
  } else {
    gr->SetLineStyle(line_ix);
    line_ix++;
    if (line_ix==6) line_ix=0;
  }
  if (sv.size()>6) {
    gr->SetLineColor(o2scl::stoi(sv[6]));
  } else {
    gr->SetLineColor(std_colors[color_ix]);
    color_ix++;
    if (color_ix==9) color_ix=0;
  }
  gr->SetLineWidth(lwidth);
  gr->Draw();
  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_arrow(std::vector<std::string> &sv, 
			      bool itive_com) {

  // sv[1],sv[2] is the first coordinate (tail)
  // sv[3],sv[4] is the second coordinate (head)
  // sv[5] (optional) is the line style
  // sv[6] (optional) is the line color

  if (int_canvas==0 || int_th1==0) {
    cerr << "Cannot draw arrow without a current plot and axes." << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------------
  // Plot data

  double x1[2]={o2scl::stod(sv[1]),o2scl::stod(sv[3])};
  double x2[2]={o2scl::stod(sv[2]),o2scl::stod(sv[4])};

  TLine *l1;
  TPolyLine *pl1;
  o2scl_graph::arrow(x1[0],x2[0],x1[1],x2[1],l1,pl1);

  if (sv.size()>5) {
    l1->SetLineStyle(o2scl::stoi(sv[5]));
  } else {
    l1->SetLineStyle(line_ix);
    line_ix++;
    if (line_ix==6) line_ix=0;
  }
  if (sv.size()>6) {
    l1->SetLineColor(o2scl::stoi(sv[6]));
    pl1->SetFillColor(o2scl::stoi(sv[6]));
  } else {
    l1->SetLineColor(std_colors[color_ix]);
    pl1->SetLineColor(std_colors[color_ix]);
    pl1->SetFillColor(std_colors[color_ix]);
    color_ix++;
    if (color_ix==9) color_ix=0;
  }
  l1->SetLineWidth(lwidth);
  l1->Draw();
  pl1->Draw();
  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_barplot(std::vector<std::string> &sv, 
				 bool itive_com) {
  
  // sv[1] is the first column name
  // sv[2] is the second column name
  // sv[3] is the third column name
  // sv[4] (optional) is the line and fill color
  // sv[5] (optional) is the fill style
  // sv[6] (optional) is the line style

  if (int_canvas==0) {
    comm_canvas(sv,itive_com);
  }

  if (int_th1==0) {

    // Initialize to avoid warnings about uninit'ed vars
    double left=0.0, right=0.0, top=0.0, bottom=0.0;
    
    if (xset) {
      
      left=xlo;
      right=xhi;

      if (yset) {

	/// Both X and Y ranges specified, so use them
	if (verbose>1) {
	  cout << "Both x and y ranges were specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {

	/// Only X-range is specified, automatically compute y-range
	if (verbose>1) {
	  cout << "The x range was specified. Computing y range." << endl;
	}

	bottom=tabp->min(sv[3]);
	top=tabp->max(sv[3]);

      }


    } else {

      // X-range was not specified, compute it.
      if (verbose>1) {
	cout << "The x range was not specified. Computing x range." << endl;
      }
    
      left=tabp->min(sv[1]);
      if (tabp->min(sv[2])<left) left=tabp->min(sv[2]);
      right=tabp->max(sv[1]);
      if (tabp->max(sv[2])<right) right=tabp->max(sv[2]);
      
      if (yset) {
      
	// Y-range was specified
	if (verbose>1) {
	  cout << "The y range was specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {
      
	// Y-range was not specified
	if (verbose>1) {
	  cout << "Computing y range." << endl;
	}

	bottom=tabp->min(sv[3]);
	top=tabp->max(sv[3]);

      }
    
    }

    if (verbose>1) {
      cout << "X range: " << left << " " << right << endl;
      cout << "Y range: " << bottom << " " << top << endl;
    }

    // -----------------------------------------------------------------
    // Draw axes

    internal_axis(left,bottom,right,top);
  }

  // -----------------------------------------------------------------
  // Plot data

  for(size_t k=1;k<tabp->get_nlines()-1;k++) {

    // First do fill
    TBox *b1=new TBox(tabp->get(sv[1],k),0.0,
		      tabp->get(sv[2],k),tabp->get(sv[3],k));
    if (sv.size()>=5) {
      b1->SetFillColor(o2scl::stoi(sv[4]));
    } else {
      b1->SetFillColor(std_colors[color_ix]);
    }
    if (sv.size()>=6) {
      b1->SetFillStyle(o2scl::stoi(sv[5]));
    } else {
      b1->SetFillStyle(1000);
    }
    b1->Draw();

    // Then do outer lines
    TBox *b2=new TBox(tabp->get(sv[1],k),0.0,
		      tabp->get(sv[2],k),tabp->get(sv[3],k));
    if (sv.size()>=5) {
      b2->SetLineColor(o2scl::stoi(sv[4]));
    } else {
      b2->SetLineColor(std_colors[color_ix]);
    }
    if (sv.size()>=7) {
      b2->SetLineStyle(o2scl::stoi(sv[6]));
    } else {
      b2->SetLineStyle(line_ix);
      line_ix++;
      if (line_ix==6) line_ix=0;
    }
    b2->Draw();

  }

  // Increment color index if necessary
  if (sv.size()<5) {
    color_ix++;
    if (color_ix==9) color_ix=0;
  }

  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_points(std::vector<std::string> &sv, 
				bool itive_com) {

  // sv[1] is the first column name
  // sv[2] is the second column name
  // sv[3] (optional) is the marker type
  // sv[4] (optional) is the marker color

  if (int_canvas==0) {
    comm_canvas(sv,itive_com);
  }

  if (int_th1==0) {

    // Initialize to avoid warnings about uninit'ed vars
    double left=0.0, right=0.0, top=0.0, bottom=0.0;
    
    if (xset) {
      
      left=xlo;
      right=xhi;

      if (yset) {

	/// Both X and Y ranges specified, so use them
	if (verbose>1) {
	  cout << "Both x and y ranges were specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {

	/// Only X-range is specified, automatically compute y-range
	if (verbose>1) {
	  cout << "The x range was specified. Computing y range." << endl;
	}

	bottom=tabp->min(sv[2]);
	top=tabp->max(sv[2]);

      }


    } else {

      // X-range was not specified, compute it.
      if (verbose>1) {
	cout << "The x range was not specified. Computing x range." << endl;
      }
    
      left=tabp->min(sv[1]);
      right=tabp->max(sv[1]);
      
      if (yset) {
      
	// Y-range was specified
	if (verbose>1) {
	  cout << "The y range was specified." << endl;
	}

	bottom=ylo;
	top=yhi;
      
      } else {
      
	// Y-range was not specified
	if (verbose>1) {
	  cout << "Computing y range." << endl;
	}

	bottom=tabp->min(sv[2]);
	top=tabp->max(sv[2]);

      }
    
    }

    if (verbose>1) {
      cout << "X range: " << left << " " << right << endl;
      cout << "Y range: " << bottom << " " << top << endl;
    }

    // -----------------------------------------------------------------
    // Draw axes

    internal_axis(left,bottom,right,top);
  }

  // -----------------------------------------------------------------
  // Plot data

  ubvector x1, x2;
  tabp->column_to_vector(sv[1],x1);
  tabp->column_to_vector(sv[2],x2);

  for(size_t jx=0;jx<tabp->get_nlines();jx++) {
    TMarker *m1;
    if (sv.size()>=4) {
      m1=new TMarker(x1[jx],x2[jx],o2scl::stoi(sv[3]));
    } else {
      m1=new TMarker(x1[jx],x2[jx],marker_ix);
    }
    if (sv.size()>=5) {
      m1->SetMarkerColor(o2scl::stoi(sv[4]));
    } else {
      m1->SetMarkerColor(std_colors[color_ix]);
    }
    m1->Draw();
  }
  int_canvas->Update();

  if (sv.size()<4) {
    marker_ix++;
    if (marker_ix==6) marker_ix=8;
    if (marker_ix==9) marker_ix=20;
    if (marker_ix==35) marker_ix=2;
  }
  if (sv.size()<5) {
    color_ix++;
    if (color_ix==9) color_ix=0;
  }

  return 0;
}

int agraph_manager::comm_plotrun(std::vector<std::string> &sv, 
				 bool itive_com) {
  theApp->Run(kTRUE);
  return 0;
}

int agraph_manager::comm_den_plot(std::vector<std::string> &sv, 
				  bool itive_com) {

  if (!threed || t3p==0) {
    cerr << "No table3d to insert columns into." << endl;
    return exc_efailed;
  }
    
  if (sv.size()<2) {
    cerr << "Not enough parameters in 'den-plot'." << endl;
    return exc_efailed;
  }
  
  size_t slix;
  if (!t3p->is_slice(sv[1],slix)) {
    cerr << "No slice named '" << sv[1] << "' in table3d object." << endl;
    return exc_efailed;
  }

  // -----------------------------------------------------------------
  // Create new canvas each time
    
  if (int_th1!=0) {
    cout << "Deleting axis." << endl;
    delete int_th1;
    int_th1=0;
  }
  if (top_axis!=0) {
    cout << "Deleting top axis." << endl;
    delete top_axis;
    top_axis=0;
  }
  if (right_axis!=0) {
    cout << "Deleting right axis." << endl;
    delete right_axis;
    right_axis=0;
  }
  if (int_pad!=0) {
    cout << "Deleting pad." << endl;
    delete int_pad;
    int_pad=0;
  }
  if (int_canvas!=0) {
    cout << "Deleting canvas." << endl;
    delete int_canvas;
    int_canvas=0;
  }

  // Add to the right margin for the density scale
  prmar+=0.14;
  comm_canvas(sv,itive_com);
  prmar-=0.14;
    
  // -----------------------------------------------------------------
  // Determine the plot limits, and store them in the variables
  // left, right, bottom, and top
    
  size_t nx, ny;
  t3p->get_size(nx,ny);
    
  // Initialize to avoid warnings about uninit'ed vars
  double left=0.0, right=0.0, top=0.0, bottom=0.0;
  
  if (xset) {
    
    left=xlo;
    right=xhi;

    if (yset) {

      /// Both X and Y ranges specified, so use them
      if (verbose>1) {
	cout << "Both x and y ranges were specified." << endl;
      }

      bottom=ylo;
      top=yhi;
      
    } else {

      /// Only X-range is specified, automatically compute y-range
      if (verbose>1) {
	cout << "The x range was specified. Computing y range." << endl;
      }

      top=t3p->get_grid_y(0);
      bottom=t3p->get_grid_y(0);
      for(size_t i=1;i<ny;i++) {
	if (bottom>t3p->get_grid_y(i)) bottom=t3p->get_grid_y(i);
	if (top<t3p->get_grid_y(i)) top=t3p->get_grid_y(i);
      }

    }


  } else {

    // X-range was not specified, compute it.
    if (verbose>1) {
      cout << "The x range was not specified. Computing x range." << endl;
    }
      
    right=t3p->get_grid_x(0);
    left=t3p->get_grid_x(0);
    for(size_t i=1;i<nx;i++) {
      if (left>t3p->get_grid_x(i)) left=t3p->get_grid_x(i);
      if (right<t3p->get_grid_x(i)) right=t3p->get_grid_x(i);
    }
    
    if (yset) {
      
      // Y-range was specified
      if (verbose>1) {
	cout << "The y range was specified." << endl;
      }

      bottom=ylo;
      top=yhi;
      
    } else {
      
      // Y-range was not specified
      if (verbose>1) {
	cout << "Computing y range." << endl;
      }

      top=t3p->get_grid_y(0);
      bottom=t3p->get_grid_y(0);
      for(size_t i=1;i<ny;i++) {
	if (bottom>t3p->get_grid_y(i)) bottom=t3p->get_grid_y(i);
	if (top<t3p->get_grid_y(i)) top=t3p->get_grid_y(i);
      }

    }
    
  }

  internal_axis(left,bottom,right,top);

  int_pad->cd();
  
  tmdp.logx=logx;
  tmdp.logy=logy;
  tmdp.logz=logz;
  
  tmdp.plot(int_pad,*t3p,sv[1],rcm);
  
  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_add_density(std::vector<std::string> &sv, 
				     bool itive_com) {
  if (sv.size()<4) {
    cerr << "Not enough parameters in 'add-density'." << endl;
    return exc_efailed;
  }
  if (t3p==0) {
    cerr << "No table to add density from." << endl;
    return exc_efailed;
  }

  string slice, label, color;
  slice=sv[1];
  label=sv[2];
  color=sv[3];
  tmdp.add(*t3p,slice,label,color);

  return 0;
}

int agraph_manager::comm_mden_plot(std::vector<std::string> &sv, 
				   bool itive_com) {

  bool debug=false;

  // -----------------------------------------------------------------
  // Create new canvas each time
    
  if (int_th1!=0) {
    cout << "Deleting axis." << endl;
    delete int_th1;
    int_th1=0;
  }
  if (top_axis!=0) {
    cout << "Deleting top axis." << endl;
    delete top_axis;
    top_axis=0;
  }
  if (right_axis!=0) {
    cout << "Deleting right axis." << endl;
    delete right_axis;
    right_axis=0;
  }
  if (int_pad!=0) {
    cout << "Deleting pad." << endl;
    delete int_pad;
    int_pad=0;
  }
  if (int_canvas!=0) {
    cout << "Deleting canvas." << endl;
    delete int_canvas;
    int_canvas=0;
  }

  // Add to the right margin for the density scale
  prmar+=0.30;
  comm_canvas(sv,itive_com);
  prmar-=0.30;

  // -----------------------------------------------------------------
  // Get a pointer to the first table3d slice
  
  table3d *t3ploc=&(tmdp.get_table(0));
  if (debug) {
    cout << "pointer: " << t3ploc->get_nx() << " " << t3ploc->get_ny()
	 << endl;
  }
  
  // -----------------------------------------------------------------
  // Determine the plot limits, and store them in the variables
  // left, right, bottom, and top
    
  size_t nx, ny;
  t3ploc->get_size(nx,ny);
    
  // Initialize to avoid warnings about uninit'ed vars
  double left=0.0, right=0.0, top=0.0, bottom=0.0;
  
  if (xset) {
    
    left=xlo;
    right=xhi;

    if (yset) {

      /// Both X and Y ranges specified, so use them
      if (verbose>1) {
	cout << "Both x and y ranges were specified." << endl;
      }

      bottom=ylo;
      top=yhi;
      
    } else {

      /// Only X-range is specified, automatically compute y-range
      if (verbose>1) {
	cout << "The x range was specified. Computing y range." << endl;
      }

      top=t3ploc->get_grid_y(0);
      bottom=t3ploc->get_grid_y(0);
      for(size_t i=1;i<ny;i++) {
	if (bottom>t3ploc->get_grid_y(i)) bottom=t3ploc->get_grid_y(i);
	if (top<t3ploc->get_grid_y(i)) top=t3ploc->get_grid_y(i);
      }

    }


  } else {

    // X-range was not specified, compute it.
    if (verbose>1) {
      cout << "The x range was not specified. Computing x range." << endl;
    }
      
    right=t3ploc->get_grid_x(0);
    left=t3ploc->get_grid_x(0);
    for(size_t i=1;i<nx;i++) {
      if (left>t3ploc->get_grid_x(i)) left=t3ploc->get_grid_x(i);
      if (right<t3ploc->get_grid_x(i)) right=t3ploc->get_grid_x(i);
    }
    
    if (yset) {
      
      // Y-range was specified
      if (verbose>1) {
	cout << "The y range was specified." << endl;
      }

      bottom=ylo;
      top=yhi;
      
    } else {
      
      // Y-range was not specified
      if (verbose>1) {
	cout << "Computing y range." << endl;
      }

      top=t3ploc->get_grid_y(0);
      bottom=t3ploc->get_grid_y(0);
      for(size_t i=1;i<ny;i++) {
	if (bottom>t3ploc->get_grid_y(i)) bottom=t3ploc->get_grid_y(i);
	if (top<t3ploc->get_grid_y(i)) top=t3ploc->get_grid_y(i);
      }

    }
    
  }

  if (debug) cout << "Going to int axis: " << endl;
  internal_axis(left,bottom,right,top);

  int_pad->cd();

  if (multi_bins>0) {
    tmdp.n_bins=((size_t)multi_bins);
  } else {
    tmdp.n_bins=30;
  }
  tmdp.logx=logx;
  tmdp.logy=logy;
  tmdp.logz=logz;
  
  if (debug) cout << "Going to multi_plot: " << endl;
  tmdp.multi_plot(rcm);
  
  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_surf_plot(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (t3p==0 || !threed) {
    cerr << "No table3d to plot in 'surf-plot'." << endl;
    return exc_efailed;
  }
    
  if (sv.size()<3) {
    cout << "Not enough parameters in 'surf-plot'." << endl;
    return exc_efailed;
  }
    
  size_t slix=t3p->lookup_slice(sv[1]);
  
  size_t nx, ny;
  t3p->get_size(nx,ny);
    
  // -----------------------------------------------------------------
  // Create canvas if necessary
    
  if (int_canvas==0) {
    comm_canvas(sv,itive_com);
  }
    
  // -----------------------------------------------------------------
  // Create ROOT 2-D histogram from table3d object

  if (surf_hist!=0) delete surf_hist;
  ubvector xbins(nx+1), ybins(ny+1);
  const ubvector &tx=t3p->get_x_data();
  const ubvector &ty=t3p->get_y_data();
  for(size_t i=0;i<nx+1;i++) {
    if (i==0) {
      xbins[i]=tx[0]-(tx[1]-tx[0])/2.0;
    } else if (i==nx) {
      xbins[i]=tx[nx-1]+(tx[nx-1]-tx[nx-2])/2.0;
    } else {
      xbins[i]=(tx[i]+tx[i-1])/2.0;
    }
  }
  for(size_t i=0;i<ny+1;i++) {
    if (i==0) {
      ybins[i]=ty[0]-(ty[1]-ty[0])/2.0;
    } else if (i==ny) {
      ybins[i]=ty[ny-1]+(ty[ny-1]-ty[ny-2])/2.0;
    } else {
      ybins[i]=(ty[i]+ty[i-1])/2.0;
    }
  }
  gStyle->SetOptStat("");
  gStyle->SetPalette(1);
  surf_hist=new TH2D("h2","",nx,&xbins[0],ny,&ybins[0]);
  
  for(size_t i=0;i<nx;i++) {
    for(size_t j=0;j<ny;j++) {
      surf_hist->Fill(t3p->get_grid_x(i),t3p->get_grid_y(j),
		      t3p->get(i,j,sv[1]));
    }
  }

  // -----------------------------------------------------------------
  // Plot 

  surf_hist->Draw(sv[2].c_str());

  // -----------------------------------------------------------------

  int_canvas->Update();

  return 0;
}

int agraph_manager::comm_print(std::vector<std::string> &sv, 
			       bool itive_com) {
  int_canvas->Update();
  int_canvas->Print(sv[1].c_str());
  return 0;
}

int agraph_manager::comm_text(std::vector<std::string> &sv, 
			      bool itive_com) {

  if (int_canvas==0 || int_pad==0 || int_th1==0) {
    cerr << "No plot to add text to." << endl;
    return exc_efailed;
  }
  tt.SetTextAlign(talign);
  if (sv.size()==2) {
    if (logx) {
      if (logy) {
	tt.DrawLatex(sqrt(a_left*a_right),sqrt(a_top*a_bottom),sv[1].c_str());
      } else {
	tt.DrawLatex(sqrt(a_left*a_right),(a_top+a_bottom)/2.0,sv[1].c_str());
      }
    } else if (logy) {
      tt.DrawLatex((a_left+a_right)/2.0,sqrt(a_top*a_bottom),sv[1].c_str());
    } else {
      tt.DrawLatex((a_left+a_right)/2.0,(a_top+a_bottom)/2.0,sv[1].c_str());
    }
  } else {
    if (sv.size()>=5) {
      double fact=o2scl::stod(sv[4]);
      tt.SetTextSize(tt.GetTextSize()*fact);
    }
    tt.DrawLatex(o2scl::stod(sv[2]),o2scl::stod(sv[3]),sv[1].c_str());
    if (sv.size()>=5) {
      double fact=o2scl::stod(sv[4]);
      tt.SetTextSize(tt.GetTextSize()/fact);
    }
  }
  int_canvas->Update();

  return 0;
}

