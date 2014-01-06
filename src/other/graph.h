/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_GRAPH_H
#define O2SCL_GRAPH_H

// Standard headers
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// ROOT headers
#include <Gtypes.h>
#include <TApplication.h>
#include <TArrow.h>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMarker.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPolyLine.h>

// O2scl headers
#include <o2scl/shared_ptr.h>
#include <o2scl/constants.h>
#include <o2scl/table_units.h>
#include <o2scl/table3d.h>

/** \file graph.h
    \brief Experimental functions for use with Root

    All of the functions in \ref graph.h are documented in
    the \ref o2scl_graph namespace.
*/

/** \brief Namespace for experimental functions for use with Root

    This namespace contains several functions, defined in graph.cpp,
    which simplify plotting with Root. They require that the Root
    libraries have been installed and are not compiled in with the
    normal \o2 libraries. Some of the example plots created for 
    the documentation are based on these functions.
*/
namespace o2scl_graph {
  
  /** \brief Draw a pretty arrow from (x1,y1) to (x2,y2)
      
      The parameter \c size determines the size of the head relative
      to the size of the arrow, \c size2 determines the position of
      the niche (1.0 is no niche), and alpha1 determines the angle of
      the head relative to the line (default is about 20 degrees).
      
      The objects \c line and \c poly are created with \c new
      and the user is reponsible for calling \c delete as needed.
  */
  void arrow(double x1, double y1, double x2, double y2,
	     TLine *&line, TPolyLine *&poly, 
	     double size=0.1, double size2=0.8, 
	     double alpha1=0.35);
  
  /** \brief Create x- and y-labels for previously plotted axes
   */
  void axis_labels(double left, double bottom, double right, double top,
		   int talign, std::string xtitle, std::string ytitle, 
		   bool logx=false, bool logy=false, 
		   double xfactor=12.5, double yfactor=9.5);

  /** \brief Create a new canvas and pad for a generic plot
      
      The CanvasName is the default output filename for graphs and
      macros also. The PadName is the name that comes up when you
      right-click on the pad. When making .eps files, the Title field
      is set to "CanvasName.eps: Windowname". This Title field is the
      title that appears, e.g., in a ghostview window. The x-axis
      limits are given by \c lleft and \c lright and the y-axis limits
      are given in \c lbottom and \c ltop. The upper left corner of
      the canvas window is at <tt>(left,top)</tt> and the lower right
      corner is <tt>(right,bottom)</tt>. The parameters \c logx and \c
      logy determine whether or not the x- or y-axes are logarithmic
      instead of linear (the default). 

      The canvas and pad objects are created using \c new (the user
      must call delete as necessary), and the histogram object \c th1
      is created with TPad::DrawFrame().
  */
  void new_graph(TCanvas *&c1, TPad *&p1, TH1 *&th1, std::string CanvasName, 
		 std::string WindowName, std::string PadName, 
		 double lleft, double lbottom, double lright, double ltop,
		 int left=0, int top=0, int right=700, int bottom=700, 
		 bool logx=false, bool logy=false);
  
  /** \brief Create a new canvas and pad for a generic plot with
      tick marks on the upper and right-hand-side edges
      
      This is the same as \ref new_graph(), except it adds tick marks
      to the right and top axes, which are created with \c new and
      returned in \c a1 and \c a2 respectively().
  */
  void new_graph_ticks(TCanvas *&c1, TPad *&p1, TH1 *&th1, 
		       std::string CanvasName, std::string WindowName, 
		       std::string PadName, 
		       double lleft, double lbottom, double lright, 
		       double ltop, TGaxis *&a1, TGaxis *&a2,
		       int left=0, int top=0, int right=700, 
		       int bottom=700, bool logx=false, bool logy=false);
  
  /** \brief Graph two columns from a data table

      This function plots the function defined by
      <tt>(scolx,scoly)</tt> given in table \c at with line style \c
      style and line color \c color. A \c TGraph object created with 
      \c new is returned. The \c name of the graph object
      is automatically set to be the same as the \c scoly, but this
      can always be changed afterwards, i.e.
      \code
      g1=table_graph(at,scolx,scoly);
      g1->SetName("Curve name");
      g1->Draw();
      \endcode
  */
  TGraph *table_graph(const o2scl::table_units<> &at, std::string scolx, 
		      std::string scoly, int style=1, int color=1);

  /** \brief Create a \c TGraph object from two vectors of equal length
   */
  template<class vec_t, class vec2_t> 
    TGraph *vector_graph(size_t n, const vec_t &x, const vec2_t &y) {
    TGraph *gr1=new TGraph(n);
    for(size_t i=0;i<n;i++) {
      gr1->SetPoint(i,x[i],y[i]);
    }
    gr1->Draw();
    return gr1;
  }

  /** \brief Create a \c TGraph object from two vectors of equal length
   */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t> 
    TPolyLine *vector_graph_region
    (size_t n, const vec_t &x, const vec2_t &x2, const vec3_t &y, 
     const vec4_t &y2, std::string option="") {
     
    std::vector<double> ox, oy;
    for(size_t i=0;i<n;i++) {
      ox.push_back(x[i]);
      oy.push_back(y[i]);
    }
    for(size_t i=0;i<n;i++) {
      ox.push_back(x2[n-1-i]);
      oy.push_back(y2[n-1-i]);
    }
    ox.push_back(x[0]);
    oy.push_back(y[0]);

    TPolyLine *pl1=new TPolyLine(2*n+1);
    for(size_t i=0;i<2*n+1;i++) {
      pl1->SetPoint(i,ox[i],oy[i]);
    }
    pl1->Draw(option.c_str());

    return pl1;
  }
  
  /** \brief Plot colums from a data table with error bars
   */
  TGraphErrors *table_graph_errors
    (const o2scl::table_units<> &at, std::string scolx, std::string scoly, 
     std::string xerr, std::string yerr, int style=1, int color=1);
  
  /** \brief Make a canvas with a two-up graph, side-by-side

      This function creates a canvas and two pads for a side-by-side
      plot with two different x-axes and the same y-axis. The x-axis
      of the left-hand plot has limits <tt>(lowx1,highx1)</tt> and
      the the x-axis for the right-hand plot has limits
      <tt>(lowx2,highx2)</tt>.
  */
  void two_up_graph(TCanvas *&c1, TPad *&p1, TPad *&p2, TH1 *&th1, TH1 *&th2,
		    std::string CanvasName, std::string WindowName, 
		    std::string Pad1Name, std::string Pad2Name,
		    double lowx1, double highx1, double lowx2, double highx2,
		    double lowy, double highy,
		    int left=0, int top=0, int right=1000, int bottom=700,
		    bool logx1=false, bool logx2=false, bool logy=false,
		    double alpha=1.3, double margin=0.1);

  /** \brief Make a canvas with two plots, one on top of the other
      
      The variable \c p1 is on bottom and \c p2 is on top
  */
  void two_up_graphy(TCanvas *&c1, TPad *&p1, TPad *&p2, 
		     TH1 *&th1, TH1 *&th2,
		     std::string CanvasName, std::string WindowName, 
		     std::string Pad1Name, std::string Pad2Name,
		     double lowx, double highx,
		     double lowy1, double highy1, 
		     double lowy2, double highy2,
		     int left=0, int top=0, int right=1000, int bottom=700,
		     bool logx=false, bool logy1=false, bool logy2=false,
		     double alpha=1.3, double margin=0.1);

  /** \brief ROOT color manager

      \warning It's not clear if multiple calls to \ref allocate()
      work gracefully or not.
  */
  class root_color_manager {
    
  protected:

    /// Current number of colors 
    size_t n_colors;

    /// Number of colors allocated
    size_t n_allocated;

    /// Index of first color
    int col_start;

    /// Minimum number of colors to allocate (default 500)
    size_t min_allocate;

  public:

    root_color_manager();
  
    /// Allocate colors
    void allocate(size_t na);

    /// Set minimum allocation
    void set_min_allocate(size_t n);
    
    /// Return first color index
    int get_col_start() const {
      return col_start;
    }

    /// Return number of colors
    size_t get_n_colors() const {
      return n_colors;
    }

    /// Return number of colors
    size_t get_min_allocate() const {
      return min_allocate;
    }

    /** \brief Create a palette based on a RGB specification

	The functions <tt>red_func, green_func</tt>, and
	<tt>blue_func</tt> should take an index number \f$ \in [0,1]
	\f$ to color components \f$ \in [0,1] \f$.
     */
    void colors_rgb
      (size_t n, std::string red_func, std::string green_func, 
       std::string blue_func);

    /** \brief Create a palette based on a HSV specification

	The functions <tt>hue_func, sat_func</tt>, and
	<tt>val_func</tt> should take as input an index number \f$ \in
	[0,1] \f$ . The output for <tt>hue_func</tt> should be \f$ \in
	[0,360] \f$ , and the outputs for <tt>sat_func</tt> and
	<tt>val_func</tt> should be \f$ \in [0,1] \f$ .
     */
    void colors_hsv
      (size_t n, std::string hue_func, std::string sat_func, 
       std::string val_func);

    /** \brief Create a simple rainbow palette

	This creates a palette corresponding to
	\f$ (0-300,1,1) \f$ in HSV color space \f$ .
     */
    void colors_rainbow(size_t n);
  
    /** \brief Rescale a set of vectors defining colors to ensure that
	\f$ r,g,b \in [0,1] \f$ 
    */
    template<class vec_t>
      void rescale_colors(size_t n, vec_t &red, vec_t &green, vec_t &blue) {

      double r_max=vector_max(n,red);
      double g_max=vector_max(n,green);
      double b_max=vector_max(n,blue);
      double r_min=vector_min(n,red);
      double g_min=vector_min(n,green);
      double b_min=vector_min(n,blue);

      double r_scale, g_scale, b_scale;
      double r_shift, g_shift, b_shift;
      
      r_scale=fabs(r_max-r_min);
      if (r_max>r_min) r_shift=-r_min/r_scale;
      else r_shift=-r_max/r_scale;
      g_scale=fabs(g_max-g_min);
      if (g_max>g_min) g_shift=-g_min/g_scale;
      else g_shift=-g_max/g_scale;
      b_scale=fabs(b_max-b_min);
      if (b_max>b_min) b_shift=-b_min/b_scale;
      else b_shift=-b_max/b_scale;

      for(size_t i=0;i<n;i++) {
	red[i]=red[i]/r_scale+r_shift;
	green[i]=green[i]/g_scale+g_shift;
	blue[i]=blue[i]/b_scale+b_shift;
      }

      return;
    }

    /** \brief Set colors based on a specification in vectors

	\note Any values less than zero are assumed to be zero and any
	values greater than one are assumed to be one.
     */
    template<class vec_t>
      void set_colors(size_t n, const vec_t &red, 
		      const vec_t &green, const vec_t &blue) {
      allocate(n);
      for (size_t i=0;i<n;i++) {
	TColor *c=gROOT->GetColor(col_start+i);
	double rtmp=red[i];
	if (rtmp<0.0) rtmp=0.0;
	if (rtmp>1.0) rtmp=1.0;
	double gtmp=green[i];
	if (gtmp<0.0) gtmp=0.0;
	if (gtmp>1.0) gtmp=1.0;
	double btmp=blue[i];
	if (btmp<0.0) btmp=0.0;
	if (btmp>1.0) btmp=1.0;
	c->SetRGB(rtmp,gtmp,btmp);
      }
      n_colors=n;
      return;
    }

    /** \brief Create a color palette based on a gradients between a
	set of index colors

	Create a palette based on gradients linearly interpolated in
	RGB color space among a list of user-specified colors. The
	user specified \c nbase colors which forms the index points in
	a color palette of site <tt>n >= nbase</tt>. The \c nbase
	index colors are specified in \c red, \c green, and \c blue,
	all of which should be vectors with at least \c nbase
	elements.

	This function is similar in concept to
	<tt>TColor::CreateGradientColorTable()</tt>, but for generic
	vector types.
    */
    template<class vec_t>
      void colors_gradients(size_t n, size_t nbase, const vec_t &red, 
			    const vec_t &green, const vec_t &blue) {

      if (n<nbase) {
	O2SCL_ERR2("Input 'n' must be >= 'nbase' in root_color_manager",
		   "::colors_gradients().",o2scl::exc_einval);
      }
      
      bool debug=false;

      allocate(n);

      std::vector<int> index;
      std::vector<double> red_full, green_full, blue_full;
      for(size_t i=0;i<nbase;i++) {
	if (i==0) {
	  index.push_back(0);
	} else if (i==nbase-1) {
	  index.push_back(n);
	} else {
	  index.push_back(((int)(((double)i)/((double)(nbase-1))*
				 ((double)n))));
	}
      }
      if (debug) {
	for(size_t i=0;i<nbase;i++) {
	  std::cout << index[i] << " " << red[i] << " " << green[i] << " " 
		    << blue[i] << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Full list: " << std::endl;
      }
      for(size_t i=0;i<nbase-1;i++) {
	for(int j=index[i];j<index[i+1];j++) {
	  if (j==index[i]) {
	    // First color in every section is exactly equal to
	    // user-specified color
	    red_full.push_back(red[i]);
	    green_full.push_back(green[i]);
	    blue_full.push_back(blue[i]);
	  } else if (i==nbase-2 && j==index[i+1]-1) {
	    // Last color in last section is exactly equal to last
	    // user-specified color
	    red_full.push_back(red[nbase-1]);
	    green_full.push_back(green[nbase-1]);
	    blue_full.push_back(blue[nbase-1]);
	  } else {
	    // In between, we linearly interpolate
	    double rt=red[i]+((double)(j-index[i]))/
	      ((double)(index[i+1]-index[i]))*(red[i+1]-red[i]);
	    double gt=green[i]+((double)(j-index[i]))/
	      ((double)(index[i+1]-index[i]))*(green[i+1]-green[i]);
	    double bt=blue[i]+((double)(j-index[i]))/
	      ((double)(index[i+1]-index[i]))*(blue[i+1]-blue[i]);
	    red_full.push_back(rt);
	    green_full.push_back(gt);
	    blue_full.push_back(bt);
	  }
	  if (debug) {
	    std::cout << i << " " << j << " " << red_full.size()-1 << " "
		      << red_full[red_full.size()-1] << " "
		      << green_full[green_full.size()-1] << " "
		      << blue_full[blue_full.size()-1] << std::endl;
	  }
	}
      }
      
      // Redefine that palette with the user-specified colors
      for (size_t i=0;i<n;i++) {
	TColor *c=gROOT->GetColor(col_start+i);
	c->SetRGB(red_full[i],green_full[i],blue_full[i]);
      }
      
      n_colors=n;
      
      return;
    }
  
  };

  /** \brief Colors for ROOT based on the W3C list

      \todo This class isn't quite finished. I need to finish the
      coordination with root_color_manager.
   */
  class html_colors {

  protected:

    /// The index of the first color
    int col_start;

    /** \brief Color name and RGB values 
	[protected subclass of \ref html_colors ]
    */
    struct color_s {
      std::string name;
      size_t index;
      short int r;
      short int g;
      short int b;
    };

    /** \brief Storage for sorted list
     */
    std::map<std::string,struct color_s,o2scl::string_comp> cmap;
  
  public:

    html_colors();

    /// Get the color index [0,146]
    int get_color_index(std::string s) const;

    /// Get color values 
    void get_color_rgb(std::string s, double &r, double &g, double &b) const;

    /** \brief Add the color list to a \ref root_color_manager object
     */
    void add_colors(root_color_manager &rcm);

    /** \brief Get the index in the ROOT palette after a call to 
	\ref add_colors()
     */
    int operator[](std::string s) const;

  };

  /** \brief Create a density plot
      
      \todo Values left, right, top and bottom are computed twice,
      once to make the canvas and once in the plot. I'm not sure
      that's necessary. It might allow for different sizes between the
      table and the pad, but it might be complicated to keep track of
      both the table and the pad sizes.

      \todo Consider separating the density plot and the z-axis into
      two separate functions?
  */
  class table3d_density_plot {

  protected:
    
    /// \name x, y, and z limits
    //@{
    bool xset;
    double xleft;
    double xright;
    bool yset;
    double ybottom;
    double ytop;
    bool zset;
    double zbottom;
    double ztop;
    //@}

    /// \name Axis objects
    //@{
    TGaxis *aright, *aleft, *atop, *abottom, *ascale;
    //@}

    /// Box for density scale
    TBox *box1;

    /// Lines for density scale 
    //@{
    TLine *l1, *l2, *l3;
    //@}

    /// Contour lines
    std::vector<TGraph *> gt;
    
  public:
    
    table3d_density_plot();

    virtual ~table3d_density_plot() {}

    /// Right margin for pad (default 0.14)
    double prmar;
    
    /// Logarithmic x-axis (default false)
    bool logx; 
    
    /// Logarithmic y-axis (default false)
    bool logy;

    /// Logarithmic z-axis (default false)
    bool logz;

    /// \name Other properties used by plot_canvas()
    //@{
    /// Canvas name (default "Canvas 1")
    std::string canvas_name; 

    /// Window name (default "Window 1")
    std::string window_name; 
    
    /// Pad name (default "Pad 1")
    std::string pad_name; 

    /// Pad name (default "")
    std::string xtitle; 

    /// Pad name (default "")
    std::string ytitle;
    //@}

    /// \name Window coordinates (default values are 0, 0, 680, 640)
    //@{
    int wleft;
    int wtop;
    int wright;
    int wbottom;
    //@}

    /// Contour lines
    std::vector<o2scl::contour_line> conts;

    /// Set x limits
    void set_x(double xl, double xr) {
      xleft=xl;
      xright=xr;
      xset=true;
      return;
    }

    /// Set y limits
    void set_y(double yb, double yt) {
      ybottom=yb;
      ytop=yt;
      yset=true;
      return;
    }

    /// Set z limits
    void set_z(double zb, double zt) {
      zbottom=zb;
      ztop=zt;
      zset=true;
      return;
    }

    /** \brief Create canvas and pad
     */
    void plot_canvas(TCanvas *&c1, TPad *&p1, TH1 *&th1, o2scl::table3d &t);

    /** \brief Create a density plot on pad \c pad based on 
	slice \c slice from table \c t
    */
    void plot(TPad *pad, o2scl::table3d &t, std::string slice, 
	      root_color_manager &rcm);

    /// Free internal ROOT objects
    virtual void free() {
      delete aright, aleft, atop, abottom, box1, ascale, l1, l2, l3;
      for(size_t i=0;i<gt.size();i++) {
	delete gt[i];
      }
      gt.clear();
      return;
    }
    
  };
  
  /** \brief Create a multiple density plot

      \todo Current code doesn't correctly normalize the
      table3d density and doesn't do anything with the
      <tt>origin</tt> parameter. 
   */
  class table3d_multi_density : public table3d_density_plot {
    
  protected:

    /// Tables to plot
    std::vector<o2scl::table3d> tables;

    /// Slice names
    std::vector<std::string> slices;

    /// Labels
    std::vector<std::string> labels;
    
    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Patterns
    ubvector patterns;

    /// Starting colors
    std::vector<double> colors;

    /// To associate colors with numbers
    html_colors hcol;

  public:

    table3d_multi_density();
    
    virtual ~table3d_multi_density() {}

    /// Number of bins in the x and y directions
    size_t n_bins;

    /// \name Method of combining colors
    //@{
    size_t combine;
    static const size_t add_colors=0;
    static const size_t overplot=1;
    //@}
    
    /** \brief Add a table to the list to plot
     */
    virtual void add(o2scl::table3d &table, std::string slice, 
		     std::string label, std::string color) {
      tables.push_back(table);
      slices.push_back(slice);
      labels.push_back(label);
      double r,g,b;
      hcol.get_color_rgb(color,r,g,b);
      colors.push_back(r);
      colors.push_back(g);
      colors.push_back(b);
      return;
    }

    /** \brief Clear the list
     */
    virtual void clear_list() {
      tables.clear();
      slices.clear();
      labels.clear();
      colors.clear();
      return;
    }

    /** \brief Set patterns to use
     */
    template<class vec_t> void set_patterns(size_t n, const vec_t &pat) {
      patterns.resize(n);
      for(size_t i=0;i<n;i++) {
	patterns[i]=pat[i];
      }
      return;
    }

    /** \brief Get a table from the list by index

	This is required to be able to plot the canvas using the
	limits from the first table in the list.
    */
    o2scl::table3d &get_table(size_t ix) {
      if (ix>=tables.size()) {
	O2SCL_ERR("Not enough tables for get_table().",o2scl::exc_efailed);
      }
      return tables[ix];
    }
    
    /** \brief Create a multi-density plot
     */
    virtual void multi_plot(root_color_manager &rcm);
      
  };
  
  /// \name Some useful colors
  //@{
  const int kGray20=kGray;
  const int kGray40=kGray+1;
  const int kGray60=kGray+2;
  const int kGray80=kGray+3;
  //@}

  /// \name Markers
  //@{
  const int m_small_dot=1;
  const int m_plus=2;
  const int m_asterisk=3;
  const int m_circle=4;
  const int m_times=5;
  const int m_med_dot=6;
  const int m_large_dot=7;
  const int m_fill_circle=8;
  const int m_fill_square=21;
  const int m_fill_up_triangle=22;
  const int m_fill_dn_triangle=23;
  const int m_open_circle=24;
  const int m_open_square=25;
  const int m_open_up_triangle=26;
  const int m_open_diamond=27;
  const int m_open_plus=28;
  const int m_fill_star=29;
  const int m_open_star=30;
  const int m_asterisk2=31;
  const int m_open_dn_triangle=32;
  const int m_fill_diamond=33;
  const int m_fill_plus=34;
  //@}

}

#endif

