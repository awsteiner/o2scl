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
#include "../../config.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/acolm.h>
#include <o2scl/graph.h>

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
#include <TLine.h>
#include <TList.h>
#include <TMarker.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TPolyLine.h>

//#ifdef O2SCL_READLINE
#include <readline/readline.h>
#include <readline/history.h>
//#endif

#include <o2scl/cli_readline.h>

/// A namespace for objects associated with the command-line utility 'agraph'
namespace o2scl_agraph {
  
  /** \brief The driver for 'agraph' command-line utility
   */
  class agraph_manager : public o2scl_acol::acol_manager {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Object for density plots
    o2scl_graph::table3d_multi_density tmdp;

    /// Title of the X-axis
    std::string xtitle;

    /// Title of the Y-axis
    std::string ytitle;

    /// Window title
    std::string wtitle;

    /// Second window title
    std::string wtitle2;

    /// Pad title
    std::string ptitle;

    /// If true, the x-axis limit has been set
    bool xset;

    /// If true, the y-axis limit has been set
    bool yset;

    /// If true, the z-axis limit has been set
    bool zset;

    /// Desc
    int multi_bins;
    
    /// Lower limit for X axis
    double xlo;
    
    /// Upper limit for X axis
    double xhi;

    /// Lower limit for Y axis
    double ylo;
    
    /// Upper limit for Y axis
    double yhi;

    /// Lower limit for Z axis
    double zlo;
    
    /// Upper limit for Z axis
    double zhi;

    /// \name Current axis limits
    //@{
    double a_left, a_bottom, a_right, a_top;
    //@}

    /// Left pad margin
    double plmar;
    /// Right pad margin
    double prmar;
    /// Top pad margin
    double ptmar;
    /// Bottom pad margin
    double pbmar;
    /// Line width
    int lwidth;
    /// Text alignment
    int talign;

    /// The ROOT Application object
    TApplication *theApp;

    /// Canvas
    TCanvas *int_canvas;
    /// Pad
    TPad *int_pad;
    /// TH1 object from Pad
    TH1 *int_th1;
    /// Top axis
    TGaxis *top_axis;
    /// Right axis
    TGaxis *right_axis;
    /// TH2 object for surface plots
    TH2D *surf_hist;
    /// Line index
    size_t line_ix;
    /// Color index
    size_t color_ix;
    /// Marker index
    size_t marker_ix;
    /// Text object
    TLatex tt;
    /// Standard line colors
    ubvector_int std_colors;

    /// If true, use a log scale for the x-axis
    bool logx;

    /// If true, use a log scale for the y-axis
    bool logy;

    /// If true, use a log scale for the y-axis
    bool logz;
    
    /// \name Parameter objects for \ref o2scl::cli
    //@{
    o2scl::cli::parameter_double p_xlo;
    o2scl::cli::parameter_double p_xhi;
    o2scl::cli::parameter_double p_ylo;
    o2scl::cli::parameter_double p_yhi;
    o2scl::cli::parameter_double p_zlo;
    o2scl::cli::parameter_double p_zhi;

    o2scl::cli::parameter_double p_lmar;
    o2scl::cli::parameter_double p_tmar;
    o2scl::cli::parameter_double p_rmar;
    o2scl::cli::parameter_double p_bmar;

    o2scl::cli::parameter_string p_xtitle;
    o2scl::cli::parameter_string p_ytitle;
    o2scl::cli::parameter_string p_wtitle;
    o2scl::cli::parameter_string p_wtitle2;
    o2scl::cli::parameter_string p_ptitle;

    o2scl::cli::parameter_bool p_logx;
    o2scl::cli::parameter_bool p_logy;
    o2scl::cli::parameter_bool p_logz;

    o2scl::cli::parameter_bool p_xset;
    o2scl::cli::parameter_bool p_yset;
    o2scl::cli::parameter_bool p_zset;

    o2scl::cli::parameter_int p_multi_bins;
    o2scl::cli::parameter_int p_lwidth;
    o2scl::cli::parameter_int p_talign;
    //@}

    /// Color manager
    o2scl_graph::root_color_manager rcm;

#endif
    
  public:

    agraph_manager();
    
    virtual ~agraph_manager();

    /** \brief Setup the \ref o2scl::cli object

	This function determines the <tt>HOME</tt> directory for
	the history file and instantiates the \ref o2scl::cli object.
     */
    virtual int setup_cli();
    
    /** \brief Setup options

	Setup the additional options and/or interactive commands
     */
    virtual int setup_options();

    /** \brief Setup help

	Adds additional information to the help description and also
	sets up the agraph command name and description.
     */
    virtual int setup_help();

    /** \brief Setup parameters

	Create all the new 'set' parameters.
     */
    virtual int setup_parameters();

    /** \brief Change color scheme for 'den-plot'
     */
    virtual int comm_density_colors(std::vector<std::string> &sv, 
				    bool itive_com);

    /** \brief Post-processing for the 'set' command
     */
    virtual int comm_set(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot the canvas

	Create the canvas with titles in \ref wtitle and \ref wtitle2,
	and the pad with title \ref ptitle and margins specified in
	\ref plmar, \ref prmar, \ref ptmar, and \ref pbmar. Also,
	the pad is set in log mode according to \ref logx and 
	\ref logy. 
     */
    virtual int comm_canvas(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot a function with a line
     */
    virtual int comm_plot(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot a function with a line
     */
    virtual int comm_hist_plot(std::vector<std::string> &sv, bool itive_com);

    /** \brief Surface plot of a \ref o2scl::table3d object
     */
    virtual int comm_surf_plot(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot a line between two points
     */
    virtual int comm_line(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot an arrow
     */
    virtual int comm_arrow(std::vector<std::string> &sv, bool itive_com);

    /** \brief A bar plot of a \ref o2scl::table object
     */
    virtual int comm_barplot(std::vector<std::string> &sv, bool itive_com);

    /** \brief Plot a set of tabulated values with a set of points
     */
    virtual int comm_points(std::vector<std::string> &sv, bool itive_com);

    /** \brief Run the ROOT application
     */
    virtual int comm_plotrun(std::vector<std::string> &sv, bool itive_com);

    /** \brief Internal function to plot the axis

	If no canvas has been created, this creates one with \ref
	comm_canvas() first.

	If an axis is already stored in \ref int_th1, then it is
	deleted before creating a new axis with the
	<tt>TPad::DrawFrame()</tt> function. The axis limits are
	additionally stored in <tt>a_left, a_bottom, a_right</tt>, and
	<tt>a_top</tt>.
     */
    int internal_axis(double left, double bottom, double right,
		      double top);
    
    /** \brief Create an axis with specified limits
     */
    virtual int comm_axis(std::vector<std::string> &sv, bool itive_com);

    /** \brief Add text
     */
    virtual int comm_text(std::vector<std::string> &sv, bool itive_com);

    /** \brief Store plot in a file
     */
    virtual int comm_print(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a density plot from a \ref o2scl::table3d object
     */
    virtual int comm_den_plot(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Create a density plot from more than one 
	\ref o2scl::table3d object
     */
    virtual int comm_mden_plot(std::vector<std::string> &sv, bool itive_com);
    virtual int comm_add_density(std::vector<std::string> &sv, bool itive_com);

  };

}
