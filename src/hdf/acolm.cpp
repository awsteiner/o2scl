/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>

/*
  Todos: 
  - sum/max/min/output/interp/deriv/integ/deriv2 for hist, hist_2d, and v<c>
  - xindex and yindex for table3d
  - create table3d output that can be read by gen3-list?
  - fix fit for table
  - use swap instead of copy in 'select' for table objects
  - Make sure get_input() is used more consistently
  - Make sure preview, output, internal, generic, and create work for 
  all types.
  - Add matrix types?
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

void acol_manager::command_add(std::string new_type) {

  const int both=cli::comm_option_both;
  
  if (new_type=="int") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","",
       0,1,"[value]","Get or set the value of the int object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="double") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","",
       0,1,"[value]","Get or set the value of the double object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="char") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","",
       0,1,"[value]","Get or set the value of the char object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="size_t") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","",
       0,1,"[value]","Get or set the value of the size_t object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="string") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","",
       0,1,"[value]","Get or set the value of the string object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="table") {
    static const size_t narr=33;
    comm_option_s options_arr[narr]={
      {'a',"assign","Assign a constant, e.g. assign pi acos(-1) .",
       0,2,"<name> [val]",
       ((string)"Assign a constant value to a name for the present table. ")+
       "Valid constant values are things like 1.618 or acos(-1.0) or sin(4^5). "
       "To remove an assignment, call assign with a blank value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_assign),
       both},
      {0,"delete-col","Delete a column.",0,1,"<name>",
       "Delete the entire column named <name>.",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_delete_col),both},
      {'d',"delete-rows","Delete rows selected by a function.",
       0,1,"<func>",((string)"Delete the set of rows for ")+
       "which a function evaluates to a number greater than 0.5. "+
       "For example, 'delete-rows if(col1+col2>10,1,0)' will delete "+
       "all columns where the sum of the entries in 'col1' and 'col2' "+
       "is larger than 10. See also 'select-rows'.",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_delete_rows),both},
      {'D',"deriv",
       "Derivative of a function defined by two columns.",
       0,3,"<x> <y> <name>",
       ((string)"Create a new column named <name> filled with the ")+
       "derivative of the function y(x) obtained from columns <x> and <y>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
       both},
      {0,"deriv2","Second derivative.",0,3,"<name> <x> <y>",
       ((string)"Create a new column named <name> filled with the second ")+
       "derivative of the function y(x) obtained from columns <x> and <y>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv2),
       both},
      {0,"cat",
       "Concatenate data from a second table object onto current table.",0,2,
       "<file> [name]",((string)"For table objects, add a ")+
       "second table to the end of the first, creating new columns "+
       "if necessary.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_cat),
       both},
      {0,"convert-unit","Convert a column to a new unit.",0,2,
       "<column> <new_unit>",((string)"(This command only works if ")+
       "the GNU 'units' command is installed and available in the current "+
       "path.) Convert the units of a column to <new unit>, multipliying "+
       "all entries in that column by the appropriate factor.",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_convert_unit),both},
      {0,"find-row","Find a row which maximizes a function (table only).",
       0,2,"<func> or find-row <col> <val>",
       ((string)"If one argument is given, then find-row finds the row ")+
       "which maximizes the value of the "+
       "expression given in <func>, and then output the entire row. "+
       "Otherwise find-row finds the row for which the value in "+
       "column named <col> is as close as possible to the value <val>. "+
       "See command 'get-row' to get a row by it's index.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_row),
       both},
      {0,"fit","Fit two columns to a function (experimental, table only).",0,7,
       "<x> <y> <yerr> <ynew> <par names> <func> <vals>","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_fit),
       both},
      {'f',"function","Create a new column from a function.",0,2,
       "<func> <name>",
       ((string)"Create a new column named <name> from a function (in ")+
       "<func>) in terms of the other columns. For example, for "+
       "a table containing columns named 'c1' and 'c2', 'function "+
       "c1-c2 c3' would create a new column c3 which contains the "+
       "difference of columns 'c1' and 'c2'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
       both},
      {0,"get-row","Get a row by index.",
       0,1,"<index>",((string)"Get a row by index. The first row ")+
       "has index 0, and the last row has index n-1, where n "+
       "is the total number of rows as returned by the 'list' command. "+
       "The 'index' command creates a column of row indexes. "+
       "To find a row which contains a particular value or maximizes "+
       "a specified function, use 'find-row'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_row),
       both},
      {0,"get-unit","Get the units for a specified column.",0,1,"<column>","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_unit),
       both},
      {0,"entry","Get a single entry in a table.",0,3,
       "<column> <index> [index2]","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_entry),
       both},
      {'N',"index","Add a column containing the row numbers.",0,1,
       "[column name]",
       ((string)"Define a new column named [column name] and fill ")+
       "the column with the row indexes, beginning with zero. If "+
       "no argument is given, the new column is named 'N'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_index),
       both},
      {0,"insert","Interpolate a column from another file.",0,6,
       "<file> <table> <oldx> <oldy> <newx> [newy]",
       ((string)"Insert a column from file <fname> interpolating it ")+
       "into the current table. The column <oldy> is the "+
       "columns in the file which is to be inserted into the table, "+
       "using the column <oldx> in the file and <newx> in the table. "+
       "The new column in the table is named <oldy>, or it is named "+
       "[newy] if the additional argument is given. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_insert),
       both},
      {0,"insert-full",
       "Interpolate a table from another file.",0,3,
       "<fname> <oldx> <newx>",
       ((string)"Insert all columns from file <fname> interpolating it ")+
       "into the current table. The column <oldy> is the "+
       "columns in the file which is to be inserted into the table, "+
       "using the column <oldx> in the file and <newx> in the table. "+
       "The new column are given the same names as in the file. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {'I',"integ",
       "Integrate a function specified by two columns.",
       0,3,"<x> <y> <name>",
       ((string)"Create a new column named <name> filled with the ")+
       "integral of the function y(x) obtained from columns <x> and <y>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_integ),
       both},
      {0,"interp","Interpolate a number into a column.",0,3,
       "<x name> <x value> <y name>",
       ((string)"Interpolate <x value> from column ")+
       "named <x name> into column named <y name>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
       both},
      {'l',"list","List the constants, column names and other info.",
       0,0,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"max","Find the maximum value of a column.",0,1,"<col>",
       "Compute the maximum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of a column.",0,1,"<col>",
       "Compute the minimum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {0,"rename","Rename a column.",0,2,"<old> <new>",
       "Rename a column from <old> to <new>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rename),
       both},
      {'s',"select","Select columns for a new table.",-1,-1,"<cols>",
       ((string)"Select creates a new table from the present table, ")+
       "including only the columns specified in <cols>. The column "+
       "specification is a list of column names, functions, or patterns "+
       "which match "+
       "the column names. Patterns must be preceeded by a colon ':' "+
       "and can use wildcards like '*' and '?'. All of the rows of data "+
       "are copied over. If functions are specified, the result can be "+
       "named using '='. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_select),
       both},
      {0,"select-rows","Select rows for a new table (table only).",
       0,1,"<row_spec>",((std::string)"Select the rows from a table for ")+
       "which the row specification in <row_spec> evaluates to a number "+
       "greater than 0.5",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_select_rows),both},
      {0,"select-rows2",
       "Select rows, with explicit column specification (table only).",
       0,-1,"<row_spec> [col1] [col2] ...",
       ((std::string)"Select the rows from a table for ")+
       "which the row specification in <row_spec> evaluates to a number "+
       "greater than 0.5 . All of the columns required to compute the row "+
       "specification must be given in [col1] [col2] ... This can be "+
       "faster than 'select-rows' for tables with many columns.",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_select_rows2),both},
      {0,"set-data","Set the entries of a column.",3,4,
       "<row_spec> <col> <val_spec>",
       ((string)"Set the value of rows specifed by the ")+
       "'row_spec' function in column 'col' to the value given by the "+
       "'val_spec' function. Rows are chosen if row_spec evaluates to a "+
       "number greater than 0.5.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_data),
       both},
      {0,"set-unit","Set the units for a specified column.",0,2,
       "<column> <unit>","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_unit),
       both},
      {'S',"sort","Sort the entire table by a column.",0,2,
       "<col> [unique]",
       ((string)"Sorts the entire table by the column specified in <col>. ")+
       "If the word \"unique\" is specified as the second argument, then "+
       "delete duplicate rows after sorting.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"stats","Show column statistics.",0,1,"<col>",
       "Output the average, std. dev, max and min of <col>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_stats),
       both},
      {0,"sum","Add data from a second table object to current table.",
       0,2,"<file> [name]",((string)"Add all columns ")+
       "from the second table to their corresponding columns in the "+
       "current table, creating new columns if necessary.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
       both},
      {0,"nlines","Add 'nlines' as a constant to a table object.",0,0,"",
       "",new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_nlines),
       both},
      {0,"to-hist","Convert a table to a histogram.",0,6,
       "[\"2d\"] <col 1> <col_2 (2d only)> <x_bins> <y_bins (2d only)> [wgts]",
       ((std::string)"The 'to-hist' command takes two forms: 'to-hist ")+
       "col1 N1 wgts' or 'to-hist col2 N1 N2 wgts'. The first creates "+
       "a 1D histogram from 'col1' using exactly 'N1' bins and "+
       "(optionally) weighting the entries by the values in column 'wgts'. "+
       "The second form creates a 2D histogram from 'col1' and 'col2' "+
       "using N1 bins in the x direction and N2 bins in the y direction, "+
       "optionally weighting the entries by the column 'wgts'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_hist),
       both},
      {0,"autocorr","Compute the autocorrelation vectors.",0,3,
       "<col> <ac> <ftom>",
       ((std::string)"Given a column <col>, this stores a vector of ")+
       "autocorrelation coefficients in column <ac> and the quantity "+
       "'5*tau/M' in column <ftom>. Columns <ac> and <ftom> are created "+
       "if they are not already present and overwritten if they "+
       "already contain data. Also, the autocorrelation length and "+
       "estimated sample size are output to the screen.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);

    if (o2graph_mode) {
      static const size_t narr2=5;
      comm_option_s options_arr2[narr2]={
	{0,"plot","Plot two columns from the table.",0,3,"<x> <y> [kwargs]",
	 ((std::string)"If the current object is a table, then plot ")+
	 "column <y> versus column <x>. If the current object is a "+
	 "one-dimensional histogram, then plot the histogram weights "+
	 "as a function of the bin representative valus. If the "+
	 "current object is a set of contour lines, then plot the "+
	 "full set of contour lines. "+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms). For example: o2graph -create x 0 10 0.2 -function \"sin(x)\" "+
	 "y -plot x y lw=0,marker='+' -show",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both},
	{0,"errorbar","Plot the specified columns with errobars.",4,5,
	 "<x> <y> <xerr> <yerr> [kwargs]",
	 ((std::string)"Plot column <y> versus column <x> with symmetric ")+
	 "error bars given in column <xerr> and <yerr>. For no uncertainty "
	 "in either the x or y direction, just use \"0\" for <xerr> "+
	 "or <yerr>, respectively. New kwargs for the errorbar command are "+
	 "ecolor=None, elinewidth=None, capsize=None, barsabove=False, "+
	 "lolims=False, uplims=False, xlolims=False, xuplims=False, "+
	 "errorevery=1, capthick=None, hold=None",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both},
	{0,"plot1","Plot the specified column.",1,2,"<y> [kwargs]",
	 ((std::string)"Plot column <y> versus row number. Some useful ")+
	 "kwargs are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms). For example: o2 -create x 0 10 0.2 -function \"sin(x)\" "+
	 "y -plot1 y ls='--',marker='o' -show",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both},
	{0,"histplot","Create a histogram plot from a table.",0,1,"<col>",
	 ((std::string)"For a table, create a histogram plot from the ")+
	 "specified column.",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both},
	{0,"hist2dplot","Create a 2-D histogram plot.",0,2,"<col x> <col y>",
	 ((std::string)"For a table, create a 2D histogram plot from the ")+
	 "specified columns.",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),
	 both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="table3d") {
    
    static const size_t narr=15;
    comm_option_s options_arr[narr]={
      {0,"cat",
       "Concatenate data from a second table3d onto current table3d.",0,2,
       "<file> [name]",((string)"Add all slices from the ")+
       "second table3d object which aren't already present in the "+
       "current table3d object.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_cat),
       both},
      {0,"contours","Create contour lines from a table3d slice.",
       0,4,"[\"frac\"] <value> <slice_name> [output file] [output name]","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_contours),
       both},
      {0,"deriv-x","Derivative with respect to x.",0,2,
       "<f> <dfdx>",
       ((string)"Create a new slice named <dfdx> filled with the ")+
       "derivative of the function from the x grid and slice named <f>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv_x),
       both},
      {0,"deriv-y","Derivative with respect to y.",0,2,
       "<f> <dfdy>",
       ((string)"Create a new slice named <dfdy> filled with the ")+
       "derivative of the function from the y grid and slice named <f>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv_y),
       both},
      {'f',"function","Create a new slice from a function.",0,2,
       "<func> <name>",
       ((string)"Create a new slice named <name> from a function (in ")+
       "<func>) in terms of the other slices. For example, for "+
       "a table3d containing slices named 's1' and 's2', 'function "+
       "s1-s2 s3' would create a new column 's3' which contains the "+
       "difference of columns 's1' and 's2'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
       both},
      {0,"entry","Get a single entry in a table3d.",0,3,
       "<slice> <x index> <y index>","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_entry),
       both},
      {0,"insert","Interpolate a slice from another file.",0,6,
       "<file> <table> <old> [new]","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_insert),
       both},
      {0,"interp","Interpolate a number into a slice.",0,3,
       "<z name> <x value> <y value> ",
       "Interpolate (<x value>,<y value>) into the slice named <z name>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
       both},
      {'l',"list","List the slice names and print out grid info.",
       0,0,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"max","Find the maximum value of a slice.",0,1,"<slice>",
       "Compute the maximum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of a slice.",0,1,"<slice>",
       "Compute the minimum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {0,"rename","Rename a slice.",0,2,"<old> <new>",
       "Rename a slice from <old> to <new>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rename),
       both},
      {0,"set-data","Set the entries of a column.",3,4,
       "<x value> <y value> <z name> <val>",
       ((string)"Set the value of ")+
       "the slice named 'z name' at the grid point closest to "+
       "(<x value>,<y value>) to the value <val>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_data),
       both},
      {0,"slice","Construct a slice.",2,2,
       "<\"x\" or \"y\"> <value>",
       ((string)"Extract a slice of a table3d object at fixed x or fixed y ")+
       "to create a new table object. This function uses interpolation "+
       "with the current interpolation type to interpolate all of the "+
       "slices in the table3d object to create a table with a column "+
       "for each slice.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_slice),
       both},
      {0,"sum","Add data from a second table3d object to current table3d.",
       0,2,"<file> [name]",((string)"Add all slides from the ")+
       "second table3d to their "+
       "corresponding slices in the current table3d, creating new slices "+
       "if necessary.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);

    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"den-plot","Create a density plot from a table3d object.",
	 0,1,"<slice name for table3d>",
	 ((std::string)"Creates a density plot from the specified slice. ")+
	 "Logarithmic x- or y-axes are handled by taking the base 10 log "+
	 "of the x- or y-grids specified in the table3d object before "+
	 "plotting. A z-axis density legend is print on the RHS if "+
	 "colbar is set to 1 before plotting. If z-axis limits are "+
	 "specified, then values larger than the upper limit are set "+
	 "equal to the upper limit and values smaller than the lower "+
	 "limit are set equal to the lower limit before plotting.",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),
	 both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }

  } else if (new_type=="hist") {

    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"plot","Plot the histogram.",0,3,"[kwargs]",
	 ((std::string)"Plot the histogram weights ")+
	 "as a function of the bin representative values. "+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms). For example: o2graph -create x 0 10 0.2 -function \"sin(x)\" "+
	 "y -plot x y lw=0,marker='+' -show",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both},
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="double[]") {

    static const size_t narr=7;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {'D',"deriv",
       "Replace the array with its derivative.",0,0,"",
       ((string)"Replace the array with its derivative using the ")+
       "current interpolation type.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
       both},
      {0,"interp","Interpolate an index into the array.",0,1,
       "<x value>",
       ((string)"Interpolate <x value> in the array."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
       both},
      {0,"autocorr","Compute the autocorrelation vector and length.",0,3,
       "<col> <ac> <ftom>",
       "Compute the autocorrelation vector and length.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both},      
      {0,"to-table","Convert to a table given a column name",0,1,
       "<column name>",
       "Convert to a table given a column name.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_table),
       both}      
    };
    cl->set_comm_option_vec(narr,options_arr);
    
    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"plot1",
	 "Plot the array",0,1,"[kwargs]",
	 ((std::string)"Plot the array. ")+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms).",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="int[]") {

    static const size_t narr=7;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {'D',"deriv",
       "Replace the array with its derivative.",0,0,"",
       ((string)"Replace the array with its derivative using the ")+
       "current interpolation type.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
       both},
      {0,"interp","Interpolate an index into the array.",0,1,
       "<x value>",
       ((string)"Interpolate <x value> in the array."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
       both},
      {0,"autocorr","Compute the autocorrelation vector and length.",0,3,
       "<col> <ac> <ftom>",
       "Compute the autocorrelation vector and length.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both},
      {0,"to-table","Convert to a table given a column name",0,1,
       "<column name>",
       "Convert to a table given a column name.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_table),
       both}      
    };
    cl->set_comm_option_vec(narr,options_arr);
    
    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"plot1",
	 "Plot the array",0,1,"[kwargs]",
	 ((std::string)"Plot the array. ")+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms).",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="size_t[]") {

    static const size_t narr=7;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value of column <col>.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {'D',"deriv",
       "Replace the array with its derivative.",0,0,"",
       ((string)"Replace the array with its derivative using the ")+
       "current interpolation type.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
       both},
      {0,"interp","Interpolate an index into the array.",0,1,
       "<x value>",
       ((string)"Interpolate <x value> in the array."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
       both},
      {0,"autocorr","Compute the autocorrelation vector and length.",0,3,
       "<col> <ac> <ftom>",
       "Compute the autocorrelation vector and length.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
       both},
      {0,"to-table","Convert to a table given a column name",0,1,
       "<column name>",
       "Convert to a table given a column name.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_table),
       both}      
    };
    cl->set_comm_option_vec(narr,options_arr);
    
    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"plot1",
	 "Plot the array",0,1,"[kwargs]",
	 ((std::string)"Plot the array. ")+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms).",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="vector<contour_line>") {

    if (o2graph_mode) {
      static const size_t narr2=1;
      comm_option_s options_arr2[narr2]={
	{0,"plot",
	 "Plot the contour lines",0,1,"[kwargs]",
	 ((std::string)"Plot the set of contour lines. ")+
	 "Some useful kwargs (which apply for all three object types) "+
	 "are color (c), dashes, linestyle (ls), linewidth (lw), "+
	 "marker, markeredgecolor (mec), markeredgewidth (mew), "+
	 "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
	 "(ms). For example: o2graph -create x 0 10 0.2 -function \"sin(x)\" "+
	 "y -plot x y lw=0,marker='+' -show",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }
    
  } else if (new_type=="hist_2d") {

    static const size_t narr=2;
    comm_option_s options_arr[narr]={
      {0,"max","Find the maximum weight.",0,0,"",
       "Find the maximum weight",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum weight.",0,0,"",
       "Find the minimum weight.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);

    if (o2graph_mode) {
      static const size_t narr2=2;
      comm_option_s options_arr2[narr2]={
	{0,"den-plot","Create a density plot from a hist_2d object.",
	 0,1,"<slice name for table3d>",
	 ((std::string)"Create ")+
	 "a density plot from the current histogram (assuming equally-"+
	 "spaced bins). "+
	 "Logarithmic x- or y-axes are handled by taking the base 10 log "+
	 "of the x- or y-grids specified in the table3d object before "+
	 "plotting. A z-axis density legend is print on the RHS if "+
	 "colbar is set to 1 before plotting. If z-axis limits are "+
	 "specified, then values larger than the upper limit are set "+
	 "equal to the upper limit and values smaller than the lower "+
	 "limit are set equal to the lower limit before plotting.",
	 new o2scl::comm_option_mfptr<acol_manager>
	 (this,&acol_manager::comm_none),
	 both},
	{0,"contours","Create contour lines from a table3d or hist_2d object.",
	 0,4,"[\"frac\"] <value> [output file] [output name]","",
	 new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_contours),
	 both}
      };
      cl->set_comm_option_vec(narr2,options_arr2);
    }

  }
  
  return;
}

acol_manager::acol_manager() : cset(this,&acol_manager::comm_set),
			       cng(o2scl_settings.get_convert_units()) {
  
  obj_name="acol";
  verbose=1;
  pretty=true;
  names_out=true;
  scientific=true;
  prec=6;
  user_ncols=-1;
  unit_fname="";
  def_args="";
  
  post_interactive=false;

  ffl.html_mode();

  type="";

  env_var_name="ACOL_DEFAULTS";
  interp_type=1;

  o2graph_mode=false;

#ifdef O2SCL_HDF5_COMP
  compress=1;
#else
  compress=0;
#endif

  cng.err_on_fail=false;
}

void acol_manager::command_del() {

  if (type=="int") {
    cl->remove_comm_option("value");
  } else if (type=="double") {
    cl->remove_comm_option("value");
  } else if (type=="char") {
    cl->remove_comm_option("value");
  } else if (type=="size_t") {
    cl->remove_comm_option("value");
  } else if (type=="string") {
    cl->remove_comm_option("value");
  } else if (type=="table") {
    cl->remove_comm_option("assign");
    cl->remove_comm_option("delete-col");
    cl->remove_comm_option("delete-rows");
    cl->remove_comm_option("deriv");
    cl->remove_comm_option("deriv2");
    cl->remove_comm_option("cat");
    cl->remove_comm_option("convert-unit");
    cl->remove_comm_option("find-row");
    cl->remove_comm_option("fit");
    cl->remove_comm_option("function");
    cl->remove_comm_option("get-row");
    cl->remove_comm_option("get-unit");
    cl->remove_comm_option("entry");
    cl->remove_comm_option("index");
    cl->remove_comm_option("insert");
    cl->remove_comm_option("insert-full");
    cl->remove_comm_option("integ");
    cl->remove_comm_option("interp");
    cl->remove_comm_option("list");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    cl->remove_comm_option("rename");
    cl->remove_comm_option("select");
    cl->remove_comm_option("select-rows");
    cl->remove_comm_option("select-rows2");
    cl->remove_comm_option("set-data");
    cl->remove_comm_option("set-unit");
    cl->remove_comm_option("sort");
    cl->remove_comm_option("stats");
    cl->remove_comm_option("sum");
    cl->remove_comm_option("nlines");
    cl->remove_comm_option("to-hist");
    cl->remove_comm_option("autocorr");

    if (o2graph_mode) {
      cl->remove_comm_option("plot");
      cl->remove_comm_option("plot1");
      cl->remove_comm_option("errorbar");
      cl->remove_comm_option("histplot");
      cl->remove_comm_option("hist2dplot");
    }
    
  } else if (type=="table3d") {
    
    cl->remove_comm_option("cat");
    cl->remove_comm_option("contours");
    cl->remove_comm_option("deriv-x");
    cl->remove_comm_option("deriv-y");
    cl->remove_comm_option("function");
    cl->remove_comm_option("entry");
    cl->remove_comm_option("insert");
    cl->remove_comm_option("interp");
    cl->remove_comm_option("list");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    cl->remove_comm_option("rename");
    cl->remove_comm_option("set-data");
    cl->remove_comm_option("slice");
    cl->remove_comm_option("sum");

    if (o2graph_mode) {
      cl->remove_comm_option("den-plot");
    }
    
  } else if (type=="hist_2d") {

    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    
    /*
      cl->remove_comm_option("deriv-x");
      cl->remove_comm_option("deriv-y");
      cl->remove_comm_option("interp");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("sum");
      
      //cl->remove_comm_option("hist2d");
      */
    if (o2graph_mode) {
      cl->remove_comm_option("den-plot");
      cl->remove_comm_option("contours");
    }
    
  } else if (type=="hist") {
    
    /*
      cl->remove_comm_option("assign");
      cl->remove_comm_option("deriv");
      cl->remove_comm_option("deriv2");
      cl->remove_comm_option("find-row");
      cl->remove_comm_option("fit");
      cl->remove_comm_option("entry");
      cl->remove_comm_option("integ");
      cl->remove_comm_option("max");
      cl->remove_comm_option("min");
      cl->remove_comm_option("rename");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("stats");
      cl->remove_comm_option("sum");
      cl->remove_comm_option("to-table");

      cl->remove_comm_option("plot");
      cl->remove_comm_option("plot1");
      cl->remove_comm_option("hist");
    */
    
    if (o2graph_mode) {
      cl->remove_comm_option("plot");
    }
    
  } else if (type=="vector<contour_line>") {

    if (o2graph_mode) {
      cl->remove_comm_option("plot");
    }
    
  } else if (type=="double[]" || type=="int[]" || type=="size_t[]") {

    cl->remove_comm_option("deriv");
    cl->remove_comm_option("interp");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    cl->remove_comm_option("sort");
    cl->remove_comm_option("autocorr");
    cl->remove_comm_option("to-table");
    
    if (o2graph_mode) {
      cl->remove_comm_option("plot1");
    }

    /*
      cl->remove_comm_option("integ");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("stats");
      cl->remove_comm_option("sum");
      cl->remove_comm_option("to-hist");
    */

    /*
      } else if (type=="int[]" || type=="size_t[]") {
      
      //cl->remove_comm_option("plot1");
      cl->remove_comm_option("deriv");
      cl->remove_comm_option("deriv2");
      cl->remove_comm_option("integ");
      cl->remove_comm_option("max");
      cl->remove_comm_option("min");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("sort");
      cl->remove_comm_option("stats");
      cl->remove_comm_option("sum");
      cl->remove_comm_option("to-hist");
    */

  }
  
  return;
}

void acol_manager::clear_obj() {

  if (type=="table") {
    table_obj.clear();
  } else if (type=="table3d") {
    table3d_obj.clear();
  } else if (type=="hist") {
    hist_obj.clear();
  } else if (type=="hist_2d") {
    hist_2d_obj.clear();
  } else if (type=="vector<contour_line>") {
    cont_obj.clear();
  } else if (type!="string[]") {
    stringv_obj.clear();
  } else if (type!="int[]") {
    intv_obj.clear();
  } else if (type!="double[]") {
    doublev_obj.clear();
  } else if (type!="string") {
    string_obj.clear();
  } else if (type!="size_t[]") {
    size_tv_obj.clear();
  }
  
  type="";
  
  return;
}

int acol_manager::setup_options() {

  const int cl_param=cli::comm_option_cl_param;
  const int both=cli::comm_option_both;

  static const int narr=14;

  // Options, sorted by long name. We allow 0 parameters in many of these
  // options so they can be requested from the user in interactive mode. 
  comm_option_s options_arr[narr]={
    {0,"calc","Compute the value of a constant expression.",0,1,"<expr>",
     ((string)"This computes the value of the constant expression ")+
     " <expr>. Examples are 'calc acos(-1)' or 'calc 2+1/sqrt(2.0e4)'. "+
     "Results are given at the current precision.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_calc),
     both},
    {'c',"create","Create an object",0,-1,"<type> [...]",
     ((string)"Create a new object of type <type>. For types char, ")+
     "double, int, size_t, and string, this takes one additional "+
     "argument which holds the value. For type table, "+
     "create a new table with one column whose entries "+
     "are an evenly-spaced grid. This takes four arguments, the name of "+
     "the column, the first value, the increment between successive values "+
     "and the maximum possible value. Note that finite precision "+
     "arithmetic may cause small deviations from the expected result. "+
     "If a table is currently in memory, it is deallocated beforehand. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
     both},
    {0,"download","Download file from specified URL.",
     0,3,"<file> <hash file> <URL>",
     ((string)"Check if a file matches a specified hash, and if not,")+
     "attempt to download a fresh copy from the specified URL.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_download),
     both},
    {0,"filelist","List objects in a HDF5 file.",0,1,"<file>",
     ((string)"This lists all the top-level datasets and groups in a ")+
     "HDF5 file and, for those groups which are in the O2scl format, "+
     "gives the type and name of the object stored in that HDF5 group.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_filelist),
     both},
    {'g',"generic","Read in a generic text file.",0,2,"<type> <file>",
     ((string)"Read an object of type <type> from a text file named ")+
     "<file>. The first line of the file must either contain numeric "+
     "data or column names "+
     "separated by white space, without carriage returns, except for "+
     "the one at the end of the line. All remaining lines are assumed "+
     "to contain data with the same number of columns as the first line. "+
     "The data must be stored in columns "+
     "where the first entry in each column is the x-axis grid point, "+
     "the second, entry in each column is the y-axis grid point, "+
     "and the remaining columns give the data for each slice at that "+
     "point. Each grid point must correspond to a line in the file, but "+
     "the lines need not be in any particular order. The columns may "+
     "have one header line at top which specifies the names of the x- "+
     "and y-grids and the names of each slice (in order).",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
     both},
    {0,"get-conv","Get a unit conversion factor.",0,2,
     "<old unit> <new unit>",((string)"This command gets a unit ")+
     "conversion factor. It only works if the conversion is one of the "
     "hard-coded O2scl conversions or if "+
     "the GNU 'units' command is installed and available in the current "+
     "path. For example, 'get-conv MeV erg' returns 1.782662e-23 and 1 MeV "+
     "is equivalent to 1.782662e-23 erg. The conversion factor is output "+
     "at the current precision, but is always internally stored with "+
     "full double precision. (Note that the 'units' command shipped on OSX "+
     "is not the same as the GNU version.)",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_conv),
     both},
    /*    
	  {'H',"html","Create a file in HTML (table3d only).",0,1,"<file>",
	  "Output the current table in HTML mode to the specified file. ",
	  new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_html),
	  both},
    */
    {'q',"interactive","Toggle the interactive interface.",
     0,0,"",((string)"If given as a command-line parameter, 'interactive' ")+
     "toggles the execution of the interactive mode after the "+
     "command-line parameters are processed. If zero arguments are given "+
     "to 'acol' on the command-line then the interactive interface is "+
     "automatically turned on.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interactive),
     cl_param},
    {'i',"internal","Output current object in the internal HDF5 format.",
     0,1,"[file]",
     ((string)"Output the current object in the internal HDF5 format. ")+
     "If no argument is given, then output is sent to the screen, "+
     "otherwise, output is sent to the specified file. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_internal),
     both},
    {'o',"output","Output the current object as text.",0,1,"[file]",
     ((string)"Output the object to the screen, or if the [file] ")+
     "argument is specified, to a file. This is the same format as "+
     "can be read using the 'generic' command.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_output),
     both},
    {'P',"preview","Preview the current object.",0,2,
     "[number of lines] [number of columns]",
     ((string)"Print out all or part of the current object in format ")+
     "suitable for the screen.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_preview),
     both},
    {'r',"read","Read an object from an O2scl-style HDF5 file.",0,2,
     "<file> [object name]",
     ((string)"Read an HDF5 file with the specified filename. ")+
     "If the [object name] argument is specified, then read the object "+
     "with the specified name. Otherwise, look for the first table object, "+
     "and if not found, look for the first table3d object, and so on, "+
     "attempting to find a readable O2scl object.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_read),
     both},
    {0,"show-units","Show the unit conversion table.",0,0,"",
     ((string)"This command does not show all possible conversions, only ")+
     "the conversions which have been previously used and are now stored "+
     "in the unit cache.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_show_units),
     both},
    {0,"type","Show current object type.",0,0,"",
     ((string)"Show the current object type, either table, ")+
     "table3d, hist, hist_2d, vector<contour_line>, int, double, "+
     "char, string, int[], double[], string[], size_t, size_t[], or "+
     "uniform_grid<double>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_type),
     both},
    {'v',"version","Print version information and O2scl settings.",0,0,"",
     "",new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_version),
     both}
  };
  /*
    {0,"find-x","Find an entry in the x-grid (3d only)",1,1,"<value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_x),
    comm_option::both},
    {0,"find-y","Find an entry in the y-grid (3d only)",1,1,"<value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_y),
    comm_option::both},
    {0,"find-xy",
    "Find the closest grid point to a pair of values (3d only)",
    1,1,"<x value> <y value>","",
    new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_y),
    comm_option::both},
  */

  cl->set_comm_option_vec(narr,options_arr);

  if (o2graph_mode) {
    
    const int cl_param=o2scl::cli::comm_option_cl_param;
    const int both=o2scl::cli::comm_option_both;
    
    static const size_t narr2=18;
    o2scl::comm_option_s options_arr2[narr2]={
      {0,"line","Plot a line.",4,5,"<x1> <y1> <x2> <y2> [kwargs]",
       ((std::string)"Plot a line from (x1,y1) to (xy,y2). Some useful ")+
       "kwargs are color (c), dashes, linestyle (ls), linewidth (lw), "+
       "marker, markeredgecolor (mec), markeredgewidth (mew), "+
       "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
       "(ms). For example: o2graph -line 0.05 0.05 0.95 0.95 "+
       "lw=0,marker='+' -show",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"rect","Plot a rectangle.",5,6,"<x1> <y1> <x2> <y2> <angle> [kwargs]",
       ((std::string)"Plot a rectange from (x1,y1) to (xy,y2) with ")+
       "rotation angle <angle>.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"new-cmaps","Define new color maps.",0,0,"",
       ((std::string)"Define new color maps, 'jet2', 'pastel2' ")+
       "'reds2', 'greens2', and 'blues2'.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"plotm","Plot the specified columns from tables in multiple files.",
       2,3,"<x> <y>",((std::string)"After using -plot-files to specify ")+
       "a list of files, plot column <y> versus column <x> for the first "+
       "table object in all of the specified files. Some useful "+
       "kwargs are color (c), dashes, linestyle (ls), linewidth (lw), "+
       "marker, markeredgecolor (mec), markeredgewidth (mew), "+
       "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
       "(ms). For example: o2graph -plot-files file1.o2 file2.o "+
       "-plotm xcol ycol lw=0,marker='+' -show",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"plot1m","Plot the specified column from tables in multiple files.",
       1,1,"<y>",((std::string)"After using -plot-files to specify ")+
       "a list of files, plot column <y> versus row number for the first "+
       "table object in all of the specified files. Some useful "+
       "kwargs are color (c), dashes, linestyle (ls), linewidth (lw), "+
       "marker, markeredgecolor (mec), markeredgewidth (mew), "+
       "markerfacecolor (mfc), markerfacecoloralt (mfcalt), markersize "+
       "(ms). For example: o2graph -plot-files file1.o2 file2.o "+
       "-plot1m ycol lw=0,marker='+' -show",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"xlimits","Set the x-axis limits.",2,2,"<low> <high>",
       ((std::string)"Set 'xlo' and 'xhi' to the specified limits, ")+
       "and set 'xset' to true. If a plotting canvas is currently open, then "+
       "the x-limits on that plot are modified. Future plots are also "+
       "set with the specified x-limits.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"reset-xlim","Reset the x-axis limits.",0,0,"",
       ((std::string)"This is an alias for 'set xset 0', and indicates ")+
       "that the values of 'xlo' and 'xhi' are to be ignored until the "
       "next call to 'xlimits'.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"ylimits","Set the y-axis limits.",2,2,"<low> <high>",
       ((std::string)"Set 'ylo' and 'yhi' to the specified limits, ")+
       "and set 'yset' to true. If a plotting canvas is currently open, then "+
       "the y-limits on that plot are modified. Future plots are also "+
       "set with the specified y-limits.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"reset-ylim","Reset the y-axis limits.",0,0,"",
       ((std::string)"This is an alias for 'set yset 0', and indicates ")+
       "that the values of 'ylo' and 'yhi' are to be ignored until the "
       "neyt call to 'ylimits'.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"canvas","Create a plotting canvas.",0,0,"",
       "Create an empty plotting canvas. For example 'o2graph -canvas -show'.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"clf","Clear the current figure.",0,0,"",
       "Clear the current figure.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"backend","Select the matplotlib backend to use.",1,1,"<backend>",
       ((std::string)"This selects the matplotlib backend. ")+
       "Typical values are 'Agg', 'TkAgg', 'WX', 'QTAgg', 'QT4Agg'. "+
       "Use -backend Agg to save the plot to a file without "+
       "opening a window.",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"move-labels","Move the labels.",0,0,"","",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"show","Show the current plot.",0,0,"",
       ((std::string)"Show the current plot on the screen and begin ")+
       "the graphical user interface. This is similar to plot.show().",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"text","Plot text in the axis coordinate system.",3,4,
       "<x> <y> <text> [kwargs]","",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"ttext","Plot text in the canvas default coordinate system.",3,4,
       "<x> <y> <text> [kwargs]","",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both},
      {0,"save","Save the current plot in a file.",1,1,"<filename>",
       ((std::string)"Save the current plot in a file similar ")+
       "to plot.savefig(). The action of this command depends on "+
       "which backend was selected. File type depends on the "+
       "extension, typically either .png, .pdf, .eps, .jpg, .raw, .svg, "+
       "and .tif .",new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),both},
      {0,"plot-files","Specify a list of files for 'plotm' and 'plot1m'.",1,-1,
       "<file 1> [file 2] ...","",
       new o2scl::comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_none),
       both}
    };

    cl->set_comm_option_vec(narr2,options_arr2);
  }
  
  
  return 0;
}

int acol_manager::setup_cli() {

  //---------------------------------------------------------------------
  // Get HOME directory and command history
  
  char *hd=getenv("HOME");
  std::string histfile;
  if (hd) {
    histfile=hd;
    histfile+="/.acol_hist";
  }
    
  //--------------------------------------------------------------------
  // Specify command-line option object
  
#ifdef O2SCL_READLINE
  cl=new cli_readline(histfile);
#else
  cl=new cli;
#endif

  return 0;
}

int acol_manager::setup_help() {

  if (o2graph_mode) {
    
    cl->cmd_name="o2graph";
    
    this->cl->desc=((string)"o2graph: A data viewing and ")+
      "processing program for O2scl.\n";

  } else {

  cl->cmd_name="acol";

  cl->desc=((string)"acol: A data viewing and ")+
    "processing program for O2scl.\n";

  }
  
  string dsc="\nNotes:\n\n";
  dsc+="1. Help for individual commands may be obtained with 'help ";
  dsc+="command'.\n   Required arguments are surrounded by <>'s and ";
  dsc+="optional arguments are\n   surrounded by []'s.\n";
  dsc+="2. Options may also be specified in the environment variable ";
  dsc+="ACOL_DEFAULTS.\n";
  dsc+="3. Long options may be preceeded by two dashes.\n";
  dsc+="4. In order to avoid confusion between arguments and functions,\n";
  dsc+="   use \"(-x*2)\" not \"-x*2\"\n\n";
  dsc+="Known operators:\n() ^ * / % + - == != < > && || << >> >= <=\n\n";
  dsc+="Known functions:\n";
  dsc+="exp(x) log(x) log10(x) sin(x) cos(x) tan(x) sqrt(x) abs(x)\n";
  dsc+="asin(x) acos(x) atan(x) sinh(x) cosh(x) tanh(x)\n";
  dsc+="asinh(x) acosh(x) atanh(x)\n\n";
  /*
    dsc+="atan2(x,y) if(x,y,z)\n";
    dsc+="cot(x) csc(x) sec(x)\n";
    dsc+="ceil(x) floor(x) int(x) max(x,y) min(x,y)\n";
  */
  
#ifdef O2SCL_UBUNTU_PKG
  dsc+=((string)"Compiled at ")+((string)__TIME__)+" on "+
    ((string)__DATE__)+" for "+((string)PACKAGE)+", version "+
    ((string)VERSION)+".\n";
#else
  dsc+=((string)"Compiled for ")+((string)PACKAGE)+", version "+
    ((string)VERSION)+".\n";
#endif
  
  cl->addl_help_cmd=dsc;
  cl->addl_help_cli=dsc;

  return 0;
}

int acol_manager::setup_parameters() {
  
  p_obj_name.str=&obj_name;
  p_unit_fname.str=&unit_fname;
  p_def_args.str=&def_args;
  p_verbose.i=&verbose;
  p_compress.i=&compress;
  p_prec.i=&prec;
  p_ncols.i=&ncols;
  p_interp_type.i=&interp_type;
  p_scientific.b=&scientific;
  p_pretty.b=&pretty;
  p_names_out.b=&names_out;
  
  p_obj_name.help="The current object name.";
  p_unit_fname.help="The unit filename.";
  p_def_args.help="The default arguments from the environment.";
  p_prec.help="The numerical precision.";
  p_verbose.help="Control the amount of output.";
  p_compress.help="If true, enable compression.";
  p_ncols.help="The number of output columns.";
  p_interp_type.help=((std::string)"The interpolation type ")+
    "(1=linear, 2=cubic spline, 3=periodic cubic spline, 4=Akima, "+
    "5=periodic Akima, 6=monotonic, 7=Steffen's monotonic).";
  p_names_out.help="If true, output column names at top.";
  p_pretty.help="If true, align the columns using spaces.";
  p_scientific.help="If true, output in scientific mode.";
  
  cl->par_list.insert(make_pair("obj_name",&p_obj_name));
  cl->par_list.insert(make_pair("unit_fname",&p_unit_fname));
  cl->par_list.insert(make_pair("def_args",&p_def_args));
  cl->par_list.insert(make_pair("precision",&p_prec));
  cl->par_list.insert(make_pair("verbose",&p_verbose));
  cl->par_list.insert(make_pair("compress",&p_compress));
  cl->par_list.insert(make_pair("ncols",&p_ncols));
  cl->par_list.insert(make_pair("interp_type",&p_interp_type));
  cl->par_list.insert(make_pair("names_out",&p_names_out));
  cl->par_list.insert(make_pair("pretty",&p_pretty));
  cl->par_list.insert(make_pair("scientific",&p_scientific));

  if (o2graph_mode) {
    xtitle="";
    ytitle="";
    xlo=0.0;
    xhi=0.0;
    ylo=0.0;
    yhi=0.0;
    zlo=0.0;
    zhi=0.0;
    xset=false;
    yset=false;
    zset=false;
    colbar=false;
    logx=false;
    logy=false;
    right_margin=0.04;
    left_margin=0.014;
    top_margin=0.04;
    bottom_margin=0.012;
    font=16.0;

    p_xtitle.str=&xtitle;
    p_xtitle.help=((std::string)"X-axis title. Latex ")+
      "works, e.g. '$\\phi$' and '$\\hat{x}$'";
    cl->par_list.insert(make_pair("xtitle",&p_xtitle));
    
    p_ytitle.str=&ytitle;
    p_ytitle.help=((std::string)"Y-axis title. Latex ")+
      "works, e.g. '$\\phi$' and '$\\hat{x}$'.";
    cl->par_list.insert(make_pair("ytitle",&p_ytitle));
    
    p_xlo.d=&xlo;
    p_xlo.help="Lower limit for x-axis (function if starts with '(').";
    cl->par_list.insert(make_pair("xlo",&p_xlo));
    
    p_xhi.d=&xhi;
    p_xhi.help="Upper limit for x-axis (function if starts with '(').";
    cl->par_list.insert(make_pair("xhi",&p_xhi));
    
    p_ylo.d=&ylo;
    p_ylo.help="Lower limit for y-axis (function if starts with '(').";
    cl->par_list.insert(make_pair("ylo",&p_ylo));
    
    p_yhi.d=&yhi;
    p_yhi.help="Upper limit for y-axis (function if starts with '(').";
    cl->par_list.insert(make_pair("yhi",&p_yhi));
    
    p_zlo.d=&zlo;
    p_zlo.help="Lower limit for z-axis (function if starts with '(').";
    cl->par_list.insert(make_pair("zlo",&p_zlo));
    
    p_zhi.d=&zhi;
    p_zhi.help="Upper limit for z-axis.";
    cl->par_list.insert(make_pair("zhi",&p_zhi));
    
    p_xset.b=&xset;
    p_xset.help="If true, x-axis limits have bion if starts with '(').";
    cl->par_list.insert(make_pair("zlo",&p_zlo));
    
    p_zhi.d=&zhi;
    p_zhi.help="Upper limit for z-axis.";
    cl->par_list.insert(make_pair("zhi",&p_zhi));
    
    p_xset.b=&xset;
    p_xset.help="If true, x-axis limits have been set.";
    cl->par_list.insert(make_pair("xset",&p_xset));
    
    p_yset.b=&yset;
    p_yset.help="If true, y-axis limits have been set.";
    cl->par_list.insert(make_pair("yset",&p_yset));
    
    p_zset.b=&zset;
    p_zset.help="If true, z-axis limits have been set.";
    cl->par_list.insert(make_pair("zset",&p_zset));
    
    p_colbar.b=&colbar;
    p_colbar.help="If true, den-plot adds a color legend.";
    cl->par_list.insert(make_pair("colbar",&p_colbar));
    
    p_logx.b=&logx;
    p_logx.help="If true, use a logarithmic x-axis.";
    cl->par_list.insert(make_pair("logx",&p_logx));
    
    p_logy.b=&logy;
    p_logy.help="If true, use a logarithmic y-axis.";
    cl->par_list.insert(make_pair("logy",&p_logy));

    p_font.i=&font;
    p_font.help="Font scaling for text objects (default 16).";
    cl->par_list.insert(make_pair("font",&p_font));

    p_left_margin.d=&left_margin;
    p_left_margin.help="Size of left margin for a new canvas.";
    cl->par_list.insert(make_pair("left_margin",&p_left_margin));

    p_right_margin.d=&right_margin;
    p_right_margin.help="Size of right margin for a new canvas.";
    cl->par_list.insert(make_pair("right_margin",&p_right_margin));

    p_bottom_margin.d=&bottom_margin;
    p_bottom_margin.help="Size of bottom margin for a new canvas.";
    cl->par_list.insert(make_pair("bottom_margin",&p_bottom_margin));

    p_top_margin.d=&top_margin;
    p_top_margin.help="Size of top margin for a new canvas.";
    cl->par_list.insert(make_pair("top_margin",&p_top_margin));

  }    
  
  return 0;
}

int acol_manager::run(int argc, char *argv[], bool set_o2graph_mode) {

  o2graph_mode=set_o2graph_mode;
  
  //--------------------------------------------------------------------
  // Default to scientific mode

  cout.setf(ios::scientific);

  // 
  if (verbose>2) {
    cout << "Setup cli class." << endl;
  }
  setup_cli();

  // 
  if (verbose>2) {
    cout << "Setup options." << endl;
  }
  setup_options();

  // 
  if (verbose>2) {
    cout << "Setup help." << endl;
  }
  setup_help();

  //-------------------------------------------------------------------

  cl->set_function(cset);

  static const int both=cli::comm_option_both;
  cl->remove_comm_option("help");
  cl->remove_comm_option("commands");
  static const size_t narr2=2;
  comm_option_s options_arr2[narr2]={
    {'h',"help","Show help information.",0,2,"[type] <command>",
     ((std::string)"Show generic help information, or, if an ")+
     "argument is given "+
     "give the documentation for the specified command. "+
     "Note that required arguments are typically given inside "+
     "angled brackes <> while optional arguments are given "+
     "inside square brackets [].",
     new o2scl::comm_option_mfptr<acol_manager>
     (this,&acol_manager::comm_help),both},
    {0,"commands","List commands",0,1,"[type]","",
     new o2scl::comm_option_mfptr<acol_manager>
     (this,&acol_manager::comm_commands),both}
  };
  cl->set_comm_option_vec(narr2,options_arr2);

  //-------------------------------------------------------------------
  // Try to get screen width
  
  int ncol=80;
  char *ncstring=getenv("COLUMNS");
  if (ncstring) {
    int nc2;
    int sret=o2scl::stoi_nothrow(ncstring,nc2);
    if (sret==0 && nc2>0) {
      ncol=nc2;
    } else {
      cerr << "Failed to interpret COLUMNS value " << ncstring
	   << " as a positive number of columns." << endl;
    }
  }
  
  set_swidth(ncol);

  //-------------------------------------------------------------------
  // Setup parameters modified by 'set' and 'get'

  if (verbose>2) {
    cout << "Setup parameters: " << endl;
  }
  setup_parameters();

  //-------------------------------------------------------------------
  // Process default options and call

  if (verbose>2) {
    cout << "Process default options" << endl;
  }
  std::vector<cmd_line_arg> ca;
  
  char *dc=getenv(env_var_name.c_str());
  if (dc) {
    def_args=dc;
    if (verbose>2) {
      cl->process_args(def_args,ca,1,true);
    } else {
      cl->process_args(def_args,ca,0,true);
    }
  }
  
  if (!o2graph_mode) {
    
    //----------------------------------------------------------------
    // Process command-line options
    
    // Note that it's ok that this appears early in the code because it
    // just processes the arguments, it doesn't do any execution based
    // on those arguments until later.
    
    if (verbose>2) {
      cout << "Process command-line options" << endl;
      cl->process_args(argc,argv,ca,1,true);
    } else {
      cl->process_args(argc,argv,ca,0,true);
    }
    if (argc<2) {
      post_interactive=true;
    }
    
    //------------------------------------------------------------------
    // Post interactive
    
    int ret2=0;
    
    if (verbose>2) {
      cout << "Post_interactive: " << post_interactive << endl;
    }
    
    if (post_interactive) {
      if (verbose>2) {
	cout << "Run interactive mode." << endl;
      }
      ret2=cl->run_interactive();
    }
    
    //--------------------------------------------------------------------
    // Notify user if error occurred
    
    if (ret2!=0) {
      cout << "An error occured." << endl;
    }
    
    delete cl;
  
  }

  return 0;

}

#ifdef O2SCL_NEVER_DEFINED
int acol_manager::comm_comment(std::vector<std::string> &sv, 
			       bool itive_com) {

  if (sv.size()==2) {
    hdf_file hf;
    hf.open(sv[1]);
    std::string def, s;
    hf.gets_def("comment",def,s);
    if (s==def) {
      cout << "No comment in file " << sv[1] << endl;
    } else {
      cout << "Comment in file " << sv[1] << " :" << endl;
      cout << s << endl;
    }
    hf.close();
    return 0;
  }
  
  hdf_file hf;
  // Make sure to open with write access
  hf.open(sv[1],1);
  
  // If it's already present as a fixed length string,
  // then we need to double check
  std::string def, s;
  int iret=hf.gets_def_fixed("comment",def,s);
  if (s!=def) {
    if (iret==1) {
      size_t len=s.length();
      if (sv[2].length()>len) {
	cerr << "Size of new comment (" << sv[2].length()
	     << ") longer than size of current "
	     << "fixed length string " << len << "." << endl;
	hf.close();
	return 1;
      } else {
	while (sv[2].length()<len) sv[2]+=' ';
      }
      hf.sets_fixed("comment",sv[2]);
    } else {
      hf.sets("comment",sv[2]);
    }
  } else {
    // String is not present so just set
    hf.sets("comment",sv[2]);
  }
  cout << "Set comment in file " << sv[1] << " to " << endl;
  cout << sv[2] << endl;
  hf.close();
  return 0;
}
#endif

int acol_manager::comm_to_hist(std::vector<std::string> &sv, 
			       bool itive_com) {

  std::string i1;

  if (type=="table") {

    bool twod_mode=false;
    
    if (sv.size()>=2 && sv[1]=="2d") {
      twod_mode=true;
      std::vector<std::string>::iterator it=sv.begin();
      it++;
      sv.erase(it);
    }
    if (sv.size()<2 && itive_com) {
      int ret=get_input_one(sv,((string)"Enter \"2d\" for 2d histogram ")+
			    +"and \"1d\" for 1d histogram",i1,"to-hist",
			    itive_com);
      if (ret!=0) return ret;
      if (i1=="2d") twod_mode=true;
    }
    
    if (twod_mode==false) {
      
      vector<string> in, pr;
      pr.push_back("Column name");
      pr.push_back("Number of bins");
      int ret=get_input(sv,pr,in,"to-hist",itive_com);
      if (ret!=0) return ret;
      
      std::string col2;
      if (sv.size()>3) {
	col2=sv[3];
      } else if (itive_com) {
	col2=cl->cli_gets("Column for weights (or blank for none): ");
      }

      size_t nbins;
      int sret=o2scl::stoszt_nothrow(in[1],nbins);
      if (sret!=0 || nbins==0) {
	cerr << "Failed to interpret " << in[1]
	     << " as a positive number of bins." << endl;
	return exc_einval;
      }
      if (col2.length()==0) {
	hist_obj.from_table(table_obj,in[0],nbins);
      } else {
	hist_obj.from_table(table_obj,in[0],col2,nbins);
      }
      command_add("hist");
      type="hist";

    } else {

      vector<string> in, pr;
      pr.push_back("Column name for x-axis");
      pr.push_back("Column name for y-axis");
      pr.push_back("Number of bins in x direction");
      pr.push_back("Number of bins in y direction");
      int ret=get_input(sv,pr,in,"to-hist-2d",itive_com);
      if (ret!=0) return ret;
      
      std::string col2;
      if (sv.size()>5) {
	col2=sv[5];
	// We don't want to prompt for weights if the user has given
	// enough arguments to proceed, so we test for sv.size()<5
      } else if (itive_com && sv.size()<5) {
	col2=cl->cli_gets("Column for weights (or blank for none): ");
      }

      size_t nbinsx, nbinsy;
      int sret=o2scl::stoszt_nothrow(in[2],nbinsx);
      if (sret!=0 || nbinsx==0) {
	cerr << "Failed to interpret " << in[2]
	     << " as a positive number of bins." << endl;
	return exc_einval;
      }
      sret=o2scl::stoszt_nothrow(in[3],nbinsy);
      if (sret!=0 || nbinsy==0) {
	cerr << "Failed to interpret " << in[2]
	     << " as a positive number of bins." << endl;
	return exc_einval;
      }

      if (col2.length()==0) {
	hist_2d_obj.from_table(table_obj,in[0],in[1],nbinsx,nbinsy);
      } else {
	hist_2d_obj.from_table(table_obj,in[0],in[1],col2,
			       nbinsx,nbinsy);
      }
      command_add("hist_2d");
      type="hist_2d";
      
    }

    return 0;
  } 

  cerr << "Cannot convert object of type " << type << " to histogram."
       << endl;
  
  return 1;
}

int acol_manager::comm_type(std::vector<std::string> &sv, 
			    bool itive_com) {
  if (type.length()==0) {
    cerr << "No current object to display type of." << endl;
    return 1;
  }
  cout << "Type is " << type << endl;
  return 0;
}

int acol_manager::comm_nlines(std::vector<std::string> &sv, 
			      bool itive_com) {
  if (type!="table") {
    cerr << "No table in 'nlines'." << endl;
    return 1;
  }

  if (table_obj.is_constant("nlines")) {
    cerr << "Constant 'nlines' already exists." << endl;
    return 2;
  }

  table_obj.add_constant("nlines",table_obj.get_nlines());
  
  return 0;
}

/*
  int acol_manager::comm_interp_type(std::vector<std::string> &sv, 
  bool itive_com) {

  if (type=="table3d") {
    
  if (type!="table3d") {
  cout << "No table to get interpolation type of in 'interp-type'." << endl;
  return exc_efailed;
  }

  if (sv.size()>1) {
  if (o2scl::stoi(sv[1])>7 || o2scl::stoi(sv[1])<0) {
  cout << "Invalid interpolation type in 'interp-type'." << endl;
  return exc_efailed;
  }
  table3d_obj.set_interp_type(o2scl::stoi(sv[1]));
  }
    
  if (sv.size()==1 || verbose>0) {
  size_t itype=table3d_obj.get_interp_type();
  cout << "Current interpolation type is " << itype;
      
  if (itype<4) {
  if (itype<3) {
  if (itype==1) {
  cout << " (linear)." << endl;
  } else {
  cout << " (cubic spline)." << endl;
  }
  } else {
  cout << " (cubic spline, periodic)." << endl;
  }
  } else {
  if (itype<6) {
  if (itype<5) {
  cout << " (Akima spline)." << endl;
  } else {
  cout << " (Akima spline, periodic)." << endl;
  }
  } else {
  if (itype<7) {
  cout << " (Monotonicity-preserving)." << endl;
  } else {
  cout << " (Steffen's monotonic)." << endl;
  }
  }
  }
  }
    
  return 0;
  }
  
  if (table_obj.get_nlines()==0) {
  cout << "No table to get interpolation type of." << endl;
  return exc_efailed;
  }

  if (sv.size()>1) {
  if (o2scl::stoi(sv[1])>7 || o2scl::stoi(sv[1])<0) {
  cout << "Invalid interpolation type in interp-type." << endl;
  return exc_efailed;
  }
  table_obj.set_interp_type(o2scl::stoi(sv[1]));
  }

  if (sv.size()==1 || verbose>0) {
  size_t itype=table_obj.get_interp_type();
  cout << "Current interpolation type is " << itype;
    
  if (itype<4) {
  if (itype<3) {
  if (itype==1) {
  cout << " (linear)." << endl;
  } else {
  cout << " (cubic spline)." << endl;
  }
  } else {
  cout << " (cubic spline, periodic)." << endl;
  }
  } else {
  if (itype<6) {
  if (itype<5) {
  cout << " (Akima spline)." << endl;
  } else {
  cout << " (Akima spline, periodic)." << endl;
  }
  } else {
  if (itype<7) {
  cout << " (Monotonicity-preserving)." << endl;
  } else {
  cout << " (Steffen's monotonic)." << endl;
  }
  }
  }
  }
  
  return 0;
  }
*/

int acol_manager::comm_output(std::vector<std::string> &sv, bool itive_com) {

  //--------------------------------------------------------------------
  // Output formatting
  
  ostream *fout;
  ofstream ffout;
  
  if (sv.size()==1) {
    fout=&cout;
  } else {
    ffout.open(sv[1].c_str());
    fout=&ffout;
  }
  
  if (scientific) fout->setf(ios::scientific);
  else fout->unsetf(ios::scientific);
  fout->precision(prec);

  if (type=="table3d") {
    
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    if (nx!=0 && ny!=0) {

      (*fout) << "Grid x: ";
      for(size_t i=0;i<((size_t)nx);i++) {
	(*fout) << table3d_obj.get_grid_x(i) << " ";
      }
      (*fout) << endl;
      (*fout) << "Grid y: ";
      for(size_t i=0;i<((size_t)ny);i++) {
	(*fout) << table3d_obj.get_grid_y(i) << " ";
      }
      (*fout) << endl;
      
      size_t nt=table3d_obj.get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  (*fout) << "Slice " << k << ": ";
	  (*fout) << table3d_obj.get_slice_name(k) << endl;
	  if (k==0) {
	    (*fout) << "Outer loops over x grid, inner loop over y grid." 
		    << endl;
	  }
	  for(size_t i=0;i<nx;i++) {
	    for(size_t j=0;j<ny;j++) {
	      (*fout) << table3d_obj.get(i,j,k) << " ";
	    }
	    (*fout) << endl;
	  }
	  fout->unsetf(ios::showpos);
	}
      }
    }

    return 0;

  } else if (type=="table") {
    
    if (table_obj.get_ncolumns()>0) {

      //--------------------------------------------------------------------
      // Count column widths

      vector<size_t> col_wids(table_obj.get_ncolumns());

      for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	col_wids[i]=prec+6;
      }

      if (names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	  if (table_obj.get_column_name(i).size()>col_wids[i]) {
	    col_wids[i]=table_obj.get_column_name(i).size();
	  }
	  std::string tunit=table_obj.get_unit(table_obj.get_column_name(i));
	  if (tunit.size()+2>col_wids[i]) {
	    col_wids[i]=tunit.size()+2;
	  }
	}
      }

      //--------------------------------------------------------------------
      // Output column names
      
      if (names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	  
	  // Preceeding space
	  if (pretty==true) {
	    (*fout) << ' ';
	  }
	  
	  // Column name
	  (*fout) << table_obj.get_column_name(i) << " ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=col_wids[i]-table_obj.get_column_name(i).size();
	    if (nsp<0) {
	      O2SCL_ERR("Column size anomaly (col names).",exc_efailed);
	    }
	    for(int j=0;j<nsp;j++) (*fout) << ' ';
	  }
	
	}
	(*fout) << endl;
      }
      
      //--------------------------------------------------------------------
      // Output units
    
      if (table_obj.get_nunits()>0 && names_out==true) {
	for(size_t i=0;i<table_obj.get_ncolumns();i++) {

	  // Preceeding space
	  if (pretty==true) {
	    (*fout) << ' ';
	  }
	
	  // Unit name
	  (*fout) << '['
		  << table_obj.get_unit(table_obj.get_column_name(i)) << "] ";
	
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=col_wids[i]-
	      table_obj.get_unit(table_obj.get_column_name(i)).size()-2;
	    if (nsp<0) {
	      O2SCL_ERR("Column size anomaly (units).",exc_efailed);
	    }
	    for(int j=0;j<nsp;j++) (*fout) << ' ';
	  }
	
	}
	(*fout) << endl;
      }
      
      //--------------------------------------------------------------------
      // Output data
      
      for(int i=0;i<((int)table_obj.get_nlines());i++) {
	
	for(size_t j=0;j<table_obj.get_ncolumns();j++) {
	  
	  // Otherwise, for normal output
	  if (pretty==true && table_obj.get(j,i)>=0.0) {
	    (*fout) << ' ';
	  }
	  (*fout) << table_obj.get(j,i) << ' ';
	  if (pretty==true) {
	    int nsp=((int)(table_obj.get_column_name(j).size()-prec-6));
	    for(int kk=0;kk<nsp;kk++) {
	      (*fout) << ' ';
	    }
	  }
	
	}
	(*fout) << endl;
      
	//--------------------------------------------------------------------
	// Continue to next row
      }
    
    }

  } else if (type=="hist") {
    
    for(size_t k=0;k<hist_obj.size();k++) {
      (*fout) << hist_obj.get_bin_low_i(k) << " ";
    }
    (*fout) << hist_obj.get_bin_high_i(hist_obj.size()-1) << endl;
    for(size_t k=0;k<hist_obj.size();k++) {
      (*fout) << hist_obj.get_wgt_i(k) << endl;
    }
    
  } else if (type=="int") {
    
    (*fout) << int_obj << endl;
    
  } else if (type=="char") {
    
    (*fout) << char_obj << endl;
    
  } else if (type=="double") {
    
    (*fout) << double_obj << endl;
    
   } else if (type=="size_t") {

    (*fout) << size_t_obj << endl;

  } else if (type=="string") {

    (*fout) << string_obj << endl;

  } else if (type=="int[]") {

    for(size_t k=0;k<intv_obj.size();k++) {
      (*fout) << intv_obj[k] << " ";
    }
    (*fout) << endl;

  } else if (type=="double[]") {

    for(size_t k=0;k<doublev_obj.size();k++) {
      (*fout) << doublev_obj[k] << " ";
    }
    (*fout) << endl;

  } else if (type=="size_t[]") {

    for(size_t k=0;k<size_tv_obj.size();k++) {
      (*fout) << size_tv_obj[k] << " ";
    }
    (*fout) << endl;

  } else if (type=="string[]") {

    (*fout) << stringv_obj.size() << endl;
    for(size_t k=0;k<stringv_obj.size();k++) {
      (*fout) << stringv_obj[k] << endl;
    }
    (*fout) << endl;

  } else if (type=="hist_2d") {
    
    for(size_t k=0;k<hist_2d_obj.size_x();k++) {
      (*fout) << hist_2d_obj.get_x_low_i(k) << " ";
    }
    (*fout) << hist_2d_obj.get_x_high_i(hist_2d_obj.size_x()-1) << endl;

    for(size_t k=0;k<hist_2d_obj.size_y();k++) {
      (*fout) << hist_2d_obj.get_y_low_i(k) << " ";
    }
    (*fout) << hist_2d_obj.get_y_high_i(hist_2d_obj.size_y()-1) << endl;

    for(size_t ki=0;ki<hist_2d_obj.size_x();ki++) {
      for(size_t kj=0;kj<hist_2d_obj.size_y();kj++) {
	(*fout) << hist_2d_obj.get_wgt_i(ki,kj) << " ";
      }
      (*fout) << endl;
    }

  } else if (type=="vector<contour_line>") {

    (*fout) << cont_obj.size() << endl;
    for(size_t k=0;k<cont_obj.size();k++) {
      (*fout) << cont_obj[k].level << " " << cont_obj[k].x.size() << endl;
      for(size_t kk=0;kk<cont_obj[k].x.size();kk++) {
	(*fout) << cont_obj[k].x[kk] << " ";
	(*fout) << cont_obj[k].y[kk] << endl;
      }
    }

  } else if (type=="uniform_grid<double>") {

    (*fout) << ug_obj.get_nbins() << " ";
    (*fout) << ug_obj.get_start() << " ";
    (*fout) << ug_obj.get_end() << " ";
    (*fout) << ug_obj.get_width() << endl;

  } else {

    cerr << "Cannot output type " << type << endl;
    return 2;
    
  }

  if (sv.size()!=1) {
    ffout.close();
  }
  
  return 0;
}

int acol_manager::comm_cat(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Not enough arguments to cat." << endl;
    return exc_efailed;
  }
  string file2=sv[1];
  
  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d to add to in command 'cat'." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Read the new table3d
    
    table3d tab2;

    hdf_file hf;
    std::string name2;
    if (sv.size()>=3) name2=sv[2];

    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    // ---------------------------------------------------------------------
    // Copy constants from the new table3d

    for(size_t i=0;i<tab2.get_nconsts();i++) {
      string tnam;
      double tval;
      tab2.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      table3d_obj.add_constant(tnam,tval);
    }

    // ---------------------------------------------------------------------
    // Copy slices over if not already present in the current table3d

    const ubvector &xg=tab2.get_x_data();
    const ubvector &yg=tab2.get_y_data();

    for(size_t k=0;k<tab2.get_nslices();k++) {
      std::string sl_name=tab2.get_slice_name(k);
      size_t slix;
      if (!table3d_obj.is_slice(sl_name,slix)) {
	table3d_obj.new_slice(sl_name);
	for(size_t i=0;i<tab2.get_nx();i++) {
	  for(size_t j=0;j<tab2.get_ny();j++) {
	    double x=xg[i];
	    double y=yg[j];
	    table3d_obj.set_val(x,y,sl_name,tab2.get(i,j,sl_name));
	  }
	}
      }
    }

  } else if (type=="table") {

    if (table_obj.get_nlines()==0) {
      cerr << "No table to add to in command 'cat'." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Read the new table 

    table_units<> tab2;

    hdf_file hf;
    std::string name2;
    if (sv.size()>=3) name2=sv[2];

    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    // ---------------------------------------------------------------------
    // Copy constants from the new table

    for(size_t i=0;i<tab2.get_nconsts();i++) {
      string tnam;
      double tval;
      tab2.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      table_obj.add_constant(tnam,tval);
    }

    // ---------------------------------------------------------------------

    size_t n1=table_obj.get_nlines();
    size_t n2=tab2.get_nlines();
    table_obj.set_nlines(n1+n2);
    for(size_t j=0;j<tab2.get_ncolumns();j++) {
      std::string col_name=tab2.get_column_name(j);
      if (!table_obj.is_column(col_name)) {
	table_obj.new_column(col_name);
	for(size_t i=0;i<n1+n2;i++) table_obj.set(col_name,i,0.0);
      }
      for(size_t i=0;i<n2;i++) {
	table_obj.set(col_name,i+n1,tab2.get(col_name,i));
      }
    }
    
  } else {

    cerr << "Cannot 'cat' with object of type " << type << endl;
    return exc_efailed;
    
  }
  
  return 0;
}

int acol_manager::comm_sum(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    if (sv.size()<2) {
      cerr << "Not enough arguments to sum." << endl;
      return exc_efailed;
    }

    string s2=sv[1];
    string name2;
    if (sv.size()>=3) name2=sv[2];

    table3d t2;

    hdf_file hf;
    int hfret=hf.open(s2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << s2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,t2,name2);
    hf.close();
  
    size_t nx=t2.get_nx();
    size_t ny=t2.get_ny();
    const ubvector &xg=t2.get_x_data();
    const ubvector &yg=t2.get_y_data();
  
    for(size_t k=0;k<t2.get_nslices();k++) {
      string slname=t2.get_slice_name(k);
      size_t slix;
      if (table3d_obj.is_slice(slname,slix)==false) {
	table3d_obj.new_slice(slname);
	table3d_obj.set_slice_all(slname,0.0);
      }
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  double x=xg[i];
	  double y=yg[j];
	  double value=table3d_obj.get_val(x,y,slname)+t2.get(i,j,slname);
	  table3d_obj.set_val(x,y,slname,value);
	}
      }
    }

  } else if (type=="table") {

    if (sv.size()<2) {
      cerr << "Not enough arguments to add." << endl;
      return exc_efailed;
    }
    if (table_obj.get_nlines()==0) {
      cerr << "No table to add to." << endl;
      return exc_efailed;
    }

    string file2=sv[1];
    std::string name2;
    if (sv.size()>=3) name2=sv[2];
    table_units<> tab2;

    hdf_file hf;
    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    size_t n1=table_obj.get_nlines();
    size_t n2=tab2.get_nlines();
    if (n2>n1) {
      table_obj.set_nlines(n1+n2);
      for(size_t j=0;j<table_obj.get_ncolumns();j++) {
	for(size_t i=n1;i<n1+n2;i++) {
	  table_obj.set(j,i,0.0);
	}
      }
    }
    for(size_t j=0;j<tab2.get_ncolumns();j++) {
      std::string col_name=tab2.get_column_name(j);
      if (!table_obj.is_column(col_name)) {
	table_obj.new_column(col_name);
	for(size_t i=0;i<table_obj.get_nlines();i++) {
	  table_obj.set(col_name,i,0.0);
	}
      }
      for(size_t i=0;i<n2;i++) {
	table_obj.set(col_name,i,tab2.get(col_name,i)+
		      table_obj.get(col_name,i));
      }
    }
    
  } else {

    cerr << "Cannot 'sum' with object of type " << type << endl;
    return exc_efailed;
    
  }
  
  return 0;
}

herr_t acol_manager::iterate_func(hid_t loc, const char *name, 
				  const H5L_info_t *inf, void *op_data) {

  // Arrange parameters
  iter_parms *ip=(iter_parms *)op_data;
  string tname=ip->tname;
  hdf_file &hf=*(ip->hf);
  int verbose=ip->verbose;

  hid_t top=hf.get_current_id();

  H5O_info_t infobuf;
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5P_DEFAULT);
  
  // If it's a group
  if (infobuf.type==H5O_TYPE_GROUP) {

    // Open the group and see if it's an O2scl object
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    string otype;
    hf.gets_def_fixed("o2scl_type","",otype);
    hf.close_group(group);
    hf.set_current_id(top);

    if (otype.length()!=0 && 
	(otype=="table" || otype=="table3d")) {
      
      if (verbose>1) {
	cout << "Found O2scl group named: " << name << " of type: " 
	     << otype << endl;
      }

      if (tname.size()>0 && tname==((string)name)) {
	if (verbose>0) {
	  cout << "Found table named '" << name << "'." << endl;
	}
	ip->found=true;
	ip->type=otype;
	return 1;
      }
      if (tname.size()==0) {
	if (verbose>0) {
	  cout << "Found table named '" << name << "'." << endl;
	}
	ip->tname=name;
	ip->found=true;
	ip->type=otype;
	return 1;
      }
      
    } else {
      if (verbose>1) {
	cout << "Non o2scl group named: " << name << endl;
      }
    }

  } else if (infobuf.type==H5O_TYPE_DATASET) {
    if (verbose>1) {
      cout << "Dataset: " << name << endl;
    }
  } else if (infobuf.type==H5O_TYPE_NAMED_DATATYPE) {
    if (verbose>1) {
      cout << "Named type: " << name << endl;
    }
  } else {
    cout << "Unexpected HDF type. " << endl;
  }
  return 0;
}

herr_t acol_manager::filelist_func(hid_t loc, const char *name, 
				   const H5L_info_t *inf, void *op_data) {

  // Arrange parameters
  iter_parms *ip=(iter_parms *)op_data;
  string tname=ip->tname;
  hdf_file &hf=*(ip->hf);
  int loc_verbose=ip->verbose;
  int mode=ip->mode;

  hid_t top=hf.get_current_id();

  H5O_info_t infobuf;
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5P_DEFAULT);

  ip->found=false;
  
  // If it's a group
  if (infobuf.type==H5O_TYPE_GROUP) {

    // Open the group and see if it's an O2scl object
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    string otype;
    hf.gets_def_fixed("o2scl_type","",otype);
    hf.close_group(group);
    hf.set_current_id(top);

    if (otype.length()!=0) {
      if (mode==0) {
	if (otype==((string)"string[]")) {
	  cout << "O2scl object \"" << name << "\" of type " 
	       << otype << "." << endl;
	} else {
	  cout << "O2scl object \"" << name << "\" of type " 
	       << otype << "." << endl;
	}
      }
      if (mode==1 && name==tname) {
	ip->type=otype;
	ip->found=true;
	return 1;
      }
    } else {
      if (mode==0) {
	cout << "Group \"" << name << "\"." << endl;
      }
    }
    
  } else if (infobuf.type==H5O_TYPE_DATASET) {
    if (mode==0) {
      cout << "Dataset \"" << name << "\" of type ";
    }

    // Open data set
    hid_t dset=H5Dopen(loc,name,H5P_DEFAULT);

    // Get type information
    hid_t type_id=H5Dget_type(dset);
    hid_t nat_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);

    // Get filter information
    int compr=0;
    hid_t plist_id=H5Dget_create_plist(dset);
    int num_filters=H5Pget_nfilters(plist_id);
    for(int i=0;i<num_filters;i++) {
      size_t n_elements=0;
      unsigned flags, filter_info;
      H5Z_filter_t filter_type=H5Pget_filter2
	(plist_id,0,&flags,&n_elements,0,0,0,&filter_info);
      if (filter_type==H5Z_FILTER_DEFLATE) {
	compr=1;
      } else if (filter_type==H5Z_FILTER_SZIP) {
	compr=2;
      }
    }

    if (mode==0) {
      if (compr==1) {
	cout << "(1) ";
      } else if (compr==2) {
	cout << "(2) ";
      }
    }
    
    // Get dataspace information
    hid_t space_id = H5Dget_space(dset);
    hsize_t dims[100];
    hsize_t max_dims[100];
    int ndims=H5Sget_simple_extent_dims(space_id,dims,max_dims);
    
    if (H5Tequal(nat_id,H5T_NATIVE_CHAR)) {
      if (ndims==1) {
	if (max_dims[0]==H5S_UNLIMITED) {
	  if (mode==0) {
	    cout << "char[" << dims[0] << "/inf] with value=\"";
	  }
	  if (mode==1 && name==tname) {
	    ip->type="string";
	    ip->found=true;
	    return 1;
	  }
	} else {
	  if (mode==0) {
	    cout << "char[" << dims[0] << "/" << max_dims[0]
		 << "] with value=\"";
	  }
	  if (mode==1 && name==tname) {
	    if (dims[0]==1) {
	      ip->type="char";
	      ip->found=true;
	      return 1;
	    } else {
	      ip->type="string(fixed)";
	      ip->found=true;
	      return 1;
	    }
	  }
	}
	if (mode==0) {
	  std::string s;
	  hf.gets(name,s);
	  if (dims[0]<20) {
	    cout << s;
	  } else {
	    cout << s[0] << s[1] << s[2] << s[3] << s[4] << s[5]
		 << " ... "
		 << s[dims[0]-6] << s[dims[0]-5] << s[dims[0]-4]
		 << s[dims[0]-3] << s[dims[0]-2] << s[dims[0]-1];
	  }
	  cout << "\".";
	}
      } else {
	if (mode==0) {
	  cout << "char[";
	  for(int i=0;i<ndims-1;i++) {
	    if (max_dims[i]==H5S_UNLIMITED) {
	      cout << dims[i] << "/inf,";
	    } else {
	      cout << dims[i] << "/" << max_dims[i] << ",";
	    }
	  }
	  if (max_dims[ndims-1]==H5S_UNLIMITED) {
	    cout << dims[ndims-1] << "/inf].";
	  } else {
	    cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	  }
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_SHORT)) {
      if (mode==0) {
	cout << "short[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_INT)) {
      if (ndims==1 && dims[0]>0) {
	if (max_dims[0]==H5S_UNLIMITED) {
	  if (mode==0) {
	    cout << "int[" << dims[0] << "/inf] with value=";
	  }
	  if (mode==1 && name==tname) {
	    ip->type="int[]";
	    ip->found=true;
	    return 1;
	  }
	} else {
	  if (mode==0) {
	    cout << "int[" << dims[0] << "/" << max_dims[0]
		 << "] with value=";
	  }
	  if (mode==1 && name==tname) {
	    if (dims[0]==1) {
	      ip->type="int";
	      ip->found=true;
	      return 1;
	    } else {
	      ip->type="int[fixed]";
	      ip->found=true;
	      return 1;
	    }
	  }
	}
	if (mode==0) {
	  std::vector<int> iarr;
	  hf.geti_vec(name,iarr);
	  if (dims[0]==1) {
	    cout << iarr[0];
	  } else if (dims[0]==2) {
	    cout << iarr[0] << ", " << iarr[1];
	  } else if (dims[0]==3) {
	    cout << iarr[0] << ", " << iarr[1] << ", " << iarr[2];
	  } else {
	    cout << iarr[0] << ", " << iarr[1] << ", ..., " << iarr[dims[0]-1];
	  }
	  cout << ".";
	}
      } else {
	if (mode==0) {
	  cout << "int[";
	  for(int i=0;i<ndims-1;i++) {
	    if (max_dims[i]==H5S_UNLIMITED) {
	      cout << dims[i] << "/inf,";
	    } else {
	      cout << dims[i] << "/" << max_dims[i] << ",";
	    }
	  }
	  if (max_dims[ndims-1]==H5S_UNLIMITED) {
	    cout << dims[ndims-1] << "/inf].";
	  } else {
	    cout << dims[ndims-1] << "/" << max_dims[ndims-1] << ").";
	  }
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_LONG)) {
      if (mode==0) {
	cout << "long[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_LLONG)) {
      if (mode==0) {
	cout << "llong[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_UCHAR)) {
      if (mode==0) {
	cout << "uchar[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_USHORT)) {
      if (mode==0) {
	cout << "ushort[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_UINT)) {
      if (mode==0) {
	cout << "uint[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULONG)) {
      if (ndims==1 && dims[0]>0) {
	if (max_dims[0]==H5S_UNLIMITED) {
	  if (mode==0) {
	    cout << "size_t[" << dims[0] << "/inf] with value=";
	  }
	  if (mode==1 && name==tname) {
	    ip->type="size_t[]";
	    ip->found=true;
	    return 1;
	  }
	} else {
	  if (mode==0) {
	    cout << "size_t[" << dims[0] << "/"
		 << max_dims[0] << "] with value=";
	  }
	  if (mode==1 && name==tname) {
	    if (dims[0]==1) {
	      ip->type="size_t";
	      ip->found=true;
	      return 1;
	    } else {
	      ip->type="size_t[fixed]";
	      ip->found=true;
	      return 1;
	    }
	  }
	}
	if (mode==0) {
	  std::vector<size_t> sarr;
	  hf.get_szt_vec(name,sarr);
	  if (dims[0]==1) {
	    cout << sarr[0];
	  } else if (dims[0]==2) {
	    cout << sarr[0] << ", " << sarr[1];
	  } else if (dims[0]==3) {
	    cout << sarr[0] << ", " << sarr[1] << ", " << sarr[2];
	  } else {
	    cout << sarr[0] << ", " << sarr[1] << ", ..., " << sarr[dims[0]-1];
	  }
	  cout << ".";
	}
      } else {
	if (mode==0) {
	  cout << "size_t[";
	  for(int i=0;i<ndims-1;i++) {
	    if (max_dims[i]==H5S_UNLIMITED) {
	      cout << dims[i] << "/inf,";
	    } else {
	      cout << dims[i] << "/" << max_dims[i] << ",";
	    }
	  }
	  if (max_dims[ndims-1]==H5S_UNLIMITED) {
	    cout << dims[ndims-1] << "/inf].";
	  } else {
	    cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	  }
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULLONG)) {
      if (mode==0) {
	cout << "ullong[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_FLOAT)) {
      if (mode==0) {
	cout << "float[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf],";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_DOUBLE)) {
      if (ndims==1 && dims[0]>0) {
	if (max_dims[0]==H5S_UNLIMITED) {
	  if (mode==0) {
	    cout << "double[" << dims[0] << "/inf] with value=";
	  }
	  if (mode==1 && name==tname) {
	    ip->type="double[]";
	    ip->found=true;
	    return 1;
	  }
	} else {
	  if (mode==0) {
	    cout << "double[" << dims[0] << "/"
		 << max_dims[0] << "] with value=";
	  }
	  if (mode==1 && name==tname) {
	    if (dims[0]==1) {
	      ip->type="double";
	      ip->found=true;
	      return 1;
	    } else {
	      ip->type="double[fixed]";
	      ip->found=true;
	      return 1;
	    }
	  }
	}
	if (mode==0) {
	  std::vector<double> darr;
	  hf.getd_vec(name,darr);
	  if (dims[0]==1) {
	    cout << darr[0];
	  } else if (dims[0]==2) {
	    cout << darr[0] << ", " << darr[1];
	  } else {
	    cout << "\n\t" << darr[0] << ", ..., " << darr[dims[0]-1];
	  }
	  cout << ".";
	}
      } else {
	if (mode==0) {
	  cout << "double[";
	  for(int i=0;i<ndims-1;i++) {
	    if (max_dims[i]==H5S_UNLIMITED) {
	      cout << dims[i] << "/inf,";
	    } else {
	      cout << dims[i] << "/" << max_dims[i] << ",";
	    }
	  }
	  if (max_dims[ndims-1]==H5S_UNLIMITED) {
	    cout << dims[ndims-1] << "/inf].";
	  } else {
	    cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	  }
	}
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_LDOUBLE)) {
      if (mode==0) {
	cout << "ldouble[";
	for(int i=0;i<ndims-1;i++) {
	  if (max_dims[i]==H5S_UNLIMITED) {
	    cout << dims[i] << "/inf]";
	  } else {
	    cout << dims[i] << "/" << max_dims[i] << ",";
	  }
	}
	if (max_dims[ndims-1]==H5S_UNLIMITED) {
	  cout << dims[ndims-1] << "/inf].";
	} else {
	  cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "].";
	}
      }
    } else {
      hid_t filetype=H5Dget_type(dset);
      size_t str_size=H5Tget_size(filetype);
      
      hsize_t dims[3];
      hid_t space=H5Dget_space(dset);
      int ndims=H5Sget_simple_extent_dims(space,dims,0);
      hid_t memtype=-1;
      if (ndims==1 && dims[0]==1) {
	memtype=H5Tcopy(H5T_C_S1);
      }
      if (memtype>0) {
	status=H5Tclose(memtype);
	status=H5Sclose(space);
	status=H5Tclose(filetype);
	if (mode==1 && name==tname) {
	  ip->type="char[fixed]";
	  ip->found=true;
	  return 1;
	}
	if (mode==0) {
	  if (str_size==0) {
	    cout << "Error. Fixed length-string with no space for null.";
	  } else {
	    cout << "fixed-length (" << str_size-1 << ") string with value \"";
	    std::string s;
	    hf.gets_fixed(name,s);
	    cout << s << "\".";
	  }
	}
      } else {
	if (mode==0) {
	  cout << "<unknown>.";
	}
	status=H5Sclose(space);
	status=H5Tclose(filetype);
      }
    }
    if (mode==0) {
      cout << endl;
    }

    H5Sclose(space_id);
    H5Tclose(nat_id);
    H5Tclose(type_id);
    H5Dclose(dset);
    
  } else if (infobuf.type==H5O_TYPE_NAMED_DATATYPE) {
    if (mode==0) {
      cout << "Named type \"" << name << "\"." << endl;
    }
  } else {
    if (mode==0) {
      cout << "Unexpected HDF type. " << endl;
    }
  }
  return 0;
}

int acol_manager::comm_filelist(std::vector<std::string> &sv, 
				bool itive_com) {
  
  std::string i1;
  if (sv.size()==1) {
    if (itive_com) {
      i1=cl->cli_gets("Enter filename (or blank to quit): ");
      if (i1.length()==0) {
        if (verbose>0) cout << "Command 'filelist' cancelled." << endl;
        return 0;
      } 
    } else {
      cout << "Filename not given for command 'filelist'." << endl;
      return exc_efailed;
    }
  } else {
    i1=sv[1];
  }

  // Use hdf_file to open the file
  hdf_file hf;
  int hfret=hf.open(i1.c_str(),false,false);
  if (hfret!=0) {
    cerr << "Failed to read file named " << i1.c_str() << endl;
    return exc_efailed;
  }

  // Set the proper output precision and mode
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  
  iter_parms ip={"",&hf,false,"",verbose,0};

  H5Literate(hf.get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	     0,filelist_func,&ip);
  
  return 0;
}

int acol_manager::comm_read(std::vector<std::string> &sv, 
			    bool itive_com) {

  std::string i1, i2;
  if (sv.size()==1) {
    if (itive_com) {
      i1=cl->cli_gets("Enter filename (or blank to quit): ");
      if (i1.length()==0) {
        if (verbose>0) cout << "Command 'read' cancelled." << endl;
        return 0;
      } 
      i2=cl->cli_gets("Enter object name (or blank for first table): ");
    } else {
      cout << "Filename not given for command 'read'." << endl;
      return exc_efailed;
    }
  } else {
    i1=sv[1];
    if (sv.size()>2) {
      i2=sv[2];
    }
  }

  // Delete previous object
  command_del();
  clear_obj();

  // Use hdf_file to open the file
  hdf_file hf;
  string type2;
  int ret;

  ret=hf.open(i1.c_str(),false,false);
  if (ret!=0) {
    cerr << "Could not find file named '" << i1 << "'. Wrong file name?" 
	 << endl;
    return exc_efailed;
  }

  if (i2.length()!=0) {
    
    //cout << "read2 name: " << i2 << endl;
    
    iter_parms ip={i2,&hf,false,type,verbose,1};
    
    H5Literate(hf.get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	       0,filelist_func,&ip);
    
    //cout << "read2 type: " << ip.type << endl;

    if (ip.found==false) {
      cerr << "Could not find readable object named " << i2
	   << " in file " << i1 << endl;
      return 1;
    }
    
    if (ip.type=="table") {
      if (verbose>2) {
	cout << "Reading table." << endl;
      }
      hdf_input(hf,table_obj,i2);
      obj_name=i2;
      interp_type=table_obj.get_interp_type();
      command_add("table");
      type="table";
      return 0;
    } else if (ip.type=="table3d") {
      if (verbose>2) {
	cout << "Reading table3d." << endl;
      }
      hdf_input(hf,table3d_obj,i2);
      obj_name=i2;
      interp_type=table3d_obj.get_interp_type();
      command_add("table3d");
      type="table3d";
      return 0;
    } else if (ip.type=="hist") {
      if (verbose>2) {
	cout << "Reading hist." << endl;
      }
      hdf_input(hf,hist_obj,i2);
      obj_name=i2;
      command_add("hist");
      type="hist";
      return 0;
    } else if (ip.type=="hist_2d") {
      if (verbose>2) {
	cout << "Reading hist_2d." << endl;
      }
      hdf_input(hf,hist_2d_obj,i2);
      obj_name=i2;
      command_add("hist_2d");
      type="hist_2d";
      return 0;
    } else if (ip.type=="vector<contour_lines>") {
      if (verbose>2) {
	cout << "Reading vector<contour_lines>." << endl;
      }
      hdf_input(hf,cont_obj,i2);
      obj_name=i2;
      command_add("vector<contour_lines>");
      type="vector<contour_lines>";
      return 0;
    } else if (ip.type=="uniform_grid<double>") {
      if (verbose>2) {
	cout << "Reading uniform_grid<double>." << endl;
      }
      hdf_input(hf,ug_obj,i2);
      obj_name=i2;
      command_add("uniform_grid<double>");
      type="uniform_grid<double>";
      return 0;
    } else if (ip.type=="string[]") {
      if (verbose>2) {
	cout << "Reading string[]." << endl;
      }
      hf.gets_vec(i2,stringv_obj);
      obj_name=i2;
      command_add("string[]");
      type="string[]";
      return 0;
    } else if (ip.type=="int") {
      if (verbose>2) {
	cout << "Reading int." << endl;
      }
      hf.geti(i2,int_obj);
      obj_name=i2;
      command_add("int");
      type="int";
      return 0;
    } else if (ip.type=="char") {
      if (verbose>2) {
	cout << "Reading char." << endl;
      }
      hf.getc(i2,char_obj);
      obj_name=i2;
      command_add("char");
      type="char";
      return 0;
    } else if (ip.type=="string") {
      if (verbose>2) {
	cout << "Reading string." << endl;
      }
      hf.gets(i2,string_obj);
      obj_name=i2;
      command_add("string");
      type="string";
      return 0;
    } else if (ip.type=="double") {
      if (verbose>2) {
	cout << "Reading double." << endl;
      }
      hf.getd(i2,double_obj);
      obj_name=i2;
      command_add("double");
      type="double";
      return 0;
    } else if (ip.type=="size_t") {
      if (verbose>2) {
	cout << "Reading size_t." << endl;
      }
      hf.get_szt(i2,size_t_obj);
      obj_name=i2;
      command_add("size_t");
      type="size_t";
      return 0;
    } else if (ip.type=="int[]") {
      if (verbose>2) {
	cout << "Reading int[]." << endl;
      }
      hf.geti_vec(i2,intv_obj);
      obj_name=i2;
      command_add("int[]");
      type="int[]";
      return 0;
    } else if (ip.type=="double[]") {
      if (verbose>2) {
	cout << "Reading double[]." << endl;
      }
      hf.getd_vec(i2,doublev_obj);
      obj_name=i2;
      command_add("double[]");
      type="double[]";
      return 0;
    } else if (ip.type=="size_t[]") {
      if (verbose>2) {
	cout << "Reading size_t[]." << endl;
      }
      hf.get_szt_vec(i2,size_tv_obj);
      obj_name=i2;
      command_add("size_t[]");
      type="size_t[]";
      return 0;
    }

    cerr << "Found object with name " << i2
	 << " in file " << i1 << " but type " << ip.type
	 << " is not readable." << endl;
    return 2;
  }
  
  ret=hf.find_group_by_type("table",i2,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, reading first table object." << endl;
    }
    hdf_input(hf,table_obj,i2);
    obj_name=i2;
    command_add("table");
    type="table";
    interp_type=table_obj.get_interp_type();
    return 0;
  }
    
  ret=hf.find_group_by_type("table3d",i2,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, reading first table3d object." << endl;
    }
    hdf_input(hf,table3d_obj,i2);
    obj_name=i2;
    interp_type=table3d_obj.get_interp_type();
      
    command_add("table3d");
    type="table3d";
      
    return 0;
  }
    
  ret=hf.find_group_by_type("hist",i2,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, reading first hist object." << endl;
    }
    hdf_input(hf,hist_obj,i2);
    obj_name=i2;
    command_add("hist");
    type="hist";
    return 0;
  }
    
  ret=hf.find_group_by_type("hist_2d",i2,verbose);
  if (ret==success) {
    if (verbose>0) {
      cout << "No name specified, reading first hist_2d object." << endl;
    }
    hdf_input(hf,hist_2d_obj,i2);
    obj_name=i2;
    command_add("hist_2d");
    type="hist_2d";
    return 0;
  }
  
  cout << "Could not find object of any readable type in file '" << i1
       << "'." << endl;
  
  return exc_efailed;
}

int acol_manager::get_input_one(vector<string> &sv, string directions,
				string &in, string comm_name,
				bool itive_com) {

  // If there are enough arguments, then just fill 'in' with the
  // correct values from 'sv'
  if (sv.size()>1) {
    in=sv[1];
    return 0;
  }
  if (itive_com) {
    string temp=directions+" (or blank to stop): ";
    in=cl->cli_gets(temp.c_str());
    if (in.length()==0 || o2scl::count_words(in)==0) {
      if (verbose>0) {
	cout << "Command '" << comm_name << "' cancelled." << endl;
      }
      return exc_efailed;
    }
  } else {
    cerr << "Not enough arguments to '" << comm_name << "'." << endl;
    return exc_efailed;
  }

  return success;
}

int acol_manager::get_input(vector<string> &sv, vector<string> &directions,
			    vector<string> &in, string comm_name,
			    bool itive_com) {

  size_t ni=directions.size();

  // If there are enough arguments, then just fill the vector
  // 'in' with the correct values from 'sv'
  if (sv.size()>ni) {
    for(size_t i=0;i<ni;i++) {
      in.push_back(sv[i+1]);
    }
    return 0;
  }

  // Otherwise, if we're in interactive mode
  if (itive_com) {
    // Prompt the user for the correct arguments
    for(size_t i=0;i<ni;i++) {
      string temp=directions[i]+" (or blank to stop): ";
      in.push_back(cl->cli_gets(temp.c_str()));
      // If the user just pressed 'enter', then cancel
      if (in[i].length()==0 || count_words(in[i])==0) {
	cout << "Command '" << comm_name << "' cancelled." << endl;
	return exc_efailed;
      }
    }
  } else {
    // We don't have enough arguments and we're not in interactive
    // mode, so we fail
    cerr << "Not enough arguments to '" << comm_name << "'." << endl;
    return exc_efailed;
  }

  return success;
}

int acol_manager::comm_max(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (table3d_obj.get_nslices()==0) {
      cerr << "No slices to find the maximum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find maximum value of",
			  i1,"max",itive_com);
    if (ret!=0) return ret;

    size_t ix;
    if (!table3d_obj.is_slice(i1,ix)) {
      cerr << "No slice named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const ubmatrix &mat=table3d_obj.get_slice(ix);
    size_t i, j;
    double max;
    matrix_max_index(mat,i,j,max);

    cout << "Maximum value of slice '" << i1 << "' is: " 
	 << max << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="hist_2d") {
    
    const ubmatrix &mat=hist_2d_obj.get_wgts();
    size_t i, j;
    double max;
    matrix_max_index(mat,i,j,max);
    
    cout << "Maximum weight is: "
	 << max << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to find the maximum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column to find maximum value of",
			  i1,"max",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    double max;
    size_t ix;
    vector_max(table_obj.get_nlines(),(table_obj)[i1],ix,max);
    cout << "Maximum value of column '" << i1 << "' is: " 
	 << max << " at row with index " << ix << "." << endl;

  } else if (type=="double[]") {

    double val;
    size_t loc;
    o2scl::vector_max<vector<double>,double>(doublev_obj.size(),
					     doublev_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="int[]") {
    
    int val;
    size_t loc;
    o2scl::vector_max<vector<int>,int>(intv_obj.size(),
				       intv_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="size_t[]") {
    
    size_t val;
    size_t loc;
    o2scl::vector_max<vector<size_t>,size_t>(size_tv_obj.size(),
					     size_tv_obj,loc,val);
    cout << "Maximum value is " << val << " at index "
	 << loc << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_autocorr(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to compute autocorrelations with." << endl;
      return exc_efailed;
    }

    vector<string> in, pr;
    pr.push_back("Enter data column name");
    pr.push_back("Enter output for autocorrelations");
    pr.push_back("Enter output for 5*tau/m");
    int ret=get_input(sv,pr,in,"autocorr",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }

    if (!table_obj.is_column(in[1])) {
      table_obj.new_column(in[1]);
    }
    if (!table_obj.is_column(in[2])) {
      table_obj.new_column(in[2]);
    }

    // Compute autocorrelation length and sample size
    vector<double> ac_vec, ftom;
    vector_autocorr_vector(table_obj[in[0]],ac_vec);
    size_t len=vector_autocorr_tau(table_obj[in[0]],ac_vec,ftom);
    cout << "Autocorrelation length: " << len << " sample size: "
	 << table_obj.get_nlines()/len << endl;

    // Add autocorrelation and ftom data to table
    for(size_t i=0;i<table_obj.get_nlines();i++) {
      if (i<ac_vec.size()) {
	table_obj.set(in[1],i,ac_vec[i]);
      } else {
	table_obj.set(in[1],i,0.0);
      }
      if (i<ftom.size()) {
	table_obj.set(in[2],i,ftom[i]);
      } else {
	table_obj.set(in[2],i,0.0);
      }
    }

  } else if (type=="double[]") {

    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(doublev_obj,ac_vec,ftom);
    cout << "Autocorrelation length: " << len << " sample size: "
	 << doublev_obj.size()/len << endl;

    doublev_obj=ac_vec;

  } else if (type=="int[]") {

    vector_copy(intv_obj,doublev_obj);
    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(doublev_obj,ac_vec,ftom);
    cout << "Autocorrelation length: " << len << " sample size: "
	 << doublev_obj.size()/len << endl;

    command_del();
    clear_obj();
    doublev_obj=ac_vec;
    command_add("double[]");
    type="double[]";
    
  } else if (type=="size_t[]") {
    
    vector_copy(size_tv_obj,doublev_obj);
    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(doublev_obj,ac_vec,ftom);
    cout << "Autocorrelation length: " << len << " sample size: "
	 << doublev_obj.size()/len << endl;

    command_del();
    clear_obj();
    doublev_obj=ac_vec;
    command_add("double[]");
    type="double[]";
    
  }
  
  return 0;
}

int acol_manager::comm_to_table(std::vector<std::string> &sv, bool itive_com) {

  if (type=="double[]") {

    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(doublev_obj.size());
    for(size_t i=0;i<doublev_obj.size();i++) {
      table_obj.set(i1,i,doublev_obj[i]);
    }

    command_del();
    clear_obj();
    command_add("table");
    type="table";
    
  } else if (type=="int[]") {

    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(intv_obj.size());
    for(size_t i=0;i<intv_obj.size();i++) {
      table_obj.set(i1,i,intv_obj[i]);
    }

    command_del();
    clear_obj();
    command_add("table");
    type="table";

  } else if (type=="size_t[]") {
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column name",i1,"to-table",itive_com);
    if (ret!=0) return ret;
    
    table_obj.clear();
    table_obj.new_column(i1);
    table_obj.set_nlines(size_tv_obj.size());
    for(size_t i=0;i<size_tv_obj.size();i++) {
      table_obj.set(i1,i,size_tv_obj[i]);
    }

    command_del();
    clear_obj();
    command_add("table");
    type="table";
    
  }
  
  return 0;
}

int acol_manager::comm_download(std::vector<std::string> &sv, bool itive_com) {

  vector<string> in, pr;
  pr.push_back("Destination filename");
  pr.push_back("Hash");
  pr.push_back("URL");
  int ret=get_input(sv,pr,in,"download",itive_com);
  if (ret!=0) return ret;

  std::string file=in[0];
  std::string hash=in[1];
  std::string url=in[2];
  std::string fname;
  
  cloud_file cf;
  cf.get_file_hash(file,hash,url,fname);
  
  return 0;
}

int acol_manager::comm_slice(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    cerr << "Slice does not work with table objects." << endl;
  }
  
  if (type!="table3d" || table3d_obj.get_nslices()==0) {
    cerr << "No table3d with slices to find the maximum value of." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Slice in 'x' or 'y' direction");
  pr.push_back("Value to interpolate in for slice");
  int ret=get_input(sv,pr,in,"slice",itive_com);
  if (ret!=0) return ret;
  
  if (sv[1]=="x") {
    table3d_obj.extract_x(std::stod(sv[2]),table_obj);
    table3d_obj.clear();
    command_add("table");
    type="table";
  } else if (sv[1]=="y") {
    table3d_obj.extract_y(std::stod(sv[2]),table_obj);
    table3d_obj.clear();
    command_add("table");
    type="table";
  } else {
    cerr << "Invalid first argument to 'slice'." << endl;
  }
  
  return 0;
}

int acol_manager::comm_min(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d" || table3d_obj.get_nslices()==0) {
      cerr << "No table3d with slices to find the minimum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find minimum of",
			  i1,"min",itive_com);
    if (ret!=0) return ret;

    size_t ix;
    if (!table3d_obj.is_slice(i1,ix)) {
      cerr << "No slice named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    const ubmatrix &mat=table3d_obj.get_slice(ix);
    size_t i, j;
    double min;
    matrix_min_index(mat,i,j,min);

    cout << "Minimum value of slice '" << i1 << "' is: " 
	 << min << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="hist_2d") {
    
    const ubmatrix &mat=hist_2d_obj.get_wgts();
    size_t i, j;
    double min;
    matrix_min_index(mat,i,j,min);
    
    cout << "Minimum weight is: "
	 << min << " at indices (" << i << "," << j << ")\n  and grid "
	 << "point (" << table3d_obj.get_grid_x(i) << ","
	 << table3d_obj.get_grid_y(j) << ")." << endl;

  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to find the minimum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter column to find minimum of",
			  i1,"min",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    double min;
    size_t ix;
    vector_min(table_obj.get_nlines(),(table_obj)[i1],ix,min);
    cout << "Minimum value of column '" << i1 << "' is: " 
	 << min << " at row with index " << ix << "." << endl;
    
  } else if (type=="double[]") {

    double val;
    size_t loc;
    o2scl::vector_min<vector<double>,double>(doublev_obj.size(),
					     doublev_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="int[]") {
    
    int val;
    size_t loc;
    o2scl::vector_min<vector<int>,int>(intv_obj.size(),
				       intv_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  } else if (type=="size_t[]") {
    
    size_t val;
    size_t loc;
    o2scl::vector_min<vector<size_t>,size_t>(size_tv_obj.size(),
					     size_tv_obj,loc,val);
    cout << "Minimum value is " << val << " at index "
	 << loc << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_set_data(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d with data to set." << endl;
      return exc_efailed;
    }
    
    vector<string> in, pr;
    pr.push_back(table3d_obj.get_x_name()+" value of point to set");
    pr.push_back(table3d_obj.get_y_name()+" value of point to set");
    pr.push_back("Slice name to set");
    pr.push_back("New value");
    int ret=get_input(sv,pr,in,"set-data",itive_com);
    if (ret!=0) return ret;
    
    table3d_obj.set_val(o2scl::stod(in[0]),o2scl::stod(in[1]),in[2],
			o2scl::stod(in[3]));
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Function describing rows to modify");
  pr.push_back("Name of column to modify");
  pr.push_back("Function describing new values");
  int ret=get_input(sv,pr,in,"set-data",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[1])==false) {
    cerr << "Could not find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  for(size_t i=0;i<table_obj.get_nlines();i++) {
    if (table_obj.row_function(in[0],i)>0.5) {
      table_obj.set(in[1],i,table_obj.row_function(in[2],i));
    }
  }

  return 0;
}

int acol_manager::comm_set_unit(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to set units of." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to set units of");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"set-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  table_obj.set_unit(in[0],in[1]);

  return 0;
}

int acol_manager::comm_contours(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table3d" && type!="hist_2d") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  bool frac_mode=false;
  
  if (sv.size()>=2 && sv[1]=="frac") {
    cout << "Fraction mode is true." << endl;
    frac_mode=true;
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
  } else {
    cout << "Fraction mode is false." << endl;
  }
  if (sv.size()<2 && itive_com) {
    string temp=((string)"Enter \"frac\" for fractions of total sum and ")
      +"\"abs\" for absolute scale", i1;
    int ret=get_input_one(sv,temp,i1,"contours",itive_com);
    if (ret!=0) return ret;
    if (i1=="frac") frac_mode=true;
  }
    
  std::string svalue, file, name="contours";

  if (type=="table3d") {

    std::string slice;

    if (sv.size()<3) {
      svalue=cl->cli_gets("Contour value (or blank to cancel): ");
      slice=cl->cli_gets("Slice (or blank to cancel): ");
      if (svalue.length()==0) return 1;
      file=cl->cli_gets("Filename (or blank to keep): ");
      if (file.length()>0) {
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==3) {
      svalue=sv[1];
      slice=sv[2];
    } else if (sv.size()==4) {
      svalue=sv[1];
      slice=sv[2];
      file=sv[3];
    } else {
      svalue=sv[1];
      slice=sv[2];
      file=sv[3];
      name=sv[4];
    }

    ubvector levs(1);
    levs[0]=o2scl::stod(svalue);
    size_t nlev=1;
    
    if (frac_mode) {
      cout << "Fraction mode not implemented with table3d objects." << endl;
    } else {
      if (file.length()>0) {
	std::vector<contour_line> clines;
	table3d_obj.slice_contours(slice,1,levs,clines);
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,clines,name);
	hf.close();
      } else {
	table3d_obj.slice_contours(slice,1,levs,cont_obj);
	command_del();
	clear_obj();
	command_add("vector<contour_line>");
	type="vector<contour_line>";
      }
    }
    
  } else if (type=="hist_2d") {

    if (sv.size()<2) {
      svalue=cl->cli_gets("Contour value (or blank to cancel): ");
      if (svalue.length()==0) return 1;
      file=cl->cli_gets("Filename (or blank to keep): ");
      if (file.length()>0) {
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==2) {
      svalue=sv[1];
    } else if (sv.size()==3) {
      svalue=sv[1];
      file=sv[2];
    } else {
      svalue=sv[1];
      file=sv[2];
      name=sv[3];
    }

    
    ubvector levs(1);
    levs[0]=o2scl::stod(svalue);
    size_t nlev=1;

    if (frac_mode) {

      // Get references to the histogram data
      size_t nx=hist_2d_obj.size_x();
      size_t ny=hist_2d_obj.size_y();
      const ubmatrix &m=hist_2d_obj.get_wgts();
      const ubvector &xbins=hist_2d_obj.get_x_bins();
      const ubvector &ybins=hist_2d_obj.get_y_bins();

      // Compute the total integral and the target fraction
      double min, max;
      o2scl::matrix_minmax(m,min,max);
      double sum=hist_2d_obj.integ_wgts();
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  sum-=min*(xbins[i+1]-xbins[i])*(ybins[j+1]-ybins[j]);
	}
      }
      double target=levs[0]*sum;
      if (verbose>1) {
	cout << "sum,target: " << sum << " " << target << endl;
      }

      // Setup the vectors to interpolate the target integral
      uniform_grid_end<double> ug(min,max,100);
      ubvector integx, integy;
      ug.vector(integx);
      size_t N=integx.size();
      if (verbose>1) {
	cout << "N integx[0] integx[1]: " << N << " "
	     << integx[0] << " " << integx[1] << endl;
      }
      integy.resize(N);

      // Fill the interpolation vectors
      for(size_t k=0;k<N;k++) {
	integy[k]=0.0;
	for(size_t i=0;i<nx;i++) {
	  for(size_t j=0;j<ny;j++) {
	    if (m(i,j)>integx[k]) {
	      integy[k]+=(m(i,j)-min)*(xbins[i+1]-xbins[i])*
		(ybins[j+1]-ybins[j]);
	    }
	  }
	}
	if (verbose>1) {
	  cout << k << " " << integx[k] << " " << integy[k] << endl;
	}
      }

      // Perform the interpolation
      bool found=false;
      double level;
      for(size_t k=0;k<N-1;k++) {
	if (integy[k]>target && integy[k+1]<target) {
	  found=true;
	  level=integx[k]+(integx[k+1]-integx[k])*(target-integy[k])/
	    (integy[k+1]-integy[k]);
	}
      }
      if (verbose>1) {
	cout << "found: " << level << endl;
      }

      // Return if the interpolation failed
      if (found==false) {
	cerr << "Failed to find a level matching requested fraction."
	     << endl;
	return 2;
      }

      // Set level from interpolated value
      levs[0]=level;
      
    }

    contour co;
    co.set_levels(nlev,levs);
    
    ubvector xreps(hist_2d_obj.size_x());
    for (size_t i=0;i<hist_2d_obj.size_x();i++) {
      xreps[i]=hist_2d_obj.get_x_rep_i(i);
    }
    ubvector yreps(hist_2d_obj.size_y());
    for (size_t i=0;i<hist_2d_obj.size_y();i++) {
      yreps[i]=hist_2d_obj.get_y_rep_i(i);
    }
    co.set_data(hist_2d_obj.size_x(),hist_2d_obj.size_y(),xreps,yreps,
		hist_2d_obj.get_wgts());

    if (file.length()>0) {
      std::vector<contour_line> clines;
      co.calc_contours(clines);
      
      hdf_file hf;
      hf.open_or_create(file);
      hdf_output(hf,clines,name);
      hf.close();
    } else {
      command_del();
      clear_obj();
      co.calc_contours(cont_obj);
      command_add("vector<contour_line>");
      type="vector<contour_line>";
    }
    
  }
  
  return 0;
}

int acol_manager::comm_show_units(std::vector<std::string> &sv, 
				  bool itive_com) {
  if (cng.use_gnu_units) {
    cout << "Using GNU units? Yes." << endl;
    cout << "Variable 'unit_fname': " << unit_fname << endl;
  } else {
    cout << "Using GNU units? No." << endl;
  }
  cng.print_cache();
  return 0;
}

int acol_manager::comm_convert_unit
(std::vector<std::string> &sv, bool itive_com) {
  
  if (type!="table") {
    cerr << "Not implemented for " << type << " objects." << endl;
    return exc_efailed;
  }
  
  if (table_obj.get_nlines()==0) {
    cerr << "No table to convert units in." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column in which to convert units");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"convert-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  if (unit_fname.length()>0) {
    cng.units_cmd_string=((string)"units -f ")+unit_fname;
  }
  ret=table_obj.convert_to_unit(in[0],in[1],false);
  if (ret!=0) {
    cerr << "Could not find column or column does not have unit." << endl;
  }

  return 0;
}

int acol_manager::comm_get_conv
(std::vector<std::string> &sv, bool itive_com) {
  
  vector<string> in, pr;
  pr.push_back("Old unit");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"convert-unit",itive_com);
  if (ret!=0) return ret;
  
  if (unit_fname.length()>0) {
    cng.units_cmd_string=((string)"units -f ")+unit_fname;
    if (verbose>=2) {
      cout << "Units command string: " << cng.units_cmd_string
	   << endl;
    }
  }
  
  // Set the proper output precision and mode
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);

  if (verbose>=2) {
    cng.verbose=1;
  } else {
    cng.verbose=0;
  }

  // If cng.verbose is non-zero, then cng.convert may output
  // verbose information to cout
  double val;
  int cret=cng.convert_ret(in[0],in[1],1.0,val);
  if (cret!=0) {
    cerr << "Conversion failed." << endl;
    return 1;
  }
  
  cout << "Conversion factor is: " << val << endl;
  
  return 0;
}

int acol_manager::comm_get_unit(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to get units for." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to get units of");
  int ret=get_input(sv,pr,in,"get-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  cout << "Units of column " << in[0] << " are: " 
       << table_obj.get_unit(in[0]) << endl;

  return 0;
}

int acol_manager::comm_find_row(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0 || table_obj.get_nlines()==0) {
    cerr << "No table or empty table in find-row." << endl;
    return exc_efailed;
  }

  // If they didn't specify any parameters
  if (sv.size()==1) {
    
    string stmp=cl->cli_gets
      ("Enter function to be maximized or column name and value: ");
    vector<string> in;
    split_string(stmp,in);
    if (in.size()==2) {
      sv.push_back(in[0]);
      sv.push_back(in[1]);
    } else if (in.size()==1) {
      sv.push_back(in[0]);
    } else {
      cerr << "Gave three parameters but expected either one or two "
	   << "in find-row. Aborting." << endl;
      return exc_efailed;
    }
  
  } 
  
  if (sv.size()>=3) {
    
    // If they specified two parameters, then it is presumed that they
    // specified a column and a value
    
    if (table_obj.is_column(sv[1])) {
      size_t row=table_obj.lookup(sv[1],o2scl::function_to_double(sv[2]));
      
      // Call get_row() for the row that was found
      std::vector<std::string> sc;
      sc.push_back("-get-row");
      sc.push_back(itos(row));
      comm_get_row(sc,itive_com);
    }

    return 0;
  }

  // Otherwise, they specified a function to be maximized.

  // Find the appropriate row
  size_t row=table_obj.function_find_row(sv[1]);
  
  // Call get_row() for the row that was found
  std::vector<std::string> sc;
  sc.push_back("get-row");
  sc.push_back(itos(row));
  comm_get_row(sc,itive_com);
  
  return 0;
}

int acol_manager::comm_get_row(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  if (table_obj.get_nlines()==0 || table_obj.get_nlines()==0) {
    cerr << "No table or empty table in get-row." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter row number to get",
			i1,"get-row",itive_com);
  if (ret!=0) return ret;

  int ix=((int)(o2scl::function_to_double(i1)));
  
  // If negative, view as distance from end
  if (ix<0) ix+=table_obj.get_nlines();

  if (ix<0) {
    cerr << "Requested negative row in 'get-row'." << endl;
    return exc_efailed;
  }
  if (ix>((int)table_obj.get_nlines())-1) {
    cerr << "Requested row beyond end of table in get-row." << endl;
    return exc_efailed;
  }
  
  //--------------------------------------------------------------------
  // Compute number of screen columns

  if (user_ncols<=0) {
    char *ncstring=getenv("COLUMNS");
    if (ncstring) {
      int nc2;
      int sret=o2scl::stoi_nothrow(ncstring,nc2);
      if (sret==0 && nc2>0) {
	ncols=nc2;
      } else {
	cerr << "Failed to interpret COLUMNS value " << ncstring
	     << " as a positive number of columns." << endl;
      }
    }
  } else {
    ncols=user_ncols;
  }

  //--------------------------------------------------------------------
  // Temporary storage strings for names and data

  vector<string> row_names, row_data;

  //--------------------------------------------------------------------
  // Process and/or output names

  if (names_out==true) {

    if (pretty) {

      size_t running_width=0;
      ostringstream str;
      // Clear ostringstream with str.str(""); and str.clear();
      str.setf(ios::scientific);
      str.precision(prec);
      
      for(size_t i=0;i<table_obj.get_ncolumns();i++) {

	// Count for space between columns and sign
	size_t this_col=2;
	// Count column name
	this_col+=table_obj.get_column_name(i).size();
	// Count extra spaces to format number
	int num_spaces=prec+6-((int)(table_obj.get_column_name(i).size()));
	if (num_spaces>0) this_col+=num_spaces;
	// See if there will be space
	if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	  row_names.push_back(str.str());
	  str.str("");
	  str.clear();
	  str.setf(ios::scientific);
	  str.precision(prec);
	  running_width=0;
	}
	// Output this column name
	str << ' ' << table_obj.get_column_name(i) << ' ';
	for(int j=0;j<num_spaces;j++) {
	  str << ' ';
	}
	running_width+=this_col;
      }
      row_names.push_back(str.str());
      str.str("");
      str.clear();
      
    } else {
      
      cout.precision(prec);
  
      for(size_t i=0;i<table_obj.get_ncolumns();i++) {
	cout << table_obj.get_column_name(i) << ' ';
      }
      cout << endl;

    }
  }
  
  //--------------------------------------------------------------------
  // Process and/or output data
  
  if (pretty) {
    
    size_t running_width=0;
    ostringstream str;
    str.setf(ios::scientific);
    str.precision(prec);
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      
      // Count space for number
      size_t this_col=prec+8;
      // Count extra spaces if necessary
      int num_spaces=((int)(table_obj.get_column_name(i).size())-prec-6);
      if (num_spaces>0) this_col+=num_spaces;
      // See if there will be space
      if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	row_data.push_back(str.str());
	str.str("");
	str.clear();
	str.setf(ios::scientific);
	str.precision(prec);
	running_width=0;
      }
      // Output the data
      if (table_obj.get(i,ix)>=0.0) {
	str << ' ' << table_obj.get(i,ix) << ' ';
      } else {
	str << table_obj.get(i,ix) << ' ';
      }
      for(int j=0;j<num_spaces;j++) {
	str << ' ';
      }
      running_width+=this_col;
    }
    row_data.push_back(str.str());
    str.str("");
    str.clear();
    
    //--------------------------------------------------------------------
    // Now output both names and data to cout

    if (row_names.size()!=row_data.size()) {
      O2SCL_ERR("Names and data size don't match in get-row.",
		exc_esanity);
    }
    for(size_t k=0;k<row_names.size();k++) {
      cout << row_names[k] << endl;
      cout << row_data[k] << endl;
    }
    
  } else {
    
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      cout << table_obj.get(i,ix) << ' ';
    }
    cout << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_rename(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {

    vector<string> pr, in;
    pr.push_back("Enter slice to be renamed");
    pr.push_back("Enter new name");
    
    int ret=get_input(sv,pr,in,"rename",itive_com);
    if (ret!=0) return ret;
    
    table3d_obj.rename_slice(in[0],in[1]);
    
    return 0;

  } else if (type=="table") {

    vector<string> pr, in;
    pr.push_back("Enter column to be renamed");
    pr.push_back("Enter new name");
    
    int ret=get_input(sv,pr,in,"rename",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.new_column(in[1]);
    
    table_obj.copy_column(in[0],in[1]);
    table_obj.delete_column(in[0]);

  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }    

  return 0;
}

int acol_manager::comm_deriv(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to take derivatives of." << endl;
      return exc_efailed;
    }
    
    vector<string> pr, in;
    pr.push_back("Enter 'x' column");
    pr.push_back("Enter 'y' column");
    pr.push_back("Enter name of new column");
    int ret=get_input(sv,pr,in,"deriv",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    if (table_obj.is_column(in[1])==false) {
      cerr << "Could not find column named '" << in[1] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.deriv(in[0],in[1],in[2]);
    
  } else if (type=="double[]") {
    
    std::vector<double> vderiv(doublev_obj.size());
    o2scl::vector_deriv_interp(doublev_obj.size(),doublev_obj,vderiv,
			       interp_type);
    doublev_obj=vderiv;
    
  } else if (type=="int[]") {

    doublev_obj.resize(intv_obj.size());
    o2scl::vector_copy(intv_obj,doublev_obj);
    std::vector<double> vderiv(doublev_obj.size());
    o2scl::vector_deriv_interp(doublev_obj.size(),doublev_obj,vderiv,
			       interp_type);
    o2scl::vector_copy(vderiv,intv_obj);
    doublev_obj.resize(0);
    
  } else if (type=="size_t[]") {
    
    doublev_obj.resize(size_tv_obj.size());
    o2scl::vector_copy(size_tv_obj,doublev_obj);
    std::vector<double> vderiv(doublev_obj.size());
    o2scl::vector_deriv_interp(doublev_obj.size(),doublev_obj,vderiv,
			       interp_type);
    o2scl::vector_copy(vderiv,size_tv_obj);
    doublev_obj.resize(0);
    
  } else {
    
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
    
  }

  return 0;
}

int acol_manager::comm_deriv_x(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table3d") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  if (table3d_obj.get_nslices()==0) {
    cerr << "No table3d with slices to take derivatives of." << endl;
    return exc_efailed;
  }
  
  vector<string> pr, in;
  pr.push_back("Enter slice containing function");
  pr.push_back("Enter name of new slice");
  int ret=get_input(sv,pr,in,"deriv",itive_com);
  if (ret!=0) return ret;

  size_t iz;
  if (table3d_obj.is_slice(in[0],iz)==false) {
    cerr << "Could not find slice named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  
  table3d_obj.deriv_x(in[0],in[1]);

  return 0;
}

int acol_manager::comm_deriv_y(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table3d") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  if (table3d_obj.get_nslices()==0) {
    cerr << "No table3d with slices to take derivatives of." << endl;
    return exc_efailed;
  }
  
  vector<string> pr, in;
  pr.push_back("Enter slice containing function");
  pr.push_back("Enter name of new slice");
  int ret=get_input(sv,pr,in,"deriv",itive_com);
  if (ret!=0) return ret;
  
  size_t iz;
  if (table3d_obj.is_slice(in[0],iz)==false) {
    cerr << "Could not find slice named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  
  table3d_obj.deriv_y(in[0],in[1]);

  return 0;
}

int acol_manager::comm_deriv2(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cout << "No table with columns to take derivatives of." << endl;
      return exc_efailed;
    }
    vector<string> pr, in;
    pr.push_back("Enter 'x' column");
    pr.push_back("Enter 'y' column");
    pr.push_back("Enter name of new column");
    int ret=get_input(sv,pr,in,"deriv2",itive_com);
    if (ret!=0) return ret;
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    if (table_obj.is_column(in[1])==false) {
      cerr << "Could not find column named '" << in[1] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.deriv2(in[0],in[1],in[2]);
    
  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }    

  return 0;
}

int acol_manager::comm_integ(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table with columns to integrate." << endl;
    return exc_efailed;
  }
  vector<string> pr, in;
  pr.push_back("Enter 'x' column");
  pr.push_back("Enter 'y' column");
  pr.push_back("Enter name of new column");
  int ret=get_input(sv,pr,in,"integ",itive_com);
  if (ret!=0) return ret;

  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  if (table_obj.is_column(in[1])==false) {
    cerr << "Could not find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  table_obj.integ(in[0],in[1],in[2]);

  return 0;
}

int acol_manager::comm_internal(std::vector<std::string> &sv, bool itive_com) {
  
  std::string i1;
  int ret=get_input_one(sv,"Enter filename",i1,"internal",itive_com);
  if (ret!=0) return ret;
  
  if (verbose>2) {
    cout << "Opening O2scl file: " << i1 << endl;
  }
  
  hdf_file hf;
  hf.compr_type=compress;
  hf.open_or_create(i1);
  
  if (type=="int") {

    hf.seti(obj_name,int_obj);
    
  } else if (type=="double") {

    hf.setd(obj_name,double_obj);
    
  } else if (type=="char") {

    hf.setc(obj_name,char_obj);
    
  } else if (type=="string") {

    hf.sets(obj_name,string_obj);
    
  } else if (type=="size_t") {

    hf.set_szt(obj_name,size_t_obj);
    
  } else if (type=="double[]") {

    hf.setd_vec(obj_name,doublev_obj);
    
  } else if (type=="int[]") {

    hf.seti_vec(obj_name,intv_obj);
    
  } else if (type=="size_t[]") {

    hf.set_szt_vec(obj_name,size_tv_obj);
    
  } else if (type=="string[]") {

    hf.sets_vec(obj_name,stringv_obj);
    
  } else if (type=="table3d") {
    
    hdf_output(hf,table3d_obj,obj_name);
    
  } else if (type=="table") {
    
    hdf_output(hf,table_obj,obj_name);

  } else if (type=="hist") {

    hdf_output(hf,hist_obj,obj_name);
    
  } else if (type=="hist_2d") {

    hdf_output(hf,hist_2d_obj,obj_name);
    
  } else if (type=="vector<contour_line>") {

    hdf_output(hf,cont_obj,obj_name);
    
  } else if (type=="uniform_grid<double>") {

    hdf_output(hf,ug_obj,obj_name);
    
  }

  hf.close();
    
  return 0;
}

int acol_manager::comm_entry(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {

    vector<string> pr, in;
    pr.push_back("Enter column name");
    pr.push_back("Enter row index");
    int ret=get_input(sv,pr,in,"function",itive_com);
    if (ret!=0) return ret;

    int row;
    int ret2=o2scl::stoi_nothrow(in[1],row);
    if (ret2!=0) {
      std::cerr << "Failed to convert " << in[1]
		<< " to a number." << endl;
      return exc_efailed;
    }
    if (row<0) {
      std::cerr << "Conversion of " << in[1]
		<< " resulted in, " << row << ", a negative number." << endl;
      return exc_efailed;
    }
    
    cout << "Entry for column " << in[0] << " at row " << in[1] << " is "
	 << table_obj.get(in[0],row) << endl;
    
  } else {
    cerr << "Command 'entry' not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_function(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    vector<string> pr, in;
    pr.push_back("Enter function for new slice");
    pr.push_back("Enter name for new slice");
    int ret=get_input(sv,pr,in,"function",itive_com);
    if (ret!=0) return ret;
    
    // Remove single or double quotes just in case
    if (in[0].size()>=3 && ((in[0][0]=='\'' && in[0][in[0].size()-1]=='\'') ||
			    (in[0][0]=='\"' && in[0][in[0].size()-1]=='\"'))) {
      in[0]=in[0].substr(1,in[0].size()-2);
    }
    
    table3d_obj.function_slice(in[0],in[1]);
    if (ret!=0) {
      cerr << "Function make_slice() failed." << endl;
      return exc_efailed;
    }
    
    return 0;

  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to add a column to." << endl;
      return exc_efailed;
    }
    
    vector<string> pr, in;
    pr.push_back("Enter function for new column");
    pr.push_back("Enter name for new column");
    int ret=get_input(sv,pr,in,"function",itive_com);
    if (ret!=0) return ret;
    
    // Remove single or double quotes just in case
    if (in[0].size()>=3 && ((in[0][0]=='\'' && in[0][in[0].size()-1]=='\'') ||
			    (in[0][0]=='\"' && in[0][in[0].size()-1]=='\"'))) {
      in[0]=in[0].substr(1,in[0].size()-2);
    }
    
    if (table_obj.is_column(in[1])==true) {
      cerr << "Already a column named '" << in[1] << "'." << endl;
      return exc_efailed;
    }
    
    table_obj.function_column(in[0],in[1]);
    
    if (ret!=0) {
      cerr << "Function make_col() failed." << endl;
      return exc_efailed;
    }

  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_calc(std::vector<std::string> &sv, bool itive_com) {

  std::string i1;
  if (sv.size()>1) {
    i1=sv[1];
  } else if (itive_com) {
    i1=cl->cli_gets("Enter expression to compute (or blank to stop): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'calc' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No expression to compute in 'calc'." << endl;
    return exc_efailed;
  }
  double d=o2scl::function_to_double(i1);
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  if (verbose>0) cout << "Result: ";
  cout << d << endl;
  return 0;
}

int acol_manager::comm_index(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to add line numbers to." << endl;
    return exc_efailed;
  }

  std::string i1="N";
  if (sv.size()>1) i1=sv[1];
  table_obj.new_column(i1);
  for(size_t i=0;i<table_obj.get_nlines();i++) table_obj.set(i1,i,((double)i));

  return 0;
}

#ifdef O2SCL_NEVER_DEFINED
int acol_manager::comm_html(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to output." << endl;
    return exc_efailed;
  }

  ostream *fout;
  ofstream ffout;

  std::string i1;
  if (sv.size()==1) {
    if (itive_com) {
      i1=cl->cli_gets("Enter filename (or blank to quit): ");
      if (i1.length()==0) {
	if (verbose>0) cout << "Command 'html' cancelled." << endl;
	return 0;
      }
    } else {
      cout << "Filename not given for 'html'." << endl;
      return exc_efailed;
    }
  } else {
    i1=sv[1];
  }
   
  ffout.open(i1.c_str());
  fout=&ffout;
  
  if (scientific) fout->setf(ios::scientific);
  else fout->unsetf(ios::scientific);

  string nlast;
  
  (*fout) << "<html><head></head><body>" << endl;
  
  //--------------------------------------------------------------------
  // Output constants
  
  for(size_t i=0;i<table_obj.get_nconsts();i++) {
    if (i==0) (*fout) << "<b>Constants:</b><br />" << endl;
    string tnam;
    double tval;
    table_obj.get_constant(i,tnam,tval);
    (*fout) << tnam << " = " << tval << "<br />" << endl;
  }
  
  if (table_obj.get_ncolumns()>0) {
    
    (*fout) << "<table border=\"0\">" << endl;
    
    //--------------------------------------------------------------------
    // Output column names
    
    (*fout) << "<tr bgcolor=\"#dddddd\">";
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      (*fout) << "<td><b>" << table_obj.get_column_name(i) << "</b></td>";
    }
    (*fout) << "</tr>" << endl;
      
    //--------------------------------------------------------------------
    // Output data
      
    for(int i=0;i<((int)table_obj.get_nlines());i++) {
      if (i%5==4) (*fout) << "<tr bgcolor=\"#dddddd\">";
      else (*fout) << "<tr>";
      for(size_t j=0;j<table_obj.get_ncolumns();j++) {
	(*fout) << "<td>" << ffl.convert(table_obj.get(j,i)) << "</td>";
      }
      (*fout) << "</tr>" << endl;
    }

    (*fout) << "</table>" << endl;

  } else {

    (*fout) << "(no data)" << endl;

  }
  
  (*fout) << "</body></html>" << endl;

  ffout.close();
  
  return 0;
}
#endif

int acol_manager::comm_value(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()>1) {
    if (type=="int") {
      int_obj=o2scl::stoi(sv[1]);
    } else if (type=="double") {
      double_obj=o2scl::function_to_double(sv[1]);
    } else if (type=="char") {
      char_obj=sv[1][0];
    } else if (type=="size_t") {
      size_t_obj=o2scl::stoszt(sv[1]);
    } else if (type=="string") {
      string_obj=sv[1];
    }
  }
  
  if (type=="int") {
    cout << "Value of " << obj_name << " is " << int_obj << endl;
  } else if (type=="double") {
    cout << "Value of " << obj_name << " is " << double_obj << endl;
  } else if (type=="char") {
    cout << "Value of " << obj_name << " is " << char_obj << endl;
  } else if (type=="size_t") {
    cout << "Value of " << obj_name << " is " << size_t_obj << endl;
  } else if (type=="string") {
    cout << "Value of " << obj_name << " is " << string_obj << endl;
  }
  
  return 0;
}
  
int acol_manager::comm_preview(std::vector<std::string> &sv, bool itive_com) {

  if (type.length()==0) {
    cerr << "No object to preview." << endl;
    return 3;
  }
  
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  
  cout.precision(prec);
  
  if (user_ncols<=0) {
    char *ncstring=getenv("COLUMNS");
    if (ncstring) {
      int nc2;
      int sret=o2scl::stoi_nothrow(ncstring,nc2);
      if (sret==0 && nc2>0) {
	ncols=nc2;
      } else {
	cerr << "Failed to interpret COLUMNS value " << ncstring
	     << " as a positive number of columns." << endl;
      }
    }
  } else {
    ncols=user_ncols;
  }

  cout << "Type: " << type << " Name: " << obj_name << endl;
  
  if (type=="table3d") {
    
    int lmar=table3d_obj.get_y_name().length()+1;
    
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    if (nx==0 || ny==0) {
      cout << "No size set. Blank table3d." << endl;
    } else {
      int nrows, ncls;
      if (sv.size()>=2) {
	int sret=o2scl::stoi_nothrow(sv[1],nrows);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[1]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	nrows=10;
      }
      if (((size_t)nrows)>nx) nrows=nx;
      if (sv.size()>=3) {
	int sret=o2scl::stoi_nothrow(sv[2],ncls);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[2]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	// 8+prec for the grid point, 4 for extra spacing,
	// and lmar for the left margin which has the x label
	if (ncols<=prec+lmar+12) ncls=1;
	else ncls=(ncols-prec-12-lmar)/(prec+8);
	std::cout << "Screen width: " << ncols << " prec: " << prec
		  << " lmar: " << lmar << " flag: " << (ncols<=prec+lmar+12)
		  << " ncols: " << ncls << endl;
      }
      if (((size_t)ncls)>ny) ncls=ny;
      int dx=nx/ncls;
      int dy=ny/nrows;
      if (dx==0) dx=1;
      if (dy==0) dy=1;

      cout << "x: " << table3d_obj.get_x_name() << " [";
      if (nx<4) {
	for(size_t i=0;i<nx-1;i++) {
	  cout << table3d_obj.get_grid_x(i) << " ";
	}
	cout << table3d_obj.get_grid_x(nx-1) << "] ";
      } else {
	cout << table3d_obj.get_grid_x(0) << " ";
	cout << table3d_obj.get_grid_x(1) << " ... ";
	cout << table3d_obj.get_grid_x(nx-2) << " ";
	cout << table3d_obj.get_grid_x(nx-1) << "] ";
      }
      cout << endl;
      cout << "y: " << table3d_obj.get_y_name() << " [";
      if (ny<4) {
	for(size_t i=0;i<ny-1;i++) {
	  cout << table3d_obj.get_grid_y(i) << " ";
	}
	cout << table3d_obj.get_grid_y(ny-1) << "] ";
      } else {
	cout << table3d_obj.get_grid_y(0) << " ";
	cout << table3d_obj.get_grid_y(1) << " ... ";
	cout << table3d_obj.get_grid_y(ny-2) << " ";
	cout << table3d_obj.get_grid_y(ny-1) << "] ";
      }
      cout << endl;
      
      size_t nt=table3d_obj.get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  // Slice name
	  cout << "Slice " << k << ": "
	       << table3d_obj.get_slice_name(k) << endl;

	  // Label row
	  for(int i=0;i<lmar+14;i++) cout << " ";
	  cout << "   " << table3d_obj.get_x_name() << endl;

	  // Set showpos
	  cout.setf(ios::showpos);

	  // Grid row
	  for(size_t i=0;i<((size_t)prec)+8+lmar;i++) cout << " ";
	  cout << "| ";
	  
	  for(size_t i=0;i<((size_t)ncls);i++) {
	    cout << table3d_obj.get_grid_x(i*dx) << " ";
	  }
	  cout << endl;

	  // Divider row
	  for(int i=0;i<lmar;i++) cout << " ";
	  for(size_t i=0;i<((size_t)prec)+8;i++) cout << "-";
	  cout << "|";
	  for(size_t i=0;i<((size_t)ncls)*(prec+8);i++) {
	    cout << "-";
	  }
	  cout << endl;

	  // Data output
	  for(int j=((int)nrows)-1;j>=0;j--) {
	    if (j==((int)nrows)-1) {
	      cout << table3d_obj.get_y_name() << " ";
	    } else {
	      for(int i=0;i<lmar;i++) cout << " ";
	    }
	    cout << table3d_obj.get_grid_y(j*dy) << " | ";
	    for(size_t i=0;i<((size_t)ncls);i++) {
	      cout << table3d_obj.get(i*dx,j*dy,k) << " ";
	    }
	    cout << endl;
	  }

	  // Unset showpos
	  cout.unsetf(ios::showpos);

	  // Newline between slices
	  cout << endl;
	}
      }
    }
    return 0;

  } else if (type=="hist") {

    int nrows=10;
    if (sv.size()>=2) {
      int sret=o2scl::stoi_nothrow(sv[1],nrows);
      if (sret!=0 || nrows<=0) {
	cerr << "Failed to interpret " << sv[1]
	     << "as a positive number of rows." << endl;
	return exc_einval;
      }
    }
    int inr=(hist_obj.size()+(nrows-1))/(nrows);
    if (inr<1) inr=1;

    cout.precision(prec);

    cout.setf(ios::left);
    cout.width(prec+8);
    cout << " low";
    cout.width(prec+8);
    cout << " high";
    cout.width(prec+8);
    cout << " weight" << endl;
    cout.unsetf(ios::left);

    for(int i=0;i<((int)hist_obj.size());i+=inr) {
      double val=hist_obj.get_bin_low_i(i);
      if (val<0.0) {
	cout << val << " ";
      } else {
	cout << " " << val << " ";
      }
      val=hist_obj.get_bin_high_i(i);
      if (val<0.0) {
	cout << val << " ";
      } else {
	cout << " " << val << " ";
      }
      val=hist_obj.get_wgt_i(i);
      if (val<0.0) {
	cout << val << endl;
      } else {
	cout << " " << val << endl;
      }
    }
    
    return 0;

  } else if (type=="hist_2d") {

    size_t nx=hist_2d_obj.size_x(), ny=hist_2d_obj.size_y();
    if (nx==0 || ny==0) {
      cout << "No size set. Blank hist_2d." << endl;
    } else {

      int nrows, ncls;
      if (sv.size()>=2) {
	int sret=o2scl::stoi_nothrow(sv[1],nrows);
	if (sret!=0 || nrows<=0) {
	  cerr << "Failed to interpret " << sv[1]
	       << "as a positive number of rows." << endl;
	  return exc_einval;
	}
      } else {
	nrows=10;
      }
      if (((size_t)nrows)>nx) nrows=nx;
      if (sv.size()>=3) {
	int sret=o2scl::stoi_nothrow(sv[2],ncls);
	if (sret!=0 || ncls<=0) {
	  cerr << "Failed to interpret " << sv[2]
	       << "as a positive number of columns.." << endl;
	  return exc_einval;
	}
      } else {
	// 8+prec for the grid point, 3 for extra spacing,
	if (ncols<=prec+11) ncls=1;
	else ncls=(ncols-prec-11)/(prec+8);
	std::cout << "Screen width: " << ncols << " prec: " << prec
		  << " flag: " << (ncols<=prec+11)
		  << " ncols: " << ncls << endl;
      }
      if (((size_t)ncls)>ny) ncls=ny;
      size_t dx=nx/ncls;
      size_t dy=ny/nrows;
      if (dx==0) dx=1;
      if (dy==0) dy=1;

      cout << "nx: " << nx << " ny: " << ny << endl;
      if (nx<=4) {
	cout << "x: [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_low_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_low_i(nx-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_rep_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_rep_i(nx-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<nx-1;i++) {
	  cout << hist_2d_obj.get_x_high_i(i) << " ";
	}
	cout << hist_2d_obj.get_x_high_i(nx-1) << "]" << endl;
      } else {
	cout << "x: [";
	cout << hist_2d_obj.get_x_low_i(0) << " ";
	cout << hist_2d_obj.get_x_low_i(1) << " ... ";
	cout << hist_2d_obj.get_x_low_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_low_i(nx-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_x_rep_i(0) << " ";
	cout << hist_2d_obj.get_x_rep_i(1) << " ... ";
	cout << hist_2d_obj.get_x_rep_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_rep_i(nx-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_x_high_i(0) << " ";
	cout << hist_2d_obj.get_x_high_i(1) << " ... ";
	cout << hist_2d_obj.get_x_high_i(nx-2) << " ";
	cout << hist_2d_obj.get_x_high_i(nx-1) << "]" << endl;
      }
      cout << endl;
      if (ny<=4) {
	cout << "y: [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_low_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_low_i(ny-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_rep_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_rep_i(ny-1) << "]" << endl;
	cout << "   [";
	for(size_t i=0;i<ny-1;i++) {
	  cout << hist_2d_obj.get_y_high_i(i) << " ";
	}
	cout << hist_2d_obj.get_y_high_i(ny-1) << "]" << endl;
      } else {
	cout << "y: [";
	cout << hist_2d_obj.get_y_low_i(0) << " ";
	cout << hist_2d_obj.get_y_low_i(1) << " ... ";
	cout << hist_2d_obj.get_y_low_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_low_i(ny-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_y_rep_i(0) << " ";
	cout << hist_2d_obj.get_y_rep_i(1) << " ... ";
	cout << hist_2d_obj.get_y_rep_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_rep_i(ny-1) << "]" << endl;
	cout << "   [";
	cout << hist_2d_obj.get_y_high_i(0) << " ";
	cout << hist_2d_obj.get_y_high_i(1) << " ... ";
	cout << hist_2d_obj.get_y_high_i(ny-2) << " ";
	cout << hist_2d_obj.get_y_high_i(ny-1) << "]" << endl;
      }
      cout << endl;
      
      // Set showpos
      cout.setf(ios::showpos);
      
      // Grid row
      for(size_t i=0;i<((size_t)prec)+8;i++) cout << " ";
      cout << "| ";
      
      for(size_t i=0;i<((size_t)ncls);i++) {
	cout << hist_2d_obj.get_x_rep_i(i*dx) << " ";
      }
      cout << endl;
      
      // Divider row
      for(size_t i=0;i<((size_t)prec)+8;i++) cout << "-";
      cout << "|";
      for(size_t i=0;i<((size_t)ncls)*(prec+8);i++) {
	cout << "-";
      }
      cout << endl;
      
      // Data output
      for(int j=((int)nrows)-1;j>=0;j--) {
	cout << hist_2d_obj.get_y_rep_i(j*dy) << " | ";
	for(size_t i=0;i<((size_t)ncls);i++) {
	  cout << hist_2d_obj.get_wgt_i(i*dx,j*dy) << " ";
	}
	cout << endl;
      }
      
      // Unset showpos
      cout.unsetf(ios::showpos);
      
    }
    
    return 0;
    
  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to preview." << endl;
      return exc_efailed;
    }
    
    if (scientific) cout.setf(ios::scientific);
    else cout.unsetf(ios::scientific);
    
    if (user_ncols<=0) {
      char *ncstring=getenv("COLUMNS");
      if (ncstring) {
	int nc2;
	int sret=o2scl::stoi_nothrow(ncstring,nc2);
	if (sret==0 && nc2>0) {
	  ncols=nc2;
	} else {
	  cerr << "Failed to interpret COLUMNS value " << ncstring
	       << " as a positive number of columns." << endl;
	}
      }
    } else {
      ncols=user_ncols;
    }
    
    cout.precision(prec);
    
    if (table_obj.get_ncolumns()>0) {
      
      string nlast;
      int inr;
      int inc;
      
      //----------------------------------------------------------------------
      // Compute number of columns which will fit
      
      size_t max_cols=(ncols)/(8+prec);
      if (max_cols>table_obj.get_ncolumns()) max_cols=table_obj.get_ncolumns();
      
      //--------------------------------------------------------------------
      // Compute column and row increment
      
      if (sv.size()==2) {
	int nrows;
	int ret=o2scl::stoi_nothrow(sv[1],nrows);
	if (ret!=0 || nrows<=0) {
	  std::cerr << "Could not interpret " << sv[1]
		    << " as a positive and nonzero number of rows." << endl;
	  return exc_efailed;
	}
	inr=(table_obj.get_nlines()+(nrows-1))/(nrows);
	if (inr<1) inr=1;
      } else {
	inr=(table_obj.get_nlines()+9)/10;
	if (inr<1) inr=1;
      }
      inc=(table_obj.get_ncolumns()+(1))/max_cols;
      if (inc<1) inc=1;
      
      //--------------------------------------------------------------------
      // Get last row number if necessary
      
      if (pretty==true) {
	nlast=itos(table_obj.get_nlines()-1);
      }
      
      //--------------------------------------------------------------------
      // Output column names
      
      if (names_out==true) {
	
	for(size_t ki=0;ki<max_cols;ki++) {
	  
	  size_t i=ki*inc;
	  if (i>=table_obj.get_ncolumns()) i=table_obj.get_ncolumns()-1;
	  
	  // Preceeding space
	  if (pretty==true) {
	    cout << ' ';
	  }
	  
	  // Column name
	  cout << table_obj.get_column_name(i) << " ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=prec+6-((int)(table_obj.get_column_name(i).size()));
	    for(int j=0;j<nsp;j++) cout << ' ';
	  } else {
	    for(size_t kk=1;kk<nlast.length();kk++) cout << ' ';
	  }
	  
	}
	cout << endl;
      }
      
      //--------------------------------------------------------------------
      // Output units
      
      if (names_out==true && table_obj.get_nunits()>0) {
	
	for(size_t ki=0;ki<max_cols;ki++) {
	  
	  size_t i=ki*inc;
	  if (i>=table_obj.get_ncolumns()) i=table_obj.get_ncolumns()-1;
	  
	  // Preceeding space
	  if (pretty==true) {
	    cout << ' ';
	  }
	  
	  // Column name
	  string cunit=table_obj.get_unit(table_obj.get_column_name(i));
	  cout << '[' << cunit << "] ";
	  
	  // Trailing spaces
	  if (pretty==true) {
	    int nsp=prec+6-cunit.size()-2;
	    if (nsp<0) nsp=0;
	    for(int j=0;j<nsp;j++) cout << ' ';
	  } else {
	    for(size_t kk=1;kk<nlast.length();kk++) cout << ' ';
	  }
	  
	}
	cout << endl;
      }
      
      //--------------------------------------------------------------------
      // Output data
      
      for(size_t i=0;i<table_obj.get_nlines();i+=inr) {
	
	for(size_t kj=0;kj<max_cols;kj++) {
	  
	  size_t j=kj*inc;
	  if (j>=table_obj.get_ncolumns()) j=table_obj.get_ncolumns()-1;
	  
	  if (pretty==true) {
	    double d=table_obj.get(j,i);
	    if (!has_minus_sign(&d)) {
	      cout << ' ';
	    }
	  }
	  cout << table_obj.get(j,i) << ' ';
	  if (pretty==true) {
	    for(int kk=0;kk<((int)(table_obj.get_column_name(j).size()-
				   prec-6));kk++) {
	      cout << ' ';
	    }
	  }
	  
	}
	cout << endl;
	
	//--------------------------------------------------------------------
	// Continue to next row
      }
      
    }
    
    return 0;

  } else if (type=="vector<contour_line>") {

    for(size_t i=0;i<cont_obj.size();i++) {
      cout << "Value: " << cont_obj[i].level << endl;
      for(size_t j=0;j<cont_obj[i].x.size();j++) {
	cout << cont_obj[i].x[j] << " ";
	cout << cont_obj[i].y[j] << endl;
      }
      cout << endl;
    }
    
    return 0;
  } else if (type=="int") {
    cout << "Value is " << int_obj << endl;
    return 0;
  } else if (type=="double") {
    cout << "Value is " << double_obj << endl;
    return 0;
  } else if (type=="char") {
    cout << "Value is " << char_obj << endl;
    return 0;
  } else if (type=="string") {
    cout << "Value is " << string_obj << endl;
    return 0;
  } else if (type=="size_t") {
    cout << "Value is " << size_t_obj << endl;
    return 0;
  } else if (type=="int[]") {
    
    vector<string> inc, outc;
    for(size_t i=0;i<intv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". "+o2scl::itos(intv_obj[i]);
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="double[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<doublev_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". "+o2scl::dtos(doublev_obj[i]);
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="size_t[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<size_tv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". "+o2scl::szttos(size_tv_obj[i]);
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="string[]") {

    vector<string> inc, outc;
    for(size_t i=0;i<stringv_obj.size();i++) {
      string tmp=o2scl::szttos(i)+". \""+stringv_obj[i]+"\"";
      inc.push_back(tmp);
    }
    screenify(inc.size(),inc,outc);
    
    int inr;
    if (sv.size()==2) {
      int nrows;
      int ret=o2scl::stoi_nothrow(sv[1],nrows);
      if (ret!=0 || nrows<=0) {
	std::cerr << "Could not interpret " << sv[1]
		  << " as a positive and nonzero number of rows." << endl;
	return exc_efailed;
      }
      inr=(outc.size()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(outc.size()+9)/10;
      if (inr<1) inr=1;
    }

    for(size_t i=0;i<outc.size();i+=inr) {
      cout << outc[i] << endl;
    }
    return 0;

  } else if (type=="uniform_grid<double>") {
    cout << "Uniform grid " << obj_name << endl;
    cout << "Number of bins: " << ug_obj.get_nbins() << endl;
    cout << "Start: " << ug_obj.get_start() << endl;
    cout << "End: " << ug_obj.get_end() << endl;
    cout << "Width: " << ug_obj.get_width() << endl;
    return 0;
  }

  cerr << "Cannot preview type " << type << "." << endl;
  return 0;
}

int acol_manager::comm_interactive(std::vector<std::string> &sv, 
				   bool itive_com) {

  post_interactive=!post_interactive;
  if (verbose>0) {
    if (post_interactive) {
      cout << "Interactive mode will run after command-line is parsed." 
	   << endl;
    } else {
      cout << "Interactive mode will not run after command-line is parsed." 
	   << endl;
    }
  }
  return 0;
}

#ifdef O2SCL_NEVER_DEFINED
int acol_manager::make_unique_name(string &col, std::vector<string> &cnames) {
  bool done;
  do {
    done=true;
    for(size_t i=0;i<cnames.size();i++) {
      if (col==cnames[i]) {
        done=false;
        i=cnames.size();
      }
    }
    if (done==false) {
      col+='_';
    }
  } while (done==false);
  return 0;
}
#endif

int acol_manager::comm_generic(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype;

  // Delete previous object
  command_del();
  clear_obj();
  
  int ret=get_input_one(sv,"Enter type of object to create",ctype,"create",
			itive_com);
  if (ret!=0) return ret;

  vector<string> sv2=sv;
  vector<string>::iterator it=sv2.begin();
  sv2.erase(it+1);
  
  // Open the file 
  ifstream ifs;
  if (sv2[1]!=((std::string)"cin")) {
    ifs.open(sv2[1].c_str());
    if (!(ifs)) {
      cerr << "Read of file named '" << sv2[1]
	   << "' failed. Non-existent file?" << endl;
      return exc_efailed;
    }
  }
    
  if (ctype=="table") {
    
    if (sv2[1]!=((std::string)"cin")) {
      table_obj.read_generic(ifs,verbose);
    } else {
      table_obj.read_generic(std::cin,verbose);
    }

  } else if (ctype=="table3d") {
    
    if (sv2[1]!=((std::string)"cin")) {
      table3d_obj.read_gen3_list(ifs,verbose);
    } else {
      table3d_obj.read_gen3_list(std::cin,verbose);
    }
    
  } else if (ctype=="int") {

    if (sv2[1]!=((std::string)"cin")) {
      ifs >> int_obj;
    } else {
      cin >> int_obj;
    }
    
  } else if (ctype=="char") {

    if (sv2[1]!=((std::string)"cin")) {
      ifs >> char_obj;
    } else {
      cin >> char_obj;
    }
    
  } else if (ctype=="double") {

    if (sv2[1]!=((std::string)"cin")) {
      ifs >> double_obj;
    } else {
      cin >> double_obj;
    }
    
  } else if (ctype=="size_t") {

    if (sv2[1]!=((std::string)"cin")) {
      ifs >> size_t_obj;
    } else {
      cin >> size_t_obj;
    }
    
  } else if (ctype=="string") {

    if (sv2[1]!=((std::string)"cin")) {
      getline(ifs,string_obj);
    } else {
      getline(cin,string_obj);
    }
    
  } else if (ctype=="int[]") {

    if (sv2[1]!=((std::string)"cin")) {
      getline(ifs,string_obj);
      std::istringstream is(string_obj);
      int itmp;
      intv_obj.clear();
      while (is >> itmp) {
	intv_obj.push_back(itmp);
      }
      string_obj.clear();
    } else {
      getline(cin,string_obj);
      std::istringstream is(string_obj);
      int itmp;
      intv_obj.clear();
      while (is >> itmp) {
	intv_obj.push_back(itmp);
      }
      string_obj.clear();
    }
    
  } else if (ctype=="double[]") {

    if (sv2[1]!=((std::string)"cin")) {
      getline(ifs,string_obj);
      std::istringstream is(string_obj);
      double dtmp;
      doublev_obj.clear();
      while (is >> dtmp) {
	doublev_obj.push_back(dtmp);
      }
      string_obj.clear();
    } else {
      getline(cin,string_obj);
      std::istringstream is(string_obj);
      double dtmp;
      doublev_obj.clear();
      while (is >> dtmp) {
	doublev_obj.push_back(dtmp);
      }
      string_obj.clear();
    }
    
  } else if (ctype=="size_t[]") {
    
    if (sv2[1]!=((std::string)"cin")) {
      getline(ifs,string_obj);
      std::istringstream is(string_obj);
      size_t sttmp;
      size_tv_obj.clear();
      while (is >> sttmp) {
	size_tv_obj.push_back(sttmp);
      }
      string_obj.clear();
    } else {
      getline(cin,string_obj);
      std::istringstream is(string_obj);
      size_t sttmp;
      size_tv_obj.clear();
      while (is >> sttmp) {
	size_tv_obj.push_back(sttmp);
      }
      string_obj.clear();
    }
    
  } else if (ctype=="string[]") {

    size_t n;
    if (sv2[1]!=((std::string)"cin")) {
      ifs >> n;
      stringv_obj.resize(n);
      for(size_t i=0;i<n;i++) {
	getline(ifs,stringv_obj[i]);
      }
    } else {
      cin >> n;
      stringv_obj.resize(n);
      for(size_t i=0;i<n;i++) {
	getline(cin,stringv_obj[i]);
      }
    }
    
  } else {

    cerr << "Cannot read generic text file for object of type "
	 << ctype << endl;
    return 1;
    
  }

  command_add(ctype);
  type=ctype;
  
  if (sv2[1]!=((std::string)"cin")) {
    ifs.close();
  }

  return 0;
}

int acol_manager::comm_version(std::vector<std::string> &sv, bool itive_com) {

  cout << "\nacol: A data table viewing and processing program for O2scl."
       << endl;
  cout << (((string)" Compiled at ")+((string)__TIME__)+" on "+
	   ((string)__DATE__)+" for "+((string)PACKAGE)+", version "+
	   ((string)VERSION)+".\n") << endl;

  cout << "O2scl version: " << o2scl_settings.o2scl_version() << endl;
  cout << "Range checking: " << o2scl_settings.range_check() << endl;
  cout << "EOS library: " << o2scl_settings.eos_installed() << endl;
  cout << "Particle library: " << o2scl_settings.part_installed() << endl;
  cout << "HDF support: " << o2scl_settings.hdf_support() << endl;
  cout << "HDF5 compression support: "
       << o2scl_settings.hdf5_compression_support() << endl;
  cout << "Armadillo support: " << o2scl_settings.armadillo_support() << endl;
  cout << "Eigen support: " << o2scl_settings.eigen_support() << endl;
  cout << "GSL V2.0+ support: " << o2scl_settings.gsl2_support() << endl;
  cout << "OpenMP support: " << o2scl_settings.openmp_support() << endl;
  cout << "Data directory: " << o2scl_settings.get_data_dir() << endl;
  cout << endl;
  cout << "o2scl_name: " << o2scl_settings.o2scl_name() << endl;
  cout << "o2scl_package: " << o2scl_settings.o2scl_package() << endl;
  cout << "o2scl_bugreport: " << o2scl_settings.o2scl_bugreport() << endl;
  cout << "o2scl_string: " << o2scl_settings.o2scl_string() << endl;
  cout << "o2scl_tarname: " << o2scl_settings.o2scl_tarname() << endl;
  cout << endl;
  cout << "config.h: " << endl;
  o2scl_settings.config_h_report();
  return 0;
}

int acol_manager::comm_assign(std::vector<std::string> &sv, bool itive_com) {
  if (type!="table" && type!="table3d") {
    cerr << "No table/table3d object to add a constant to." << endl;
    return exc_efailed;
  }

  if (sv.size()==2) {
    if (verbose>0) {
      cout << "Removing constant named '" << sv[1] << "'." << endl;
    }
    if (type=="table3d") {
      table3d_obj.remove_constant(sv[1]);
    } else {
      table_obj.remove_constant(sv[1]);
    }
    return 0;
  }

  vector<string> pr, in;
  pr.push_back("Name of constant");
  pr.push_back("Value");
  int ret=get_input(sv,pr,in,"assign",itive_com);
  if (ret!=0) return ret;
  
  if (type=="table3d") {
    table3d_obj.add_constant(sv[1],function_to_double(sv[2]));
  } else {
    table_obj.add_constant(sv[1],function_to_double(sv[2]));
  }

  return ret;
}

int acol_manager::comm_list(std::vector<std::string> &sv, bool itive_com) {
  cout.precision(prec);
  if (type=="table3d") {
    cout << "Table3d name: " << obj_name << endl;
    table3d_obj.summary(&cout,ncols);
  } else if (type=="table") {
    cout << "Table name: " << obj_name << endl;
    if (table_obj.get_nunits()>0) {
      table_obj.summary(&cout,ncols);
    } else {
      table_obj.table<std::vector<double> >::summary(&cout,ncols);
    }
  } else {
    cerr << "No table to list columns for." << endl;
    return exc_efailed;
  }
  return 0;
}

int acol_manager::comm_sort(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
  
    std::string i1, i2;
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to sort." << endl;
      return exc_efailed;
    }
    
    if (sv.size()>2) {
      i1=sv[1];
      i2=sv[2];
    } else if (sv.size()>1) {
      i1=sv[1];
    } else {
      if (itive_com) {
	i1=cl->cli_gets("Enter column to sort by (or blank to stop): ");
	if (i1.length()==0) {
	  cout << "Command 'sort' cancelled." << endl;
	  return 0;
	}
      } else {
	cerr << "Not enough arguments for 'sort'." << endl;
	return exc_efailed;
      }
    }
    
    bool unique=false;
    if (i2==((std::string)"unique")) unique=true;
    
    if (table_obj.is_column(i1)==false) {
      cerr << "Could not find column named '" << i1 << "'." << endl;
      return exc_efailed;
    }
    
    if (verbose>1) {
      cout << "Sorting by column " << i1 << endl; 
    }
    table_obj.sort_table(i1);
    
    if (unique) {
      std::cout << "Going to dir." << std::endl;
      table_obj.delete_idadj_rows();
      std::cout << "Done in dir." << std::endl;
    }

  } else if (type=="double[]") {

    vector_sort<vector<double>,double>(doublev_obj.size(),doublev_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }
    
  } else if (type=="int[]") {

    vector_sort<vector<int>,int>(intv_obj.size(),intv_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }
    
  } else if (type=="size_t[]") {

    vector_sort<vector<size_t>,size_t>(size_tv_obj.size(),size_tv_obj);
    if (verbose>0) {
      cout << "Object of type double[] sorted." << endl;
    }

  } else {

    cout << "Not implemented for type " << type << endl;
    return 1;
    
  }
  
  return 0;
}

int acol_manager::comm_stats(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }
  
  std::string i1;

  if (table_obj.get_nlines()==0) {
    cerr << "No table to analyze." << endl;
    return exc_efailed;
  }
  
  if (sv.size()>1) {
    i1=sv[1];
  } else {
    if (itive_com) {
      i1=cl->cli_gets("Enter column to get info on (or blank to stop): ");
      if (i1.length()==0) {
	cout << "Command 'stats' cancelled." << endl;
	return 0;
      }
    } else {
      cerr << "Not enough arguments for 'stats'." << endl;
      return exc_efailed;
    }
  }

  if (table_obj.is_column(i1)==false) {
    cerr << "Could not find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  const vector<double> &cref=table_obj.get_column(i1);
  cout << "N        : " << table_obj.get_nlines() << endl;
  cout << "Sum      : " << vector_mean(table_obj.get_nlines(),cref)*
    table_obj.get_nlines() << endl;
  cout << "Mean     : " << vector_mean(table_obj.get_nlines(),cref) << endl;
  cout << "Std. dev.: " << vector_stddev(table_obj.get_nlines(),cref) << endl;
  size_t ix;
  double val;
  vector_min(table_obj.get_nlines(),cref,ix,val);
  cout << "Min      : " << val << " at index: " << ix << endl;
  vector_max(table_obj.get_nlines(),cref,ix,val);
  cout << "Max      : " << val << " at index: " << ix << endl;

  size_t dup=0, inc=0, dec=0;
  for(size_t i=0;i<table_obj.get_nlines()-1;i++) {
    if (cref[i+1]==cref[i]) dup++;
    if (cref[i]<cref[i+1]) inc++;
    if (cref[i]>cref[i+1]) dec++;
  }
  if (inc>0 && dec==0) {
    if (dup>0) {
      cout << "Increasing (" << dup << " duplicates)." << endl;
    } else {
      cout << "Strictly increasing. No duplicates." << endl;
    }
  } else if (dec>0 && inc==0) {
    if (dup>0) {
      cout << "Decreasing (" << dup << " duplicates)." << endl;
    } else {
      cout << "Strictly decreasing. No duplicates." << endl;
    }
  } else if (dec==0 && inc==0) {
    cout << "Constant (" << dup << " duplicates)." << endl;
  } else {
    cout << "Non-monotonic (" << inc << " increasing, " << dec 
	 << " decreasing, and " << dup << " duplicates)." << endl;
  }
  if ((dup+inc+dec)!=(table_obj.get_nlines()-1)) {
    cout << "Counting mismatch from non-finite values or signed zeros." << endl;
  }
  
  
  return 0;
}

int acol_manager::comm_set(std::vector<std::string> &sv, bool itive_com) {

  // This is currently taken care of inside cli

  return 0;
}

int acol_manager::comm_help(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()==3) {
    string temp_type=sv[1];
    string cur_type=type;

    command_del();
    command_add(temp_type);
    
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
    
    int ret=cl->comm_option_help(sv,itive_com);

    command_del();
    command_add(cur_type);
    return ret;
  }
  
  return cl->comm_option_help(sv,itive_com);
}

int acol_manager::comm_commands(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()==2) {
    cout << "Commands argument: " << sv[1] << endl;
    string temp_type=sv[1];
    string cur_type=type;

    command_del();
    command_add(temp_type);
    
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
    int ret=cl->comm_option_commands(sv,itive_com);

    command_del();
    command_add(cur_type);
    return ret;
  }
  return cl->comm_option_commands(sv,itive_com);
}

int acol_manager::comm_select(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table" && type!="table3d") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }
  
  // Remove outermost single or double quotes from all arguments
  for(size_t i=0;i<sv.size();i++) {
    size_t svs=sv[i].size();
    if (svs>2 && ((sv[i][0]=='\'' && sv[i][svs-1]=='\'') ||
		  (sv[i][0]=='\"' && sv[i][svs-1]=='\"'))) {
      sv[i]=sv[i].substr(1,sv[i].size()-2);
    }
  }

  // ----------------------------------------------------------------
  // Parse arguments into names and values
  // ----------------------------------------------------------------

  // The number of column arguments
  int nargs=((int)sv.size())-1;
  // The name part of the argument
  std::vector<std::string> names;
  // The function part of the argument
  std::vector<std::string> args;
  // Each element is true if the corresponding argument is a pattern
  std::vector<bool> is_pattern;
  
  // ----------------------------------------------------------------
  
  if (nargs==0) {
    cerr << "No arguments to 'select'. Command 'select' aborted." << endl;
    return exc_efailed;
  }
  
  // Reserve enough space
  names.reserve(nargs);
  args.reserve(nargs);
  is_pattern.reserve(nargs);
  
  // Go through each column
  for(int i=0;i<nargs;i++) {

    args.push_back(sv[i+1]);
    
    if (args[i][0]==':') {
      
      // If it has a colon, add to the patterns list
      args[i]=args[i].substr(1,args[i].length()-1);
      is_pattern.push_back(true);
      names.push_back("");
      
    } else {
      
      // If there's no colon, look for an '=' sign
      int tmpindx;
      tmpindx=args[i].find('=',0);
      
      if (tmpindx!=((int)string::npos)) {
	
	// Parse column name and function
	names.push_back(args[i].substr(0,tmpindx));
	args[i]=args[i].substr(tmpindx+1,args[i].size()-tmpindx-1);
	
      } else {
	
	// There's no equals sign, just assume its a column name
	names.push_back(args[i]);
      }
      
      is_pattern.push_back(false);
      
    }
  }
  
  // Output the results of the parsing
  if (verbose>2) {
    cout << "Parsed arguments: " << endl;
    for(int i=0;i<nargs;i++) {
      cout.width(2);
      cout << i << " ";
      if (is_pattern[i]==true) {
	cout << "Pattern: " << args[i] << endl;
      } else {
	cout << names[i] << " = " << args[i] << endl;
      }
    }
  }

  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d to select from." << endl;
      return exc_efailed;
    }
    
    // ---------------------------------------------------------------------
    // Create new table3d and copy grid over
    // ---------------------------------------------------------------------

    table3d *new_table3d=new table3d;
    size_t nx, ny;
    table3d_obj.get_size(nx,ny);
    new_table3d->set_xy(table3d_obj.get_x_name(),nx,table3d_obj.get_x_data(),
			table3d_obj.get_y_name(),ny,table3d_obj.get_y_data());
	
    // ---------------------------------------------------------------------
    // Copy constants from old to new table3d
    // ---------------------------------------------------------------------

    for(size_t i=0;i<table3d_obj.get_nconsts();i++) {
      string tnam;
      double tval;
      table3d_obj.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      new_table3d->add_constant(tnam,tval);
    }
  
    // ---------------------------------------------------------------------
    // Add slides and data to new table3d
    // ---------------------------------------------------------------------

    std::vector<bool> matched(table3d_obj.get_nslices());
    for(size_t i=0;i<table3d_obj.get_nslices();i++) {
      matched[i]=false;
    }

    // In this next loop, we need fix the code to ensure that when a
    // non-pattern column is given, it sets matched to true
    // appropriately.

    // Copy the data over
    for(int i=0;i<nargs;i++) {

      if (is_pattern[i]==false) {
	
        // Return an error if the slice doesn't exist
	size_t ix;
        if (names[i]==args[i] && table3d_obj.is_slice(args[i],ix)==false) {
          cerr << "Slice '" << args[i] << "' is not in the table." << endl;
          return exc_einval;
        }

        // Add the new slice to the new table
        new_table3d->new_slice(names[i]);

        // Fill slice with the new data
        ubmatrix mat(nx,ny);
        table3d_obj.function_matrix(args[i],mat,false);
        new_table3d->copy_to_slice(mat,names[i]);

      } else {
	
        // Find the matching slices
        for(size_t j=0;j<table3d_obj.get_nslices();j++) {
	  
          if (matched[j]==false &&  
              fnmatch(args[i].c_str(),
                      table3d_obj.get_slice_name(j).c_str(),0)==0) {
	    
            // If we've found a match, add it to the new table
            matched[j]=true;
	    
            // Add the new slice to the new table
            new_table3d->new_slice(table3d_obj.get_slice_name(j));
	    
            // Fill it with the new data
	    ubmatrix mat(nx,ny);
	    table3d_obj.function_matrix(args[i],mat,false);
	    new_table3d->copy_to_slice(mat,table3d_obj.get_slice_name(j));
          }
        }
      }
    }
    
    // Todo: Replace this copy with std::swap
    table3d_obj=*new_table3d;
    
  } else {

    if (table_obj.get_nlines()==0) {
      cerr << "No table to select from." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Create new table
    // ---------------------------------------------------------------------
  
    table_units<> new_table;
  
    new_table.set_nlines(table_obj.get_nlines());

    // ---------------------------------------------------------------------
    // Copy constants from old to new table
    // ---------------------------------------------------------------------

    for(size_t i=0;i<table_obj.get_nconsts();i++) {
      string tnam;
      double tval;
      table_obj.get_constant(i,tnam,tval);
      new_table.add_constant(tnam,tval);
    }
  
    // ---------------------------------------------------------------------
    // Add columns and data to new table
    // ---------------------------------------------------------------------

    std::vector<bool> matched(table_obj.get_ncolumns());
    for(size_t i=0;i<table_obj.get_ncolumns();i++) {
      matched[i]=false;
    }

    // In this next loop, we need fix the code to ensure that when a
    // non-pattern column is given, it sets matched to true
    // appropriately.

    for(int i=0;i<nargs;i++) {

      if (is_pattern[i]==false) {
      
	// Return an error if the column doesn't exist
	if (names[i]==args[i] && table_obj.is_column(args[i])==false) {
	  cerr << "Column '" << args[i] << "' is not in the table." << endl;
	  return exc_einval;
	}

	// Add the new column to the new table
	new_table.new_column(names[i]);

	// If necessary, set units
	if (names[i]==args[i] && table_obj.get_unit(args[i]).length()!=0) {
	  new_table.set_unit(names[i],table_obj.get_unit(args[i]));
	}

	// Fill column with the new data
	ubvector vec(table_obj.get_nlines());
	table_obj.function_vector(args[i],vec,false);
	new_table.copy_to_column(vec,names[i]);

      } else {

	// Find the matching columns
	for(size_t j=0;j<table_obj.get_ncolumns();j++) {

	  if (matched[j]==false &&  
	      fnmatch(args[i].c_str(),
		      table_obj.get_column_name(j).c_str(),0)==0) {

	    // If we've found a match, add it to the new table
	    matched[j]=true;

	    // Add the new column to the new table
	    new_table.new_column(table_obj.get_column_name(j));

	    // If necessary, set units
	    string tmp=table_obj.get_column_name(j);
	    if (table_obj.get_unit(tmp).length()!=0) {
	      new_table.set_unit(tmp,table_obj.get_unit(tmp));
	    }

	    // Fill it with the new data
	    ubvector vec(table_obj.get_nlines());
	    table_obj.function_vector(table_obj.get_column_name(j),vec,false);
	    new_table.copy_to_column(vec,table_obj.get_column_name(j));
	  }
	}
      }
    }

    // Todo: Replace this copy with std::swap
    table_obj=new_table;

  }

  // Call list command
  if (verbose>0) {
    comm_list(sv,itive_com);
  }

  return 0;
}

int acol_manager::comm_delete_rows(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to delete rows from." << endl;
    return exc_efailed;
  }

  std::string i1;
  if (sv.size()>=2) {
    i1=sv[1];
  } else if (itive_com==true) {
    i1=cl->cli_gets("Function to specify rows (or blank to quit): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'delete_rows' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No rows to delete." << endl;
    return exc_efailed;
  }
  
  table_obj.delete_rows(i1);

  return 0;
}

int acol_manager::comm_select_rows(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to select rows from." << endl;
    return exc_efailed;
  }

  std::string i1;
  if (sv.size()>=2) {
    i1=sv[1];
  } else if (itive_com==true) {
    i1=cl->cli_gets("Function to specify rows (or blank to quit): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'select_rows' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No rows to delete." << endl;
    return exc_efailed;
  }
    
  // ---------------------------------------------------------------------
  // Create new table
  // ---------------------------------------------------------------------
  
  table_units<> *new_table=new table_units<>;
  
  // ---------------------------------------------------------------------
  // Copy constants from old to new table
  // ---------------------------------------------------------------------

  for(size_t i=0;i<table_obj.get_nconsts();i++) {
    string tnam;
    double tval;
    table_obj.get_constant(i,tnam,tval);
    new_table->add_constant(tnam,tval);
  }
  
  // ---------------------------------------------------------------------
  // Add column names to new table
  // ---------------------------------------------------------------------

  for(int i=0;i<((int)table_obj.get_ncolumns());i++) {
    new_table->new_column(table_obj.get_column_name(i));
  }

  // ---------------------------------------------------------------------
  // Copy data from selected rows
  // ---------------------------------------------------------------------

  int new_lines=0;
  for(int i=0;i<((int)table_obj.get_nlines());i++) {
    if (table_obj.row_function(i1,i)>0.5) {
      new_table->set_nlines(new_lines+1);
      for(int j=0;j<((int)table_obj.get_ncolumns());j++) {
	new_table->set(j,new_lines,table_obj.get(j,i));
      }
      new_lines++;
    }
  }
  
  // Replace the old table with the new one
  table_obj.clear();
  table_obj=*new_table;
  delete new_table;

  return 0;
}

int acol_manager::comm_select_rows2(std::vector<std::string> &sv, 
				    bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to select rows from." << endl;
    return exc_efailed;
  }

  std::string i1;
  if (sv.size()>=2) {
    i1=sv[1];
  } else if (itive_com==true) {
    i1=cl->cli_gets("Function to specify rows (or blank to quit): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'select_rows' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No rows to delete." << endl;
    return exc_efailed;
  }
    
  // ---------------------------------------------------------------------
  // Create new table
  // ---------------------------------------------------------------------
  
  table_units<> *new_table=new table_units<>;
  
  // ---------------------------------------------------------------------
  // Copy constants from old to new table
  // ---------------------------------------------------------------------

  for(size_t i=0;i<table_obj.get_nconsts();i++) {
    string tnam;
    double tval;
    table_obj.get_constant(i,tnam,tval);
    new_table->add_constant(tnam,tval);
  }
  
  // ---------------------------------------------------------------------
  // Add column names to new table
  // ---------------------------------------------------------------------

  for(int i=0;i<((int)table_obj.get_ncolumns());i++) {
    new_table->new_column(table_obj.get_column_name(i));
  }

  // ---------------------------------------------------------------------
  // Copy data from selected rows
  // ---------------------------------------------------------------------

  calculator calc;
  std::map<std::string,double> vars;

  vector<string> cols;
  for(size_t i=2;i<sv.size();i++) {
    cols.push_back(sv[i]);
  }
  
  int new_lines=0;
  for(int i=0;i<((int)table_obj.get_nlines());i++) {
    if (i%10000==0) std::cout << "I: " << i << endl;
    for(size_t j=0;j<cols.size();j++) {
      vars[cols[j]]=table_obj.get(cols[j],i);
    }
    calc.compile(i1.c_str(),&vars);
    if (calc.eval(&vars)>0.5) {
      new_table->set_nlines(new_lines+1);
      for(int j=0;j<((int)table_obj.get_ncolumns());j++) {
	new_table->set(j,new_lines,table_obj.get(j,i));
      }
      new_lines++;
    }
  }
  
  // Replace the old table with the new one
  table_obj.clear();
  table_obj=*new_table;
  delete new_table;
  
  return 0;
}

int acol_manager::comm_delete_col(std::vector<std::string> &sv, 
				  bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to delete columns from." << endl;
    return exc_efailed;
  }

  std::string i1;
  if (sv.size()>=2) {
    i1=sv[1];
  } else if (itive_com==true) {
    i1=cl->cli_gets("Column to delete: ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'delete_col' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No rows to delete." << endl;
    return exc_efailed;
  }
    
  if (table_obj.is_column(i1)==false) {
    cerr << "Could not find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  table_obj.delete_column(i1);

  return 0;
}

int acol_manager::comm_create(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype, tval;

  // Delete previous object
  command_del();
  clear_obj();
  
  int ret=get_input_one(sv,"Enter type of object to create",ctype,"create",
			itive_com);
  if (ret!=0) return ret;

  vector<string> sv2=sv;
  vector<string>::iterator it=sv2.begin();
  sv2.erase(it+1);
  
  if (ctype=="int") {

    int ret=get_input_one(sv2,"Enter integer",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    int_obj=o2scl::stoi(tval);
    type="int";
    command_add("int");
    obj_name="int";
    
  } else if (ctype=="size_t") {

    int ret=get_input_one(sv2,"Enter size_t",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    size_t_obj=o2scl::stoszt(tval);
    type="size_t";
    command_add("size_t");
    obj_name="size_t";
    
  } else if (ctype=="char") {

    int ret=get_input_one(sv2,"Enter char",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    char_obj=tval[0];
    type="char";
    command_add("char");
    obj_name="char";
    
  } else if (ctype=="double") {

    int ret=get_input_one(sv2,"Enter double",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    double_obj=o2scl::function_to_double(tval);
    type="double";
    command_add("double");
    obj_name="double";

  } else if (ctype=="string") {

    int ret=get_input_one(sv2,"Enter string",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    string_obj=tval;
    type="string";
    command_add("string");
    obj_name="string";
    
  } else if (ctype=="double[]") {
    
    vector<string> in, pr;
    pr.push_back("Size");
    pr.push_back("Function of i (starting with zero)");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    calculator calc;
    std::map<std::string,double> vars;
    std::map<std::string,double>::const_iterator mit;
    size_t nn=o2scl::stoszt(in[0]);
    doublev_obj.clear();
    calc.compile(in[1].c_str(),&vars);
    for(size_t i=0;i<nn;i++) {
      vars["i"]=((double)i);
      doublev_obj.push_back(calc.eval(&vars));
    }
    command_add("double[]");
    type="double[]";
    
  } else if (ctype=="int[]") {
    
    vector<string> in, pr;
    pr.push_back("Size");
    pr.push_back("Function of i (starting with zero)");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    calculator calc;
    std::map<std::string,double> vars;
    std::map<std::string,double>::const_iterator mit;
    size_t nn=o2scl::stoszt(in[0]);
    intv_obj.clear();
    calc.compile(in[1].c_str(),&vars);
    for(size_t i=0;i<nn;i++) {
      vars["i"]=((double)i);
      intv_obj.push_back(((int)(calc.eval(&vars))));
    }
    command_add("int[]");
    type="int[]";
    
  } else if (ctype=="size_t[]") {
    
    vector<string> in, pr;
    pr.push_back("Size");
    pr.push_back("Function of i (starting with zero)");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    calculator calc;
    std::map<std::string,double> vars;
    std::map<std::string,double>::const_iterator mit;
    size_t nn=o2scl::stoszt(in[0]);
    size_tv_obj.clear();
    calc.compile(in[1].c_str(),&vars);
    for(size_t i=0;i<nn;i++) {
      vars["i"]=((double)i);
      size_tv_obj.push_back(((size_t)(calc.eval(&vars))));
    }
    command_add("size_t[]");
    type="size_t[]";
    
  } else if (ctype=="table") {
    
    vector<string> in, pr;
    pr.push_back("Name of new column");
    pr.push_back("Value for first row");
    pr.push_back("Maximum value");
    pr.push_back("Increment");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;
    
    double d2=function_to_double(in[1]);
    double d3=function_to_double(in[2]);
    double d4=function_to_double(in[3]);
    d3+=d4/1.0e4;
    int cnl=((int)((d3-d2)/d4))+1;
    
    table_obj.clear();
    table_obj.line_of_names(in[0]);
    table_obj.set_nlines(cnl);
    
    for(int li=0;li<cnl;li++) {
      table_obj.set(in[0],li,d2+((double)li)*d4);
    }
    command_add("table");
    type="table";
    
  } else if (ctype=="table3d") {
    
    vector<string> in;
    vector<string> pr=
      {"Enter x-axis name (or blank to stop): ",
       "Enter x-axis lower limit (or blank to stop): ",
       "Enter x-axis upper limit (or blank to stop): ",
       "Enter x-axis step size (or blank to stop): ",
       "Enter y-axis name (or blank to stop): ",
       "Enter y-axis lower limit (or blank to stop): ",
       "Enter y-axis upper limit (or blank to stop): ",
       "Enter y-axis step size (or blank to stop): ",
       "Enter slice name (or blank to stop): ",
       "Enter slice function (or blank to stop): "};
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;
    
    std::string xname=in[0];
    double x0=function_to_double(in[1]);
    double x1=function_to_double(in[2]);
    double dx=function_to_double(in[3]);
    uniform_grid_end<double> ugx(x0,x1,((size_t)(x1-x0)/(dx*(1.0-1.0e-14))));
    
    std::string yname=in[4];
    double y0=function_to_double(in[5]);
    double y1=function_to_double(in[6]);
    double dy=function_to_double(in[7]);
    uniform_grid_end<double> ugy(y0,y1,((size_t)(y1-y0)/(dy*(1.0-1.0e-14))));
    
    std::string zname=in[8];
    std::string zfunc=in[9];
    
    table3d_obj.set_xy(xname,ugx,yname,ugy);
    
    table3d_obj.function_slice(zfunc,zname);
    
    command_add("table3d");
    type="table3d";

  } else {

    cerr << "Cannot create object of type " << ctype << endl;
    return 1;
      
  }

  return 0;
}

int acol_manager::comm_insert(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table" && type!="table3d") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (type=="table3d") {

    std::string in[4], pr[4]=
      {"Enter filename of external table (or blank to stop): ",
       "Enter name of table in file (or blank for first table): ",
       "Enter slice in external table (or blank to stop): ",
       "Enter name of new slice in present table (blank to keep old name): "};

    if (sv.size()>3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
      if (sv.size()>4) {
	in[3]=sv[4];
      } else {
	in[3]="";
      }
    } else {
      if (itive_com) {
	for(size_t is=0;is<4;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (is!=3 && is!=1 && in[is].length()==0) {
	    cout << "Command 'insert' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to 'insert'" << endl;
	return exc_efailed;
      }
    }

    if (true || verbose>2) {
      cout << "Read table3d  named " << in[1] << " from file " << in[0] << endl;
      cout << "old slice, new slice: " << in[2] << " " << in[3] << endl;
    }

    // Read table from file
    hdf_file hf;
    table3d tmp;
    int hfret=hf.open(in[0],false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << in[0] << endl;
      return exc_efailed;
    }
    
    std::string tmp_name;
    if (in[1].length()>0) tmp_name=in[1];
    hdf_input(hf,tmp,tmp_name);
    hf.close();

    table3d_obj.add_slice_from_table(tmp,in[2],in[3]);

    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into." << endl;
    return exc_efailed;
  }

  std::string in[6], pr[6]=
    {"Enter filename of external table (or blank to stop): ",
     "Enter name of table in file (or blank for first table): ",
     "Enter index column in external table (or blank to stop): ",
     "Enter data column in external table (or blank to stop): ",
     "Enter index column in present table (or blank to stop): ",
     "Enter name of new column in present table (or blank to keep old name): "};
  if (sv.size()>5) {
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
    in[3]=sv[4];
    in[4]=sv[5];
    if (sv.size()>6) {
      in[5]=sv[6];
    } else {
      in[5]="";
    }
  } else {
    if (itive_com) {
      for(size_t is=0;is<6;is++) {
	in[is]=cl->cli_gets(pr[is].c_str());
	if (is!=5 && is!=1 && in[is].length()==0) {
	  cout << "Command 'insert' cancelled." << endl;
	  return 0;
	}
      }
    } else {
      cerr << "Not enough arguments to 'insert'" << endl;
      return exc_efailed;
    }
  }

  if (true || verbose>2) {
    cout << "Read table named " << in[1] << " from file " << in[0] << endl;
    cout << "oldx,oldy,newx,newy: " << in[2] << " " << in[3] << " " 
	 << in[4] << " " << in[5] << endl;
    cout << endl;
  }

  // Read table from file
  hdf_file hf;
  table_units<> tmp;
  int hfret=hf.open(in[0],false,false);
  if (hfret!=0) {
    cerr << "Failed to read file named " << in[0] << endl;
    return exc_efailed;
  }
  std::string tmp_name;
  if (in[1].length()>0) tmp_name=in[1];
  hdf_input(hf,tmp,tmp_name);
  hf.close();

  table_obj.add_col_from_table(tmp,in[2],in[3],in[4],in[5]);

  return 0;

}

int acol_manager::comm_insert_full(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (type!="table") {
    cout << "Not implemented for type " << type << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table to insert columns into in command 'insert-full'." << endl;
    return exc_efailed;
  }

  std::string in[3], pr[3]=
    {"Enter filename of external table (or blank to stop): ",
     "Enter index column in present table (or blank to stop): ",
     "Enter index column in external table (or blank to stop): "};
  if (sv.size()>=3) {
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
  } else {
    if (itive_com) {
      for(size_t is=0;is<3;is++) {
	in[is]=cl->cli_gets(pr[is].c_str());
	if (in[is].length()==0) {
	  cout << "Command 'insert' cancelled." << endl;
	  return 0;
	}
      }
    } else {
      cerr << "Not enough arguments to command 'insert'" << endl;
      return exc_efailed;
    }
  }

  cout << "Unimplemented." << endl;

  return 1;

#ifdef O2SCL_NEVER_DEFINED

  for (size_t j=0;j<tmp->get_ncolumns();j++) {
    string ty=tmp->get_column_name(j);
    int tret=table_obj.add_col_from_table(in[1],*tmp,in[2],ty,ty);
    if (tret!=0) {
      cerr << "Adding column " << ty << " failed." << endl;
      ret=tret;
      // We don't return here so that "delete tmp;" is called below
    }
  }
  
#endif
}

int acol_manager::comm_interp(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table3d") {
    
    // --------------------------------------------------------------
    // 3d table interpolation
    
    if (type!="table3d") {
      cerr << "No table to interpolate into." << endl;
      return exc_efailed;
    }

    std::string in[3], pr[3]=
      {"Enter slice name (or blank to stop): ",
       "Enter x value (or blank to stop): ",
       "Enter y value (or blank to stop): "};
    if (sv.size()>=3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
    } else {
      if (itive_com) {
	for(size_t is=0;is<3;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (in[is].length()==0) {
	    cout << "Command 'interp' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to command 'interp'" << endl;
	return exc_efailed;
      }
    }
  
    double ret=table3d_obj.interp(function_to_double(in[1]),
				  function_to_double(in[2]),in[0]);
    if (err_hnd->get_errno()!=0) {
      cerr << "Interpolation failed." << endl;
      return exc_efailed;
    } else {
      cout << "Interpolation result: " << ret << endl;
    }

    return 0;

  } else if (type=="table") {

    // --------------------------------------------------------------
    // 2d table interpolation
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to interpolate into." << endl;
      return exc_efailed;
    }
    
    std::string in[3], pr[3]=
      {"Enter column name of independent variable (or blank to stop): ",
       "Enter value of independent variable (or blank to stop): ",
       "Enter column name of dependent variable (or blank to stop): "};
    if (sv.size()>=3) {
      in[0]=sv[1];
      in[1]=sv[2];
      in[2]=sv[3];
    } else {
      if (itive_com) {
	for(size_t is=0;is<3;is++) {
	  in[is]=cl->cli_gets(pr[is].c_str());
	  if (in[is].length()==0) {
	    cout << "Command 'interp' cancelled." << endl;
	    return 0;
	  }
	}
      } else {
	cerr << "Not enough arguments to 'interp'" << endl;
	return exc_efailed;
      }
    }
    
    if (table_obj.is_column(in[0])==false) {
      cerr << "Could not find column named '" << in[0] << "'." << endl;
      return exc_efailed;
    }
    if (table_obj.is_column(in[2])==false) {
      cerr << "Could not find column named '" << in[2] << "'." << endl;
      return exc_efailed;
    }
    
    double ret=table_obj.interp(in[0],function_to_double(in[1]),in[2]);
    if (err_hnd->get_errno()!=0) {
      cerr << "Interpolation failed." << endl;
      return exc_efailed;
    } else {
      cout << "Interpolation result: " << ret << endl;
    }

    // --------------------------------------------------------------
    
  } else if (type=="double[]") {

    size_t n=doublev_obj.size();
    std::vector<double> index(n);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,doublev_obj) << endl;

  } else if (type=="int[]") {

    size_t n=intv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(intv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,value) << endl;

  } else if (type=="size_t[]") {

    size_t n=size_tv_obj.size();
    std::vector<double> index(n), value(n);
    o2scl::vector_copy(size_tv_obj,value);
    for(size_t i=0;i<n;i++) index[i]=((double)i);

    o2scl::interp<vector<double> > it(interp_type);
    double x=o2scl::stod(sv[1]);
    cout << "Interpolation result: "
	 << it.eval(x,n,index,value) << endl;

  } else {
    cout << "Not implemented for type " << type << endl;
    return 1;
  }    
  
  return 0;
}

int acol_manager::comm_fit(std::vector<std::string> &sv, bool itive_com) {

  cout << "Not implemented." << endl;
  return 0;

#ifdef O2SCL_NEVER_DEFINED
  if (type=="table3d") {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (table_obj.get_nlines()==0) {
    cerr << "No table with data to fit." << endl;
    return exc_efailed;
  }

  const size_t nargs=7;
  std::string in[nargs], pr[nargs]=
    {"Enter column name of independent variable (or blank to stop): ",
     "Enter column name of dependent variable (or blank to stop): ",
     "Enter column name of errors on dependent variable (or blank to stop): ",
     "Enter column name of fitted values (or blank to stop): ",
     "Enter comma-delimited parameter name list (or blank to stop): ",
     "Enter function to fit (or blank to stop): ",
     "Enter space-delimited list of initial values (or blank to stop): "};
  if (sv.size()>=nargs) {
    in[0]=sv[1];
    in[1]=sv[2];
    in[2]=sv[3];
    in[3]=sv[4];
    in[4]=sv[5];
    in[5]=sv[6];
    in[6]=sv[7];
  } else {
    if (itive_com) {
      for(size_t is=0;is<nargs;is++) {
	in[is]=cl->cli_gets(pr[is].c_str());
	if (in[is].length()==0) {
	  cout << "Command 'fit' cancelled." << endl;
	  return 0;
	}
      }
    } else {
      cerr << "Not enough arguments to 'fit'" << endl;
      return exc_efailed;
    }
  }
  
  // Create data to fit to
  size_t ndat=table_obj.get_nlines();
  ubvector xdat(ndat), ydat(ndat), yerr(ndat);
  for(size_t i=0;i<ndat;i++) {
    xdat[i]=table_obj.get(in[0],i);
    ydat[i]=table_obj.get(in[1],i);
    yerr[i]=table_obj.get(in[2],i);
  }
  if (verbose>=2) {
    cout << "Data summary:" << endl;
    for(size_t i=0;i<ndat;i+=ndat/10) {
      cout.width(4);
      cout << i << " ";
      cout << xdat[i] << " " << ydat[i] << " " << yerr[i] << endl;
    }
  }

  // Parse initial parameter values
  std::vector<std::string> param_list;
  split_string(in[6],param_list);
  size_t n_parms=param_list.size();
  ubvector params(n_parms);
  for(size_t k=0;k<n_parms;k++) {
    params[k]=o2scl::stod(param_list[k]);
  }
  if (verbose>=1) {
    cout << "Initial parameters: ";
    for(size_t k=0;k<n_parms;k++) {
      cout << params[k] << " ";
    }
    cout << endl;
    cout << "Function: " << in[5] << endl;
  }

  // Set up fitting function
  fit_funct_strings ffs(in[5],in[4],in[0]);
  ffs.set_aux_parms(params);

  // Fitting function object
  chi_fit_funct<ubvector,ubmatrix,fit_funct_strings<> > 
    cff(ndat,xdat,ydat,yerr,ffs);
  
  // Perform fit
  double chi2;
  ubmatrix covar(n_parms,n_parms);
  fit_nonlin<> gf;
  gf.err_nonconv=false;
  int ret=gf.fit(n_parms,params,covar,chi2,cff);

  if (ret!=0) {
    cout << "Fit failed." << endl;
    return exc_einval;
  }

  // Create and fill the new fitted value column
  table_obj.new_column(in[3]);
  for(int k=0;k<((int)table_obj.get_nlines());k++) {
    table_obj.set(in[3],k,ffs(n_parms,params,(xdat)[k]));
  }
  
  // Output results to cout
  cout << "\nFit results: " << endl;
  for(int k=0;k<((int)n_parms);k++) {
    cout << "Parameter " << k+1 << ": ";
    if (params[k]<0.0) cout << params[k] << endl;
    else cout << " " << params[k] << endl;
  }
  cout << "Covariance matrix:" << endl;
  matrix_out(cout,n_parms,n_parms,covar);
  cout << "Chi-squared: " << chi2 << endl;
  cout << endl;
  
#endif
  
  return 0;
}
