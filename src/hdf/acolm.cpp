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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/*
  This source file contains the constructor and the functions
  command_add(), get_input(), get_input_one(), command_del(),
  clear_obj(), run(), setup_parameters(), setup_help(), setup_cli(),
  and setup_options().
*/

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

#ifdef O2SCL_HDF5_COMP
  compress=1;
#else
  compress=0;
#endif

  cng.err_on_fail=false;

  {
    vector<std::string> itmp={"value"};
    type_comm_list.insert(std::make_pair("int",itmp));
    type_comm_list.insert(std::make_pair("double",itmp));
    type_comm_list.insert(std::make_pair("char",itmp));
    type_comm_list.insert(std::make_pair("size_t",itmp));
    type_comm_list.insert(std::make_pair("string",itmp));
  }
  {
    vector<std::string> itmp={"assign","autocorr","delete-col","delete-rows",
			      "delete-rows-tol","deriv","deriv2","cat",
			      "convert-unit","find-row","fit","function",
			      "get-row","get-unit","entry","index",
			      "insert","insert-full","integ","interp",
			      "list","max","min","nlines","rename",
			      "select","select-rows","select-rows2",
			      "set-data","set-unit","sort","stats","sum",
			      "to-hist","to-hist-2d","wstats"};
    type_comm_list.insert(std::make_pair("table",itmp));
  }
  {
    vector<std::string> itmp={"cat","contours","deriv-x","deriv-y",
			      "function","entry","insert","interp",
			      "list","max","min","rename","set-data",
			      "slice","sum","x-name","y-name"};
    type_comm_list.insert(std::make_pair("table3d",itmp));
  }
  {
    vector<std::string> itmp={"list","min","max","to-table3d",
			      "rearrange"};
    type_comm_list.insert(std::make_pair("tensor<int>",itmp));
    type_comm_list.insert(std::make_pair("tensor<size_t>",itmp));
  }
  {
    vector<std::string> itmp={"list","diag","to-table3d","to-table3d-sum",
			      "max","min","to-tensor-grid","rearrange"};
    type_comm_list.insert(std::make_pair("tensor",itmp));
  }
  {
    vector<std::string> itmp={"to-table3d"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_amr",itmp));
  }
  {
    vector<std::string> itmp={"list","to-table3d","slice","to-table",
			      "set-grid","max","min"};
    type_comm_list.insert(std::make_pair("tensor_grid",itmp));
  }
  {
    vector<std::string> itmp={"max","min"};
    type_comm_list.insert(std::make_pair("hist_2d",itmp));
  }
  {
    vector<std::string> itmp={"deriv","interp","max","min","sort",
			      "autocorr","to-table","function","sum"};
    type_comm_list.insert(std::make_pair("double[]",itmp));
    type_comm_list.insert(std::make_pair("int[]",itmp));
    type_comm_list.insert(std::make_pair("size_t[]",itmp));
  }

}

void acol_manager::command_add(std::string new_type) {

  const int both=cli::comm_option_both;
  
  if (new_type=="int") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","Get or set the value of the int object",
       0,1,"[value]","Get or set the value of the int object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="double") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","Get or set the value of the double object",
       0,1,"[value]","Get or set the value of the double object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="char") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","Get or set the value of the char object",
       0,1,"[value]","Get or set the value of the char object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="size_t") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","Get or set the value of the size_t object",
       0,1,"[value]","Get or set the value of the size_t object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="string") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"value","Get or set the value of the string object",
       0,1,"[value]","Get or set the value of the string object",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="table") {
    static const size_t narr=36;
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
      {0,"delete-rows-tol","Delete rows within a tolerance.",
       0,2,"[relative tol.] [absolute tol.]",
       ((std::string)("This command deletes all rows which match "))+
       "within the specified tolerances. If verbose is larger than zero "+
       "then information about how many rows were deleted is provided.",
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
      {0,"find-row","Find a row which maximizes a function.",
       0,2,"<func> or find-row <col> <val>",
       ((string)"If one argument is given, then find-row finds the row ")+
       "which maximizes the value of the "+
       "expression given in <func>, and then output the entire row. "+
       "Otherwise find-row finds the row for which the value in "+
       "column named <col> is as close as possible to the value <val>. "+
       "See command 'get-row' to get a row by it's index.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_row),
       both},
      {0,"fit","Fit two columns to a function (experimental).",0,7,
       "<x> <y> <yerr> <ynew> <par names> <func> <vals>","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_fit),
       both},
      {'f',"function","Set a column from a function.",0,2,
       "<func> <name>",
       ((string)"Set the column named <name> to the result of a function, ")+
       "<func>, in terms of the other columns. If the column does not "+
       "already exist, a new one is added to the table. For example, for "+
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
      {0,"get-unit","Get the units for a specified column.",0,1,"<column>",
       "Obtains the units for the specified column.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_unit),
       both},
      {0,"entry","Get a single entry in a table.",0,3,
       "<column> <row index> [value]",((std::string)"This command ")+
       "gets or sets the value in the specified column and row.",
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
       0,0,"","List the constants, column names and other info.",
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
      {0,"select-rows","Select rows for a new table.",
       0,1,"<row_spec>",((std::string)"Select the rows from a table for ")+
       "which the row specification in <row_spec> evaluates to a number "+
       "greater than 0.5",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_select_rows),both},
      {0,"select-rows2",
       "Select rows, with explicit column specification.",
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
       "<column> <unit>","Set the units for a specified column.",
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
      {0,"wstats","Show weighted column statistics.",0,2,"<col> <weights>",
       "Output the average, std. dev, max and min of <col>. ",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_wstats),
       both},
      {0,"sum","Add data from a second table object to current table.",
       0,2,"<file> [name]",((string)"Add all columns ")+
       "from the second table to their corresponding columns in the "+
       "current table, creating new columns if necessary.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
       both},
      {0,"nlines","Add 'nlines' as a constant to a table object.",0,0,
       "",((std::string)"Add a constant called 'nlines' to the table and ")+
       "set it equal to the number of lines (rows) in the table",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_nlines),
       both},
      {0,"to-hist","Convert a table to a histogram.",0,3,
       "<col> <n_bins> [wgts]",
       ((std::string)"The 'to-hist' command creates ")+
       "a 1D histogram from 'col' using exactly 'n_bins' bins and "+
       "(optionally) weighting the entries by the values in column 'wgts'. "+
       "The second form creates a 2D histogram from 'col1' and 'col2' "+
       "using N1 bins in the x direction and N2 bins in the y direction, "+
       "optionally weighting the entries by the column 'wgts'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_hist),
       both},
      {0,"to-hist-2d","Convert a table to a 2d histogram.",0,5,
       "<col x> <col y> <n_x_bins> <n_y_bins> [wgts]",
       ((std::string)"The 'to-hist-3d' creates a 2D histogram ")+
       "from 'col x' and 'col y' using 'n_x_bins' bins in the x "+
       "direction and 'n_y_bins' bins in the y direction, "+
       "optionally weighting the entries by the column 'wgts'.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_to_hist_2d),
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

  } else if (new_type=="table3d") {
    
    static const size_t narr=17;
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
       ((string)"Set the slice named <name> to the result of a function, ")+
       "<func>, in terms of the other slices. If the slice does not "+
       "already exist, a new one is created. For example, for "+
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
       0,0,"","List the slice names and print out grid info.",
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
       both},
      {0,"x-name","Get or set the 'x' grid name",
       0,1,"[name]","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_x_name),
       both},
      {0,"y-name","Get or set the 'y' grid name",
       0,1,"[name]","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_y_name),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor") {
    
    static const size_t narr=8;
    comm_option_s options_arr[narr]={
      {'l',"list","List the tensor rank and index sizes.",
       0,0,"","List the tensor rank and index sizes.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"rearrange","Rearrange the tensor.",
       -1,-1,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rearrange),
       both},
      {0,"to-table3d","Select two indices and convert to a table3d object.",
       -1,-1,"<x index> <y index> <slice name> [fixed 1] [fixed 2] ...",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d),both},
      {0,"to-table3d-sum",
       "Select two indices and convert to a table3d object.",
       -1,-1,"<x name> <y name> <slice name> [fixed 1] [fixed 2] ...",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d_sum),both},
      {0,"diag","Get diagonal elements.",-1,-1,"",
       ((string)"Extract only the elements on the main diagonal ")+
       "to create a double[] object.",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_diag),both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both},
      {0,"to-tensor-grid","Convert the tensor to a tensor_grid object.",
       -1,-1,"[function 1] [function 2] ...",
       ((string)"Convert a tensor to a tensor_grid object, using ")+
       "functions to specify the grid for each index. The functions "+
       "should be specified as functions of the variable 'i', which "+
       "runs from 0 to size-1 for each index. Any user-specified "+
       "functions are used up to the rank of the tensor, and if "+
       "not enough functions are specified, then the function 'i' is "+
       "used.",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_tensor_grid),both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="prob_dens_mdim_amr") {
    
    static const size_t narr=1;
    comm_option_s options_arr[narr]={
      {0,"to-table3d","Select two indices and convert to a table3d object.",
       -1,-1,((std::string)"<x index> <y index> ")+
       "<x name> <x points> <y name> <y points> <slice name>",
       "Select two indices and convert to a table3d object.",
       new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d),both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor<int>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]={
      {0,"rearrange","Rearrange the tensor.",
       -1,-1,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rearrange),
       both},
      {'l',"list","List the rank and sizes.",
       0,0,"","List the rank and sizes.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"to-table3d","Select two indices and convert to a table3d object.",
       -1,-1,"<x name> <y name> <slice name>",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d),both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor<size_t>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]={
      {0,"rearrange","Rearrange the tensor.",
       -1,-1,"","",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rearrange),
       both},
      {'l',"list","List the rank and sizes.",
       0,0,"","List the rank and sizes.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"to-table3d","Select two indices and convert to a table3d object.",
       -1,-1,"<x name> <y name> <slice name>",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d),both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor_grid") {
    
    static const size_t narr=7;
    comm_option_s options_arr[narr]={
      {'l',"list","List the slice names and print out grid info.",
       0,0,"","List the slice names and print out grid info.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
       both},
      {0,"to-table3d","Select two indices and convert to a table3d object.",
       -1,-1,"<x index> <y index> <new slice name> [values of fixed indices]",
       ((string)"To set the x- or y- grid names afterwards, ")+
       "use 'x-name' or 'y-name'.",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table3d),both},
      {0,"to-table","Convert to a two-column table object.",
       -1,-1,"<index> <grid name> <data name> [values of fixed indices]",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_to_table),both},
      {0,"slice","Slice to a smaller rank tensor_grid object.",
       -1,-1,"<index 1> <value 1> <index 2> <value 2> ...",
       "",new comm_option_mfptr<acol_manager>
       (this,&acol_manager::comm_slice),both},
      {0,"set-grid","Set the tensor grid.",-1,-1,
       ((std::string)"<function for rank 0> ")+
       "<function for rank 1> ... <function for rank n-1>",
       ((std::string)"Given a function which specifies the grid ")+
       "value as a function of the variable 'i' for each rank, "+
       "this command sets the tensor grid. The value of 'i' ranges "+
       "from 0 to m-1, where 'm' is the tensor size for each rank.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_grid),
       both},
      {0,"max","Find the maximum value and index.",0,0,"",
       "Compute the maximum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum value of and index.",0,0,"",
       "Compute the minimum value.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="hist") {

  } else if (new_type=="double[]") {

    static const size_t narr=9;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"sum","Compute the vector sum.",0,0,"",
       ((string)"Compute the vector sum."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
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
       both},      
      {0,"function","Set the values of the array given a function",0,1,
       "<function.",((string)"Set the values of the array ")+
       "given a user-specified function of 'i'. For example, "+
       "\"(sin(i)>1)*4\".",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
       both}      
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="int[]") {

    static const size_t narr=9;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"sum","Compute the vector sum.",0,0,"",
       ((string)"Compute the vector sum."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
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
       both},
      {0,"function","Set the values of the array given a function",0,1,
       "<function.",((string)"Set the values of the array ")+
       "given a user-specified function of 'i'. For example, "+
       "\"(sin(i)>1)*4\".",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
       both}     
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="size_t[]") {

    static const size_t narr=9;
    comm_option_s options_arr[narr]={
      {0,"sort","Sort the vector.",0,0,"",
       ((string)"Sorts the vector."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
       both},
      {0,"sum","Compute the vector sum.",0,0,"",
       ((string)"Compute the vector sum."),
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
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
       both},
      {0,"function","Set the values of the array given a function",0,1,
       "<function.",((string)"Set the values of the array ")+
       "given a user-specified function of 'i'. For example, "+
       "\"(sin(i)>1)*4\".",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vector<contour_line>") {

  } else if (new_type=="hist_2d") {

    static const size_t narr=2;
    comm_option_s options_arr[narr]={
      {0,"max","Find the maximum weight.",0,0,"",
       "Find the maximum weight and print out the location.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
       both},
      {0,"min","Find the minimum weight.",0,0,"",
       "Find the minimum weight and print out the location.",
       new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
       both}
    };
    cl->set_comm_option_vec(narr,options_arr);

  }
  
  return;
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
    cl->remove_comm_option("autocorr");
    cl->remove_comm_option("delete-col");
    cl->remove_comm_option("delete-rows");
    cl->remove_comm_option("delete-rows-tol");
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
    cl->remove_comm_option("nlines");
    cl->remove_comm_option("rename");
    cl->remove_comm_option("select");
    cl->remove_comm_option("select-rows");
    cl->remove_comm_option("select-rows2");
    cl->remove_comm_option("set-data");
    cl->remove_comm_option("set-unit");
    cl->remove_comm_option("sort");
    cl->remove_comm_option("stats");
    cl->remove_comm_option("sum");
    cl->remove_comm_option("to-hist");
    cl->remove_comm_option("to-hist-2d");
    cl->remove_comm_option("wstats");

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
    cl->remove_comm_option("x-name");
    cl->remove_comm_option("y-name");

  } else if (type=="tensor<int>") {
    
    cl->remove_comm_option("list");
    cl->remove_comm_option("to-table3d");
    cl->remove_comm_option("min");
    cl->remove_comm_option("max");
    cl->remove_comm_option("rearrange");

  } else if (type=="tensor<size_t>") {
    
    cl->remove_comm_option("list");
    cl->remove_comm_option("to-table3d");
    cl->remove_comm_option("min");
    cl->remove_comm_option("max");
    cl->remove_comm_option("rearrange");

  } else if (type=="tensor") {
    
    cl->remove_comm_option("list");
    cl->remove_comm_option("diag");
    cl->remove_comm_option("to-table3d");
    cl->remove_comm_option("to-table3d-sum");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    cl->remove_comm_option("to-tensor-grid");
    cl->remove_comm_option("rearrange");

  } else if (type=="prob_dens_mdim_amr") {
    
    cl->remove_comm_option("to-table3d");

  } else if (type=="tensor_grid") {
    
    cl->remove_comm_option("list");
    cl->remove_comm_option("to-table3d");
    cl->remove_comm_option("slice");
    cl->remove_comm_option("to-table");
    cl->remove_comm_option("set-grid");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");

  } else if (type=="hist_2d") {

    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    
    /*
      cl->remove_comm_option("deriv-x");
      cl->remove_comm_option("deriv-y");
      cl->remove_comm_option("interp");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("sum");
    */
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
    */
    
  } else if (type=="vector<contour_line>") {

  } else if (type=="double[]" || type=="int[]" || type=="size_t[]") {

    cl->remove_comm_option("deriv");
    cl->remove_comm_option("interp");
    cl->remove_comm_option("max");
    cl->remove_comm_option("min");
    cl->remove_comm_option("sort");
    cl->remove_comm_option("autocorr");
    cl->remove_comm_option("to-table");
    cl->remove_comm_option("function");
    cl->remove_comm_option("sum");
    
    /*
      cl->remove_comm_option("integ");
      cl->remove_comm_option("set-data");
      cl->remove_comm_option("stats");
      cl->remove_comm_option("sum");
      cl->remove_comm_option("to-hist");
    */

    /*
      } else if (type=="int[]" || type=="size_t[]") {
      
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
  } else if (type=="tensor") {
    tensor_obj.clear();
  } else if (type=="tensor<int>") {
    tensor_int_obj.clear();
  } else if (type=="tensor<size_t>") {
    tensor_size_t_obj.clear();
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
    {'c',"create","Create an object.",0,-1,"<type> [...]",
     ((string)"Create a new object of type <type>. For types char, ")+
     "double, int, size_t, and string, this takes one additional "+
     "argument which holds the value. For type table, "+
     "this option creates a new table with one column whose entries "+
     "are an evenly-spaced grid. In this case four additional arguments "+
     "are needed: the name of "+
     "the column, the first value, the maximum possible value, and the "+
     "increment between successive values. For array types double[] "+
     "int[], and size_t[], the user must specify the size of the array "+
     "and a function of the array index 'i' to fill the array. "+
     "If an object is currently in memory, it is deallocated before "+
     "creating the new object.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
     both},
    {0,"download","Download file from specified URL.",0,4,
     "<file> <URL> [hash, \"file:\"hash_filename, or \"none\"] [directory]",
     ((string)"Check if a file matches a specified hash, and if not, ")+
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
     "<file>. The allowed text file formats depend on the particular "+
     "type specified. For int, char, double, or size_t objects, "+
     "the file is assumed to begin with the descired object and it is "+
     "read using operator>>(). For string objects, the first line is "
     "read using std::getline(). For array objects, it is assumed "+
     "that all array entries are on the first line of the file and no "+
     "carriage returns are present between entries. For table objects, "
     "the first line of the file must either contain numeric "+
     "data or column names "+
     "separated by white space, without carriage returns, except for "+
     "the one at the end of the line. If the first line contains "+
     "column names, the second line may optionally contain unit "+
     "expressions for each column, enclosed by square brackets. "+
     "All remaining lines are assumed "+
     "to contain data with the same number of columns as the first line. "+
     "For table3d objects, the data must be stored in columns "+
     "with the first column specifying the x-axis grid point and "+
     "the second column specifying the y-axis grid point. The "+
     "remaining columns give the data for each slice at that point. "+
     "Each grid point must correspond to a row in the file, but "+
     "the lines need not be in any particular order. The columns may "+
     "have one header line at top which specifies the names of the x- "+
     "and y-grids and the names of each slice (in order).",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
     both},
    {0,"get-conv","Get a unit conversion factor.",0,2,
     "<old unit> <new unit>",((string)"This command gets a unit ")+
     "conversion factor. It only works if the conversion is one of the "
     "hard-coded O2scl conversions or if HAVE_POPEN is defined and "+
     "the 'units' command is available in the current "+
     "path. For example, 'get-conv MeV erg' returns 1.602177e-6 and 1 MeV "+
     "is equivalent to 1.602177e-6 erg. The conversion factor is output "+
     "at the current precision, but is always internally stored with "+
     "full double precision. O2scl has several unit conversions which "+
     "implicitly assume hbar=c=1.",
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
     "char, string, int[], double[], string[], size_t, size_t[], "+
     "uniform_grid<double>, tensor_grid, tensor, tensor<int>, "+
     "tensor<size_t>, or prob_dens_mdim_amr.",
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

  cl->cmd_name="acol";
  
  cl->desc=((string)"acol: A data viewing and ")+
    "processing program for O2scl.\n";
  
  string dsc="\nNotes:\n\n";
  dsc+="1. Help for general commands may be obtained with 'help ";
  dsc+="<command>'. Help for \n   type-specific commands can be obtained ";
  dsc+="by 'help <type> <command>'. A \n   list of commands for each type ";
  dsc+="can be obtained with 'commands <type>'.\n   Required arguments ";
  dsc+="are surrounded by ";
  dsc+="<>'s and optional arguments are\n   surrounded by []'s.\n";
  dsc+="2. Options may also be specified in the environment variable ";
  dsc+="ACOL_DEFAULTS.\n";
  dsc+="3. Long options may be preceeded by two dashes.\n";
  dsc+="4. In order to avoid confusion between arguments and functions,\n";
  dsc+="   use parenthesis and quotes, i.e. \"(-x*2)\" instead of -x*2.\n";
  dsc+="5. Also, do not use a unary minus next to a binary operator,";
  dsc+=" i.e. use\n   \"a>(-1)\" instead of \"a>-1\".\n\n";
  dsc+="Known operators:\n\n() ^ * / % + - == != < > && || << >> >= <=\n\n";
  dsc+="Known functions:\n\n";
  dsc+="exp(x) log(x) log10(x) sin(x) cos(x) tan(x) sqrt(x) abs(x)\n";
  dsc+="asin(x) acos(x) atan(x) sinh(x) cosh(x) tanh(x)\n";
  dsc+="asinh(x) acosh(x) atanh(x)\n\n";
  /*
    dsc+="atan2(x,y) if(x,y,z)\n";
    dsc+="cot(x) csc(x) sec(x)\n";
    dsc+="ceil(x) floor(x) int(x) max(x,y) min(x,y)\n";
  */
  
  dsc+="List of additional type-specific commands\n";
  dsc+="  (use 'help <type> <command>' for more info):\n\n";
  std::map<std::string,std::vector<std::string> >::iterator it;
  for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
    std::vector<std::string> &clist=it->second;
    string stempx=it->first+": ";
    for(size_t j=0;j<clist.size()-1;j++) {
      stempx+=clist[j]+", ";
    }
    stempx+=clist[clist.size()-1]+"\n";
    std::vector<std::string> stempy;
    rewrap(stempx,stempy,77);
    for(size_t j=0;j<stempy.size();j++) {
      if (j>0) {
	dsc+=((std::string)"  ")+stempy[j]+"\n";
      } else {
	dsc+=stempy[j]+"\n";
      }
    }
  }
  dsc+="\n";
  
#ifndef O2SCL_UBUNTU_PKG
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
  p_def_args.help=((std::string)"The default arguments from the ")+
    "environment varable ACOL_DEFAULTS.";
  p_prec.help="The numerical precision.";
  p_verbose.help="Control the amount of output.";
  p_compress.help=((std::string)"If true, enable compression ")+
    "(defaults to true if compression was enabled in O2scl).";
  p_ncols.help="The number of columns for screen output.";
  p_interp_type.help=((std::string)"The interpolation type ")+
    "(1=linear, 2=cubic spline, 3=periodic cubic spline, 4=Akima, "+
    "5=periodic Akima, 6=monotonic, 7=Steffen's monotonic).";
  p_names_out.help="If true, output column names at top.";
  p_pretty.help="If true, make the output more readable.";
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

  return 0;
}

int acol_manager::run(int argc, char *argv[], bool full_process) {

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
    {0,"commands","List available commands.",0,1,"[type]","",
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
  
  if (full_process) {
    
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

