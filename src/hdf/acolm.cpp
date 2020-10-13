/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
#include <o2scl/cursesw.h>

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

  type_list.push_back("table");
  type_list.push_back("table3d");
  type_list.push_back("hist");
  type_list.push_back("hist_2d");
  type_list.push_back("vector<contour_line>");
  type_list.push_back("int");
  type_list.push_back("double");
  type_list.push_back("char");
  type_list.push_back("string");
  type_list.push_back("int[]");
  type_list.push_back("double[]");
  type_list.push_back("string[]");
  type_list.push_back("size_t");
  type_list.push_back("size_t[]");
  type_list.push_back("uniform_grid<double>");
  type_list.push_back("tensor_grid");
  type_list.push_back("tensor");
  type_list.push_back("tensor<int>");
  type_list.push_back("tensor<size_t>");
  type_list.push_back("prob_dens_mdim_amr");
  vector_sort<vector<string>,string>(type_list.size(),type_list);
  
  {
    vector<std::string> itmp={"value"};
    type_comm_list.insert(std::make_pair("int",itmp));
    type_comm_list.insert(std::make_pair("double",itmp));
    type_comm_list.insert(std::make_pair("char",itmp));
    type_comm_list.insert(std::make_pair("size_t",itmp));
    type_comm_list.insert(std::make_pair("string",itmp));
  }
  {
    vector<std::string> itmp={"ac-len","average-rows",
			      "assign","cat","convert-unit",
			      "correl","delete-col",
			      "delete-rows",
			      "delete-rows-tol","deriv","deriv2",
			      "entry-grid",
			      "find-row","fit","function",
			      "get-row","get-unit","entry","index",
			      "insert","insert-full","integ","interp",
			      "list","max","min","nlines","refine","rename",
			      "select","select-rows","select-rows2",
			      "set-data","set-unit","sort","stats","sum",
			      "to-hist","to-hist-2d","to-table3d","wstats",
    };
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("table",itmp));
  }
  {
    vector<std::string> itmp={"cat","contours","deriv-x","deriv-y",
			      "function","entry","entry-grid","get-grid",
			      "insert","interp","stats","select",
			      "list","max","min","rename","set-data",
			      "slice","slice-hist","sum",
			      "to-tensor-grid","x-name","y-name"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("table3d",itmp));
  }
  {
    vector<std::string> itmp={"list","min","max","to-table3d",
			      "rearrange"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor<int>",itmp));
    type_comm_list.insert(std::make_pair("tensor<size_t>",itmp));
  }
  {
    vector<std::string> itmp={"list","diag","to-table3d","to-table3d-sum",
			      "max","min","to-tensor-grid","rearrange",
			      "entry","function","sum","stats"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor",itmp));
  }
  {
    vector<std::string> itmp={"to-table3d"};
    type_comm_list.insert(std::make_pair("prob_dens_mdim_amr",itmp));
  }
  {
    vector<std::string> itmp={"list","to-table3d","slice","to-table",
			      "set-grid","max","min","rearrange",
			      "get-grid","interp","entry","to-tensor",
			      "entry-grid","function","sum","stats",
			      "binary"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("tensor_grid",itmp));
  }
  {
    vector<std::string> itmp={"max","min","contours","to-table3d"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist_2d",itmp));
  }
  {
    vector<std::string> itmp={"to-table","function"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("hist",itmp));
  }
  {
    vector<std::string> itmp={"deriv","interp","max","min","sort",
			      "to-table","function","sum"};
    vector_sort<vector<string>,string>(itmp.size(),itmp);
    type_comm_list.insert(std::make_pair("double[]",itmp));
    type_comm_list.insert(std::make_pair("int[]",itmp));
    type_comm_list.insert(std::make_pair("size_t[]",itmp));
  }

}

void acol_manager::command_add(std::string new_type) {

  const int both=cli::comm_option_both;
  
  if (new_type=="int") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"value","Get or set the value of the int object.",
	0,1,"[value]","Get or set the value of the int object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="double") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"value","Get or set the value of the double object",
	0,1,"[value-spec]","Get or set the value of the double object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="char") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"value","Get or set the value of the char object",
	0,1,"[value]","Get or set the value of the char object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="size_t") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"value","Get or set the value of the size_t object.",
	0,1,"[value]","Get or set the value of the size_t object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="string") {
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"value","Get or set the value of the string object.",
	0,1,"[value]","Get or set the value of the string object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_value),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
  } else if (new_type=="table") {
    static const size_t narr=41;
    comm_option_s options_arr[narr]=
      {{0,"ac-len","Autocorrelation length using 'acor'.",0,1,"<colum>","",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_ac_len),
	both},
       {0,"average-rows","Average rows of some or all columns together",
	0,3,"<column or '*' for all> <window> [block averages]",
	((string)"The first argument is the column to be ")+
	"modified. If the first argument is '*', then all columns are "+
	"averaged. The second argument is the size "+
	"of the window. If the third argument evaluates to false, then "+
	"block averages instead of rolling averages are computed, and then "+
	"the number of rows is divided by the window parameter. If block "+
	"averages are requested, then the first argument must be '*'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_average_rows),both},
       {'a',"assign","Assign a constant, e.g. assign pi acos(-1) .",
	0,2,"<name> [val]",
	((string)"Assign a constant value to a name for the present table. ")+
	"Valid constant values are things like 1.618 or "+
	"acos(-1.0) or sin(4^5). "+
	"To remove an assignment, call assign with a blank value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_assign),
	both},
       {0,"cat",
	"Concatenate a second table object onto current table.",0,2,
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
       {0,"correl","Compute correlation between two columns",0,2,
	"<column 1> <column 2>",((string)"Compute the correlation ")+
	"coefficient between two columns, or, if no arguments are "+
	"given, then compute the correlation coefficients between all "+
	"pairs of columns.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_correl),
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
       {0,"deriv",
	"Derivative of a function defined by two columns.",
	0,3,"<x> <y> <name>",
	((string)"Create a new column named <name> filled with the ")+
	"derivative of the function y(x) obtained from columns <x> and <y>. ",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
	both},
       {0,"deriv2","Second derivative.",0,3,"<name> <x> <y>",
	((string)"Create a new column named <name> filled with the second ")+
	"derivative of the function y(x) obtained from columns <x> and <y>. ",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_deriv2),
	both},
       {0,"entry","Get or set a single entry in a table.",0,3,
	"<column> <row> [value or \"none\"]",
	((std::string)"This command ")+
	"gets or sets the value in the specified column and row. If "+
	"\"none\" is specified as the third argument, then \"entry\" "+
	"just prints out the specified entry as if the third argument "+
	"was not specified.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry),both},
       {0,"entry-grid","Get or set a single entry in a table.",0,4,
	"<index column> <index value> <target column> [value or \"none\"]",
	((std::string)"The \"entry-grid\" command ")+
	"first looks for the value closest to <index value> in the column "+
	"<index column> to determine a row in the table. "+
	"Next \"entry-grid\" gets or sets the value of the "+
	"target column in that row. If "+
	"\"none\" is specified as the fourth argument, then \"entry\" "+
	"just prints out the specified entry as if the third argument "+
	"was not specified.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry_grid),both},
       {0,"find-row","Find a row which maximizes a function.",
	0,2,"<func> or find-row <col> <val>",
	((string)"If one argument is given, then find-row finds the row ")+
	"which maximizes the value of the "+
	"expression given in <func>, and then output the entire row. "+
	"Otherwise find-row finds the row for which the value in "+
	"column named <col> is as close as possible to the value <val>. "+
	"See command 'get-row' to get a row by it's index.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_find_row),both},
       {0,"fit","Fit two columns to a function (experimental).",0,7,
	"<x> <y> <yerr> <ynew> <par names> <func> <vals>","",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_fit),both},
       {0,"function","Set a column from a function.",0,2,
	"<func> <name>",
	((string)"Set the column named <name> to the result of a function, ")+
	"<func>, in terms of the other columns. If the column does not "+
	"already exist, a new one is added to the table. For example, for "+
	"a table containing columns named 'c1' and 'c2', 'function "+
	"c1-c2 c3' would create a new column c3 which contains the "+
	"difference of columns 'c1' and 'c2'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_function),both},
       {0,"get-row","Get a row by index.",
	0,1,"<index>",((string)"Get a row by index. The first row ")+
	"has index 0, and the last row has index n-1, where n "+
	"is the total number of rows as returned by the 'list' command. "+
	"The 'index' command creates a column of row indexes. "+
	"To find a row which contains a particular value or maximizes "+
	"a specified function, use 'find-row'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_get_row),both},
       {0,"get-unit","Get the units for a specified column.",0,1,"<column>",
	"Obtains the units for the specified column.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_get_unit),both},
       {'N',"index","Add a column containing the row numbers.",0,1,
	"[column name]",
	((string)"Define a new column named [column name] and fill ")+
	"the column with the row indexes, beginning with zero. If "+
	"no argument is given, the new column is named 'N'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_index),both},
       {0,"insert","Interpolate a column from another file.",0,6,
	"<file> <table> <oldx> <oldy> <newx> [newy]",
	((string)"Insert a column from file <fname> interpolating it ")+
	"into the current table. The column <oldy> is the "+
	"columns in the file which is to be inserted into the table, "+
	"using the column <oldx> in the file and <newx> in the table. "+
	"The new column in the table is named <oldy>, or it is named "+
	"[newy] if the additional argument is given. ",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_insert),both},
       {0,"insert-full",
	"Insert a table from another file.",0,4,
	"<fname> [table name] [old_x new_x]",
	((string)"Insert all columns from file <fname> ")+
	"into the current table. The first table is used or the "+
	"table object named table_name, if specified. "+
	"If index columns old_x and new_x "+
	"are not specified, then the insert requires both the current "+
	"and the source table to have the same number of rows. If they "+
	"are specified, then interpolation using those index columns "+
	"is used. If columns in the new table are not present in the "+
	"current table, then they are added automatically. If a column "+
	"in the current table has the same name as one in the new table "+
	"then it is rewritten with new data, with one exception. If "+
	"a column in the new table has the same name as old_x, then it "+
	"is left unmodified.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_insert_full),both},
       {0,"integ",
	"Integrate a function specified by two columns.",
	0,3,"<x> <y> <name>",
	((string)"Create a new column named <name> filled with the ")+
	"integral of the function y(x) obtained from columns <x> and <y>. ",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_integ),both},
       {0,"interp","Interpolate a number into a column.",0,3,
	"<x name> <x value> <y name>",
	((string)"Interpolate <x value> from column ")+
	"named <x name> into column named <y name>.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_interp),both},
       {'l',"list","List the constants, column names and other info.",
	0,0,"","List the constants, column names and other info.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_list),both},
       {0,"max","Find the maximum value of a column.",0,1,"<col>",
	"Compute the maximum value of column <col>.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_max),both},
       {0,"min","Find the minimum value of a column.",0,1,"<col>",
	"Compute the minimum value of column <col>.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_min),both},
       {0,"nlines","Add 'nlines' as a constant to a table object.",0,0,
	"",((std::string)"Add a constant called 'nlines' to the table and ")+
	"set it equal to the number of lines (rows) in the table",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_nlines),both},
       {0,"refine","Refine the table.",0,2,"<index column> <factor>","",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_refine),both},
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
	"and use the rules from the fnmatch() function. A '*' matches "+
	"everything, and '?' matches any single character. The pattern "+
	"[list] matches any character in list, and the pattern "+
	"[!list] matches all characters not in the list. The pattern "+
	"[*] matches the character '*' and [?] matches the character '?'. "+
	"All of the rows of data "+
	"are copied over. If functions are specified, the result can be "+
	"named using '='.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_select),
	both},
       {0,"select-rows","Select rows for a new table.",
	0,1,"<row_spec>",((std::string)"Select the rows from a table for ")+
	"which the row specification in <row_spec> evaluates to a number "+
	"greater than 0.5",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_select_rows2),both},
       {0,"select-rows2",
	"Select rows, with column specification (deprecated.",
	0,-1,"<row_spec>",
	((std::string)"Select the rows from a table for ")+
	"which the row specification in <row_spec> evaluates to a number "+
	"greater than 0.5 . This command should be equivalent to but "+
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
	((std::string)"The 'to-hist-2d' command creates a 2D histogram ")+
	"from 'col x' and 'col y' using 'n_x_bins' bins in the x "+
	"direction and 'n_y_bins' bins in the y direction, "+
	"optionally weighting the entries by the column 'wgts'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_hist_2d),both},
       {0,"to-table3d","Convert a table to a table3d object.",0,4,
	"<x column> <y column> [empty value] [eps]",
	((std::string)"The 'to-table3d' creates a table3d object using ")+
	"'x column' and 'y column' as the data for the x and y grids. "+
	"If 'empty value', then this value is used for points not given "+
	"by the table. If 'eps' is specified, then use that value as the "+
	"minimum value between grid points.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both}
      };
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="table3d") {
    
    static const size_t narr=23;
    comm_option_s options_arr[narr]=
      {{0,"to-tensor-grid",
	"Convert a slice of the table3d to a tensor_grid object.",
	1,1,"<slice>","",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_tensor_grid),both},
       {0,"get-grid","Get the table3d grid.",0,0,"",
	"Output the table3d grid as a series of columns.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_grid),
	both},
       {'s',"select","Select slices for a new table3d object.",-1,-1,
	"<slice spec.>",
	((string)"Select creates a new table3d from the present table3d, ")+
	"including only the slices specified in <slice spec.>. The slice "+
	"specification is a list of slice names, functions, or patterns "+
	"which match "+
	"the slice names. Patterns must be preceeded by a colon ':' "+
	"and can use wildcards like '*' and '?'. All of the rows of data "+
	"are copied over. If functions are specified, the result can be "+
	"named using '='. ",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_select),
	both},
       {0,"cat",
	"Concatenate data from a second table3d onto current table3d.",0,2,
	"<file> [name]",((string)"Add all slices from the ")+
	"second table3d object which aren't already present in the "+
	"current table3d object.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_cat),
	both},
       {0,"contours","Create contour lines from a table3d slice.",
	0,5,"[\"frac\"] <value> <slice_name> [output_filename object_name]",
	((string)"If the argument \"frac\" is not present, the ")+
	"\"contours\" command constructs a set of contour lines using "+
	"the data in slice named <slice> at the fixed value given in "+
	"<value>. If two additional arguments are given, then the "+
	"contour lines are stored in the file named output_filename "+
	"and the object is named object_name. If the file does not "+
	"exist, it is created. If no contours are found, then no file "+
	"I/O is performed and the current table3d object is unmodified. "+
	"If the argument \"frac\" is present, then the operation is "+
	"the same except that <value> is interpreted as a fraction of "+
	"the total integral under the data.",
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
       {0,"entry","Get a single entry in a table3d object.",0,4,
	"<slice> <x index> <y index> [value or \"none\"]","",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_entry),
	both},
       {0,"entry-grid","Get a single entry in a table3d object.",0,4,
	"<slice> <x value> <y value> [value or \"none\"]","",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry_grid),both},
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
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_set_data),both},
       {0,"slice","Construct a slice.",2,2,
	"<\"x\" or \"y\"> <value>",
	((string)"Extract a slice of a table3d object at fixed x or fixed y ")+
	"to create a new table object. This function uses interpolation "+
	"with the current interpolation type to interpolate all of the "+
	"slices in the table3d object to create a table with a column "+
	"for each slice.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_slice),
	both},
       {0,"slice-hist","Construct a histogram from a slice.",1,1,
	"<slice>","",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_slice_hist),both},
       {0,"stats","Show slice statistics.",0,1,"<slice>",
	"Output the size, sum, max and min of <slice>. ",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_stats),
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
    
    static const size_t narr=12;
    comm_option_s options_arr[narr]=
      {
       {0,"sum","Output the sum of all the tensor entries.",0,0,"",
	((string)"The \"sum\" command outputs the total tensor size ")+
	"and the sum over all entries.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
	both},
       {0,"stats","Show stats for the data in the tensor.",0,0,"",
	((string)"The 'stats' command outputs the number of entries, ")+
	"their mean, standard deviation, minimum and maximum. It also "+
	"counts the number of infinite or NaN values.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_stats),
	both},
       {0,"diag","Get diagonal elements.",-1,-1,"",
	((string)"Extract only the elements on the main diagonal ")+
	"to create a double[] object.",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_diag),both},
       {0,"entry","Get a single entry in a tensor object.",
	-1,-1,"<value 1> <value 2> <value 3> ... [value or \"none\"]",
	"",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry),both},
       {'f',"function","Set tensor value from a function.",0,1,
	"[cond. function] <function of v, i0, i1, ...>",
	((string)"The \"function\" command ")+
	"sets all entries in a tensor equal to a user-specified "+
	"mathematical function of the indices. When the conditional "+
	"function evaluates to a number "+
	"less than or equal to 0.5, then the tensor entry will be "
	"unchanged. (For more help with "+
	"functions, type \""+cl->cmd_name+" -help functions\".)",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
	both},
       {'l',"list","List the tensor rank and index sizes.",
	0,0,"","List the tensor rank and index sizes.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
	both},
       {0,"max","Find the maximum value and index.",0,0,"",
	"Compute the maximum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
	both},
       {0,"min","Find the minimum value of and index.",0,0,"",
	"Compute the minimum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
	both},
       {0,"rearrange","Rearrange the tensor.",
	-1,-1,"<index spec. 1> [index spec. 2] ...",
	((std::string)"Index specifications are: index(ix), fixed(ix), ")+
	"sum(ix), trace(ix1,ix2), reverse(ix), and range(ix,start,end).",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_rearrange),both},
       {0,"to-table3d","Select two indices and convert to a table3d object.",
	-1,-1,"<x index> <y index> <slice name> [fixed 1] [fixed 2] ...",
	((string)"This command uses two indices in the current ")+
	"tensor object to create a table3d object. The values for "+
	"the remaining indices fixed to [fixed 1], "+
	"[fixed 2], etc. in that order. For example, \"to-table3d 3 1 "+
	"z 5 3\" uses index 3 for the "+
	"x coordinate of the new table3d object, uses index 1 for "+
	"the y coordinate of the new table3d object, uses 5 for index "+
	"0, and uses 3 for index 2."+
	"The x- and y-grids in "+
	"the table3d object are named \"x\" and \"y\" and filled with "+
	"the grid index by default."+
	"To set the x- or y-grid names afterwards, "+
	"use commands 'x-name' and 'y-name'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both},
       {0,"to-table3d-sum",
	"Select two indices and convert to a table3d object.",
	-1,-1,"<x name> <y name> <slice name> [fixed 1] [fixed 2] ...",
	"",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d_sum),both},
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
    
  } else if (new_type=="tensor<int>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]=
      {
       {'l',"list","List the rank and sizes.",
	0,0,"","List the rank and sizes.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
	both},
       {0,"max","Find the maximum value and index.",0,0,"",
	"Compute the maximum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
	both},
       {0,"min","Find the minimum value of and index.",0,0,"",
	"Compute the minimum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
	both},
       {0,"rearrange","Rearrange the tensor.",
	-1,-1,"<index spec. 1> [index spec. 2] ...",
	((std::string)"Index specifications are: index(ix), fixed(ix), ")+
	"sum(ix), trace(ix1,ix2), reverse(ix), and range(ix,start,end).",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_rearrange),both},
       {0,"to-table3d","Select two indices and convert to a table3d object.",
	-1,-1,"<x name> <y name> <slice name>",
	((string)"This command uses two indices in the current ")+
	"tensor<int> object to create a table3d object. The values for "+
	"the remaining indices fixed to [fixed 1], "+
	"[fixed 2], etc. in that order. For example, \"to-table3d 3 1 "+
	"z 5 3\" uses index 3 for the "+
	"x coordinate of the new table3d object, uses index 1 for "+
	"the y coordinate of the new table3d object, uses 5 for index "+
	"0, and uses 3 for index 2."+
	"The x- and y-grids in "+
	"the table3d object are named \"x\" and \"y\" and filled with "+
	"the grid index by default."+
	"To set the x- or y-grid names afterwards, "+
	"use commands 'x-name' and 'y-name'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both}
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="tensor<size_t>") {
    
    static const size_t narr=5;
    comm_option_s options_arr[narr]=
      {
       {0,"rearrange","Rearrange the tensor.",
	-1,-1,"<index spec. 1> [index spec. 2] ...",
	((std::string)"Index specifications are: index(ix), fixed(ix), ")+
	"sum(ix), trace(ix1,ix2), reverse(ix), and range(ix,start,end).",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_rearrange),both},
       {'l',"list","List the rank and sizes.",
	0,0,"","List the rank and sizes.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
	both},
       {0,"to-table3d","Select two indices and convert to a table3d object.",
	-1,-1,"<x name> <y name> <slice name>",
	((string)"This command uses two indices in the current ")+
	"tensor<size_t> object to create a table3d object. The values for "+
	"the remaining indices fixed to [fixed 1], "+
	"[fixed 2], etc. in that order. For example, \"to-table3d 3 1 "+
	"z 5 3\" uses index 3 for the "+
	"x coordinate of the new table3d object, uses index 1 for "+
	"the y coordinate of the new table3d object, uses 5 for index "+
	"0, and uses 3 for index 2."+
	"The x- and y-grids in "+
	"the table3d object are named \"x\" and \"y\" and filled with "+
	"the grid index by default."+
	"To set the x- or y-grid names afterwards, "+
	"use commands 'x-name' and 'y-name'.",
	new comm_option_mfptr<acol_manager>
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
    
    static const size_t narr=17;
    comm_option_s options_arr[narr]=
      {
       {0,"sum","Output the sum of all the tensor entries.",0,0,"",
	((string)"The \"sum\" command outputs the total tensor size ")+
	"and the sum over all entries.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sum),
	both},
       {0,"binary","Apply a binary function to two tensor_grid objects.",
	-1,-1,"<file> <object name> <function>",
	((string)"Long ")+"desc.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_binary),
	both},
       {0,"stats","Show stats for the data in the tensor.",0,0,"",
	((string)"The 'stats' command outputs the number of entries, ")+
	"their mean, standard deviation, minimum and maximum. It also "+
	"counts the number of infinite or NaN values.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_stats),
	both},
       {0,"entry","Get a single entry in a tensor_grid object.",
	-1,-1,"<index 1> <index 2> <index 3> ... [value or \"none\"]",
	((string)"The \"entry\" command gets or sets a value in the ")+
	"tensor_grid object. The arguments are a list of indices and "+
	"(optionally) a new value to store in that location.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry),both},
       {0,"entry-grid","Get a single entry in a tensor_grid object.",
	-1,-1,"<value 1> <value 2> <value 3> ... [value or \"none\"]",
	((string)"The \"entry-grid\" command gets or sets a value in the ")+
	"tensor_grid object. The arguments are a list of grid values and "+
	"(optionally) a new value to store in the location closest to "+
	"the specified grid values.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_entry_grid),both},
       {'f',"function","Set tensor value from a function.",0,1,
	"[conditional func.] <func. of v, i0, i1, ... and x0, x1, ...>",
	((string)"The \"function\" command sets ")+
	"all the data entries in a tensor_grid equal to a user-specified "+
	"mathematical function of the indices (i0, i1, ...) or grid "+
	"points (x0, x1, ...). If two function arguments are given and "+
	"if the first function argument is not \"none\", then "+
	"the first function specifies which tensor entries are to be "+
	"modified. When the conditional function evaluates to a number "+
	"less than or equal to 0.5, then the tensor entry will be "
	"unchanged. (For more help with "+
	"functions, type \""+cl->cmd_name+" -help functions.\".)",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_function),both},
       {0,"get-grid","Get the tensor grid.",0,0,"",
	"Output the tensor grid as a series of columns.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_get_grid),both},
       {0,"interp","Linearly interpolate in the grid.",
	-1,-1,"<value 1> <value 2> <value 3> ...",
	((string)"The command \"interp\" uses linear interpolation to ")+
	"interpolate an array with size equal to the tensor rank "+
	"into the tensor grid and outputs the result.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_interp),both},
       {'l',"list","List the slice names and print out grid info.",
	0,0,"","List the slice names and print out grid info.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
	both},
       {0,"max","Find the maximum value and indices.",0,0,"",
	"Compute the maximum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
	both},
       {0,"min","Find the minimum value and indices.",0,0,"",
	"Compute the minimum value.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
	both},
       {0,"rearrange","Rearrange the tensor_grid object.",
	-1,-1,"<index spec. 1> [index spec. 2] ...",
	((std::string)"Index specifications are: index(ix), fixed(ix), ")+
	"sum(ix), trace(ix1,ix2), reverse(ix), range(ix,start,end), "+
	"interp(ix,value), grid(ix,begin,end,n_bins,log), and "+
	"gridw(ix,begin,end,bin_width,log).",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_rearrange),both},
       {0,"slice","Slice to a smaller rank tensor_grid object.",
	-1,-1,"<index 1> <value 1> <index 2> <value 2> ...",
	"",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_slice),both},
       {0,"set-grid","Set the tensor grid.",0,2,
	"<index> <func. or vector spec> ",
	((std::string)"The first argument for the \"set-grid\" command ")+
	"specifies the index for which grid to set. The second argument "+
	"specifies the grid. If it contains a ':', it is assumed "+
	"to be a vector specification (see 'help vector-spec'). Otherwise, "+
	"the argument is assumed to be a function which specifies the grid "+
	"value as a function of the variables 'i' and 'x'. The value of "+
	"'i' ranges from 0 to m-1, where 'm' is the tensor size for each "+
	"rank and the value of 'x' is equal to the previous grid value.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_set_grid),both},
       {0,"to-table3d","Select two indices and convert to a table3d object.",
	-1,-1,"<x index> <y index> <new slice> [value 1] [value 2] ...",
	((string)"This command uses two indices in the current ")+
	"tensor_grid object to create a table3d object. The values for "+
	"the remaining indices are by interpolation to [value 1], "+
	"[value 2], etc. in that order. For example, \"to-table3d 3 1 "+
	"z 0.5 2.0\" uses index 3 for the "+
	"x coordinate of the new table3d object, uses index 1 for "+
	"the y coordinate of the new table3d object, uses interpolation "+
	"to set the value of the index 0 to 0.5, and uses interpolation "+
	"to set the value of index 2 to to 2.0. The x- and y-grids in "+
	"the table3d object are named \"x\" and \"y\" by default. "+
	"To set the x- or y-grid names afterwards, "+
	"use commands 'x-name' and 'y-name'.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both},
       {0,"to-table","Convert to a two-column table object.",
	-1,-1,"<index> <grid name> <data name> [values of fixed indices]",
	"",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table),both},
       {0,"to-tensor","Convert to a tensor object.",
	0,0,"Convert the tensor-grid object to a tensor object.",
	"",new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_tensor),both}
      };
    cl->set_comm_option_vec(narr,options_arr);

  } else if (new_type=="prob_dens_mdim_amr") {
    
    static const size_t narr=1;
    comm_option_s options_arr[narr]=
      {
       {0,"to-table3d","Select two indices and convert to a table3d object.",
	-1,-1,((std::string)"<x index> <y index> ")+
	"<x name> <x points> <y name> <y points> <slice name>",
	"Select two indices and convert to a table3d object.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both}
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="hist") {

    static const size_t narr=2;
    comm_option_s options_arr[narr]=
      {
       {0,"function","Apply a function to the weights.",0,1,"",
	((string)"Apply a function to the weights."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_function),both},
       {0,"to-table","Convert to a table object.",0,0,"",
	((string)"Convert to a table object."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table),both}
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="double[]") {
    
    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
       {0,"deriv",
	"Replace the array with its derivative.",0,0,"",
	((string)"Replace the array with its derivative using the ")+
	"current interpolation type.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_deriv),both},
       {0,"interp","Interpolate an index into the array.",0,1,
	"<x value>",
	((string)"Interpolate <x value> in the array."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_interp),both},
       {0,"max","Find the maximum value and index.",0,0,"",
	"Compute the maximum value of column <col>.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_max),both},
       {0,"min","Find the minimum value of and index.",0,0,"",
	"Compute the minimum value of column <col>.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_min),both},
       {0,"sort","Sort the vector.",0,0,"",
	((string)"Sorts the vector."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_sort),both},
       {0,"sum","Compute the vector sum.",0,0,"",
	((string)"Compute the vector sum."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_sum),both},
       {0,"to-table","Convert to a table given a column name",0,1,
	"<column name>",
	"Convert to a table given a column name.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table),both},      
       {0,"function","Set the values of the array given a function",0,1,
	"<function>.",((string)"Set the values of the array ")+
	"given a user-specified function of 'i'. For example, "+
	"\"(sin(i)>1)*4\".",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_function),both}      
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="int[]") {

    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
       {0,"sort","Sort the vector.",0,0,"",
	((string)"Sorts the vector."),
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_sort),both},
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
       {0,"deriv",
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
       {0,"to-table","Convert to a table given a column name",0,1,
	"<column name>",
	"Convert to a table given a column name.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table),both},
       {0,"function","Set the values of the array given a function",0,1,
	"<function>.",((string)"Set the values of the array ")+
	"given a user-specified function of 'i'. For example, "+
	"\"(sin(i)>1)*4\".",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_function),both}     
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="size_t[]") {

    static const size_t narr=8;
    comm_option_s options_arr[narr]=
      {
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
       {0,"deriv",
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
       {0,"to-table","Convert to a table given a column name",0,1,
	"<column name>",
	"Convert to a table given a column name.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table),both},
       {0,"function","Set the values of the array given a function",0,1,
	"<function>.",((string)"Set the values of the array ")+
	"given a user-specified function of 'i'. For example, "+
	"\"(sin(i)>1)*4\".",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
	both}
      };
    cl->set_comm_option_vec(narr,options_arr);
    
  } else if (new_type=="vector<contour_line>") {

  } else if (new_type=="hist_2d") {

    static const size_t narr=4;
    comm_option_s options_arr[narr]=
      {
       {0,"max","Find the maximum weight.",0,0,"",
	"Find the maximum weight and print out the location.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
	both},
       {0,"min","Find the minimum weight.",0,0,"",
	"Find the minimum weight and print out the location.",
	new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
	both},
       {0,"to-table3d","Convert to a table3d object",-1,-1,
	"<x name> <y name> <weight name>",
	"Convert to a table3d object using the specified names.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_to_table3d),both},
       {0,"contours","Create contour lines from a table3d slice.",
	0,4,"[\"frac\"] <value> [output file] [output name]",
	((string)"If the argument \"frac\" is not present, the ")+
	"\"contours\" command constructs a set of contour lines using "+
	"at the fixed value given in "+
	"<value>. If two additional arguments are given, then the "+
	"contour lines are stored in the file named output_filename "+
	"and the object is named object_name. If the file does not "+
	"exist, it is created. If no contours are found, then no file "+
	"I/O is performed and the current table3d object is unmodified."+
	"If the argument \"frac\" is present, then the operation is "+
	"the same except that <value> is interpreted as a fraction of "+
	"the total integral under the data.",
	new comm_option_mfptr<acol_manager>
	(this,&acol_manager::comm_contours),both},
      };
    cl->set_comm_option_vec(narr,options_arr);
  }
  
  return;
}

void acol_manager::command_del(std::string loc_type) {

  std::map<std::string,std::vector<std::string> >::iterator it;
  for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
    if (it->first==loc_type) {
      std::vector<std::string> &clist=it->second;
      for(size_t j=0;j<clist.size();j++) {
	cl->remove_comm_option(clist[j]);
      }
    }
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

  cl->cmd_name="acol";

  const int cl_param=cli::comm_option_cl_param;
  const int both=cli::comm_option_both;

  static const int narr=21;

  string type_list_str;
  for(size_t i=0;i<type_list.size()-1;i++) {
    type_list_str+=type_list[i]+", ";
  }
  type_list_str+="or "+type_list[type_list.size()-1]+'.';
  
  // Options, sorted by long name. We allow 0 parameters in many of these
  // options so they can be requested from the user in interactive mode. 
  comm_option_s options_arr[narr]=
    {{0,"autocorr","Compute the autocorrelation coefficients.",0,-1,
      "[arguments depend on current object type.]",
      ((std::string)"The behavior of the autocorr command depends on ")+
      "the type of the current object.\n\nNumerical array: "+
      "(no arguments)\n\nReplace the current "+
      "object with a vector of doubles which contains the autocorrelation "+
      "coefficient as a function of the step size.\n\n"+
      "table: <ac> <ftom> <column or vec. spec> [column or vec. spec. 2]"+
      "\n\nThree arguments are required "+
      "A column name <ac>, a column name <ftom>, and arguments which "+
      "specify the data. The "
      "autocorrelation coefficients are stored in column <ac> and "+
      "the quantity '5*tau/M' is stored in "+
      "column <ftom>. The data may be either a column "+
      "in the table or a vector specification. "+
      "Columns <ac> and <ftom> are created "+
      "if they are not already present and overwritten if they "+
      "already contain data. Also, the autocorrelation length and "+
      "estimated sample size are output to the screen. If multiple "+
      "data sources are given, then the autocorrelation coefficients "+
      "are averaged together. See "+cl->cmd_name+" -help vector-spec"+
      " for help on multiple vector specifications.\n\n"+
      "no object: <mult. vec. spec. 1> "+
      "[mult. vec. spec 2] "+
      "...\n\n Compute the autocorrelation coefficient for all vectors "+
      "specified as arguments then average those autocorrelation "+
      "coefficients together. If there is no current object "+
      "then the averaged autocorrelation coefficients are kept "+
      "as a double[] object. See "+cl->cmd_name+" -help mult-vector-spec"+
      " for help on multiple vector specifications.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_autocorr),
      both},
     {0,"calc","Compute the value of a constant expression.",0,1,"<expr>",
      ((string)"This computes the value of the constant expression ")+
      " <expr>. Examples are \"calc acos(-1)\" or \"calc 2+1/sqrt(2.0e4)\". "+
      "Results are given at the current precision. To see valid "+
      "expressions type \""+cl->cmd_name+" -help functions\".",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_calc),
      both},
     {0,"clear","Clear the current object.",0,0,"",
      "Deallocate the memory associated with the current object.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_clear),
      both},
     {'c',"create","Create an object.",0,-1,"<type> [...]",
      ((string)"Create a new object of type <type>. ")+
      "If an object is currently in memory, it is deallocated before "+
      "creating the new object.\n\n"+
      "\"create <type> <value>\": For types char, "+
      "int, size_t, and string, create an object and give it the "+
      "initial value specified.\n\n"+
      "\"create double <value spec.>\": Create a double object and set "+
      "it equal to the value specified by <value spec.>. "+
      "(See \"acol -help "+
      "functions\" for help on specifying functions and \"acol -help "+
      "value-spec\" for help on value specifications.)\n\n"+
      "\"create <type> <size> <function of \"i\">\": For array types "+
      "int[] and size_t[], the user must specify the size of "+
      "the array and a function of the array index 'i' to fill the array.\n\n"+
      "\"create double[] [<size> <function of \"i\">] or [vector spec.]\": "+
      "For double[] the user must either "+
      "give a vector specification, or specify the size of "+
      "the array and a function of the array index 'i'.\n\n"+
      "\"create table <name> <vector spec.>\": "+
      "Create a new table object with one column named <name> "+
      "from a vector specification (see \"acol -help vector-spec\" "+
      "for the syntax).\n\n"+
      "\"create tensor <rank> <size 0> <size 1> ...\": Create a tensor "+
      "object with the specified rank and sizes. All tensor entries "+
      "are initialized to zero.\n\n"+
      "\"create tensor_grid <rank> <size 0> <size 1> ...\": Create a "+
      "tensor_grid object with the specified rank and sizes. The tensor "+
      "grid is initialized to count each index (beginning with zero) and "+
      "the entries of the tensor "+
      "are initialized to zero. The grid can be specified afterwards "+
      "using \"set-grid\".\n\n"+
      "\"create table3d <x name> <x vector spec.> <y name> "+
      "<y vector spec.>\n  <slice name> <slice function>\": Create "+
      "a new table3d object which has one slice. The x and y grids "+
      "are given as vector specifications (see \"acol -help vector-spec\" "+
      "for the syntax). The slice function can be written in terms "+
      "of the x- and y-grid values which are referred to by name.\n\n",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
      both},
     {0,"docs","Open local HTML docs for acol or an O2scl topic.",0,1,
      "[topic]",
      ((string)"If [topic] is unspecified, this command opens up the ")+
      "local HTML documentation for acol in the default web browser "+
      "using 'open' on OSX and "+
      "'xdg-open' on other systems. If a topic is specified, then "+
      "the closest O2scl documentation web page is opened. In order "+
      "to open the remote version of the documentation instead of "+
      "the local copy, use 'wdocs' instead.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_docs),
      both},
     {0,"wdocs","Open remote HTML docs for acol or an O2scl topic.",0,2,
      "[search_term], [topic] or [section search_term]",
      ((string)"If no arguments are given, this command opens up the ")+
      "remote HTML documentation for acol in the default web browser "+
      "using 'open' on OSX and "+
      "'xdg-open' on other systems. If a [topic] is specified, then "+
      "the associated O2slcl web page is opened. If the argument "+
      "does not match an already known topic, then the search feature "+
      "on the O2scl web page is opened using the specified search term. "+
      "Note that, for search terms, spaces can be included using e.g. "+
      "'-wdocs \"Simulated annealing\"'. Valid sections are either "+
      "\"eos\" or \"part\". In order "+
      "to open the local version of the documentation instead of "+
      "the remote copy, use 'docs' instead of 'wdocs'.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_wdocs),
      both},
     {0,"download","Download file from specified URL.",0,4,
      "<file> <URL> [hash, \"file:\"hash_filename, or \"none\"] [directory]",
      ((string)"Check if a file matches a specified hash, and if not, ")+
      "attempt to download a fresh copy from the specified URL. If the "+
      "filename is \"_\", then the file is extracted from the end of "+
      "the URL.",
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
      "type specified.\n\nFor int, char, double, or size_t objects, "+
      "the file is assumed to begin with the desired object and it is "+
      "read using operator>>().\n\nFor string objects, the first line is "
      "read using std::getline().\n\nFor array objects, it is assumed "+
      "that all array entries are on the first line of the file and no "+
      "carriage returns are present between entries.\n\nFor table objects, "
      "the first line of the file must either contain numeric "+
      "data or column names "+
      "separated by white space, without carriage returns, except for "+
      "the one at the end of the line. If the first line contains "+
      "column names, the second line may optionally contain unit "+
      "expressions for each column, enclosed by square brackets. "+
      "All remaining lines are assumed "+
      "to contain data with the same number of columns as the first line. "+
      "\n\nFor table3d objects, the data must be stored in columns "+
      "with the first column specifying the x-axis grid point and "+
      "the second column specifying the y-axis grid point. The "+
      "remaining columns give the data for each slice at that point. "+
      "Each grid point must correspond to a row in the file, but "+
      "the lines need not be in any particular order. The columns may "+
      "have one header line at top which specifies the names of the x- "+
      "and y-grids and the names of each slice (in order).",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
      both},
     {0,"get-conv","Get a unit conversion factor.",0,3,
      "<old unit> <new unit> [value to convert]",
      ((string)"This command gets a unit ")+
      "conversion factor and optionally applies than conversion factor "+
      "to a user-specified value. It only works if the conversion is one of "
      "the hard-coded O2scl conversions or if HAVE_POPEN is defined and "+
      "the 'units' command is available in the current "+
      "path. For example, 'get-conv MeV erg' returns '1.000000e+00 MeV = "+
      "1.602177e-6 erg. The conversion factor is output "+
      "at the current precision, but is always internally stored with "+
      "full double precision. O2scl has several unit conversions which "+
      "implicitly assume hbar=c=1.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_conv),
      both},
     {0,"h5-copy","Copy hdf5 file (experimental).",-1,-1,
      "<source> <destination>",((string)"Copy all O2scl objects from ")+
      "one HDF5 file to another. This may not work for HDF5 files "+
      "generated outside of O2scl. The source and destination filenames "+
      "may not be identical. The destination file may not be the same "+
      "size as the source, but will contain the same information.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_h5_copy),
      both},
     {0,"constant","Get constant.",0,2,"<name or pattern> [unit]",
      "Get a hard-coded constant from the library.",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_constant),
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
     {0,"slack","Send a slack message.",0,2,"[\"#channel\"] <strings-spec>",
      ((string)"Send a message to slack, using the specified channel. ")+
      "If the channel is not specified, it is taken from the "+
      "environment variable O2SCL_SLACK_CHANNEL. The '#' sign "+
      "should be included with the channel name. "+
      "The Slack webhook URL is taken "+
      "from the environment variable O2SCL_SLACK_URL and the username "+
      "is taken from the environment variable O2SCL_SLACK_USERNAME. "+
      "The message is constructed from the string list specification "+
      "in <strings-spec> (see 'acol -help strings-spec' for more "+
      "information).",
      new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_slack),
      both},
     {0,"type","Show current object type.",0,0,"",
      ((string)"Show the current object type, either table, ")+
      type_list_str,
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

  terminal ter;
  cl->desc=((string)"acol: A data viewing and processing ")+
    "program for "+ter.bold()+"O2scl"+ter.default_fg()+".\n";

#ifdef O2SCL_NEVER_DEFINED

  // This code is now moved to acolm_gi.cpp
  
  ostringstream oss;
  oss << ((char)27) << '(' << '0';
  for(size_t i=0;i<78;i++) oss << 'q';
  oss << ((char)27) << '(' << 'B';
  string line=oss.str();
  
  string stemp;
  string dsc=line+"\nNotes:\n\n";
  vector<std::string> sv;
  
  stemp="1. Help for general commands may be obtained with '-help ";
  stemp+="<command>'. Help for type-specific commands can be obtained ";
  stemp+="by '-help <type> <command>'. A list of commands for each type ";
  stemp+="can be obtained with '-commands <type>', or for all commands ";
  stemp+=" use '-commands all'. Required arguments ";
  stemp+="are surrounded by ";
  stemp+="<>'s and optional arguments are surrounded by []'s.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }
  
  stemp="2. Options may also be specified in the environment variable ";
  stemp+=env_var_name+".\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="3. Long options may be preceeded by two dashes.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="4. In order to avoid confusion between arguments and functions, ";
  stemp+="use parenthesis and quotes, i.e. \"(-x*2)\" instead of -x*2.\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  stemp="5. Also, do not use a unary minus next to a binary operator, ";
  stemp+="i.e. use \"a>(-1)\" instead of \"a>-1\".\n\n";
  rewrap_ignore_vt100(stemp,sv,76);
  dsc+=sv[0]+"\n";
  for(size_t j=1;j<sv.size();j++) {
    dsc+="   "+sv[j]+"\n";
  }

  dsc+=line+"\n";
  
  dsc+="List of additional help topics (e.g. \"acol -help <topic>\"): ";
  dsc+="functions, mult-vector-spec, types, value-spec, and vector-spec.\n\n";
  
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

#endif
  
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
    "environment varable "+env_var_name+".";
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
  comm_option_s options_arr2[narr2]=
    {
     {'h',"help","Show help information.",0,2,
      "[type command] or [command] or [topic]",
      ((std::string)"Show generic help information, or, if an ")+
      "argument is given "+
      "give the documentation for the specified command or topic. "+
      "If two arguments are given, show the help for a type-specific "+
      "command. "+
      "Note that required arguments are typically given inside "+
      "angled brackes <> while optional arguments are given "+
      "inside square brackets [].",
      new o2scl::comm_option_mfptr<acol_manager>
      (this,&acol_manager::comm_help),both},
     {0,"commands",
      "List available commands for current or specified type.",
      0,1,"[type]",((string)"Output the commands available for ")+
      "the current type, or, if the optional type argument is given "+
      "then output the commands available for that type.",
      new o2scl::comm_option_mfptr<acol_manager>
      (this,&acol_manager::comm_commands),both}
    };
  cl->set_comm_option_vec(narr2,options_arr2);

  //-------------------------------------------------------------------
  // Try to get screen width
  
  int nrow, ncol=80;
  // Use curses
  // AWS: this prints out a bunch of crazy characters on Ubuntu
  //get_screen_size_curses(nrow,ncol);
  //#else
  
  // If not, attempt to obtain the result from the environment
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
  
  //#endif
  
  set_swidth(ncol);

  //-------------------------------------------------------------------
  // Setup parameters modified by 'set' and 'get'

  if (verbose>2) {
    cout << "Setup parameters: " << endl;
  }
  setup_parameters();

  //-------------------------------------------------------------------
  // Process default options and call

  std::vector<cmd_line_arg> ca;
  
  char *dc=getenv(env_var_name.c_str());
  if (dc) {
    if (verbose>2) {
      cout << "Process default options: " << dc << endl;
    }
    def_args=dc;
    if (verbose>2) {
      cl->process_args_str(def_args,ca,1,true);
    } else {
      cl->process_args_str(def_args,ca,0,true);
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
      cl->process_args_c(argc,argv,ca,1,true);
    } else {
      cl->process_args_c(argc,argv,ca,0,true);
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
    for(size_t i=0;i<sv.size()-1;i++) {
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

