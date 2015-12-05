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
#include "acolm.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

acol_manager::acol_manager() : cng(o2scl_settings.get_convert_units()) {
  table_name="acol";
  tabp=0;
  t3p=0;
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

  threed=false;

  env_var_name="ACOL_DEFAULTS";
}

int acol_manager::setup_options() {

  const int cl_param=cli::comm_option_cl_param;
  const int both=cli::comm_option_both;

  static const int narr=41;

  // Options, sorted by long name. We allow 0 parameters in many of these
  // options so they can be requested from the user in interactive mode. 
  comm_option_s options_arr[narr]={
    {'a',"assign","Assign a constant, e.g. assign pi acos(-1)",
     0,2,"<name> [val]",
     ((string)"Assign a constant value to a name for the present table. ")+
     "Valid constant values are things like 1.618 or acos(-1.0) or sin(4^5). "
     "To remove an assignment, call assign with a blank value.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_assign),
     both},
    {0,"calc","Compute the value of a constant expression.",0,1,"<expr>",
     ((string)"This computes the value of the constant expression ")+
     " <expr>. Examples are 'calc acos(-1)' or 'calc 2+1/sqrt(2.0e4)'. "+
     "Results are given at the current precision.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_calc),
     both},
    {'c',"create","Create a table from uniform grid.",
     0,4,"<name> <low> <hi> <step>",
     ((string)"Create a new table with one column whose entries ")+
     "are an evenly-spaced grid. This takes four arguments, the name of "+
     "the column, the first value, the increment between successive values "+
     "and the maximum possible value. Note that finite precision "+
     "arithmetic may cause small deviations from the expected result. "+
     "If a table is currently in memory, it is deallocated beforehand. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_create),
     both},
    {0,"delete-col","Delete a column (2d only).",0,1,"<name>",
     "Delete the entire column named <name>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_delete_col),
     both},
    {'d',"delete-rows","Delete rows selected by a function.",
     0,1,"<func>",((string)"Delete the set of rows for ")+
     "which a function evaluates to a number greater than 0.5. "+
     "For example, 'delete-rows if(col1+col2>10,1,0)' will delete "+
     "all columns where the sum of the entries in 'col1' and 'col2' "+
     "is larger than 10 (2d only).",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_delete_rows),
     both},
    {'D',"deriv",
     "Derivative of a function defined by two columns (2d only).",
     0,3,"<x> <y> <name>",
     ((string)"Create a new column named <name> filled with the ")+
     "derivative of the function y(x) obtained from columns <x> and <y>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv),
     both},
    {0,"deriv2","Second derivative (2d only).",0,3,"<name> <x> <y>",
     ((string)"Create a new column named <name> filled with the second ")+
     "derivative of the function y(x) obtained from columns <x> and <y>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_deriv2),
     both},
    {0,"filelist","List objects in a HDF5 file.",0,1,"<file>",
     ((string)"This lists all the top-level datasets and groups in a ")+
     "HDF5 file and, for those groups which are in the O2scl format, "+
     "gives the type and name of the object stored in that HDF5 group.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_filelist),
     both},
    {0,"find-row","Find a row which maximizes a function (2d only).",
     0,2,"<func> or find-row <col> <val>",
     ((string)"If one argument is given, then find-row finds the row ")+
     "which maximizes the value of the "+
     "expression given in <func>, and then output the entire row. "+
     "Otherwise find-row finds the row for which the value in "+
     "column named <col> is as close as possible to the value <val>. "+
     "See command 'get-row' to get a row by it's index.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_find_row),
     both},
    {0,"fit","Fit two columns to a function (experimental, 2d only).",0,7,
     "<x> <y> <yerr> <ynew> <par names> <func> <vals>","",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_fit),
     both},
    {'f',"function","Create a new column or slice from a function.",0,2,
     "<func> <name>",
     ((string)"Create a new column named <name> from a function (in ")+
     "<func>) in terms of the other columns. For example, for "+
     "a table containing columns named 'c1' and 'c2', 'function "+
     "c1-c2 c3' would create a new column c3 which contains the "+
     "difference of columns 'c1' and 'c2'.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_function),
     both},
    {'g',"generic","Read in a generic data file (2d only).",0,1,"<file>",
     ((string)"Read a generic data file with the given filename. ")+
     "The first line of the file is assumed to contain column names "+
     "separated by white space, without carriage returns, except for "+
     "the one at the end of the line. All remaining lines are assumed "+
     "to contain data. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_generic),
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
    {'H',"html","Create a file in HTML (2d only).",0,1,"<file>",
     "Output the current table in HTML mode to the specified file. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_html),
     both},
    {'N',"index","Add a column containing the row numbers (2d only).",0,1,
     "[column name]",
     ((string)"Define a new column named [column name] and fill ")+
     "the column with the row indexes, beginning with zero. If "+
     "no argument is given, the new column is named 'N'.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_index),
     both},
    {0,"insert","Interpolate a column from another file.",0,6,
     ((string)"2D: <file> <table> <oldx> <oldy> <newx> [newy],\n\t\t")+
     "3D: <file> <table> <old> [new]",
     ((string)"Insert a column from file <fname> interpolating it ")+
     "into the current table. The column <oldy> is the "+
     "columns in the file which is to be inserted into the table, "+
     "using the column <oldx> in the file and <newx> in the table. "+
     "The new column in the table is named <oldy>, or it is named "+
     "[newy] if the additional argument is given. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_insert),
     both},
    {0,"insert-full","Interpolate a table from another file (2d only).",0,3,
     "<fname> <oldx> <newx>",
     ((string)"Insert all columns from file <fname> interpolating it ")+
     "into the current table. The column <oldy> is the "+
     "columns in the file which is to be inserted into the table, "+
     "using the column <oldx> in the file and <newx> in the table. "+
     "The new column are given the same names as in the file. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
     both},
    {'I',"integ","Integrate a function specified by two columns (2d only).",
     0,3,"<x> <y> <name>",
     ((string)"Create a new column named <name> filled with the ")+
     "integral of the function y(x) obtained from columns <x> and <y>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_integ),
     both},
    {'q',"interactive","Load the interactive interface.",
     0,0,"",((string)"If given as a command-line parameter, 'interactive' ")+
     "toggles the execution of the interactive mode after the "+
     "command-line parameters are processed. If zero arguments are given "+
     "to 'acol' on the command-line then the interactive interface is "+
     "automatically turned on.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interactive),
     cl_param},
    {'i',"internal","Output in the internal HDF5 format.",0,1,"[file]",
     ((string)"Output the current table in the internal HDF5 format. ")+
     "If no argument is given, then output is sent to the screen, "+
     "otherwise, output is sent to the specified file. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_internal),
     both},
    {0,"interp","Interpolate a number into a column or slice.",0,3,
     "2d: <x name> <x value> <y name>, 3d: <z name> <x value> <y value> ",
     ((string)"For a 2d table, interpolate <x value> from column ")+
     "named <x name> into column named <y name>. For a 3d table "+
     "interpolate (<x value>,<y value>) into the slice named <z name>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp),
     both},
    {'l',"list","List the constants, column/slice names and other info.",
     0,0,"","",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_list),
     both},
    {0,"max","Find the maximum value of a column or slice.",0,1,"<col>",
     "Compute the maximum value of column <col>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_max),
     both},
    {0,"min","Find the minimum value of a column or slice.",0,1,"<col>",
     "Compute the minimum value of column <col>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_min),
     both},
    {'o',"output","Output the current table.",0,1,"[file]",
     ((string)"Output the table to the screen, or if the [file] ")+
     "argument is specified, to a file. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_output),
     both},
    {0,"interp-type","Get/set the current interpolation type.",0,1,"[type]",
     ((string)"Get or set the current object's interpolation type. ")+
     "Values are 0 (linear), 1 (cubic spline), 2 (cubic spline, periodic) "+
     "3 (Akima spline), and 4 (Akima spline, periodic).",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_interp_type),
     both},
    {'P',"preview","Preview the current table",0,1,"[nlines]",
     ((string)"Print out [nlines] lines of data for as many columns as ")+
     "will fit on the screen. The value of [nlines] defaults to 10. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_preview),
     both},
    {'r',"read","Read a table or table3d from a file.",0,2,
     "<file> [table name]",
     ((string)"Read the internally-formatted (either text or binary) ")+
     "file with the specified filename and make it the current table. " +
     "If the [table name] argument is specified, then read the table "+
     "with the specified name in a file with more than one table.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_read),
     both},
    {0,"rename","Rename a column or slice.",0,2,"<old> <new>",
     "Rename a column from <old> to <new>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_rename),
     both},
    {'s',"select","Select columns or slices for a new table.",-1,-1,"<cols>",
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
    {0,"select-rows","Select rows for a new table (2d only).",
     0,1,"<row_spec>","",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_select_rows),
     both},
    {0,"set-data","Set the entries of a column.",3,4,
     "2d: <row_spec> <col> <val_spec> 3d: <x value> <y value> <z name> <val>",
     ((string)"For a 2d table, sfet the value of rows specifed by the ")+
     "'row_spec' function in column 'col' to the value given by the "+
     "'val_spec' function. Rows are chosen if row_spec evaluates to a "+
     "number greater than 0.5. For a 3d table, just set the value of "+
     "the slice named 'z name' at the grid point closest to "+
     "(<x value>,<y value>) to the value <val>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_data),
     both},
    {0,"get-unit","Get the units for a specified column.",0,1,"<column>","",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_unit),
     both},
    {0,"set-unit","Set the units for a specified column.",0,2,
     "<column> <unit>","",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_set_unit),
     both},
    {0,"show-units","Show the unit conversion table.",0,0,"",
     ((string)"(This doesn't show all possible conversions, only ")+
     "the conversions which have been previously used and are now stored "+
     "in the unit cache.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_show_units),
     both},
    {0,"convert-unit","Convert a column to a new unit",0,2,
     "<column> <new_unit>",((string)"(This command only works if ")+
     "the GNU 'units' command is installed and available in the current "+
     "path.) Convert the units of a column to <new unit>, multipliying "+
     "all entries in that column by the appropriate factor.",
     new comm_option_mfptr<acol_manager>
     (this,&acol_manager::comm_convert_unit),both},
    {0,"get-conv","Get a unit conversion factor.",0,2,
     "<old unit> <new unit>",((string)"(This command only works if ")+
     "the GNU 'units' command is installed and available in the current "+
     "path.) For example, 'get-conv MeV erg' returns 1.602e-6 and 1 MeV "+
     "is equivalent to 1.602e-6 erg. The conversion factor is output "+
     "at the current precision, but is always internally stored with "+
     "full double precision.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_get_conv),
     both},
    {'S',"sort","Sort the entire table by a column (2d only).",0,1,"<col>",
     "Sorts the entire table by the column specified in <col>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_sort),
     both},
    {0,"stats","Show column statistics (2d only).",0,1,"<col>",
     "Output the average, std. dev, max and min of <col>. ",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_stats),
     both},
    {'v',"version","Print version information.",0,0,"",
     "",new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_version),
     both},
    {0,"add","Add the data from two table3d objects.",0,3,
     "<file 1> <file 2> <sum file>",((string)"Read the two equally-")+
     "sized table3d objects in <file 1> and <file 2> and construct a "+
     "new table where each data slice represents the sum of the data "+
     "from each of the original files. Store this sum in <sum file>.",
     new comm_option_mfptr<acol_manager>(this,&acol_manager::comm_add),
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

  //---------------------------------------
  // Get HOME directory and command history
  
  char *hd=getenv("HOME");
  std::string histfile;
  histfile=hd;
  histfile+="/.acol_hist";
  
  //---------------------------------------
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
  
  cl->desc=((string)"acol: A data table viewing and ")+
    "processing program for O2scl.\n";
  
  string dsc="\nNotes:\n\n";
  dsc+="1. Help for individual commands may be obtained with 'help ";
  dsc+="command'.\n   Required arguments are surrounded by <>'s and ";
  dsc+="optional arguments are\n   surrounded by []'s.\n\n";
  dsc+="2. Options may also be specified in the environment variable ";
  dsc+="ACOL_DEFAULTS.\n\n";
  dsc+="3. Long options may be preceeded by two dashes.\n\n";
  dsc+="4. In order to avoid confusion between arguments and functions,\n";
  dsc+="   use \"(-x*2)\" not \"-x*2\"\n\n";
  dsc+="5. The select command assumes that the first equals sign indicates\n";
  dsc+="   renaming, so use \"eq=if(x=3,1,0)\" not \"if(x=3,1,0)\".\n\n";
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

int acol_manager::setup_parameters() {
  
  p_table_name.str=&table_name;
  p_unit_fname.str=&unit_fname;
  p_def_args.str=&def_args;
  p_verbose.i=&verbose;
  p_prec.i=&prec;
  p_ncols.i=&ncols;
  p_scientific.b=&scientific;
  p_pretty.b=&pretty;
  p_names_out.b=&names_out;
  
  p_table_name.help="The current table name";
  p_unit_fname.help="The unit filename";
  p_def_args.help="The default arguments from the environment";
  p_prec.help="The numerical precision";
  p_verbose.help="Control the amount of output";
  p_ncols.help="The number of output columns";
  p_names_out.help="If true, output column names at top";
  p_pretty.help="If true, align the columns using spaces";
  p_scientific.help="If true, output in scientific mode";
  
  cl->par_list.insert(make_pair("table_name",&p_table_name));
  cl->par_list.insert(make_pair("unit_fname",&p_unit_fname));
  cl->par_list.insert(make_pair("def_args",&p_def_args));
  cl->par_list.insert(make_pair("precision",&p_prec));
  cl->par_list.insert(make_pair("verbose",&p_verbose));
  cl->par_list.insert(make_pair("ncols",&p_ncols));
  cl->par_list.insert(make_pair("names_out",&p_names_out));
  cl->par_list.insert(make_pair("pretty",&p_pretty));
  cl->par_list.insert(make_pair("scientific",&p_scientific));

  return 0;
}

int acol_manager::run(int argv, char *argc[]) {
  
  //---------------------------------------
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

  // 
  comm_option_mfptr<acol_manager> cset(this,&acol_manager::comm_set);
  cl->set_function(cset);

  //----------------------------------------------------------------
  // Process default options

  if (verbose>2) {
    cout << "Process default options" << endl;
  }
  std::vector<cmd_line_arg> ca;
  
  char *dc=getenv(env_var_name.c_str());
  if (dc) {
    def_args=dc;
    if (verbose>2) {
      cl->process_args(dc,ca,1);
    } else {
      cl->process_args(dc,ca);
    }
  }
  
  //----------------------------------------------------------------
  // Process command-line options

  // Note that it's ok that this appears early in the code because it
  // just processes the arguments, it doesn't really try to do
  // any execution based on those arguments

  if (verbose>2) {
    cout << "Process command-line options" << endl;
    cl->process_args(argv,argc,ca,1);
  } else {
    cl->process_args(argv,argc,ca);
  }
  if (argv<2) {
    post_interactive=true;
  }

  //----------------------------------------------------------------
  // Try to get screen width
  
  int ncol=80;
  char *ncstring=getenv("COLUMNS");
  if (ncstring) ncol=o2scl::stoi(ncstring);
  
  set_swidth(ncol);

  if (verbose>2) {
    cout << "Setup parameters: " << endl;
  }
  //
  setup_parameters();

  //----------------------------------------------------------------
  // Main execution
  
  int ret=0, ret2=0;

  if (ca.size()>0) {
    if (verbose>2) {
      cout << "Process default options and command-line arguments." << endl;
    }
    ret=cl->call_args(ca);
  }

  if (verbose>2) {
    cout << "Post_interactive: " << post_interactive << endl;
  }
  
  if (post_interactive) {
    if (verbose>2) {
      cout << "Run interactive mode." << endl;
    }
    ret2=cl->run_interactive();
  }

  //---------------------------------------
  // Notify user if error occurred

  if (ret!=0 || ret2!=0) {
    cout << "An error occured." << endl;
  }

  delete cl;
  
  return 0;

}

int acol_manager::comm_interp_type(std::vector<std::string> &sv, 
				   bool itive_com) {

  if (threed) {
    
    if (t3p==0) {
      cout << "No table to get interpolation type of." << endl;
      return exc_efailed;
    }

    if (sv.size()>1) {
      if (o2scl::stoi(sv[1])>4 || o2scl::stoi(sv[1])<0) {
	cout << "Invalid interpolation type in interp-type." << endl;
	return exc_efailed;
      }
      t3p->set_interp_type(o2scl::stoi(sv[1]));
    }
    
    if (sv.size()==1 || verbose>0) {
      size_t itype=t3p->get_interp_type();
      cout << "Current interpolation type is " << itype;
      if (itype<2) {
	if (itype==0) {
	  cout << " (linear)." << endl;
	} else {
	  cout << " (cubic spline)." << endl;
	}
      } else {
	if (itype==2) {
	  cout << " (cubic spline, periodic)." << endl;
	} else if (itype==3) {
	  cout << " (Akima spline)." << endl;
	} else {
	  cout << " (Akima spline, periodic)." << endl;
	}
      }
    }
    
    return 0;
  }
  
  if (tabp==0) {
    cout << "No table to get interpolation type of." << endl;
    return exc_efailed;
  }

  if (sv.size()>1) {
    if (o2scl::stoi(sv[1])>4 || o2scl::stoi(sv[1])<0) {
      cout << "Invalid interpolation type in interp-type." << endl;
      return exc_efailed;
    }
    tabp->set_interp_type(o2scl::stoi(sv[1]));
  }

  if (sv.size()==1 || verbose>0) {
    size_t itype=tabp->get_interp_type();
    cout << "Current interpolation type is " << itype;
    if (itype<2) {
      if (itype==0) {
	cout << " (linear)." << endl;
      } else {
	cout << " (cubic spline)." << endl;
      }
    } else {
      if (itype==2) {
	cout << " (cubic spline, periodic)." << endl;
      } else if (itype==3) {
	cout << " (Akima spline)." << endl;
      } else {
	cout << " (Akima spline, periodic)." << endl;
      }
    }
  }
  
  return 0;
}

int acol_manager::comm_output(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    
    if (t3p==0) {
      cerr << "No table3d to output." << endl;
      return exc_efailed;
    }

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

    size_t nx, ny;
    t3p->get_size(nx,ny);
    if (nx!=0 && ny!=0) {

      (*fout) << "Grid x: ";
      for(size_t i=0;i<((size_t)nx);i++) {
	(*fout) << t3p->get_grid_x(i) << " ";
      }
      (*fout) << endl;
      (*fout) << "Grid y: ";
      for(size_t i=0;i<((size_t)ny);i++) {
	(*fout) << t3p->get_grid_y(i) << " ";
      }
      (*fout) << endl;

      size_t nt=t3p->get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  (*fout) << "Slice " << k << ": " << t3p->get_slice_name(k) << endl;
	  if (k==0) {
	    (*fout) << "Outer loops over x grid, inner loop over y grid." 
		    << endl;
	  }
	  for(size_t i=0;i<nx;i++) {
	    for(size_t j=0;j<ny;j++) {
	      (*fout) << t3p->get(i,j,k) << " ";
	    }
	    (*fout) << endl;
	  }
	  fout->unsetf(ios::showpos);
	}
      }
    }

    ffout.close();

    return 0;
  }

  if (tabp==0) {
    cerr << "No table to output." << endl;
    return exc_efailed;
  }

  //---------------------------------------
  // Create stream

  ostream *fout;
  ofstream ffout;

  if (sv.size()==1) {
    fout=&cout;
  } else {
    ffout.open(sv[1].c_str());
    fout=&ffout;
  }

  //---------------------------------------
  // Output formatting

  if (scientific) fout->setf(ios::scientific);
  else fout->unsetf(ios::scientific);
  fout->precision(prec);

  if (tabp->get_ncolumns()>0) {

    //---------------------------------------
    // Count column widths

    size_t *col_wids=new size_t[tabp->get_ncolumns()];

    for(size_t i=0;i<tabp->get_ncolumns();i++) {
      col_wids[i]=prec+6;
    }

    if (names_out==true) {
      for(size_t i=0;i<tabp->get_ncolumns();i++) {
	if (tabp->get_column_name(i).size()>col_wids[i]) {
	  col_wids[i]=tabp->get_column_name(i).size();
	}
	if (tabp->get_unit(tabp->get_column_name(i)).size()+2>col_wids[i]) {
	  col_wids[i]=tabp->get_unit(tabp->get_column_name(i)).size()+2;
	}
      }
    }

    //---------------------------------------
    // Output column names
      
    if (names_out==true) {
      for(size_t i=0;i<tabp->get_ncolumns();i++) {
	  
	// Preceeding space
	if (pretty==true) {
	  (*fout) << ' ';
	}
	  
	// Column name
	(*fout) << tabp->get_column_name(i) << " ";
	  
	// Trailing spaces
	if (pretty==true) {
	  int nsp=col_wids[i]-tabp->get_column_name(i).size();
	  if (nsp<0) {
	    O2SCL_ERR("Column size anomaly (col names).",exc_efailed);
	  }
	  for(int j=0;j<nsp;j++) (*fout) << ' ';
	}
	
      }
      (*fout) << endl;
    }
      
    //---------------------------------------
    // Output units
    
    if (tabp->get_nunits()>0 && names_out==true) {
      for(size_t i=0;i<tabp->get_ncolumns();i++) {

	// Preceeding space
	if (pretty==true) {
	  (*fout) << ' ';
	}
	
	// Unit name
	(*fout) << '[' << tabp->get_unit(tabp->get_column_name(i)) << "] ";
	
	// Trailing spaces
	if (pretty==true) {
	  int nsp=col_wids[i]-
	    tabp->get_unit(tabp->get_column_name(i)).size()-2;
	  if (nsp<0) {
	    O2SCL_ERR("Column size anomaly (units).",exc_efailed);
	  }
	  for(int j=0;j<nsp;j++) (*fout) << ' ';
	}
	
      }
      (*fout) << endl;
    }
      
    //---------------------------------------
    // Output data
      
    for(int i=0;i<((int)tabp->get_nlines());i++) {
	
      for(size_t j=0;j<tabp->get_ncolumns();j++) {
	  
	// Otherwise, for normal output
	if (pretty==true && tabp->get(j,i)>=0.0) {
	  (*fout) << ' ';
	}
	(*fout) << tabp->get(j,i) << ' ';
	if (pretty==true) {
	  int nsp=((int)(tabp->get_column_name(j).size()-prec-6));
	  for(int kk=0;kk<nsp;kk++) {
	    (*fout) << ' ';
	  }
	}
	
      }
      (*fout) << endl;
      
      //---------------------------------------
      // Continue to next row
    }
    
  }

  ffout.close();
  
  return 0;
}

int acol_manager::comm_add(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<6) {
    cout << "Not enough arguments to add." << endl;
    return -1;
  }

  string s1=sv[1];
  string name1=sv[2];
  string s2=sv[3];
  string name2=sv[4];
  string sum=sv[5];

  table3d *t1=new table3d;
  table3d *t2=new table3d;
  table3d tsum;

  hdf_file hf;
  hf.open(s1);
  hdf_input(hf,*t1,name1);
  hf.close();
  hf.open(s2);
  hdf_input(hf,*t2,name2);
  hf.close();
  
  if (t1->get_nx()!=t2->get_nx() ||
      t1->get_ny()!=t2->get_ny()) {
    cout << "First table3d object. nx,ny: " 
	 << t1->get_nx() << " " << t1->get_ny() << endl;
    cout << "Second table3d object. nx,ny: " 
	 << t2->get_nx() << " " << t2->get_ny() << endl;
    cout << "Tables not compatible." << endl;
    return 0;
  }
  
  size_t nx=t1->get_nx();
  size_t ny=t1->get_ny();
  string sx=t1->get_x_name();
  string sy=t1->get_y_name();
  const ubvector &xg=t1->get_x_data();
  const ubvector &yg=t1->get_y_data();
  
  tsum.set_xy(sx,nx,xg,sy,ny,yg);
  for(size_t k=0;k<t1->get_nslices();k++) {
    string slname=t1->get_slice_name(k);
    tsum.new_slice(slname);
    for(size_t i=0;i<nx;i++) {
      for(size_t j=0;j<ny;j++) {
	tsum.set(i,j,slname,t1->get(i,j,slname)+t2->get(i,j,slname));
      }
    }
  }
  
  hf.open(sum);
  hdf_output(hf,tsum,name1);
  hf.close();
  
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

int acol_manager::comm_read(std::vector<std::string> &sv, bool itive_com) {
  
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

  // Delete previous table
  if (threed && t3p!=0) {
    threed=0;
    delete t3p;
    t3p=0;
  } else if (tabp!=0) {
    delete tabp;
    tabp=0;
  }

  // Use hdf_file to open the file
  hdf_file hf;
  string type;
  int ret;

  ret=hf.open(i1.c_str(),false);
  if (ret!=0) {
    cerr << "Couldn't find file named '" << i1 << "'. Wrong file name?" 
	 << endl;
    return exc_efailed;
  }

  // If a name was specified, look for an object with that name
  if (i2.length()>0) {
    if (verbose>2) {
      cout << "Looking for object with name '" << i2 << "'." << endl;
    }
    ret=hf.find_group_by_name(hf,i2,type,verbose);
    if (ret==exc_enotfound) {
      cout << "Could not find object named '" << i2 
	   << "' in file '" << i1 << "'." << endl;
      return exc_efailed;
    }
    if (verbose>2) {
      cout << "Found object with type '" << type << "'." << endl;
    }
    if (type=="table") {
      if (verbose>2) {
	cout << "Reading table." << endl;
      }
      tabp=new table_units<>;
      hdf_input(hf,*tabp,i2);
      table_name=i2;
      threed=false;
      return 0;
    } else if (type=="table3d") {
      if (verbose>2) {
	cout << "Reading table3d." << endl;
      }
      t3p=new table3d;
      hdf_input(hf,*t3p,i2);
      table_name=i2;
      threed=true;
      return 0;
    } else if (type=="hist") {
      if (verbose>2) {
	cout << "Reading hist." << endl;
      }
      threed=false;
      tabp=new table_units<>;
      hist h;
      hdf_input(hf,h,i2);
      if (verbose>0) {
	cout << "Creating a table from the histogram with columns named\n";
	cout << "'bins', 'low', 'high', and 'weights'." << endl;
      }
      h.copy_to_table(*tabp,"bins","low","high","weights");
      return 0;
    } else if (type=="hist_2d") {
      if (verbose>2) {
	cout << "Reading hist_2d." << endl;
      }
      threed=true;
      t3p=new table3d;
      hist_2d h;
      hdf_input(hf,h,i2);
      if (verbose>0) {
	cout << "Creating a table3d from the histogram with slice named\n";
	cout << "'weights'." << endl;
      }
      h.copy_to_table(*t3p,"x","y","weights");
      return 0;
    }
    cout << "Incompatible type '" << i2 << "'." << endl;
    return exc_efailed;
  } 

  if (verbose>2) {
    cout << "Looking for table." << endl;
  }
  ret=hf.find_group_by_type(hf,"table",i2,verbose);
  if (ret==success) {
    tabp=new table_units<>;
    hdf_input(hf,*tabp,i2);
    table_name=i2;
    threed=false;
    return 0;
  }

  if (verbose>2) {
    cout << "Looking for table3d." << endl;
  }
  ret=hf.find_group_by_type(hf,"table3d",i2,verbose);
  if (ret==success) {
    t3p=new table3d;
    hdf_input(hf,*t3p,i2);
    table_name=i2;
    threed=true;
    return 0;
  }

  if (verbose>2) {
    cout << "Looking for hist." << endl;
  }
  ret=hf.find_group_by_type(hf,"hist",i2,verbose);
  if (ret==success) {
    hist h;
    hdf_input(hf,h,i2);
    tabp=new table_units<>;
    if (verbose>0) {
      cout << "Creating a table from the histogram with columns named\n";
      cout << "'bins', 'low', 'high', and 'weights'." << endl;
    }
    h.copy_to_table(*tabp,"bins","low","high","weights");
    table_name=i2;
    threed=false;
    return 0;
  }
  
  if (verbose>2) {
    cout << "Looking for hist_2d." << endl;
  }
  ret=hf.find_group_by_type(hf,"hist_2d",i2,verbose);
  if (ret==success) {
    hist_2d h;
    hdf_input(hf,h,i2);
    t3p=new table3d;
    if (verbose>0) {
      cout << "Creating a table3d from the histogram with slice named\n";
      cout << "'weights'." << endl;
    }
    h.copy_to_table(*t3p,"x","y","weights");
    table_name=i2;
    threed=true;
    return 0;
  }

  cout << "Could not find object of any readable type in file '" << i1
       << "'." << endl;
  
  return exc_efailed;
}

herr_t acol_manager::filelist_func(hid_t loc, const char *name, 
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

    if (otype.length()!=0) {
      cout << "O2scl object named '" << name << "' of type " 
	   << otype << "." << endl;
    } else {
      cout << "Group named '" << name << "'." << endl;
    }

  } else if (infobuf.type==H5O_TYPE_DATASET) {
    cout << "Dataset named '" << name << "'." << endl;
  } else if (infobuf.type==H5O_TYPE_NAMED_DATATYPE) {
    cout << "Named type with name '" << name << "'." << endl;
  } else {
    cout << "Unexpected HDF type. " << endl;
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
  hf.open(i1.c_str());

  iter_parms ip={"",&hf,false,"",verbose};

  H5Literate(hf.get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	     0,filelist_func,&ip);
  
  return 0;
}

int acol_manager::get_input_one(vector<string> &sv, string directions,
				string &in, string comm_name,
				bool itive_com) {
  if (sv.size()>1) {
    in=sv[1];
    return 0;
  }
  if (itive_com) {
    string temp=directions+" (or blank to stop): ";
    in=cl->cli_gets(temp.c_str());
    if (in.length()==0) {
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

  if (sv.size()>ni) {
    for(size_t i=0;i<ni;i++) {
      in.push_back(sv[i+1]);
    }
    return 0;
  }
  if (itive_com) {
    for(size_t i=0;i<ni;i++) {
      string temp=directions[i]+" (or blank to stop): ";
      in.push_back(cl->cli_gets(temp.c_str()));
      if (in[i].length()==0) {
	cout << "Command '" << comm_name << "' cancelled." << endl;
	return exc_efailed;
      }
    }
  } else {
    cerr << "Not enough arguments to '" << comm_name << "'." << endl;
    return exc_efailed;
  }

  return success;
}

int acol_manager::comm_max(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {

    if (t3p==0 || t3p->get_nslices()==0) {
      cerr << "No table3d with slices to find the maximum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find maximum value of",
			  i1,"max",itive_com);
    if (ret!=0) return ret;

    double max=t3p->get(0,0,i1);
    for(size_t i=0;i<t3p->get_nx();i++) {
      for(size_t j=0;j<t3p->get_ny();j++) {
	if (t3p->get(i,j,i1)>max) max=t3p->get(i,j,i1);
      }
    }
    cout << "Maximum value of slice '" << i1 << "' is: " 
	 << max << endl;

    return 0;
  }
  
  if (tabp==0) {
    cerr << "No table with columns to find the maximum value of." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter column to find maximum value of",
			i1,"max",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(i1)==false) {
    cerr << "Couldn't find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  cout << "Maximum value of column '" << i1 << "' is: " 
       << tabp->max(i1) << endl;
  
  return 0;
}

int acol_manager::comm_min(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {

    if (t3p==0 || t3p->get_nslices()==0) {
      cerr << "No table3d with slices to find the minimum value of." << endl;
      return exc_efailed;
    }
    
    std::string i1;
    int ret=get_input_one(sv,"Enter slice to find minimum of",
			  i1,"min",itive_com);
    if (ret!=0) return ret;

    double min=t3p->get(0,0,i1);
    for(size_t i=0;i<t3p->get_nx();i++) {
      for(size_t j=0;j<t3p->get_ny();j++) {
	if (t3p->get(i,j,i1)<min) min=t3p->get(i,j,i1);
      }
    }
    cout << "Minimum value of slice '" << i1 << "' is: " 
	 << min << endl;

    return 0;
  }
  
  if (tabp==0) {
    cerr << "No table with columns to find the minimum value of." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter column to find minimum of",
			i1,"min",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(i1)==false) {
    cerr << "Couldn't find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  cout << "Minimum value of column '" << i1 << "' is: " 
       << tabp->min(i1) << endl;
  
  return 0;
}

int acol_manager::comm_set_data(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {

    if (t3p==0) {
      cerr << "No table3d with data to set." << endl;
      return exc_efailed;
    }
    
    vector<string> in, pr;
    pr.push_back(t3p->get_x_name()+" value of point to set");
    pr.push_back(t3p->get_y_name()+" value of point to set");
    pr.push_back("Slice name to set");
    pr.push_back("New value");
    int ret=get_input(sv,pr,in,"set-data",itive_com);
    if (ret!=0) return ret;
    
    t3p->set_val(o2scl::stod(in[0]),o2scl::stod(in[1]),in[2],
		 o2scl::stod(in[3]));
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to insert columns into." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Function describing rows to modify");
  pr.push_back("Name of column to modify");
  pr.push_back("Function describing new values");
  int ret=get_input(sv,pr,in,"set-data",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(in[1])==false) {
    cerr << "Couldn't find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  for(size_t i=0;i<tabp->get_nlines();i++) {
    if (tabp->row_function(in[0],i)>0.5) {
      tabp->set(in[1],i,tabp->row_function(in[2],i));
    }
  }

  return 0;
}

int acol_manager::comm_set_unit(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cerr << "Not implemented for 3d." << endl;
    return exc_efailed;
  }

  if (tabp==0) {
    cerr << "No table to set units of." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to set units of");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"set-unit",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  tabp->set_unit(in[0],in[1]);

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
  
  if (threed) {
    cerr << "Not implemented for 3d." << endl;
    return exc_efailed;
  }
  
  if (tabp==0) {
    cerr << "No table to convert units in." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column in which to convert units");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"convert-unit",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  if (unit_fname.length()>0) {
    cng.units_cmd_string=((string)"units -f ")+unit_fname;
  }
  ret=tabp->convert_to_unit(in[0],in[1],false);
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
  }
  
  // Set the proper output precision and mode
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  
  cout << "Conversion factor is: " << cng.convert(in[0],in[1],1.0) << endl;
  
  return 0;
}

int acol_manager::comm_get_unit(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cerr << "Not implemented for 3d." << endl;
    return exc_efailed;
  }

  if (tabp==0) {
    cerr << "No table to get units for." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column to get units of");
  int ret=get_input(sv,pr,in,"get-unit",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  cout << "Units of column " << in[0] << " are: " 
       << tabp->get_unit(in[0]) << endl;

  return 0;
}

int acol_manager::comm_find_row(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cerr << "Not implemented for table3d." << endl;
    return exc_efailed;
  }

  if (tabp==0 || tabp->get_nlines()==0) {
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
    
    if (tabp->is_column(sv[1])) {
      size_t row=tabp->lookup(sv[1],o2scl::stod(sv[2]));
      
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
  size_t row=tabp->function_find_row(sv[1]);
  
  // Call get_row() for the row that was found
  std::vector<std::string> sc;
  sc.push_back("get-row");
  sc.push_back(itos(row));
  comm_get_row(sc,itive_com);
  
  return 0;
}

int acol_manager::comm_get_row(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }
  
  if (tabp==0 || tabp->get_nlines()==0) {
    cerr << "No table or empty table in get-row." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter row number to get",
			i1,"get-row",itive_com);
  if (ret!=0) return ret;
  size_t ix=o2scl::stoi(i1);
  
  if (ix>tabp->get_nlines()-1) {
    cerr << "Requested row beyond end of table in get-row." << endl;
    return exc_efailed;
  }
  
  //---------------------------------------
  // Compute number of screen columns

  if (user_ncols<=0) {
    char *ncstring=getenv("COLUMNS");
    if (ncstring) ncols=o2scl::stoi(ncstring);
  } else {
    ncols=user_ncols;
  }

  //---------------------------------------
  // Temporary storage strings for names and data

  vector<string> row_names, row_data;

  //---------------------------------------
  // Process and/or output names

  if (names_out==true) {

    if (pretty) {

      size_t running_width=0;
      ostringstream *str=new ostringstream;
      str->setf(ios::scientific);
      str->precision(prec);
      
      for(size_t i=0;i<tabp->get_ncolumns();i++) {

	// Count for space between columns and sign
	size_t this_col=2;
	// Count column name
	this_col+=tabp->get_column_name(i).size();
	// Count extra spaces to format number
	int num_spaces=prec+6-((int)(tabp->get_column_name(i).size()));
	if (num_spaces>0) this_col+=num_spaces;
	// See if there will be space
	if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	  row_names.push_back(str->str());
	  delete str;
	  str=new ostringstream;
	  str->setf(ios::scientific);
	  str->precision(prec);
	  running_width=0;
	}
	// Output this column name
	(*str) << ' ' << tabp->get_column_name(i) << ' ';
	for(int j=0;j<num_spaces;j++) {
	  (*str) << ' ';
	}
	running_width+=this_col;
      }
      row_names.push_back(str->str());
      delete str;
      
    } else {
      
      cout.precision(prec);
  
      for(size_t i=0;i<tabp->get_ncolumns();i++) {
	cout << tabp->get_column_name(i) << ' ';
      }
      cout << endl;

    }
  }
  
  //---------------------------------------
  // Process and/or output data
  
  if (pretty) {
    
    size_t running_width=0;
    ostringstream *str=new ostringstream;
    str->setf(ios::scientific);
    str->precision(prec);
    
    for(size_t i=0;i<tabp->get_ncolumns();i++) {
      
      // Count space for number
      size_t this_col=prec+8;
      // Count extra spaces if necessary
      int num_spaces=((int)(tabp->get_column_name(i).size())-prec-6);
      if (num_spaces>0) this_col+=num_spaces;
      // See if there will be space
      if (running_width>0 && ((int)(running_width+this_col))>=ncols) {
	row_data.push_back(str->str());
	delete str;
	str=new ostringstream;
	str->setf(ios::scientific);
	str->precision(prec);
	running_width=0;
      }
      // Output the data
      if (tabp->get(i,ix)>=0.0) {
	(*str) << ' ' << tabp->get(i,ix) << ' ';
      } else {
	(*str) << tabp->get(i,ix) << ' ';
      }
      for(int j=0;j<num_spaces;j++) {
	(*str) << ' ';
      }
      running_width+=this_col;
    }
    row_data.push_back(str->str());
    delete str;
    
    //---------------------------------------
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
    
    for(size_t i=0;i<tabp->get_ncolumns();i++) {
      cout << tabp->get(i,ix) << ' ';
    }
    cout << endl;
    
  }
  
  return 0;
}

int acol_manager::comm_rename(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {

    vector<string> pr, in;
    pr.push_back("Enter slice to be renamed");
    pr.push_back("Enter new name");
    
    int ret=get_input(sv,pr,in,"rename",itive_com);
    if (ret!=0) return ret;
    
    t3p->rename_slice(in[0],in[1]);
    
    return 0;
  }

  vector<string> pr, in;
  pr.push_back("Enter column to be renamed");
  pr.push_back("Enter new name");
  
  int ret=get_input(sv,pr,in,"rename",itive_com);
  if (ret!=0) return ret;
  
  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  tabp->new_column(in[1]);

  tabp->copy_column(in[0],in[1]);
  tabp->delete_column(in[0]);

  return 0;
}

int acol_manager::comm_deriv(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
    cerr << "No table with columns to take derivatives of." << endl;
    return exc_efailed;
  }

  vector<string> pr, in;
  pr.push_back("Enter 'x' column");
  pr.push_back("Enter 'y' column");
  pr.push_back("Enter name of new column");
  int ret=get_input(sv,pr,in,"deriv",itive_com);
  if (ret!=0) return ret;

  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  if (tabp->is_column(in[1])==false) {
    cerr << "Couldn't find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  tabp->deriv(in[0],in[1],in[2]);

  return 0;
}

int acol_manager::comm_deriv2(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
    cout << "No table with columns to take derivatives of." << endl;
    return exc_efailed;
  }
  vector<string> pr, in;
  pr.push_back("Enter 'x' column");
  pr.push_back("Enter 'y' column");
  pr.push_back("Enter name of new column");
  int ret=get_input(sv,pr,in,"deriv2",itive_com);
  if (ret!=0) return ret;

  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  if (tabp->is_column(in[1])==false) {
    cerr << "Couldn't find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  tabp->deriv2(in[0],in[1],in[2]);

  return 0;
}

int acol_manager::comm_integ(std::vector<std::string> &sv, bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
    cerr << "No table with columns to integrate." << endl;
    return exc_efailed;
  }
  vector<string> pr, in;
  pr.push_back("Enter 'x' column");
  pr.push_back("Enter 'y' column");
  pr.push_back("Enter name of new column");
  int ret=get_input(sv,pr,in,"integ",itive_com);
  if (ret!=0) return ret;

  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  if (tabp->is_column(in[1])==false) {
    cerr << "Couldn't find column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }

  tabp->integ(in[0],in[1],in[2]);

  return 0;
}

int acol_manager::comm_internal(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    if (t3p==0) {
      cerr << "No table3d to write to a file." << endl;
      return exc_efailed;
    }

    std::string i1;
    int ret=get_input_one(sv,"Enter filename",i1,"internal",itive_com);
    if (ret!=0) return ret;

    if (verbose>2) {
      cout << "Creating O2scl file: " <<  i1 << endl;
    }

    hdf_file hf;
    hf.open_or_create(i1);
    table3d *tp=(table3d *)t3p;
    hdf_output(hf,*tp,table_name);
    hf.close();
  
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to write to a file." << endl;
    return exc_efailed;
  }

  std::string i1;
  int ret=get_input_one(sv,"Enter filename",i1,"internal",itive_com);
  if (ret!=0) return ret;

  if (verbose>2) {
    cout << "Creating O2scl file: " <<  i1 << endl;
  }

  hdf_file hf;
  hf.open_or_create(i1);
  hdf_output(hf,*tabp,table_name);
  hf.close();
  
  return 0;
}

int acol_manager::comm_function(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    
    if (t3p==0) {
      cerr << "No table3d to add a slice to." << endl;
      return exc_efailed;
    }
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
    
    t3p->function_slice(in[0],in[1]);
    if (ret!=0) {
      cerr << "Function make_slice() failed." << endl;
      return exc_efailed;
    }
    
    return 0;
  }

  if (tabp==0) {
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

  if (tabp->is_column(in[1])==true) {
    cerr << "Already a column named '" << in[1] << "'." << endl;
    return exc_efailed;
  }
  
  tabp->function_column(in[0],in[1]);
  
  if (ret!=0) {
    cerr << "Function make_col() failed." << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_calc(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

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
  double d=o2scl::function_to_double(i1,false);
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  if (verbose>0) cout << "Result: ";
  cout << d << endl;
  return 0;
}

int acol_manager::comm_index(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to add line numbers to." << endl;
    return exc_efailed;
  }

  std::string i1="N";
  if (sv.size()>1) i1=sv[1];
  tabp->new_column(i1);
  for(size_t i=0;i<tabp->get_nlines();i++) tabp->set(i1,i,((double)i));

  return 0;
}

int acol_manager::comm_html(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
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
  
  //---------------------------------------
  // Output constants
  
  for(size_t i=0;i<tabp->get_nconsts();i++) {
    if (i==0) (*fout) << "<b>Constants:</b><br />" << endl;
    string tnam;
    double tval;
    tabp->get_constant(i,tnam,tval);
    (*fout) << tnam << " = " << tval << "<br />" << endl;
  }
  
  if (tabp->get_ncolumns()>0) {
    
    (*fout) << "<table border=\"0\">" << endl;
    
    //---------------------------------------
    // Output column names
    
    (*fout) << "<tr bgcolor=\"#dddddd\">";
    for(size_t i=0;i<tabp->get_ncolumns();i++) {
      (*fout) << "<td><b>" << tabp->get_column_name(i) << "</b></td>";
    }
    (*fout) << "</tr>" << endl;
      
    //---------------------------------------
    // Output data
      
    for(int i=0;i<((int)tabp->get_nlines());i++) {
      if (i%5==4) (*fout) << "<tr bgcolor=\"#dddddd\">";
      else (*fout) << "<tr>";
      for(size_t j=0;j<tabp->get_ncolumns();j++) {
	(*fout) << "<td>" << ffl.convert(tabp->get(j,i)) << "</td>";
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

int acol_manager::comm_preview(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    if (t3p==0) {
      cerr << "No table3d to preview." << endl;
      return exc_efailed;
    }

    if (scientific) cout.setf(ios::scientific);
    else cout.unsetf(ios::scientific);
    
    if (user_ncols<=0) {
      char *ncstring=getenv("COLUMNS");
      if (ncstring) ncols=o2scl::stoi(ncstring);
    } else {
      ncols=user_ncols;
    }
    
    cout.precision(prec);

    size_t nx, ny;
    t3p->get_size(nx,ny);
    if (nx==0 || ny==0) {
      cout << "No size set. Blank table3d." << endl;
    } else {
      int nrows, ncls;
      if (sv.size()>=2) {
	nrows=o2scl::stoi(sv[1]);
      } else {
	nrows=nx;
      }
      if (((size_t)nrows)>nx) nrows=nx;
      if (sv.size()>=3) {
	ncls=o2scl::stoi(sv[2]);
      } else {
	ncls=ny;
      }
      if (((size_t)ncls)>ny) ncls=ny;
      size_t dx=nx/nrows;
      size_t dy=ny/ncls;
      if (dx==0) dx=1;
      if (dy==0) dy=1;

      size_t nt=t3p->get_nslices(); 
      if (nt!=0) {
	for(size_t k=0;k<nt;k++) {
	  cout << "Slice " << k << ": " << t3p->get_slice_name(k) << endl;
	  
	  cout.setf(ios::showpos);
	  for(size_t i=0;i<((size_t)prec)+8;i++) cout << " ";
	  
	  for(size_t i=0;i<((size_t)ncls);i++) {
	    cout << t3p->get_grid_y(i*dy) << " ";
	  }
	  cout << endl;

	  for(size_t j=0;j<((size_t)nrows);j++) {
	    cout << t3p->get_grid_x(j*dx) << " ";
	    for(size_t i=0;i<((size_t)ncls);i++) {
	      cout << t3p->get(j*dx,i*dy,k) << " ";
	    }
	    cout << endl;
	  }
	  cout.unsetf(ios::showpos);
	}
      }
    }
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to preview." << endl;
    return exc_efailed;
  }

  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);

  if (user_ncols<=0) {
    char *ncstring=getenv("COLUMNS");
    if (ncstring) ncols=o2scl::stoi(ncstring);
  } else {
    ncols=user_ncols;
  }

  cout.precision(prec);

  if (tabp->get_ncolumns()>0) {

    string nlast;
    int inr;
    int inc;

    //-----------------------------------------
    // Compute number of columns which will fit

    size_t max_cols=(ncols)/(8+prec);
    if (max_cols>tabp->get_ncolumns()) max_cols=tabp->get_ncolumns();
    
    //---------------------------------------
    // Compute column and row increment
    
    if (sv.size()==2) {
      int nrows=o2scl::stoi(sv[1]);
      inr=(tabp->get_nlines()+(nrows-1))/(nrows);
      if (inr<1) inr=1;
    } else {
      inr=(tabp->get_nlines()+9)/10;
      if (inr<1) inr=1;
    }
    inc=(tabp->get_ncolumns()+(1))/max_cols;
    if (inc<1) inc=1;
    
    //---------------------------------------
    // Get last row number if necessary
      
    if (pretty==true) {
      nlast=itos(tabp->get_nlines()-1);
    }

    //---------------------------------------
    // Output column names
    
    if (names_out==true) {

      for(size_t ki=0;ki<max_cols;ki++) {

	size_t i=ki*inc;
	if (i>=tabp->get_ncolumns()) i=tabp->get_ncolumns()-1;
	  
	// Preceeding space
	if (pretty==true) {
	  cout << ' ';
	}
	  
	// Column name
	cout << tabp->get_column_name(i) << " ";
	  
	// Trailing spaces
	if (pretty==true) {
	  int nsp=prec+6-((int)(tabp->get_column_name(i).size()));
	  for(int j=0;j<nsp;j++) cout << ' ';
	} else {
	  for(size_t kk=1;kk<nlast.length();kk++) cout << ' ';
	}
	
      }
      cout << endl;
    }
      
    //---------------------------------------
    // Output units
    
    if (names_out==true && tabp->get_nunits()>0) {

      for(size_t ki=0;ki<max_cols;ki++) {

	size_t i=ki*inc;
	if (i>=tabp->get_ncolumns()) i=tabp->get_ncolumns()-1;
	  
	// Preceeding space
	if (pretty==true) {
	  cout << ' ';
	}
	  
	// Column name
	string cunit=tabp->get_unit(tabp->get_column_name(i));
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
      
    //---------------------------------------
    // Output data
      
    for(size_t i=0;i<tabp->get_nlines();i+=inr) {
      
      for(size_t kj=0;kj<max_cols;kj++) {

	size_t j=kj*inc;
	if (j>=tabp->get_ncolumns()) j=tabp->get_ncolumns()-1;

	if (pretty==true) {
	  double d=tabp->get(j,i);
	  if (!has_minus_sign(&d)) {
	    cout << ' ';
	  }
	}
	cout << tabp->get(j,i) << ' ';
	if (pretty==true) {
	  for(int kk=0;kk<((int)(tabp->get_column_name(j).size()-
				 prec-6));kk++) {
	    cout << ' ';
	  }
	}
	
      }
      cout << endl;
      
      //---------------------------------------
      // Continue to next row
    }
    
  }

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

/*
int acol_manager::comm_gen3_list(std::vector<std::string> &sv,
				 bool itive_com) {

  // Delete previous table
  if (threed && t3p!=0) {
    threed=0;
    delete t3p;
    t3p=0;
  } else if (tabp!=0) {
    delete tabp;
    tabp=0;
  }
  
  t3p=new table3d;

  return 0;
}
*/

int acol_manager::comm_generic(std::vector<std::string> &sv, bool itive_com) {

  // Input a generic file
  ifstream ifs;
  ifs.open(sv[1].c_str());
  if (!(ifs)) {
    cerr << "Read failed. Non-existent file?" << endl;
    return exc_efailed;
  }

  // Delete previous table
  if (threed && t3p!=0) {
    threed=0;
    delete t3p;
    t3p=0;
  } else if (tabp!=0) {
    delete tabp;
    tabp=0;
  }

  tabp=new table_units<>(100);
  tabp->read_generic(ifs,verbose);
  
  ifs.close();

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
  cout << "Armadillo support: " << o2scl_settings.armadillo_support() << endl;
  cout << "Eigen support: " << o2scl_settings.eigen_support() << endl;
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
  if (!threed && tabp==0) {
    cerr << "No table to add a constant to." << endl;
    return exc_efailed;
  }
  if (threed && t3p==0) {
    cerr << "No table3d to add a constant to." << endl;
    return exc_efailed;
  }

  if (sv.size()==2) {
    if (verbose>0) {
      cout << "Removing constant named '" << sv[1] << "'." << endl;
    }
    if (threed) {
      t3p->remove_constant(sv[1]);
    } else {
      tabp->remove_constant(sv[1]);
    }
    return 0;
  }

  vector<string> pr, in;
  pr.push_back("Name of constant");
  pr.push_back("Value");
  int ret=get_input(sv,pr,in,"assign",itive_com);
  if (ret!=0) return ret;
  
  if (threed) {
    t3p->add_constant(sv[1],function_to_double(sv[2]));
  } else {
    tabp->add_constant(sv[1],function_to_double(sv[2]));
  }

  return ret;
}

int acol_manager::comm_list(std::vector<std::string> &sv, bool itive_com) {
  if (threed && t3p!=0) {
    cout.precision(prec);
    cout << "Table3d name: " << table_name << endl;
    t3p->summary(&cout,ncols);
  } else if (tabp!=0) {
    cout.precision(prec);
    cout << "Table name: " << table_name << endl;
    if (tabp->get_nunits()>0) {
      tabp->summary(&cout,ncols);
    } else {
      tabp->table<std::vector<double> >::summary(&cout,ncols);
    }
  } else {
    cerr << "No table to list columns for." << endl;
    return exc_efailed;
  }
  return 0;
}

int acol_manager::comm_sort(std::vector<std::string> &sv, bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }
  
  std::string i1;

  if (tabp==0) {
    cerr << "No table to sort." << endl;
    return exc_efailed;
  }
  
  if (sv.size()>1) {
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

  if (tabp->is_column(i1)==false) {
    cerr << "Couldn't find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  if (verbose>1) {
    cout << "Sorting by column " << i1 << endl; 
  }
  tabp->sort_table(i1);
  
  return 0;
}

int acol_manager::comm_stats(std::vector<std::string> &sv, bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }
  
  std::string i1;

  if (tabp==0) {
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

  if (tabp->is_column(i1)==false) {
    cerr << "Couldn't find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  const vector<double> &cref=tabp->get_column(i1);
  cout << "N        : " << tabp->get_nlines() << endl;
  cout << "Sum      : " << vector_mean(tabp->get_nlines(),cref)*
    tabp->get_nlines() << endl;
  cout << "Mean     : " << vector_mean(tabp->get_nlines(),cref) << endl;
  cout << "Std. dev.: " << vector_stddev(tabp->get_nlines(),cref) << endl;
  size_t ix;
  double val;
  vector_min(tabp->get_nlines(),cref,ix,val);
  cout << "Min      : " << val << " at index: " << ix << endl;
  vector_max(tabp->get_nlines(),cref,ix,val);
  cout << "Max      : " << val << " at index: " << ix << endl;

  size_t dup=0, inc=0, dec=0;
  for(size_t i=0;i<tabp->get_nlines()-1;i++) {
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
  if ((dup+inc+dec)!=(tabp->get_nlines()-1)) {
    cout << "Counting mismatch from non-finite values or signed zeros." << endl;
  }
  
  
  return 0;
}

int acol_manager::comm_set(std::vector<std::string> &sv, bool itive_com) {

  // This is taken care of inside cli

  return 0;
}

int acol_manager::comm_select(std::vector<std::string> &sv, bool itive_com) {

  // Remove outermost single or double quotes from all arguments
  for(size_t i=0;i<sv.size();i++) {
    size_t svs=sv[i].size();
    if (svs>2 && ((sv[i][0]=='\'' && sv[i][svs-1]=='\'') ||
		  (sv[i][0]=='\"' && sv[i][svs-1]=='\"'))) {
      sv[i]=sv[i].substr(1,sv[i].size()-2);
    }
  }

  if (threed) {

    if (t3p==0) {
      cerr << "No table3d to select from." << endl;
      return exc_efailed;
    }
    
    size_t nargs=sv.size()-1;
    vector<string> names;
    
    if (nargs==0) {
      cerr << "No columns selected. Command 'select' aborted." << endl;
      return exc_efailed;
    }
    
    for(size_t i=0;i<nargs;i++) {
      names.push_back(sv[i+1]);
    }

    // Make new table3d
    table3d *newt=new table3d;
    
    // Copy grid over
    size_t nx, ny;
    t3p->get_size(nx,ny);
    double *xv=new double[nx];
    double *yv=new double[ny];
    
    for(size_t i=0;i<nx;i++) {
      xv[i]=t3p->get_grid_x(i);
    }
    for(size_t i=0;i<ny;i++) {
      yv[i]=t3p->get_grid_y(i);
    }
    newt->set_xy(t3p->get_x_name(),nx,xv,t3p->get_y_name(),ny,yv);
    delete[] xv;
    delete[] yv;
	
    // Create the new slices
    for(size_t i=0;i<nargs;i++) {
      newt->new_slice(names[i]);
    }

    // Copy the data over
    for(size_t k=0;k<nargs;k++) {

      // Return an error if the slice doesn't exist
      size_t k2;
      if (t3p->is_slice(names[k],k2)==false) {
	cerr << "Slice '" << names[k] << "' is not in the table." << endl;
	return exc_einval;
      }

      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  newt->set(i,j,names[k],t3p->get(i,j,names[k]));
	}
      }
    }
    
    // Delete the old table3d and copy the new one over
    delete t3p;
    t3p=newt;
    
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to select from." << endl;
    return exc_efailed;
  }

  // ----------------------------------------------------------------
  // Parse arguments into names and values
  // ----------------------------------------------------------------

  // The number of column arguments
  int nargs;
  // The name part of the argument
  std::vector<std::string> names;
  // The function part of the argument
  std::vector<std::string> args;
  // Each element is true if the corresponding argument is a pattern
  std::vector<bool> is_pattern;
  
  // ----------------------------------------------------------------

  nargs=sv.size()-1;
  
  if (nargs==0) {
    cerr << "No columns selected. Command 'select' aborted." << endl;
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

  // ----------------------------------------
  // Create new table
  // ----------------------------------------
  
  table_units<> *new_table;
  new_table=new table_units<>(tabp->get_nlines());
  
  new_table->set_nlines(tabp->get_nlines());

  // ----------------------------------------
  // Copy constants from old to new table
  // ----------------------------------------

  for(size_t i=0;i<tabp->get_nconsts();i++) {
    string tnam;
    double tval;
    tabp->get_constant(i,tnam,tval);
    new_table->add_constant(tnam,tval);
  }
  
  // ----------------------------------------
  // Add columns and data to new table
  // ----------------------------------------

  bool *matched=new bool[tabp->get_ncolumns()];
  for(size_t i=0;i<tabp->get_ncolumns();i++) {
    matched[i]=false;
  }

  // In this next loop, we need fix the code to ensure that when a
  // non-pattern column is given, it sets matched to true
  // appropriately.

  for(int i=0;i<nargs;i++) {

    if (is_pattern[i]==false) {
      
      // Return an error if the column doesn't exist
      if (names[i]==args[i] && tabp->is_column(args[i])==false) {
	cerr << "Column '" << args[i] << "' is not in the table." << endl;
	return exc_einval;
      }

      // Add the new column to the new table
      new_table->new_column(names[i]);

      // If necessary, set units
      if (names[i]==args[i] && tabp->get_unit(args[i]).length()!=0) {
	new_table->set_unit(names[i],tabp->get_unit(args[i]));
      }

      // Fill column with the new data
      ubvector vec(tabp->get_nlines());
      tabp->function_vector(args[i],vec,false);
      new_table->copy_to_column(vec,names[i]);

    } else {

      // Find the matching columns
      for(size_t j=0;j<tabp->get_ncolumns();j++) {

	if (matched[j]==false &&  
	    fnmatch(args[i].c_str(),
		    tabp->get_column_name(j).c_str(),0)==0) {

	  // If we've found a match, add it to the new table
	  matched[j]=true;

	  // Add the new column to the new table
	  new_table->new_column(tabp->get_column_name(j));

	  // If necessary, set units
	  string tmp=tabp->get_column_name(j);
	  if (tabp!=0 && tabp->get_unit(tmp).length()!=0) {
	    new_table->set_unit(tmp,tabp->get_unit(tmp));
	  }

	  // Fill it with the new data
	  ubvector vec(tabp->get_nlines());
	  tabp->function_vector(tabp->get_column_name(j),vec,false);
	  new_table->copy_to_column(vec,tabp->get_column_name(j));
	}
      }
    }
  }

  delete[] matched;

  // Replace the old table with the new one
  delete tabp;
  tabp=new_table;

  // Call list command
  if (verbose>0) {
    comm_list(sv,itive_com);
  }

  return 0;
}

int acol_manager::comm_delete_rows(std::vector<std::string> &sv, 
				   bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
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
  
  tabp->delete_rows(i1);

  return 0;
}

int acol_manager::comm_select_rows(std::vector<std::string> &sv, 
				   bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
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
    
  // ----------------------------------------
  // Create new table
  // ----------------------------------------
  
  table_units<> *new_table=new table_units<>;
  
  // ----------------------------------------
  // Copy constants from old to new table
  // ----------------------------------------

  for(size_t i=0;i<tabp->get_nconsts();i++) {
    string tnam;
    double tval;
    tabp->get_constant(i,tnam,tval);
    new_table->add_constant(tnam,tval);
  }
  
  // ----------------------------------------
  // Add column names to new table
  // ----------------------------------------

  for(int i=0;i<((int)tabp->get_ncolumns());i++) {
    new_table->new_column(tabp->get_column_name(i));
  }

  // ----------------------------------------
  // Copy data from selected rows
  // ----------------------------------------

  int new_lines=0;
  for(int i=0;i<((int)tabp->get_nlines());i++) {
    if (tabp->row_function(i1,i)>0.5) {
      new_table->set_nlines(new_lines+1);
      for(int j=0;j<((int)tabp->get_ncolumns());j++) {
	new_table->set(j,new_lines,tabp->get(j,i));
      }
      new_lines++;
    }
  }
  
  // Replace the old table with the new one
  delete tabp;
  tabp=new_table;

  return 0;
}

int acol_manager::comm_delete_col(std::vector<std::string> &sv, 
				  bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
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
    
  if (tabp->is_column(i1)==false) {
    cerr << "Couldn't find column named '" << i1 << "'." << endl;
    return exc_efailed;
  }

  tabp->delete_column(i1);

  return 0;
}

int acol_manager::comm_create(std::vector<std::string> &sv, bool itive_com) {
  std::string i1, i2, i3, i4;
  
  if (sv.size()>=5) {
    i1=sv[1];
    i2=sv[2];
    i3=sv[3];
    i4=sv[4];
  } else if (itive_com) {
    i1=cl->cli_gets("Name of new column (or blank to quit): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'create' cancelled." << endl;
      return exc_efailed;
    }
    i2=cl->cli_gets("Value for first row (or blank to quit): ");
    if (i2.length()==0) {
      if (verbose>0) cout << "Command 'create' cancelled." << endl;
      return exc_efailed;
    }
    i3=cl->cli_gets("Maximum value (or blank to quit): ");
    if (i3.length()==0) {
      if (verbose>0) cout << "Command 'create' cancelled." << endl;
      return exc_efailed;
    }
    i4=cl->cli_gets("Increment (or blank to quit): ");
    if (i4.length()==0) {
      if (verbose>0) cout << "Command 'create' cancelled." << endl;
      return exc_efailed;
    }
  } else {
    cerr << "No arguments specified to 'create'." << endl;
    return exc_efailed;
  }
  
  double d2=function_to_double(i2);
  double d3=function_to_double(i3);
  double d4=function_to_double(i4);
  d3+=d4/1.0e4;
  int cnl=((int)((d3-d2)/d4))+1;
  tabp=new table_units<>(cnl);
  
  tabp->line_of_names(i1);

  for(int li=0;li<cnl;li++) {
    tabp->set(i1,li,o2scl::stod(i2)+((double)li)*o2scl::stod(i4));
  }

  return 0;
}

int acol_manager::comm_insert(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {

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
    hf.open(in[0]);
    std::string tmp_name;
    if (in[1].length()>0) tmp_name=in[1];
    hdf_input(hf,tmp,tmp_name);
    hf.close();

    t3p->add_slice_from_table(tmp,in[2],in[3]);

    return 0;
  }

  if (tabp==0) {
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
  hf.open(in[0]);
  std::string tmp_name;
  if (in[1].length()>0) tmp_name=in[1];
  hdf_input(hf,tmp,tmp_name);
  hf.close();

  tabp->add_col_from_table(tmp,in[2],in[3],in[4],in[5]);

  return 0;

}

int acol_manager::comm_insert_full(std::vector<std::string> &sv, 
				   bool itive_com) {
  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
    cerr << "No table to insert columns into." << endl;
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
      cerr << "Not enough arguments to 'insert'" << endl;
      return exc_efailed;
    }
  }

  cout << "Unimplemented." << endl;

  return 1;

#ifdef O2SCL_NEVER_DEFINED

  for (size_t j=0;j<tmp->get_ncolumns();j++) {
    string ty=tmp->get_column_name(j);
    int tret=tabp->add_col_from_table(in[1],*tmp,in[2],ty,ty);
    if (tret!=0) {
      cerr << "Adding column " << ty << " failed." << endl;
      ret=tret;
      // We don't return here so that "delete tmp;" is called below
    }
  }
  
#endif
}

int acol_manager::comm_interp(std::vector<std::string> &sv, bool itive_com) {
  if (threed) {

    // --------------------------------------------------------------
    // 3d table interpolation
    
    if (t3p==0) {
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
	cerr << "Not enough arguments to 'interp'" << endl;
	return exc_efailed;
      }
    }
  
    double ret=t3p->interp(function_to_double(in[1]),
			   function_to_double(in[2]),in[0]);
    if (err_hnd->get_errno()!=0) {
      cerr << "Interpolation failed." << endl;
      return exc_efailed;
    } else {
      cout << "Interpolation result: " << ret << endl;
    }

    return 0;
  }

  // --------------------------------------------------------------
  // 2d table interpolation

  if (tabp==0) {
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
  
  if (tabp->is_column(in[0])==false) {
    cerr << "Couldn't find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }
  if (tabp->is_column(in[2])==false) {
    cerr << "Couldn't find column named '" << in[2] << "'." << endl;
    return exc_efailed;
  }

  double ret=tabp->interp(in[0],function_to_double(in[1]),in[2]);
  if (err_hnd->get_errno()!=0) {
    cerr << "Interpolation failed." << endl;
    return exc_efailed;
  } else {
    cout << "Interpolation result: " << ret << endl;
  }
  
  return 0;
}

int acol_manager::separate(std::string str, std::vector<std::string> &sv) {
  string tmp, tmp2;

  istringstream *is=new istringstream(str.c_str());
  while ((*is) >> tmp) {
    
    // If it begins with a quote, add more words accordingly
    if (tmp[0]=='\"') {
      tmp2=tmp.substr(1,tmp.size()-1);
      bool done=false;
      while (done==false) {
	if (!((*is) >> tmp)) {
	  done=true;
	} else if (tmp[tmp.size()-1]=='\"') {
	  tmp=tmp.substr(0,tmp.size()-1);
	  done=true;
	}
	tmp2+=" ";
	tmp2+=tmp;
      }
      tmp=tmp2;
    }

    // Add to the list
    sv.push_back(tmp);
  }
  delete is;

  return 0;
}

int acol_manager::comm_fit(std::vector<std::string> &sv, bool itive_com) {

  if (threed) {
    cout << "Not implemented for table3d." << endl;
    return 0;
  }

  if (tabp==0) {
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
  size_t ndat=tabp->get_nlines();
  ubvector xdat(ndat), ydat(ndat), yerr(ndat);
  for(size_t i=0;i<ndat;i++) {
    xdat[i]=tabp->get(in[0],i);
    ydat[i]=tabp->get(in[1],i);
    yerr[i]=tabp->get(in[2],i);
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
  separate(in[6],param_list);
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
  fit_funct11_strings<> ffs(in[5],in[4],in[0],0);
  ffs.set_aux_parms(params);

  // Fitting function object
  chi_fit_funct<ubvector,ubmatrix,fit_funct11_strings<> > 
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
  tabp->new_column(in[3]);
  for(int k=0;k<((int)tabp->get_nlines());k++) {
    tabp->set(in[3],k,ffs(n_parms,params,(xdat)[k]));
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
  
  return 0;
}
