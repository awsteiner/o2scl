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
#ifndef O2SCL_ACOLM_H
#define O2SCL_ACOLM_H

/** \file acolm.h
    \brief The \ref o2scl_acol::acol_manager class header
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <fnmatch.h>
#include <o2scl/misc.h>
#include <o2scl/cli.h>
#include <o2scl/fit_nonlin.h>
#include <o2scl/table_units.h>
#include <o2scl/table3d.h>
#include <o2scl/format_float.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/contour.h>
#include <o2scl/tensor_grid.h>

#ifdef O2SCL_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

/// A namespace for objects associated with the command-line utility 'acol'
namespace o2scl_acol {

  /** \brief The driver for 'acol' command-line utility
      \nothing

      \todo Find a way to reorganize the source code into
      smaller files. Separate out the global functions for
      o2graph into their own file.

      \comment
      There was some concern about the confusion b/w the commands
      "get row by index" (get-row) and "get row by function"
      (find-row) and "select rows by function" (select-rows),
      but it might be ok.
      \endcomment

      \future Stack-like operations (push, pop, swap, 
      stack-list, etc.)?

      \future Use get_input() in comm_create?

      \future Add functionality to ensure that three digit exponents
      are still handled gracefully (do this by creating a new boolean
      setting which, if true, always makes three spaces for
      exponents?)

      \future Fix insert and insert_full so that it automatically
      renames columns
      
      \future Allow "insert" commands to be restrictive, avoiding
      extrapolation
      
      \future Replace ~ with $HOME in filenames (this might be best
      done inside the \ref o2scl::cli class). (Some progress made.
      Function cli::expand_tilde() is written but not yet
      implemented.)
      
      \future New 3d commands: transpose (for table3d), contours, 
      find_grid_x, find_grid_y.

      \hline
  */
  class acol_manager {

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief If true, then run in o2graph mode
     */
    bool o2graph_mode;
    
    /** \brief The object for the set function
     */
    o2scl::comm_option_mfptr<acol_manager> cset;
    
    /** \brief Add new commands for type \c new_type
     */
    void command_add(std::string new_type);
    
#ifdef DOXYGEN
    /// The number formatter for html output
    format_float ffl;
#else
    o2scl::format_float ffl;
#endif

    /// Convert units object (initialized by constructor to global object)
    o2scl::convert_units &cng;

    /// \name Parameters modifiable by the user
    //@{
    /// The output precision (default 6)
    int prec;

    /// The verbosity level (default 1)
    int verbose;

    /// True if we should make the output into neat columns (default true)
    bool pretty;
    
    /// True if we should output column names
    bool names_out;

    /// The name of the table
    std::string obj_name;
  
    /// Filename for units command
    std::string unit_fname;

    /// Default arguments from environment
    std::string def_args;

    /// The number of columns requested by the user
    int user_ncols;

    /// The interpolation type
    int interp_type;

    /// If set, try to compress
    int compress;
    
    /// True for scientific output mode
    bool scientific;
    //@}

    /// \name The parameter objects
    //@{
    o2scl::cli::parameter_string p_obj_name;
    o2scl::cli::parameter_string p_unit_fname;
    o2scl::cli::parameter_string p_def_args;
    o2scl::cli::parameter_int p_verbose;
    o2scl::cli::parameter_int p_compress;
    o2scl::cli::parameter_int p_prec;
    o2scl::cli::parameter_int p_ncols;
    o2scl::cli::parameter_int p_interp_type;
    o2scl::cli::parameter_bool p_scientific;
    o2scl::cli::parameter_bool p_pretty;
    o2scl::cli::parameter_bool p_names_out;
    //@}

    /// Number of columns in screen
    int ncols;

#endif

  public:

    acol_manager();

    virtual ~acol_manager() {}

    /// String designating the current type
    std::string type;

    /// Dummy cli object for cli::cli_gets()
#ifdef DOXYGEN
    cli *cl;
#else
    o2scl::cli *cl;
#endif

    /// \name Object storage
    //@{
    o2scl::table_units<> table_obj;
    o2scl::table3d table3d_obj;
    o2scl::hist hist_obj;
    o2scl::hist_2d hist_2d_obj;

    int int_obj;
    char char_obj;
    double double_obj;
    size_t size_t_obj;
    std::string string_obj;

    std::vector<o2scl::contour_line> cont_obj;
    o2scl::uniform_grid<double> ug_obj;
    
    std::vector<int> intv_obj;
    std::vector<double> doublev_obj;
    std::vector<size_t> size_tv_obj;
    std::vector<std::string> stringv_obj;

    o2scl::tensor<> tensor_obj;
    o2scl::tensor<int> tensor_int_obj;
    o2scl::tensor<size_t> tensor_size_t_obj;
    o2scl::tensor_grid<> tensor_grid_obj;

    o2scl::prob_dens_mdim_amr<> pdma_obj;
    //@}
    
    /** \brief True if we should run interactive mode after parsing
	the command-line
    */
    bool post_interactive;

    /// The environment variable to read from 
    std::string env_var_name;

    /// String parameters
    std::map<std::string,std::string *> str_params;

    /// Integer parameters
    std::map<std::string,int *> int_params;

  protected:

    /** \brief Clear memory associated with the current object and set
	type to ""
    */
    void clear_obj();

    /** \brief Remove the type-specific commands
     */
    void command_del();
    
    // Ensure \c col is unique from entries in \c cnames
    //int make_unique_name(std::string &col, std::vector<std::string> &cnames);

  public:
    
    /** \brief Main run function

	Process command-line options using cli object, interface with
	the operating system via getenv(), instantiate and call the
	acol_manager object.
    */
    virtual int run(int argv, char *argc[], bool full_process=true);

    /// Create the cli object (with readline support if available)
    virtual int setup_cli();

    /// Add the options to the cli object
    virtual int setup_options();

    /// Add the help text to the cli object
    virtual int setup_help();

    /// Add the parameters for 'set' to the cli object
    virtual int setup_parameters();

    /// Assign a constant
    virtual int comm_assign(std::vector<std::string> &sv, bool itive_com);

    /// Compute a scalar value
    virtual int comm_calc(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output the help text
     */
    virtual int comm_help(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief List commands, with an optional type argument
     */
    virtual int comm_commands(std::vector<std::string> &sv, bool itive_com);
    
    /// Create a table from a column of equally spaced values
    virtual int comm_create(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the grid for a \ref o2scl::tensor_grid object
     */
    virtual int comm_set_grid(std::vector<std::string> &sv, bool itive_com);

    /// Download a file from a specified URL
    virtual int comm_download(std::vector<std::string> &sv, bool itive_com);

    /// Delete a column
    virtual int comm_delete_col(std::vector<std::string> &sv, bool itive_com);

    /// Delete rows specified by a function
    virtual int comm_delete_rows(std::vector<std::string> &sv, bool itive_com);
    
    /// Delete rows which match to within a specified tolerance
    virtual int comm_delete_rows_tol(std::vector<std::string> &sv,
				     bool itive_com);

    /// Create a column which is the derivative of another
    virtual int comm_deriv(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table object
     */
    virtual int comm_to_table(std::vector<std::string> &sv, bool itive_com);

    /** \brief For tensor object, get entries along the main diagonal
     */
    virtual int comm_diag(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object
     */
    virtual int comm_to_table3d(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the autocorrelation coefficient
     */
    virtual int comm_autocorr(std::vector<std::string> &sv, bool itive_com);

    /// Create a slice which is the derivative wrt x of another
    virtual int comm_deriv_x(std::vector<std::string> &sv, bool itive_com);

    /// Create a slice which is the derivative wrt y of another
    virtual int comm_deriv_y(std::vector<std::string> &sv, bool itive_com);

    /// Create a column which is the second derivative of another
    virtual int comm_deriv2(std::vector<std::string> &sv, bool itive_com);

    /// Parameters for iterate_func()
    typedef struct {
      std::string tname;
      o2scl_hdf::hdf_file *hf;
      bool found;
      std::string type;
      int verbose;
      int mode;
    } iter_parms;

    /// HDF object iteration function
    static herr_t iterate_func(hid_t loc, const char *name, 
			       const H5L_info_t *inf, void *op_data);
    
    /// Read a file and list the O2scl objects
    virtual int comm_filelist(std::vector<std::string> &sv, bool itive_com);

    /// Read an object from a file
    virtual int comm_read(std::vector<std::string> &sv, bool itive_com);

    /// Add 'nlines' as a constant to a \ref o2scl::table object
    virtual int comm_nlines(std::vector<std::string> &sv, bool itive_com);

    /// Convert a \ref o2scl::table object to a \ref o2scl::hist object
    virtual int comm_to_hist(std::vector<std::string> &sv, bool itive_com);

    /// Convert a \ref o2scl::table object to a \ref o2scl::hist object
    virtual int comm_to_hist_2d(std::vector<std::string> &sv, bool itive_com);

    /// Output the type of the current object to the screen
    virtual int comm_type(std::vector<std::string> &sv, bool itive_com);
    
    /// Find a row from a function
    virtual int comm_find_row(std::vector<std::string> &sv, bool itive_com);
    
    /// Create a column from a function
    virtual int comm_function(std::vector<std::string> &sv, bool itive_com);

    /// Read a generic data file
    virtual int comm_generic(std::vector<std::string> &sv, bool itive_com);

    /// Print out an entire row
    virtual int comm_get_row(std::vector<std::string> &sv, bool itive_com);

    /** \brief Extract a slice from a table3d object to generate a 
	\ref o2scl::table object
    */
    virtual int comm_slice(std::vector<std::string> &sv, bool itive_com);

    /// Fit two columns to a function
    virtual int comm_fit(std::vector<std::string> &sv, bool itive_com);

    /// Insert a column from an external table using interpolation
    virtual int comm_insert(std::vector<std::string> &sv, bool itive_com);

    /// Insert an external table using interpolation
    virtual int comm_insert_full(std::vector<std::string> &sv, bool itive_com);

    /// Create a column which is the integral of another
    virtual int comm_integ(std::vector<std::string> &sv, bool itive_com);

    /// Toggle interactive mode
    virtual int comm_interactive(std::vector<std::string> &sv, bool itive_com);

    /// Output to a file in internal format
    virtual int comm_internal(std::vector<std::string> &sv, bool itive_com);

    /// Create an html file
    virtual int comm_interp(std::vector<std::string> &sv, bool itive_com);

    /// List columns in table 'tp' named 'tname' assuming screen size 'ncol'
    virtual int comm_list(std::vector<std::string> &sv, bool itive_com);

    /// Compute the maximum value of a colum
    virtual int comm_max(std::vector<std::string> &sv, bool itive_com);
    
    /// Compute the minimum value of a colum
    virtual int comm_min(std::vector<std::string> &sv, bool itive_com);
    
    /// Add a column for line numbers
    virtual int comm_index(std::vector<std::string> &sv, bool itive_com);

    /// Output to screen or file
    virtual int comm_output(std::vector<std::string> &sv, bool itive_com);

    /// Preview the table
    virtual int comm_preview(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get or set the value 
     */
    virtual int comm_value(std::vector<std::string> &sv, bool itive_com);

    /// Concatenate two table/table3d objects
    virtual int comm_cat(std::vector<std::string> &sv, bool itive_com);

    /// Sum two table/table3d objects
    virtual int comm_sum(std::vector<std::string> &sv, bool itive_com);

    /// Rename a column
    virtual int comm_rename(std::vector<std::string> &sv, bool itive_com);

    /// Select several columns for a new table
    virtual int comm_select(std::vector<std::string> &sv, bool itive_com);

    /// Select several rows for a new table
    virtual int comm_select_rows(std::vector<std::string> &sv,
				 bool itive_com);

    /** \brief A faster form of select rows which requires one to specify
	the columns needed for the selection criteria first
    */
    virtual int comm_select_rows2(std::vector<std::string> &sv,
				  bool itive_com);

    /// Post-processing for setting a value
    virtual int comm_set(std::vector<std::string> &sv, bool itive_com);

    /// Set an individual data point at a specified row and column
    virtual int comm_set_data(std::vector<std::string> &sv, bool itive_com);
    
    /// Set units of a column
    virtual int comm_set_unit(std::vector<std::string> &sv, bool itive_com);
    
    /// Compute contour lines
    virtual int comm_contours(std::vector<std::string> &sv, bool itive_com);
    
    /// Set units of a column
    virtual int comm_show_units(std::vector<std::string> &sv, bool itive_com);
    
    /// Get units of a column
    virtual int comm_get_unit(std::vector<std::string> &sv, bool itive_com);
    
    /// Get an entry
    virtual int comm_entry(std::vector<std::string> &sv, bool itive_com);
    
    /// Convert units of a column
    virtual int comm_convert_unit(std::vector<std::string> &sv, 
				  bool itive_com);
    
    /// Sort the table by a column
    virtual int comm_sort(std::vector<std::string> &sv, bool itive_com);

    /// Get column stats
    virtual int comm_stats(std::vector<std::string> &sv, bool itive_com);

    /// Get column stats with weights specified in a second column
    virtual int comm_wstats(std::vector<std::string> &sv, bool itive_com);

    /// Print version
    virtual int comm_version(std::vector<std::string> &sv, bool itive_com);

    /// Get a conversion factor
    virtual int comm_get_conv(std::vector<std::string> &sv, bool itive_com);

    /// Set screen width
    int set_swidth(int ncol) {
      ncols=ncol;
      return 0;
    }

  protected:
    
    /// An internal command for prompting the user for command arguments
    int get_input(std::vector<std::string> &sv, 
		  std::vector<std::string> &directions,
		  std::vector<std::string> &in, std::string comm_name,
		  bool itive_com);

    /// An internal command for prompting the user for one command argument
    int get_input_one(std::vector<std::string> &sv, std::string directions,
		      std::string &in, std::string comm_name,
		      bool itive_com);
    
  public:
    
    /// \name Temporary storage for \ref o2scl_acol_get_slice()
    //@{
    std::vector<double> xtemp;
    std::vector<double> ytemp;
    std::vector<double> stemp;
    //@}
    
  };
  
}

extern "C" {

  // (remember that \ref's don't work in \name groups)
  /// \name Functions to integrate o2scl_acol::acol_manager with python
  //@{
  /** \brief Create an \ref o2scl_acol::acol_manager object
      
      This function creates an object of type
      \ref o2scl_acol::acol_manager with the <tt>new</tt>
      operator and then calls the function
      \ref o2scl_acol::acol_manager::run() .
   */
  void *o2scl_create_acol_manager();
  
  /** \brief Free memory associated with a \ref
      o2scl_acol::acol_manager object

      This function uses <tt>delete</tt> to free the
      memory associated with an object of type
      \ref o2scl_acol::acol_manager .
  */
  void o2scl_free_acol_manager(void *vp);

  /** \brief Set the command name, the short description,
      and the environment variable name
   */
  void o2scl_acol_set_names(void *vp, int n1, char *cmd_name,
			    int n2, char *short_desc, int n3,
			    char *env_var);
  
  /** \brief Construct a string vector from the data in 
      \c n_entries, \c sizes, and \c str
  */
  std::vector<std::string> o2scl_acol_parse_arrays
  (int n_entries, int *sizes, char *str);
  
  /** \brief Parse the set of commands in \c n_entries, \c sizes
      and \c str
   */
  void o2scl_acol_parse(void *vp, int n_entries, int *sizes, 
			char *str);

  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object
   */
  int o2scl_acol_get_column(void *vp, char *col_name,
			    int &n, double *&ptr);

  /** \brief Return the size and a pointer to the row
      with index \c row_index in a \ref o2scl::table object
   */
  int o2scl_acol_get_row_ser(void *vp, char *parttern, int row_index,
			     int &n, double *&ptr);
  
  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object
   */
  int o2scl_acol_get_double_arr(void *vp, int &n, double *&ptr);
  
  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object
   */
  int o2scl_acol_get_hist_reps(void *vp, int &n, double *&ptr);

  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object
   */
  int o2scl_acol_get_hist_wgts(void *vp, int &n, double *&ptr);

  /** \brief Return the dimensionality, mesh size, and 
      lower and upper limits for a \ref o2scl::prob_dens_mdim_amr 
      object.
   */
  int o2scl_acol_pdma_get_base(void *vp, int &ndim, int &n, 
			       double *&low, double *&high);

  /** \brief Return the lower and upper limits, fractional volume, and
      weight for the \ref o2scl::prob_dens_mdim_amr::hypercube object
      of index \c ix
   */
  int o2scl_acol_pdma_get_cube(void *vp, int ix, 
			       double *&low, double *&high,
			       double &frac_vol, double &weight);

  /** \brief Return the number of contour lines associated with
      the current contour line vector object
   */
  int o2scl_acol_contours_n(void *vp);
  
  /** \brief For the current contour line vector object, set the
      pointers to the x- and y-values in the contour lines and return
      the contour level
  */
  double o2scl_acol_contours_line(void *vp, int i, int &n, double *&ptrx,
				  double *&ptry);

  /** \brief Return the type of the current object 
   */
  void o2scl_acol_get_type(void *vp, int &n, char *&str);

  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object
   */
  int o2scl_acol_get_slice(void *vp, char *slice_name,
			   int &nx, double *&xptr,
			   int &ny, double *&yptr,
			   double *&data);
  
  /** \brief For a two-dimensional histogram, return the bin edges,
      number of bins in both directions, and the weights in each bin
   */
  int o2scl_acol_get_hist_2d(void *vp, 
			     int &nx, double *&xptr,
			     int &ny, double *&yptr,
			     double *&data);
  //@}
  
}

#endif
