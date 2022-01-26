/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/uniform_grid.h>
#include <o2scl/slack_messenger.h>
#include <o2scl/rng.h>

#ifdef O2SCL_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

/// A namespace for objects associated with the command-line utility 'acol'
namespace o2scl_acol {

  /** \brief The driver for 'acol' command-line utility
      \nothing

      \verbatim embed:rst

      .. todo::

         In class acol_manager:

         - (Future) There is quite a bit of code duplication in
           comm_autocorr() between the "table" and "other" types. 
           This could be streamlined.

         - (Future) sum/max/min/output/interp/deriv/integ/deriv2 
           for hist, hist_2d, and v<c>

         - (Future) Commands xindex and yindex for table3d.

         - (Future) Enable set_grid() for table3d similar to tensor_grid.

         - (Future) Fix fit for table.

         - (Future) Use swap instead of copy in 'select' for table objects.

         - (Future) Make sure get_input() is used more consistently.

         - (Future) Make sure preview, output, internal, generic, and create
           work consistently across all types.

         - (Future) Stack-like operations (push, pop, swap, 
           stack-list, etc.)?

         - (Future) Add functionality to ensure that three digit exponents
           are still handled gracefully (do this by creating a new boolean
           setting which, if true, always makes three spaces for
           exponents?)

         - (Future) Fix insert and insert_full so that it automatically
           renames columns
      
         - (Future) Allow "insert" commands to be restrictive, avoiding
           extrapolation

      \endverbatim
  */
  class acol_manager {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Random number generator
    o2scl::rng<> rng;
    
    /// The object which sends Slack messages
    o2scl::slack_messenger smess;
    
    /** \brief A list of all type-specific commands for each type
     */
    std::map<std::string,std::vector<std::string> > type_comm_list;

    /** \brief A list of all types
     */
    std::vector<std::string> type_list;
    
    /** \brief If true, then run in o2graph mode
     */
    bool o2graph_mode;
    
    /** \brief The object for the set function
     */
    o2scl::comm_option_mfptr<acol_manager> cset;
    
#ifdef DOXYGEN
    /// The number formatter for html output
    format_float ffl;
#else
    o2scl::format_float ffl;
#endif

    /// Convert units object (initialized by constructor to global object)
    o2scl::convert_units<double> &cng;

    /// \name Parameters modifiable by the user
    //@{
    /// The output precision (default 6)
    int prec;

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

    /// \name Other data [protected]
    //@{
    /// Number of columns in screen
    int ncols;
    //@}

#endif

  public:

    acol_manager();

    virtual ~acol_manager() {}

    /// The verbosity level (default 1)
    int verbose;

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

    /// Integer parameters
    std::map<std::string,int *> int_params;

  protected:

    /** \brief Clear memory associated with the current object and set
        type to ""
    */
    void clear_obj();

  public:

    /** \brief Add new commands for type \c new_type
     */
    void command_add(std::string new_type);
    
    /** \brief Remove the type-specific commands

        \note This needs to be public for the o2graph interface
     */
    void command_del(std::string ltype);
    
    /** \brief Get the verbose parameter
        
        This function is used in \ref o2scl_acol_mult_vectors_to_conts() .
     */
    int get_verbose() {
      return verbose;
    }
    
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

    /// \name Functions for the interface
    //@{
    /// Assign a constant
    virtual int comm_assign(std::vector<std::string> &sv, bool itive_com);

    /// Convert a series of histograms to a table3d object
    virtual int comm_ser_hist_t3d(std::vector<std::string> &sv,
                                  bool itive_com);

    /// Binary function for tensors
    virtual int comm_binary(std::vector<std::string> &sv, bool itive_com);

    /// Average rows together in a table
    virtual int comm_average_rows(std::vector<std::string> &sv,
                                  bool itive_com);

    /// Compute correlation between table columns
    virtual int comm_correl(std::vector<std::string> &sv, bool itive_com);

    /// Refine an object
    virtual int comm_refine(std::vector<std::string> &sv, bool itive_com);

    /// Compute a scalar value
    virtual int comm_calc(std::vector<std::string> &sv, bool itive_com);

    /// Clear the current object
    virtual int comm_clear(std::vector<std::string> &sv, bool itive_com);

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

    /** \brief Get the grid for a \ref o2scl::tensor_grid object
     */
    virtual int comm_get_grid(std::vector<std::string> &sv, bool itive_com);

    /// Download a file from a specified URL
    virtual int comm_download(std::vector<std::string> &sv, bool itive_com);

    /// Open the local HTML documentation
    virtual int comm_docs(std::vector<std::string> &sv, bool itive_com);

    /// Open the HTML documentation
    virtual int comm_wdocs(std::vector<std::string> &sv, bool itive_com);

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

    /** \brief Convert object to a \ref o2scl::tensor_grid object
     */
    virtual int comm_to_tensor_grid(std::vector<std::string> &sv,
                                    bool itive_com);
    
    /** \brief Convert object to a \ref o2scl::tensor object
     */
    virtual int comm_to_tensor(std::vector<std::string> &sv,
                               bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object
        by summing over tensor indices
     */
    virtual int comm_to_table3d_sum(std::vector<std::string> &sv,
                                    bool itive_com);

    /** \brief Compute the autocorrelation coefficients
     */
    virtual int comm_autocorr(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the autocorrelation coefficient
     */
    virtual int comm_ac_len(std::vector<std::string> &sv,
                            bool itive_com);

    /// Create a slice which is the derivative wrt x of another
    virtual int comm_deriv_x(std::vector<std::string> &sv, bool itive_com);

    /// Create a slice which is the derivative wrt y of another
    virtual int comm_deriv_y(std::vector<std::string> &sv, bool itive_com);

    /// Create a column which is the second derivative of another
    virtual int comm_deriv2(std::vector<std::string> &sv, bool itive_com);

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

    /// Add a column from a vector_specification
    virtual int comm_add_vec(std::vector<std::string> &sv, bool itive_com);

    /// Read a generic data file
    virtual int comm_generic(std::vector<std::string> &sv, bool itive_com);

    /// Print out an entire row
    virtual int comm_get_row(std::vector<std::string> &sv, bool itive_com);

    /** \brief Extract a slice from a table3d object to generate a 
        \ref o2scl::table object
    */
    virtual int comm_slice(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert a slice to a histogram
     */
    virtual int comm_slice_hist(std::vector<std::string> &sv, bool itive_com);

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

    /// Set the name of the x grid
    virtual int comm_x_name(std::vector<std::string> &sv, bool itive_com);

    /// Set the name of the y grid
    virtual int comm_y_name(std::vector<std::string> &sv, bool itive_com);
    
    /// Add a column for line numbers
    virtual int comm_index(std::vector<std::string> &sv, bool itive_com);

    /// Output to screen or file
    virtual int comm_output(std::vector<std::string> &sv, bool itive_com);

    /// Rearrange a tensor
    virtual int comm_rearrange(std::vector<std::string> &sv, bool itive_com);

    /// Preview the table
    virtual int comm_preview(std::vector<std::string> &sv, bool itive_com);

    /// Send a slack message
    virtual int comm_slack(std::vector<std::string> &sv, bool itive_com);

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

    /** \brief A faster form of select rows which requires one to specify
        the columns needed for the selection criteria first
    */
    virtual int comm_select_rows(std::vector<std::string> &sv,
                                  bool itive_com);

    /// Post-processing for setting a value
    virtual int comm_set(std::vector<std::string> &sv, bool itive_com);

    /// Set an individual data point at a specified row and column
    virtual int comm_set_data(std::vector<std::string> &sv, bool itive_com);
    
    /// Set units of a column
    virtual int comm_set_unit(std::vector<std::string> &sv, bool itive_com);
    
    /// Compute contour lines
    virtual int comm_contours(std::vector<std::string> &sv, bool itive_com);
    
    /// Get units of a column
    virtual int comm_get_unit(std::vector<std::string> &sv, bool itive_com);
    
    /// Get an entry
    virtual int comm_entry(std::vector<std::string> &sv, bool itive_com);

    /// Get an entry by grid point
    virtual int comm_entry_grid(std::vector<std::string> &sv, bool itive_com);
    
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
    virtual int comm_convert(std::vector<std::string> &sv, bool itive_com);

    /// Copy an HDF5 file
    virtual int comm_h5_copy(std::vector<std::string> &sv, 
                             bool itive_com);

    /// Search for, or add or delete a constant
    virtual int comm_constant(std::vector<std::string> &sv, bool itive_com);
    //@}
    
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
    
    /** \name Temporary storage 

        These are used in \ref o2scl_acol_get_slice(), \ref
        o2scl_acol_get_hist_reps(), \ref o2scl_acol_get_hist_wgts(),
        \ref o2scl_acol_get_hist_bins(), and \ref
        o2scl_acol_get_hist_2d(),
    */
    //@{
    std::vector<double> xtemp;
    std::vector<double> ytemp;
    std::vector<double> stemp;
    std::vector<int> itemp;
    std::vector<char> ctemp;
    //@}

    // End of class acol_manager
  };

  /** \brief Construct a string vector from the data in 
      \c n_entries, \c sizes, and \c str

      This function operates on an integer \c n_entries, an array \c
      sizes (which has length \c n_entries) and an array of characters
      \c str which has a length equal to the sum of the entries in the
      array \c sizes. The \c sizes array contains the length of each
      string, and the \c str array contains the characters in multiple
      strings, concatenated together to form a single combined string.
      This function takes the data in these three objects and creates
      an object of type <tt>vector&lt;string&gt;</tt> from it, similar
      to the way that \ref o2scl_hdf::hdf_file::gets_vec() reads a string
      array from an HDF5 file.

      This function is used in \ref o2scl_acol_parse(), \ref
      o2scl_acol_alias_counts() and \ref o2scl_acol_apply_aliases() .
  */
  std::vector<std::string> parse_arrays
    (int n_entries, int *sizes, char *str);
  
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

  /** \brief Using the commands stored in
      <tt>(n_entries,sizes,str)</tt>, apply the aliases stored in the
      \ref o2scl::cli object and return the counts
      <tt>(n_new,s_new)</tt> for memory allocation

      This function is used in \o2y in
      <tt>o2graph_plotter::parse_argv()</tt> to process aliases in the
      \c o2graph executable. It converts the input data to a
      <tt>vector&lt;string&gt;</tt> object, applies any aliases stored
      in the \ref o2scl::cli (or \ref o2scl::cli_readline) object, and
      then counts the new number of arguments (\c n_new) and string
      length for all arguments (\c s_new). These counts are used to
      allocate memory in Python in
      <tt>o2graph_plotter::parse_argv()</tt> in order to prepare for a
      call to \ref o2scl_acol_apply_aliases() .
   */
  void o2scl_acol_alias_counts(void *vp, int n_entries, int *sizes, 
                               char *str, int &n_new, int &s_new);

  /** \brief \brief Using the commands stored in
      <tt>(n_entries,sizes,str)</tt>, apply the aliases stored in the
      \ref o2scl::cli object and place them in the pre-allocated 
      arrays \c sizes_new and \c str_new

      This function is used in \o2y in
      <tt>o2graph_plotter::parse_argv()</tt> to process aliases in the
      \c o2graph executable. It converts the input data to a
      <tt>vector&lt;string&gt;</tt> object, applies any aliases stored
      in the \ref o2scl::cli (or \ref o2scl::cli_readline) object, and
      then stores the results in \c sizes_new and \c str_new (which
      are allocated beforehand in Python.
  */
  void o2scl_acol_apply_aliases(void *vp, int n_entries, int *sizes, 
                                char *str, int *sizes_new, char *str_new);
  
  /** \brief Set the command name, the short description,
      and the environment variable name
      
      This function is used in \o2y in
      <tt>o2graph_plotter::parse_argv()</tt> to communicate three
      strings which are used in the \ref o2scl_acol::acol_manager class.
  */
  void o2scl_acol_set_names(void *vp, int n1, char *cmd_name,
                            int n2, char *short_desc, int n3,
                            char *env_var);

  /** \brief Convert a rank 2 \ref o2scl::tensor (with data types
      \c double, \c int, or \c size_t) or \ref o2scl::tensor_grid
      object to a \ref o2scl::table3d object

      There are two sets of values for \c i1 and \c i2 which are
      allowed, either <tt>i1=0, i2=1</tt> or <tt>i1=1, i2=0</tt>, the
      latter of which corresponds to transposing the two indices.

      This function is used in o2graph_plotter::den_plot().
   */
  int o2scl_acol_tensor_to_table3d(void *vp, int i1, int i2);
  
  /** \brief Parse the set of commands in \c n_entries, \c sizes
      and \c str
      
      This function uses the executes the commands stored \c
      n_entries, \c sizes, and \c str using the \ref o2scl::cli object
      in \ref o2scl_acol::acol_manager as if they were typed on the command
      line.

      This function is used in \o2y in o2graph_plotter::set_wrapper(),
      o2graph_plotter::get_wrapper(), and o2graph_plotter::gen_acol().
   */
  void o2scl_acol_parse(void *vp, int n_entries, int *sizes, 
                        char *str);

  /** \brief Return the size and a pointer to the column
      named \c col_name in a \ref o2scl::table object

      This function is used in o2graph_plotter::plot(),
      o2graph_plotter::plot1(), o2graph_plotter::rplot(),
      o2graph_plotter::scatter(), o2graph_plotter::histplot(),
      o2graph_plotter::hist2dplot(), and o2graph_plotter::errorbar().
   */
  int o2scl_acol_get_column(void *vp, char *col_name,
                            int &n, double *&ptr);

  /** \brief Return the size and a pointer to the row
      with index \c row_index in a \ref o2scl::table object
      
      \note This function is currently unused. It may have been a
      precursor for a mult-vector-spec?

      \note Deprecating this for now.
   */
  //int o2scl_acol_get_row_ser(void *vp, char *parttern, int row_index,
  //int &n, double *&ptr);
  
  /** \brief Return the size and a pointer to a double array
      corresponding to a <tt>int[]</tt>, <tt>size_t[]</tt>, or
      <tt>double[]</tt> object

      This function is used in o2graph_plotter::plot1().
   */
  int o2scl_acol_get_double_arr(void *vp, int &n, double *&ptr);

  /** \brief Return the sizes, grid, and data pointer for 
      a rank 3 \ref o2scl::tensor_grid object

      This function is used in <tt>o2graph_plotter</tt> for 
      <tt>yt-add-vol</tt>.
   */
  int o2scl_acol_get_tensor_grid3(void *vp, int &nx, int &ny,
                                  int &nz, const double *&xg,
                                  const double *&yg,
                                  const double *&zg, const double *&data);
  
  /** \brief Return the size and a pointer to the histogram
      representative x values in a \ref o2scl::hist object

      This function is used in o2graph_plotter::plot() and
      o2graph_plotter::hist_plot().
   */
  int o2scl_acol_get_hist_reps(void *vp, int &n, double *&ptr);

  /** \brief Return the size and a pointer to the histogram bin edges
      in a \ref o2scl::hist object

      This function is used in o2graph_plotter::hist_plot().
   */
  int o2scl_acol_get_hist_bins(void *vp, int &n, double *&ptr);

  /** \brief Return the size and a pointer to the histogram weights in
      a \ref o2scl::hist object

      This function is used in o2graph_plotter::plot() and
      o2graph_plotter::hist_plot().
   */
  int o2scl_acol_get_hist_wgts(void *vp, int &n, double *&ptr);

  /** \brief Return the dimensionality, mesh size, and 
      lower and upper limits for a \ref o2scl::prob_dens_mdim_amr 
      object.

      This function is used in o2graph_plotter::plot().
   */
  int o2scl_acol_pdma_get_base(void *vp, int &ndim, int &n, 
                               double *&low, double *&high);

  /** \brief Return the lower and upper limits, fractional volume, and
      weight for the \ref o2scl::prob_dens_mdim_amr::hypercube object
      of index \c ix

      This function is used in o2graph_plotter::plot().
   */
  int o2scl_acol_pdma_get_cube(void *vp, int ix, 
                               double *&low, double *&high,
                               double &frac_vol, double &weight);

  /** \brief Return the number of contour lines associated with
      the current contour line vector object

      This function is used in o2graph_plotter::plot() and 
      o2graph_plotter::plotv().
   */
  int o2scl_acol_contours_n(void *vp);
  
  /** \brief For the current contour line vector object, set the
      pointers to the x- and y-values in the contour lines and return
      the contour level
  */
  double o2scl_acol_contours_line(void *vp, int i, int &n, double *&ptrx,
                                  double *&ptry);

  /** \brief Return the type of the current object 

      This function is used in o2graph_plotter::get_type(), 
      o2graph_plotter::den_plot(), o2graph_plotter::plot(), 
      o2graph_plotter::rplot(), o2graph_plotter::scatter(), 
      o2graph_plotter::histplot(), o2graph_plotter::hist2dplot(), 
      o2graph_plotter::errorbar(), o2graph_plotter::plot1(), 
      and o2graph_plotter::parse_string_list().
   */
  void o2scl_acol_get_type(void *vp, int &n, char *&str);

  /** \brief Return the size and a pointer to the slice
      named \c sl_name in a \ref o2scl::table object

      This function is used in o2graph_plotter::den_plot().
   */
  int o2scl_acol_get_slice(void *vp, char *slice_name,
                           int &nx, double *&xptr,
                           int &ny, double *&yptr,
                           double *&data);
  
  /** \brief For a two-dimensional histogram, return the bin edges,
      number of bins in both directions, and the weights in each bin

      This function is used in o2graph_plotter::den_plot().
   */
  int o2scl_acol_get_hist_2d(void *vp, 
                             int &nx, double *&xptr,
                             int &ny, double *&yptr,
                             double *&data);

  /** \brief Convert two multiple vector specifications to
      the a list of \ref o2scl::contour_line objects

      This function is used in o2graph_plotter::plotv().
  */
  int o2scl_acol_mult_vectors_to_conts(void *vp, char *str1,
                                       char *str2);

  /** \brief Get the list of parameters from the acol_manager cli
      object
   */
  int o2scl_acol_get_cli_parameters(void *vp, int &n, int *&sizes,
                                    char *&chlist);

  /** \brief Get the list of options/commands from the acol_manager
      cli object
   */
  int o2scl_acol_get_cli_options(void *vp, int &n, int *&sizes,
                                 char *&chlist);

  /** \brief Get the list of options/commands from the acol_manager
      cli object
   */
  int o2scl_acol_get_cli_options_type(void *vp, char *type,
                                      int &n, int *&sizes,
                                      char *&chlist);
  
  /** \brief Obtain the description of a parameter from the 
      acol_manager cli object
   */
  int o2scl_acol_cli_param_desc(void *vp, char *name, int &ndesc, 
                                    char *&chlist);
  
  /** \brief Obtain the short description of an option/command from
      the acol_manager cli object
   */
  int o2scl_acol_cli_option_desc(void *vp, char *name, int &ndesc, 
                                     char *&chlist);
  
  //@}
  
}

#endif
