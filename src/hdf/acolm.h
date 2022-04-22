/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
  This file is part of O₂scl.
  
  O₂scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O₂scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O₂scl. If not, see <http://www.gnu.org/licenses/>.

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
    
    /// The number formatter for html output
    o2scl::format_float ffl;

    /// Document strings for commands
    std::vector<std::vector<std::string>> cmd_doc_strings;
    
    /// Document strings for parameters
    std::vector<std::vector<std::string>> param_doc_strings;

    /// Document strings for help topics
    std::vector<std::vector<std::string>> help_doc_strings;
    
    /// Convert units object (initialized by constructor to global object)
    o2scl::convert_units<double> &cng;

    /// \name Parameters modifiable by the user
    //@{
    /// The output precision (default 6)
    int precision;

    /// True if we should make the output into neat columns (default true)
    bool pretty;
    
    /** \brief If true, output names at the top of some objects
     */
    bool names_out;

    /// True to use regex (false)
    bool use_regex;

    /// The name of the table
    std::string obj_name;
  
    /// Default arguments from environment
    std::string def_args;

    /** \brief The number of columns requested by the user 
        (default 0 for autodetect)
    */
    int ncols;
    
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
    o2scl::cli::parameter_string p_def_args;
    o2scl::cli::parameter_int p_verbose;
    o2scl::cli::parameter_int p_compress;
    o2scl::cli::parameter_int p_precision;
    o2scl::cli::parameter_int p_ncols;
    o2scl::cli::parameter_int p_interp_type;
    o2scl::cli::parameter_bool p_scientific;
    o2scl::cli::parameter_bool p_pretty;
    o2scl::cli::parameter_bool p_names_out;
    o2scl::cli::parameter_bool p_use_regex;
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
    o2scl::cli *cl;
    
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
    std::vector<std::vector<std::string>> vvstring_obj;

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

    /** \brief Desc
     */
    void xml_replacements(std::string &s,
                          std::vector<std::string> &clist);

    /** \brief Add new commands for type \c new_type
     */
    void command_add(std::string new_type);

    /** \brief Desc
     */
    void update_o2_docs(size_t narr,
                        o2scl::comm_option_s *options_arr,
                        std::string new_type="");
    
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

    /// \name Functions for the command-line interface
    //@{
    /** \brief Assign a constant

        For objects of type table:

        Assign a constant to the table, e.g. <tt>-assign pi
        "acos(-1)"</tt>.

        <name> val

        Assign a constant value to a name for the present table. Valid
        constant values are things like 1.618 or acos(-1.0) or
        sin(4^5). To remove an assignment, call assign with a blank
        value.
     */
    virtual int comm_assign(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert a series of histograms to a table3d object

        For objects of type table:

        Histogram series in a table3d object.

        <grid vector spec.> <direction (\"x\" or \"y\")> <grid name>
        <bin edges vector spec.> "<bin grid vector spec.> <bin name>
        <pattern> <new slice>

        Detailed desc.
     */
    virtual int comm_ser_hist_t3d(std::vector<std::string> &sv,
                                  bool itive_com);

    /** \brief Binary function for tensors

        For objects of type tensor_grid:

        Apply a binary function to two tensor_grid objects.

        <file> <object name> <function>

        Read tensor_grid named <object name> from file <file> and use
        it along with the function <function> to modify the current
        tensor_grid object. The <function> parameter should be a
        mathematical function of the value in the current tensor (v),
        the value in the tensor named <object name> (w), the indices
        (i0, i1, ...) or the grid points (x0, x1, ...).
     */
    virtual int comm_binary(std::vector<std::string> &sv, bool itive_com);

    /** \brief Average rows together 

        For objects of type table:

        Average rows of some or all columns together.

        <column or '*' for all> <window> [block averages]

        The first argument is the column to be modified. If the first
        argument is '*', then all columns are averaged. The second
        argument is the size of the window. If the third argument
        evaluates to false, then block averages instead of rolling
        averages are computed, and then the number of rows is divided
        by the window parameter. If block averages are requested, then
        the first argument must be '*'.
     */
    virtual int comm_average_rows(std::vector<std::string> &sv,
                                  bool itive_com);

    /** \brief Compute correlation

        For objects of type table:

        Compute the correlation coefficient between two columns.

        <column 1> <column 2>

        Compute the correlation coefficient between two columns, or,
        if no arguments are given, then compute the correlation
        coefficients between all pairs of columns.
     */
    virtual int comm_correl(std::vector<std::string> &sv, bool itive_com);

    /** \brief Refine an object

        For objects of type table:

        Refine the table.

        <index column> <factor>

        Detailed desc.
     */
    virtual int comm_refine(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the value of a constant expression.

        <expr>

        This computes the value of the constant expression <expr>.
        Examples are "calc acos(-1)" or "calc 2+1/sqrt(2.0e4)". To see
        valid expressions type 'acol -help <tt>functions</tt>'.

        Results are given at the current value of <tt>precision</tt>.
        Values of precision up to 50 are allowed, and multiprecision
        (rather than double precision) arithmetic is used if
        necessary. For example, try 'acol -set precision 45 -calc
        "acos(-1)"'.

        Constant values from the constant library (see 'acol -help
        <tt>constant</tt>') will automatically be used, so long as
        they have a unique value in MKS units. However, constant
        values are currently only stored to double precision, so they
        will result in an error if the value of <tt>precision</tt> is
        not larger than 15. Unicode is also supported for constants,
        so try, e.g. 'acol -set precision 15 -calc π'.
    */
    virtual int comm_calc(std::vector<std::string> &sv, bool itive_com);

    /** \brief Clear the current object

        (No arguments.)

        Deallocate the memory associated with the current object. This
        command does not clear the object name stored in \c obj_name.
     */
    virtual int comm_clear(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get help

        [command or parameter or type or topic]

        If no argument is specified, this outputs all of the commands
        which are valid for the current type, and then some generic
        help for using acol. If a command is specified as the
        argument, then help for that command is printed. If a
        parameter which is modifiable by the \c get or \c set commands
        is specified as the argument, then the help for that parameter
        is printed. If a type is specified as the argument, then a
        list of commands which operate on objects of that type is
        printed. Finally, if a help topic is specified as the
        argument, then that help information is printed.
     */
    virtual int comm_help(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief List commands, with an optional type argument

        [type]

        If no argument is specified, list all valid commands for the
        current type (including those commands which do not require a
        current object). If a type argument is given, then list all
        valid commands for the specified type.
     */
    virtual int comm_commands(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Create an object.

        <type> [...]

        Create a new object of type <type>. If an object is currently
        in memory, it is deallocated before creating the new object.

        "<tt>create</tt> <type> <value>": For types <tt>char</tt>,
        <tt>int</tt>, <tt>size_t</tt>, and <tt>string</tt>, create an
        object and give it the initial value specified.

        "<tt>create</tt> <tt>double</tt> <value spec.>": Create a
        <tt>double</tt> object and set it equal to the value specified by
        <value spec.>. (See "acol -help functions" for help on
        specifying functions and "acol -help value-spec" for help on
        value specifications.)

        "<tt>create</tt> <type> <size> <function of "i">": For array
        types <tt>int[]</tt> and <tt>size_t[]</tt>, the user must
        specify the size of the array and a function of the array
        index <tt>i</tt> to fill the array.

        "<tt>create</tt> <tt>double[]</tt> [<size> <function of "i">]
        or [vector spec.]": For <tt>double[]</tt> the user must either
        give a vector specification, or specify the size of the array
        and a function of the array index <tt>i</tt>.

        "<tt>create</tt> <tt>table</tt> <name> <vector spec.>":
        Create a new <tt>table</tt> object with one column named <name>
        from a vector specification (see ``Vector specifications``
        for the syntax).
        
        "<tt>create</tt> <tt>tensor</tt> <rank> <size 0> <size 1>
        ...": Create a <tt>tensor</tt> object with the specified rank and
        sizes. All tensor entries are initialized to zero.

        "<tt>create</tt> <tt>tensor_grid</tt> <rank> <size 0> <size 1>
        ...": <tt>Create a tensor_grid object</tt> with the specified
        rank and sizes. The tensor grid is initialized to count each
        index (beginning with zero) and the entries of the tensor are
        initialized to zero. The grid can be specified afterwards
        using <tt>set-grid</tt>.

        "<tt>create</tt> <tt>table3d</tt> <x name> <x vector spec.> <y
        name> <y vector spec.>\n <slice name> <slice function>":
        Create a new table3d object which has one slice. The x and y
        grids are given as vector specifications (see "acol -help
        vector-spec" for the syntax). The slice function can be
        written in terms of the x- and y-grid values which are
        referred to by name.

        For example, using <tt>o2graph</tt> from o2sclpy:

        o2graph -create table3d x func:100:i/200 y func:100:i/200 z
        "sin(1/(x+0.01))* sin(1/(y+0.01))" -den-plot z -xtitle x
        -ytitle y -show
    */
    virtual int comm_create(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the grid for a \ref o2scl::tensor_grid object

        For objects of type tensor_grid:

        Set the tensor grid.

        <index> <func. or vector spec>

        The first argument for the \c set-grid command specifies the
        index for which grid to set. The second argument specifies the
        grid. If it contains a ':', it is assumed to be a vector
        specification (see 'help vector-spec'). Otherwise, the
        argument is assumed to be a function which specifies the grid
        value as a function of the variables 'i' and 'x'. The value of
        'i' ranges from 0 to m-1, where 'm' is the tensor size for
        each rank and the value of 'x' is equal to the previous grid
        value.
     */
    virtual int comm_set_grid(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get the grid for a \ref o2scl::tensor_grid object

        For objects of type table3d:

        Print out the table3d grid.

        (No arguments.)

        Output the table3d grid as a series of columns.

        For objects of type tensor_grid:

        Get the tensor grid.

        (No arguments.)

        Output the tensor grid as a series of columns.
     */
    virtual int comm_get_grid(std::vector<std::string> &sv, bool itive_com);

    /** \brief Download a file from the specified URL.

        <file> <URL> [hash, \"file:\"hash_filename, or \"none\"] 
        [directory]

        Check if a file matches a specified hash, and if not, attempt
        to download a fresh copy from the specified URL. If the
        filename is "_", then the file is extracted from the end of
        the URL.
     */
    virtual int comm_download(std::vector<std::string> &sv, bool itive_com);

    /** \brief Parse doxygen XML to generate runtime docs.
        
        (No arguments.)

        When pugixml is enabled, this function reads the
        doxygen XML output and generates an HDF5 file which acol reads
        to generate the runtime documentation.
    */
    virtual int comm_xml_to_o2(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Open local HTML docs for acol or an O₂scl topic.

        [topic]

        If [topic] is unspecified, this command opens up the local
        HTML documentation for acol in the default web browser using
        'open' on OSX and 'xdg-open' on other systems. If a topic is
        specified, then the closest O₂scl documentation web page is
        opened. In order to open the remote version of the
        documentation instead of the local copy, use the
        <tt>wdocs</tt> command instead.
    */
    virtual int comm_docs(std::vector<std::string> &sv, bool itive_com);

    /** \brief Open remote HTML docs for acol or an O₂scl topic.

        [search_term], [topic] or [section search_term]
        
        If no arguments are given, this command opens up the remote
        HTML documentation for acol in the default web browser using
        'open' on OSX and 'xdg-open' on other systems. If a [topic] is
        specified, then the associated O₂scl web page is opened. If
        the argument does not match an already known topic, then the
        search feature on the O₂scl web page is opened using the
        specified search term. Note that, for search terms, spaces can
        be included using e.g. '-wdocs \"Simulated annealing\"'. Valid
        sections are either \"eos\" or \"part\". In order to open
        the local version of the documentation instead of the remote
        copy, use <tt>docs</tt> instead of <tt>wdocs</tt>.
    */
    virtual int comm_wdocs(std::vector<std::string> &sv, bool itive_com);

    /** \brief Delete a column

        For objects of type table:

        Delete a table column.

        <name>

        Delete the entire column named <name>.
     */
    virtual int comm_delete_col(std::vector<std::string> &sv, bool itive_com);

    /** \brief Delete rows

        For objects of type table:

        Delete rows selected by a function.

        <function>

        Delete the set of rows for which a function evaluates to a
        number greater than 0.5. For example, <tt>-delete-rows
        if(col1+col2>10,1,0)</tt> will delete all columns where the
        sum of the entries in \c col1 and \c col2 is larger than 10.
        See also the \c select-rows command.
     */
    virtual int comm_delete_rows(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Delete rows which match to within a specified tolerance

        For objects of type table:

        Delete rows which match to within a specified tolerance

        [relative tol.] [absolute tol.]

        This command deletes all rows which match within the specified
        tolerances. If verbose is larger than zero then information
        about how many rows were deleted is provided.
     */
    virtual int comm_delete_rows_tol(std::vector<std::string> &sv,
                                     bool itive_com);

    /** \brief Compute a derivative

        For objects of type table:

        Derivative of a function defined by two columns.

        <x> <y> <name>

        Create a new column named <name> filled with the derivative
        of the function y(x) obtained from columns <x> and <y>.

        For objects of type tensor:

        Compute the derivative of the tensor object w.r.t. an index

        <index>

        The <tt>deriv</tt> command differentiates the tensor object
        with respect to one of the indices.

        For objects of type tensor_grid:

        Compute the derivative of the tensor object w.r.t. an index

        <index>

        The <tt>deriv</tt> command differentiates the tensor object
        with respect to one of the indices.

        For objects of type double[]:

        Replace the array with its derivative.

        (No arguments.)

        Replace the array with its derivative using the current
        interpolation type.

        For objects of type int[]:

        Replace the array with its derivative.

        (No arguments.)

        Replace the array with its derivative using the current
        interpolation type, converting it to a double[].

        For objects of type size_t[]:

        Replace the array with its derivative.

        (No arguments.)

        Replace the array with its derivative using the current
        interpolation type, converting it to a double[].
     */
    virtual int comm_deriv(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table object

        For objects of type double[]:

        Convert to a table object.

        <column name>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type int[]:

        Convert to a table object.

        <column name>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type size_t[]:

        Convert to a table object.

        <column name>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type tensor_grid:

        Convert to a table object.

        <index> <grid name> <data name> [values of fixed indices]
        
        Detailed desc.

        For objects of type hist:

        Convert to a table object.

        (No arguments.)
        
        Convert to a table object.
     */
    virtual int comm_to_table(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get entries along the main diagonal

        For objects of type tensor:

        Get diagonal elements.

        (No arguments.)

        Extract only the elements on the main diagonal to create a
        double[] object.
     */
    virtual int comm_diag(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object

        For objects of type table:

        Convert a table to a table3d object.

        <x column> <y column> [empty value] [eps]

        The 'to-table3d' creates a table3d object using 'x column' and
        'y column' as the data for the x and y grids. If 'empty
        value', then this value is used for points not given by the
        table. If 'eps' is specified, then use that value as the
        minimum value between grid points.

        For objects of type tensor:

        Select two indices and convert to a table3d object.

        <x index> <y index> <slice name> [fixed 1] [fixed 2] ...

        This command uses two indices in the current tensor object to
        create a table3d object. The values for the remaining indices
        fixed to [fixed 1], [fixed 2], etc. in that order. For
        example, "to-table3d 3 1 z 5 3" uses index 3 for the x
        coordinate of the new table3d object, uses index 1 for the y
        coordinate of the new table3d object, uses 5 for index 0, and
        uses 3 for index 2."+ The x- and y-grids in he table3d object
        are named "x" and "y" and filled with the grid index by
        default."+ To set the x- or y-grid names afterwards, use
        commands 'x-name' and 'y-name'.

        For objects of type tensor<int>:

        Select two indices and convert to a table3d object.

        <x index> <y index> <slice name> [fixed 1] [fixed 2] ...

        This command uses two indices in the current tensor object to
        create a table3d object. The values for the remaining indices
        fixed to [fixed 1], [fixed 2], etc. in that order. For
        example, "to-table3d 3 1 z 5 3" uses index 3 for the x
        coordinate of the new table3d object, uses index 1 for the y
        coordinate of the new table3d object, uses 5 for index 0, and
        uses 3 for index 2."+ The x- and y-grids in he table3d object
        are named "x" and "y" and filled with the grid index by
        default."+ To set the x- or y-grid names afterwards, use
        commands 'x-name' and 'y-name'.

        For objects of type tensor<size_t>:

        Select two indices and convert to a table3d object.

        <x index> <y index> <slice name> [fixed 1] [fixed 2] ...

        This command uses two indices in the current tensor object to
        create a table3d object. The values for the remaining indices
        fixed to [fixed 1], [fixed 2], etc. in that order. For
        example, "to-table3d 3 1 z 5 3" uses index 3 for the x
        coordinate of the new table3d object, uses index 1 for the y
        coordinate of the new table3d object, uses 5 for index 0, and
        uses 3 for index 2."+ The x- and y-grids in he table3d object
        are named "x" and "y" and filled with the grid index by
        default."+ To set the x- or y-grid names afterwards, use
        commands 'x-name' and 'y-name'.

        For objects of type tensor_grid:

        Select two indices and convert to a table3d object.

        <x index> <y index> <slice name> [value 1] [value 2] ...

        This command uses two indices in the current tensor_grid
        object to create a table3d object. The values for the
        remaining indices are by interpolation to [value 1], [value
        2], etc. in that order. For example, "to-table3d 3 1 z 0.5
        2.0" uses index 3 for the x coordinate of the new table3d
        object, uses index 1 for the y coordinate of the new table3d
        object, uses interpolation to set the value of the index 0 to
        0.5, and uses interpolation to set the value of index 2 to to
        2.0. The x- and y-grids in the table3d object are named "x"
        and "y" by default. To set the x- or y-grid names
        afterwards, use commands 'x-name' and 'y-name'.

        For objects of type hist_2d:

        Convert to a hist_2d object.

        <x name> <y name> <weight name>

        Convert to a hist_2d object using the specified names.
        
        For objects of type prob_dens_mdim_amr:

        Select two indices and convert to a table3d object.

        <x index> <y index> <x name> <x points> <y name> <y points>
        <slice name>

        Select two indices and convert to a table3d object.
     */
    virtual int comm_to_table3d(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::tensor_grid object

        For objects of type table3d:

        Convert a slice of the table3d object to a tensor_grid object.

        <slice>

        Detailed desc.

        For objects of type tensor:

        Convert the tensor to a tensor_grid object.

        [function 1] [function 2] ...

        Convert a tensor to a tensor_grid object, using functions to
        specify the grid for each index. The functions should be
        specified as functions of the variable 'i', which runs from 0
        to size-1 for each index. Any user-specified functions are
        used up to the rank of the tensor, and if not enough functions
        are specified, then the function 'i' is used.
     */
    virtual int comm_to_tensor_grid(std::vector<std::string> &sv,
                                    bool itive_com);
    
    /** \brief Convert object to a \ref o2scl::tensor object

        For objects of type tensor_grid:

        Convert to a tensor object.

        (No arguments.)
        
        Convert to a tensor object, ignoring the grid.
     */
    virtual int comm_to_tensor(std::vector<std::string> &sv,
                               bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object
        by summing over tensor indices

        For objects of type tensor:

        Select two indices and convert to a table3d object.

        <x name> <y name> <slice name> [fixed 1] [fixed 2] ...

        Detailed desc.
    */
    virtual int comm_to_table3d_sum(std::vector<std::string> &sv,
                                    bool itive_com);

    /** \brief Compute the autocorrelation coefficients

        If there is no current object:

        Compute autocorrelation coefficients from a set of vectors

        <mult. vec. spec. 1> [mult. vec. spec. 2]
        
        This command computes the autocorrelation coefficients for all
        vectors specified as multiple vector specifications in the
        arguments, then averages those autocorrelation coefficients
        together. The averaged autocorrelation coefficients are kept
        as a new <tt>double[]</tt> object. See ``Multiple vector
        specifications`` for more information.

        For objects of type int[]:

        (No arguments.)

        Replace the current object with a <tt>double[]</tt> object which
        contains the autocorrelation coefficient as a function of the
        step size.

        For objects of type double[]:

        (No arguments.)

        Replace the current object with a <tt>double[]</tt> object which
        contains the autocorrelation coefficient as a function of the
        step size.

        For objects of type table:

        <ac> <ftom> <column or vector specification> [second column or vector
        specification] ... 

        Compute autocorrelation coefficients from a column of a table.
        The first argument, <ac>, is the name of the column in which
        the autocorrelation coefficients will be stored. The second
        argument, <ftom>, is the name of the column in which the
        quantity '5*tau/M' will be stored. The data may be either a
        column in the table or a vector specification. Columns <ac>
        and <ftom> are created if they are not already present and
        overwritten if they already contain data. Also, the
        autocorrelation length and estimated sample size are output to
        the screen. If multiple data sources are given, then the
        autocorrelation coefficients are averaged together. See also
        ``Vector specifications`` for more information on the third
        and fourth arguments.

    */
    virtual int comm_autocorr(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the autocorrelation coefficient using \c acor

        For objects of type table:

        Compute the autocorrelation coefficient using \c acor

        <column>

        Detailed desc.
     */
    virtual int comm_ac_len(std::vector<std::string> &sv,
                            bool itive_com);

    /** \brief Create a slice which is the derivative wrt x of another

        For objects of type table3d:

        Derivative with respect to x.

        <f> <dfdx>

        Create a new slice named <dfdx> filled with the derivative of
        the function from the x grid and slice named <f>.
     */
    virtual int comm_deriv_x(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a slice which is the derivative wrt y of another

        For objects of type table3d:

        Derivative with respect to y.

        <f> <dfdy>

        Create a new slice named <dfdy> filled with the derivative of
        the function from the y grid and slice named <f>.
     */
    virtual int comm_deriv_y(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute a second derivative

        For objects of type table:

        Second derivative of a function defined by two columns.

        <x> <y> <name>

        Create a new column named <name> filled with the second
        derivative of the function y(x) obtained from columns <x> and
        <y>.
     */
    virtual int comm_deriv2(std::vector<std::string> &sv, bool itive_com);

    /** \brief List objects in a HDF5 file

        <file>

        This lists all the top-level datasets and groups in a HDF5
        file and, for those groups which are in the O₂scl format,
        gives the type and name of the object stored in that HDF5
        group.
     */
    virtual int comm_filelist(std::vector<std::string> &sv, bool itive_com);

    /** \brief Read an object from an O₂scl-style HDF5 file.

        <file> [object name]

        Read an HDF5 file with the specified filename. If the [object
        name] argument is specified, then read the object with the
        specified name. Otherwise, look for the first <tt>table</tt>
        object, and if not found, look for the first <tt>table3d</tt>
        object, and so on, attempting to find a readable O₂scl object.
    */
    virtual int comm_read(std::vector<std::string> &sv, bool itive_com);

    /** \brief Add 'nlines' as a constant to a \ref o2scl::table object

        For objects of type table:

        Add 'nlines' as a constant to a table object.

        (No arguments.)

        Add a constant called 'nlines' to the table and set it equal
        to the number of lines (rows) in the table.
     */
    virtual int comm_nlines(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert to a \ref o2scl::hist object

        For objects of type table:

        Convert a table to a histogram.

        <col> <n_bins> [wgts]

        The 'to-hist' command creates a 1D histogram from 'col' using
        exactly 'n_bins' bins and (optionally) weighting the entries
        by the values in column 'wgts'. The second form creates a 2D
        histogram from 'col1' and 'col2' using N1 bins in the x
        direction and N2 bins in the y direction, optionally weighting
        the entries by the column 'wgts'.
     */
    virtual int comm_to_hist(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert to a \ref o2scl::hist_2d object

        For objects of type table:

        Convert a table to a 2d histogram.

        <col x> <col y> <n_x_bins> <n_y_bins> [wgts]

        The 'to-hist-2d' command creates a 2D histogram from 'col x'
        and 'col y' using 'n_x_bins' bins in the x direction and
        'n_y_bins' bins in the y direction, optionally weighting the
        entries by the column 'wgts'.

        For objects of type table3d:

        Convert a table3d slice to a 2d histogram.

        <slice>

        The 'to-hist-2d' command creates a 2D histogram from slice
        <slice>.
     */
    virtual int comm_to_hist_2d(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output the type of the current object

        (No arguments.)

        Show the current object type, either <tt>table</tt>,
        <tt>table3d</tt>, <tt>hist</tt>, <tt>hist_2d</tt>,
        <tt>vector<contour_line></tt>, <tt>int</tt>, <tt>double</tt>,
        <tt>char</tt>, <tt>string</tt>, <tt>int[]</tt>,
        <tt>double[]</tt>, <tt>string[]</tt>, <tt>size_t</tt>,
        <tt>size_t[]</tt>, <tt>uniform_grid<double></tt>,
        <tt>tensor_grid</tt>, <tt>tensor</tt>, <tt>tensor<int></tt>,
        <tt>tensor<size_t></tt> or <tt>prob_dens_mdim_amr</tt>.
     */
    virtual int comm_type(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Find a row

        For objects of type table:

        Find a row which maximizes a function.

        <func> or find-row <col> <val>

        If one argument is given, then find-row finds the row which
        maximizes the value of the expression given in <func>, and
        then output the entire row. Otherwise find-row finds the row
        for which the value in column named <col> is as close as
        possible to the value <val>. See command 'get-row' to get a
        row by it's index.
     */
    virtual int comm_find_row(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Create a column from a function

        For objects of type table:

        Create a column from a function

        <func> <name>

        Set the column named <name> to the result of a function,
        <func>, in terms of the other columns. If the column does not
        already exist, a new one is added to the table. For example,
        for a table containing columns named 'c1' and 'c2', 'function
        c1-c2 c3' would create a new column c3 which contains the
        difference of columns 'c1' and 'c2'.

        For objects of type double[]:

        Set the values of the array given a function.

        <function>

        Set the values of the array given a user-specified function of
        'i'. For example, "(sin(i)>1)*4".

        For objects of type int[]:

        Set the values of the array given a function.

        <function>

        Set the values of the array given a user-specified function of
        'i'. For example, "(sin(i)>1)*4".

        For objects of type size_t[]:

        Set the values of the array given a function.

        <function>

        Set the values of the array given a user-specified function of
        'i'. For example, "(sin(i)>1)*4".

        For objects of type hist:

        Apply a function to the weights.

        <function>

        Apply a function to the weights.

        For objects of type table3d:

        Create a new slice from a function.

        <func> <name>

        Set the slice named <name> to the result of a function,
        <func>, in terms of the other slices. If the slice does not
        already exist, a new one is created. For example, for a
        table3d containing slices named 's1' and 's2', 'function s1-s2
        s3' would create a new column 's3' which contains the
        difference of columns 's1' and 's2'.

        For objects of type tensor:

        Set the tensor values from a function.

        [cond. function] <function of v, i0, i1, ...>

        The \c function command sets all entries in a tensor equal to
        a user-specified mathematical function of the indices. When
        the conditional function evaluates to a number less than or
        equal to 0.5, then the tensor entry will be unchanged. (For
        more help with functions, type acol -help functions)

        For objects of type tensor_grid:

        Set the tensor values from a function.

        [conditional func.] <func. of v, i0, i1, ... and x0, x1, ...>

        The "function" command sets all the data entries in a
        tensor_grid equal to a user-specified mathematical function of
        the value in the tensor (v), the indices (i0, i1, ...) or the
        grid points (x0, x1, ...). If two function arguments are given
        and if the first function argument is not "none", then the
        first function specifies which tensor entries are to be
        modified. When the conditional function evaluates to a number
        less than or equal to 0.5, then the tensor entry will be "
        unchanged. (For more help with functions, type "acol -help
        functions.".)
    */
    virtual int comm_function(std::vector<std::string> &sv, bool itive_com);

    /** \brief Add a vector_specification

        For objects of type table:

        Add column from a vector specification to the table.

        <vec. spec.> <column>

        Detailed spec.
     */
    virtual int comm_add_vec(std::vector<std::string> &sv, bool itive_com);

    /** \brief Read an object generic text file
        
        <type> <file>
        
        Read an object of type <type> from a text file named <file>.
        The allowed text file formats depend on the particular type
        specified.

        For <tt>int</tt>, <tt>char</tt>, <tt>double</tt>, or 
        <tt>size_t</tt> objects, the
        file is assumed to begin with the desired object and it is
        read using operator>>().

        For <tt>string</tt> objects, the first line is read using
        <tt>std::getline()</tt>.

        For array objects, it is assumed that all array entries are on
        the first line of the file and no carriage returns are present
        between entries.

        For <tt>table</tt> objects, the first line of the file must
        either contain numeric data or column names separated by white
        space, without carriage returns, except for the one at the end
        of the line. If the first line contains column names, the
        second line may optionally contain unit expressions for each
        column, enclosed by square brackets. All remaining lines are
        assumed to contain data with the same number of columns as the
        first line.
        
        For <tt>table3d</tt> objects, the data must be stored in
        columns with the first column specifying the x-axis grid point
        and the second column specifying the y-axis grid point. The
        remaining columns give the data for each slice at that point.
        Each grid point must correspond to a row in the file, but the
        lines need not be in any particular order. The columns may
        have one header line at top which specifies the names of the
        x- and y-grids and the names of each slice (in order).
     */
    virtual int comm_generic(std::vector<std::string> &sv, bool itive_com);

    /** \brief Print out an entire row

        For objects of type table:

        Get a row by index.

        <index>

        Get a row by index. The first row has index 0, and the last
        row has index n-1, where n is the total number of rows as
        returned by the 'list' command. The 'index' command creates a
        column of row indexes. To find a row which contains a
        particular value or maximizes a specified function, use
        'find-row'.
     */
    virtual int comm_get_row(std::vector<std::string> &sv, bool itive_com);

    /** \brief Extract a slice from a table3d object to generate a 
        \ref o2scl::table object

        For objects of type table3d:

        Convert a slice to a table object.

        <"x" or "y"> <value>

        Extract a slice of a table3d object at fixed x or fixed y to
        create a new table object. This function uses interpolation
        with the current interpolation type to interpolate all of the
        slices in the table3d object to create a table with a column
        for each slice.

        For objects of type tensor_grid:

        Slice to a smaller rank tensor_grid object.

        <index 1> <value 1> <index 2> <value 2> ...

        Detailed desc.
    */
    virtual int comm_slice(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert a slice to a histogram

        For objects of type table3d:

        Construct a histogram from a slice.

        <slice>

        Detailed desc.
     */
    virtual int comm_slice_hist(std::vector<std::string> &sv, bool itive_com);

    /** \brief Fit data to a function

        For objects of type table:

        Fit two columns to a function (experimental).

        <x> <y> <yerr> <ynew> <par names> <func> <vals>

        Detailed desc.
    */
    virtual int comm_fit(std::vector<std::string> &sv, bool itive_com);

    /** \brief Insert a column from an external table using 
        interpolation

        For objects of type table:

        Interpolate a column from another file.

        <file> <table> <oldx> <oldy> <newx> [newy]

        Insert a column from file <fname> interpolating it into the
        current table. The column <oldy> is the columns in the file
        which is to be inserted into the table, using the column
        <oldx> in the file and <newx> in the table. The new column in
        the table is named <oldy>, or it is named [newy] if the
        additional argument is given.

        For objects of type table3d:

        Interpolate a slice from another file.

        <file> <table> <old> [new]
    */
    virtual int comm_insert(std::vector<std::string> &sv, bool itive_com);

    /** \brief Insert an external table using interpolation

        For objects of type table:

        Insert a table from another file.

        <fname> [table name] [old_x new_x]

        Insert all columns from file <fname> into the current table.
        The first table is used or the table object named table_name,
        if specified. If index columns old_x and new_x are not
        specified, then the insert requires both the current and the
        source table to have the same number of rows. If they are
        specified, then interpolation using those index columns is
        used. If columns in the new table are not present in the
        current table, then they are added automatically. If a column
        in the current table has the same name as one in the new table
        then it is rewritten with new data, with one exception. If a
        column in the new table has the same name as old_x, then it is
        left unmodified.
     */
    virtual int comm_insert_full(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a column which is the integral of another

        For objects of type table:

        Integrate a function specified by two columns.

        <x> <y> <name>

        Create a new column named <name> filled with the integral of
        the function y(x) obtained from columns <x> and <y>.
     */
    virtual int comm_integ(std::vector<std::string> &sv, bool itive_com);

    /** \brief Toggle interactive mode

        (No arguments.)

        If given as a command-line parameter, 'interactive' toggles
        the execution of the interactive mode after the command-line
        parameters are processed. If zero arguments are given to
        'acol' on the command-line then the interactive interface is
        automatically turned on.
     */
    virtual int comm_interactive(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output current object in the internal HDF5 format.
        
        [file]

        Output the current object in the internal HDF5 format.
        If no argument is given, then output is sent to the screen,
        otherwise, output is sent to the specified file. 
     */
    virtual int comm_internal(std::vector<std::string> &sv, bool itive_com);

    /** \brief Perform an interpolation using the current object
        
        For objects of type table:

        Interpolate a number into a column.

        <x name> <x value> <y name>

        Interpolate <x value> from column named <x name> into column
        named <y name>.

        For objects of type double[]:

        Interpolate an index into the array

        <x value>

        Interpolate <x value> in the array.

        For objects of type int[]:

        Interpolate an index into the array

        <x value>

        Interpolate <x value> in the array and print out the result
        as a double.

        For objects of type size_t[]:

        Interpolate an index into the array

        <x value>

        Interpolate <x value> in the array and print out the result
        as a double.

        For objects of type table3d:

        Interpolate x and y values into a slice.

        <z name> <x value> <y value>

        Interpolate (<x value>,<y value>) into the slice named <z
        name>.

        For objects of type tensor_grid:

        Linearly interpolate in the grid.

        <value 1> <value 2> <value 3> ...

        The command "interp" uses linear interpolation to interpolate
        an array with size equal to the tensor rank into the tensor
        grid and outputs the result.
    */
    virtual int comm_interp(std::vector<std::string> &sv, bool itive_com);

    /** \brief List properties of an object

        For objects of type table:

        List the constants, column names and other info.

        (No arguments.)

        List the constants, column names and other info.

        For objects of type table3d:

        List the slice names and print out grid info.

        (No arguments.)

        List the slice names and print out grid info.

        For objects of type tensor:

        List the tensor rank and index sizes.

        (No arguments.)

        List the tensor rank and index sizes.

        For objects of type tensor<int>:

        List the tensor rank and index sizes.

        (No arguments.)

        List the tensor rank and index sizes.

        For objects of type tensor<size_t>:

        List the tensor rank and index sizes.

        (No arguments.)

        List the tensor rank and index sizes.

        For objects of type tensor_grid:

        List the tensor rank and index sizes.

        (No arguments.)

        List the tensor rank and index sizes.

        For objects of type hist_2d:

        List the bin edges.

        (No arguments.)

        For objects of type hist_2d:
    */
    virtual int comm_list(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the maximum value of a column

        For objects of type table:

        Compute the maximum value of a column.

        <column name>

        Compute the maximum value of a column.

        For objects of type double[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type int[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type size_t[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type table3d:

        Compute the maximum value of a slice.

        <slice name>

        Compute the maximum value of a slice.

        For objects of type tensor:

        Compute the maximum value and the associated index.

        (No arguments.)
        
        Compute the maximum value and the associated index.

        For objects of type tensor<int>:

        Compute the maximum value and the associated index.

        (No arguments.)
        
        Compute the maximum value and the associated index.

        For objects of type tensor<size_t>:

        Compute the maximum value and the associated index.

        (No arguments.)
        
        Compute the maximum value and the associated index.

        For objects of type tensor_grid:

        Compute the maximum value and the associated index and grid point.

        (No arguments.)
        
        Compute the maximum value and the associated index and grid point.

        For objects of type hist_2d:

        Find the maximum weight.

        (No arguments.)

        Find the maximum weight and print out the location.
    */
    virtual int comm_max(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Compute the minimum value of a column

        For objects of type table:

        Compute the minimum value of a column.

        <column name>

        Compute the minimum value of a column.

        For objects of type double[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type int[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type size_t[]:

        Compute the maximum value and the associated index.

        (No arguments.)

        Compute the maximum value and the associated index.

        For objects of type table3d:

        Compute the minimum value of a slice.

        <slice name>

        Compute the minimum value of a slice.

        For objects of type tensor:

        Compute the minimum value and the associated index.

        (No arguments.)
        
        Compute the minimum value and the associated index.

        For objects of type tensor<int>:

        Compute the minimum value and the associated index.

        (No arguments.)
        
        Compute the minimum value and the associated index.

        For objects of type tensor<size_t>:

        Compute the minimum value and the associated index.

        (No arguments.)
        
        Compute the minimum value and the associated index.

        For objects of type tensor_grid:

        Compute the minimum value and the associated index and grid point.

        (No arguments.)
        
        Compute the minimum value and the associated index and grid point.

        For objects of type hist_2d:

        Find the minimum weight.

        (No arguments.)

        Find the minimum weight and print out the location.
    */
    virtual int comm_min(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the name of the x grid

        For objects of type table3d:

        Get or set the name of the x grid.

        [name]

        Get or set the name of the x grid.
     */
    virtual int comm_x_name(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the name of the y grid

        For objects of type table3d:

        Get or set the name of the y grid.

        [name]

        Get or set the name of the y grid.
     */
    virtual int comm_y_name(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Add a column for line numbers

        For objects of type table:

        Add a column containing the row numbers.

        [column name]
        
        Define a new column named [column name] and fill the column
        with the row indexes, beginning with zero. If no argument is
        given, the new column is named 'N'.
    */
    virtual int comm_index(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output the current object to screen or text file.

        [file]

        Output the object to the screen, or if the [file] argument is
        specified, to a file. This is the same format as can be read
        using the 'generic' command.
     */
    virtual int comm_output(std::vector<std::string> &sv, bool itive_com);

    /** \brief Rearrange a tensor

        For objects of type tensor:

        Rearrange the tensor

        <index spec. 1> [index spec. 2] ...

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". See
        '-help index-spec for more information on the tensor index
        specifications.

        For objects of type tensor<int>:

        Rearrange the tensor.

        <index spec. 1> [index spec. 2] ...

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". See
        '-help index-spec for more information on the tensor index
        specifications.

        For objects of type tensor<size_t>:

        Rearrange the tensor.

        <index spec. 1> [index spec. 2] ...

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". See
        '-help index-spec for more information on the tensor index
        specifications.

        For objects of type tensor_grid:

        Rearrange the tensor_grid object.

        <index spec. 1> [index spec. 2] ...

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), range(ix,start,end),
        interp(ix,value), grid(ix,begin,end,n_bins,log), and
        gridw(ix,begin,end,bin_width,log). Index specifications may be
        specified as separate arguments e.g. "index(1)" "fixed(2,10)"
        or multiple index specifications may be given in a single
        argument separated by spaces or commas, e.g. "index(1)
        fixed(2,10)" or "index(1),fixed(2,10)". See '-help index-spec'
        for more information on the tensor index specifications.
    */
    virtual int comm_rearrange(std::vector<std::string> &sv, bool itive_com);

    /** \brief Preview the current object

        [number of lines] [number of columns]

        Print out all or part of the current object in format suitable
        for the screen.
    */
    virtual int comm_preview(std::vector<std::string> &sv, bool itive_com);

    /** \brief Send a slack message

        ["#channel"] <strings-spec>

        Send a message to slack, using the specified channel. If the
        channel is not specified, it is taken from the environment
        variable O2SCL_SLACK_CHANNEL. The '#' sign should be included
        with the channel name. The Slack webhook URL is taken from the
        environment variable O2SCL_SLACK_URL and the username is taken
        from the environment variable O2SCL_SLACK_USERNAME. The
        message is constructed from the string list specification in
        <strings-spec> (see 'acol -help strings-spec' for more
        information).
     */
    virtual int comm_slack(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get or set the value of an object

        For objects of type int:

        [value]

        Get or set the integer.

        For objects of type size_t:

        [value]

        Get or set the size_t object.

        For objects of type string:

        [value]

        Get or set the string.

        For objects of type double:

        [value spec.]

        Get or set the value of the double object. See ``Value
        specifications`` for more information.

        For objects of type char:

        [value]

        Get or set the character.

     */
    virtual int comm_value(std::vector<std::string> &sv, bool itive_com);

    /** \brief Concatenate two objects

        For objects of type table:

        Concatenate a second table object onto current table.

        <file> [name]

        Add a second table to the end of the first, creating new
        columns if necessary.

        For objects of type table3d:

        Concatenate data from a second table3d onto current table3d.

        <file> [name]

        Add all slices from the second table3d object which aren't
        already present in the current table3d object.

     */
    virtual int comm_cat(std::vector<std::string> &sv, bool itive_com);

    /** \brief Sum two objects

        For objects of type table:

        Add data from a second table object to current table.

        <file> [name]

        Add all columns from the second table to their corresponding
        columns in the current table, creating new columns if
        necessary.

        For objects of type double[]:

        Compute the vector sum.

        (No arguments.)

        Compute the vector sum.

        For objects of type int[]:

        Compute the vector sum.

        (No arguments.)

        Compute the vector sum.

        For objects of type size_t[]:

        Compute the vector sum.

        (No arguments.)

        Compute the vector sum.

        For objects of type table3d:

        Add data from a second table3d object to current table3d.
        
        <file> [name]

        Add all slides from the second table3d to their corresponding
        slices in the current table3d, creating new slices if
        necessary.

        For objects of type tensor:

        Output the sum of all the tensor entries.

        (No arguments.)

        The "sum" command outputs the total tensor size
        and the sum over all entries. Note, to perform a partial
        sum over sum of the tensor indices, use the 
        <tt>rearrange</tt> command.

        For objects of type tensor_grid:

        Output the sum of all the tensor entries.

        (No arguments.)

        The "sum" command outputs the total tensor size
        and the sum over all entries. Note, to perform a partial
        sum over sum of the tensor indices, use the 
        <tt>rearrange</tt> command.
    */
    virtual int comm_sum(std::vector<std::string> &sv, bool itive_com);

    /** \brief Rename a part of an object

        For objects of type table:

        Rename a column.

        <old> <new>

        Rename a column from <old> to <new>. Note that to rename the
        entire object, you should use <tt>-set obj_name new_name</tt>.

        For objects of type table3d:

        Rename a slice.

        <old> <new>

        Rename a slice from <old> to <new>. Note that to rename the
        entire object, you should use <tt>-set obj_name new_name</tt>.
     */
    virtual int comm_rename(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select part of an object

        For objects of type table:

        Select columns for a new table.

        <column 1> [column 2] ...

        Select creates a new table from the present table, including
        only the columns specified in <cols>. The column specification
        is a list of column names, functions, or patterns which match
        the column names. Patterns must be preceeded by a colon ':'
        and use ECMAScript regular expressions. All of the rows of
        data are copied over. If functions are specified, the result
        can be named using '='.

        For objects of type table3d:

        Select columns for a new table3d.

        <slice 1> [slice 2] ...

        Select creates a new table3d from the present table3d,
        including only the slices specified in <slice spec.>. The
        slice specification is a list of slice names, functions, or
        patterns which match the slice names. Patterns must be
        preceeded by a colon ':' and can use wildcards like '*' and
        '?'. All of the rows of data are copied over. If functions are
        specified, the result can be named using '='.

     */
    virtual int comm_select(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select rows from an object

        For objects of type table:
        
        <row specification>

        Select the rows from a table for which the row specification
        in <row_spec> evaluates to a number greater than 0.5.
    */
    virtual int comm_select_rows(std::vector<std::string> &sv,
                                  bool itive_com);

    /// Post-processing for setting a value
    virtual int comm_set(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set an individual data point at a specified row and
        column

        For objects of type table:

        Set the entries of a column.

        <row_spec> <col> <val_spec>

        Set the value of rows specifed by the 'row_spec' function in
        column 'col' to the value given by the 'val_spec' function.
        Rows are chosen if row_spec evaluates to a number greater than
        0.5.

        For objects of type table3d:

        Set the entries of a slice.

        <x value> <y value> <z name> <val>

        Set the value of the slice named 'z name' at the grid point
        closest to (<x value>,<y value>) to the value <val>.
     */
    virtual int comm_set_data(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Set units of a column

        For objects of type table:

        Set the units for a specified column.

        <column> <unit>

        Detailed desc.
     */
    virtual int comm_set_unit(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Compute contour lines

        For objects of type table3d:

        Create contour lines from a table3d slice.

        <value> <slice_name> [output_filename object_name]

        The "contours" command constructs a set of contour lines using
        the data in slice named <slice> at the fixed value given in
        <value>. If two additional arguments are given, then the
        contour lines are stored in the file named output_filename and
        the object is named object_name. If the file does not exist,
        it is created. If no contours are found, then no file I/O is
        performed and the current table3d object is unmodified.

        For objects of type hist_2d:

        Create contour lines from a table3d slice.

        ["frac"] <value> [output file] [output name]

        If the argument "frac" is not present, the "contours" command
        constructs a set of contour lines using at the fixed value
        given in <value>. If two additional arguments are given, then
        the contour lines are stored in the file named output_filename
        and the object is named object_name. If the file does not
        exist, it is created. If no contours are found, then no file
        I/O is performed and the current table3d object is
        unmodified."+ If the argument "frac" is present, then the
        operation is the same except that <value> is interpreted as a
        fraction of the total integral under the data.
     */
    virtual int comm_contours(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Get units of a column

        For objects of type table:

        Get the units for a specified column.

        <column>

        Obtains the units for the specified column.
     */
    virtual int comm_get_unit(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Get or set an entry

        For objects of type table:

        Get or set a single entry in a table.

        <column> <row> [value or "none"]

        This command gets or sets the value in the specified column
        and row. If "none" is specified as the third argument, then
        "entry" just prints out the specified entry as if the third
        argument was not specified.

        For objects of type table3d:

        Get or set a single entry in a table3d object.

        <slice> <x index> <y index> [value or "none"]

        Detailed desc.

        For objects of type tensor:

        Get or set a single entry in a tensor object.

        <index 1> <index 2> <index 3> ... [value or "none"]

        Detailed desc.

        For objects of type tensor_grid:

        Get or set a single entry in a tensor_grid object.

        <index 1> <index 2> <index 3> ... [value or "none"]

        The \"entry\" command gets or sets a value in the tensor_grid
        object. The arguments are a list of indices and (optionally) a
        new value to store in that location.
    */
    virtual int comm_entry(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get an entry by grid point

        For objects of type table:

        Get or set a single entry in a table.
        
        <index column> <index value> <target column> [value or "none"]
        
        The "entry-grid" command first looks for the value closest to
        <index value> in the column <index column> to determine a row
        in the table. Next "entry-grid" gets or sets the value of the
        target column in that row. If "none" is specified as the
        fourth argument, then "entry" just prints out the specified
        entry as if the third argument was not specified.

        For objects of type table3d:

        Get a single entry in a table3d object.

        <slice> <x value> <y value> [value or "none"]

        Detailed desc.

        For objects of type tensor_grid:

        Get a single entry in a tensor_grid object.

        <value 1> <value 2> <value 3> ... [value or "none"]

        The "entry-grid" command gets or sets a value in the
        tensor_grid object. The arguments are a list of grid values
        and (optionally) a new value to store in the location closest
        to the specified grid values.
    */
    virtual int comm_entry_grid(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Convert units

        For objects of type table:

        Convert a column to a new unit.

        <column> <new_unit>

        Convert the units of a column to <new unit>, multipliying all
        entries in that column by the appropriate factor.
    */
    virtual int comm_convert_unit(std::vector<std::string> &sv, 
                                  bool itive_com);
    
    /** \brief Sort data

        For objects of type int[]:

        Sort the vector.

        (No arguments.)

        Sort the vector (in-place).

        For objects of type size_t[]:

        Sort the vector.

        (No arguments.)

        Sort the vector (in-place).

        For objects of type double[]:

        Sort the vector.

        (No arguments.)

        Sort the vector (in-place).

        For objects of type table:

        Sort the entire table by one column.

        <col> [unique]

        Sorts the entire table by the column specified in <col>.
        If the word "unique" is specified as the second argument, then
        delete duplicate rows after sorting.
     */
    virtual int comm_sort(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get object statistics

        For objects of type table:

        Show column statistics.

        <column>

        Output the average, std. dev, max and min of <column>.

        For objects of type table3d:

        Show slice statistics.

        <slice>

        Output the average, std. dev, max and min of <slice>.

        For objects of type tensor:

        Show tensor statistics.

        (No arguments.)

        The 'stats' command outputs the number of entries, their mean,
        standard deviation, minimum and maximum. It also counts the
        number of infinite or NaN values.

        For objects of type tensor_grid:

        Show tensor statistics.

        (No arguments.)

        The 'stats' command outputs the number of entries, their mean,
        standard deviation, minimum and maximum. It also counts the
        number of infinite or NaN values.
     */
    virtual int comm_stats(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get weighted statistics

        For objects of type table:

        Show weighted column statistics.

        <column> <weights>

        Output the average, std. dev, max and min of <column>, using
        the weights specified in <weights>.

        For objects of type tensor:
        
        Show stats for the data in the tensor.

        (No arguments.)

        The 'stats' command outputs the number of entries, their mean,
        standard deviation, minimum and maximum. It also counts the
        number of infinite or NaN values.
     */
    virtual int comm_wstats(std::vector<std::string> &sv, bool itive_com);

    /** \brief Print version information and O₂scl settings.
     */
    virtual int comm_version(std::vector<std::string> &sv, bool itive_com);

    /** \brief Manipulate or use a unit conversion factor

        <old unit (or "list", "add", "del", or "nat")> <new unit>
        [value to convert]

        The <tt>convert</tt> command handles unit conversions. To
        compute a unit conversion factor and then optionally apply
        than conversion factor to a user-specified value. use the form
        <tt>'acol -convert <old unit> <new unit> [value]'</tt>.
        Conversions which presume ħ=c=kB=1 are allowed by default. For
        example, <tt>'acol -convert MeV 1/fm'</tt> returns
        '1.000000e+00 MeV = 5.067731e-03 1/fm'. The conversion factor
        is output at the current value of <tt>precision</tt>, but is
        always internally stored with full double precision. 

        If no value to convert is specified, then a value of 1.0
        is assumed.

        Conversions are cached, so that if the user requests 
        an identical conversion with a different numerical value,
        then obtaining the conversion from the cache is faster than
        looking it up and processing it each time.

        The <tt>convert</tt> command attempts to handle arbitrary
        combinations of powers of SI base units to automatically
        compute new unit conversion. For example, <tt>'acol -convert
        "fm^10/g^30" "m^10/kg^30"'</tt> reports <tt>1.000000e+00
        fm^10/g^30 = 1.000000e-60 m^10/kg^30</tt>. 

        Unit conversions containing constants stored in the
        <tt>constant</tt> library are also allowed. For example,
        <tt>'acol -convert "Msun^2" "g^2"'</tt> gives <tt>1.000000e+00
        Msun^2 = 3.953774e+66 g^2</tt>. SI units are also understood,
        and both μ and "mu" are interpreted as the "micro" prefix. For
        example, <tt>'acol -convert "μm" "pc"'</tt> or <tt>'acol
        -convert "mum" "pc"'</tt> both report the conversion between
        micrometers and parsecs.

        To print the list of known units, SI prefixes, and the unit
        conversion cache, use <tt>'acol -convert list'</tt>. 

        To add a unit (only MKS is supported) the format is:

        -convert add <unit> <power of meters> <power of kg>
        <power of seconds> <power of Kelvin> <power of amps>
        <power of moles> <power of candelas> <value> <long name>

        To delete a unit, the format is:

        -convert del <unit>

        However, note that deleting a unit does not delete its
        occurences in the unit conversion cache. 

        While ħ=c=kB=1 is assumed by default, the user can disable
        conversions taking advantage of these assignments. To modify
        the use of natural units, use:

        -convert nat <boolean for c=1> <boolean for ħ=1> <boolean for kB=1>
     */
    virtual int comm_convert(std::vector<std::string> &sv, bool itive_com);

    /** \brief Copy an O₂scl-generated HDF5 file

        <source> <destination>

        Copy all O₂scl objects from one HDF5 file to another. This may
        not work for HDF5 files generated outside of O₂scl. The source
        and destination filenames may not be identical. The
        destination file may not be the same size as the source, but
        will contain the same information.
     */
    virtual int comm_h5_copy(std::vector<std::string> &sv, 
                             bool itive_com);

    /** \brief Get or modify a physical or numerical constant.

        <name, pattern, "add", "del", "list", "list-full"> [unit]

        If the constant has no units, like the Euler-Mascheroni
        constant, then e.g. <tt>acol -constant euler</tt> will report
        the value 5.772157e-1 (at the default precision which is 6).
        If the user requests a <tt>precision</tt> larger than 15
        (double precision), then the <tt>constant</tt> command fails
        and prints an error message.

        If the constant has units but no units are specified as
        arguments to the <tt>constant</tt> command, then all values of
        the constant in the library are written to the screen. If a
        unit is specified, then the <tt>constant</tt> command tries to
        find the unique value with the specified unit. The user can
        specify, <tt>mks</tt>, <tt>cgs</tt>, or the exact unit string
        of the constant. For example, <tt>'acol -constant hbar cgs'</tt>
        and <tt>'acol -constant hbar g*cm^2/s'</tt> both work and return
        the same value.

        Note that some constants in the library are not known to 
        full double precision and acol currently has no way of
        reporting this.

        Search patterns are also allowed, for example <tt>'acol
        -constant "satur*"'</tt> returns all the constants related to
        saturn in both MKS and CGS units. If <tt>use_regex</tt> is set
        to true, then regex is used to do the pattern matching, and
        otherwise <tt>fnmatch()</tt> is used. Unicode is allowed, but
        pattern matching and unicode is not fully functional.

        To list all the constants in the library, use <tt>'acol
        -constant list'</tt>. Alternatively, <tt>acol -constant
        list-full</tt> gives all information for all constants,
        including all aliases, the source, and all the decompositions
        into base units.

        One can delete a constant with, e.g. <tt>'acol -del
        pi'</tt> (this doesn't quite work yet for constants with
        different values in different unit systems).

        To add a constant, one must specify the name of the constant,
        the value, the unit system, the unit flag, the source, and
        then the 7 powers of the SI base units (in order
        m,kg,s,K,A,mol,cd).
     */
    virtual int comm_constant(std::vector<std::string> &sv, bool itive_com);
    //@}
    
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
