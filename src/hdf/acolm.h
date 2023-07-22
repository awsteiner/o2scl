/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
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
#include <o2scl/kde_python.h>
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

         - (Future) Fix documentation for value-grid command.

         - (Future) There is quite a bit of code duplication in
           comm_autocorr() between the "table" and "other" types. 
           This could be streamlined.

         - (Future) sum/max/min/output/interp/deriv/integ/deriv2 
           for hist, hist_2d, and v<c>

         - (Future) Commands xindex and yindex for table3d.

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

         - (Future) For strings, allow selection of words or substrings

      \endverbatim
  */
  class acol_manager {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// \name Typedefs for multiprecision commands
    //@{
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> >
    cpp_dec_float_100;
    
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> >
    cpp_dec_float_50;
    
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
    cpp_dec_float_35;
    
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
    cpp_dec_float_25;

#ifdef O2SCL_MPFR
    
    typedef boost::multiprecision::number<
      boost::multiprecision::mpfr_float_backend<100> > mpfr_100;
    
    typedef boost::multiprecision::number<
      boost::multiprecision::mpfr_float_backend<50> > mpfr_50;
    
    typedef boost::multiprecision::number<
      boost::multiprecision::mpfr_float_backend<35> > mpfr_35;
    
    typedef boost::multiprecision::number<
      boost::multiprecision::mpfr_float_backend<25> > mpfr_25;
    
#endif
    
    //@}
    
    /// Random number generator
    o2scl::rng<> rng;
    
    /// Random number generator for long double types
    o2scl::rng<long double> rng_ld;
    
    /// The object which sends Slack messages
    o2scl::slack_messenger smess;
    
    /** \brief A list of all type-specific commands for each type
     */
    std::map<std::string,std::vector<std::string> > type_comm_list;

    /** \brief A list of all types
     */
    std::vector<std::string> type_list;
    
    /** \brief The object for the set function
     */
    o2scl::comm_option_mfptr<acol_manager> cset;
    
    /// Document strings for commands
    std::vector<std::vector<std::string>> cmd_doc_strings;
    
    /// Document strings for parameters
    std::vector<std::vector<std::string>> param_doc_strings;
    
    /// Document strings for help topics
    std::vector<std::vector<std::string>> help_doc_strings;
    
    /// Convert units object (initialized by constructor to global object)
    o2scl::convert_units<double> &cng;

    /// The number formatter for unicode output
    o2scl::format_float ff;

    /// \name Parameters modifiable by the user
    //@{
    /// The output precision (default 6)
    int precision;

    /** \brief Output data into neat columns where possible (default
        true)
     */
    bool pretty;
    
    /** \brief If true, output names at the top of some objects
     */
    bool names_out;

    /// If true, use regex (false)
    bool use_regex;

    /// The name of the table
    std::string obj_name;
  
    /** \brief The number of columns requested by the user 
        (default 0 for autodetect)
    */
    int ncols;
    
    /** \brief The interpolation type (default is 1; linear)

        The valid interpolation types are 1: linear, 2: cubic, 3:
        periodic cubic, 4: akima, 5: periodic akima, 6: monotonic
        [experimental], 7: Steffen's monotonic, 8: nearest neighbor
        [experimental].

        Note that when this value is modified and the current object
        is a table, table3d or a histogram object, it also modifies
        the interpolation type of the object.
    */
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
    o2scl::cli::parameter_string p_color_spec;
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

    /// Default arguments from environment
    std::string def_args;

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
    std::vector<std::vector<double>> vvdouble_obj;

    o2scl::tensor<> tensor_obj;
    o2scl::tensor<int> tensor_int_obj;
    o2scl::tensor<size_t> tensor_size_t_obj;
    o2scl::tensor_grid<> tensor_grid_obj;

    o2scl::prob_dens_mdim_amr<> pdma_obj;
    o2scl::prob_dens_mdim_gaussian<> pdmg_obj;
    o2scl::prob_dens_mdim_gmm<> pgmm_obj;
    o2scl::kde_python<> pkde_obj;
    //@}
    
    /** \brief True if we should run interactive mode after parsing
        the command-line
    */
    bool post_interactive;

    /// The environment variable to read from 
    std::string env_var_name;

    /// \name Colors
    //@{
    /// Color for commands
    std::string command_color;
    /// Color for types
    std::string type_color;
    /// Color for parameters
    std::string param_color;
    /// Color for help topics
    std::string help_color;
    /// Color for executable strings
    std::string exec_color;
    /// Color for URLSs
    std::string url_color;
    /// Default color
    std::string default_color;
    //@}

    /// Default acol options for any type
    std::vector<o2scl::comm_option_s> opts_new;
    
    /** \brief Color specification for terminal output

        Note that full 3 byte colors are not supported in OSX terminal
        and not yet supported here.

        The 16 basic VT100 colors are 1-15 except for black which is
        assigned value 256. For 8 bit colors, use values greater than
        15 and smaller than 256. The thousands digit is used for the
        underline attribute, the tens of thousands digit is used for
        the intensity attribute (0 is normal, 1 is high intensity, and
        2 is low intensity). Background colors are specified similarly
        to foreground colors, but multipled by 100,000.

        By default, commands are cyan, types are magenta, parameters
        are red, help topics are green, executable commands are white,
        and URLS are underlined.
    */
    std::string color_spec;
    
  protected:

    /** \brief Clear memory associated with the current object and set
        type to ""
    */
    void clear_obj();

  public:

    /** \brief Return true if acol can provide help on the
        specified arguments

        This function is used by o2graph to determine if
        the acol help command should be called.
     */
    bool help_found(std::string arg1, std::string arg2="");
    
    /** \brief Make all of the XML replacements in string \c s
        based on the command list \c clist
    */
    void xml_replacements(std::string &s,
                          std::vector<std::string> &clist);

    /** \brief Add new commands for type \c new_type
     */
    void command_add(std::string new_type);

    /** \brief Perform the current color replacements on string \c s

        This function replaces "[c]", "[d]", etc. with the
        user-specified color strings.
     */
    void color_replacements(std::string &s);
  
    /** \brief Update the command documentation from the o2scl
        data file

        This function iterates through all of the command line
        interface options and sets their documentation from the o2scl
        file.

        This function is called in \ref command_add() to update
        the documentation for a new type and in \ref setup_options()
        to set up the initial documentation.
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

    /** \brief Call \ref run() with no arguments and with 
        \c full_process equal to false
     */
    int run_empty() {
      return run(0,0,false);
    }

    /** \brief Process and call the arguments specified in the
        string vector \c args 
     */
    void parse_vec_string(std::vector<std::string> &args) {
      std::vector<o2scl::cmd_line_arg> ca;
      cl->process_args(args,ca,0);
      cl->call_args(ca);
      return;
    }      
    
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
    /** \brief Add a vector_specification

        For objects of type table:

        Add column from a vector specification to the table.

        Arguments: <tt><vec. spec.> <column></tt>

        Add the vector specification <vec. spec.> to the table and
        store in the column named <column>.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::vector_spec()` for help
        on vector specifications.
        \endverbatim
    */
    virtual int comm_add_vec(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Compute autocorrelation coefficients

        Arguments: <tt><options> ...</tt>
        
        This command computes autocorrelation coefficients for the
        specified vector(s) and then estimates the associated 
        autocorrelation length(s). 

        There are three algorithms: a brute force algorithm (the
        default specified with the string "def"), the algorithm from
        "acor", and a method using FFTW ("fft"). The choice of
        algorithm can be specified in the options argument. Note
        that the FFT method produces a mirror image of the 
        autocorrelations on the RHS. The autocorrelation computation
        from thesethree methods is relatively good, but the estimate
        of the autocorrelation length is only approximate. 

        If multiple vectors are given, then there are two options.
        Option "max" (the default) means that autocorrelation
        coefficients are computed separately for each vector and then
        the maximum autocorrelation length is specified at the end.
        Option "avg" means that the autocorrelation coefficients are
        averaged over all vectors and a single autocorrelation length
        for the averaged data is reported at the end.

        Finally, the "store" option, if specified, means that the
        autocorrelation coefficients are stored afterwards in 
        a <tt>vec_vec_double</tt> object (and the current object,
        if present, is cleared).

        Options can be combined with commas (but no spaces), for
        example "fft,max,store" or "def,avg".

        If there is no current object, then the user may specify one
        or more multiple vector specifications. If the current object
        is of type <tt>int[]</tt>, <tt>double[]</tt>, or
        <tt>size_t[]</tt>, then the autocorrelation coefficients of
        the current vector are computed. If the current object is of
        type <tt>table</tt>, then the user may specify either columns
        of the table or a multiple vector specification as additional
        arguments.

        For example, using <tt>o2graph</tt> from o2sclpy: 

        <tt>o2graph -create table x grid:0,4000,1 -function
        "exp(x/400)*sin(floor(x/10)/5)+rand" z -autocorr store z -plot
        0 -show</tt>

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::mult_vector_spec()` for help
        on multiple vector specifications.
        \endverbatim
    */
    virtual int comm_autocorr(std::vector<std::string> &sv, bool itive_com);

    /** \brief Perform a forward or reverse FFT

        For objects of type table:

        Perform an FFT on a pair of columns

        Arguments: <tt><real part in> <complex part in>
        <real part out> <complex part out> ["backward"]</tt>
        
        Perform an FFT.

        For objects of type table3d:

        Perform an FFT on a pair of slices

        Arguments: <tt><real part in> <complex part in>
        <real part out> <complex part out> ["backward"]</tt>
        
        Perform an FFT.
    */
    virtual int comm_fft(std::vector<std::string> &sv, bool itive_com);

    /** \brief Assign a constant

        For objects of type table:

        Assign a constant to the table, e.g. <tt>-assign pi
        "acos(-1)"</tt>.

        Arguments: <tt><name> val</tt>

        Assign a constant value to a name for the present table. Valid
        constant values are things like <tt>1.618</tt>,
        <tt>acos(-1.0)</tt> or <tt>sin(4^5)</tt>. To remove an
        assignment, call assign with a blank value.

    */
    virtual int comm_assign(std::vector<std::string> &sv, bool itive_com);

    /** \brief Average rows together 

        For objects of type table:

        Average rows of some or all columns together.

        Arguments: <tt><column or '*' for all> <window> [block
        averages]</tt>

        The first argument is the column to be modified. If the first
        argument is '*', then all columns are averaged. The second
        argument is the size of the window. By default, rolling
        averages are computed and the size of the column is unchanged.
        If the third argument evaluates to true and the first argument
        is '*', then block averages instead of rolling averages are
        computed, and then the number of rows is divided by the window
        parameter. 
    */
    virtual int comm_average_rows(std::vector<std::string> &sv,
                                  bool itive_com);

    /** \brief Binary function for tensors

        For objects of type tensor_grid:

        Apply a binary function to two tensor_grid objects.
        
        Arguments: <tt><file> <object name> <function></tt>

        Read a <tt>tensor_grid</tt> named <object name> from file
        <file> and use it along with the function <function> to modify
        the current <tt>tensor_grid</tt> object. The <function>
        parameter should be a mathematical function of the value in
        the current tensor (v), the value in the tensor named <object
        name> (w), the indices (i0, i1, ...) or the grid points (x0,
        x1, ...).
    */
    virtual int comm_binary(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the value of a constant expression.

        Arguments: <tt><expr> ["1" or "true" for adaptive
        multiprecision]</tt>

        This computes the value of the constant mathematical
        expression <expr>. Examples are <tt>-calc acos(-1)</tt> or
        <tt>-calc 2+1/sqrt(2.0e4)</tt>. To see which operators and
        functions can be used, use '-help <tt>functions</tt>'.

        Results are given at the current value of <tt>precision</tt>.
        Values of precision up to 50 are allowed, and multiprecision
        (rather than double precision) arithmetic is used if
        necessary. For example, try <tt>-set precision 45 -calc
        "acos(-1)"</tt>. When adaptive multiprecision is not used, the
        calculations are relatively fast, but small errors due to the
        finite precision may give incorrect results at the requested
        precision. If the second optional argument evaluates to true,
        the calc command uses multiprecision to attempt to ensure the
        result is exact to within the requested precision. However,
        this option also makes the calculation slower by at least a
        factor of two.

        Note that adaptive multiprecision is only available for OSX at
        the moment.

        Constant values from the constant library (see e.g. '-help
        <tt>constant</tt>') will automatically be used, so long as
        they have a unique value in MKS units. However, some constant
        values are currently only stored to double precision and will
        be arbitrarily promoted to higher-precision without warning.
        Unicode is also supported for constants, so try, e.g. 
        <tt>-set precision 15 -calc π</tt>.

        Note that the variable <tt>precision</tt> is used for the
        argument to the <tt>cout.precision()</tt> function, so a
        precision of 10 is actually 11 significant figures. 
        When adaptive multiprecision is enabled, the
        value is computed to within a relative tolerance of \f$
        10^{-\mathrm{precision}-1} \f$.
    */
    virtual int comm_calc(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Concatenate two objects

        For objects of type table:

        Concatenate a second table object onto current table.

        Arguments: <tt><file> [name]</tt>

        This command adds the rows in second table to the end of the
        current table. If [name] is not provided, then the second
        table will be first object of type <tt>table</tt> in the the
        specified file. The resulting table will always have a number
        of rows equal to the sum of the rows in the original two
        tables. If the two tables have column names which are the
        same, then the data from these two columns will be
        concatenated even if the column ordering is different. If the
        second table has a column not present in the current table,
        then a new column in the current table is created. 

        Columns which are present in only one of the two tables will
        result in columns which have multiple zero entries in the new
        resulting table.

        For example:

        <tt>-create table x grid:0,5,1 -function "sin(x)" y -internal
        temp.o2 -create table x grid:0,3,1 -function "cos(x)" z -cat
        temp.o2 -output</tt>

        For objects of type table3d:

        Concatenate data from a second table3d onto current table3d.

        Arguments: <tt><file> [name]</tt>

        The <tt>cat</tt> command adds all slices from the first \c
        table3d object in the specified file which are not already
        present in the current \c table3d object. The x and y grids in
        the current table3d object are unmodified. If the x and y
        grids in the two tables are the same, then the values are
        simply copied over. If the grids are different, then the x and
        y grids are interpolated into the new table3d object (using
        it's associated interpolation type) to fill the new slice in
        the current table3d object.
    */
    virtual int comm_cat(std::vector<std::string> &sv, bool itive_com);

    /** \brief Clear the current object

        Arguments: (No arguments.)

        Deallocate the memory associated with the current object. This
        command does not clear the object name stored in \c obj_name.
    */
    virtual int comm_clear(std::vector<std::string> &sv, bool itive_com);

    /** \brief List commands, with an optional type argument

        Arguments: <tt>[type]</tt>

        If no argument is specified, list all valid commands for the
        current type (including those commands which do not require a
        current object). If a type argument is given, then list all
        valid commands for the specified type.
    */
    virtual int comm_commands(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Get or modify a physical or numerical constant.

        Arguments: <tt><name, pattern, "add", "del", "list", or 
        "list-full"> [unit] or [value unit flag source m kg s K A mol 
        cd]</tt>

        If the constant has no units, like the Euler-Mascheroni
        constant, then e.g. <tt>-constant euler</tt> will report
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
        of the constant. For example, <tt>-constant hbar cgs</tt>
        and <tt>-constant hbar g*cm^2/s</tt> both work and return
        the same value.

        Note that some constants in the library are not known to 
        full double precision and \c acol currently has no way of
        reporting this.

        Search patterns are also allowed, for example 
        <tt>-constant "satur*"</tt> returns all the constants related to
        Saturn in both MKS and CGS units. If <tt>use_regex</tt> is set
        to true, then regex is used to do the pattern matching, and
        otherwise <tt>fnmatch()</tt> is used. Unicode is allowed, but
        pattern matching and Unicode do not work well together.

        To list all the constants in the library, use 
        <tt>-constant list</tt>. Alternatively, <tt>-constant
        list-full</tt> gives all information for all constants,
        including all aliases, the source, and all the decompositions
        into base units.

        One can delete a constant with, e.g. <tt>-constant del
        pi</tt> (this doesn't quite work yet for constants with
        different values in different unit systems), but this
        deletion only lasts for the current invocation of the acol
        command. 

        To add a constant, one must specify the name of the constant,
        the value, the unit system, the unit flag, the source, and
        then the 7 powers of the SI base units (in order
        m,kg,s,K,A,mol,cd).
    */
    virtual int comm_constant(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute contour lines

        For objects of type table3d:

        Create contour lines from a table3d slice.

        Arguments: <tt><val> <slice_name> ["output file"
        "object name"]</tt>
        
        The \c contours command constructs a set of contour lines
        using the data in slice named <slice> at the fixed value given
        in <val>. If no additional arguments are given, then the \c
        table3d object is deleted from memory and the contour lines
        which were computed become the new current object of type
        <tt>vector<contour_line></tt>. If two additional arguments are
        given, then the contour lines are stored in the file named
        [output file] and the object is named object_name and the
        current \c table3d object is retained. If the file does not
        exist, it is created. If no contours are found, then a message
        is output to the screen, no file I/O is performed, and the
        current \c table3d object is unmodified.

        Countours are computed by the \ref o2scl::contours class, by
        piecing together line segments across grid lines. The best way
        to obtain more accurate contours is to compute them using data
        with a smaller grid spacing. The \c refine command can be used
        to refine the data using a simple interpolation scheme, but
        the accuracy of the contours is then limited by the accuracy
        of the interpolation.

        For objects of type hist_2d:

        Create contour lines from a table3d slice.

        Arguments: <tt>["frac" or "frac2"] <val> ["output file"
        "output name"]</tt>

        If the optional arguments "frac" or "frac2" are not present,
        the \c contours command constructs a set of contour lines
        using at the fixed value given in <val>. If no additional
        arguments are given, then the \c table3d object is deleted
        from memory and the contour lines which were computed become
        the new current object of type <tt>vector<contour_line></tt>
        If two additional arguments are given, then the contour lines
        are stored in the file named [output file], the object is
        named [object name], and the current \c hist_2d object
        retained. If the file does not exist, it is created. If no
        contours are found, then a message is output to the screen, no
        file I/O is performed, and the current table3d object is
        unmodified.

        If the argument "frac" is present, then the operation of the
        \c contours command is the same except that <val> is
        interpreted as a fraction of the total integral under the
        data. The integral is computed as the sum of (w-min) Δx Δy .
        where "min" is the minimum weight over all bins. This
        definition ensures that the contours exist for all values
        between 0 and 1.

        If the argument "frac2" is present, then the <val> is
        interpreted as a fraction of the total integral, but the
        integral is computed as the sum of w Δx Δy and thus there are
        no guarantees that the requested contours exist.

        Countours are computed by the \ref o2scl::contours class, by
        piecing together line segments across grid lines. The best way
        to obtain more accurate contours is to compute them using data
        with a smaller grid spacing. The \c refine command can be used
        to refine the data using a simple interpolation scheme, but
        the accuracy of the contours is then limited by the accuracy
        of the interpolation.
    */
    virtual int comm_contours(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Manipulate or use a unit conversion factor

        Arguments: <tt><old unit (or "list", "add", "del", or "nat")> 
        <new unit> [value to convert]</tt>

        The <tt>convert</tt> command handles unit conversions. To
        compute a unit conversion factor and then optionally apply
        than conversion factor to a user-specified value. use the form
        <tt>-convert <old unit> <new unit> [value]</tt>.
        Conversions which presume ħ=c=kB=1 are allowed by default. For
        example, <tt>-convert MeV 1/fm</tt> returns
        <tt>1.000000e+00 MeV = 5.067731e-03 1/fm</tt>. The conversion
        factor is output at the current value of <tt>precision</tt>,
        but is always internally stored with full double precision.

        If no value to convert is specified, then a value of 1.0
        is assumed.

        Conversions are cached, so that if the user requests 
        an identical conversion with a different numerical value,
        then obtaining the conversion from the cache is faster than
        looking it up and processing it each time.

        The <tt>convert</tt> command attempts to handle arbitrary
        combinations of powers of SI base units to automatically
        compute new unit conversion. For example, <tt>-convert
        "fm^10/g^30" "m^10/kg^30"</tt> reports <tt>1.000000e+00
        fm^10/g^30 = 1.000000e-60 m^10/kg^30</tt>. 

        Unit conversions containing constants stored in the
        <tt>constant</tt> library are also allowed. For example,
        <tt>-convert "Msun^2" "g^2"</tt> gives <tt>1.000000e+00
        Msun^2 = 3.953774e+66 g^2</tt>. SI units are also understood,
        and both μ and "mu" are interpreted as the "micro" prefix. For
        example, <tt>-convert "μm" "pc"</tt> or 
        <tt>-convert "mum" "pc"</tt> both report the conversion between
        micrometers and parsecs.

        To print the list of known units, SI prefixes, and the unit
        conversion cache, use <tt>acol -convert list</tt>. 

        To add a unit (only MKS is supported) the format is:

        <tt>-convert add <unit> <power of meters> <power of kg>
        <power of seconds> <power of Kelvin> <power of amps>
        <power of moles> <power of candelas> <val> <long name></tt>

        To delete a unit, the format is:

        <tt>-convert del <unit></tt>

        However, note that deleting a unit does not delete its
        occurrences in the unit conversion cache. 

        While ħ=c=kB=1 is assumed by default, the user can disable
        conversions taking advantage of these assignments. To modify
        the use of natural units, use e.g. <tt>-convert nat 0 1 0</tt>.
        The second argument is 1 for c=1 and 0 otherwise, the third
        argument is 1 for ħ=1 and 0 otherwise, and finally the
        fourth argument is 1 for kB=1 and 0 otherwise.
    */
    virtual int comm_convert(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert units

        For objects of type table:

        Convert a column to a new unit.

        Arguments: <tt><column> <new unit></tt>

        Convert the units of column <column> to <new unit>,
        multipliying all entries in that column by the appropriate
        factor. See <tt>acol -convert list</tt> for a list of units
        and prefixes which can be automatically converted.
    */
    virtual int comm_convert_unit(std::vector<std::string> &sv, 
                                  bool itive_com);
    
    /** \brief Compute correlation coefficients.

        For objects of type table:

        Compute the correlation coefficient between two columns.

        Arguments: <tt>[column 1, column 2]</tt>

        Compute the correlation coefficient between two columns, or,
        if no arguments are given, then compute the correlation
        coefficients between all pairs of columns. Results are output
        to the screen.
    */
    virtual int comm_correl(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create an object.

        Arguments: <tt><type> [...]</tt>

        Create a new object of type <type>. If an object is currently
        in memory, it is deallocated before creating the new object.

        <tt>create <type> <val></tt>: For types <tt>char</tt>,
        <tt>int</tt>, <tt>size_t</tt>, and <tt>string</tt>, create an
        object and give it the initial value specified.

        <tt>create double <value spec.></tt>: Create a <tt>double</tt>
        object and set it equal to the value specified by <value
        spec.>.

        <tt>create <type> <size> <function of "i"></tt>: For array
        types <tt>int[]</tt> and <tt>size_t[]</tt>, the user must
        specify the size of the array and a function of the array
        index <tt>i</tt> to fill the array.

        <tt>create double[] [<size> <function of "i">]
        or [vector spec.]</tt>: For <tt>double[]</tt> the user must either
        give a vector specification, or specify the size of the array
        and a function of the array index <tt>i</tt>.

        <tt>create table <name> <vector spec.></tt>:
        Create a new <tt>table</tt> object with one column named <name>
        from a vector specification.
        
        <tt>create table-mv <mult. string spec.> <mult. vector spec.></tt>:
        Create a new <tt>table</tt> object with several columns with
        names taken from the multiple string specification and 
        data taken from the multiple vector specification.
        
        <tt>create tensor <rank> <size 0> <size 1> ...</tt>: Create a
        <tt>tensor</tt> object with the specified rank and sizes. All
        tensor entries are initialized to zero.

        <tt>create tensor_grid <rank> <size 0> <size 1> ...</tt>:
        Create a <tt>tensor_grid</tt> object with the specified rank
        and sizes. The tensor grid is initialized to count each index
        (beginning with zero) and the entries of the tensor are
        initialized to zero. The grid can be specified afterwards
        using <tt>set-grid</tt>.

        <tt>create table3d <x name> <x vector spec.> <y name> <y
        vector spec.> <slice name> <slice func.></tt>: Create a
        new <tt>table3d</tt> object which has one slice. The x and y
        grids are given as vector specifications (see "acol -help
        vector-spec" for the syntax). The slice function can be
        written in terms of the x- and y-grid values which are
        referred to by name.

        For example, using <tt>o2graph</tt> from o2sclpy:

        <tt>o2graph -create table3d x func:100:i/200 y func:100:i/200
        z "sin(1/(x+0.01))* sin(1/(y+0.01))" -den-plot z -xtitle x
        -ytitle y -show</tt>

        <tt>create vec_vec_double <mult. vector spec.></tt>: Create a
        <tt>vec_vec_double</tt> object using the given multiple 
        vector specification.
        
        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::value_spec()` for help on value
        specifications, :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications, and 
        :cpp:func:`o2scl_hdf::mult_vector_spec()` for help
        on multiple vector specifications.

        \endverbatim
    */
    virtual int comm_create(std::vector<std::string> &sv, bool itive_com);

    /** \brief Delete a column

        For objects of type table:

        Delete a table column.

        Arguments: <tt><name></tt>

        Delete the entire column named <name>.
    */
    virtual int comm_delete_col(std::vector<std::string> &sv, bool itive_com);

    /** \brief Delete rows

        For objects of type table:

        Delete rows selected by a function.

        Arguments: <tt><function></tt>

        Delete the set of rows for which a function evaluates to a
        number greater than 0.5. For example, <tt>-delete-rows
        if(col1+col2>10,1,0)</tt> will delete all columns where the
        sum of the entries in col1 and col2 is larger than 10. See
        also the \c select-rows command.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim
    */
    virtual int comm_delete_rows(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Delete rows which match to within a specified tolerance

        For objects of type table:

        Delete rows which match to within a specified tolerance

        Arguments: <tt>[relative tol.] [absolute tol.]</tt>

        for every column, the entries in the two rows do not match if
        either of their absolute values are greater the absolute
        tolerance and their relative deviation is greater than the
        relative tolerance. if all columns match, then the rows match.

        The \c delete-rows-tol command deletes all rows which
        match a row earlier in the table. if verbose is larger than
        zero then information about how many rows were deleted is
        provided. 
    */
    virtual int comm_delete_rows_tol(std::vector<std::string> &sv,
                                     bool itive_com);

    /** \brief Compute a derivative

        For objects of type table:

        Derivative of a function defined by two columns.

        Arguments: <tt><x> <y> <name></tt>

        Create a new column named <name> filled with the derivative
        of the function \f$ y(x) \f$ obtained from columns <x> and <y>.

        For objects of type tensor:

        Compute the derivative of the tensor object w.r.t. an index

        Arguments: <tt><index></tt>

        The <tt>deriv</tt> command differentiates the tensor object
        with respect to one of the indices.

        For objects of type tensor_grid:

        Compute the derivative of the tensor object w.r.t. an index

        Arguments: <tt><index></tt>

        The <tt>deriv</tt> command differentiates the tensor object
        with respect to one of the indices.

        For objects of type double[]:

        Replace the array with its derivative.

        Arguments: (No arguments.)

        Replace the array with its derivative using the current
        interpolation type.

        For objects of type int[]:

        Replace the array with its derivative.

        Arguments: (No arguments.)

        Replace the array with its derivative using the current
        interpolation type, converting it to a double[].

        For objects of type size_t[]:

        Replace the array with its derivative.

        Arguments: (No arguments.)

        Replace the array with its derivative using the current
        interpolation type, converting it to a double[].
    */
    virtual int comm_deriv(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a slice which is the derivative wrt x of another

        For objects of type table3d:

        Derivative with respect to x.

        Arguments: <tt><f> <dfdx></tt>

        Create a new slice named <dfdx> filled with the derivative of
        the function from the x grid and slice named <f>.
    */
    virtual int comm_deriv_x(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a slice which is the derivative wrt y of another

        For objects of type table3d:

        Derivative with respect to y.

        Arguments: <tt><f> <dfdy></tt>

        Create a new slice named <dfdy> filled with the derivative of
        the function from the y grid and slice named <f>.
    */
    virtual int comm_deriv_y(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute a second derivative

        For objects of type table:

        Second derivative of a function defined by two columns.

        Arguments: <tt><x> <y> <name></tt>

        Create a new column named <name> filled with the second
        derivative of the function y(x) obtained from columns <x> and
        <y>.
    */
    virtual int comm_deriv2(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get entries along the main diagonal

        For objects of type tensor:

        Get diagonal elements.

        Arguments: (No arguments.)

        Extract only the elements on the main diagonal to create a
        <tt>double[]</tt> object.
    */
    virtual int comm_diag(std::vector<std::string> &sv, bool itive_com);

    /** \brief Open local HTML docs for O₂scl documentation.

        Arguments: <tt>[section, class, or function]</tt>

        If [topic] is unspecified, this command opens up the local
        HTML documentation for O₂scl in the default web browser using
        'open' on OSX and 'xdg-open' on other systems. If a topic is
        specified, then the <tt>docs</tt> command looks for the O₂scl
        local HTML documentation for the specified section, class, or
        function. If the documentation is found, then that web page is
        opened instead. 

        The operation of the <tt>docs</tt> command depends on the
        assumption that the HTML documentation files have not been
        moved after the installation process. In order to open the
        remote version of the documentation instead of the local copy,
        use the <tt>wdocs</tt> command instead.
    */
    virtual int comm_docs(std::vector<std::string> &sv, bool itive_com);

    /** \brief Download a file from the specified URL.

        Arguments: <tt><file> <URL> [hash, \"file:\"hash_filename, or
        \"none\"] [directory]</tt>

        First, look for the file named <file> (in directory given
        [directory] if specified). If a hash is not specified and the
        file exists, then return success. If a hash is specified, then
        compare the file with the specified hash. If they match, then
        return success. If the file doesn't exist or doesn't match
        the specified hash, then use <tt>curl</tt> to
        download the file. Again compare with the hash if specified.
        If the download files or the file doesn't match the hatch,
        then call the error handler. 

        If the filename is "_", then the file is extracted from the
        end of the URL. 

        This function exits immediately if it fails, preventing the
        user from reading a data file which is corrupted.

        This function uses \ref o2scl::cloud_file to handle the 
        file acquisition.
    */
    virtual int comm_download(std::vector<std::string> &sv, bool itive_com);

    /** \brief List objects in a HDF5 file

        Arguments: <tt><file> [group]</tt>

        This lists all the top-level datasets and groups in a HDF5
        file and, for those groups which are in the O₂scl format,
        gives the type and name of the object stored in that HDF5
        group.

        If a group is specified, then all of the top-level datasets
        and groups inside that specified group are listed.
    */
    virtual int comm_filelist(std::vector<std::string> &sv, bool itive_com);

    /** \brief Find a value in an object

        For objects of type double[]:

        Find a value in the array

        Arguments: <tt><value></tt>

        Find the closest value to <value> in the array and print 
        out the associated index.
        
        For objects of type int[]:

        Find a value in the array

        Arguments: <tt><value></tt>

        Find the closest value to <value> in the array and print 
        out the associated index.
        
        For objects of type size_t[]:

        Find a value in the array

        Arguments: <tt><value></tt>

        Find the closest value to <value> in the array and print 
        out the associated index.
        
     */
    virtual int comm_find(std::vector<std::string> &sv, bool itive_com);

    /** \brief Find a row

        For objects of type table:

        Find a row which has a specified value or  maximizes a function.

        Arguments: <tt><func> or find-row <col> <val></tt>

        If one argument is given, then <tt>find-row</tt> finds the row
        which maximizes the value of the expression given in <func>,
        and then output the entire row. Otherwise, <tt>find-row</tt>
        finds the row for which the value in column named <col> is as
        close as possible to the value <val>. See command
        <tt>get-row</tt> to get a row by its index.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim
    */
    virtual int comm_find_row(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Fit data to a function

        For objects of type table:

        Fit two columns to a function (experimental).

        Arguments: <tt><x> <y> <yerr> <ynew> <par names> <func> <vals></tt>

        Detailed desc.
    */
    virtual int comm_fit(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a column from a function

        For objects of type table:

        Create a column from a function

        Arguments: <tt><func> <name></tt>

        Set the column named <name> to the result of a function,
        <func>, in terms of the other columns. If the column does not
        already exist, a new one is added to the table. For example,
        for a table containing columns named 'c1' and 'c2', <tt>function
        c1-c2 c3</tt> would create a new column c3 which contains the
        difference of columns 'c1' and 'c2'. 

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type double[]:

        Set the values of the array given a function.

        Arguments: <tt><function></tt>

        Set the values of the array given a user-specified function of
        'i'. For example, <tt>(sin(i)>1)*4</tt>.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type int[]:

        Set the values of the array given a function.

        Arguments: <tt><function></tt>

        Set the values of the array given a user-specified function of
        'i'. For example, <tt>(sin(i)>1)*4</tt>.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type size_t[]:

        Set the values of the array given a function.

        Arguments: <tt><function></tt>

        Set the values of the array given a user-specified function of
        'i'. For example, <tt>(sin(i)>1)*4</tt>.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type hist:

        Apply a function to the weights.
        
        Arguments: <tt><function></tt>

        Apply a function to the weights.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type table3d:

        Create a new slice from a function.

        Arguments: <tt><func> <name></tt>

        Set the slice named <name> to the result of a function,
        <func>, in terms of the other slices. If the slice does not
        already exist, a new one is created. For example, for a
        table3d containing slices named 's1' and 's2', 'function s1-s2
        s3' would create a new column 's3' which contains the
        difference of columns 's1' and 's2'.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type tensor:

        Set the tensor values from a function.

        Arguments: <tt>[cond. function] <function of v, i0, i1, ...></tt>

        The \c function command sets all entries in a tensor equal to
        a user-specified mathematical function of the indices. When
        the conditional function evaluates to a number less than or
        equal to 0.5, then the tensor entry will be unchanged. (For
        more help with functions, type acol -help functions)

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type tensor_grid:

        Set the tensor values from a function.

        Arguments: <tt>[conditional func.] <func. of v, i0, i1, ...
        and x0, x1, ...></tt>

        The \c function command sets all the data entries in a
        tensor_grid equal to a user-specified mathematical function of
        the value in the tensor (v), the indices (i0, i1, ...) or the
        grid points (x0, x1, ...). If two function arguments are given
        and if the first function argument is not "none", then the
        first function specifies which tensor entries are to be
        modified. When the conditional function evaluates to a number
        less than or equal to 0.5, then the tensor entry will be "
        unchanged. (For more help with functions, type "acol -help
        functions.".)

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim
    */
    virtual int comm_function(std::vector<std::string> &sv, bool itive_com);

    /** \brief Read an object from a generic text file
        
        Arguments: <tt><type> <file or "cin"></tt>
        
        Read an object of type <type> from a text file named <file>,
        or from <tt>std::cin</tt>.

        For <tt>int</tt>, <tt>char</tt>, <tt>double</tt>, or
        <tt>size_t</tt> objects, the input is assumed to begin
        with the object and it is read using <tt>operator>>()</tt>.

        For <tt>int[]</tt>, <tt>double[]</tt>, or
        <tt>size_t[]</tt> objects, the input file is assumed to begin
        with the desired object and it is read using operator>>().
        The array is presumed to end at the end of the file or the
        end of the input. 

        To read a numeric array stored in a single line,
        use types <tt>int[]-line</tt>, <tt>double[]-line</tt>, or
        <tt>size_t[]-line</tt>. The array is assumed to end when it 
        reaches the end of the first line. 

        If you want to read a column of numbers where the first number
        is the number of entries in the vector and subsequent numbers
        contain data (whitespace, including carriage returns is
        ignored), then use types <tt>int[]-n</tt>,
        <tt>double[]-n</tt>, or <tt>size_t[]-n</tt>.

        For <tt>string</tt> or <tt>string[]</tt> objects, the strings
        are read using <tt>std::getline()</tt>. The strings can contain
        whitespace, and they are presumed to end at the end of the line
        (i.e. at the carriage return).

        To read a list of strings stored in a single line,
        use type <tt>string[]-line</tt>. The strings cannot contain
        whitespace because the whitespace is presumed to separate 
        individual strings. The array is assumed to end when it 
        reaches the end of the first line. 

        To read a list of strings where the first input is the number
        of strings, then use type <tt>string[]-n</tt>. The strings can
        contain whitespace, and they are presumed to end at the end of
        the line (i.e. at the carriage return).

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

        Objects of type <tt>prob_dens_mdim_gaussians</tt> and
        <tt>prob_dens_mdim_gmm</tt> read a generic format which is
        defined in the O₂scl documentation for the read_generic()
        function for these classes.
    */
    virtual int comm_generic(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get the grid 

        For objects of type table3d:

        Print out the table3d grid.

        Arguments: (No arguments.)

        Output the table3d grid as a series of columns.

        For objects of type tensor_grid:

        Get the tensor grid.

        Arguments: (No arguments.)

        Output the tensor grid as a series of columns.
    */
    virtual int comm_get_grid(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get a table row by index

        For objects of type table:

        Get a table row by index.

        Arguments: <tt><index></tt>
        
        Get a row by index. The first row has index 0, and the last
        row has index n-1, where n is the total number of rows (which
        can be determined by the \c nlines or \c list commands). The
        <tt>index</tt> command creates a column of row indexes (which
        is called 'N' by default). To find a row which contains a
        particular value or maximizes a specified function, use
        <tt>find-row</tt>.
    */
    virtual int comm_get_row(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get units of a column

        For objects of type table:

        Get the units for a specified column.

        Arguments: <tt><column></tt>

        Obtains the units for the specified column. 
    */
    virtual int comm_get_unit(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Output help

        Arguments: <tt>[command or parameter or type or topic]</tt>

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
    
    /** \brief Copy an O₂scl-generated HDF5 file

        Arguments: <tt><source> <destination></tt>

        Copy all O₂scl objects from one HDF5 file to another. This
        command may not work for HDF5 files generated outside of
        O₂scl. The source and destination filenames may not be
        identical. The destination file may not be the same size as
        the source, but will contain the same information.
    */
    virtual int comm_h5_copy(std::vector<std::string> &sv, 
                             bool itive_com);

    /** \brief Add a column for line numbers

        For objects of type table:

        Add a column containing the row numbers.

        Arguments: <tt>[column name]</tt>
        
        Define a new column named [column name] and fill the column
        with the row indexes, beginning with zero. If no argument is
        given, the new column is named 'N'. If a column named 'N' is
        already present, then \c index adds underscores (up to a
        maximum of 10) in an attempt to create a new unique column
        name.
    */
    virtual int comm_index(std::vector<std::string> &sv, bool itive_com);

    /** \brief Insert a column from an external table using 
        interpolation

        For objects of type table:

        Interpolate a column from another file.

        Arguments: <tt><file> <table> <oldx> <oldy> <newx> [newy]</tt>

        Insert a column from file <fname>, interpolating it into the
        current table. The column <oldy> is the column in the file
        which is to be inserted into the table, using the column
        <oldx> in the file and <newx> in the table. The new column in
        the table is named <oldy>, or it is named [newy] if the
        additional argument is given. Note that extrapolation is
        allowed, and the operation of the \c insert command is not
        well-defined if some of the column entries are not finite.

        For objects of type table3d:

        Interpolate a slice from another file.

        Arguments: <tt><file> <table> <old> [new]</tt>

        Insert a slice from file <fname>, interpolating it into the
        current table. The slice <old> is the slice in the file which
        is to be inserted into the \c table3d object, using the slice
        <old> in the file and <new> in the current \c table3d object.
        The new slice in the table is named <old>, or [new] if the
        additional argument is given. Note that extrapolation is
        allowed, and the operation of the \c insert command is not
        well-defined if some of the slice entries are not finite.
    */
    virtual int comm_insert(std::vector<std::string> &sv, bool itive_com);

    /** \brief Insert an external table using interpolation

        For objects of type table:

        Insert a table from another file.

        Arguments: <tt><fname> [table name] [old_x new_x]</tt>

        Insert all columns from file <fname> into the current table.
        The first table is used or the table object named table_name,
        if specified. If index columns [old_x] and [new_x] are not
        specified, then the insert requires both the current and the
        source table to have the same number of rows. If they are
        specified, then interpolation using those index columns is
        used. If columns in the new table are not present in the
        current table, then they are added automatically. If a column
        in the current table has the same name as one in the new table
        then it is rewritten with new data, with one exception. If a
        column in the new table has the same name as [old_x], then it is
        left unmodified.
    */
    virtual int comm_insert_full(std::vector<std::string> &sv, bool itive_com);

    /** \brief Create a column which is the integral of another

        For objects of type table:

        Integrate a function specified by two columns.

        Arguments: <tt><x> <y> <name></tt>
        
        Create a new column named <name> filled with the integral of
        the function y(x) obtained from columns <x> and <y>. The lower
        limit of the integral is the value of the column <x> in the
        first row and the upper limit of the integral for each
        specified row is the value of the column <x> at that row.
    */
    virtual int comm_integ(std::vector<std::string> &sv, bool itive_com);

    /** \brief Toggle interactive mode

        Arguments: (No arguments.)

        If given as a command-line parameter, the \c interactive
        command toggles the execution of the interactive mode after
        the command-line parameters are processed. If zero arguments
        are given to \c acol on the command-line then the interactive
        interface is automatically turned on.
    */
    virtual int comm_interactive(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output current object in the internal HDF5 format.
        
        Arguments: <tt><file></tt>

        Output the current object to the specified file in the
        internal HDF5 format. If there is an object of the same
        name in the HDF5 file of the same type, then that object
        is overwritten. If there is an object of the same name
        in the HDF5 file with a different time, then the file
        write will fail and \c acol will exit.
    */
    virtual int comm_internal(std::vector<std::string> &sv, bool itive_com);

    /** \brief Perform an interpolation using the current object
        
        For objects of type table:

        Interpolate a number into a column.

        Arguments: <tt><x name> <x value> <y name></tt>

        Interpolate <x value> from column named <x name> into column
        named <y name> using the current interpolation type 
        as specified in \c interp_type .

        For objects of type double[]:

        Interpolate an index into the array using the current
        interpolation type as specified in \c interp_type .

        Arguments: <tt><x value></tt>

        Interpolate <x value> in the array using the current
        interpolation type as specified in \c interp_type .

        For objects of type int[]:

        Interpolate an index into the array using the current
        interpolation type as specified in \c interp_type .

        Arguments: <tt><x value></tt>

        Interpolate <x value> in the array using the current
        interpolation type as specified in \c interp_type and print
        out the result as a double.

        For objects of type size_t[]:

        Interpolate an index into the array

        Arguments: <tt><x value></tt>

        Interpolate <x value> in the array using the current
        interpolation type as specified in \c interp_type and print
        out the result as a double.

        For objects of type table3d:

        Interpolate x and y values into a slice.

        Arguments: <tt><z name> <x value> <y value></tt>

        Interpolate (<x value>,<y value>) into the slice named <z
        name> using the current interpolation type 
        as specified in \c interp_type .

        For objects of type tensor_grid:

        Linearly interpolate in the grid.

        Arguments: <tt><value 1> <value 2> <value 3> ...</tt>

        When the current object is a \c tensor_grid object, the \c
        interp uses linear interpolation (independent of the value of
        \c interp_type) to interpolate an array with size equal to the
        tensor rank into the tensor grid and outputs the result.
    */
    virtual int comm_interp(std::vector<std::string> &sv, bool itive_com);

    /** \brief Interpolate into a table3d object
        
        For objects of type table:

        Interpolate the values to create a table3d object

        Arguments: <tt><method> <options> 
        <x column> <y column> <x grid or "auto">
        <y grid or "auto"> <z column 1> <z column 2> ...</tt>
    */
    virtual int comm_interp_table3d(std::vector<std::string> &sv,
                                    bool itive_com);

    /** \brief List properties of an object

        For objects of type table:

        List the constants, column names and other info.

        Arguments: (No arguments.)

        List the constants, column names and other info.

        For objects of type table3d:

        List the slice names and print out grid info.

        Arguments: (No arguments.)

        List the slice names and print out grid info.

        For objects of type hist:

        List the ...

        Arguments: (No arguments.)

        List the ...

        For objects of type tensor:

        List the tensor rank and index sizes.

        Arguments: (No arguments.)

        List the tensor rank and index sizes.

        For objects of type vector<contour_line>:

        List the contour line levels and sizes.

        Arguments: (No arguments.)

        List the contour line levels and sizes.

        For objects of type vector<vector<double>>:

        List the number of entries and their respective sizes.

        Arguments: (No arguments.)

        List the number of entries and their respective sizes.

        For objects of type vector<vector<string>>:

        List the number of entries and their respective sizes.

        Arguments: (No arguments.)

        List the number of entries and their respective sizes.

        For objects of type tensor<int>:

        List the tensor rank and index sizes.

        Arguments: (No arguments.)

        List the tensor rank and index sizes.

        For objects of type tensor<size_t>:

        List the tensor rank and index sizes.

        Arguments: (No arguments.)

        List the tensor rank and index sizes.

        For objects of type tensor_grid:

        List the tensor rank and index sizes.

        Arguments: (No arguments.)

        List the tensor rank and index sizes.

        For objects of type hist_2d:

        List the bin edges.

        Arguments: (No arguments.)

        For objects of type hist_2d:

        List the bin edges.
    */
    virtual int comm_list(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the maximum value of a column

        For objects of type table:

        Compute the maximum value of a column.

        Arguments: <tt><column name></tt>

        Compute the maximum value of a column.

        For objects of type double[]:

        Compute the maximum value and the associated index.

        Arguments: (No arguments.)

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

        Arguments: <tt><slice name></tt>

        Compute the maximum value of <slice name>.

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

        Arguments: <tt><column name></tt>

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

        Arguments: <tt><slice name></tt>

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

    /** \brief Numerically integrate a user-specified function 

        Arguments: <tt><function> <variable> <lower limit> 
        <upper limit> [multip=false,method=kb]</tt>

        This command numerically integrates <function> with respect to
        <variable> from <lower limit> to <upper limit>. 

        The fifth argument is a set of keyword arguments. If multip is
        set to  either \c "1" or \c "true", then multiprecision is
        used to attempt to ensure the result is accurate to within the
        requested precision (multiprecision only works on OSX right now
        possibly because the Ubuntu release of boost is still behind). 

        There are three methods, kb (Kronrod from boost), deb 
        (double exponential from boost) or ac (adaptive integration
        based on CERNLIB). 

        Infinite upper or lower limits are supported (for all three
        methods), use "-infty" or "infty", respectively.

        Note that the variable \c precision is used for the argument to
        the <tt>cout.precision()</tt> function, so precision of 10 is
        actually 11 significant figures. Thus in multiprecision mode,
        the integral is computed to within a relative tolerance of \f$
        10^{-11} \f$.

        Here is an example demonstrating a dilogarithm ladder. First,
        the exact result, computed using the \c calc command.

        <tt>acol -set verbose 2 -set precision 30 -calc
        "pi^2/10-(log((sqrt(5)-1)/2)^2)" 1</tt>

        <tt>Result (cpp_dec_float_35): 7.553956195317414693865200287561e-01
        </tt>

        In this example, the \c calc command begins by
        obtaining the value of π to 35 digits. It then starts with
        35-digit precision and then compares that result to that
        obtained with 50-digit precision and finds that those two are
        equal to within the requested precision. Now, using 
        the \c ninteg command to achieve the same result:

        <tt>acol -set verbose 2 -set precision 30 -ninteg
        "(-log(1-t)/t)" t 0 "(sqrt(5)-1)/2" multip=true</tt>

        <tt>Result (cpp_dec_float_35): 7.553956195317414693865200287561e-01
        </tt>

        The <tt>ninteg</tt> command computes this same value using
        numerical integration (and obtains the same result), using
        42-digit precision internally to evaluate the integrand.
    */
    virtual int comm_ninteg(std::vector<std::string> &sv, bool itive_com);

    /** \brief Add 'nlines' as a constant to a \ref o2scl::table object

        For objects of type table:

        Add 'nlines' as a constant to a <tt>table</tt> object.

        Arguments: (No arguments.)

        Add a constant called 'nlines' to the table and set it equal
        to the number of lines (rows) in the table. The number of
        lines is also output to the screen.
    */
    virtual int comm_nlines(std::vector<std::string> &sv, bool itive_com);

    /** \brief Output the current object to screen or text file.

        Arguments: <tt>[file]</tt>

        Output the object to the screen, or if the [file] argument is
        specified, to a file. This is (supposed to be) the same format
        as can be read using the <tt>generic</tt> command (independent
        of whether or not the pretty flag is true or false).
    */
    virtual int comm_output(std::vector<std::string> &sv, bool itive_com);

    /** \brief Read an object from an O₂scl-style HDF5 file.

        Arguments: <tt><file> [object name]</tt>

        Read an HDF5 file with the specified filename. If the [object
        name] argument is specified, then read the object with the
        specified name. Otherwise, look for the first <tt>table</tt>
        object, and if not found, look for the first <tt>table3d</tt>
        object, and so on, attempting to find a readable O₂scl object.
    */
    virtual int comm_read(std::vector<std::string> &sv, bool itive_com);

    /** \brief Rearrange a tensor

        For objects of type tensor:

        Rearrange the tensor

        Arguments: <tt><index spec. 1> [index spec. 2] ...</tt>

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". 

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::index_spec()` for help
        on index specifications.
        \endverbatim

        For objects of type tensor<int>:

        Rearrange the tensor.

        Arguments: <tt><index spec. 1> [index spec. 2] ...</tt>

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". 

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::index_spec()` for help
        on index specifications.
        \endverbatim

        For objects of type tensor<size_t>:

        Rearrange the tensor.

        Arguments: <tt><index spec. 1> [index spec. 2] ...</tt>

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), and range(ix,start,end). Index
        specifications may be specified as separate arguments e.g.
        "index(1)" "fixed(2,10)" or multiple index specifications may
        be given in a single argument separated by spaces or commas,
        e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". 

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::index_spec()` for help
        on index specifications.
        \endverbatim

        For objects of type tensor_grid:

        Rearrange the tensor_grid object.

        Arguments: <tt><index spec. 1> [index spec. 2] ...</tt>

        Index specifications are: index(ix), fixed(ix), sum(ix),
        trace(ix1,ix2), reverse(ix), range(ix,start,end),
        interp(ix,value), grid(ix,begin,end,n_bins,log), and
        gridw(ix,begin,end,bin_width,log). Index specifications may be
        specified as separate arguments e.g. "index(1)" "fixed(2,10)"
        or multiple index specifications may be given in a single
        argument separated by spaces or commas, e.g. "index(1)
        fixed(2,10)" or "index(1),fixed(2,10)". 

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::index_spec()` for help
        on index specifications.
        \endverbatim
    */
    virtual int comm_rearrange(std::vector<std::string> &sv, bool itive_com);

    /** \brief Refine an object

        For objects of type table:

        Refine the table.

        Arguments: <tt><index column> <factor></tt>

        Refine the data by interpolating a fixed number of rows in
        between every row of the original table. The type of
        interpolation is determined by the value of
        <tt>interp_type</tt>. If the initial number of rows is N, then
        the final number of rows is 1+(N-1)*<factor>.

        For objects of type table3d:

        Refine the table3d.

        Arguments: <tt><factor> [log mode]</tt>
        
        Refine the data by interpolating. The type of interpolation is
        determined by the value of <tt>interp_type</tt>. If the
        initial number of rows is N, then the final number of rows is
        1+(N-1)*<factor>. If [log mode] is unspecified or "auto",
        then O₂scl attempts to automatically determine if the 
        x or y grids are more linear or logarithmic. Other options
        for [log mode] include "none" (for linear refinement),
        "x" (logarithmic in x and linear in y), "y", or "xy".

        For objects of type hist_2d:

        Refine the hist_2d.

        Arguments: <tt><factor> [log mode]</tt>
        
        Refine the data by interpolating. The type of interpolation is
        determined by the value of <tt>interp_type</tt>. If the
        initial number of rows is N, then the final number of rows is
        1+(N-1)*<factor>. If [log mode] is unspecified or "auto",
        then O₂scl attempts to automatically determine if the 
        x or y grids are more linear or logarithmic. Other options
        for [log mode] include "none" (for linear refinement),
        "x" (logarithmic in x and linear in y), "y", or "xy".
    */
    virtual int comm_refine(std::vector<std::string> &sv, bool itive_com);

    /** \brief Rename a part of an object

        For objects of type table:

        Rename a column.

        Arguments: <tt><old> <new></tt>

        Rename a column from <old> to <new>. Note that to rename the
        entire object, you should use <tt>-set obj_name new_name</tt>.

        For objects of type table3d:

        Rename a slice.

        Arguments: <tt><old> <new></tt>

        Rename a slice from <old> to <new>. Note that to rename the
        entire object, you should use <tt>-set obj_name new_name</tt>.
    */
    virtual int comm_rename(std::vector<std::string> &sv, bool itive_com);

    /** \brief Resize an object

        For objects of type double[]:

        Resize the vector.

        Arguments: <tt><new size></tt>

        Resize the vector. If the new size is larger than the old
        size, then the old entries are retained and the new entries
        are set to zero. If the old size is larger, then the first new
        size entries are left unchanged.

        For objects of type int[]:

        Resize the vector.

        Arguments: <tt><new size></tt>

        Resize the vector. If the new size is larger than the old
        size, then the old entries are retained and the new entries
        are set to zero. If the old size is larger, then the first new
        size entries are left unchanged.

        For objects of type size_t[]:

        Resize the vector.

        Arguments: <tt><new size></tt>

        Resize the vector. If the new size is larger than the old
        size, then the old entries are retained and the new entries
        are set to zero. If the old size is larger, then the first new
        size entries are left unchanged.

        For objects of type string[]:

        Resize the vector.

        Arguments: <tt><new size></tt>

        Resize the vector. If the new size is larger than the old
        size, then the old entries are retained and the new entries
        are set to zero. If the old size is larger, then the first new
        size entries are left unchanged.
    */
    virtual int comm_resize(std::vector<std::string> &sv, bool itive_com);

    /** \brief Preview the current object

        Arguments: <tt>[number of lines] [number of columns]</tt>

        Print out all or part of the current object in format suitable
        for the screen.
    */
    virtual int comm_preview(std::vector<std::string> &sv, bool itive_com);

    /** \brief Sample a distribution

        For objects of type prob_dens_mdim_gaussian:

        Arguments: <number of samples>

        Sample the distribution to create a <tt>table</tt> object.

        For objects of type prob_dens_mdim_gmm:

        Arguments: <number of samples>

        Sample the Gaussian mixture to create a <tt>table</tt> object.

        For objects of type prob_dens_mdim_amr:
        
        Arguments: <number of samples>

        Sample the distribution to create a <tt>table</tt> object.
     */
    virtual int comm_sample(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select part of an object

        For objects of type table:

        Select columns for a new table.

        Arguments: <tt><column, pattern, or function 1> 
        [column, pattern, or function 2] ...</tt>

        The \c select command creates a new table from the present
        table, including some or all of the columns based on the
        arguments. Each argument can be a column, a pattern, or a
        assignment and a new function, If the argument is a column,
        the specified column is added to the new table. If the
        argument begins with a colon ':', then the remainder of the
        string is interpreted as a pattern. All columns which match
        the pattern are added to the new table. If a column in the old
        table is specified multiple times in the arguments, then it is
        only included in the new table once. Finally, if the argument
        contains an equals sign, '=', then the string to the left of
        '=' is interpreted as a new column name and the string to the
        right of the column is interpreted as a function (built from
        the columns in the old table).

        Depending on the value of \c use_regex, patterns are built
        upon the rules of fnmatch() or regex. For fnmatch(), '*'
        represents multiple and '?' represents a single character.

        For example, given a table with columns 'c1dab c2dxy c13d
        c13def' the command <tt>-select c1dab :c?d*
        d2=c13d+c13def</tt> creates a new table with columns 'c1dab
        c2dxy d2' where the third column in the new table is given by
        the sum of the third and fourth columns in the old table.

        For objects of type table3d:

        Select columns for a new table3d object.

        Arguments: <tt><slice, pattern, or function 1> 
        [slice, pattern, or function 2] ...</tt>

        The \c select command creates a new \c table3d object from the
        present \c table3d object, including some or all of the slices
        based on the arguments. Each argument can be a slice, a
        pattern, or a assignment and a new function, If the argument
        is a slice, the specified slice is added to the new \c table3d
        object. If the argument begins with a colon ':', then the
        remainder of the string is interpreted as a pattern. All
        slices which match the pattern are added to the new \c table3d
        object. If a slice in the old \c table3d object is specified
        multiple times in the arguments, then it is only included in
        the new \c table3d object once. Finally, if the argument
        contains an equals sign, '=', then the string to the left of
        '=' is interpreted as a new slice name and the string to the
        right of the slice is interpreted as a function (built from
        the slices in the old \c table3d object).

        Depending on the value of \c use_regex, patterns are built
        upon the rules of fnmatch() or regex. For fnmatch(), '*'
        represents multiple and '?' represents a single character.

        For example, given a \c table3d object with slices 'c1dab
        c2dxy c13d c13def' the command <tt>-select c1dab :c?d*
        d2=c13d+c13def</tt> creates a new \c table3d object with
        slices 'c1dab c2dxy d2' where the third slice in the new \c
        table3d object is given by the sum of the third and fourth
        slices in the old \c table3d object.
    */
    virtual int comm_select(std::vector<std::string> &sv, bool itive_com);

    /** \brief Select rows from an object

        For objects of type table:

        Select rows from the table
        
        Arguments: <tt><row specification></tt>

        Select the rows from a table for which the function (based on
        columns in the current \c table) in <row_spec> evaluates to a
        number greater than 0.5.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

    */
    virtual int comm_select_rows(std::vector<std::string> &sv,
                                 bool itive_com);

    /** \brief Convert a series of histograms to a table3d object

        For objects of type table:

        Store a histogram series in a table3d object.

        Arguments: <tt><grid vector spec.> <direction: "x" or
        "y"> <grid name> <bin spec.> <bin name> <pattern> <new slice></tt>

        Construct a series of histograms from a series of columns in a
        <tt>table</tt> object and then store them in a
        <tt>table3d</tt> object. The <tt>ser-hist-t3d</tt> command
        begins by creating a set of histograms, one for each entry of
        the grid specified in <grid vector spec.>. This grid is used
        as either the x or the y-grid in the new table3d object,
        depending on which is given for the <direction> argument. 

        The bin edges for the histogram are specified in a "bin
        specification", <bin spec.>, which consists of two arguments.
        The first argument is the bin edges for the histograms, the
        second is the grid values of the table3d object. The grid
        values in the table3d object can be computed automatically as
        the average of the bin edges if the word "auto" is given for
        the second argument of the bin specification. The 
        linear_or_log() function is used to attempt to automatically
        determine logarithmic bin edges. 

        The data for the histograms is pulled from the columns
        specified in the search pattern <pattern>. There must be one
        column for each histogram, and thus each element in <grid
        vector spec.>. The <bin name> argument is used as the name of
        the bin grid in the table3d object. Finally, the <new slice>
        argument is the name of the new slice in the table3d object
        which stores the histogram data. 

        Finally, if the bin specification is the string "auto" and
        then a number of bins, then the bin edges and the associated
        grid are automatically computed by the minimum and maximum
        value over all the columns specified in <pattern>.
    */
    virtual int comm_ser_hist_t3d(std::vector<std::string> &sv,
                                  bool itive_com);

    /// Post-processing for setting a value
    virtual int comm_set(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set an individual data point at a specified row and
        column

        For objects of type table:

        Set the entries of a column.

        Arguments: <tt><row_spec> <col> <val_spec></tt>

        Set the value of rows specifed by the 'row_spec' function in
        column 'col' to the value given by the 'val_spec' function.
        Rows are chosen if row_spec evaluates to a number greater than
        0.5.

        For objects of type table3d:

        Set the entries of a slice.

        Arguments: <tt><x value> <y value> <z name> <val></tt>

        Set the value of the slice named 'z name' at the grid point
        closest to (<x value>,<y value>) to the value <val>.
    */
    virtual int comm_set_data(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Set the grid for a \ref o2scl::tensor_grid object

        For objects of type tensor_grid:

        Set the tensor grid.

        Arguments: <tt><index> <func. or vector spec></tt>

        The first argument for the \c set-grid command specifies the
        index for which grid to set. The second argument specifies the
        grid. If it contains a ':', it is assumed to be a vector
        specification. Otherwise, the
        argument is assumed to be a function which specifies the grid
        value as a function of the variables 'i', 'm', and 'x'. The value of
        'i' ranges from 0 to m-1, where 'm' is the tensor size for
        each rank and the value of 'x' is equal to the previous grid
        value.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::vector_spec()` for help on vector
        specifications and :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim

        For objects of type table3d:

        Set either the x or y grid of the table3d object.

        Arguments: <tt><index or 'x' or 'y'> <func. or vector spec></tt>

        The first argument for the \c set-grid command specifies the
        index or label for which grid to set. The second argument
        specifies the grid. If it contains a ':', it is assumed to be
        a vector specification. Otherwise, the argument is assumed to
        be a function which specifies the grid value as a function of
        the variables 'i', 'j', 'm', 'n', 'x', and 'y'. The value of
        'i' ranges from 0 to 'm'-1, where 'm' is the size of the x
        index. Similarly, the value lf 'j' ranges from 0 to 'n'-1,
        where 'n' is the size of the y index. The variables 'x' and
        'y' hold the x and y grid values.

        \verbatim embed:rst
        See :cpp:func:`o2scl_hdf::vector_spec()` for help on vector
        specifications and :cpp:func:`o2scl_hdf::functions()` for help
        on function specifications.
        \endverbatim
    */
    virtual int comm_set_grid(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set units of a column

        For objects of type table:

        Set the units for a specified column.

        Arguments: <tt><column> <unit></tt>

        Set the unit string of <column>. Any string is allowed, but
        only those based on the output of <tt>acol -convert list</tt>
        can be automatically converted.
    */
    virtual int comm_set_unit(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Send a slack message

        Arguments: <tt>["#channel"] <strings-spec></tt>

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

    /** \brief Extract a slice from a table3d object to generate a 
        \ref o2scl::table object

        For objects of type table3d:

        Convert a slice to a table object.

        Arguments: <tt><"x" or "y"> <val></tt>

        Extract a slice of a table3d object at fixed x or fixed y to
        create a new table object. This function uses interpolation
        with the current interpolation type to interpolate all of the
        slices in the table3d object to create a table with a column
        for each slice.

        TODO: explain how this is different from the \c to-table
        command.

        For objects of type tensor_grid:

        Slice to a smaller rank tensor_grid object.

        Arguments: <tt><index 1> <value 1> <index 2> <value 2> ...</tt>

        Detailed desc.

        TODO: explain how this is different from to-table
    */
    virtual int comm_slice(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert a slice to a histogram

        For objects of type table3d:

        Construct a histogram from a slice.

        Arguments: <tt><slice></tt>

        Detailed desc.
    */
    virtual int comm_slice_hist(std::vector<std::string> &sv, bool itive_com);

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

        For objects of type string[]:

        Sort the vector.

        (No arguments.)

        Sort the vector (in-place).

        For objects of type table:

        Sort the entire table by one column.

        Arguments: <tt><column> [unique]</tt>

        Sorts the entire table by the column specified in <column>.
        If the word "unique" is specified as the second argument, then
        delete duplicate rows after sorting.
    */
    virtual int comm_sort(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get object statistics

        For objects of type table:

        Show column statistics.

        Arguments: <tt><column></tt>

        Output the average, standard deviation, max and min of
        <column>.

        For objects of type double[]:

        Show vector statistics.

        Arguments: (None.)

        Output the average, standard deviation, max and min.

        For objects of type table3d:

        Show slice statistics.

        Arguments: <tt><slice></tt>

        Output the average, standard deviation, max and min of
        <slice>.

        For objects of type tensor:

        Show tensor statistics.

        Arguments: (None.)

        The <tt>stats</tt> command outputs the number of entries,
        their mean, standard deviation, minimum and maximum. It also
        counts the number of infinite or NaN values.

        For objects of type tensor_grid:

        Show tensor statistics.

        Arguments: (None.)

        The <tt>stats</tt> command outputs the number of entries,
        their mean, standard deviation, minimum and maximum. It also
        counts the number of infinite or NaN values.
    */
    virtual int comm_stats(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute the sum of two objects

        For objects of type table:

        Add data from a second table object to current table.

        Arguments: <tt><file> [name]</tt>

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
        
        Arguments: <tt><file> [name]</tt>

        Add all slides from the second table3d to their corresponding
        slices in the current table3d, creating new slices if
        necessary.

        For objects of type tensor:

        Output the sum of all the tensor entries.

        (No arguments.)

        The \c sum command outputs the total tensor size
        and the sum over all entries. Note, to perform a partial
        sum over sum of the tensor indices, use the 
        <tt>rearrange</tt> command.

        For objects of type tensor_grid:

        Output the sum of all the tensor entries.

        (No arguments.)

        The \c sum command outputs the total tensor size and the sum
        over all entries. Note, to perform a partial sum over sum of
        the tensor indices, use the <tt>rearrange</tt> command.
    */
    virtual int comm_sum(std::vector<std::string> &sv, bool itive_com);

    /** \brief Construct a multivariate Gaussian distribution

        For objects of type table:
        
        Construct a multivariate Gaussian distribution

        Arguments: <column 1> [column 2] ...

        This creates an object of type <tt>prob_dens_mdim_gaussian</tt>
        based on the given columns of data in the table.
     */
    virtual int comm_to_gaussian(std::vector<std::string> &sv,
                                 bool itive_com);

    /** \brief Construct a Gaussian mixture model

        For objects of type table:
        
        Construct a multivariate Gaussian distribution

        Arguments: <number of Gaussians> <column 1> [column 2] ...
        
        This creates an object of type <tt>prob_dens_mdim_gmm</tt>
        based on the given columns of data in the table.
     */
    virtual int comm_to_gmm(std::vector<std::string> &sv,
                                 bool itive_com);
    
    /** \brief Construct a KDE

        For objects of type table:
        
        Construct a KDE

        Arguments: <options or 'none'> <column 1> [column 2] ...
        
        This creates an object of type <tt>prob_dens_mdim_kde</tt>
        based on the given columns of data in the table.
    */
    virtual int comm_to_kde(std::vector<std::string> &sv,
                                 bool itive_com);
    
    /** \brief Construct a AMR-based probability distribution

        For objects of type table:

        Construct a AMR-based probability distribution

        Arguments: <column 1> [column 2] ...

        This creates an object of type <tt>prob_dens_mdim_amr</tt>
        based on the given columns of data in the table.
     */
    virtual int comm_to_pdma(std::vector<std::string> &sv,
                             bool itive_com);

    /** \brief Convert to a \ref o2scl::hist object

        For objects of type table:

        Convert a table column to a histogram.

        Arguments: <tt><col> <n_bins> [wgts]</tt>

        The <tt>to-hist</tt> command creates a 1D histogram from
        column <col> using exactly <n_bins> bins and (optionally)
        weighting the entries by the values in column [wgts]. 

        For objects of type table3d:

        Convert a table3d slice to a histogram.

        Arguments: <tt><slice> <n_bins></tt>

        The <tt>to-hist</tt> command creates a 1D histogram from
        slice <slice> using exactly <n_bins> bins.
    */
    virtual int comm_to_hist(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert to a \ref o2scl::hist_2d object

        For objects of type table:

        Convert a table to a 2d histogram.

        Arguments: <tt><col x> <col y> <n_x_bins> <n_y_bins>
        [wgts]</tt>

        The <tt>to-hist-2d</tt> command creates a 2D histogram from
        <col x> and <col y> using <n_x_bins> bins in the x direction
        and <n_y_bins> bins in the y direction, optionally weighting
        the entries by the column [wgts].

        For objects of type table3d:

        Convert a table3d slice to a 2d histogram.

        <slice>

        The <tt>to-hist-2d</tt> command creates a 2D histogram from
        slice <slice>.
    */
    virtual int comm_to_hist_2d(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table object

        For objects of type double[]:

        Convert to a table object.

        Arguments: <tt><column name></tt>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type int[]:

        Convert to a table object.

        Arguments: <tt><column name></tt>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type size_t[]:

        Convert to a table object.

        Arguments: <tt><column name></tt>
        
        Convert the vector to a table with a single column named
        <column name>.

        For objects of type tensor_grid:

        Convert to a table object.

        Arguments: <tt><index> <grid name> <data name> [values of
        fixed indices]</tt>
        
        Convert the \c tensor_grid object to a \c table object by
        choosing an index to vary and fixing the remaining indices.
        The resulting table has two columns. The grid associated with
        the specified index is stored in a new column named <grid
        name> and the values of the tensor are stored in a column
        named <data name>. Linear interpolation is used, so the values
        for the fixed indices need not lie on grid points.

        For objects of type table3d:

        Convert to a table object.

        Arguments: (No arguments.)
        
        Given a \c table3d object with C slices which each have N
        values in the x direction and M values in the y direction,
        this command creates a table with C+2 columns and N times M
        rows. The x name from the \c table3d object is the first
        column, the y name from the \c table3d object is the 
        second column, and the remaining columns correspond to the
        slices. Then the \c to-table command loops through all
        of the x and y values, creating a row from the values in
        each slice. 

        For objects of type hist:

        Convert to a \c hist object to a \c table object.

        Arguments: (No arguments.)
        
        This function creates a new table with four columns from
        the current \c hist object. The columns are named
        "rep", "low", "high", "wgt". In order, for each bin,
        each row stores the bin representative, the lower edge,
        the upper edge, and the bin weight. 

        For objects of type prob_dens_mdim_kde:

        Convert to a \c prob_dens_mdim_kde object to a \c table object.

        Arguments: (No arguments.)
        
        This function creates a new table...
    */
    virtual int comm_to_table(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object

        For objects of type table:

        Convert a <tt>table</tt> to a <tt>table3d</tt> object.

        Arguments: <tt><x column> <y column> [empty value] [eps]</tt>

        The \c to-table3d command converts the entire \c table to a a
        \c table3d object assuming that the \c table is consists of a
        list of entries in the new \c table3d object. The <x column>
        and <y column> arguments are used as the data for the x and y
        grids in the new \c table3d object. If [empty value] is given,
        then this value is used for points not given by the table. If
        [eps] is specified, then that value (instead of the default of
        10^{-12}) is used as the minimum value between grid points.

        For objects of type tensor:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x index> <y index> <slice name> [fixed 1]
        [fixed 2] ...</tt>

        The \c to-table3d command uses two indices in the current
        tensor object to create a \c table3d object. The values for
        the remaining indices fixed to [fixed 1], [fixed 2], etc. in
        that order. For example, <tt>to-table3d 3 1 z 5 3</tt> uses
        index 3 for the x coordinate of the new table3d object, uses
        index 1 for the y coordinate of the new table3d object, uses 5
        for index 0, and uses 3 for index 2. The x- and y-grids in
        he table3d object are named "x" and "y" and filled with the
        grid index by default. To set the x- or y-grid names
        afterwards, use commands 'x-name' and 'y-name'.

        For objects of type tensor<int>:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x index> <y index> <slice name> [fixed 1]
        [fixed 2] ...</tt>

        This command uses two indices in the current tensor object to
        create a table3d object. The values for the remaining indices
        fixed to [fixed 1], [fixed 2], etc. in that order. For
        example, <tt>to-table3d 3 1 z 5 3</tt> uses index 3 for the x
        coordinate of the new table3d object, uses index 1 for the y
        coordinate of the new table3d object, uses 5 for index 0, and
        uses 3 for index 2. The x- and y-grids in he table3d object
        are named "x" and "y" and filled with the grid index by
        default. To set the x- or y-grid names afterwards, use
        commands 'x-name' and 'y-name'.

        For objects of type tensor<size_t>:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x index> <y index> <slice name> [fixed 1]
        [fixed 2] ...</tt>

        This command uses two indices in the current tensor object to
        create a table3d object. The values for the remaining indices
        fixed to [fixed 1], [fixed 2], etc. in that order. For
        example, <tt>to-table3d 3 1 z 5 3</tt> uses index 3 for the x
        coordinate of the new table3d object, uses index 1 for the y
        coordinate of the new table3d object, uses 5 for index 0, and
        uses 3 for index 2. The x- and y-grids in he table3d object
        are named "x" and "y" and filled with the grid index by
        default. To set the x- or y-grid names afterwards, use
        commands 'x-name' and 'y-name'.

        For objects of type tensor_grid:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x index> <y index> <slice name> [value 1]
        [value 2] ...</tt>

        This command uses two indices in the current tensor_grid
        object to create a table3d object. The values for the
        remaining indices are by interpolation to [value 1], [value
        2], etc. in that order. For example, <tt>to-table3d 3 1 z 0.5
        2.0</tt> uses index 3 for the x coordinate of the new table3d
        object, uses index 1 for the y coordinate of the new table3d
        object, uses interpolation to set the value of the index 0 to
        0.5, and uses interpolation to set the value of index 2 to to
        2.0. The x- and y-grids in the table3d object are named "x"
        and "y" by default. To set the x- or y-grid names
        afterwards, use commands 'x-name' and 'y-name'.

        For objects of type hist_2d:

        Convert to a hist_2d object.

        Arguments: <tt><x name> <y name> <weight name></tt>

        Convert to a hist_2d object using the specified names.
        
        For objects of type prob_dens_mdim_amr:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x index> <y index> <x name> <x points> <y
        name> <y points> <slice name></tt>

        Select two indices and convert to a table3d object.

        For objects of type prob_dens_mdim_kde:

        Convert to a \c prob_dens_mdim_kde object to a \c table3d object.

        Arguments: (No arguments.)
        
        This function creates a new table3d...
    */
    virtual int comm_to_table3d(std::vector<std::string> &sv, bool itive_com);

    /** \brief Convert object to a \ref o2scl::table3d object
        by summing over tensor indices

        For objects of type tensor:

        Select two indices and convert to a table3d object.

        Arguments: <tt><x name> <x index> <y name> <y index> 
        <slice name></tt>

        Convert the \c tensor object to a \c table3d object by
        assigning index <x index> to the x value, index <y index> to
        the <y value>, and summing over all remaining indices to
        create a slice of named <slice name>.
    */
    virtual int comm_to_table3d_sum(std::vector<std::string> &sv,
                                    bool itive_com);

    /** \brief Convert object to a \ref o2scl::tensor object

        For objects of type tensor_grid:

        Convert to a tensor object.

        Arguments: (No arguments.)
        
        Convert to a tensor object, removing the grid information.
    */
    virtual int comm_to_tensor(std::vector<std::string> &sv,
                               bool itive_com);

    /** \brief Convert object to a \ref o2scl::tensor_grid object

        For objects of type table3d:

        Convert a slice of the table3d object to a tensor_grid object.

        Arguments: <tt><slice></tt>

        Convert a \c table3d object to a rank 2 \c tensor_grid object,
        directly copying the grid over.

        For objects of type tensor:

        Convert the tensor to a tensor_grid object.

        Arguments: <tt>[function 1] [function 2] ...</tt>

        Convert a \c tensor to a \c tensor_grid object, using
        functions to specify the grid for each index. The functions
        should be specified as functions of the variable 'i', which
        runs from 0 to size-1 for each index. Any user-specified
        functions are used up to the rank of the tensor, and if not
        enough functions are specified, then the function 'i' is used.
    */
    virtual int comm_to_tensor_grid(std::vector<std::string> &sv,
                                    bool itive_com);

    /** \brief Convert object to a \ref o2scl::tensor_grid object

        For objects of type table3d:

        Convert a slice of the table3d object to a tensor_grid object.

        Arguments: <tt><slice> <n points></tt>

        Detailed desc.

        For objects of type tensor:

        Convert the tensor to a tensor_grid object.

        Arguments: <tt>[function 1] [function 2] ...</tt>

        Convert a tensor to a tensor_grid object, using functions to
        specify the grid for each index. The functions should be
        specified as functions of the variable 'i', which runs from 0
        to size-1 for each index. Any user-specified functions are
        used up to the rank of the tensor, and if not enough functions
        are specified, then the function 'i' is used.
    */
    virtual int comm_to_tg_fermi(std::vector<std::string> &sv,
                                 bool itive_com);
    
    /** \brief Output the type of the current object

        Arguments: (No arguments.)
        
        Show the current object type. Type can be <tt>char</tt>,
        <tt>double</tt>, <tt>double[]</tt>, <tt>hist</tt>,
        <tt>hist_2d</tt>, <tt>int</tt>, <tt>int[]</tt>,
        <tt>prob_dens_mdim_amr</tt>, <tt>prob_dens_mdim_gaussian</tt>,
        <tt>size_t</tt>, <tt>size_t[]</tt>, <tt>string</tt>,
        <tt>string[]</tt>, <tt>table</tt>, <tt>table3d</tt>,
        <tt>tensor</tt>, <tt>tensor<int></tt>,
        <tt>tensor<size_t></tt>, <tt>tensor_grid</tt>,
        <tt>uniform_grid<double></tt>, <tt>vec_vec_double</tt>,
        <tt>vec_vec_string</tt>, or <tt>vector<contour_line></tt>.
    */
    virtual int comm_type(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Get or set the value of an object

        For objects of type int:

        Get or set the integer.

        Arguments: <tt>[value]</tt>

        Get or set the integer.

        For objects of type size_t:

        Get or set the size_t object.

        Arguments: <tt>[value]</tt>

        Get or set the size_t object.

        For objects of type string:

        Get or set the string.

        Arguments: <tt>[value]</tt>

        Get or set the string.

        For objects of type double:

        Get or set the value of the double object. 

        Arguments: <tt>[value spec.]</tt>

        Get or set the value of the double object. See ``Value
        specifications`` for more information.

        For objects of type char:

        Get or set the character.

        Arguments: <tt>[value]</tt>

        Get or set the character.

        For objects of type double[]:

        Get or set an entry in the array

        Arguments: <tt><index> [value]</tt>

        Get or set entry at index <index>.

        For objects of type int[]:

        Get or set an entry in the array

        Arguments: <tt><index> [value]</tt>

        Get or set entry at index <index>.

        For objects of type size_t[]:

        Get or set an entry in the array

        Arguments: <tt><index> [value]</tt>

        Get or set entry at index <index>.

        For objects of type string[]:

        Get or set an entry in the array

        Arguments: <tt><index> [value]</tt>

        Get or set entry at index <index>.

        For objects of type table:

        Get or set a single entry in a table.

        Arguments: <tt><column> <row> [value or "none"]</tt>

        This command gets or sets the value in the specified column
        and row. If "none" is specified as the third argument, then
        the \c value command just prints out the specified entry as if
        the third argument was not specified.

        To refer to a location by the x and y values instead of 
        indices, use the \c value-grid command.

        For objects of type table3d:

        Get or set a single entry in a table3d object.

        Arguments: <tt><slice> <x index> <y index> [value or "none"]</tt>

        This command gets or sets the value in the specified slice at
        the location specified by <x index> and <y index>. If "none"
        is specified as the fourth argument, or if only three
        arguments are given, then the \c value command just prints out
        the specified value.

        To refer to a location by the x and y values instead of 
        indices, use the \c value-grid command.

        For objects of type tensor:

        Get or set a single entry in a tensor object.

        Arguments: <tt><index 1> <index 2> <index 3> ... [value or
        "none"]</tt>

        This command gets or sets the value in the tensor at the
        location given by the specified indices. If an an argument is
        given at the end and that arguemnt is not "none", then it is
        used to set the new value. If only the indices are given, then
        the \c value command just prints out the specified value.

        For objects of type tensor_grid:

        Get or set a single entry in a tensor_grid object.

        Arguments: <tt><index 1> <index 2> <index 3> ... [value or "none"]</tt>

        The \c value command gets or sets a value in the \c
        tensor_grid object. The arguments are a list of indices and
        (optionally) a new value to store in that location.

        See the \c value-grid command to specify a grid location
        rather than specifying indices.
    */
    virtual int comm_value(std::vector<std::string> &sv, bool itive_com);
    
    /** \brief Get an entry by grid point

        For objects of type table:

        Get or set a single entry in a table.
        
        Arguments: <tt><index column> <index value> <target column>
        [value or "none"]</tt>
        
        The \c value-grid command first looks for the value closest to
        <index value> in the column <index column> to determine a row
        in the table. Next \c value-grid gets or sets the value of the
        target column in that row. If "none" is specified as the
        fourth argument, then \c value-grid just prints out the
        specified entry as if the third argument was not specified.

        For objects of type table3d:

        Get a single entry in a table3d object.

        Arguments: <tt><slice> <x value> <y value> [value or "none"]</tt>

        The \c value-grid command first looks for the value closest to
        <x value> and <y value> in the slice <index slice> to
        determine a row in the table. Next, \c value-grid gets or sets
        the value of of the specified slice in that location. If
        "none" is specified as the fourth argument, then \c value-grid
        just prints out the specified entry as if the third argument
        was not specified.

        For objects of type tensor_grid:

        Get a single entry in a \c tensor_grid object.

        Arguments: <tt><value 1> <value 2> <value 3> ... [value or
        "none"]</tt>

        The \c value-grid command gets or sets a value in the
        \c tensor_grid object. The arguments are a list of grid values
        and (optionally) a new value to store in the location closest
        to the specified grid values.
    */
    virtual int comm_value_grid(std::vector<std::string> &sv, bool itive_com);

    /** \brief Print version information and O₂scl settings.

        Arguments: (No arguments.)

        This prints the O₂scl version, when it was compiled, HDF5
        version information, what packages were enabled during
        installation, the data and documentation directories, and
        other miscellaneous information about O₂scl.
     */
    virtual int comm_version(std::vector<std::string> &sv, bool itive_com);

    /** \brief Open remote HTML docs for acol or an O₂scl topic.

        Arguments: <tt>["dev"] [search_term], [topic] or [section
        search_term]</tt>
        
        If no arguments are given, this command opens up the remote
        HTML documentation for acol in the default web browser using
        \c open on OSX and \c xdg-open on other systems. If a help
        topic, [topic] is specified, then the associated O₂scl web
        page is opened. If the argument does not match an already
        known topic, then the search feature on the O₂scl web page is
        opened using the specified search term. Note that, for search
        terms, spaces can be included using e.g.

        <tt>-wdocs \"Simulated annealing\"</tt>

        If the optional argument "dev" is given, then the development
        rather than release documentation is used. In order to open
        the local version of the documentation instead of the remote
        copy, use <tt>docs</tt> instead of <tt>wdocs</tt>.
    */
    virtual int comm_wdocs(std::vector<std::string> &sv, bool itive_com);

    /** \brief Get weighted statistics

        For objects of type table:

        Show weighted column statistics.

        Arguments: <tt><column> <weights></tt>

        Output the average, standard deviation, max and min of
        <column>, using the weights specified in <weights>.
    */
    virtual int comm_wstats(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the name of the x grid

        For objects of type table3d:

        Get or set the name of the x grid.

        Arguments: <tt>[name]</tt>

        Get the name of the x grid, or, if [name] is specified, change
        the x grid name to [name].
    */
    virtual int comm_x_name(std::vector<std::string> &sv, bool itive_com);

    /** \brief Parse doxygen XML to generate runtime docs.
        
        Arguments: (No arguments.)

        When pugixml is enabled, this function reads the doxygen XML
        output and generates an HDF5 file which acol reads to generate
        the runtime documentation. This command is principally
        designed for developers and requires several additional tools
        not otherwise required during installation.
    */
    virtual int comm_xml_to_o2(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set the name of the y grid

        For objects of type table3d:

        Get or set the name of the y grid.

        Arguments: <tt>[name]</tt>

        Get the name of the y grid, or, if [name] is specified, change
        the y grid name to [name].
    */
    virtual int comm_y_name(std::vector<std::string> &sv, bool itive_com);
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
    
    /** \brief Validate the setting of \ref interp_type
        
        Used in \ref comm_set().
     */
    int validate_interp_type();
    
    // End of class acol_manager
  };

}

#endif
