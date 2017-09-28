/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
/** \file cli.h
    \brief File defining command-line interface in \ref o2scl::cli
*/
#ifndef O2SCL_CLI_H
#define O2SCL_CLI_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>

#include <o2scl/columnify.h>
#include <o2scl/vector.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Base for \ref o2scl::cli command function
      
      See the \ref o2scl::cli class for more details.
  */
  class comm_option_funct {

  public:

    comm_option_funct() {}

    virtual ~comm_option_funct() {}

    /// The basic function called by \ref o2scl::cli
    virtual int operator()(std::vector<std::string> &cstr, bool itive_com)=0;

  };

  /// Function pointer for \ref o2scl::cli command function
  class comm_option_fptr : public comm_option_funct {

  public:

    /// Create from a member function pointer from the specified class
    comm_option_fptr(int (*fp)(std::vector<std::string> &, bool)) {
      fptr=fp;
    }

    virtual ~comm_option_fptr() {}

    /// The basic function called by \ref o2scl::cli
    virtual int operator()(std::vector<std::string> &cstr, bool itive_com) {
      return (*fptr)(cstr,itive_com);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The pointer to the member function
    int (*fptr)(std::vector<std::string> &cstr, bool itive_com);

    // Copy constructor
    //comm_option_fptr(const comm_option_fptr &f) {
    //fptr=f.fptr;
    //}
    
    /// Copy constructor
    comm_option_fptr& operator=(const comm_option_fptr &f) {
      fptr=f.fptr;
      return *this;
    }
    
#endif
    
  };

  /// Member function pointer for \ref o2scl::cli command function
  template<class tclass> class comm_option_mfptr : public comm_option_funct {

  public:

    /// Create from a member function pointer from the specified class
    comm_option_mfptr(tclass *tp, int (tclass::*fp)(std::vector<std::string> &,
						bool)) {
      tptr=tp;
      fptr=fp;
    }

    virtual ~comm_option_mfptr() {}

    /// The basic function called by \ref o2scl::cli
    virtual int operator()(std::vector<std::string> &cstr, bool itive_com) {
      return (*tptr.*fptr)(cstr,itive_com);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The pointer to the member function
    int (tclass::*fptr)(std::vector<std::string> &cstr, bool itive_com);

    /// The pointer to the class 
    tclass *tptr;

    /// Copy constructor
    comm_option_mfptr(const comm_option_mfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
    }
    
    /// Copy constructor
    comm_option_mfptr& operator=(const comm_option_mfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      return *this;
    }

#endif
    
  };

  /** \brief Command for interactive mode in \ref o2scl::cli

      See the \ref o2scl::cli class for more details.

      \comment
      This was at one point converted into a class, but it wasn't that
      easy to use, in comparison to structs which are easy to
      initialize in aggregate.
      \endcomment
  */
  typedef struct {

    /// Short option (\c '\\0' for none, must be unique if present)
    char shrt;
    /// Long option (must be specified and must be unique)
    std::string lng;
    /// Description for help
    std::string desc;
    /// Minimum number of parameters (0 for none, -1 for variable)
    int min_parms;
    /// Maximum number of parameters (0 for none, -1 for variable)
    int max_parms;
    /// Description of parameters
    std::string parm_desc;
    /// The help description
    std::string help;
    /// The pointer to the function to be called (or 0 for no function)
    comm_option_funct *func;
    /// Type: command-line parameter, command, or both
    int type;

  } comm_option_s;
  
  /** \brief A command-line argument for \ref o2scl::cli

      This is the internal structure that \ref o2scl::cli uses to package
      command-line arguments.
   */
  typedef struct {
    /// The argument
    std::string arg;
    /// Is an option?
    bool is_option;
    /// Is a properly formatted option
    bool is_valid;
    /// List of parameters (empty, unless it's an option)
    std::vector<std::string> parms;
    /// A pointer to the appropriate option (0, unless it's an option)
    comm_option_s *cop;
  } cmd_line_arg;

  /** \brief Configurable command-line interface

      This class is experimental.

      Default commands: help, get/set, quit, exit, '!', verbose, license,
      warranty, alias, run.

      Note that if the shell command is allowed (as it is by default)
      there are some potential security issues which are not solved
      here.
      
      Commands which begin with a '#' character are ignored.

      \note In interactive mode, commands are limited to 300 characters.

      \future Warn in run_interactive() when extra parameters are given
      \future Include a "remove command" function
      \future A replace command function, there's already some code
      in cli.cpp for this.
      \future There's some code duplication between comm_option_run()
      and run_interactive()
      \future Allow the user to set the tilde string
      \future Disallow direct access to \ref o2scl::cli::par_list in order to
      ensure parameter names do not contain whitespace

      <b>Concepts</b>

      As a matter of definition, the command-line arguments are simply
      called arguments. These arguments may be options (in which case
      they begin with either one dash or two) or parameters to these
      options.  When run in interactive mode, these options are also
      commands.
      
  */
  class cli {

  public:
    
    /// Parameter for \ref o2scl::cli
    class parameter {
      
    public:

      virtual ~parameter() {}
      
      /// Help description
      std::string help;

      /// Set from string
      virtual int set(std::string s)=0;

      /// Convert to string
      virtual std::string get()=0;
      
    };

    /// String parameter for \ref o2scl::cli
    class parameter_string : public parameter {

    public:

      virtual ~parameter_string() {}

      /// Parameter
      std::string *str;

      /// Set from string
      virtual int set(std::string s) {
	*str=s;
	return 0;
      }

      /// Convert to string
      virtual std::string get() {
	return *str;
      }
      
    };

    /// String parameter for \ref o2scl::cli
    class parameter_bool : public parameter {

    public:

      virtual ~parameter_bool() {}

      /// Parameter
      bool *b;

      /// Set from string
      virtual int set(std::string s) {
	*b=o2scl::stob(s);
	return 0;
      }

      /// Convert to string
      virtual std::string get() {
	return btos(*b);
      }
      
    };

    /// Double parameter for \ref o2scl::cli
    class parameter_double : public parameter {

    public:
      
      parameter_double() {
	parse_strings=true;
      }

      virtual ~parameter_double() {}

      /// Parameter
      double *d;

      /// If true, 
      bool parse_strings;

      /// Set from string
      virtual int set(std::string s) {
	if (parse_strings) {
	  *d=function_to_double(s);
	} else {
	  *d=o2scl::stod(s);
	}
	return 0;
      }

      /// Convert to string
      virtual std::string get() {
	return dtos(*d);
      }
      
    };

    /// Integer parameter for \ref o2scl::cli
    class parameter_int : public parameter {

    public:
      
      virtual ~parameter_int() {}

      /// Parameter
      int *i;

      /// Set from string
      virtual int set(std::string s) {
	*i=o2scl::stoi(s);
	return 0;
      }

      /// Convert to string
      virtual std::string get() {
	return itos(*i);
      }
      
    };

    /// Integer parameter for \ref o2scl::cli
    class parameter_size_t : public parameter {

    public:
      
      virtual ~parameter_size_t() {}

      /// Parameter
      size_t *s;

      /// Set from string
      virtual int set(std::string st) {
	return o2scl::stoszt_nothrow(st,*s);
      }

      /// Convert to string
      virtual std::string get() {
	return szttos(*s);
      }
      
    };

    /// \name Parameter storage and associated iterator type
    //@{
    /// Parameter list
    std::map<std::string,parameter *,std::less<std::string> > par_list;
    /// List iterator
    typedef std::map<std::string,parameter *,
      std::greater<std::string> >::iterator par_t;
    //@}
    
#ifndef DOXYGEN_NO_O2NS_INTERNAL
    
  protected:
    
    /// Output the parameter list
    int output_param_list();

    /** \brief Attempt to expand a tilde to a user's home directory

        Experimental and currently unused.
     */
    int expand_tilde(std::vector<std::string> &sv);

    /// Replace all occurences of \c sold with \c snew in \c sv
    int apply_alias(std::vector<std::string> &sv, 
		    std::string sold, std::string snew);

    /// Control screen output
    int verbose;

    /// \name The hard-coded command functions
    //@{
    int comm_option_alias(std::vector<std::string> &sv, bool itive_com);
    int comm_option_commands(std::vector<std::string> &sv, bool itive_com);
    int comm_option_get(std::vector<std::string> &sv, bool itive_com);
    int comm_option_help(std::vector<std::string> &sv, bool itive_com);
    int comm_option_license(std::vector<std::string> &sv, bool itive_com);
    int comm_option_no_intro(std::vector<std::string> &sv, bool itive_com);
    int comm_option_run(std::vector<std::string> &sv, bool itive_com);
    int comm_option_set(std::vector<std::string> &sv, bool itive_com);
    int comm_option_warranty(std::vector<std::string> &sv, bool itive_com);
    //@}

    /// Storage for getline
    char buf[300];

    /// Storage for the function to call after setting a parameter
    comm_option_funct *user_set_func;
    
    /// List of commands
    std::vector<comm_option_s> clist;
  
    /// \name Help for parameters
    //@{
    std::vector<std::string> ph_name, ph_desc;
    //@}

    /// \name Aliases
    //@{
    std::map<std::string,std::string,std::greater<std::string> > als;
    typedef std::map<std::string,std::string,
      std::greater<std::string> >::iterator al_it;
    //@}
    
    /// Compare two strings, treating dashes and underscores as equivalent
    bool string_equal_dash(std::string s1, std::string s2);

#endif
  
  public:

    cli();

    virtual ~cli();

    /// String to replace tildes with
    std::string tilde_string;

    /** \brief If true, output the usual GNU intro when run_interactive() 
	is called (default true).

	In order to conform to GNU standards, this ought not be set to
	false by default.
    */
    bool gnu_intro;

    /// Function to call when a \c set command is issued
    int set_function(comm_option_funct &usf) {
      user_set_func=&usf;
      return 0;
    }

    /// \name Value to indicate whether commands are also command-line options
    //@{
    static const int comm_option_command=0;
    static const int comm_option_cl_param=1;
    static const int comm_option_both=2;
    //@}

    /// \name The default command objects
    //@{
    comm_option_s c_commands;
    comm_option_s c_help;
    comm_option_s c_quit;
    comm_option_s c_exit;
    comm_option_s c_license;
    comm_option_s c_warranty;
    comm_option_s c_set;
    comm_option_s c_get;
    comm_option_s c_run;
    comm_option_s c_no_intro;
    comm_option_s c_alias;
    //@}
    
    /// If true, then sync cli::verbose, with a parameter of the same name
    bool sync_verbose;

    /** \brief If true, allow the user to use ! to execute a shell command 
	(default true)
    */
    bool shell_cmd_allowed;

    /// The prompt (default <tt>"> "</tt>)
    std::string prompt;
    
    /// A one- or two-line description (default is empty string)
    std::string desc;

    /// The name of the command
    std::string cmd_name;

    /// Additional help text for interactive mode (default is empty string)
    std::string addl_help_cmd;

    /// Additional help text for command-line (default is empty string)
    std::string addl_help_cli;

    /// \name Basic operation
    //@{
    /** \brief Add a new command
	
	Each command/option must have either a short form in
	comm_option_s::shrt or a long from in comm_option_s::lng,
	which is unique from the other commands/options already
	present. You cannot add two commands/options with the same
	short form, even if they have different long forms, and vice
	versa.

    */
    int set_comm_option(comm_option_s &ic);

    /** \brief Remove a command with long name \cmd
     */
    void remove_comm_option(std::string cmd);
    
    /// Add a vector containing new commands/options
    template<class vec_t> int set_comm_option_vec
      (size_t list_size, vec_t &option_list) {

      for(size_t k=0;k<list_size;k++) {

	if (option_list[k].lng.length()<2) {
	  std::string str=((std::string)"Long option '")+option_list[k].lng+
	    "' does not have at "+
	    "least two characters in cli::set_comm_option().";
	  O2SCL_ERR(str.c_str(),exc_efailed);
	}
	
	bool found=false;
	for(size_t i=0;found==false && i<clist.size();i++) {
	  // If short or long options match
	  if ((option_list[k].shrt!=0 && 
	       clist[i].shrt==option_list[k].shrt) || 
	      (option_list[k].lng.length()>0 && 
	       clist[i].lng==option_list[k].lng)) {
	    found=true;
	  }
	}
	if (found==true) {
	  // Call the error handler
	  if (option_list[k].shrt!=0) {
	    std::string err="Option ";
	    err+=option_list[k].shrt;
	    err+=((std::string)" , ")+option_list[k].lng+" already present.";
	    O2SCL_ERR(err.c_str(),exc_einval);
	  } else {
	    std::string err="Option ";
	    err+=option_list[k].lng+" already present.";
	    O2SCL_ERR(err.c_str(),exc_einval);
	  }
	}
	// Add the option to the option list
	clist.push_back(option_list[k]);
      }

      return 0;
    }

    /// Set one-line help text for a parameter named \c param
    int set_param_help(std::string param, std::string help);

    /** \brief Automatically parse arguments to main and 
	call interactive mode if required
    */
    int run_auto(int argc, char *argv[], int debug=0);
    //@}
    
    /** \brief The function which obtains input from the user

	\future Think about whether or not this should be protected?
	(Possibly not, as it's extensively used by acolm.cpp)
    */
    virtual char *cli_gets(const char *c);
    
    /// Call functions corresponding to command-line args
    int call_args(std::vector<cmd_line_arg> &ca);

    /** \brief Process command-line arguments from a const char array
	
        This doesn't actually execute the functions for the
	corresponding options, but simply processes the parameters \c
	argv and \c argv and packs the information into \c ca.

	This function assumes that <tt>argc[0]</tt> just contains
	the name of the command, and should thus be ignored.
    */
    int process_args(int argc, char *argv[], 
		     std::vector<cmd_line_arg> &ca, int debug=0,
		     bool also_call_args=false);

    /** \brief Process command-line arguments from a vector of strings
	
        This doesn't actually execute the functions for the
	corresponding options, but simply processes the arguments in
	\c sv and packs the information into \c ca.
    */
    int process_args(std::vector<std::string> &sv,
		     std::vector<cmd_line_arg> &ca, int debug);
    
    /** \brief Process command-line arguments from a string

	\future There's a typecast in this function to (char *)
	from (const char *) which needs reworking.
     */
    int process_args(std::string s, std::vector<cmd_line_arg> &ca, 
		     int debug=0, bool also_call_args=false);

    /** \brief Set verbosity
	
	Most errors are output to the screen even if verbose is zero.
    */
    int set_verbose(int v);

    /// Run the interactive mode
    int run_interactive();

    // Create a new command
    //    int replace_command(comm_option &ic);

    /** \brief Set an alias \c alias for the string \c str

	Aliases can also be set using the command \c 'alias', but 
	that version allows only one-word aliases. 
     */
    int set_alias(std::string alias, std::string str);

    /** \brief Set an alias \c alias for the string \c str

	Aliases can also be set using the command \c 'alias', but 
	that version allows only one-word aliases. 
     */
    std::string get_alias(std::string alias);
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
