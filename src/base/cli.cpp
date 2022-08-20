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
#include <cstdlib>
#include <o2scl/cli.h>
#include <o2scl/xml.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/cursesw.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_NEVER_DEFINED
/*
  This code used to print out a short usage string, e.g.
  
  Usage: file [-bcikLhnNsvz] [-f namefile] [-F separator] 
  [-m magicfiles] file...
  
  It's not used in acol and is commented out here but it might be
  useful in the future. It could be added to the run_auto() function
  (i.e. if no parameters to the command are given, print out the
  usage).
*/
string cmd_line::usage() {
  string s;
  s+=cmd_name+" ";
  int pos=cmd_name.length()+1;
  for(size_t is=0;is<opts.size();is++) {
    string tmp;
    if (opts[is].shrt!='\0') {
      if (opts[is].parm_desc.length()>0) {
	tmp=((string)"[-")+opts[is].shrt+"|-"+opts[is].lng+" "+
	  opts[is].parm_desc+"] ";
      } else {
	tmp=((string)"[-")+opts[is].shrt+"|-"+opts[is].lng+"] ";
      }
    } else {
      if (opts[is].parm_desc.length()>0) {
	tmp=((string)"[--")+opts[is].lng+" "+opts[is].parm_desc+"] ";
      } else {
	tmp=((string)"[--")+opts[is].lng+"] ";
      }
    }
    if (pos+tmp.length()>80) {
      s+="\n  "+tmp;
      pos=2+tmp.length();
    } else {
      s+=tmp;
      pos+=tmp.length();
    }
  }
  s+="\n";
  return s;
}

#endif

int cli::parameter_double::set(std::string s) {
  if (parse_strings) {
    
    // This doesn't work because value_spec is in o2scl_hdf, which is
    // currently a separate library.
    //
    // o2scl_hdf::value_spec(s,*d);
    
    int iret=function_to_double_nothrow(s,*d);
    return iret;
    
  } else {
    *d=o2scl::stod(s);
  }
  return 0;
}

int cli::parameter_int::set(std::string s) {
  if (parse_strings) {
    double d;
    int iret=function_to_double_nothrow(s,d);
    if (iret==0) *i=((int)d);
    return iret;
  }
  return o2scl::stoi_nothrow(s,*i);
}

int cli::parameter_size_t::set(std::string st) {
  if (parse_strings) {
    double d;
    int iret=function_to_double_nothrow(st,d);
    if (iret==0) *s=((size_t)d);
    return iret;
  }
  return o2scl::stoszt_nothrow(st,*s);
}

cli::cli() {

  verbose=1;
  sync_verbose=true;
  prompt="> ";
  
  c_help.shrt='h';
  c_help.lng="help";
  c_help.min_parms=0;
  c_help.max_parms=1;
  c_help.desc="Show help information.";
  c_help.help="Show generic help information, or, if an argument is given ";
  c_help.help+="give the documentation for the specified command. ";
  c_help.help+="Note that required arguments are typically given inside ";
  c_help.help+="angled brackes <> while optional arguments are given ";
  c_help.help+="inside square brackets []. If an argument is given ";
  c_help.help+="and it is not a command name but it is a parameter name ";
  c_help.help+="then output the name, value, and description of that ";
  c_help.help+="parameter.";
  c_help.parm_desc="[command name or parameter name]";
  c_help.func=new comm_option_mfptr<cli>(this,&cli::comm_option_help);
  c_help.type=comm_option_both;
  
  c_no_intro.shrt=0;
  c_no_intro.lng="no-intro";
  c_no_intro.min_parms=0;
  c_no_intro.max_parms=0;
  c_no_intro.desc="Do not print introductory text.";
  c_no_intro.help=((string)"If this command is used, then no ")+
    "introductory text (including the program name and the \"no "+
    "warranty\" statement) is output in interactive mode.";
  c_no_intro.parm_desc="";
  c_no_intro.func=new comm_option_mfptr<cli>(this,&cli::comm_option_no_intro);
  c_no_intro.type=comm_option_cl_param;

  c_commands.shrt=0;
  c_commands.lng="commands";
  c_commands.min_parms=0;
  c_commands.max_parms=0;
  c_commands.desc="List all available commands.";
  c_commands.help="This command lists all the available commands";
  c_commands.parm_desc="";
  c_commands.func=new comm_option_mfptr<cli>(this,&cli::comm_option_commands);
  c_commands.type=comm_option_both;

  c_alias.shrt=0;
  c_alias.lng="alias";
  // The value of min_parms has to be two so that it is possible
  // to parse aliases which have definitions which begin with
  // dashes.
  c_alias.min_parms=2;
  c_alias.max_parms=-1;
  c_alias.desc="Create a command alias.";
  c_alias.help=((string)"This creates an alias for a command ")+
    "or arguments or a combination of both. To list all aliases, use "+
    "either \"-alias list list\" or \"alias list list\" in "+
    "interactive mode. In general, aliases are applied only after "+
    "they are defined. However, aliases given anywhere as arguments on "+
    "the command line will be applied to the rest of the command line, "+
    "and if multiple definitions for the same alias are given on the "+
    "command line, then only the last one applies.";
  c_alias.parm_desc="<alias name> <alias definition>";
  c_alias.func=new comm_option_mfptr<cli>(this,&cli::comm_option_alias);
  c_alias.type=comm_option_both;
  
  c_get.shrt=0;
  c_get.lng="get";
  c_get.min_parms=1;
  c_get.max_parms=1;
  c_get.desc="Get the value of a parameter.";
  c_get.help=((string)"This command gets the value of a parameter. ")+
    "A list of parameters can be obtained by e.g. -help get "+
    "on the command-line or \"help get\" in interactive mode.";
  c_get.parm_desc="<parameter name>";
  c_get.func=new comm_option_mfptr<cli>(this,&cli::comm_option_get);
  c_get.type=comm_option_both;

  c_run.shrt=0;
  c_run.lng="run";
  c_run.min_parms=1;
  c_run.max_parms=1;
  c_run.parm_desc="<file name>";
  c_run.desc="Run a file containing a list of commands.";
  c_run.help=((string)"This command executes a series of commands ")+
    "specified in a file, line by line. Quotes can be used "+
    "to specify parameters which have single spaces. All whitespace "+
    "is compressed into a single space. Quotes inside a word with "+
    "no whitespace are ignored.";
  c_run.func=new comm_option_mfptr<cli>(this,&cli::comm_option_run);
  c_run.type=comm_option_both;
      
  c_shell.shrt=0;
  c_shell.lng="shell";
  c_shell.min_parms=1;
  c_shell.max_parms=-1;
  c_shell.parm_desc="<commands>";
  c_shell.desc="Run a shell command.";
  c_shell.help=((string)"This command excutes the specified shell ")+
    "command. If verbose>0 (the default), then the shell output "+
    "is preceeded by \"Executing system command: <command>\" and "+
    "followed by \"Done with system command (returned <return "+
    "value>).\". In interactive mode, shell commands can also "+
    "be run using ! instead of \"shell \", i.e. \"!ls\".";
  c_shell.func=new comm_option_mfptr<cli>(this,&cli::comm_option_shell);
  c_shell.type=comm_option_both;
      
  c_set.shrt=0;
  c_set.lng="set";
  c_set.min_parms=2;
  c_set.max_parms=2;
  c_set.parm_desc="<parameter name> <value>";
  c_set.desc="Set the value of a parameter.";
  c_set.help=((string)"This command sets the value of a parameter. ")+
    "A list of parameters can be obtained by e.g. -help set "+
    "on the command-line or \"help set \" in interactive mode."+
    "For floating point or integer parameters, constant expressions "+
    "can be used, e.g. acos(-1). "+
    "For string-like parameters, quotes can be used "+
    "to specify parameters which have single spaces. All whitespace "+
    "is compressed into a single space. Quotes inside a word with "+
    "no whitespace are ignored.";
  c_set.func=new comm_option_mfptr<cli>(this,&cli::comm_option_set);
  c_set.type=comm_option_both;
      
  c_quit.shrt=0;
  c_quit.lng="quit";
  c_quit.min_parms=0;
  c_quit.max_parms=0;
  c_quit.desc="Quit (synonymous with 'exit').";
  c_quit.help="Quit the current program.";
  c_quit.parm_desc="";
  c_quit.type=comm_option_both;
      
  c_exit.shrt=0;
  c_exit.lng="exit";
  c_exit.min_parms=0;
  c_exit.max_parms=0;
  c_exit.desc="Exit (synonymous with 'quit').";
  c_exit.help="Exit the current program.";
  c_exit.parm_desc="";
  c_exit.type=comm_option_both;
  
  c_license.shrt=0;
  c_license.lng="license";
  c_license.min_parms=0;
  // If license has a parameter, it is the file to output the license to
  c_license.max_parms=1;
  c_license.desc="Show license information.";
  c_license.help=((string)"This command outputs the GPLv3 license ")+
    "information. If the optional [filename] argument is given, then "+
    "the license will be written to the specified file.";
  c_license.parm_desc="[filename]";
  c_license.func=new comm_option_mfptr<cli>(this,&cli::comm_option_license);
  c_license.type=comm_option_both;
      
  c_warranty.shrt=0;
  c_warranty.lng="warranty";
  c_warranty.min_parms=0;
  c_warranty.max_parms=0;
  c_warranty.desc="Show warranty information.";
  c_warranty.help="This command shows the GPLv3 no warranty statement.";
  c_warranty.parm_desc="";
  c_warranty.func=new comm_option_mfptr<cli>(this,&cli::comm_option_warranty);
  c_warranty.type=comm_option_both;
  
  c_xml_to_o2.shrt=0;
  c_xml_to_o2.lng="xml-to-o2";
  c_xml_to_o2.min_parms=0;
  c_xml_to_o2.max_parms=0;
  c_xml_to_o2.desc="Read doxygen XML files to generate run-time documentation";
  c_xml_to_o2.help="X";
  c_xml_to_o2.parm_desc="(no parameters).";
  c_xml_to_o2.func=new comm_option_mfptr<cli>(this,&cli::comm_option_xml_to_o2);
  c_xml_to_o2.type=comm_option_both;
  
  clist.push_back(c_alias);
  clist.push_back(c_commands);
  clist.push_back(c_exit);
  clist.push_back(c_get);
  clist.push_back(c_help);
  clist.push_back(c_license);
  clist.push_back(c_no_intro);
  clist.push_back(c_quit);
  clist.push_back(c_run);
  clist.push_back(c_set);
  clist.push_back(c_shell);
  clist.push_back(c_warranty);
  clist.push_back(c_xml_to_o2);

  user_set_func=0;

  gnu_intro=true;
  
  addl_help_cmd="";
  addl_help_cli="";
  shell_cmd_allowed=true;
}

cli::~cli() {
  delete c_help.func;
  delete c_no_intro.func;
  delete c_commands.func;
  delete c_alias.func;
  delete c_get.func;
  delete c_run.func;
  delete c_shell.func;
  delete c_set.func;
  delete c_license.func;
  delete c_warranty.func;
}

int cli::apply_aliases(vector<string> &sv, size_t istart, bool debug) {

  // No new alias definitions are allowed, this function
  // only makes the relevant replacements

  // Apply each alias in order
  for(al_it it=als.begin();it!=als.end();it++) {

    for(size_t i=istart;i<sv.size();i++) {

      // Skip new alias definitions while applying aliases
      if (sv[i]==((string)"-alias") ||
	  sv[i]==((string)"alias") ||
	  sv[i]==((string)"--alias")) {
	i+=2;
	
      } else if (sv[i]==it->first) {
	
	if (debug) {
	  cout << "Replacing " << sv[i] << " with " << it->second << endl;
	  cout << "Before: " << endl;
	  for(size_t j=0;j<sv.size();j++) {
	    cout << j << ": " << sv[j] << endl;
	  }
	}
	
	// ---------------------------------------------------------
	// Refactor the replacement into an array and insert into sv
	
	// First, split the alias definition into sv3
	vector<std::string> sv3;
	o2scl::split_string(it->second,sv3);
	
	// Then, construct a new vector, sv2, which contains the
	// correct entries in order
	vector<string> sv2;
	for(size_t j=0;j<i;j++) {
	  sv2.push_back(sv[j]);
	}
	for(size_t j=0;j<sv3.size();j++) {
	  sv2.push_back(sv3[j]);
	}
	for(size_t j=i+1;j<sv.size();j++) {
	  sv2.push_back(sv[j]);
	}
	
	if (debug) {
	  cout << "After: " << endl;
	  for(size_t j=0;j<sv2.size();j++) {
	    cout << j << ": " << sv2[j] << endl;
	  }
	}
	
	// Set sv equal to the new vector
	sv=sv2;
      }
    }
  }
  return 0;
}

bool cli::string_equal_dash(string s1, string s2) {
  for(size_t i=0;i<s1.size();i++) {
    if (s1[i]=='-') s1[i]='_';
  }
  for(size_t i=0;i<s2.size();i++) {
    if (s2[i]=='-') s2[i]='_';
  }
  if (s1==s2) return true;
  return false;
}

int cli::set_verbose(int v) {
  verbose=v;
  return 0;
}

int cli::comm_option_alias(vector<string> &sv, bool itive_com) {

  // This function is only a placeholder. Aliases are parsed separately
  // by process_args().
  
  return 0;
}

int cli::comm_option_no_intro(vector<string> &sv, bool itive_com) {
  gnu_intro=false;
  return 0;
}

int cli::comm_option_commands(vector<string> &sv, bool itive_com) {

  terminal ter;

  // Form a list of the commands
  std::string *slist=new string[clist.size()];
  for(size_t i=0;i<clist.size();i++) {
    slist[i]=ter.bold()+ter.cyan_fg()+clist[i].lng+ter.default_fg();
  }

  // Sort
  vector_sort<string *,string>(clist.size(),slist);

  // Reformat using screenify
  int nout;
  vector<string> outs;
  screenify<string *>(clist.size(),slist,outs);
  nout=outs.size();

  // Output the command list
  //cout << "Command list:\n" << endl;
  for(int i=0;i<nout;i++) {
    cout << outs[i] << endl;
  }
  cout << endl;

  return 0;
}

int cli::process_args_str(string s, vector<cmd_line_arg> &ca, 
			  int debug, bool also_call_args) {
  
  // Reformat string s into the (argc,argv) format
  vector<string> sv;
  split_string(s,sv);

  int ret=process_args(sv,ca,debug,also_call_args);

  return ret;
}

int cli::parse_for_aliases(std::vector<std::string> &svsv,
			   bool allow_undashed) {
  
  for (int c2=0;c2<((int)svsv.size());c2++) {
    if ((svsv[c2]==((string)"alias") || svsv[c2]==((string)"-alias") ||
	 svsv[c2]==((string)"--alias")) && c2+2<((int)svsv.size())) {
      if (svsv[c2+1]==((string)"list") &&
	  svsv[c2+2]==((string)"list")) {
	cout << "Current aliases: " << endl;
	for (al_it it=als.begin();it!=als.end();it++) {
	  cout << it->first << " = " << it->second << endl;
	}
      } else {
	if (als.find(svsv[c2+1])==als.end()) {
	  // Add alias
	  als.insert(std::make_pair(svsv[c2+1],svsv[c2+2]));
	  if (verbose>1) {
	    cout << "New alias \"" << svsv[c2+1] << "\" = \"" << svsv[c2+2]
		 << "\"" << endl;
	  }
	} else {
	  if (verbose>1) {
	    cout << "Redefining alias \"" << svsv[c2+1]
		 << "\" from \"" << als[svsv[c2+1]]
		 << "\" to \"" << svsv[c2+2] << "\"" << endl;
	  }
	  als[svsv[c2+1]]=svsv[c2+2];
	}
      }
      c2+=2;
    }
  }
  return 0;
}

int cli::process_args(std::vector<std::string> &svsv,
		      std::vector<cmd_line_arg> &ca, int debug,
		      bool also_call_args) {
  
  int retval=0;
  
  // Temporary storage for a command-line argument
  cmd_line_arg c;
  
  if (svsv.size()==0) {
    if (debug>0) cout << "No arguments. Returning." << endl;
    return 0;
  }

  // Index of current argument
  int current=0;

  // Take care of aliases. We don't want this to be recursive, so we
  // loop through to look for new alias definitions, and then
  // afterwards we make all of the alias replacements

  // Apply all new alias definitions in \c svsv
  parse_for_aliases(svsv,false);

  // Now apply aliases to \c svsv
  apply_aliases(svsv,current);
  
  bool done=false;
  while(done==false) {

    string s=svsv[current];

    if (debug>0) cout << "Processing: " << s << endl;
    
    // If it's an option, it begins with '-'. Otherwise, just skip it
    if (s[0]=='-') {

      // The corresponding index in clist
      int option_index=0;
      bool found=false;
      
      // Find the corresponding short option
      if (s.length()==2) {

	if (debug>0) cout << "Short option." << endl;

	for(size_t i=0;found==false && i<clist.size();i++) {
	  if (clist[i].shrt==s[1] && 
	      clist[i].type!=comm_option_command) {

	    if (debug>0) cout << "Found option letter: " << s[1] << endl;
	    
	    c.arg=s;
	    c.is_option=true;
	    c.is_valid=true;
	    c.parms.clear();
	    c.cop=&clist[i];

	    found=true;
	    option_index=i;
	  }
	}

      } else {

	// Find the corresponding long option

	if (debug>0) cout << "Long option." << endl;

	// Remove the leading dash
	string s2;
	if (s[1]=='-') s2=s.substr(2,s.length()-2);
	else s2=s.substr(1,s.length()-1);
	if (debug>0) cout << "Looking for: " << s2 << endl;

	for(size_t i=0;found==false && i<clist.size();i++) {
	  if (string_equal_dash(clist[i].lng,s2) && 
	      clist[i].type!=comm_option_command) {
	    if (debug>0) cout << "Found option name: " << s2 << endl;

	    c.arg=s;
	    c.is_option=true;
	    c.is_valid=true;
	    c.parms.clear();
	    c.cop=&clist[i];

	    found=true;
	    option_index=i;
	  }
	}

      }

      if (found==false) {

	string stmp=s.substr(1,s.length()-1);
	par_t it=par_list.find(stmp);
	if (it!=par_list.end()) {

	  // There are no more arguments left, so assume a 'get'
	  if (current+1>=((int)svsv.size())) {

	    bool found_get_cmd=false;
	    for(size_t i=0;found_get_cmd==false && i<clist.size();i++) {

	      if (clist[i].lng=="get") {
		// We've found a get command so arrange the cmd_line_arg
		// object accordingly
		found_get_cmd=true;
		c.arg="get";
		c.is_option=true;
		c.is_valid=true;
		c.parms.clear();
		c.parms.push_back(stmp);
		c.cop=&clist[i];
		current++;
	      }
	    }

	    // No 'get' command was found in the list
	    if (found_get_cmd==false) {
	      
	      if (verbose>0) {
		cerr << "Option '" << stmp << "' not found." << endl;
	      }
	      
	      c.is_option=true;
	      c.is_valid=false;
	      c.arg=stmp;
	      c.parms.clear();
	      c.cop=0;
	      current++;
	      retval=exc_efailed;

	    }

	  } else {
	    
	    // The next string on the command-line
	    string s2=svsv[current+1];

	    // The next string is an argument, so assume 'get'
	    if (s2[0]=='-') {

	      bool found_get_cmd=false;
	      for(size_t i=0;found_get_cmd==false && i<clist.size();i++) {
		
		if (clist[i].lng=="get") {
		  // We've found a get command so arrange the cmd_line_arg
		  // object accordingly
		  found_get_cmd=true;
		  c.arg="get";
		  c.is_option=true;
		  c.is_valid=true;
		  c.parms.clear();
		  c.parms.push_back(stmp);
		  c.cop=&clist[i];
		  current++;
		}
	      }
	      
	      // No 'get' command was found in the list
	      if (found_get_cmd==false) {

		if (verbose>0) {
		  cerr << "Option '" << s << "' not found." << endl;
		}
		
		c.is_option=true;
		c.is_valid=false;
		c.arg=stmp;
		c.parms.clear();
		c.cop=0;
		current++;
		retval=exc_efailed;
	      }

	    } else {

	      // The next string is not an argument, so assume 'set'

	      bool found_set_cmd=false;
	      for(size_t i=0;found_set_cmd==false && i<clist.size();i++) {
		
		if (clist[i].lng=="set") {
		  // We've found a set command so arrange the cmd_line_arg
		  // object accordingly
		  found_set_cmd=true;
		  c.arg="set";
		  c.is_option=true;
		  c.is_valid=true;
		  c.parms.clear();
		  c.parms.push_back(stmp);
		  c.parms.push_back(s2);
		  c.cop=&clist[i];
		  current+=2;
		  cout << "Assuming \"" << s << " " << s2 
		       << "\" implies \"-set " << stmp << " " << s2
		       << "\"" << endl;
		}
	      }
	      
	      // No 'set' command was found in the list
	      if (found_set_cmd==false) {

		if (verbose>0) {
		  cerr << "Option '" << s << "' not found." << endl;
		}
		
		c.is_option=true;
		c.is_valid=false;
		c.arg=stmp;
		c.parms.clear();
		c.cop=0;
		current++;
		retval=exc_efailed;
	      }

	    }

	  }
	  
	  // End of "if (it!=par_list.end())"
	} else {

	  if (verbose>0) {
	    cerr << "Option '" << s << "' not found." << endl;
	  }
	  
	  c.is_option=true;
	  c.is_valid=false;
	  c.arg=s;
	  c.parms.clear();
	  c.cop=0;
	  current++;
	  
	  retval=exc_efailed;

	}
	
      } else {
	  
	// Now, having found the option, look for arguments

	if (debug>0) cout << "Processing possible arguments." << endl;
	
	bool option_done=false;
	
	// Proceed until either we've gone through the maximum number
	// of parameters or the current entry is a new argument
	current++;
	for(int j=0;((j<clist[option_index].max_parms) || 
		     (clist[option_index].max_parms==-1)) && 
	      option_done==false;j++) {
	  
	  // If we're past the end of the argument list, then
	  // we must be done
	  if (current>=((int)svsv.size())) {

	    option_done=true;

	  } else {

	    s=svsv[current];
	    if (debug>0) cout << "Possible argument: " << s << endl;

	    // If we still haven't got the minimum number of parameters
	    // or if the next entry is not an argument
	    if (j<clist[option_index].min_parms-1 || s[0]!='-') {
		
	      if (debug>0) {
		cout << "Adding argument (before min) to list." << endl;
	      }
	      c.parms.push_back(s);
	      current++;

	    } else {

	      // If the entry begins with a dash or and we've got at
	      // least the minimum number of parameters, then we've
	      // got a new option.
	      if (s[0]=='-' && ((j>clist[option_index].min_parms-1) || 
				(clist[option_index].min_parms==-1))) {
		if (debug>0) cout << "Parsed an option" << endl;
		option_done=true;
	      } else {
		// Otherwise, we have a new parameter for the option.
		if (debug>0) {
		  cout << "Adding argument (after min) to list." << endl;
		}
		c.parms.push_back(s);
		current++;
	      }
	    } 
	  }
	}
	
	// Check to see if we have got enough parameters, if not,
	// report an error.
	if (((int)c.parms.size())<clist[option_index].min_parms &&
	    clist[option_index].min_parms!=-1) {
	  if (verbose>0) {
	    cerr << "Too few arguments for option " 
		 << clist[option_index].lng << "." << endl;
	    retval=exc_efailed;
	  }
	  c.is_valid=false;
	}
	
      }
     
      // End of loop "if (s[0]=='-')" loop
    } else {

      // AWS 7/9/20: Previously this was debug>0, but it's surprising
      // that acol, for example, silently skips non-options. This
      // default is good for now, but we could create a new bool class
      // variable which turned off this message if necessary.
      if (verbose>0) {
	cout << "Argument '" << s << "' is not an option. Skipping." << endl;
      }
      
      c.is_option=false;
      c.is_valid=true;
      c.arg=s;
      c.parms.clear();
      c.cop=0;
      
      current++;

    }

    // Add argument to the list of arguments
    ca.push_back(c);

    if (also_call_args && c.is_option && c.is_valid) {
      vector<string> sv;
      sv.push_back(c.arg);
      for(size_t j=0;j<c.parms.size();j++) {
	sv.push_back(c.parms[j]);
      }
      if (c.arg=="-quit" || c.arg=="-exit" ||
	  c.arg=="--quit" || c.arg=="--exit") {
	return 0;
      }
      (*(c.cop->func))(sv,false);
    }
    
    if (current==((int)svsv.size())) {
      done=true;
    }

  }
  
  if (debug>0) cout << "Done with process_args()." << endl;

  return retval;
}

int cli::process_args_c(int argc, char *argv[], 
			vector<cmd_line_arg> &ca, int debug,
			bool also_call_args) {

  vector<string> sv;
  for(int i=1;i<argc;i++) {
    sv.push_back(argv[i]);
  }
  return process_args(sv,ca,debug,also_call_args);
}

int cli::call_args(vector<cmd_line_arg> &ca, int debug) {
  for(size_t i=0;i<ca.size();i++) {
    if (ca[i].is_option && ca[i].is_valid) {
      vector<string> sv;
      sv.push_back(ca[i].arg);
      for(size_t j=0;j<ca[i].parms.size();j++) {
	sv.push_back(ca[i].parms[j]);
      }
      if (ca[i].arg=="-quit" || ca[i].arg=="-exit" ||
	  ca[i].arg=="--quit" || ca[i].arg=="--exit") {
	if (debug>0) {
	  cout << "In cli::call_args(), quitting." << endl;
	}
	return 0;
      }
      if (debug>0) {
	cout << "In cli::call_args(), calling function for '"
	     << ca[i].arg << "'." << endl;
      }
      (*(ca[i].cop->func))(sv,false);
      if (debug>0) {
	cout << "In cli::call_args(), done calling function for '"
	     << ca[i].arg << "'." << endl;
      }
    }
  }
  return 0;
}

int cli::comm_option_get(vector<string> &sv, bool itive_com) {
  
  if (sv.size()>1) {
    bool found=false;
    for(par_t it=par_list.begin();it!=par_list.end() && found==false;it++) {
      if (string_equal_dash(it->first,sv[1])) {
	cout << "The value of '" << sv[1] << "' is: " << (it->second)->get() 
	     << endl;
	found=true;
      }
    }
    if (found==false) {
      cout << "Parameter named '" << sv[1] << "' not found." << endl;
    }

  } else {
    cout << "You must specify the parameter to get." << endl;
    return exc_efailed;
  }

  return 0;
}

int cli::comm_option_set(vector<string> &sv, bool itive_com) {

  // If a variable and value have been specified
  if (sv.size()>2) {
    
    if (sync_verbose && sv[1]=="verbose") verbose=o2scl::stoi(sv[2]);

    bool found=false;
    for(par_t it=par_list.begin();it!=par_list.end() && found==false;it++) {
      if (string_equal_dash(it->first,sv[1])) {
	(it->second)->set(sv[2]);
	found=true;
      }
    }
    if (found==false) {
      cout << "Parameter named '" << sv[1] << "' not found." << endl;
    }

  } else {
    cerr << "You must specify the parameter to set and it's value." 
	 << endl;
    return exc_efailed;
  }
  
  if (user_set_func!=0) {
    int iret=(*user_set_func)(sv,itive_com);
    if (iret!=success) {
      cerr << "Setting parameter " << sv[1] << " to value "
	   << sv[2] << " failed." << endl;
      return exc_efailed;
    }
  }

  return 0;
}

int cli::output_param_list() {

  terminal ter;

  int ncols_loc;
  int srow, scol;
  int iret=get_screen_size_ioctl(srow,scol);
  if (scol>10 || iret!=0) ncols_loc=scol;
  else ncols_loc=80;
  
  // One row for each parameter and an additional row for
  // the header
  size_t nr=par_list.size()+1;
  
  if (par_list.size()>0) {

    cout << "Parameter list:\n" << endl;

    // Construct a new array of strings, 'tab' containing 
    // the name and value of each parameter
    vector<vector<string> > tab_names(1);
    vector<vector<string> > tab_values(1);
    tab_names[0].resize(nr);
    tab_values[0].resize(nr);
    
    par_t it=par_list.begin();
    tab_names[0][0]="Name";
    tab_values[0][0]="Value";
    for(size_t i=0;i<par_list.size();i++) {
      tab_names[0][i+1]=ter.red_fg()+ter.bold()+it->first+ter.default_fg();
      string stmp=it->second->get();
      // FIXME, AWS 3/28/22, why 64?
      static const size_t nc2=64;
      if (stmp.length()>=nc2) stmp=stmp.substr(0,nc2-3)+"...";
      tab_values[0][i+1]=stmp;
      it++;
    }
    
    // Reformat 'tab' into columns and store in 'tab2'
    columnify c;
    vector<string> col_names(nr);
    vector<string> col_values(nr);
    int align[1]={columnify::align_left};
    c.align(tab_names,1,nr,col_names,align);
    c.align(tab_values,1,nr,col_values,align);
    
    // Output the names, values, and also the help description
    for(size_t i=0;i<nr;i++) {
      if (i==0) {
	// Remove extra space since there are no color charactersx
	cout << col_names[0].substr(0,col_names[0].length()-12)
	     << " " << col_values[i] << endl;
      } else {
	cout << col_names[i] << " " << col_values[i] << endl;
      }
      if (i==0) {
	cout << " Description" << endl;
      } else {

	// If we're on the first parameter, set the iterator
	if (i==1) it=par_list.begin();

        vector<string> desc2;
        rewrap_keep_endlines(it->second->help,desc2,ncols_loc-2);
        for(size_t j=0;j<desc2.size();j++) {
          cout << " " << desc2[j] << endl;
        }

	// Advance the iterator
	it++;
      }

      // Output a trailing endine
      cout << endl;
    }

  }

  return 0;
}

bool cli::is_valid_option(std::string str) {
  // Look for the argument in the command list
  for(size_t i=0;i<clist.size();i++) {
    if (string_equal_dash(clist[i].lng,str)) {
      return true;
    }
  }
  return false;
}

bool cli::is_parameter(std::string name) {
  for(par_t it=par_list.begin();it!=par_list.end();it++) {
    if (string_equal_dash(it->first,name)) {
      return true;
    }
  }
  return false;
}

std::string cli::parameter_desc(std::string name) {
  for(par_t it=par_list.begin();it!=par_list.end();it++) {
    if (string_equal_dash(it->first,name)) {
      return it->second->help;
    }
  }
  O2SCL_ERR("Parameter not found in cli::parameter_desc().",
	    o2scl::exc_einval);
  return "";
}

std::string cli::option_short_desc(std::string name) {
  for(size_t i=0;i<clist.size();i++) {
    if (string_equal_dash(clist[i].lng,name)) {
      return clist[i].desc;
    }
  }
  O2SCL_ERR("Option not found in cli::option_short_desc().",
	    o2scl::exc_einval);
  return "";
}

std::vector<std::string> cli::get_option_list() {
  std::vector<std::string> list;
  for(size_t i=0;i<clist.size();i++) {
    list.push_back(clist[i].lng);
  }
  return list;
}

std::vector<std::string> cli::get_parameter_list() {
  std::vector<std::string> list;
  for(par_t it=par_list.begin();it!=par_list.end();it++) {
    list.push_back(it->first);
  }
  return list;
}

int cli::print_option_list() {

  terminal ter;
  
  cout << "List of command-line options:\n" << endl;
  
  // Repackage the command list
  vector<string> c[3];
  size_t tot=0;
  for(size_t i=0;i<clist.size();i++) {
    if (clist[i].type!=comm_option_command) {
      c[0].push_back(" ");
      c[1].push_back(clist[i].lng);
      c[2].push_back(" ");
      tot++;
    }
  }
  
  // Sort
  sort(c[1].begin(),c[1].end());
  
  // Add the descriptions and the short options
  for(size_t i=0;i<clist.size();i++) {
    for(size_t j=0;j<tot;j++) {
      if (c[1][j]==clist[i].lng) {
	c[2][j]=clist[i].desc;
	if (clist[i].shrt!=0) {
	  c[0][j]=((string)"-")+clist[i].shrt;
	} else {
	  c[0][j]=" ";
	}
      }
    }
  }
  
  // Add the dashes to the long options
  for(size_t j=0;j<tot;j++) {
    c[1][j]=((string)"-")+c[1][j];
  }
  
  size_t maxlen=0;
  for(size_t i=0;i<tot;i++) {
    if (c[1][i].length()>=maxlen) maxlen=c[1][i].length();
  }
  for(size_t i=0;i<tot;i++) {
    if (c[0][i].length()>=2 && c[0][i][0]=='-') {
      cout << "-" << c[0][i][1] << " ";
    } else {
      cout << "   ";
    }
    if (c[1][i].length()>=2 && c[1][i][0]=='-') {
      size_t len=c[1][i].length();
      cout << "-" << ter.cyan_fg() << ter.bold()
	   << c[1][i].substr(1,len-1) << ter.default_fg() << " ";
      for(size_t j=c[1][i].length();j<maxlen;j++) {
	cout << " ";
      }
    } else {
      for(size_t j=c[1][i].length();j<maxlen+1;j++) {
	cout << " ";
      }
    }
    cout << c[2][i] << endl;
  }
  
  return 0;
}

int cli::comm_option_help(vector<string> &sv, bool itive_com) {

  terminal ter;
  
  int ncols_loc;
  int srow, scol;
  int iret=get_screen_size_ioctl(srow,scol);
  if (scol>10 || iret!=0) ncols_loc=scol;
  else ncols_loc=80;
  
  // If not in interactive mode and no arguments are given,
  // Just output the full command list

  if (itive_com==false && sv.size()==1) {

    cout << desc << endl;

    print_option_list();

    // 5/12: I've taken this out because it's too long most of the
    // time. In the future, maybe I should make outputting the
    // parameter list optional here

    // Output parameter list
    // cout << endl;
    //output_param_list();

    if (addl_help_cmd.length()>0) {
      vector<string> str_rewrap;
      rewrap_keep_endlines(addl_help_cmd,str_rewrap,ncols_loc-1);
      for(size_t j=0;j<str_rewrap.size();j++) {
	cout << str_rewrap[j] << endl;
      }
    }

    return 0;
  }

  // The rest of the function handles the 'help' command in
  // interactive mode

  // If this is true, the list of commands will be output
  // at the end
  bool comlist=true;

  // If there's an argument to the help command
  if (sv.size()>1) {

    comlist=false;

    // Look for the argument in the command list
    int ix=-1;
    for(size_t i=0;i<clist.size();i++) {
      if (string_equal_dash(clist[i].lng,sv[1])) {
	ix=i;
	i=clist.size();
      }
    }

    // We can't find the argument
    if (ix==-1) {

      // Look for the argument among get/set parameters
      bool getset_found=false;
      for(par_t it=par_list.begin();it!=par_list.end() &&
	    getset_found==false;it++) {
	
	if (string_equal_dash(it->first,sv[1])) {
	  
	  cout << "Parameter: " << ter.red_fg() << ter.bold() << sv[1]
	       << ter.default_fg() << " value: " << (it->second)->get()
               << endl;

	  vector<string> desc2;
	  rewrap_keep_endlines(((string)"Description: ")+it->second->help,
                               desc2,ncols_loc-1);
	  for(size_t j=0;j<desc2.size();j++) {
	    cout << desc2[j] << endl;
	  }
          
	  getset_found=true;
	}
      }

      if (getset_found==false) {
	cout << "Command '" << sv[1] << "' not found.\n" << endl;
	comlist=true;
      }

    } else {

      // The user has given a command name as an parameter, so 
      // print out usage information for that command

      bool redirected=false;
      if (!isatty(STDOUT_FILENO)) redirected=true;
      
      if (clist[ix].parm_desc.length()==0) {
	cout << "Usage: " << ter.cyan_fg() << ter.bold() 
	     << clist[ix].lng << ter.default_fg()
	     << " (no arguments)" << endl;
      } else {
	cout << "Usage: " << ter.cyan_fg() << ter.bold() 
	     << clist[ix].lng << ter.default_fg()
	     << " " << clist[ix].parm_desc << endl;
      }

      if (clist[ix].desc.length()==0) {
	cout << "(No description.)" << endl;
      } else {
	cout << "Short description: " << clist[ix].desc << endl;
      }

      // Output any additional help text specified
      if (clist[ix].help.length()>0) {

	cout << "\nLong description:" << endl;
	{
	  // Rewrap before writing additional help text
	  vector<string> desc2;
	  rewrap_keep_endlines(clist[ix].help,desc2,ncols_loc-1);
	  for(size_t j=0;j<desc2.size();j++) {
	    cout << desc2[j] << endl;
	  }
	}

      }
      
      // If 'help set' or 'help get' was requested, output
      // a list of the valid parameters
      if ((sv[1]=="set" || sv[1]=="get")) {

	// Output parameter list
	cout << endl;
	output_param_list();

      }


    }

  } 

  // If no argument or an invalid argument was given, just
  // print out the command list
  if (comlist) {

    cout << desc << endl;

    print_option_list();
    
  }

  return 0;
}

int cli::comm_option_shell(vector<string> &sv, bool itive_com) {

  if (shell_cmd_allowed) {
    if (sv.size()==1) {
      cerr << "Not enough arguments to shell." << endl;
      return 1;
    }
    string cmd;
    for(size_t i=1;i<sv.size();i++) {
      if (i==1) {
        cmd+=sv[i];
      } else {
        cmd+=" "+sv[i];
      }
    }
    if (verbose>0) {
      cout << "Executing system command: " 
	   << cmd << endl;
    }
    int sret=system(cmd.c_str());
    if (verbose>0) {
      cout << "Done with system command (returned " 
	   << sret << ")." << endl;
    }
  } else {
    cerr << "Shell commands not allowed." << endl;
    return 1;
  }
  
  return 0;
}

int cli::comm_option_run(vector<string> &sv, bool itive_com) {

  string i1;
  bool suc=false;
  if (sv.size()>=2) {
    i1=sv[1];
    suc=true;
  } else {
    i1=cli_gets("Enter filename to run (or blank to stop): ");
    if (i1.length()>0) {
      suc=true;
    }
  }
  sv.clear();

  if (suc==false) {
    return 0;
  }

  ifstream fin(i1.c_str());
  vector<string> sw;

  string entry;
  while(getline(fin,entry)) {

    if (entry.length()>0 && entry!="exit" && entry!="quit") {

      // Remove up to two preceeding dashes
      if (entry[0]=='-' && entry[1]=='-') {
	entry=entry.substr(2,entry.length()-2);
      } else if (entry[0]=='-') {
	entry=entry.substr(1,entry.length()-1);
      }
      
      split_string(entry,sw);

      if (sw[0]==((string)"alias")) {
	parse_for_aliases(sw,true);
      } else {
	// Apply any aliases
	for(al_it it=als.begin();it!=als.end();it++) {
	  apply_aliases(sw,0);
	}
      }

      if (sw[0][0]=='!') {

	if (shell_cmd_allowed) {
	  if (verbose>0) {
	    cout << "> " << entry << endl;
	  }
	  entry=entry.substr(1,entry.length()-1);
	  if (verbose>0) {
	    cout << cmd_name << ": Executing system command: " 
		 << entry << endl;
	  }
	  int sret=system(entry.c_str());
	  if (verbose>0) {
	    cout << cmd_name << ": Done with system command (returned " 
		 << sret << ")." << endl;
	  }
	}
	  
      } else if (sw[0][0]!='#') {

	bool found=false;
	
	for(size_t i=0;i<clist.size();i++) {
	  
	  // Find the command
	  if (string_equal_dash(clist[i].lng,sw[0]) && 
	      clist[i].type!=comm_option_cl_param) {

	    if (verbose>0) {
	      cout << "> " << entry << endl;
	    }
	    
	    // If there is an associated function
	    if (clist[i].func!=0) {

	      found=true;
	  
	      // Run the appropriate function
	      (*(clist[i].func))(sw,true);
	    
	      // Force exiting from the loop
	      i=clist.size();
	    }
	  }
	}
	
	// If the command wasn't found, then interpret it
	// as a request to report the value of a variable

	if (found==false && entry!="exit" && entry!="quit") {
	  suc=false;
	  
	  if (sw.size()==1) {

	    // Check if the user gave a parameter
	    par_t it=par_list.find(sw[0]);
	    if (it!=par_list.end()) {
	      suc=true;
	      cout << "The value of '" << sw[0] 
		   << "' is: " << (it->second)->get() 
		   << endl;
	    }

	    // If get fails, it could be just a typo, not a bad get 
	    // command, so we ignore it here and just leave 'suc'
	    // to false
	    
	  } else {

	    par_t it=par_list.find(sw[0]);
	    if (it!=par_list.end()) {
	      suc=true;
	      // We output here even if verbose is 0, because we want
	      // to make sure the user isn't surprised by the fact
	      // that they're setting a variable without using
	      // 'set'. To ensure no output, the user should use the
	      // full form with the 'set' command.
	      cout << "Setting variable '" << sw[0] 
		   << "' to " << sw[1] << endl;
	      
	      it->second->set(sw[1]);

	      sw.insert(sw.begin(),"set");
	      if (user_set_func!=0) {
		(*user_set_func)(sw,itive_com);
	      }

	    }

	    // If there is no parameter with that name, it could be
	    // just a typo so we ignore it here and just leave 'suc'
	    // to false

	  }
	}
	
	if (found==false && entry!="exit" && entry!="quit" && suc==false) {
	  cout << "Command '" << sw[0] << "' not found." << endl;
	}


      } else {

	// Just copy the comments to output
	cout << entry << endl;

      }

      // Clear for the next command from the file
      sw.clear();

    }

  }

  return 0;
}

int cli::run_interactive() {

  bool suc=false;
  
  if (gnu_intro) {
    cout << desc << endl;
    cout << "This program comes with ABSOLUTELY NO WARRANTY; for details\n"
	 << "type `warranty'.  This is free software, and you are welcome\n"
	 << "to redistribute it under certain conditions; type `license'\n"
	 << "for details." << endl;
    cout << "------------------------------------------------------------"
	 << endl;
  }
  
  string entry;
  vector<string> sv;
    
  do {
      
    entry=cli_gets(prompt.c_str());

    if (entry[0]=='!') {
      
#ifndef O2SCL_NO_SYSTEM_FUNC
      
      if (shell_cmd_allowed) {
	entry=entry.substr(1,entry.length()-1);
	if (verbose>0) {
	  cout << cmd_name << ": Executing system command: " 
	       << entry << endl;
	}
	int sret=system(entry.c_str());
	if (verbose>0) {
	  cout << cmd_name << ": Done with system command (returned " 
	       << sret << ")." << endl;
	}
	//return 0;

      } else {
	
	cout << "Shell commands disabled by shell_cmd_allowed=false." << endl;
	return exc_efailed;
	
      }

#else
      
      cout << "Shell commands disabled by O2SCL_NO_SYSTEM_FUNC." << endl;
      return exc_efailed;

#endif

    } else if (entry.length()>0) {

      // Remove up to two preceeding dashes
      if (entry[0]=='-' && entry[1]=='-') {
	entry=entry.substr(2,entry.length()-2);
      } else if (entry[0]=='-') {
	entry=entry.substr(1,entry.length()-1);
      }
      
      split_string(entry,sv);

      if (sv[0]==((string)"alias")) {
	parse_for_aliases(sv,true);
      } else {
	// Apply any aliases
	apply_aliases(sv,0);
      }
	  
      if (sv[0][0]!='#') {

	bool found=false;
      
	for(size_t i=0;i<clist.size();i++) {

	  // Find the command
	  if (string_equal_dash(clist[i].lng,sv[0]) && 
	      clist[i].type!=comm_option_cl_param) {
	  
	    // If there is an associated function
	    if (sv[0]!="quit" && sv[0]!="exit" && clist[i].func!=0) {

	      found=true;

	      // Double check that the number of parameters is correct
	      if (clist[i].min_parms!=-1 && 
		  clist[i].min_parms>((int)(sv.size()))-1) {
		cerr << "Minimum number of parameters is " 
		     << clist[i].min_parms << " but only " 
		     << ((int)sv.size())-1
		     << " parameters were given." << endl;
		vector<string> hc;
		hc.push_back("help");
		hc.push_back(clist[i].lng);
		comm_option_help(hc,true);
	      } else if (clist[i].max_parms!=-1 && 
			 clist[i].max_parms<((int)sv.size())-1) {
		cerr << "Maximum number of parameters is " 
		     << clist[i].max_parms << " and " << ((int)sv.size())-1
		     << " parameters were given." << endl;
		vector<string> hc;
		hc.push_back("help");
		hc.push_back(clist[i].lng);
		comm_option_help(hc,true);
	      } else {
	  
		// Run the appropriate function
		(*(clist[i].func))(sv,true);
		
	      }
	    
	      // Force exiting from the loop
	      i=clist.size();
	    }
	  }
	}

	if (found==false && entry!="exit" && entry!="quit") {
	  suc=false;

	  // If the command was not found in the list, interpret it as
	  // either gettting or setting a parameter and report the value

	  if (sv.size()==1) {
	    
	    // Check if the user gave a parameter
	    par_t it=par_list.find(sv[0]);
	    if (it!=par_list.end()) {
	      suc=true;
	      cout << "The value of '" << sv[0] 
		   << "' is: " << (it->second)->get() 
		   << endl;
	    }

	    // If get fails, it could be just a typo, not a bad get 
	    // command, so we ignore it here and just leave 'suc'
	    // to false
	    
	  } else {

	    par_t it=par_list.find(sv[0]);
	    if (it!=par_list.end()) {
	      suc=true;

	      // We output here even if verbose is 0, because we want
	      // to make sure the user isn't surprised by the fact
	      // that they're setting a variable without using
	      // 'set'. To ensure no output, the user should use the
	      // full form with the 'set' command.
	      cout << "Setting variable '" << sv[0] 
		   << "' to " << sv[1] << endl;

	      it->second->set(sv[1]);

	      sv.insert(sv.begin(),"set");
	      if (user_set_func!=0) {
		(*user_set_func)(sv,true);
	      }
	    }

	    // If there is no parameter with that name, it could be
	    // just a typo so we ignore it here and just leave 'suc'
	    // to false

	  }
	  
	}
	
	if (found==false && entry!="exit" && entry!="quit" && suc==false) {
	  cout << "Command '" << sv[0] << "' not found." << endl;
	}

      }

      sv.clear();

    }
      
  } while(entry!="exit" && entry!="quit");
  
  cout << endl;
  
  return 0;
}

void cli::remove_comm_option(std::string cmd) {
  std::vector<comm_option_s>::iterator it;
  
  for(it=clist.begin();it!=clist.end();it++) {
    if (it->lng==cmd) {
      clist.erase(it);
      return;
    }
  }
  string str=((string)"Option ")+cmd+" not found in "+
    "cli::remove_comm_option().";
  O2SCL_ERR(str.c_str(),o2scl::exc_einval);
  return;
}  

int cli::set_comm_option(comm_option_s &ic) {
  if (ic.lng.length()<2) {
    string str=((string)"Long option '")+ic.lng+"' does not have at "+
      "least two characters in cli::set_comm_option().";
    O2SCL_ERR(str.c_str(),exc_efailed);
  }
  bool found=false;
  size_t i_found=0;
  for(size_t i=0;found==false && i<clist.size();i++) {
    if ((ic.shrt!=0 && clist[i].shrt==ic.shrt) || 
	(ic.lng.length()>0 && clist[i].lng==ic.lng)) {
      found=true;
      i_found=i;
    }
  }
  if (found==true) {
    std::string err="Cannot add option ";
    err+=ic.lng+' ';
    if (ic.shrt!=0) {
      err+="('"+ic.shrt+((std::string)"') ");
    }
    err+="because option "+clist[i_found].lng+" ";
    if (clist[i_found].shrt!=0) {
      err+="('"+clist[i_found].shrt+((std::string)"') ");
    }
    err+="already present in cli::set_comm_option().";
    O2SCL_ERR(err.c_str(),exc_einval);
  }
  clist.push_back(ic);
  return 0;
}

int cli::run_auto(int argc, char *argv[], int debug) {

  int ret;

  std::vector<cmd_line_arg> ca;
  
  // ---------------------------------------
  // Process command-line options
  
  ret=process_args_c(argc,argv,ca,debug);
  if (ret!=0) {
    O2SCL_ERR("Failed to process command-line in cli::run_auto().",
              exc_efailed);
  }
  
  if (ca.size()<1) {
    
    ret=run_interactive();
    if (ret!=0) {
      O2SCL_ERR2("Interactive mode failed in ",
		 "cli::run_auto().",exc_efailed);
    }
    
  } else {
    
    ret=call_args(ca,debug);
    if (ret!=0) {
      O2SCL_ERR("Function call_args() failed in cli::run_auto().",
                exc_efailed);
    }
    
  }

  return success;
}

char *cli::cli_gets(const char *c) {
  std::cout << c << std::flush;
  std::cin.getline(buf,300);
  return buf;
}

int cli::set_alias(std::string alias, std::string str) {
  als.insert(std::make_pair(alias,str));
  return 0;
}

std::string cli::get_alias(std::string alias) {
  return als.find(alias)->second;
}

/*
  int cli::replace_command(comm_option &ic) {
  bool found=false;
  for(size_t i=0;found==false && i<clist.size();i++) {
  if ((ic.shrt!=0 && clist[i].shrt==ic.shrt) || 
  (ic.lng.length()>0 && clist[i].lng==ic.lng)) {
  clist[i]=&ic;
  found=true;
  }
  }
  if (found==false) {
  O2SCL_ERR((((string)"Option ")+ic.shrt+" , "+ic.lng+
  " not present.").c_str(),exc_einval);
  }
  return 0;
  }
*/

int cli::set_param_help(string param, string helps) {
  ph_name.push_back(param);
  ph_desc.push_back(helps);
  return 0;
}

void cli::xml_replacements(std::string &s) {

  terminal ter;
  
  // Make the manual replacements from the 'subs' list
  // above
  for(size_t i=0;i<xml_subs.size();i+=2) {
    string_replace(s,xml_subs[i],xml_subs[i+1]);
  }
  
  // Make the command replacements from the current option
  // list
  std::vector<std::string> comm_list=get_option_list();
  for(size_t i=0;i<comm_list.size();i++) {
    string_replace(s,"<computeroutput> "+comm_list[i]+
                   " </computeroutput>",
                   ter.cyan_fg()+ter.bold()+comm_list[i]+
                   ter.default_fg());
  }
  
  // Make the command replacements from the current parameter list
  for(par_t it=par_list.begin();it!=par_list.end();it++) {
    string_replace(s,"<computeroutput> "+it->first+
                   " </computeroutput>",
                   ter.red_fg()+ter.bold()+it->first+
                   ter.default_fg());
  }
  
  string_replace(s,"  "," ");
  string_replace(s," )",")");
  string_replace(s," ,",",");
  string_replace(s," .",".");
                  
  return;
}

int cli::comm_option_xml_to_o2(vector<string> &sv, bool itive_com) {

#ifdef O2SCL_PUGIXML

  terminal ter;
  
  vector<vector<std::string>> cmd_doc_strings, param_doc_strings;

  // XML walkers
  ostream_walker ow;
  vec_string_walker vsw;
  vsw.indent=false;

  std::string stmp;

  // Go through all the commands
  if (verbose>1) {
    cout << "In cli::comm_option_xml_to_o2(): " << clist.size()
         << " commands." << endl;
  }
  for(size_t j=0;j<clist.size();j++) {

    // This command name, the brief description, the parameter description,
    // and then all of the paragraphs in the long-form documentation,
    // in that order.
    vector<std::string> vs_tmp;

    if (clist[j].doc_type==1) {
      
      pugi::xml_document doc;
      pugi::xml_document doc2;
      
      if (verbose>1) {
        cout << "cmd,fn: " << clist[j].lng << " "
             << clist[j].doc_name << endl;
      }
      
      std::string fn=clist[j].doc_xml_file;

      if (fn.length()>0) {
      
      pugi::xml_node n3=doxygen_xml_member_get
        (fn,clist[j].doc_class,clist[j].doc_name,"briefdescription",doc);

      if (n3!=0) {
      
        // We found a brief description, so add the command name
        // to the list
        vs_tmp.push_back(clist[j].lng);
        
        if (verbose>2) {
          cout << "brief_desc name,value: "
               << n3.name() << " " << n3.value() << endl;
          n3.traverse(ow);
        }

        // Store the XML for the brief description in vsw.output
        n3.traverse(vsw);

        // Combine with spaces, and remove the outer paragraph
        stmp="";
        for(size_t k=0;k<vsw.output.size();k++) {
          if (stmp.length()!=0) {
            stmp+=' ';
          }
          if (vsw.output[k]!=((string)"<para>") &&
              vsw.output[k]!=((string)"</para>")) {
            stmp+=vsw.output[k];
          }
        }
        
        xml_replacements(stmp);
        
        // Add brief description to stmp
        vs_tmp.push_back(stmp);

        // Look for a detailed description
        pugi::xml_node n4=doxygen_xml_member_get
          (fn,clist[j].doc_class,clist[j].doc_name,"detaileddescription",
           doc2);
        
        if (n4!=0) {
        
          if (verbose>2) {
            cout << "desc name,value: "
                 << n4.name() << " " << n4.value() << endl;
            n4.traverse(ow);
          }
        
          bool found=false;

          // Store the XML for the detailed description in vsw.output
          n4.traverse(vsw);
          
          bool done=false;
          stmp="";
          
          for(size_t k=0;k<vsw.output.size() && done==false;k++) {
            
            if (vsw.output[k].find("End of runtime documentation.")!=
                string::npos) {
              done=true;
            } else {
              if (vsw.output[k]!=((string)"<para>")) {
                if (vsw.output[k]==((string)"</para>")) {
                  if (verbose>1) {
                    //cout << "stmp: " << stmp << endl;
                  }
                  found=true;

                  xml_replacements(stmp);

                  vs_tmp.push_back(stmp);
                  stmp.clear();
                } else {
                  if (stmp.length()==0) stmp=vsw.output[k];
                  else stmp+=' '+vsw.output[k];
                }
              }
            }
          }

          // End of if (n4!=0) 
        }
        
        if (vs_tmp.size()>=2) {
          if (verbose>1) {
            for(size_t jj=0;jj<vs_tmp.size();jj++) {
              cout << jj << ": " << vs_tmp[jj] << endl;
            }
            cout << endl;
          } else if (verbose>0) {
            cout << "Adding documention for command: " << vs_tmp[0]
                 << endl;
          }
          cmd_doc_strings.push_back(vs_tmp);
        }

        
        // End of if (n3!=0) 
      }
      
      }
    }

  }

  // Go through all the parameters
  for(par_t itp=par_list.begin();itp!=par_list.end();itp++) {

    // This parameter name and then all of the paragraphs in the
    // long-form description, in that order.
    vector<std::string> vs_tmp;
    
    pugi::xml_document doc;
    pugi::xml_document doc2;
    
    if (verbose>1) {
      cout << "parameter name, doc_name: " << itp->first << " "
           << itp->second->doc_name << endl;
    }

    std::string fn=itp->second->doc_xml_file;
    
    if (fn.length()>0) {
    
    pugi::xml_node n3=doxygen_xml_member_get
      (fn,itp->second->doc_class,itp->second->doc_name,
       "briefdescription",doc);
    if (n3==0) {
      cout << "n3 0 " << itp->second->doc_class << " "
           << itp->second->doc_name << " " << fn << endl;
    }
    
    if (n3!=0) {
      
      // We found a brief description, so add the parameter name
      // to the list
      vs_tmp.push_back(itp->first);
      
      if (verbose>2) {
        cout << "brief desc. name, value: "
             << n3.name() << " " << n3.value() << endl;
        n3.traverse(ow);
      }
      
      // Store the XML for the brief description in vsw.output
      n3.traverse(vsw);
      
      // Combine with spaces, and remove the outer paragraph
      stmp="";
      for(size_t k=0;k<vsw.output.size();k++) {
        if (stmp.length()!=0) {
          stmp+=' ';
        }
        if (vsw.output[k]!=((string)"<para>") &&
            vsw.output[k]!=((string)"</para>")) {
          stmp+=vsw.output[k];
        }
      }

      xml_replacements(stmp);

      // Remove trailing spaces
      while (stmp.length()>2 && stmp[stmp.length()-1]==' ') {
        stmp=stmp.substr(0,stmp.length()-1);
      }
      
      // Add a period at the end, because this brief description
      // is combined with the detailed description by the
      // cli 'help' command above.
      if (stmp.length()>2 && stmp[stmp.length()-1]!='.') {
        stmp+='.';
      }
      
      // Add brief description to stmp
      vs_tmp.push_back(stmp);

      // Look for a detailed description
      pugi::xml_node n4=doxygen_xml_member_get
        (fn,itp->second->doc_class,itp->second->doc_name,
         "detaileddescription",doc2);
      
      if (n4!=0) {
        
        if (verbose>2) {
          cout << "detailed desc. name, value: "
               << n4.name() << " " << n4.value() << endl;
          n4.traverse(ow);
        }
        
        // Store the XML for the brief description in vsw.output
        n4.traverse(vsw);
        
        bool done=false;
        stmp="";
        
        for(size_t k=0;k<vsw.output.size() && done==false;k++) {
          
          if (vsw.output[k].find("End of runtime documentation.")!=
              string::npos) {
            done=true;
          } else {
            if (vsw.output[k]!=((string)"<para>")) {
              if (vsw.output[k]==((string)"</para>")) {
                if (verbose>1) {
                  //cout << "stmp: " << stmp << endl;
                }
                
                xml_replacements(stmp);
                
                vs_tmp.push_back(stmp);
                stmp.clear();
              } else {
                if (stmp.length()==0) stmp=vsw.output[k];
                else stmp+=' '+vsw.output[k];
              }
            }
          }
        }

        // End of if (n4!=0)
      }
      
      if (vs_tmp.size()>=2) {
        if (verbose>1) {
          for(size_t jj=0;jj<vs_tmp.size();jj++) {
            cout << jj << ": " << vs_tmp[jj] << endl;
          }
          cout << endl;
        } else if (verbose>0) {
          cout << "Adding documentation for parameter " << vs_tmp[0] << endl;
        }
        param_doc_strings.push_back(vs_tmp);
      }

      // End of if (n3!=0)
    }

    }
    
  }
  
  if (cmd_doc_strings.size()>0 || param_doc_strings.size()>0) {
    if (doc_o2_file.length()==1) {
      O2SCL_ERR("No doc_o2_file specified in cli::comm_option_xml_to_o2().",
                o2scl::exc_einval);
    }
    o2scl_hdf::hdf_file hf;
    hf.open_or_create(doc_o2_file);
    if (cmd_doc_strings.size()>0) {
      hf.sets_vec_vec_copy("cmd_doc_strings",cmd_doc_strings);
    }
    if (param_doc_strings.size()>0) {
      hf.sets_vec_vec_copy("param_doc_strings",param_doc_strings);
    }
    hf.close();
    cout << "Created file " << doc_o2_file << endl;
  }
  
#else

  cout << "Pugixml must be enabled to create the documentation strings "
       << "from the doxygen\n XML output." << endl;
  
#endif
  
  return 0;
}

int cli::comm_option_warranty(vector<string> &sv, bool itive_com) {

  std::string stmp=desc;
  stmp+="\n";
  
  stmp+="This program is free software: you ";
  stmp+="can redistribute it and/or modify ";
  stmp+="it under the terms of the GNU General ";
  stmp+="Public License as published by ";
  stmp+="the Free Software Foundation, either ";
  stmp+="version 3 of the License, or ";
  stmp+="(at your option) any later version.\n\n";
  
  stmp+="This program is distributed in the ";
  stmp+="hope that it will be useful, ";
  stmp+="but WITHOUT ANY WARRANTY; without ";
  stmp+="even the implied warranty of ";
  stmp+="MERCHANTABILITY or FITNESS FOR A ";
  stmp+="PARTICULAR PURPOSE. See the ";
  stmp+="GNU General Public License for more details.\n\n";
  
  stmp+="The GNU General Public License can also be obtained using the ";
  stmp+="'license' command.";

  int ncols_loc;
  int srow, scol;
  int iret=get_screen_size_ioctl(srow,scol);
  if (scol>10 || iret!=0) ncols_loc=scol;
  else ncols_loc=80;
  
  vector<string> desc2;
  rewrap_keep_endlines(stmp,desc2,ncols_loc-2);
  for(size_t j=0;j<desc2.size();j++) {
    cout << desc2[j] << endl;
  }
  
  return 0;
}

int cli::comm_option_license(vector<string> &sv, bool itive_com) {
  ostream *outs;
  ofstream fout;
  if (sv.size()>1) {
    fout.open(sv[1].c_str());
    outs=&fout;
  } else {
    outs=&cout;
  }
  
  (*outs) << desc << endl;
  
  (*outs) << "--------------------------------------"
	  << "--------------------------------------\n" 
	  << endl;

  (*outs) << "                    GNU GENERAL PUBLIC "
	  << "LICENSE" << endl;
  (*outs) << "                       Version 3, 29 Ju"
	  << "ne 2007" << endl;
  (*outs) << "" << endl;
  (*outs) << " Copyright (C) 2007 Free Software Found"
	  << "ation, Inc. <http://fsf.org/>" << endl;
  (*outs) << " Everyone is permitted to copy and dist"
	  << "ribute verbatim copies" << endl;
  (*outs) << " of this license document, but changing"
	  << " it is not allowed." << endl;
  (*outs) << "" << endl;
  (*outs) << "                            Preamble" << endl;
  (*outs) << "" << endl;
  (*outs) << "  The GNU General Public License is a f"
	  << "ree, copyleft license for" << endl;
  (*outs) << "software and other kinds of works." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The licenses for most software and ot"
	  << "her practical works are designed" << endl;
  (*outs) << "to take away your freedom to share and "
	  << "change the works.  By contrast," << endl;
  (*outs) << "the GNU General Public License is inten"
	  << "ded to guarantee your freedom to" << endl;
  (*outs) << "share and change all versions of a prog"
	  << "ram--to make sure it remains free" << endl;
  (*outs) << "software for all its users.  We, the Fr"
	  << "ee Software Foundation, use the" << endl;
  (*outs) << "GNU General Public License for most of "
	  << "our software; it applies also to" << endl;
  (*outs) << "any other work released this way by its"
	  << " authors.  You can apply it to" << endl;
  (*outs) << "your programs, too." << endl;
  (*outs) << "" << endl;
  (*outs) << "  When we speak of free software, we ar"
	  << "e referring to freedom, not" << endl;
  (*outs) << "price.  Our General Public Licenses are"
	  << " designed to make sure that you" << endl;
  (*outs) << "have the freedom to distribute copies o"
	  << "f free software (and charge for" << endl;
  (*outs) << "them if you wish), that you receive sou"
	  << "rce code or can get it if you" << endl;
  (*outs) << "want it, that you can change the softwa"
	  << "re or use pieces of it in new" << endl;
  (*outs) << "free programs, and that you know you ca"
	  << "n do these things." << endl;
  (*outs) << "" << endl;
  (*outs) << "  To protect your rights, we need to pr"
	  << "event others from denying you" << endl;
  (*outs) << "these rights or asking you to surrender"
	  << " the rights.  Therefore, you have" << endl;
  (*outs) << "certain responsibilities if you distrib"
	  << "ute copies of the software, or if" << endl;
  (*outs) << "you modify it: responsibilities to resp"
	  << "ect the freedom of others." << endl;
  (*outs) << "" << endl;
  (*outs) << "  For example, if you distribute copies"
	  << " of such a program, whether" << endl;
  (*outs) << "gratis or for a fee, you must pass on t"
	  << "o the recipients the same" << endl;
  (*outs) << "freedoms that you received.  You must m"
	  << "ake sure that they, too, receive" << endl;
  (*outs) << "or can get the source code.  And you mu"
	  << "st show them these terms so they" << endl;
  (*outs) << "know their rights." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Developers that use the GNU GPL prote"
	  << "ct your rights with two steps:" << endl;
  (*outs) << "(1) assert copyright on the software, a"
	  << "nd (2) offer you this License" << endl;
  (*outs) << "giving you legal permission to copy, di"
	  << "stribute and/or modify it." << endl;
  (*outs) << "" << endl;
  (*outs) << "  For the developers' and authors' prot"
	  << "ection, the GPL clearly explains" << endl;
  (*outs) << "that there is no warranty for this free"
	  << " software.  For both users' and" << endl;
  (*outs) << "authors' sake, the GPL requires that mo"
	  << "dified versions be marked as" << endl;
  (*outs) << "changed, so that their problems will no"
	  << "t be attributed erroneously to" << endl;
  (*outs) << "authors of previous versions." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Some devices are designed to deny use"
	  << "rs access to install or run" << endl;
  (*outs) << "modified versions of the software insid"
	  << "e them, although the manufacturer" << endl;
  (*outs) << "can do so.  This is fundamentally incom"
	  << "patible with the aim of" << endl;
  (*outs) << "protecting users' freedom to change the"
	  << " software.  The systematic" << endl;
  (*outs) << "pattern of such abuse occurs in the are"
	  << "a of products for individuals to" << endl;
  (*outs) << "use, which is precisely where it is mos"
	  << "t unacceptable.  Therefore, we" << endl;
  (*outs) << "have designed this version of the GPL t"
	  << "o prohibit the practice for those" << endl;
  (*outs) << "products.  If such problems arise subst"
	  << "antially in other domains, we" << endl;
  (*outs) << "stand ready to extend this provision to"
	  << " those domains in future versions" << endl;
  (*outs) << "of the GPL, as needed to protect the fr"
	  << "eedom of users." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Finally, every program is threatened "
	  << "constantly by software patents." << endl;
  (*outs) << "States should not allow patents to rest"
	  << "rict development and use of" << endl;
  (*outs) << "software on general-purpose computers, "
	  << "but in those that do, we wish to" << endl;
  (*outs) << "avoid the special danger that patents a"
	  << "pplied to a free program could" << endl;
  (*outs) << "make it effectively proprietary.  To pr"
	  << "event this, the GPL assures that" << endl;
  (*outs) << "patents cannot be used to render the pr"
	  << "ogram non-free." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The precise terms and conditions for "
	  << "copying, distribution and" << endl;
  (*outs) << "modification follow." << endl;
  (*outs) << "" << endl;
  (*outs) << "                       TERMS AND CONDIT"
	  << "IONS" << endl;
  (*outs) << "" << endl;
  (*outs) << "  0. Definitions." << endl;
  (*outs) << "" << endl;
  (*outs) << "  \"This License\" refers to version 3 "
	  << "of the GNU General Public License." << endl;
  (*outs) << "" << endl;
  (*outs) << "  \"Copyright\" also means copyright-li"
	  << "ke laws that apply to other kinds of" << endl;
  (*outs) << "works, such as semiconductor masks." << endl;
  (*outs) << "" << endl;
  (*outs) << "  \"The Program\" refers to any copyrig"
	  << "htable work licensed under this" << endl;
  (*outs) << "License.  Each licensee is addressed as"
	  << " \"you\".  \"Licensees\" and" << endl;
  (*outs) << "\"recipients\" may be individuals or or"
	  << "ganizations." << endl;
  (*outs) << "" << endl;
  (*outs) << "  To \"modify\" a work means to copy fr"
	  << "om or adapt all or part of the work" << endl;
  (*outs) << "in a fashion requiring copyright permis"
	  << "sion, other than the making of an" << endl;
  (*outs) << "exact copy.  The resulting work is call"
	  << "ed a \"modified version\" of the" << endl;
  (*outs) << "earlier work or a work \"based on\" the"
	  << " earlier work." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A \"covered work\" means either the u"
	  << "nmodified Program or a work based" << endl;
  (*outs) << "on the Program." << endl;
  (*outs) << "" << endl;
  (*outs) << "  To \"propagate\" a work means to do a"
	  << "nything with it that, without" << endl;
  (*outs) << "permission, would make you directly or "
	  << "secondarily liable for" << endl;
  (*outs) << "infringement under applicable copyright"
	  << " law, except executing it on a" << endl;
  (*outs) << "computer or modifying a private copy.  "
	  << "Propagation includes copying," << endl;
  (*outs) << "distribution (with or without modificat"
	  << "ion), making available to the" << endl;
  (*outs) << "public, and in some countries other act"
	  << "ivities as well." << endl;
  (*outs) << "" << endl;
  (*outs) << "  To \"convey\" a work means any kind o"
	  << "f propagation that enables other" << endl;
  (*outs) << "parties to make or receive copies.  Mer"
	  << "e interaction with a user through" << endl;
  (*outs) << "a computer network, with no transfer of"
	  << " a copy, is not conveying." << endl;
  (*outs) << "" << endl;
  (*outs) << "  An interactive user interface display"
	  << "s \"Appropriate Legal Notices\"" << endl;
  (*outs) << "to the extent that it includes a conven"
	  << "ient and prominently visible" << endl;
  (*outs) << "feature that (1) displays an appropriat"
	  << "e copyright notice, and (2)" << endl;
  (*outs) << "tells the user that there is no warrant"
	  << "y for the work (except to the" << endl;
  (*outs) << "extent that warranties are provided), t"
	  << "hat licensees may convey the" << endl;
  (*outs) << "work under this License, and how to vie"
	  << "w a copy of this License.  If" << endl;
  (*outs) << "the interface presents a list of user c"
	  << "ommands or options, such as a" << endl;
  (*outs) << "menu, a prominent item in the list meet"
	  << "s this criterion." << endl;
  (*outs) << "" << endl;
  (*outs) << "  1. Source Code." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The \"source code\" for a work means th"
	  << "e preferred form of the work" << endl;
  (*outs) << "for making modifications to it.  \"Objec"
	  << "t code\" means any non-source" << endl;
  (*outs) << "form of a work." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A \"Standard Interface\" means an int"
	  << "erface that either is an official" << endl;
  (*outs) << "standard defined by a recognized standa"
	  << "rds body, or, in the case of" << endl;
  (*outs) << "interfaces specified for a particular p"
	  << "rogramming language, one that" << endl;
  (*outs) << "is widely used among developers working"
	  << " in that language." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The \"System Libraries\" of an execut"
	  << "able work include anything, other" << endl;
  (*outs) << "than the work as a whole, that (a) is i"
	  << "ncluded in the normal form of" << endl;
  (*outs) << "packaging a Major Component, but which "
	  << "is not part of that Major" << endl;
  (*outs) << "Component, and (b) serves only to enabl"
	  << "e use of the work with that" << endl;
  (*outs) << "Major Component, or to implement a Stan"
	  << "dard Interface for which an" << endl;
  (*outs) << "implementation is available to the publ"
	  << "ic in source code form.  A" << endl;
  (*outs) << "\"Major Component\", in this context, mea"
	  << "ns a major essential component" << endl;
  (*outs) << "(kernel, window system, and so on) of t"
	  << "he specific operating system" << endl;
  (*outs) << "(if any) on which the executable work r"
	  << "uns, or a compiler used to" << endl;
  (*outs) << "produce the work, or an object code int"
	  << "erpreter used to run it." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The \"Corresponding Source\" for a wo"
	  << "rk in object code form means all" << endl;
  (*outs) << "the source code needed to generate, ins"
	  << "tall, and (for an executable" << endl;
  (*outs) << "work) run the object code and to modify"
	  << " the work, including scripts to" << endl;
  (*outs) << "control those activities.  However, it "
	  << "does not include the work's" << endl;
  (*outs) << "System Libraries, or general-purpose to"
	  << "ols or generally available free" << endl;
  (*outs) << "programs which are used unmodified in p"
	  << "erforming those activities but" << endl;
  (*outs) << "which are not part of the work.  For ex"
	  << "ample, Corresponding Source" << endl;
  (*outs) << "includes interface definition files ass"
	  << "ociated with source files for" << endl;
  (*outs) << "the work, and the source code for share"
	  << "d libraries and dynamically" << endl;
  (*outs) << "linked subprograms that the work is spe"
	  << "cifically designed to require," << endl;
  (*outs) << "such as by intimate data communication "
	  << "or control flow between those" << endl;
  (*outs) << "subprograms and other parts of the work"
	  << "." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The Corresponding Source need not inc"
	  << "lude anything that users" << endl;
  (*outs) << "can regenerate automatically from other"
	  << " parts of the Corresponding" << endl;
  (*outs) << "Source." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The Corresponding Source for a work i"
	  << "n source code form is that" << endl;
  (*outs) << "same work." << endl;
  (*outs) << "" << endl;
  (*outs) << "  2. Basic Permissions." << endl;
  (*outs) << "" << endl;
  (*outs) << "  All rights granted under this License"
	  << " are granted for the term of" << endl;
  (*outs) << "copyright on the Program, and are irrev"
	  << "ocable provided the stated" << endl;
  (*outs) << "conditions are met.  This License expli"
	  << "citly affirms your unlimited" << endl;
  (*outs) << "permission to run the unmodified Progra"
	  << "m.  The output from running a" << endl;
  (*outs) << "covered work is covered by this License"
	  << " only if the output, given its" << endl;
  (*outs) << "content, constitutes a covered work.  T"
	  << "his License acknowledges your" << endl;
  (*outs) << "rights of fair use or other equivalent,"
	  << " as provided by copyright law." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may make, run and propagate cover"
	  << "ed works that you do not" << endl;
  (*outs) << "convey, without conditions so long as y"
	  << "our license otherwise remains" << endl;
  (*outs) << "in force.  You may convey covered works"
	  << " to others for the sole purpose" << endl;
  (*outs) << "of having them make modifications exclu"
	  << "sively for you, or provide you" << endl;
  (*outs) << "with facilities for running those works"
	  << ", provided that you comply with" << endl;
  (*outs) << "the terms of this License in conveying "
	  << "all material for which you do" << endl;
  (*outs) << "not control copyright.  Those thus maki"
	  << "ng or running the covered works" << endl;
  (*outs) << "for you must do so exclusively on your "
	  << "behalf, under your direction" << endl;
  (*outs) << "and control, on terms that prohibit the"
	  << "m from making any copies of" << endl;
  (*outs) << "your copyrighted material outside their"
	  << " relationship with you." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Conveying under any other circumstanc"
	  << "es is permitted solely under" << endl;
  (*outs) << "the conditions stated below.  Sublicens"
	  << "ing is not allowed; section 10" << endl;
  (*outs) << "makes it unnecessary." << endl;
  (*outs) << "" << endl;
  (*outs) << "  3. Protecting Users' Legal Rights Fro"
	  << "m Anti-Circumvention Law." << endl;
  (*outs) << "" << endl;
  (*outs) << "  No covered work shall be deemed part "
	  << "of an effective technological" << endl;
  (*outs) << "measure under any applicable law fulfil"
	  << "ling obligations under article" << endl;
  (*outs) << "11 of the WIPO copyright treaty adopted"
	  << " on 20 December 1996, or" << endl;
  (*outs) << "similar laws prohibiting or restricting"
	  << " circumvention of such" << endl;
  (*outs) << "measures." << endl;
  (*outs) << "" << endl;
  (*outs) << "  When you convey a covered work, you w"
	  << "aive any legal power to forbid" << endl;
  (*outs) << "circumvention of technological measures"
	  << " to the extent such circumvention" << endl;
  (*outs) << "is effected by exercising rights under "
	  << "this License with respect to" << endl;
  (*outs) << "the covered work, and you disclaim any "
	  << "intention to limit operation or" << endl;
  (*outs) << "modification of the work as a means of "
	  << "enforcing, against the work's" << endl;
  (*outs) << "users, your or third parties' legal rig"
	  << "hts to forbid circumvention of" << endl;
  (*outs) << "technological measures." << endl;
  (*outs) << "" << endl;
  (*outs) << "  4. Conveying Verbatim Copies." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may convey verbatim copies of the"
	  << " Program's source code as you" << endl;
  (*outs) << "receive it, in any medium, provided tha"
	  << "t you conspicuously and" << endl;
  (*outs) << "appropriately publish on each copy an a"
	  << "ppropriate copyright notice;" << endl;
  (*outs) << "keep intact all notices stating that th"
	  << "is License and any" << endl;
  (*outs) << "non-permissive terms added in accord wi"
	  << "th section 7 apply to the code;" << endl;
  (*outs) << "keep intact all notices of the absence "
	  << "of any warranty; and give all" << endl;
  (*outs) << "recipients a copy of this License along"
	  << " with the Program." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may charge any price or no price "
	  << "for each copy that you convey," << endl;
  (*outs) << "and you may offer support or warranty p"
	  << "rotection for a fee." << endl;
  (*outs) << "" << endl;
  (*outs) << "  5. Conveying Modified Source Versions." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may convey a work based on the Pr"
	  << "ogram, or the modifications to" << endl;
  (*outs) << "produce it from the Program, in the for"
	  << "m of source code under the" << endl;
  (*outs) << "terms of section 4, provided that you a"
	  << "lso meet all of these conditions:" << endl;
  (*outs) << "" << endl;
  (*outs) << "    a) The work must carry prominent no"
	  << "tices stating that you modified" << endl;
  (*outs) << "    it, and giving a relevant date." << endl;
  (*outs) << "" << endl;
  (*outs) << "    b) The work must carry prominent no"
	  << "tices stating that it is" << endl;
  (*outs) << "    released under this License and any"
	  << " conditions added under section" << endl;
  (*outs) << "    7.  This requirement modifies the r"
	  << "equirement in section 4 to" << endl;
  (*outs) << "    \"keep intact all notices\"." << endl;
  (*outs) << "" << endl;
  (*outs) << "    c) You must license the entire work"
	  << ", as a whole, under this" << endl;
  (*outs) << "    License to anyone who comes into po"
	  << "ssession of a copy.  This" << endl;
  (*outs) << "    License will therefore apply, along"
	  << " with any applicable section 7" << endl;
  (*outs) << "    additional terms, to the whole of t"
	  << "he work, and all its parts," << endl;
  (*outs) << "    regardless of how they are packaged"
	  << ".  This License gives no" << endl;
  (*outs) << "    permission to license the work in a"
	  << "ny other way, but it does not" << endl;
  (*outs) << "    invalidate such permission if you h"
	  << "ave separately received it." << endl;
  (*outs) << "" << endl;
  (*outs) << "    d) If the work has interactive user"
	  << " interfaces, each must display" << endl;
  (*outs) << "    Appropriate Legal Notices; however,"
	  << " if the Program has interactive" << endl;
  (*outs) << "    interfaces that do not display Appr"
	  << "opriate Legal Notices, your" << endl;
  (*outs) << "    work need not make them do so." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A compilation of a covered work with "
	  << "other separate and independent" << endl;
  (*outs) << "works, which are not by their nature ex"
	  << "tensions of the covered work," << endl;
  (*outs) << "and which are not combined with it such"
	  << " as to form a larger program," << endl;
  (*outs) << "in or on a volume of a storage or distr"
	  << "ibution medium, is called an" << endl;
  (*outs) << "\"aggregate\" if the compilation and its "
	  << "resulting copyright are not" << endl;
  (*outs) << "used to limit the access or legal right"
	  << "s of the compilation's users" << endl;
  (*outs) << "beyond what the individual works permit"
	  << ".  Inclusion of a covered work" << endl;
  (*outs) << "in an aggregate does not cause this Lic"
	  << "ense to apply to the other" << endl;
  (*outs) << "parts of the aggregate." << endl;
  (*outs) << "" << endl;
  (*outs) << "  6. Conveying Non-Source Forms." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may convey a covered work in obje"
	  << "ct code form under the terms" << endl;
  (*outs) << "of sections 4 and 5, provided that you "
	  << "also convey the" << endl;
  (*outs) << "machine-readable Corresponding Source u"
	  << "nder the terms of this License," << endl;
  (*outs) << "in one of these ways:" << endl;
  (*outs) << "" << endl;
  (*outs) << "    a) Convey the object code in, or em"
	  << "bodied in, a physical product" << endl;
  (*outs) << "    (including a physical distribution "
	  << "medium), accompanied by the" << endl;
  (*outs) << "    Corresponding Source fixed on a dur"
	  << "able physical medium" << endl;
  (*outs) << "    customarily used for software inter"
	  << "change." << endl;
  (*outs) << "" << endl;
  (*outs) << "    b) Convey the object code in, or em"
	  << "bodied in, a physical product" << endl;
  (*outs) << "    (including a physical distribution "
	  << "medium), accompanied by a" << endl;
  (*outs) << "    written offer, valid for at least t"
	  << "hree years and valid for as" << endl;
  (*outs) << "    long as you offer spare parts or cu"
	  << "stomer support for that product" << endl;
  (*outs) << "    model, to give anyone who possesses"
	  << " the object code either (1) a" << endl;
  (*outs) << "    copy of the Corresponding Source fo"
	  << "r all the software in the" << endl;
  (*outs) << "    product that is covered by this Lic"
	  << "ense, on a durable physical" << endl;
  (*outs) << "    medium customarily used for softwar"
	  << "e interchange, for a price no" << endl;
  (*outs) << "    more than your reasonable cost of p"
	  << "hysically performing this" << endl;
  (*outs) << "    conveying of source, or (2) access "
	  << "to copy the" << endl;
  (*outs) << "    Corresponding Source from a network"
	  << " server at no charge." << endl;
  (*outs) << "" << endl;
  (*outs) << "    c) Convey individual copies of the "
	  << "object code with a copy of the" << endl;
  (*outs) << "    written offer to provide the Corres"
	  << "ponding Source.  This" << endl;
  (*outs) << "    alternative is allowed only occasio"
	  << "nally and noncommercially, and" << endl;
  (*outs) << "    only if you received the object cod"
	  << "e with such an offer, in accord" << endl;
  (*outs) << "    with subsection 6b." << endl;
  (*outs) << "" << endl;
  (*outs) << "    d) Convey the object code by offeri"
	  << "ng access from a designated" << endl;
  (*outs) << "    place (gratis or for a charge), and"
	  << " offer equivalent access to the" << endl;
  (*outs) << "    Corresponding Source in the same wa"
	  << "y through the same place at no" << endl;
  (*outs) << "    further charge.  You need not requi"
	  << "re recipients to copy the" << endl;
  (*outs) << "    Corresponding Source along with the"
	  << " object code.  If the place to" << endl;
  (*outs) << "    copy the object code is a network s"
	  << "erver, the Corresponding Source" << endl;
  (*outs) << "    may be on a different server (opera"
	  << "ted by you or a third party)" << endl;
  (*outs) << "    that supports equivalent copying fa"
	  << "cilities, provided you maintain" << endl;
  (*outs) << "    clear directions next to the object"
	  << " code saying where to find the" << endl;
  (*outs) << "    Corresponding Source.  Regardless o"
	  << "f what server hosts the" << endl;
  (*outs) << "    Corresponding Source, you remain ob"
	  << "ligated to ensure that it is" << endl;
  (*outs) << "    available for as long as needed to "
	  << "satisfy these requirements." << endl;
  (*outs) << "" << endl;
  (*outs) << "    e) Convey the object code using pee"
	  << "r-to-peer transmission, provided" << endl;
  (*outs) << "    you inform other peers where the ob"
	  << "ject code and Corresponding" << endl;
  (*outs) << "    Source of the work are being offere"
	  << "d to the general public at no" << endl;
  (*outs) << "    charge under subsection 6d." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A separable portion of the object cod"
	  << "e, whose source code is excluded" << endl;
  (*outs) << "from the Corresponding Source as a Syst"
	  << "em Library, need not be" << endl;
  (*outs) << "included in conveying the object code w"
	  << "ork." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A \"User Product\" is either (1) a \"con"
	  << "sumer product\", which means any" << endl;
  (*outs) << "tangible personal property which is nor"
	  << "mally used for personal, family," << endl;
  (*outs) << "or household purposes, or (2) anything "
	  << "designed or sold for incorporation" << endl;
  (*outs) << "into a dwelling.  In determining whethe"
	  << "r a product is a consumer product," << endl;
  (*outs) << "doubtful cases shall be resolved in fav"
	  << "or of coverage.  For a particular" << endl;
  (*outs) << "product received by a particular user, "
	  << "\"normally used\" refers to a" << endl;
  (*outs) << "typical or common use of that class of "
	  << "product, regardless of the status" << endl;
  (*outs) << "of the particular user or of the way in"
	  << " which the particular user" << endl;
  (*outs) << "actually uses, or expects or is expecte"
	  << "d to use, the product.  A product" << endl;
  (*outs) << "is a consumer product regardless of whe"
	  << "ther the product has substantial" << endl;
  (*outs) << "commercial, industrial or non-consumer "
	  << "uses, unless such uses represent" << endl;
  (*outs) << "the only significant mode of use of the"
	  << " product." << endl;
  (*outs) << "" << endl;
  (*outs) << "  \"Installation Information\" for a User"
	  << " Product means any methods," << endl;
  (*outs) << "procedures, authorization keys, or othe"
	  << "r information required to install" << endl;
  (*outs) << "and execute modified versions of a cove"
	  << "red work in that User Product from" << endl;
  (*outs) << "a modified version of its Corresponding"
	  << " Source.  The information must" << endl;
  (*outs) << "suffice to ensure that the continued fu"
	  << "nctioning of the modified object" << endl;
  (*outs) << "code is in no case prevented or interfe"
	  << "red with solely because" << endl;
  (*outs) << "modification has been made." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If you convey an object code work und"
	  << "er this section in, or with, or" << endl;
  (*outs) << "specifically for use in, a User Product"
	  << ", and the conveying occurs as" << endl;
  (*outs) << "part of a transaction in which the righ"
	  << "t of possession and use of the" << endl;
  (*outs) << "User Product is transferred to the reci"
	  << "pient in perpetuity or for a" << endl;
  (*outs) << "fixed term (regardless of how the trans"
	  << "action is characterized), the" << endl;
  (*outs) << "Corresponding Source conveyed under thi"
	  << "s section must be accompanied" << endl;
  (*outs) << "by the Installation Information.  But t"
	  << "his requirement does not apply" << endl;
  (*outs) << "if neither you nor any third party reta"
	  << "ins the ability to install" << endl;
  (*outs) << "modified object code on the User Produc"
	  << "t (for example, the work has" << endl;
  (*outs) << "been installed in ROM)." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The requirement to provide Installati"
	  << "on Information does not include a" << endl;
  (*outs) << "requirement to continue to provide supp"
	  << "ort service, warranty, or updates" << endl;
  (*outs) << "for a work that has been modified or in"
	  << "stalled by the recipient, or for" << endl;
  (*outs) << "the User Product in which it has been m"
	  << "odified or installed.  Access to a" << endl;
  (*outs) << "network may be denied when the modifica"
	  << "tion itself materially and" << endl;
  (*outs) << "adversely affects the operation of the "
	  << "network or violates the rules and" << endl;
  (*outs) << "protocols for communication across the "
	  << "network." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Corresponding Source conveyed, and In"
	  << "stallation Information provided," << endl;
  (*outs) << "in accord with this section must be in "
	  << "a format that is publicly" << endl;
  (*outs) << "documented (and with an implementation "
	  << "available to the public in" << endl;
  (*outs) << "source code form), and must require no "
	  << "special password or key for" << endl;
  (*outs) << "unpacking, reading or copying." << endl;
  (*outs) << "" << endl;
  (*outs) << "  7. Additional Terms." << endl;
  (*outs) << "" << endl;
  (*outs) << "  \"Additional permissions\" are terms th"
	  << "at supplement the terms of this" << endl;
  (*outs) << "License by making exceptions from one o"
	  << "r more of its conditions." << endl;
  (*outs) << "Additional permissions that are applica"
	  << "ble to the entire Program shall" << endl;
  (*outs) << "be treated as though they were included"
	  << " in this License, to the extent" << endl;
  (*outs) << "that they are valid under applicable la"
	  << "w.  If additional permissions" << endl;
  (*outs) << "apply only to part of the Program, that"
	  << " part may be used separately" << endl;
  (*outs) << "under those permissions, but the entire"
	  << " Program remains governed by" << endl;
  (*outs) << "this License without regard to the addi"
	  << "tional permissions." << endl;
  (*outs) << "" << endl;
  (*outs) << "  When you convey a copy of a covered w"
	  << "ork, you may at your option" << endl;
  (*outs) << "remove any additional permissions from "
	  << "that copy, or from any part of" << endl;
  (*outs) << "it.  (Additional permissions may be wri"
	  << "tten to require their own" << endl;
  (*outs) << "removal in certain cases when you modif"
	  << "y the work.)  You may place" << endl;
  (*outs) << "additional permissions on material, add"
	  << "ed by you to a covered work," << endl;
  (*outs) << "for which you have or can give appropri"
	  << "ate copyright permission." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Notwithstanding any other provision o"
	  << "f this License, for material you" << endl;
  (*outs) << "add to a covered work, you may (if auth"
	  << "orized by the copyright holders of" << endl;
  (*outs) << "that material) supplement the terms of "
	  << "this License with terms:" << endl;
  (*outs) << "" << endl;
  (*outs) << "    a) Disclaiming warranty or limiting"
	  << " liability differently from the" << endl;
  (*outs) << "    terms of sections 15 and 16 of this"
	  << " License; or" << endl;
  (*outs) << "" << endl;
  (*outs) << "    b) Requiring preservation of specif"
	  << "ied reasonable legal notices or" << endl;
  (*outs) << "    author attributions in that materia"
	  << "l or in the Appropriate Legal" << endl;
  (*outs) << "    Notices displayed by works containi"
	  << "ng it; or" << endl;
  (*outs) << "" << endl;
  (*outs) << "    c) Prohibiting misrepresentation of"
	  << " the origin of that material, or" << endl;
  (*outs) << "    requiring that modified versions of"
	  << " such material be marked in" << endl;
  (*outs) << "    reasonable ways as different from t"
	  << "he original version; or" << endl;
  (*outs) << "" << endl;
  (*outs) << "    d) Limiting the use for publicity p"
	  << "urposes of names of licensors or" << endl;
  (*outs) << "    authors of the material; or" << endl;
  (*outs) << "" << endl;
  (*outs) << "    e) Declining to grant rights under "
	  << "trademark law for use of some" << endl;
  (*outs) << "    trade names, trademarks, or service"
	  << " marks; or" << endl;
  (*outs) << "" << endl;
  (*outs) << "    f) Requiring indemnification of lic"
	  << "ensors and authors of that" << endl;
  (*outs) << "    material by anyone who conveys the "
	  << "material (or modified versions of" << endl;
  (*outs) << "    it) with contractual assumptions of"
	  << " liability to the recipient, for" << endl;
  (*outs) << "    any liability that these contractua"
	  << "l assumptions directly impose on" << endl;
  (*outs) << "    those licensors and authors." << endl;
  (*outs) << "" << endl;
  (*outs) << "  All other non-permissive additional t"
	  << "erms are considered \"further" << endl;
  (*outs) << "restrictions\" within the meaning of se"
	  << "ction 10.  If the Program as you" << endl;
  (*outs) << "received it, or any part of it, contain"
	  << "s a notice stating that it is" << endl;
  (*outs) << "governed by this License along with a t"
	  << "erm that is a further" << endl;
  (*outs) << "restriction, you may remove that term. "
	  << " If a license document contains" << endl;
  (*outs) << "a further restriction but permits relic"
	  << "ensing or conveying under this" << endl;
  (*outs) << "License, you may add to a covered work "
	  << "material governed by the terms" << endl;
  (*outs) << "of that license document, provided that"
	  << " the further restriction does" << endl;
  (*outs) << "not survive such relicensing or conveyi"
	  << "ng." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If you add terms to a covered work in"
	  << " accord with this section, you" << endl;
  (*outs) << "must place, in the relevant source file"
	  << "s, a statement of the" << endl;
  (*outs) << "additional terms that apply to those fi"
	  << "les, or a notice indicating" << endl;
  (*outs) << "where to find the applicable terms." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Additional terms, permissive or non-p"
	  << "ermissive, may be stated in the" << endl;
  (*outs) << "form of a separately written license, o"
	  << "r stated as exceptions;" << endl;
  (*outs) << "the above requirements apply either way." << endl;
  (*outs) << "" << endl;
  (*outs) << "  8. Termination." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may not propagate or modify a cov"
	  << "ered work except as expressly" << endl;
  (*outs) << "provided under this License.  Any attem"
	  << "pt otherwise to propagate or" << endl;
  (*outs) << "modify it is void, and will automatical"
	  << "ly terminate your rights under" << endl;
  (*outs) << "this License (including any patent lice"
	  << "nses granted under the third" << endl;
  (*outs) << "paragraph of section 11)." << endl;
  (*outs) << "" << endl;
  (*outs) << "  However, if you cease all violation o"
	  << "f this License, then your" << endl;
  (*outs) << "license from a particular copyright hol"
	  << "der is reinstated (a)" << endl;
  (*outs) << "provisionally, unless and until the cop"
	  << "yright holder explicitly and" << endl;
  (*outs) << "finally terminates your license, and (b"
	  << ") permanently, if the copyright" << endl;
  (*outs) << "holder fails to notify you of the viola"
	  << "tion by some reasonable means" << endl;
  (*outs) << "prior to 60 days after the cessation." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Moreover, your license from a particu"
	  << "lar copyright holder is" << endl;
  (*outs) << "reinstated permanently if the copyright"
	  << " holder notifies you of the" << endl;
  (*outs) << "violation by some reasonable means, thi"
	  << "s is the first time you have" << endl;
  (*outs) << "received notice of violation of this Li"
	  << "cense (for any work) from that" << endl;
  (*outs) << "copyright holder, and you cure the viol"
	  << "ation prior to 30 days after" << endl;
  (*outs) << "your receipt of the notice." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Termination of your rights under this"
	  << " section does not terminate the" << endl;
  (*outs) << "licenses of parties who have received c"
	  << "opies or rights from you under" << endl;
  (*outs) << "this License.  If your rights have been"
	  << " terminated and not permanently" << endl;
  (*outs) << "reinstated, you do not qualify to recei"
	  << "ve new licenses for the same" << endl;
  (*outs) << "material under section 10." << endl;
  (*outs) << "" << endl;
  (*outs) << "  9. Acceptance Not Required for Having"
	  << " Copies." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You are not required to accept this L"
	  << "icense in order to receive or" << endl;
  (*outs) << "run a copy of the Program.  Ancillary p"
	  << "ropagation of a covered work" << endl;
  (*outs) << "occurring solely as a consequence of us"
	  << "ing peer-to-peer transmission" << endl;
  (*outs) << "to receive a copy likewise does not req"
	  << "uire acceptance.  However," << endl;
  (*outs) << "nothing other than this License grants "
	  << "you permission to propagate or" << endl;
  (*outs) << "modify any covered work.  These actions"
	  << " infringe copyright if you do" << endl;
  (*outs) << "not accept this License.  Therefore, by"
	  << " modifying or propagating a" << endl;
  (*outs) << "covered work, you indicate your accepta"
	  << "nce of this License to do so." << endl;
  (*outs) << "" << endl;
  (*outs) << "  10. Automatic Licensing of Downstream"
	  << " Recipients." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Each time you convey a covered work, "
	  << "the recipient automatically" << endl;
  (*outs) << "receives a license from the original li"
	  << "censors, to run, modify and" << endl;
  (*outs) << "propagate that work, subject to this Li"
	  << "cense.  You are not responsible" << endl;
  (*outs) << "for enforcing compliance by third parti"
	  << "es with this License." << endl;
  (*outs) << "" << endl;
  (*outs) << "  An \"entity transaction\" is a transact"
	  << "ion transferring control of an" << endl;
  (*outs) << "organization, or substantially all asse"
	  << "ts of one, or subdividing an" << endl;
  (*outs) << "organization, or merging organizations."
	  << "  If propagation of a covered" << endl;
  (*outs) << "work results from an entity transaction"
	  << ", each party to that" << endl;
  (*outs) << "transaction who receives a copy of the "
	  << "work also receives whatever" << endl;
  (*outs) << "licenses to the work the party's predec"
	  << "essor in interest had or could" << endl;
  (*outs) << "give under the previous paragraph, plus"
	  << " a right to possession of the" << endl;
  (*outs) << "Corresponding Source of the work from t"
	  << "he predecessor in interest, if" << endl;
  (*outs) << "the predecessor has it or can get it wi"
	  << "th reasonable efforts." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You may not impose any further restri"
	  << "ctions on the exercise of the" << endl;
  (*outs) << "rights granted or affirmed under this L"
	  << "icense.  For example, you may" << endl;
  (*outs) << "not impose a license fee, royalty, or o"
	  << "ther charge for exercise of" << endl;
  (*outs) << "rights granted under this License, and "
	  << "you may not initiate litigation" << endl;
  (*outs) << "(including a cross-claim or counterclai"
	  << "m in a lawsuit) alleging that" << endl;
  (*outs) << "any patent claim is infringed by making"
	  << ", using, selling, offering for" << endl;
  (*outs) << "sale, or importing the Program or any p"
	  << "ortion of it." << endl;
  (*outs) << "" << endl;
  (*outs) << "  11. Patents." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A \"contributor\" is a copyright hold"
	  << "er who authorizes use under this" << endl;
  (*outs) << "License of the Program or a work on whi"
	  << "ch the Program is based.  The" << endl;
  (*outs) << "work thus licensed is called the contri"
	  << "butor's \"contributor version\"." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A contributor's \"essential patent cl"
	  << "aims\" are all patent claims" << endl;
  (*outs) << "owned or controlled by the contributor,"
	  << " whether already acquired or" << endl;
  (*outs) << "hereafter acquired, that would be infri"
	  << "nged by some manner, permitted" << endl;
  (*outs) << "by this License, of making, using, or s"
	  << "elling its contributor version," << endl;
  (*outs) << "but do not include claims that would be"
	  << " infringed only as a" << endl;
  (*outs) << "consequence of further modification of "
	  << "the contributor version.  For" << endl;
  (*outs) << "purposes of this definition, \"control\" "
	  << "includes the right to grant" << endl;
  (*outs) << "patent sublicenses in a manner consiste"
	  << "nt with the requirements of" << endl;
  (*outs) << "this License." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Each contributor grants you a non-exc"
	  << "lusive, worldwide, royalty-free" << endl;
  (*outs) << "patent license under the contributor's "
	  << "essential patent claims, to" << endl;
  (*outs) << "make, use, sell, offer for sale, import"
	  << " and otherwise run, modify and" << endl;
  (*outs) << "propagate the contents of its contribut"
	  << "or version." << endl;
  (*outs) << "" << endl;
  (*outs) << "  In the following three paragraphs, a "
	  << "\"patent license\" is any express" << endl;
  (*outs) << "agreement or commitment, however denomi"
	  << "nated, not to enforce a patent" << endl;
  (*outs) << "(such as an express permission to pract"
	  << "ice a patent or covenant not to" << endl;
  (*outs) << "sue for patent infringement).  To \"gran"
	  << "t\" such a patent license to a" << endl;
  (*outs) << "party means to make such an agreement o"
	  << "r commitment not to enforce a" << endl;
  (*outs) << "patent against the party." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If you convey a covered work, knowing"
	  << "ly relying on a patent license," << endl;
  (*outs) << "and the Corresponding Source of the wor"
	  << "k is not available for anyone" << endl;
  (*outs) << "to copy, free of charge and under the t"
	  << "erms of this License, through a" << endl;
  (*outs) << "publicly available network server or ot"
	  << "her readily accessible means," << endl;
  (*outs) << "then you must either (1) cause the Corr"
	  << "esponding Source to be so" << endl;
  (*outs) << "available, or (2) arrange to deprive yo"
	  << "urself of the benefit of the" << endl;
  (*outs) << "patent license for this particular work"
	  << ", or (3) arrange, in a manner" << endl;
  (*outs) << "consistent with the requirements of thi"
	  << "s License, to extend the patent" << endl;
  (*outs) << "license to downstream recipients.  \"Kno"
	  << "wingly relying\" means you have" << endl;
  (*outs) << "actual knowledge that, but for the pate"
	  << "nt license, your conveying the" << endl;
  (*outs) << "covered work in a country, or your reci"
	  << "pient's use of the covered work" << endl;
  (*outs) << "in a country, would infringe one or mor"
	  << "e identifiable patents in that" << endl;
  (*outs) << "country that you have reason to believe"
	  << " are valid." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If, pursuant to or in connection with"
	  << " a single transaction or" << endl;
  (*outs) << "arrangement, you convey, or propagate b"
	  << "y procuring conveyance of, a" << endl;
  (*outs) << "covered work, and grant a patent licens"
	  << "e to some of the parties" << endl;
  (*outs) << "receiving the covered work authorizing "
	  << "them to use, propagate, modify" << endl;
  (*outs) << "or convey a specific copy of the covere"
	  << "d work, then the patent license" << endl;
  (*outs) << "you grant is automatically extended to "
	  << "all recipients of the covered" << endl;
  (*outs) << "work and works based on it." << endl;
  (*outs) << "" << endl;
  (*outs) << "  A patent license is \"discriminatory\" "
	  << "if it does not include within" << endl;
  (*outs) << "the scope of its coverage, prohibits th"
	  << "e exercise of, or is" << endl;
  (*outs) << "conditioned on the non-exercise of one "
	  << "or more of the rights that are" << endl;
  (*outs) << "specifically granted under this License"
	  << ".  You may not convey a covered" << endl;
  (*outs) << "work if you are a party to an arrangeme"
	  << "nt with a third party that is" << endl;
  (*outs) << "in the business of distributing softwar"
	  << "e, under which you make payment" << endl;
  (*outs) << "to the third party based on the extent "
	  << "of your activity of conveying" << endl;
  (*outs) << "the work, and under which the third par"
	  << "ty grants, to any of the" << endl;
  (*outs) << "parties who would receive the covered w"
	  << "ork from you, a discriminatory" << endl;
  (*outs) << "patent license (a) in connection with c"
	  << "opies of the covered work" << endl;
  (*outs) << "conveyed by you (or copies made from th"
	  << "ose copies), or (b) primarily" << endl;
  (*outs) << "for and in connection with specific pro"
	  << "ducts or compilations that" << endl;
  (*outs) << "contain the covered work, unless you en"
	  << "tered into that arrangement," << endl;
  (*outs) << "or that patent license was granted, pri"
	  << "or to 28 March 2007." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Nothing in this License shall be cons"
	  << "trued as excluding or limiting" << endl;
  (*outs) << "any implied license or other defenses t"
	  << "o infringement that may" << endl;
  (*outs) << "otherwise be available to you under app"
	  << "licable patent law." << endl;
  (*outs) << "" << endl;
  (*outs) << "  12. No Surrender of Others' Freedom." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If conditions are imposed on you (whe"
	  << "ther by court order, agreement or" << endl;
  (*outs) << "otherwise) that contradict the conditio"
	  << "ns of this License, they do not" << endl;
  (*outs) << "excuse you from the conditions of this "
	  << "License.  If you cannot convey a" << endl;
  (*outs) << "covered work so as to satisfy simultane"
	  << "ously your obligations under this" << endl;
  (*outs) << "License and any other pertinent obligat"
	  << "ions, then as a consequence you may" << endl;
  (*outs) << "not convey it at all.  For example, if "
	  << "you agree to terms that obligate you" << endl;
  (*outs) << "to collect a royalty for further convey"
	  << "ing from those to whom you convey" << endl;
  (*outs) << "the Program, the only way you could sat"
	  << "isfy both those terms and this" << endl;
  (*outs) << "License would be to refrain entirely fr"
	  << "om conveying the Program." << endl;
  (*outs) << "" << endl;
  (*outs) << "  13. Use with the GNU Affero General P"
	  << "ublic License." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Notwithstanding any other provision o"
	  << "f this License, you have" << endl;
  (*outs) << "permission to link or combine any cover"
	  << "ed work with a work licensed" << endl;
  (*outs) << "under version 3 of the GNU Affero Gener"
	  << "al Public License into a single" << endl;
  (*outs) << "combined work, and to convey the result"
	  << "ing work.  The terms of this" << endl;
  (*outs) << "License will continue to apply to the p"
	  << "art which is the covered work," << endl;
  (*outs) << "but the special requirements of the GNU"
	  << " Affero General Public License," << endl;
  (*outs) << "section 13, concerning interaction thro"
	  << "ugh a network will apply to the" << endl;
  (*outs) << "combination as such." << endl;
  (*outs) << "" << endl;
  (*outs) << "  14. Revised Versions of this License." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The Free Software Foundation may publ"
	  << "ish revised and/or new versions of" << endl;
  (*outs) << "the GNU General Public License from tim"
	  << "e to time.  Such new versions will" << endl;
  (*outs) << "be similar in spirit to the present ver"
	  << "sion, but may differ in detail to" << endl;
  (*outs) << "address new problems or concerns." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Each version is given a distinguishin"
	  << "g version number.  If the" << endl;
  (*outs) << "Program specifies that a certain number"
	  << "ed version of the GNU General" << endl;
  (*outs) << "Public License \"or any later version\" a"
	  << "pplies to it, you have the" << endl;
  (*outs) << "option of following the terms and condi"
	  << "tions either of that numbered" << endl;
  (*outs) << "version or of any later version publish"
	  << "ed by the Free Software" << endl;
  (*outs) << "Foundation.  If the Program does not sp"
	  << "ecify a version number of the" << endl;
  (*outs) << "GNU General Public License, you may cho"
	  << "ose any version ever published" << endl;
  (*outs) << "by the Free Software Foundation." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If the Program specifies that a proxy"
	  << " can decide which future" << endl;
  (*outs) << "versions of the GNU General Public Lice"
	  << "nse can be used, that proxy's" << endl;
  (*outs) << "public statement of acceptance of a ver"
	  << "sion permanently authorizes you" << endl;
  (*outs) << "to choose that version for the Program." << endl;
  (*outs) << "" << endl;
  (*outs) << "  Later license versions may give you a"
	  << "dditional or different" << endl;
  (*outs) << "permissions.  However, no additional ob"
	  << "ligations are imposed on any" << endl;
  (*outs) << "author or copyright holder as a result "
	  << "of your choosing to follow a" << endl;
  (*outs) << "later version." << endl;
  (*outs) << "" << endl;
  (*outs) << "  15. Disclaimer of Warranty." << endl;
  (*outs) << "" << endl;
  (*outs) << "  THERE IS NO WARRANTY FOR THE PROGRAM,"
	  << " TO THE EXTENT PERMITTED BY" << endl;
  (*outs) << "APPLICABLE LAW.  EXCEPT WHEN OTHERWISE "
	  << "STATED IN WRITING THE COPYRIGHT" << endl;
  (*outs) << "HOLDERS AND/OR OTHER PARTIES PROVIDE TH"
	  << "E PROGRAM \"AS IS\" WITHOUT WARRANTY" << endl;
  (*outs) << "OF ANY KIND, EITHER EXPRESSED OR IMPLIE"
	  << "D, INCLUDING, BUT NOT LIMITED TO," << endl;
  (*outs) << "THE IMPLIED WARRANTIES OF MERCHANTABILI"
	  << "TY AND FITNESS FOR A PARTICULAR" << endl;
  (*outs) << "PURPOSE.  THE ENTIRE RISK AS TO THE QUA"
	  << "LITY AND PERFORMANCE OF THE PROGRAM" << endl;
  (*outs) << "IS WITH YOU.  SHOULD THE PROGRAM PROVE "
	  << "DEFECTIVE, YOU ASSUME THE COST OF" << endl;
  (*outs) << "ALL NECESSARY SERVICING, REPAIR OR CORR"
	  << "ECTION." << endl;
  (*outs) << "" << endl;
  (*outs) << "  16. Limitation of Liability." << endl;
  (*outs) << "" << endl;
  (*outs) << "  IN NO EVENT UNLESS REQUIRED BY APPLIC"
	  << "ABLE LAW OR AGREED TO IN WRITING" << endl;
  (*outs) << "WILL ANY COPYRIGHT HOLDER, OR ANY OTHER"
	  << " PARTY WHO MODIFIES AND/OR CONVEYS" << endl;
  (*outs) << "THE PROGRAM AS PERMITTED ABOVE, BE LIAB"
	  << "LE TO YOU FOR DAMAGES, INCLUDING ANY" << endl;
  (*outs) << "GENERAL, SPECIAL, INCIDENTAL OR CONSEQU"
	  << "ENTIAL DAMAGES ARISING OUT OF THE" << endl;
  (*outs) << "USE OR INABILITY TO USE THE PROGRAM (IN"
	  << "CLUDING BUT NOT LIMITED TO LOSS OF" << endl;
  (*outs) << "DATA OR DATA BEING RENDERED INACCURATE "
	  << "OR LOSSES SUSTAINED BY YOU OR THIRD" << endl;
  (*outs) << "PARTIES OR A FAILURE OF THE PROGRAM TO "
	  << "OPERATE WITH ANY OTHER PROGRAMS)," << endl;
  (*outs) << "EVEN IF SUCH HOLDER OR OTHER PARTY HAS "
	  << "BEEN ADVISED OF THE POSSIBILITY OF" << endl;
  (*outs) << "SUCH DAMAGES." << endl;
  (*outs) << "" << endl;
  (*outs) << "  17. Interpretation of Sections 15 and 16." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If the disclaimer of warranty and lim"
	  << "itation of liability provided" << endl;
  (*outs) << "above cannot be given local legal effec"
	  << "t according to their terms," << endl;
  (*outs) << "reviewing courts shall apply local law "
	  << "that most closely approximates" << endl;
  (*outs) << "an absolute waiver of all civil liabili"
	  << "ty in connection with the" << endl;
  (*outs) << "Program, unless a warranty or assumptio"
	  << "n of liability accompanies a" << endl;
  (*outs) << "copy of the Program in return for a fee"
	  << "." << endl;
  (*outs) << "" << endl;
  (*outs) << "                     END OF TERMS AND C"
	  << "ONDITIONS" << endl;
  (*outs) << "" << endl;
  (*outs) << "            How to Apply These Terms to"
	  << " Your New Programs" << endl;
  (*outs) << "" << endl;
  (*outs) << "  If you develop a new program, and you"
	  << " want it to be of the greatest" << endl;
  (*outs) << "possible use to the public, the best wa"
	  << "y to achieve this is to make it" << endl;
  (*outs) << "free software which everyone can redist"
	  << "ribute and change under these terms." << endl;
  (*outs) << "" << endl;
  (*outs) << "  To do so, attach the following notice"
	  << "s to the program.  It is safest" << endl;
  (*outs) << "to attach them to the start of each sou"
	  << "rce file to most effectively" << endl;
  (*outs) << "state the exclusion of warranty; and ea"
	  << "ch file should have at least" << endl;
  (*outs) << "the \"copyright\" line and a pointer to w"
	  << "here the full notice is found." << endl;
  (*outs) << "" << endl;
  (*outs) << "    <one line to give the program's nam"
	  << "e and a brief idea of what it does.>" << endl;
  (*outs) << "    Copyright (C) <year>  <name of author>" << endl;
  (*outs) << "" << endl;
  (*outs) << "    This program is free software: you "
	  << "can redistribute it and/or modify" << endl;
  (*outs) << "    it under the terms of the GNU Gener"
	  << "al Public License as published by" << endl;
  (*outs) << "    the Free Software Foundation, eithe"
	  << "r version 3 of the License, or" << endl;
  (*outs) << "    (at your option) any later version." << endl;
  (*outs) << "" << endl;
  (*outs) << "    This program is distributed in the "
	  << "hope that it will be useful," << endl;
  (*outs) << "    but WITHOUT ANY WARRANTY; without e"
	  << "ven the implied warranty of" << endl;
  (*outs) << "    MERCHANTABILITY or FITNESS FOR A PA"
	  << "RTICULAR PURPOSE.  See the" << endl;
  (*outs) << "    GNU General Public License for more"
	  << " details." << endl;
  (*outs) << "" << endl;
  (*outs) << "    You should have received a copy of "
	  << "the GNU General Public License" << endl;
  (*outs) << "    along with this program.  If not, s"
	  << "ee <http://www.gnu.org/licenses/>." << endl;
  (*outs) << "" << endl;
  (*outs) << "Also add information on how to contact "
	  << "you by electronic and paper mail." << endl;
  (*outs) << "" << endl;
  (*outs) << "  If the program does terminal interact"
	  << "ion, make it output a short" << endl;
  (*outs) << "notice like this when it starts in an i"
	  << "nteractive mode:" << endl;
  (*outs) << "" << endl;
  (*outs) << "    <program>  Copyright (C) <year>  <n"
	  << "ame of author>" << endl;
  (*outs) << "    This program comes with ABSOLUTELY "
	  << "NO WARRANTY; for details type `show w'." << endl;
  (*outs) << "    This is free software, and you are "
	  << "welcome to redistribute it" << endl;
  (*outs) << "    under certain conditions; type `sho"
	  << "w c' for details." << endl;
  (*outs) << "" << endl;
  (*outs) << "The hypothetical commands `show w' and "
	  << "`show c' should show the appropriate" << endl;
  (*outs) << "parts of the General Public License.  O"
	  << "f course, your program's commands" << endl;
  (*outs) << "might be different; for a GUI interface"
	  << ", you would use an \"about box\"." << endl;
  (*outs) << "" << endl;
  (*outs) << "  You should also get your employer (if"
	  << " you work as a programmer) or school," << endl;
  (*outs) << "if any, to sign a \"copyright disclaimer"
	  << "\" for the program, if necessary." << endl;
  (*outs) << "For more information on this, and how t"
	  << "o apply and follow the GNU GPL, see" << endl;
  (*outs) << "<http://www.gnu.org/licenses/>." << endl;
  (*outs) << "" << endl;
  (*outs) << "  The GNU General Public License does n"
	  << "ot permit incorporating your program" << endl;
  (*outs) << "into proprietary programs.  If your pro"
	  << "gram is a subroutine library, you" << endl;
  (*outs) << "may consider it more useful to permit l"
	  << "inking proprietary applications with" << endl;
  (*outs) << "the library.  If this is what you want "
	  << "to do, use the GNU Lesser General" << endl;
  (*outs) << "Public License instead of this License."
	  << "  But first, please read" << endl;
  (*outs) << "<http://www.gnu.org/philosophy/why-not-"
	  << "lgpl.html>.\n" << endl;
  
  (*outs) << "--------------------------------------"
	  << "--------------------------------------\n" 
	  << endl;
  
  if (sv.size()>1) {
    
    fout.close();
    if (verbose>0) cout << "License written to file: " << sv[1] << endl;
    
  } else {
    
    (*outs) << "To output this information to a file, give the filename "
	    << "as the first\nargument of the 'license' command. The file,"
	    << "if already present, will be\noverwritten." << endl;
  }

  return 0;
}

#ifdef NEVER_DEFINED

int cli::comm_option_license(vector<string> &sv, bool itive_com) {

  // For this new version, we need to create a string here of
  // the license text.

  int ncols_loc=80;
  if (sv.size()==1) {
    int srow, scol;
    int iret=get_screen_size_ioctl(srow,scol);
    if (scol>=10 || iret!=0) ncols_loc=scol;
  }
  
  vector<string> desc2;
  rewrap_keep_endlines(stmp,desc2,ncols_loc-2);
  
  ostream *outs;
  ofstream fout;
  if (sv.size()>1) {
    fout.open(sv[1].c_str());
    //for(size_t j=0;j<desc2.size();j++) {
    //fout << desc2[j] << endl;
    //}
    fout << stmp << endl;
    fout.close();
    if (verbose>0) cout << "License written to file: " << sv[1] << endl;
  } else {
    //for(size_t j=0;j<desc2.size();j++) {
    //cout << desc2[j] << endl;
    //}
    cout << stmp << endl;
    cout << "To output this information to a file, give the filename "
         << "as the first\nargument of the 'license' command. The file,"
         << "if already present, will be\noverwritten." << endl;
  }
  
  return 0;
}

#endif

int cli::read_docs() {

  if (file_exists(doc_o2_file)) {
  
    o2scl_hdf::hdf_file hf;
    hf.open(doc_o2_file);
    vector<vector<string>> cmd_doc_strings, param_doc_strings;
    hf.gets_vec_vec_copy("cmd_doc_strings",cmd_doc_strings);
    hf.gets_vec_vec_copy("param_doc_strings",param_doc_strings);
    hf.close();
    
    for(size_t j=0;j<clist.size();j++) {
      bool found=false;
      for(size_t k=0;k<cmd_doc_strings.size() && found==false;k++) {
        if (cmd_doc_strings[k][0]==clist[j].lng) {
          if (cmd_doc_strings[k].size()>=2) {
            clist[j].desc=cmd_doc_strings[k][1];
            if (cmd_doc_strings[k].size()>=3) {
              clist[j].parm_desc=cmd_doc_strings[k][2];
              if (cmd_doc_strings[k].size()>=4) {
                clist[j].help=cmd_doc_strings[k][3];
                for(size_t kk=4;kk<cmd_doc_strings[k].size();kk++) {
                  clist[j].help+="\n\n"+cmd_doc_strings[k][kk];
                }
              }
            }
          }
          found=true;
        }
      }
      if (verbose>2) {
        if (found==true) {
          cout << "Function cli::read_docs() "
               << "found documentation for command "
               << clist[j].lng << " ." << endl;
        } else {
          cout << "Function cli::read_docs() could not "
               << "find documentation for command "
               << clist[j].lng << " ." << endl;
        }
      }
    }
    
    for(par_t it=par_list.begin();it!=par_list.end();it++) {
      bool found=false;
      for(size_t k=0;k<param_doc_strings.size() && found==false;k++) {
        if (param_doc_strings[k][0]==it->first) {
          it->second->help=param_doc_strings[k][1];
          for(size_t kk=2;kk<param_doc_strings[k].size();kk++) {
            it->second->help+="\n\n"+param_doc_strings[k][kk];
          }
          found=true;
        }
      }
      if (verbose>2) {
        if (found) {
          cout << "Function cli::read_docs() found "
               << "documentation for parameter "
               << it->first << " ." << endl;
        } else {
          cout << "Function cli::read_docs() could not find "
               << "documentation for parameter "
               << it->first << " ." << endl;
        }
      }
    }
    
  } else {
    O2SCL_ERR2("File does not exist in ",
               "cli::read_docs().",o2scl::exc_efilenotfound);
  }

  return 0;
}

comm_option_s *cli::get_option_pointer(std::string name) {
  for(size_t i=0;i<clist.size();i++) {
    if (string_equal_dash(clist[i].lng,name)) {
      return &clist[i];
    }
  }
  O2SCL_ERR("Option not found in cli::get_option_pointer().",
            o2scl::exc_esanity);
  return 0;
}
