/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/string_conv.h>
#include <o2scl/err_hnd.h>
#include <o2scl/shunting_yard.h>
#include <o2scl/find_constants.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

std::string o2scl::btos(bool b) {
  if (b) return "1";
  return "0";
}

bool o2scl::has_minus_sign(double *x) {
  string s=dtos(*x);
  if (s[0]=='-') return true;
  return false;
}

string o2scl::itos(int x) {
#ifdef O2SCL_OLDER_COMPILER
  ostringstream strout;
  
  if (strout << x) {
    return strout.str();
  } 
  
  O2SCL_ERR("Conversion from int to string failed in itos(int).",
	    exc_einval);
  return "";
#else
  return std::to_string(x);
#endif
}

string o2scl::szttos(size_t x) {
#ifdef O2SCL_OLDER_COMPILER
  /*
    For some reason I'm having trouble with std::to_string(int) on
    intel compilers, so this is a fallback version
  */
  ostringstream strout;
  
  if (strout << x) {
    return strout.str();
  } 
  
  O2SCL_ERR("Conversion from size_t to string failed in szttos(size_t).",
	    exc_einval);
  return "";
#else
  return std::to_string(x);
#endif
}

string o2scl::ptos(void *p) {
  string ret;
  ostringstream strout;
  strout.setf(ios::scientific);

  if (strout << p) {
    return strout.str();
  }
  
  O2SCL_ERR("Conversion from pointer to string failed in ptos().",
	    exc_einval);
  return "";
}

std::string o2scl::dtos(double x, std::ostream &format) {
  string ret;
  ostringstream strout;
  
  strout.precision(format.precision());
  strout.fill(format.fill());
  strout.flags(format.flags());

  if (strout << x) {
    return strout.str();
  } 
  
  O2SCL_ERR("Conversion from double to string failed in dtos(double,ostream).",
	    exc_einval);
  return "";
}

size_t o2scl::size_of_exponent(double x) {
  string ret;
  ostringstream strout;
  
  strout.setf(ios::scientific);

  if (strout << x) {
    ret=strout.str();
    if (ret.length()==0 || ret.find('e')==std::string::npos) {
      O2SCL_ERR("No exponent found in size_of_exponent().",exc_einval);
      return 0;
    }
    return ret.length()-ret.find('e')-2;
  }

  O2SCL_ERR("Failed to convert to string in size_of_exponent().",exc_einval);
  return 0;
}
  
int o2scl::stoi(string s) {
  return std::stoi(s);
}

int o2scl::stoi_nothrow(string s, int &result) {
  istringstream ins(s);
  if (ins >> result) {
    return 0;
  }
  return exc_einval;
}

size_t o2scl::stoszt(string s) {
  size_t ret;
  istringstream ins(s);
  if (ins >> ret) {
    return ret;
  }
  // AWS, 3/26/19: I've tried including the string parameter 's' in
  // the error message, but it leads to badly formatted error messages
  // so this version is better.
  O2SCL_ERR("Conversion from string to size_t failed in stoszt().",
	    exc_einval);
  return 0;
}

int o2scl::stoszt_nothrow(string s, size_t &result) {
  istringstream ins(s);
  if (ins >> result) {
    return 0;
  }
  return exc_einval;
}

int o2scl::stod_nothrow(string s, double &result) {
  istringstream ins(s);
  if (ins >> result) {
    return 0;
  }
  return exc_einval;
}

bool o2scl::stob(string s, bool err_on_fail) {
  bool ret;
  // Read into a string stream to remove initial whitespace
  istringstream ins(s);
  string s2;
  if (ins >> s2) {
    if (s2.length()>0) {
      if (s2[0]=='t' || s2[0]=='T') return true;
      if (s2[0]>'0' && s2[0]<='9') return true;
    }
    return false;
  }
  // If the read into the istringstream failed
  if (err_on_fail) {
    // AWS, 3/26/19: I've tried including the string parameter 's' in
    // the error message, but it can lead to badly formatted error messages
    // so this version is better.
    O2SCL_ERR("Conversion from string to bool failed in stob().",
	      exc_einval);
  }
  return false;
}

double o2scl::stod(string s) {
  return std::stod(s);
}

bool o2scl::is_number(std::string s) {
  // Number of non-number-like characters
  size_t n_char=0;
  for(size_t i=0;i<s.length();i++) {
    if (s[i]<'0' && s[i]>'9' && s[i]!='+' && s[i]!='-' && s[i]!='.' &&
	s[i]!='e' && s[i]!='E' && s[i]!='d' && s[i]!='D') {
      n_char++;
    }
  }
    
  if (n_char>0) return false;
    
  double ret;
  std::istringstream ins(s);
  if (ins >> ret) return true;
    
  return false;
}

double o2scl::function_to_double(std::string s) {
  // Remove quotes and apostrophes
  for(size_t i=0;i<s.length();i++) {
    if (s[i]=='\"' || s[i]=='\'') {
      string t;
      if (i>0) t+=s.substr(0,i);
      if (i<s.length()-1) t+=s.substr(i+1,s.length()-i-1);
      s=t;
      i=0;
    }
  }
  calculator calc;
  calc.compile(s.c_str(),0);
  double dat=calc.eval(0);
  return dat;
}

double o2scl::find_constant(std::string name, std::string unit) {
  o2scl::find_constants &fc=o2scl_settings.get_find_constants();
  return fc.find_unique(name,unit);
}

int o2scl::function_to_double_nothrow(std::string s, double &result) {
  // Remove quotes and apostrophes
  for(size_t i=0;i<s.length();i++) {
    if (s[i]=='\"' || s[i]=='\'') {
      string t;
      if (i>0) t+=s.substr(0,i);
      if (i<s.length()-1) t+=s.substr(i+1,s.length()-i-1);
      s=t;
      i=0;
    }
  }
  calculator calc;
  int ret=calc.compile_nothrow(s.c_str(),0);
  if (ret!=0) return ret;
  int ret2=calc.eval_nothrow(0,result);
  if (ret2!=0) return ret2;
  return 0;
}

void o2scl::split_string(string str, vector<string> &sv) {

  string tmp, tmp2;
  
  istringstream is(str.c_str());

  while (is >> tmp) {
    
    // If it begins with a quote...
    if (tmp[0]=='\"') {

      // If it also ends with a quote, just remove them both
      if (tmp.length()>=2 && tmp[tmp.length()-1]=='\"') {

	// Remove the initial and final quotes
	tmp2=tmp.substr(1,tmp.size()-2);

	// Copy the reformatted string to the original 'tmp' variable
	tmp=tmp2;

	// Otherwise, look for the next word with a final quote
      } else {
	
	// Remove the initial quote
	if (tmp.length()>=2) {
	  tmp2=tmp.substr(1,tmp.size()-1);
	} else {
	  tmp2=" ";
	}
	
	// Add entries until a final quote is found
	bool done=false;
	while (done==false) {
	  
	  // If there are no more words, or if the next word ends in a
	  // quote, then we're done
	  if (!(is >> tmp)) {
	    done=true;
	  } else if (tmp[tmp.size()-1]=='\"') {
	    if (tmp.length()>=2) {
	      tmp=tmp.substr(0,tmp.size()-1);
	    } else {
	      tmp="";
	    }
	    done=true;
	  }
	  tmp2+=" ";
	  tmp2+=tmp;
	}
	
	// Copy the reformatted string to the original 'tmp' variable
	tmp=tmp2;

      }
    }

    // Add to the list
    sv.push_back(tmp);
  }
  
  return;
}

int o2scl::split_string_delim(string str, vector<string> &list,
			      char delim) {
  
  list.clear();
  size_t k=0;
  while (k<str.length()) {
    size_t loc=str.find(delim,k);
    if (loc!=std::string::npos) {
      std::string stemp=str.substr(k,loc-k);
      list.push_back(stemp);
      k+=stemp.length()+1;
    } else {
      if (k<str.length()) {
	list.push_back(str.substr(k,str.length()-k));
      }
      k=str.length();
    }
  }
  
  return 0;
}

void o2scl::rewrap(std::string str, std::vector<std::string> &sv,
		   size_t ncol) {

  if (sv.size()>0) sv.clear();
  
  std::vector<std::string> sv_tmp;
  split_string(str,sv_tmp);

  string stmp;
  if (sv.size()>0) sv.clear();
  for(size_t old_ix=0;old_ix<sv_tmp.size();old_ix++) {
    if (stmp.length()+sv_tmp[old_ix].length()+1<ncol) {
      if (stmp.length()==0) {
	stmp+=sv_tmp[old_ix];
      } else {
	stmp+=((string)" ")+sv_tmp[old_ix];
      }
    } else {
      sv.push_back(stmp);
      stmp=sv_tmp[old_ix];
    }
  }
  if (stmp.size()>0) {
    sv.push_back(stmp);
  }
  
  return;
}

void o2scl::rewrap_ignore_vt100(std::string str,
				std::vector<std::string> &sv,
				size_t ncol) {

  if (sv.size()>0) sv.clear();
  terminal ter;
  
  std::vector<std::string> sv_tmp;
  split_string(str,sv_tmp);

  string stmp;
  if (sv.size()>0) sv.clear();
  for(size_t old_ix=0;old_ix<sv_tmp.size();old_ix++) {
    if (ter.str_len(stmp)+ter.str_len(sv_tmp[old_ix])+1<ncol) {
      if (ter.str_len(stmp)==0) {
	stmp+=sv_tmp[old_ix];
      } else {
	stmp+=((string)" ")+sv_tmp[old_ix];
      }
    } else {
      sv.push_back(stmp);
      stmp=sv_tmp[old_ix];
    }
  }
  if (stmp.size()>0) {
    sv.push_back(stmp);
  }
  
  return;
}

void o2scl::rewrap_color(std::string str, std::vector<std::string> &sv,
			 size_t ncol) {

  if (sv.size()>0) sv.clear();
  
  std::vector<std::string> sv_tmp;
  split_string(str,sv_tmp);

  terminal ter;
  size_t color_len=(ter.magenta_fg()+ter.bold()+
		    ter.default_fg()).length();
  size_t extra_len=0;

  string stmp;
  if (sv.size()>0) sv.clear();
  for(size_t old_ix=0;old_ix<sv_tmp.size();old_ix++) {
    if (stmp.length()+sv_tmp[old_ix].length()+1<ncol+extra_len) {
      size_t len=sv_tmp[old_ix].length();
      if (stmp.length()==0) {
	if (len>=2 && sv_tmp[old_ix][len-1]==':') {
	  stmp+=ter.magenta_fg()+ter.bold()+
	    sv_tmp[old_ix].substr(0,len-1)+ter.default_fg()+":";
	  extra_len+=color_len;
	} else {
	  if (sv_tmp[old_ix][len-1]==',') {
	    stmp+=ter.cyan_fg()+ter.bold()+
	      sv_tmp[old_ix].substr(0,len-1)+ter.default_fg()+",";
	    extra_len+=color_len;
	  } else {
	    stmp+=ter.cyan_fg()+ter.bold()+
	      sv_tmp[old_ix]+ter.default_fg();
	    extra_len+=color_len;
	  }
	}
      } else {
	stmp+=" ";
	if (sv_tmp[old_ix][len-1]==',') {
	  stmp+=ter.cyan_fg()+ter.bold()+
	    sv_tmp[old_ix].substr(0,len-1)+ter.default_fg()+",";
	  extra_len+=color_len;
	} else {
	  stmp+=ter.cyan_fg()+ter.bold()+
	    sv_tmp[old_ix]+ter.default_fg();
	  extra_len+=color_len;
	}
      }
    } else {
      sv.push_back(stmp);
      stmp=ter.cyan_fg()+ter.bold()+sv_tmp[old_ix]+ter.default_fg();
      extra_len=color_len;
    }
  }
  if (stmp.size()>0) {
    sv.push_back(stmp);
  }
  
  return;
}

void o2scl::rewrap_keep_endlines(std::string str,
				 std::vector<std::string> &sv,
				 size_t ncol, int verbose,
				 bool ignore_vt100) {

  terminal ter;
  
  if (sv.size()>0) sv.clear();

  // First, split with the endline characters to
  // create string list "sv_endlines"
  
  std::vector<std::string> sv_endlines;
  split_string_delim(str,sv_endlines,'\n');
  
  if (verbose>1) {
    std::cout << "sv_endlines: " << sv_endlines.size() << endl;
    for(size_t k=0;k<sv_endlines.size();k++) {
      cout << "Q" << sv_endlines[k] << "Q" << endl;
    }
    
  }

  for(size_t k=0;k<sv_endlines.size();k++) {
    
    if (ter.str_len(sv_endlines[k])<ncol) {

      // If the string is not too wide, even if it's empty,
      // just add it to the final list
      sv.push_back(sv_endlines[k]);
      
    } else {

      // Otherwise, now split this entry in sv_endlines by the space
      // character
      std::vector<std::string> sv_tmp;
      split_string_delim(sv_endlines[k],sv_tmp,' ');

      string stmp;

      // Flag if stmp begins with a space
      bool preceeding_space=false;

      // Iterate through the list
      for(size_t old_ix=0;old_ix<sv_tmp.size();old_ix++) {
	  
	// Proceed the next addition will not make the string too
	// long
	if (ter.str_len(stmp)+ter.str_len(sv_tmp[old_ix])+1<ncol) {
	    
	  // If this entry in svn_tmp is empty, then that's because
	  // there were adjacent spaces at the beginning of
	  // this entry in svn_endlines. Add a space to the
	  // temporary string and flip the flag
	  if (ter.str_len(sv_tmp[old_ix])==0) {
	      
	    preceeding_space=true;
	    stmp+=' ';
	      
	  } else {
	      
	    // Add the next entry to the temporary string,
	    // with or without a space as necessary
	    if (ter.str_len(stmp)==0 || preceeding_space) {
	      stmp+=sv_tmp[old_ix];
	    } else {
	      stmp+=((string)" ")+sv_tmp[old_ix];
	    }
	      
	    // Flip the preceeding space flag since
	    // we're no longer at the beginning of the string
	    preceeding_space=false;
	  }
	    
	} else {
	    
	  // Otherwise, the next addition will make the string too
	  // long, so add the temporary string to sv, reset the
	  // preceeding_space flag, and start the temporary string
	  // with the current entry from sv_tmp
	    
	  sv.push_back(stmp);
	  preceeding_space=false;
	  stmp=sv_tmp[old_ix];
	}
      }

      // Add any remaining text in the temporary string to sv
      if (stmp.size()>0) {
	sv.push_back(stmp);
	preceeding_space=false;
      }
    }
    
  }
  
  return;
}

void o2scl::parse_fortran_format(std::string line,
                                 std::string format,
                                 vector<string> &entries) {

  bool debug=false;
  
  entries.clear();
  vector<string> format_list;
  split_string_delim(format,format_list,',');
  int index=0;
  for(size_t j=0;j<format_list.size();j++) {
    // For entries not ending in 'x'
    int size;
    if (format_list[j][format_list[j].length()-1]!='x') {
      if (format_list[j].find('.')!=std::string::npos) {
        // Get string from second character to the '.'
        size=o2scl::stoi(format_list[j].substr
                         (1,format_list[j].find('.')-1));
      } else {
        // Remove the character at the front and convert
        // to an integer
        size=o2scl::stoi(format_list[j].substr
                         (1,format_list[j].length()-1));
      }
      // If there is at least one character left in the line, then go
      // ahead and add the rest of the string to the entries array
      if (((int)line.length())>index+1) {
        entries.push_back(line.substr(index,size));
      }
    } else {
      // Remove the 'x' at the end and convert to an integer
      size=o2scl::stoi(format_list[j].substr
                       (0,format_list[j].length()-1));
    }
    if (debug) {
      cout << "format_list[j],size: " << format_list[j] << " "
           << size << " " << entries.size() << " " << index << endl;
      char ch;
      cin >> ch;
    }
    index+=size;
  }
  
  return;
}

void o2scl::string_to_char_array(std::string s, char *x, int len) {
  if (((int)s.length())+1>len) {
    cerr << "Not enough space." << endl;
    cerr << s << endl;
    exit(-1);
  }
  remove_whitespace(s);
  for(size_t j=0;j<s.length();j++) {
    x[j]=s[j];
  }
  x[s.length()]='\0';
  return;
}

