/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
  
string o2scl::dtos(double x, int prec, bool auto_prec) {
  string ret;
  ostringstream strout;

  if (prec!=0) {
    if (!auto_prec) strout.setf(ios::scientific);
    strout.precision(prec);
  }

  if (strout << x) {
    return strout.str();
  }
  
  O2SCL_ERR2("Conversion from double to string failed in ",
	     "dtos(double,int,bool).",exc_einval);
  return "";
}

int o2scl::stoi(string s, bool err_on_fail) {
  return std::stoi(s);
}

size_t o2scl::stoszt(string s, bool err_on_fail) {
  size_t ret;
  istringstream ins(s);
  if (ins >> ret) {
    return ret;
  }
  if (err_on_fail) {
    O2SCL_ERR("Conversion from string to size_t failed in stoszt().",
	      exc_einval);
  }
  return 0;
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
    O2SCL_ERR("Conversion from string to bool failed in stob().",
	      exc_einval);
  }
  return false;
}

double o2scl::stod(string s, bool err_on_fail) {
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

double o2scl::function_to_double(std::string s, bool err_on_fail) {
  calculator calc;
  calc.compile(s.c_str(),0);
  double dat=calc.eval(0);
  return dat;
}

void o2scl::split_string(string str, vector<string> &sv) {
  
  string tmp, tmp2;
  
  istringstream is(str.c_str());

  while (is >> tmp) {
    
    // If it begins with a quote...
    if (tmp[0]=='\"') {

      // If it also ends with a quote, just remove them both
      if (tmp[tmp.length()-1]=='\"') {

	// Remove the initial and final quotes
	tmp2=tmp.substr(1,tmp.size()-2);

	// Copy the reformatted string to the original 'tmp' variable
	tmp=tmp2;

	// Otherwise, look for the next word with a final quote
      } else {
	
	// Remove the initial quote
	tmp2=tmp.substr(1,tmp.size()-1);
	
	// Add entries until a final quote is found
	bool done=false;
	while (done==false) {
	  
	  // If there are no more words, or if the next word ends in a
	  // quote, then we're done
	  if (!(is >> tmp)) {
	    done=true;
	  } else if (tmp[tmp.size()-1]=='\"') {
	    tmp=tmp.substr(0,tmp.size()-1);
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

void o2scl::rewrap(std::string str, std::vector<std::string> &sv,
		   size_t ncol) {
  
  std::vector<std::string> sv_tmp;
  split_string(str,sv_tmp);

  //for(size_t i=0;i<sv_tmp.size();i++) {
  //cout << "Phase 1: " << i << " x" << sv_tmp[i] << "x" << endl;
  //}

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
  
  //for(size_t i=0;i<sv.size();i++) {
  //cout << "Phase 2: " << i << " x" << sv[i] << "x" << endl;
  //}
  //exit(-1);

  return;
}
