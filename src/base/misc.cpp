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

// unistd.h is for isatty()
#include <unistd.h>

// For gsl_finite() and gsl_hypot()
#include <gsl/gsl_sys.h>

// For glob()
#include <glob.h>

// For wordexp()
#include <wordexp.h>

// For stat() in file_exists()
#include <sys/stat.h>

#include <boost/math/special_functions/hypot.hpp>

#include <o2scl/misc.h>

using namespace std;
using namespace o2scl;

bool o2scl::file_exists(std::string fname) {
  if (fname.length()==0) return false;
  struct stat buffer;   
  int val=stat(fname.c_str(),&buffer);
  if (val==-1) return false;
  return true;
}     

int o2scl_python_test(int x) {
  return x*x;
}

float o2scl::o2abs(const float x) {
  return fabsf(x);
}

double o2scl::o2abs(const double x) {
  return fabs(x);
}

long double o2scl::o2abs(const long double x) {
  return fabsl(x);
}

bool o2scl::o2isfinite(const double x) {
  return std::isfinite(x);
}

float o2scl::o2hypot(const float x, const float y) {
  return hypotf(x,y);
}

double o2scl::o2hypot(const double x, const double y) {
  return gsl_hypot(x,y);
}

long double o2scl::o2hypot(const long double x, const long double y) {
  return hypotl(x,y);
}

#ifdef O2SCL_LD_TYPES

// cpp_dec_float_35 functions

typedef
boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
cpp_dec_float_35;

cpp_dec_float_35 o2scl::o2abs(const cpp_dec_float_35 x) {
  return boost::multiprecision::abs(x);
}

bool o2scl::o2isfinite(const cpp_dec_float_35 x) {
  return boost::math::isfinite(x);
}

cpp_dec_float_35 o2scl::o2hypot(const cpp_dec_float_35 x,
				const cpp_dec_float_35 y) {
  return boost::math::hypot(x,y);
}

// cpp_dec_float_50 functions

boost::multiprecision::cpp_dec_float_50
o2scl::o2abs(const boost::multiprecision::cpp_dec_float_50 x) {
  return boost::multiprecision::abs(x);
}

bool o2scl::o2isfinite(const boost::multiprecision::cpp_dec_float_50 x) {
  return boost::math::isfinite(x);
}

boost::multiprecision::cpp_dec_float_50
o2scl::o2hypot(const boost::multiprecision::cpp_dec_float_50 x,
	       const boost::multiprecision::cpp_dec_float_50 y) {
  return boost::math::hypot(x,y);
}

// cpp_dec_float_100 functions

boost::multiprecision::cpp_dec_float_100
o2scl::o2abs(const boost::multiprecision::cpp_dec_float_100 x) {
  return boost::multiprecision::abs(x);
}

bool o2scl::o2isfinite(const boost::multiprecision::cpp_dec_float_100 x) {
  return boost::math::isfinite(x);
}

boost::multiprecision::cpp_dec_float_100
o2scl::o2hypot(const boost::multiprecision::cpp_dec_float_100 x,
	       const boost::multiprecision::cpp_dec_float_100 y) {
  return boost::math::hypot(x,y);
}

#endif

int o2scl::pipe_cmd_string(std::string cmd, std::string &result,
			   bool err_on_fail, int nmax) {
  
#ifdef HAVE_POPEN
  
  FILE *ps_pipe=popen(cmd.c_str(),"r");
  if (!ps_pipe) {
    if (err_on_fail) {
      O2SCL_ERR("Pipe could not be opened in o2scl::pipe_cmd_string().",
		o2scl::exc_efailed);
    }
    return 1;
  }
  
  char char_arr[nmax];
  
  // Variable 'cret' is unused, but put here to avoid
  // unused return value errors
  char *cret=fgets(char_arr,nmax,ps_pipe);
  if (cret==0) {
    if (err_on_fail) {
      O2SCL_ERR("Null pointer returned by fgets in o2scl::pipe_cmd_string().",
		o2scl::exc_efailed);
    }
    return 2;
  }
  
  result=char_arr;
  
  int pret=pclose(ps_pipe);
  if (pret!=0) {
    if (err_on_fail) {
      O2SCL_ERR("Close pipe returned non-zero in o2scl::pipe_cmd_string().",
		o2scl::exc_efailed);
    }
    return 4;
  }
  
#else
  
  if (err_on_fail) {
    O2SCL_ERR("Compiled without popen support in o2scl::pipe_cmd_string().",
	      o2scl::exc_efailed);
  }
  return 3;
  
#endif
  
  return 0;
}

int o2scl::python_cmd_string(std::string cmd, std::string &result,
			     bool err_on_fail, int nmax) {
  
#ifdef HAVE_POPEN

  std::string pycmd="python -c "+cmd;
  FILE *ps_pipe=popen(pycmd.c_str(),"r");
  if (!ps_pipe) {
    if (err_on_fail) {
      O2SCL_ERR("Pipe could not be opened in o2scl::python_cmd_string().",
		o2scl::exc_efailed);
    }
    return 1;
  }
  
  char char_arr[nmax];
  
  // Variable 'cret' is unused, but put here to avoid
  // unused return value errors
  char *cret=fgets(char_arr,nmax,ps_pipe);
  if (cret==0) {
    if (err_on_fail) {
      O2SCL_ERR("Null pointer returned by fgets in o2scl::python_cmd_string().",
		o2scl::exc_efailed);
    }
    return 2;
  }
  
  result=char_arr;
  
  int pret=pclose(ps_pipe);
  if (pret!=0) {
    if (err_on_fail) {
      O2SCL_ERR("Close pipe returned non-zero in o2scl::python_cmd_string().",
		o2scl::exc_efailed);
    }
    return 4;
  }
  
#else
  
  if (err_on_fail) {
    O2SCL_ERR("Compiled without popen support in o2scl::python_cmd_string().",
	      o2scl::exc_efailed);
  }
  return 3;
  
#endif
  
  return 0;
}

std::string o2scl::pipe_cmd_string(std::string cmd, int nmax) {
  std::string result;
  pipe_cmd_string(cmd,result,true,nmax);
  return result;
}

std::string o2scl::binary_to_hex(std::string s) {
  std::string t="";
  char nums[16]={'0','1','2','3','4','5','6','7',
		 '8','9','A','B','C','D','E','F'};
  for(size_t i=0;i<s.length();i++) {
    bool found=false;
    if (s[i]=='1' || s[i]=='0') {
      if (i+1<s.length() && (s[i+1]=='1' || s[i+1]=='0')) {
	if (i+2<s.length() && (s[i+2]=='1' || s[i+2]=='0')) {
	  if (i+3<s.length() && (s[i+3]=='1' || s[i+3]=='0')) {
	    found=true;
	    int cnt=0;
	    if (s[i]=='1') cnt+=8;
	    if (s[i+1]=='1') cnt+=4;
	    if (s[i+2]=='1') cnt+=2;
	    if (s[i+3]=='1') cnt++;
	    t+=nums[cnt];
	    i+=3;
	  }
	}
      }
    }
    if (found==false) {
      t+=s[i];
    }
  }
  return t;
}

double o2scl::fermi_function(double E, double mu, double T, double limit) {
  double ret, x=(E-mu)/T;
  
  if (x>limit) {
    ret=0.0;
  } else if (x<-limit) {
    ret=1.0;
  } else {
    ret=1.0/(1.0+exp(x));
  }
  return ret;
}

double o2scl::bose_function(double E, double mu, double T, double limit) {
  double ret, x=(E-mu)/T;
  
  if (x>limit) {
    ret=0.0;
  } else if (x<-limit) {
    ret=-1.0;
  } else if (fabs(x)<1.0e-3) {
    double x2=x*x;
    double x3=x2*x;
    double x5=x3*x2;
    double x7=x5*x2;
    ret=1.0/x-0.5+x/12.0-x3/720.0+x5/30240.0-x7/1209600.0;
  } else {
    ret=1.0/(exp(x)-1.0);
  }
  return ret;
}

size_t o2scl::count_words(string str) {
  string st;
  istringstream *is=new istringstream(str.c_str());
  size_t ctr=0;
  while((*is) >> st) ctr++;
  delete is;

  return ctr;
}

void o2scl::remove_whitespace(std::string &s) {
  // The index 'i' must be an integer so we can do 'i--' below
  for(int i=0;i<((int)s.length());i++) {
    if (s[i]==9 || s[i]==10 || s[i]==11 || s[i]==12 || s[i]==13 
	|| s[i]==32) {
      s=s.substr(0,i)+s.substr(i+1,s.length()-i-1);
      i--;
    }
  }
  return;
}

void o2scl::RGBtoHSV(double r, double g, double b, 
		     double &h, double &s, double &v) {
  
  double min, max, delta;

  if (r<g && r<b) min=r;
  else if (g<b) min=g;
  else min=b;
  if (r>g && r>b) max=r;
  else if (g>b) max=g;
  else max=b;
  v=max;

  delta=max-min;

  if (max!=0) {
    s=delta/max;
  } else {
    // r=g=b=0		
    // s=0,v is undefined
    s=0;
    h=-1;
    return;
  }

  if (r==max) {
    // between yellow & magenta
    h=(g-b)/delta;		
  } else if (g==max) {
    // between cyan & yellow
    h=2+(b-r)/delta;	
  } else {
    // between magenta & cyan
    h=4+(r-g)/delta;	
  }
  
  // degrees
  h*=60;				
  if(h<0) {
    h+=360.0;
  }

  return;
}

/** \brief Convert HSV color to RGB
 */
void o2scl::HSVtoRGB(double h, double s, double v, 
		     double &r, double &g, double &b) {

  int i;
  double f,p,q,t;
  
  if (s==0.0) {
    // achromatic (grey)
    r=g=b=v;
    return;
  }

  if (h==360.0) {
    h=0.0;
  }

  // sector 0 to 5
  h/=60.0;			
  i=((int)floor(h));
  // fractional part of h
  f=h-i;			
  p=v*(1-s);
  q=v*(1-s*f);
  t=v*(1-s*(1-f));

  switch(i) {
  case 0:
    r=v;
    g=t;
    b=p;
    break;
  case 1:
    r=q;
    g=v;
    b=p;
    break;
  case 2:
    r=p;
    g=v;
    b=t;
    break;
  case 3:
    r=p;
    g=q;
    b=v;
    break;
  case 4:
    r=t;
    g=p;
    b=v;
    break;
  default:
    r=v;
    g=p;
    b=q;
    break;
  }

  return;
}

int o2scl::glob_wrapper(std::string pattern,
			std::vector<std::string> &matches) {
  glob_t pglob;
  pglob.gl_offs=0;
  pglob.gl_pathc=0;
  matches.clear();
  int ret=glob(pattern.c_str(),GLOB_MARK | GLOB_TILDE,NULL,&pglob);
  if (ret==0) {
    for(size_t i=0;i<pglob.gl_pathc;i++) {
      matches.push_back(pglob.gl_pathv[i]);
    }
  }
  globfree(&pglob);
  return ret;
}

int o2scl::wordexp_wrapper(std::string word,
			   std::vector<std::string> &matches) {
  wordexp_t pwordexp;
  char **w;
  matches.clear();
  int ret=wordexp(word.c_str(),&pwordexp,0);
  if (ret==0) {
    w=pwordexp.we_wordv;
    for(size_t i=0;i<pwordexp.we_wordc;i++) {
      matches.push_back(w[i]);
    }
  }
  wordfree(&pwordexp);
  return ret;
}

void o2scl::wordexp_single_file(std::string &fname) {
  std::vector<std::string> matches;
  int wret=wordexp_wrapper(fname,matches);
  if (wret!=0) {
    O2SCL_ERR2("Function wordexp_wrapper() failed in ",
	       "wordexp_single_file().",o2scl::exc_einval);
  }
  if (matches.size()>1) {
    O2SCL_ERR2("More than one match found for ",
	       "wordexp_single_file().",o2scl::exc_einval);
  }
  if (matches.size()==0) {
    O2SCL_ERR2("Zero matches in ",
	       "wordexp_single_file().",o2scl::exc_einval);
  }
  fname=matches[0];
  return;
}

terminal::terminal() {
  redirected=false;
  if (!isatty(STDOUT_FILENO)) {
    redirected=true;
  }
}

std::string terminal::hrule(size_t n) {
  std::ostringstream oss;
  if (redirected) {
    for(size_t i=0;i<n;i++) oss << '-';
  } else {
    oss << ((char)27) << '(' << '0';
    for(size_t i=0;i<n;i++) oss << 'q';
    oss << ((char)27) << '(' << 'B';
  }
  return oss.str();
}
  
std::string terminal::bold() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[1m";
  return oss.str();
}

std::string terminal::eight_bit_fg(short col) {
  if (col<0 || col>255) {
    O2SCL_ERR("Color out of range in vt100_eight_bit_fg().",
	      o2scl::exc_efailed);
  }
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[38;5;" << col << "m";
  return oss.str();
}

std::string terminal::three_byte_fg(short red, short green, short blue) {
  if (red<0 || red>255 ||green<0 || green>255 ||blue<0 || blue>255) {
    O2SCL_ERR("Color out of range in vt100_three_byte_bg().",
	      o2scl::exc_efailed);
  }
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[38;2;" << red << ";" << green << ";" << blue << "m";
  return oss.str();
}

std::string terminal::eight_bit_bg(short col) {
  if (col<0 || col>255) {
    O2SCL_ERR("Color out of range in vt100_eight_bit_bg().",
	      o2scl::exc_efailed);
  }
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[48;5;" << col << "m";
  return oss.str();
}

std::string terminal::three_byte_bg(short red, short green, short blue) {
  if (red<0 || red>255 ||green<0 || green>255 ||blue<0 || blue>255) {
    O2SCL_ERR("Color out of range in vt100_three_byte_bg().",
	      o2scl::exc_efailed);
  }
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[48;2;" << red << ";" << green << ";" << blue << "m";
  return oss.str();
}

std::string terminal::lowint() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[2m";
  return oss.str();
}

std::string terminal::underline() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[4m";
  return oss.str();
}

std::string terminal::reverse() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[7m";
  return oss.str();
}

std::string terminal::cyan_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[36m";
  return oss.str();
}

std::string terminal::magenta_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[35m";
  return oss.str();
}

std::string terminal::yellow_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[33m";
  return oss.str();
}

std::string terminal::red_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[31m";
  return oss.str();
}

std::string terminal::green_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[32m";
  return oss.str();
}

std::string terminal::blue_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[34m";
  return oss.str();
}

std::string terminal::cyan_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[46m";
  return oss.str();
}

std::string terminal::magenta_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[45m";
  return oss.str();
}

std::string terminal::yellow_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[43m";
  return oss.str();
}

std::string terminal::red_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[41m";
  return oss.str();
}

std::string terminal::green_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[42m";
  return oss.str();
}

std::string terminal::blue_bg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[44m";
  return oss.str();
}

std::string terminal::default_fg() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "[m";
  return oss.str();
}

std::string terminal::alt_font() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "(0";
  return oss.str();
}

std::string terminal::normal_font() {
  if (redirected) return "";
  std::ostringstream oss;
  oss << ((char)27) << "(B";
  return oss.str();
}

std::string terminal::three_byte_summ_long() {
  if (redirected) return "";
  std::ostringstream oss;
  for(size_t i=0;i<256;i+=17) {
    for(size_t j=0;j<256;j+=17) {
      for(size_t k=0;k<256;k+=17) {
	oss << ((char)27) << "[38;2;" << i << ";" << j << ";" << k << "m";
	if (i<10) {
	  oss << "  " << i;
	} else if (i<100) {
	  oss << " " << i;
	} else {
	  oss << i;
	}
	if (j<10) {
	  oss << "  " << j;
	} else if (j<100) {
	  oss << " " << j;
	} else {
	  oss << j;
	}
	if (k<10) {
	  oss << "  " << k;
	} else if (k<100) {
	  oss << " " << k;
	} else {
	  oss << k;
	}
	oss << ((char)27) << "[m";
	oss << " ";
	oss << ((char)27) << "[48;2;" << i << ";" << j << ";" << k << "m";
	if (i<10) {
	  oss << "  " << i;
	} else if (i<100) {
	  oss << " " << i;
	} else {
	  oss << i;
	}
	if (j<10) {
	  oss << "  " << j;
	} else if (j<100) {
	  oss << " " << j;
	} else {
	  oss << j;
	}
	if (k<10) {
	  oss << "  " << k;
	} else if (k<100) {
	  oss << " " << k;
	} else {
	  oss << k;
	}
	oss << ((char)27) << "[m";
	oss << " ";
	if (k==51 || k==119 || k==187 || k==255) {
	  oss << endl;
	}
      }
    }
  }
  return oss.str();
}

std::string terminal::three_byte_summ() {
  if (redirected) return "";
  std::ostringstream oss;
  for(size_t i=0;i<256;i+=51) {
    for(size_t j=0;j<256;j+=51) {
      for(size_t k=0;k<256;k+=51) {
	oss << ((char)27) << "[38;2;" << i << ";" << j << ";" << k << "m";
	if (i<10) {
	  oss << "  " << i;
	} else if (i<100) {
	  oss << " " << i;
	} else {
	  oss << i;
	}
	if (j<10) {
	  oss << "  " << j;
	} else if (j<100) {
	  oss << " " << j;
	} else {
	  oss << j;
	}
	if (k<10) {
	  oss << "  " << k;
	} else if (k<100) {
	  oss << " " << k;
	} else {
	  oss << k;
	}
	oss << ((char)27) << "[m";
	oss << " ";
	oss << ((char)27) << "[48;2;" << i << ";" << j << ";" << k << "m";
	if (i<10) {
	  oss << "  " << i;
	} else if (i<100) {
	  oss << " " << i;
	} else {
	  oss << i;
	}
	if (j<10) {
	  oss << "  " << j;
	} else if (j<100) {
	  oss << " " << j;
	} else {
	  oss << j;
	}
	if (k<10) {
	  oss << "  " << k;
	} else if (k<100) {
	  oss << " " << k;
	} else {
	  oss << k;
	}
	oss << ((char)27) << "[m";
	oss << " ";
	if (k==102 || k==255) {
	  oss << endl;
	}
      }
    }
  }
  return oss.str();
}

std::string terminal::eight_bit_summ() {
  if (redirected) return "";
  std::ostringstream oss;
  for(size_t i=0;i<256;i++) {
    oss << ((char)27) << "[38;5;" << i << "m";
    if (i<10) {
      oss << "  " << i;
    } else if (i<100) {
      oss << " " << i;
    } else {
      oss << i;
    }
    oss << ((char)27) << "[m";
    oss << " ";
    oss << ((char)27) << "[48;5;" << i << "m";
    if (i<10) {
      oss << "  " << i;
    } else if (i<100) {
      oss << " " << i;
    } else {
      oss << i;
    }
    oss << ((char)27) << "[m";
    oss << " ";
    if (i%10==9) {
      oss << endl;
    }
  }
  return oss.str();
}

size_t terminal::str_len(std::string str) {
  size_t cnt=0, len=str.length();
  for(size_t i=0;i<len;i++) {
    int ic=((int)str[i]);
    if (ic!=27) {
      cnt++;
    } else if (i+2<len && str[i+1]=='[' && str[i+2]=='m') {
      i+=2;
    } else if (i+2<len && str[i+1]=='[' && str[i+3]=='m') {
      i+=3;
    } else if (i+2<len && str[i+1]=='[' && str[i+4]=='m') {
      i+=4;
    } else if (i+2<len && str[i+1]=='(' && str[i+2]=='0') {
      i+=2;
    } else if (i+2<len && str[i+1]=='(' && str[i+2]=='B') {
      i+=2;
    }
  }
  return cnt;
}
