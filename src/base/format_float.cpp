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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "format_float.h"

// For gsl_finite()
#include <gsl/gsl_sys.h>

using namespace std;
using namespace o2scl;

format_float::format_float() {
  prefx="";
  sgn="-";
  suffx="";
  sci_prefx="";
  tmes=" x ";
  exp_prefx="^{";
  exp_sgn="-";
  sci_sgn="-";
  exp_suffx="}";
  sci_suffx="";
  not_finte="Nan";
  zeros="0";
  ex_mn=-2;
  ex_mx=3;
  sig_fgs=5;
  pad_zeros=false;
  dpt='.';
  exp_dgs=0;
  show_exp_sgn=false;
}

void format_float::html_mode() {
  prefx="";
  sgn="-";
  suffx="";
  sci_prefx="";
  tmes=" &times; ";
  exp_prefx="10<sup>";
  exp_sgn="-";
  sci_sgn="-";
  exp_suffx="</sup>";
  sci_suffx="";
  not_finte="Nan";
  zeros="0";
  exp_dgs=0;
  show_exp_sgn=false;
  return;
}

void format_float::latex_mode() {
  prefx="";
  sgn="$-$";
  suffx="";
  sci_prefx="";
  tmes=" $\\times ";
  exp_prefx="10^{";
  sci_sgn="$-$";
  exp_sgn="-";
  exp_suffx="}$";
  sci_suffx="";
  not_finte="Nan";
  zeros="0";
  exp_dgs=0;
  show_exp_sgn=false;
  return;
}      

void format_float::c_mode() {
  prefx="";
  sgn="-";
  suffx="";
  sci_prefx="";
  tmes="";
  exp_prefx="e";
  exp_sgn="-";
  sci_sgn="-";
  exp_suffx="";
  sci_suffx="";
  not_finte="NaN";
  zeros="0";
  ex_mn=-4;
  ex_mx=5;
  exp_dgs=2;
  show_exp_sgn=true;
  return;
}

int format_float::remove_zeros_dpt(string &s) {

  // If necessary, remove extra zeros
  if (!pad_zeros) {
    
    size_t sz=s.length();
    for(int i=((int)(sz-1));i>0;i--) {
      if (s[i]=='0') {
	sz=i;
      } else {
	i=0;
      }
    }
    s=s.substr(0,sz);

    // Remove extra decimal point if necessary
    string dptemp=s.substr(s.length()-dpt.size(),dpt.size());
    if (dptemp==dpt) {
      s=s.substr(0,s.length()-dpt.size());
    }

  }
  
  return 0;
}

string format_float::convert(double x, bool debug) {

  // Handle special cases
  if (!o2scl::is_finite(x)) return not_finte;
  if (x==0.0) return zeros;

  if (debug) cout.setf(ios::scientific);
  if (debug) cout << "Sig_figs: " << sig_fgs << endl;

  // Check sign
  int sign=1;
  if (x<0.0) {
    x=-x; 
    sign=-1;
  }

  // Compute exponent and mantissa separately
  int expo=((int)log10(x)), itmp;
  double mant=x/pow(10.0,expo);

  // This prevented errors with the mantissa=10 case on isospin,
  // possibly because it prevents too much optimizing the calculation
  // above with the conditionals below.
  ostringstream s2;
  s2 << mant;

  // Occasionally, finite precision errors compute the 
  // mantissa incorrectly. Fix this here.
  if (mant<1.0) {
    mant*=10.0;
    expo-=1;
  }
  if (mant>=10.0) {
    mant/=10.0;
    expo+=1;
  }

  if (debug) {
    cout << "Expo: " << expo << " Mant: " << mant << " "
	 << "Mant-10: " << mant-10.0 << endl;
  }
  
  // Pick off the digits
  int digits[16];
  string sx;
  
  if (sig_fgs>1) {

    string digits_str=dtos(mant,sig_fgs-1);

    sx+=digits_str[0];
    digits[0]=o2scl::stoi(sx);
    sx.clear();
    if (debug) {
      cout << "Digits_str: " << digits_str << endl;
      cout << "Digits: " << endl;
      cout << digits[0] << " ";
    }
    for(size_t i=1;i<sig_fgs;i++) {
      sx+=digits_str[i+1];
      digits[i]=o2scl::stoi(sx);
      sx.clear();
      if (debug) cout << digits[i] << " ";
    }
    if (debug) cout << endl;

  } else {
    
    // dtos() fails to round when sig_fgs is 1 (i.e. 1.6 doesn't easily
    // round to 2), so we manually do it here
    
    // Ensure truncation works by adding a small amount to the mantissa
    double eps=5.0*pow(10.0,-((int)sig_fgs));
    mant+=eps;
    
    // We want the string conversion to give us at least one extra digit
    string digits_str=dtos(mant,17);
    sx+=digits_str[0];
    digits[0]=o2scl::stoi(sx);
    sx.clear();
    if (debug) {
      cout << "Digits_str: " << digits_str << endl;
      cout << "Digits: " << endl;
      cout << digits[0] << " ";
    }
    // This number 16 must be one less than the size of the 
    // digits array
    for(size_t i=1;i<15;i++) {
      sx+=digits_str[i+1];
      digits[i]=o2scl::stoi(sx);
      sx.clear();
      if (debug) cout << digits[i] << " ";
    }
    if (debug) cout << endl;

  }
  
  // Begin constructing the string for output
  string s;

  if (expo<=ex_mx && expo>=ex_mn) {

    // Normal notation

    // Prefix and sign
    if (sign==-1) {
      s=prefx+sgn;
    } else {
      s=prefx;
    }

    // Digits and decimal point
    if (sig_fgs==1) {

      if (expo>=0) {
	// If exponent is positive or zero
	s+=itos(digits[0]);
	for(int i=0;i<expo;i++) {
	  s+="0";
	}
      } else {
	// If exponent is negative
	s="0"+dpt;
	for(int i=-1;i>expo;i--) {
	  s+="0";
	}
	s+=itos(digits[0]);
      }

    } else {

      if (expo>=0) {
	// If exponent is positive or zero
	int i=0;
	// Digits before decimal point
	s+=itos(digits[0]);
	for(i=0;i<expo;i++) {
	  if (i<((int)(sig_fgs))-1) {
	    s+=itos(digits[i+1]);
	  } else {
	    s+="0";
	  }
	}
	// Digits after decimal point
	if (i<((int)sig_fgs)-1) {
	  s+=dpt;
	  for(;i<((int)sig_fgs)-1;i++) {
	    s+=itos(digits[i+1]);
	  }
	  remove_zeros_dpt(s);
	}
      } else {
	// If exponent is negative
	s+="0"+dpt;
	for(int i=-1;i>expo;i--) {
	  s+="0";
	}
	for(size_t i=0;i<sig_fgs;i++) {
	  s+=itos(digits[i]);
	}
	remove_zeros_dpt(s);
      }

    }
    
    s+=suffx;
	
  } else {

    // Scientific notation

    // Compute sign of exponent
    int expo_sign=1;
    if (expo<0) {
      expo_sign=-1;
      expo*=-1;
    }
	
    // Prefix and sign
    if (sign==-1) {
      s=sci_prefx+sci_sgn;
    } else {
      s=sci_prefx;
    }

    // Digits and decimal point
    if (sig_fgs==1) {
      s+=itos(digits[0]);
    } else {
      s+=itos(digits[0])+dpt;
      for(size_t i=1;i<sig_fgs;i++) s+=itos(digits[i]);
    }

    remove_zeros_dpt(s);

    // Add zeros to exponent if necessary
    string expo_str=itos(expo);
    if (exp_dgs>0 && expo_str.length()<exp_dgs) {
      for(size_t i=0;i<exp_dgs-expo_str.length();i++) {
	expo_str='0'+expo_str;
      }
    }
    
    // Add exponent and suffix
    if (expo_sign<0) {
      s+=tmes+exp_prefx+exp_sgn+expo_str+exp_suffx+sci_suffx;
    } else {
      if (show_exp_sgn) {
	s+=tmes+exp_prefx+"+"+expo_str+exp_suffx+sci_suffx;
      } else {
	s+=tmes+exp_prefx+expo_str+exp_suffx+sci_suffx;
      }
    }
	
  }

  if (debug) {
    cout << "Final: " << s << endl;
    //char ch;
    //cin >> ch;
  }
      
  return s;

}
