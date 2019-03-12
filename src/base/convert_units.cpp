/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <o2scl/convert_units.h>

using namespace std;
using namespace o2scl;

convert_units::convert_units() {
  verbose=0;
  use_gnu_units=true;
  units_cmd_string="units";
  err_on_fail=true;
  combine_two_conv=true;
}

double convert_units::convert(std::string from, std::string to,
			      double val)  {
  double converted;
  int ret=convert_ret(from,to,val,converted);
  if (ret==exc_efilenotfound) {
    O2SCL_ERR("Pipe could not be opened in convert_units::convert().",
	      exc_efilenotfound);
  }
  if (ret==exc_efailed) {
    O2SCL_ERR("Pipe could not be closed in convert_units::convert().",
	      exc_efailed);
  }
  if (ret==exc_enotfound) {
    string str=((string)"Conversion between ")+from+" and "+to+
      " not found in convert_units::convert().";
    O2SCL_ERR(str.c_str(),exc_enotfound);
  }
  return converted;
}

double convert_units::convert_const(std::string from, std::string to,
			      double val) const {
  double converted;
  int ret=convert_ret_const(from,to,val,converted);
  if (ret==exc_enotfound) {
    string str=((string)"Conversion between ")+from+" and "+to+
      " not found in convert_units::convert().";
    O2SCL_ERR(str.c_str(),exc_enotfound);
  }
  return converted;
}

int convert_units::convert_ret(std::string from, std::string to, double val,
			       double &converted) {

  double factor;
  bool new_conv;
  
  int ret=convert_internal(from,to,val,converted,factor,new_conv);

  if (ret==0 && new_conv) {
    // Add the newly computed conversion to the table
    unit_t ut;
    ut.f=from;
    ut.t=to;
    ut.c=factor;
    std::string both=from+","+to;
    mcache.insert(make_pair(both,ut));
  }

  return ret;
}

int convert_units::convert_cache(std::string from, std::string to,
				 double val, double &converted,
				 double &factor) const {
  
  // Look in cache for conversion
  std::string both=from+","+to;
  mciter m3=mcache.find(both);
  if (m3!=mcache.end()) {
    factor=m3->second.c;
    converted=val*factor;
    return 0;
  }

  // Look in cache for reverse conversion
  std::string both2=to+","+from;
  m3=mcache.find(both2);
  if (m3!=mcache.end()) {
    factor=1.0/m3->second.c;
    converted=val*factor;
    return 0;
  }

  if (combine_two_conv) {
    
    // Look for combined conversions
    for(mciter m=mcache.begin();m!=mcache.end();m++) {
      if (m->second.f==from) {
	std::string b=m->second.t+","+to;
	mciter m2=mcache.find(b);
	if (m2!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << " " << m->second.c << std::endl;
	    std::cout << " (1)          and: " << m2->second.f << " , "
		      << m2->second.t << " " << m2->second.c << std::endl;
	  }
	  factor=m->second.c*m2->second.c;
	  converted=val*factor;
	  return 0;
	}
	std::string b2=to+","+m->second.t;
	mciter m4=mcache.find(b2);
	if (m4!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << std::endl;
	    std::cout << " (2)          and: " << m4->second.f << " , "
		      << m4->second.t << std::endl;
	  }
	  factor=m->second.c/m4->second.c;
	  converted=val*factor;
	  return 0;
	}
      } else if (m->second.t==from) {
	std::string b=m->second.f+","+to;
	mciter m2=mcache.find(b);
	if (m2!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << std::endl;
	    std::cout << " (3)          and: " << m2->second.f << " , "
		      << m2->second.t << std::endl;
	  }
	  factor=m2->second.c/m->second.c;
	  converted=val*factor;
	  return 0;
	}
      } else if (m->second.f==to) {
	std::string b=m->second.t+","+from;
	mciter m2=mcache.find(b);
	if (m2!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << std::endl;
	    std::cout << " (4)          and: " << m2->second.f << " , "
		      << m2->second.t << std::endl;
	  }
	  factor=1.0/m->second.c/m2->second.c;
	  converted=val*factor;
	  return 0;
	}
	std::string b2=from+","+m->second.t;
	mciter m4=mcache.find(b2);
	if (m4!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << std::endl;
	    std::cout << " (5)          and: " << m4->second.f << " , "
		      << m4->second.t << std::endl;
	  }
	  factor=m4->second.c/m->second.c;
	  converted=val*factor;
	  return 0;
	}
      } else if (m->second.t==to) {
	std::string b=m->second.f+","+from;
	mciter m2=mcache.find(b);
	if (m2!=mcache.end()) {
	  if (verbose>0) {
	    std::cout << "Using conversions: " << m->second.f << " , "
		      << m->second.t << std::endl;
	    std::cout << " (6)          and: " << m2->second.f << " , "
		      << m2->second.t << std::endl;
	  }
	  factor=m->second.c/m2->second.c;
	  converted=val*factor;
	  return 0;
	}
      }
    }

  }

  return exc_efailed;
}

int convert_units::convert_gnu_units(std::string from, std::string to,
				     double val, double &converted,
				     double &factor, bool &new_conv) const {
  
  // Run the GNU 'units' command
  std::string cmd=units_cmd_string+" '"+from+"' '"+to+"'";
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Units command is " << cmd << std::endl;
  }
  
#ifdef HAVE_POPEN
  
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Define constant popen is defined." << std::endl;
  }
  
  cout << "Using GNU units: " << cmd << endl;
  FILE *ps_pipe=popen(cmd.c_str(),"r");
  if (!ps_pipe) {
    if (err_on_fail) {
      O2SCL_ERR2("Pipe could not be opened in ",
		 "convert_units::convert_gnu_units().",exc_efilenotfound);
    }
    return exc_efilenotfound;
  }
  
  char line1[80];
  int size=80;
  
  // Variable 'cret' is unused, but put here to avoid
  // unused return value errors
  char *cret=fgets(line1,80,ps_pipe);
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Units output is " << line1 << std::endl;
  }
  
  // Read the output from the 'units' command and compute the 
  // conversion factor
  std::string s=line1, t1, t2;
  std::istringstream *ins=new std::istringstream(s);
  (*ins) >> t1 >> t2;
  delete ins;
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Units string to convert is "
	      << t2 << std::endl;
  }
  int sret=o2scl::stod_nothrow(t2,factor);
  if (sret!=0) {
    if (err_on_fail) {
      string str=((string)"Conversion between ")+from+" and "+to+
	" not found in convert_units::convert_gnu_units().";
      O2SCL_ERR(str.c_str(),exc_enotfound);
    }
    return exc_enotfound;
  }
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Converted value is "
	      << factor*val << std::endl;
  }
  
  // Cleanup
  if (pclose(ps_pipe)!=0) {
    if (err_on_fail) {
      O2SCL_ERR2("Pipe could not be closed in ",
		 "convert_units::convert_gnu_units().",exc_efailed);
    }
    return exc_efailed;
  }

  converted=factor*val;
  
#else
  
  if (verbose>0) {
    std::cout << "convert_units::convert_gnu_units(): "
	      << "Define constant popen is not defined." << std::endl;
  }
  
  return exc_efailed;
  
#endif
      
  return 0;
}

int convert_units::convert_internal(std::string from, std::string to,
				    double val, double &converted,
				    double &factor, bool &new_conv) const {

  // Remove whitespace
  remove_whitespace(from);
  remove_whitespace(to);

  int ret_cache=convert_cache(from,to,val,converted,factor);

  if (ret_cache==0) {
    new_conv=false;
    return 0;
  }
  
  if (use_gnu_units) {

    if (verbose>0) {
      std::cout << "convert_units::convert(): "
		<< "Value of use_gnu_units is true." << std::endl;
    }

    int ret_gnu=convert_gnu_units(from,to,val,converted,factor,new_conv);

    if (ret_gnu==0) {
      new_conv=true;
      return 0;
    }

  } else {
    
    if (verbose>0) {
      std::cout << "convert_units::convert(): "
		<< "Value of use_gnu_units is false." << std::endl;
    }
    
  }    
  
  if (err_on_fail) {
    string str=((string)"Conversion between ")+from+" and "+to+
      " not found in convert_units::convert_ret().";
    O2SCL_ERR(str.c_str(),exc_enotfound);
  }
  
  return exc_enotfound;
}

int convert_units::test_cache() {
  err_on_fail=false;
  mciter m, m2;
  cout << "units_cmd_string: " << units_cmd_string << endl;
  for (m=mcache.begin();m!=mcache.end();m++) {
    for (m2=m;m2!=mcache.end();m2++) {
      string from=m->second.f;
      string to=m2->second.t;
      if (from!=to) {
	double v=1.0, c, f1, f2;
	int cret=convert_cache(from,to,v,c,f1);
	if (cret==0) {
	  bool new_conv;
	  int gret=convert_gnu_units(from,to,v,c,f2,new_conv);
	  if (gret==0) {
	    if (fabs(f1-f2)/f1>1.0e-6) {
	      cout << "* ";
	    } else {
	      cout << "  ";
	    }
	    cout.width(10);
	    cout << from << " ";
	    cout.width(10);
	    cout << to << " " << f1 << " " << f2 << " "
		 << fabs(f1-f2)/f1 << endl;
	  }
	}
      }
      to=m2->second.f;
      if (from!=to) {
	double v=1.0, c, f1, f2;
	int cret=convert_cache(from,to,v,c,f1);
	if (cret==0) {
	  bool new_conv;
	  int gret=convert_gnu_units(from,to,v,c,f2,new_conv);
	  if (gret==0) {
	    if (fabs(f1-f2)/f1>1.0e-6) {
	      cout << "* ";
	    } else {
	      cout << "  ";
	    }
	    cout.width(10);
	    cout << from << " ";
	    cout.width(10);
	    cout << to << " " << f1 << " " << f2 << " "
		 << fabs(f1-f2)/f1 << endl;
	  }
	}
      }
    }
  }
  return 0;
}

int convert_units::convert_ret_const(std::string from, std::string to,
				     double val, double &converted) const {
				    
  double factor;
  bool new_conv;
  
  return convert_internal(from,to,val,converted,factor,new_conv);
}

void convert_units::remove_cache(std::string from, std::string to) {

  // Remove whitespace
  remove_whitespace(from);
  remove_whitespace(to);
  std::string both=from+","+to;

  miter m3=mcache.find(both);
  if (m3!=mcache.end()) {
    mcache.erase(m3);
    return;
  }
  
  if (err_on_fail) {
    O2SCL_ERR((((string)"Conversion ")+from+" -> "+to+
	       " not found in convert_units::remove_cache().").c_str(),
	      exc_enotfound);
  }

  return;
}

void convert_units::insert_cache
(std::string from, std::string to, double conv) {

  // Remove whitespace
  remove_whitespace(from);
  remove_whitespace(to);
  
  if (err_on_fail && 
      (from.find(',')!=std::string::npos || 
       to.find(',')!=std::string::npos)) {
    O2SCL_ERR("Units cannot contain comma in insert_cache()",
	      exc_efailed);
  }

  unit_t ut;
  ut.f=from;
  ut.t=to;
  ut.c=conv;
  std::string both=from+","+to;

  // If it's already inside the cache, just update the conversion
  // value
  miter m3=mcache.find(both);
  if (m3!=mcache.end()) {
    m3->second=ut;
  }

  mcache.insert(make_pair(both,ut));
  return;
}
    
void convert_units::print_cache() const {
  mciter m;
  if (mcache.size()==0) {
    cout << "No units in cache." << endl;
  } else {
    std::cout << "Unit cache: " << std::endl;
    for (m=mcache.begin();m!=mcache.end();m++) {
      std::cout.setf(std::ios::left);
      std::cout.width(25);
      std::cout << m->second.f << " ";
      std::cout.width(25);
      std::cout << m->second.t << " ";
      std::cout.unsetf(std::ios::left);
      std::cout.precision(10);
      std::cout << m->second.c << std::endl;
    }
  }
  return;
}

void convert_units::make_units_dat(std::string fname, bool c_1, 
				   bool hbar_1, bool K_1) const {
  
  std::ofstream fout(fname.c_str());
  fout.precision(14);

  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# Fundamental units" << std::endl;
  fout << "m\t!" << std::endl;
  fout << "meter\tm" << std::endl;
  if (c_1==false) {
    fout << "s\t!" << std::endl;
    fout << "second\ts" << std::endl;
  } else {
    fout << "s\t" << o2scl_mks::speed_of_light << " m" << std::endl;
    fout << "second\ts" << std::endl;
  }
  if (hbar_1==false) {
    fout << "kg\t!" << std::endl;
    fout << "kilogram\tkg" << std::endl;
  } else {
    fout << "kg\t" << 1.0/o2scl_mks::plancks_constant_hbar 
	 << " s / m^2" << std::endl;
    fout << "kilogram\tkg" << std::endl;
  }
  fout << "A\t!" << std::endl;
  fout << "ampere\tA" << std::endl;
  fout << "amp\tA" << std::endl;
  fout << "cd\t!" << std::endl;
  fout << "candela\tcd" << std::endl;
  fout << "mol\t!" << std::endl;
  fout << "mole\tmol" << std::endl;
  if (K_1==false) {
    fout << "K\t!" << std::endl;
    fout << "kelvin\tK" << std::endl;
  } else {
    fout << "K\t" << o2scl_mks::boltzmann << " kg m^2 / s^2" << std::endl;
    fout << "kelvin\tK" << std::endl;
  }
  fout << "radian\t!" << std::endl;
  fout << "sr\t!" << std::endl;
  fout << "steradian\tsr" << std::endl;
  fout << "US$\t!" << std::endl;
  fout << "bit\t!" << std::endl;
  fout << std::endl;
      
  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# SI and common prefixes" << std::endl;
  fout << "yotta-\t\t1e24" << std::endl;
  fout << "zetta-\t\t1e21" << std::endl;
  fout << "exa-\t\t1e18" << std::endl;
  fout << "peta-\t\t1e15" << std::endl;
  fout << "tera-\t\t1e12" << std::endl;
  fout << "giga-\t\t1e9" << std::endl;
  fout << "mega-\t\t1e6" << std::endl;
  fout << "myria-\t\t1e4" << std::endl;
  fout << "kilo-\t\t1e3" << std::endl;
  fout << "hecto-\t\t1e2" << std::endl;
  fout << "deca-\t\t1e1" << std::endl;
  fout << "deka-\t\tdeca" << std::endl;
  fout << "deci-\t\t1e-1" << std::endl;
  fout << "centi-\t\t1e-2" << std::endl;
  fout << "milli-\t\t1e-3" << std::endl;
  fout << "micro-\t\t1e-6" << std::endl;
  fout << "nano-\t\t1e-9" << std::endl;
  fout << "pico-\t\t1e-12" << std::endl;
  fout << "femto-\t\t1e-15" << std::endl;
  fout << "atto-\t\t1e-18" << std::endl;
  fout << "zepto-\t\t1e-21" << std::endl;
  fout << "yocto-\t\t1e-24" << std::endl;
  fout << "quarter-\t1|4" << std::endl;
  fout << "semi-\t\t0.5" << std::endl;
  fout << "demi-\t\t0.5" << std::endl;
  fout << "hemi-\t\t0.5" << std::endl;
  fout << "half-\t\t0.5" << std::endl;
  fout << "double-\t\t2" << std::endl;
  fout << "triple-\t\t3" << std::endl;
  fout << "treble-\t\t3" << std::endl;
  fout << std::endl;

  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# SI prefix abbreviations" << std::endl;
  fout << "Y-                      yotta" << std::endl;
  fout << "Z-                      zetta" << std::endl;
  fout << "E-                      exa" << std::endl;
  fout << "P-                      peta" << std::endl;
  fout << "T-                      tera" << std::endl;
  fout << "G-                      giga" << std::endl;
  fout << "M-                      mega" << std::endl;
  fout << "k-                      kilo" << std::endl;
  fout << "h-                      hecto" << std::endl;
  fout << "da-                     deka" << std::endl;
  fout << "d-                      deci" << std::endl;
  fout << "c-                      centi" << std::endl;
  fout << "m-                      milli" << std::endl;
  fout << "n-                      nano" << std::endl;
  fout << "p-                      pico" << std::endl;
  fout << "f-                      femto" << std::endl;
  fout << "a-                      atto" << std::endl;
  fout << "z-                      zepto" << std::endl;
  fout << "y-                      yocto" << std::endl;
  fout << std::endl;

  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# Basic numbers" << std::endl;
  fout << "one                     1" << std::endl;
  fout << "two                     2" << std::endl;
  fout << "double                  2" << std::endl;
  fout << "couple                  2" << std::endl;
  fout << "three                   3" << std::endl;
  fout << "triple                  3" << std::endl;
  fout << "four                    4" << std::endl;
  fout << "quadruple               4" << std::endl;
  fout << "five                    5" << std::endl;
  fout << "quintuple               5" << std::endl;
  fout << "six                     6" << std::endl;
  fout << "seven                   7" << std::endl;
  fout << "eight                   8" << std::endl;
  fout << "nine                    9" << std::endl;
  fout << "ten                     10" << std::endl;
  fout << "twenty                  20" << std::endl;
  fout << "thirty                  30" << std::endl;
  fout << "forty                   40" << std::endl;
  fout << "fifty                   50" << std::endl;
  fout << "sixty                   60" << std::endl;
  fout << "seventy                 70" << std::endl;
  fout << "eighty                  80" << std::endl;
  fout << "ninety                  90" << std::endl;
  fout << "hundred                 100" << std::endl;
  fout << "thousand                1000" << std::endl;
  fout << "million                 1e6" << std::endl;
  fout << "billion                 1e9" << std::endl;
  fout << "trillion                1e12" << std::endl;
  fout << "quadrillion             1e15" << std::endl;
  fout << "quintillion             1e18" << std::endl;
  fout << "sextillion              1e21" << std::endl;
  fout << "septillion              1e24" << std::endl;
  fout << "octillion               1e27" << std::endl;
  fout << "nonillion               1e30" << std::endl;
  fout << "noventillion            nonillion" << std::endl;
  fout << "decillion               1e33" << std::endl;
  fout << "undecillion             1e36" << std::endl;
  fout << "duodecillion            1e39" << std::endl;
  fout << "tredecillion            1e42" << std::endl;
  fout << "quattuordecillion       1e45" << std::endl;
  fout << "quindecillion           1e48" << std::endl;
  fout << "sexdecillion            1e51" << std::endl;
  fout << "septendecillion         1e54" << std::endl;
  fout << "octodecillion           1e57" << std::endl;
  fout << "novemdecillion          1e60" << std::endl;
  fout << "vigintillion            1e63" << std::endl;
  fout << "googol                  1e100" << std::endl;
  fout << "centillion              1e303" << std::endl;
  fout << std::endl;

  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# Basic SI units" << std::endl;
  fout << "newton                  kg m / s^2   " << std::endl;
  fout << "N                       newton" << std::endl;
  fout << "pascal                  N/m^2        " << std::endl;
  fout << "Pa                      pascal" << std::endl;
  fout << "joule                   N m          " << std::endl;
  fout << "J                       joule" << std::endl;
  fout << "watt                    J/s          " << std::endl;
  fout << "W                       watt" << std::endl;
  fout << "coulomb                 A s          " << std::endl;
  fout << "C                       coulomb" << std::endl;
  fout << "volt                    W/A          " << std::endl;
  fout << "V                       volt" << std::endl;
  fout << "ohm                     V/A          " << std::endl;
  fout << "siemens                 A/V          " << std::endl;
  fout << "S                       siemens" << std::endl;
  fout << "farad                   C/V          " << std::endl;
  fout << "F                       farad" << std::endl;
  fout << "weber                   V s          " << std::endl;
  fout << "Wb                      weber" << std::endl;
  fout << "henry                   Wb/A         " << std::endl;
  fout << "H                       henry" << std::endl;
  fout << "tesla                   Wb/m^2       " << std::endl;
  fout << "T                       tesla" << std::endl;
  fout << "hertz                   1/s           " << std::endl;
  fout << "Hz                      hertz" << std::endl;
  fout << "gram                    millikg       " << std::endl;
  fout << "g                       gram" << std::endl;
  fout << std::endl;

  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# Dimensional analysis units" << std::endl;
  fout << "LENGTH                  meter" << std::endl;
  fout << "AREA                    LENGTH^2" << std::endl;
  fout << "VOLUME                  LENGTH^3" << std::endl;
  fout << "MASS                    kilogram" << std::endl;
  fout << "CURRENT                 ampere" << std::endl;
  fout << "AMOUNT                  mole" << std::endl;
  fout << "ANGLE                   radian" << std::endl;
  fout << "SOLID_ANGLE             steradian" << std::endl;
  fout << "MONEY                   US$" << std::endl;
  fout << "FORCE                   newton" << std::endl;
  fout << "PRESSURE                FORCE / AREA" << std::endl;
  fout << "STRESS                  FORCE / AREA" << std::endl;
  fout << "CHARGE                  coulomb" << std::endl;
  fout << "CAPACITANCE             farad" << std::endl;
  fout << "RESISTANCE              ohm" << std::endl;
  fout << "CONDUCTANCE             siemens" << std::endl;
  fout << "INDUCTANCE              henry" << std::endl;
  fout << "FREQUENCY               hertz" << std::endl;
  fout << "VELOCITY                LENGTH / TIME" << std::endl;
  fout << "ACCELERATION            VELOCITY / TIME" << std::endl;
  fout << "DENSITY                 MASS / VOLUME" << std::endl;
  fout << "LINEAR_DENSITY          MASS / LENGTH" << std::endl;
  fout << "VISCOSITY               FORCE TIME / AREA" << std::endl;
  fout << "KINEMATIC_VISCOSITY     VISCOSITY / DENSITY" << std::endl;
  fout << std::endl;
      
  fout << "################################################" 
       << "##############" << std::endl;
  fout << "# GSL constants" << std::endl;
  fout << "schwarzchild_radius\t\t" << o2scl_mks::schwarzchild_radius
       << " m" << std::endl;
  fout << "Rschwarz\t\tschwarzchild_radius" << std::endl;
  fout << "speed_of_light\t\t" << o2scl_mks::speed_of_light
       << " m / s" << std::endl;
  fout << "c\t\tspeed_of_light" << std::endl;
  fout << "gravitational_constant\t\t" 
       << o2scl_mks::gravitational_constant
       << " m^3 / kg s^2" << std::endl;
  fout << "G\t\tgravitational_constant" << std::endl;
  fout << "plancks_constant_h\t\t" 
       << o2scl_mks::plancks_constant_h
       << " kg m^2 / s" << std::endl;
  fout << "h\t\tplancks_constant_h" << std::endl;
  fout << "plancks_constant_hbar\t\t" 
       << o2scl_mks::plancks_constant_hbar
       << " kg m^2 / s" << std::endl;
  fout << "hbar\t\tplancks_constant_hbar" << std::endl;
  fout << "astronomical_unit\t\t" << o2scl_mks::astronomical_unit
       << " m" << std::endl;
  fout << "au\t\tastronomical_unit" << std::endl;
  fout << "light_year\t\t" << o2scl_mks::light_year
       << " m" << std::endl;
  fout << "lyr\t\tlight_year" << std::endl;
  fout << "parsec\t\t" << o2scl_mks::parsec
       << " m" << std::endl;
  fout << "pc\t\tparsec" << std::endl;
  fout << "grav_accel\t\t" << o2scl_mks::grav_accel
       << " m / s^2" << std::endl;
  fout << "electron_volt\t\t" << o2scl_mks::electron_volt
       << " kg m^2 / s^2" << std::endl;
  fout << "eV\t\telectron_volt" << std::endl;
  fout << "mass_electron\t\t" << o2scl_mks::mass_electron
       << " kg" << std::endl;
  fout << "mass_muon\t\t" << o2scl_mks::mass_muon
       << " kg" << std::endl;
  fout << "mass_proton\t\t" << o2scl_mks::mass_proton
       << " kg" << std::endl;
  fout << "mass_neutron\t\t" << o2scl_mks::mass_neutron
       << " kg" << std::endl;
  fout << "rydberg\t\t" << o2scl_mks::rydberg
       << " kg m^2 / s^2" << std::endl;
  fout << "boltzmann\t\t" << o2scl_mks::boltzmann
       << " kg m^2 / K s^2" << std::endl;
  fout << "bohr_magneton\t\t" << o2scl_mks::bohr_magneton
       << " A m^2" << std::endl;
  fout << "nuclear_magneton\t\t" << o2scl_mks::nuclear_magneton
       << " A m^2" << std::endl;
  fout << "electron_magnetic_moment\t\t" 
       << o2scl_mks::electron_magnetic_moment
       << " A m^2" << std::endl;
  fout << "proton_magnetic_moment\t\t" 
       << o2scl_mks::proton_magnetic_moment
       << " A m^2" << std::endl;
  fout << "molar_gas\t\t" << o2scl_mks::molar_gas
       << " kg m^2 / K mol s^2" << std::endl;
  fout << "standard_gas_volume\t\t" << o2scl_mks::standard_gas_volume
       << " m^3 / mol" << std::endl;
  fout << "minute\t\t" << o2scl_mks::minute
       << " s" << std::endl;
  fout << "min\t\tminute" << std::endl;
  fout << "hour\t\t" << o2scl_mks::hour
       << " s" << std::endl;
  fout << "day\t\t" << o2scl_mks::day
       << " s" << std::endl;
  fout << "week\t\t" << o2scl_mks::week
       << " s" << std::endl;
  fout << "inch\t\t" << o2scl_mks::inch
       << " m" << std::endl;
  fout << "foot\t\t" << o2scl_mks::foot
       << " m" << std::endl;
  fout << "yard\t\t" << o2scl_mks::yard
       << " m" << std::endl;
  fout << "mile\t\t" << o2scl_mks::mile
       << " m" << std::endl;
  fout << "nautical_mile\t\t" << o2scl_mks::nautical_mile
       << " m" << std::endl;
  fout << "fathom\t\t" << o2scl_mks::fathom
       << " m" << std::endl;
  fout << "mil\t\t" << o2scl_mks::mil
       << " m" << std::endl;
  fout << "point\t\t" << o2scl_mks::point
       << " m" << std::endl;
  fout << "texpoint\t\t" << o2scl_mks::texpoint
       << " m" << std::endl;
  fout << "micron\t\t" << o2scl_mks::micron
       << " m" << std::endl;
  fout << "angstrom\t\t" << o2scl_mks::angstrom
       << " m" << std::endl;
  fout << "hectare\t\t" << o2scl_mks::hectare
       << " m^2" << std::endl;
  fout << "acre\t\t" << o2scl_mks::acre
       << " m^2" << std::endl;
  fout << "barn\t\t" << o2scl_mks::barn
       << " m^2" << std::endl;
  fout << "liter\t\t" << o2scl_mks::liter
       << " m^3" << std::endl;
  fout << "us_gallon\t\t" << o2scl_mks::us_gallon
       << " m^3" << std::endl;
  fout << "gallon\t\tus_gallon" << std::endl;
  fout << "quart\t\t" << o2scl_mks::quart
       << " m^3" << std::endl;
  fout << "pint\t\t" << o2scl_mks::pint
       << " m^3" << std::endl;
  fout << "cup\t\t" << o2scl_mks::cup
       << " m^3" << std::endl;
  fout << "fluid_ounce\t\t" << o2scl_mks::fluid_ounce
       << " m^3" << std::endl;
  fout << "tablespoon\t\t" << o2scl_mks::tablespoon
       << " m^3" << std::endl;
  fout << "teaspoon\t\t" << o2scl_mks::teaspoon
       << " m^3" << std::endl;
  fout << "canadian_gallon\t\t" << o2scl_mks::canadian_gallon
       << " m^3" << std::endl;
  fout << "uk_gallon\t\t" << o2scl_mks::uk_gallon
       << " m^3" << std::endl;
  fout << "miles_per_hour\t\t" << o2scl_mks::miles_per_hour
       << " m / s" << std::endl;
  fout << "mph\t\tmiles_per_hour" << std::endl;
  fout << "kilometers_per_hour\t\t" << o2scl_mks::kilometers_per_hour
       << " m / s" << std::endl;
  fout << "kph\t\tkilometers_per_hour" << std::endl;
  fout << "knot\t\t" << o2scl_mks::knot
       << " m / s" << std::endl;
  fout << "pound_mass\t\t" << o2scl_mks::pound_mass
       << " kg" << std::endl;
  fout << "ounce_mass\t\t" << o2scl_mks::ounce_mass
       << " kg" << std::endl;
  fout << "ton\t\t" << o2scl_mks::ton
       << " kg" << std::endl;
  fout << "metric_ton\t\t" << o2scl_mks::metric_ton
       << " kg" << std::endl;
  fout << "uk_ton\t\t" << o2scl_mks::uk_ton
       << " kg" << std::endl;
  fout << "troy_ounce\t\t" << o2scl_mks::troy_ounce
       << " kg" << std::endl;
  fout << "carat\t\t" << o2scl_mks::carat
       << " kg" << std::endl;
  fout << "unified_atomic_mass\t\t" << o2scl_mks::unified_atomic_mass
       << " kg" << std::endl;
  fout << "gram_force\t\t" << o2scl_mks::gram_force
       << " kg m / s^2" << std::endl;
  fout << "pound_force\t\t" << o2scl_mks::pound_force
       << " kg m / s^2" << std::endl;
  fout << "kilopound_force\t\t" << o2scl_mks::kilopound_force
       << " kg m / s^2" << std::endl;
  fout << "poundal\t\t" << o2scl_mks::poundal
       << " kg m / s^2" << std::endl;
  fout << "calorie\t\t" << o2scl_mks::calorie
       << " kg m^2 / s^2" << std::endl;
  fout << "btu\t\t" << o2scl_mks::btu
       << " kg m^2 / s^2" << std::endl;
  fout << "therm\t\t" << o2scl_mks::therm
       << " kg m^2 / s^2" << std::endl;
  fout << "horsepower\t\t" << o2scl_mks::horsepower
       << " kg m^2 / s^3" << std::endl;
  fout << "hp\t\thorsepower" << std::endl;
  fout << "bar\t\t" << o2scl_mks::bar
       << " kg / m s^2" << std::endl;
  fout << "std_atmosphere\t\t" << o2scl_mks::std_atmosphere
       << " kg / m s^2" << std::endl;
  fout << "torr\t\t" << o2scl_mks::torr
       << " kg / m s^2" << std::endl;
  fout << "meter_of_mercury\t\t" << o2scl_mks::meter_of_mercury
       << " kg / m s^2" << std::endl;
  fout << "inch_of_mercury\t\t" << o2scl_mks::inch_of_mercury
       << " kg / m s^2" << std::endl;
  fout << "inch_of_water\t\t" << o2scl_mks::inch_of_water
       << " kg / m s^2" << std::endl;
  fout << "psi\t\t" << o2scl_mks::psi
       << " kg / m s^2" << std::endl;
  fout << "poise\t\t" << o2scl_mks::poise
       << " kg / m s " << std::endl;
  fout << "stokes\t\t" << o2scl_mks::stokes
       << " m^2 / s" << std::endl;
  fout << "faraday\t\t" << o2scl_mks::faraday
       << " A s / mol" << std::endl;
  fout << "electron_charge\t\t" << o2scl_mks::electron_charge
       << " A s" << std::endl;
  fout << "gauss\t\t" << o2scl_mks::gauss
       << " kg / A s^2" << std::endl;
  fout << "stilb\t\t" << o2scl_mks::stilb
       << " cd / m^2" << std::endl;
  fout << "lumen\t\t" << o2scl_mks::lumen
       << " cd sr" << std::endl;
  fout << "lux\t\t" << o2scl_mks::lux
       << " cd sr / m^2" << std::endl;
  fout << "phot\t\t" << o2scl_mks::phot
       << " cd sr / m^2" << std::endl;
  fout << "footcandle\t\t" << o2scl_mks::footcandle
       << " cd sr / m^2" << std::endl;
  fout << "lambert\t\t" << o2scl_mks::lambert
       << " cd sr / m^2" << std::endl;
  fout << "footlambert\t\t" << o2scl_mks::footlambert
       << " cd sr / m^2" << std::endl;
  fout << "curie\t\t" << o2scl_mks::curie
       << " 1 / s" << std::endl;
  fout << "roentgen\t\t" << o2scl_mks::roentgen
       << " A s / kg" << std::endl;
  fout << "rad\t\t" << o2scl_mks::rad
       << " m^2 / s^2" << std::endl;
  fout << "solar_mass\t\t" << o2scl_mks::solar_mass
       << " kg" << std::endl;
  fout << "Msun\t\tsolar_mass" << std::endl;
  fout << "bohr_radius\t\t" << o2scl_mks::bohr_radius
       << " m" << std::endl;
  fout << "dyne\t\t" << o2scl_mks::dyne
       << " kg m / s^2" << std::endl;
  fout << "erg\t\t" << o2scl_mks::erg
       << " kg m^2 / s^2" << std::endl;
  fout << "stefan_boltzmann_constant\t\t" 
       << o2scl_mks::stefan_boltzmann_constant
       << " kg / K^4 s^3" << std::endl;
  fout << "thomson_cross_section\t\t" 
       << o2scl_mks::thomson_cross_section
       << " m^2" << std::endl;
  fout << "vacuum_permittivity\t\t" << o2scl_mks::vacuum_permittivity
       << " A^2 s^4 / kg m^3" << std::endl;
  fout << "vacuum_permeability\t\t" << o2scl_mks::vacuum_permeability
       << " kg m / A^2 s^2" << std::endl;
  fout.close();
      
  return;
}

