/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2020, Andrew W. Steiner
  
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
#ifndef O2SCL_CONVERT_UNITS_H
#define O2SCL_CONVERT_UNITS_H

/** \file convert_units.h
    \brief File defining \ref o2scl::convert_units
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>

#include <o2scl/misc.h>
#include <o2scl/constants.h>
#include <o2scl/string_conv.h>
#include <o2scl/shunting_yard.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Convert units

      Allow the user to convert between two different units after
      specifying a conversion factor. This class will also
      automatically combine two conversion factors to create a new
      unit conversion (but it cannot combine more than two).

      Conversions are performed by the \ref convert() function. The
      run-time unit cache is initially filled with hard-coded
      conversions, and \ref convert() searches this cache is searched
      for the requested conversion first. If the conversion is not
      found and if \ref use_gnu_units is true, then \ref convert()
      tries to open a pipe to open a shell to get the conversion
      factor from <tt>'units'</tt>. If this is successful, then the
      conversion factor is returned and the conversion is added to the
      cache.

      If the GNU units command is not in the local path, the user may
      modify \ref units_cmd_string to specify the full pathname. One
      can also modify \ref units_cmd_string to specify a different
      <tt>units.dat</tt> file.

      Example:
      \code
      convert_units cu;
      cout << "A solar mass is " << cu.convert("Msun","g",1.0) 
      << " g. " << endl;
      \endcode

      An object of this type is created by \ref o2scl_settings
      (of type \ref lib_settings_class) for several unit 
      conversions used internally in \o2 .

      \note Combining two conversions allows for some surprising
      apparent contradictions from numerical precision errors. If
      there are two matching unit conversion pairs which give the same
      requested conversion factor, then one can arrange a situation
      where the same conversion factor is reported with slightly
      different values after adding a related conversion to the table.
      One way to fix this is to force the class not to combine two
      conversions by setting \ref combine_two_conv to false.
      Alternatively, one can ensure that no combination is necessary
      by manually adding the desired combination conversion to the
      cache after it is first computed.

      \note \o2 uses some unit aliases which are not used the the GNU
      or OSX units commands, like "Msun" for the solar mass and adds
      some units not present like "Rschwarz" for the Schwarzchild
      radius of a 1 solar mass black hole.
      
      \note Only the const versions, \ref convert_const and
      \ref convert_ret_const are guaranteed to be thread-safe,
      since they are not allowed to update the unit cache.

      \future Add G=1. 

      \future An in_cache() function to test
      to see if a conversion is currently in the cache. 
  */
  template<class fp_t=double> class convert_units {

#ifndef DOXYGEN_INTERNAL

  protected:

  /// The type for caching unit conversions
  typedef struct {
    /// The input unit
    std::string f;
    /// The output unit
    std::string t;
    /// The conversion factor
    fp_t c;
  } unit_t;

  /// The cache where unit conversions are stored
  std::map<std::string,unit_t,std::greater<std::string> > mcache;
    
  /// The iterator type
  typedef typename std::map<std::string,unit_t,
  std::greater<std::string> >::iterator miter;
      
  /// The const iterator type
  typedef typename std::map<std::string,unit_t,
  std::greater<std::string> >::const_iterator mciter;

  /** \brief Type of a derived unit
   */
  typedef struct der_unit_s {
    /// Unit label
    std::string label;
    /// Power of length
    int m;
    /// Power of mass
    int k;
    /// Power of time
    int s;
    /// Power of temperature
    int K;
    /// Power of current
    int A;
    /// Power of moles
    int mol;
    /// Power of luminous intensity
    int cd;
    /// Value
    fp_t val;
  } der_unit;

  /// Number of SI prefixes
  static const size_t n_prefixes=24;

  /// SI prefix labels
  std::vector<std::string> prefixes;

  /// SI prefix factors
  std::vector<fp_t> prefix_facts;
  
  /** \brief SI-like units which allow prefixes
      
      Set in constructor.
  */
  std::vector<der_unit> SI;

  /** \brief Other units which do not allow prefixing
      
      Set in constructor.
  */
  std::vector<der_unit> other;

  /// \name Flags for natural units
  //@{
  bool c_is_1;
  bool hbar_is_1;
  bool kb_is_1;
  bool G_is_1;
  //@}

  /** \brief Set variables for the calculator object for 
      \ref convert_calc()
  */
  void set_vars(fp_t m, fp_t k, fp_t s, fp_t K, fp_t A, fp_t mol,
		fp_t cd, std::map<std::string, fp_t> &vars,
		bool test_vars=false) const {

    vars.clear();

    // meters
    vars.insert(std::make_pair("m",m));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"m")!=vars.end()) {
	O2SCL_ERR("m already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"m",prefix_facts[i]*m));
    }
    
    // kilograms
    if (vars.find("g")!=vars.end()) {
      O2SCL_ERR("g already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("g",1.0e-3*k));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"g")!=vars.end()) {
	O2SCL_ERR("g already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"g",
				 prefix_facts[i]*1.0e-3*k));
    }

    // seconds
    if (test_vars && vars.find("s")!=vars.end()) {
      O2SCL_ERR("s already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("s",s));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"s")!=vars.end()) {
	O2SCL_ERR("s already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"s",prefix_facts[i]*s));
    }
    
    // Kelvin
    if (vars.find("K")!=vars.end()) {
      O2SCL_ERR("K already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("K",K));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"K")!=vars.end()) {
	O2SCL_ERR("K already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"K",
				 prefix_facts[i]*K));
    }

    // Amperes
    if (vars.find("A")!=vars.end()) {
      O2SCL_ERR("A already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("A",A));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"A")!=vars.end()) {
	O2SCL_ERR("A already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"A",
				 prefix_facts[i]*A));
    }

    // moles
    if (vars.find("mol")!=vars.end()) {
      O2SCL_ERR("mol already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("mol",mol));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"mol")!=vars.end()) {
	O2SCL_ERR("mol already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"mol",
				 prefix_facts[i]*mol));
    }

    // candelas
    if (vars.find("cd")!=vars.end()) {
      O2SCL_ERR("cd already found in list.",o2scl::exc_esanity);
    }
    vars.insert(std::make_pair("cd",cd));
    for(size_t i=0;i<n_prefixes;i++) {
      if (test_vars && vars.find(prefixes[i]+"cd")!=vars.end()) {
	O2SCL_ERR("cd already found in list.",o2scl::exc_esanity);
      }
      vars.insert(std::make_pair(prefixes[i]+"cd",
				 prefix_facts[i]*cd));
    }

    // For SI-like units, add with all of the various prefixes
    for(size_t i=0;i<SI.size();i++) {
      if (vars.find(SI[i].label)!=vars.end()) {
	O2SCL_ERR((((std::string)"SI unit ")+SI[i].label+
		   " already found in list.").c_str(),o2scl::exc_esanity);
      }
      fp_t val=pow(m,SI[i].m)*pow(k,SI[i].k)*pow(s,SI[i].s);
      vars.insert(std::make_pair(SI[i].label,val));
      for(size_t j=0;j<n_prefixes;j++) {
	if (test_vars && vars.find(prefixes[j]+SI[i].label)!=vars.end()) {
	O2SCL_ERR((((std::string)"SI unit ")+prefixes[j]+SI[i].label+
		   " already found in list.").c_str(),o2scl::exc_esanity);
	}
	fp_t val=pow(m,SI[i].m)*pow(k,SI[i].k)*pow(s,SI[i].s)*SI[i].val;
	vars.insert(std::make_pair(prefixes[j]+SI[i].label,
				   val*prefix_facts[j]));
      }
    }
    
    // For other units, just add the value
    for(size_t i=0;i<other.size();i++) {
      if (vars.find(other[i].label)!=vars.end()) {
	O2SCL_ERR((((std::string)"Non SI-unit ")+other[i].label+
		   " already found in list.").c_str(),o2scl::exc_esanity);
      }
      fp_t val=pow(m,other[i].m)*pow(k,other[i].k)*pow(s,other[i].s)*
	other[i].val;
      vars.insert(std::make_pair(other[i].label,val));
    }

    if (false) {
      for (typename std::map<std::string,fp_t>::iterator it=vars.begin();
	   it!=vars.end();it++) {
	std::cout << it->first << " " << it->second << std::endl;
      }
    }
    
    return;
  }
  
  /** \brief The internal conversion function which tries the
      cache first and, if that failed, tries GNU units.

      This function returns 0 if the conversion was successful. If
      the conversion fails and \ref err_on_fail is \c true, then the
      error handler is called. If the conversion fails and \ref
      err_on_fail is \c false, then the value \ref
      o2scl::exc_enotfound is returned.

      The public conversion functions in this class are
      basically just wrappers around this internal function.
  */
  int convert_internal(std::string from, std::string to,
		       fp_t val, fp_t &converted,
		       fp_t &factor, bool &new_conv) const {

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
	std::cout << "convert_units::convert_internal(): "
	<< "Value of use_gnu_units is true." << std::endl;
      }

      int ret_gnu=convert_gnu_units(from,to,val,converted,factor);

      if (ret_gnu==0) {
	new_conv=true;
	return 0;
      }

    } else {
    
      if (verbose>0) {
	std::cout << "convert_units::convert_internal(): "
	<< "Value of use_gnu_units is false." << std::endl;
      }
    
    }    
  
    if (err_on_fail) {
      std::string str=((std::string)"Conversion between ")+from+" and "+to+
      " not found in convert_units::convert_internal().";
      O2SCL_ERR(str.c_str(),exc_enotfound);
    }
  
    return exc_enotfound;
  }

  /** \brief Attempt to use GNU units to perform a conversion

      This function attempts to open a pipe to GNU units independent
      of the value of \ref use_gnu_units. However, it will always
      return a non-zero value if \c HAVE_POPEN is not defined
      signaling that the <tt>popen()</tt> function is not available
      (but does not call the error handler in this case). The
      function returns 0 if the conversion was successful. If
      HAVE_POPEN is defined but the call to GNU units fails
      for some reason, then the error handler is called if
      \ref err_on_fail is true.

      \future Finer control on verbose=1 vs. verbose=2
  */
  int convert_gnu_units(std::string from, std::string to,
			fp_t val, fp_t &converted,
			fp_t &factor) const {

    // Run the GNU 'units' command
    std::string cmd=units_cmd_string+" '"+from+"' '"+to+"'";
    if (verbose>0) {
      std::cout << "convert_units::convert_gnu_units(): "
		<< "Units command is " << cmd << std::endl;
    }
  
#ifdef HAVE_POPEN
  
    if (verbose>0) {
      std::cout << "convert_units::convert_gnu_units(): "
      << "GNU units command is \"" << cmd
      << "\"." << std::endl;
    }
  
    FILE *ps_pipe=popen(cmd.c_str(),"r");
    if (!ps_pipe) {
      if (verbose>0) {
	std::cout << "convert_units::convert_gnu_units(): "
		  << "Pope open failed." << std::endl;
      }
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
		<< "Units output is \"" << line1 << "\"." << std::endl;
    }
  
    // Read the output from the 'units' command and compute the 
    // conversion factor
    std::string s=line1, t1, t2;
    std::istringstream *ins=new std::istringstream(s);
    (*ins) >> t1 >> t2;
    delete ins;
    if (verbose>0) {
      std::cout << "convert_units::convert_gnu_units(): "
		<< "Second word of output is \""
		<< t2 << "\"." << std::endl;
    }
    int sret=o2scl::stod_nothrow(t2,factor);
    if (sret!=0) {
      if (verbose>0) {
	std::cout << "convert_units::convert_gnu_units(): "
		  << "GNU unit conversion failed." << std::endl;
      }
      if (err_on_fail) {
	std::string str=((std::string)"Conversion between ")+from+" and "+to+
	  " not found in convert_units::convert_gnu_units().";
	O2SCL_ERR(str.c_str(),exc_enotfound);
      }
      return exc_enotfound;
    }
    if (verbose>0) {
      std::cout << "convert_units::convert_gnu_units(): "
      << "Converted value is " << factor*val << std::endl;
    }
  
    // Cleanup
    if (pclose(ps_pipe)!=0) {
      if (verbose>0) {
	std::cout << "convert_units::convert_gnu_units(): "
	<< "Close pipe failed." << std::endl;
      }
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

  /** \brief Attempt to construct a conversion from the internal
      unit cache

      This function returns 0 if the conversion was successful and
      \ref o2scl::exc_efailed otherwise. This function
      does not call the error handler.
  */
  int convert_cache(std::string from, std::string to,
		    fp_t val, fp_t &converted,
		    fp_t &factor) const {
  
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

  
#endif

  public:

  /// Create a unit-conversion object
  convert_units() {
    verbose=0;
    use_gnu_units=true;
    units_cmd_string="units";
    err_on_fail=true;
    combine_two_conv=true;

    prefixes={"Q","R","Y","Z","E","P","T","G","M",
	      "k","h","da","d","c","m","mu","n",
	      "p","f","a","z","y","r","q"};
    
    prefix_facts={1.0e30,1.0e27,1.0e24,1.0e21,1.0e18,1.0e15,
		  1.0e12,1.0e9,1.0e6,1.0e3,1.0e2,10.0,0.1,
		  1.0e-2,1.0e-3,1.0e-6,1.0e-9,1.0e-12,1.0e-15,
		  1.0e-18,1.0e-21,1.0e-24,1.0e-27,1.0e-30};
    
    // SI derived units, in order m kg s K A mol cd
    std::vector<der_unit> SI_=
    {{"J",2,1,-2,0,0,0,0,1.0},
     {"N",1,1,-2,0,0,0,0,1.0},
     {"Pa",0,1,-2,0,0,0,0,1.0},
     {"W",2,1,-3,0,0,0,0,1.0},
     {"C",0,1,0,0,1,0,0,1.0},
     {"V",2,1,-3,0,-1,0,0,1.0},
     {"ohm",2,1,-3,0,-2,0,0,1.0},
     {"S",-2,-1,3,0,2,0,0,1.0},
     {"F",-2,-1,4,0,2,0,0,1.0},
     {"Wb",2,1,-2,0,-1,0,0,1.0},
     {"H",2,1,-2,0,-2,0,0,1.0},
     {"T",0,1,-2,0,-1,0,0,1.0},
     {"Hz",0,0,-1,0,0,0,0,1.0},
     {"lx",-1,0,0,0,0,0,1,1.0},
     {"Bq",0,0,-1,0,0,0,0,1.0},
     {"Gy",2,0,-2,0,0,0,0,1.0},
     {"Sv",2,0,-2,0,0,0,0,1.0},
     {"kat",0,0,-1,0,0,1,0,1.0},
     // liters, "l" and "L"
     {"l",3,0,0,0,0,0,0,1.0e-3},
     {"L",3,0,0,0,0,0,0,1.0e-3},
     // metric tons
     {"t",0,1,0,0,0,0,0,1.0e3},
     // Daltons (atomic mass units)
     {"Da",0,1,0,0,0,0,0,o2scl_mks::unified_atomic_mass},
     // Electron volts
     {"eV",2,1,-2,0,0,0,0,o2scl_mks::electron_volt}};
    SI=SI_;
    
    // Other units, in order m kg s K A mol cd
    std::vector<der_unit> other_=
    {

     // Length units
     // We cannot use 'ft' because the metric ton is a unit "accepted
     // for use with the SI" and then ft is a "femto-ton", which is a
     // bit confusing
     {"foot",1,0,0,0,0,0,0,o2scl_mks::foot},
     {"in",1,0,0,0,0,0,0,o2scl_mks::inch},
     {"yd",1,0,0,0,0,0,0,o2scl_mks::yard},
     {"mi",1,0,0,0,0,0,0,o2scl_mks::mile},
     {"nmi",1,0,0,0,0,0,0,o2scl_mks::nautical_mile},
     {"fathom",1,0,0,0,0,0,0,o2scl_mks::fathom},
     {"angstrom",1,0,0,0,0,0,0,o2scl_mks::angstrom},
     {"mil",1,0,0,0,0,0,0,o2scl_mks::mil},
     {"point",1,0,0,0,0,0,0,o2scl_mks::point},
     {"texpoint",1,0,0,0,0,0,0,o2scl_mks::texpoint},
     {"micron",1,0,0,0,0,0,0,o2scl_mks::micron},
     // AU's, "au" and "AU"
     {"AU",1,0,0,0,0,0,0,o2scl_mks::astronomical_unit},
     {"au",1,0,0,0,0,0,0,o2scl_mks::astronomical_unit},
     // light years "ly" and "lyr"
     {"ly",1,0,0,0,0,0,0,o2scl_mks::light_year},
     {"lyr",1,0,0,0,0,0,0,o2scl_mks::light_year},
     // We add common SI-like prefixes for non-SI units
     {"Gpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e9},
     {"Mpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e6},
     {"kpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e3},
     {"pc",1,0,0,0,0,0,0,o2scl_mks::parsec},

     // Area units
     // hectares, "ha" and "hectare"
     {"hectare",2,0,0,0,0,0,0,o2scl_mks::hectare},
     {"ha",2,0,0,0,0,0,0,1.0e4},
     // acre
     {"acre",2,0,0,0,0,0,0,o2scl_mks::acre},
     // barn
     {"barn",2,0,0,0,0,0,0,o2scl_mks::barn},
     
     // Volume units
     {"us_gallon",3,0,0,0,0,0,0,o2scl_mks::us_gallon},
     {"quart",3,0,0,0,0,0,0,o2scl_mks::quart},
     {"pint",3,0,0,0,0,0,0,o2scl_mks::pint},
     {"cup",3,0,0,0,0,0,0,o2scl_mks::cup},
     {"tbsp",3,0,0,0,0,0,0,o2scl_mks::tablespoon},
     {"tsp",3,0,0,0,0,0,0,o2scl_mks::teaspoon},
     {"ca_gallon",3,0,0,0,0,0,0,o2scl_mks::canadian_gallon},
     {"uk_gallon",3,0,0,0,0,0,0,o2scl_mks::uk_gallon},

     // Mass units
     // Solar masses, "Msun", and "Msolar"
     {"Msun",0,1,0,0,0,0,0,o2scl_mks::solar_mass},
     {"Msolar",0,1,0,0,0,0,0,o2scl_mks::solar_mass},
     {"pound",0,1,0,0,0,0,0,o2scl_mks::pound_mass},
     {"ounce",0,1,0,0,0,0,0,o2scl_mks::ounce_mass},
     {"uk_ton",0,1,0,0,0,0,0,o2scl_mks::uk_ton},
     {"troy_ounce",0,1,0,0,0,0,0,o2scl_mks::troy_ounce},
     {"carat",0,1,0,0,0,0,0,o2scl_mks::carat},
     
     // Velocity units
     {"knot",1,0,-1,0,0,0,0,o2scl_mks::knot},
     {"c",1,0,-1,0,0,0,0,o2scl_const::speed_of_light_f<fp_t>()},
     
     // Energy units
     {"cal",2,1,-2,0,0,0,0,o2scl_mks::calorie},
     {"btu",2,1,-2,0,0,0,0,o2scl_mks::btu},
     {"erg",2,1,-2,0,0,0,0,o2scl_mks::erg},

     // Power units
     {"therm",2,1,-3,0,0,0,0,o2scl_mks::therm},
     {"horsepower",2,1,-3,0,0,0,0,o2scl_mks::horsepower},

     // Pressure units
     {"atm",-1,2,-2,0,0,0,0,o2scl_mks::std_atmosphere},
     {"bar",-1,1,-2,0,0,0,0,o2scl_mks::bar},
     {"torr",-1,1,-2,0,0,0,0,o2scl_mks::torr},
     {"psi",-1,1,-2,0,0,0,0,o2scl_mks::psi},

     // Time units
     {"yr",0,0,1,0,0,0,0,31556926},
     {"wk",0,0,1,0,0,0,0,o2scl_mks::week},
     {"d",0,0,1,0,0,0,0,o2scl_mks::day},

     // Hours, "hr", to avoid confusion with Planck's constant
     {"hr",0,0,1,0,0,0,0,o2scl_mks::hour},
     {"min",0,0,1,0,0,0,0,o2scl_mks::minute},

     // Inverse time units
     {"curie",0,0,-1,0,0,0,0,o2scl_mks::curie},
     
     // Force units
     {"dyne",1,1,-2,0,0,0,0,o2scl_mks::dyne},
     
     // Viscosity units
     {"poise",-1,1,-1,0,0,0,0,o2scl_mks::poise},

     // Units of heat capacity or entropy
     // Boltzmann's constant. Could be confused with kilobytes?
     {"kB",2,1,-2,-1,0,0,0,o2scl_const::boltzmann_f<fp_t>()},
     
     {"hbar",2,1,-1,0,0,0,0,o2scl_const::hbar_f<fp_t>()},
     // "Planck" instead of "h", to avoid confusing with hours
     {"Planck",2,1,-1,0,0,0,0,o2scl_const::planck_f<fp_t>()},
     
     // Gravitational constant. We cannot use "G" because of
     // confusion with "Gauss"
     {"GNewton",3,-1,-2,0,0,0,0,o2scl_mks::gravitational_constant},
     
     // Gauss, and note the possible confusion with the gravitational
     // constant
     {"G",0,1,-2,0,-1,0,0,o2scl_mks::gauss},

     {"NA",0,0,0,0,0,0,0,o2scl_const::avogadro}
     
    };
    other=other_;
    
    c_is_1=false;
    hbar_is_1=false;
    kb_is_1=false;
    G_is_1=false;
  }
    
  virtual ~convert_units() {}

  /** \brief Set natural units
   */
  void set_natural_units(bool c_is_one=true, bool hbar_is_one=true,
			 bool kb_is_one=true) {
    if (hbar_is_one && !c_is_one) {
      O2SCL_ERR2("The convert_units class does not support hbar=1,",
		 "when c!=1.",o2scl::exc_einval);
    }
    c_is_1=c_is_one;
    hbar_is_1=hbar_is_one;
    kb_is_1=kb_is_one;
    return;
  }

  /** \brief Test to make sure all units are unique
   */
  void test_unique() {
    std::map<std::string,fp_t> vars;
    set_vars(1.0,1.0,1.0,1.0,1.0,1.0,1.0,vars,true);
    return;
  }
  
  /** \brief Print the units in the data base
   */
  void print_units(std::ostream &out) {

    out << "SI-like:  label  m kg  s  K  A mol cd value" << std::endl;
    out << "--------------- -- -- -- -- -- --- -- ------------"
    << std::endl;
    for(size_t i=0;i<SI.size();i++) {
      out.width(15);
      out << SI[i].label << " ";
      out.width(2);
      out << SI[i].m << " ";
      out.width(2);
      out << SI[i].k << " ";
      out.width(2);
      out << SI[i].s << " ";
      out.width(2);
      out << SI[i].K << " ";
      out.width(2);
      out << SI[i].A << " ";
      out.width(3);
      out << SI[i].mol << " ";
      out.width(2);
      out << SI[i].cd << " ";
      out << SI[i].val << std::endl;
    }
    out << std::endl;
    
    out << "Other:    label  m kg  s  K  A mol cd value:" << std::endl;
    out << "--------------- -- -- -- -- -- --- -- ------------"
    << std::endl;
    
    for(size_t i=0;i<other.size();i++) {
      out.width(15);
      out << other[i].label << " ";
      out.width(2);
      out << other[i].m << " ";
      out.width(2);
      out << other[i].k << " ";
      out.width(2);
      out << other[i].s << " ";
      out.width(2);
      out << other[i].K << " ";
      out.width(2);
      out << other[i].A << " ";
      out.width(3);
      out << other[i].mol << " ";
      out.width(2);
      out << other[i].cd << " ";
      out << other[i].val << std::endl;
    }
    
    return;
  }
  
  /** \brief Automatic unit conversion between SI-based units
      with a \ref o2scl::calculator object
   */
  int convert_calc(std::string from, std::string to,
		   fp_t val, fp_t &converted,
		   fp_t &factor) const {
    
    // These calculator objects have to be inside this
    // function to make the function const
    o2scl::calculator calc;
    o2scl::calculator calc2;

    std::map<std::string, fp_t> vars;

    set_vars(1.0,1.0,1.0,1.0,1.0,1.0,1.0,vars);
    
    int cret1=calc.compile_nothrow(from.c_str());
    if (cret1!=0) return 1;
    int cret2=calc2.compile_nothrow(to.c_str());
    if (cret2!=0) return 2;
    
    fp_t before, after;
    int cret3=calc.eval_nothrow(&vars,before);
    if (cret3!=0) return 3;
    int cret4=calc2.eval_nothrow(&vars,after);
    if (cret4!=0) return 4;
    
    factor=before/after;
    converted=val*factor;

    // Now, having verified that a conversion is possible, we need to
    // separately verify the conversion is sensible by rescaling the
    // base units to make sure sure both sides scale the same way. The
    // prevents, e.g. conversions between m and 1/m. We initialize
    // factor_m to avoid uninitialized warnings.

    fp_t factor_m, factor_s, factor_kg, factor_K, addl=1;

    if (c_is_1) {

      if (hbar_is_1) {

	if (kb_is_1) {
	  
	  set_vars(2.0,1.0,4.0,1.0,1.0,1.0,1.0,vars);
	  fp_t factor_1=calc.eval(&vars)/calc2.eval(&vars);
	  std::cout << std::endl;
	  std::cout << "0. " << calc.eval(&vars) << " "
		    << calc2.eval(&vars) << std::endl;
	  
	  set_vars(2.0,1.0,1.0,4.0,1.0,1.0,1.0,vars);
	  fp_t factor_2=calc.eval(&vars)/calc2.eval(&vars);
	  std::cout << "0. " << calc.eval(&vars) << " "
		    << calc2.eval(&vars) << std::endl;
	  
	  set_vars(1.0,1.0,4.0,2.0,1.0,1.0,1.0,vars);
	  fp_t factor_3=calc.eval(&vars)/calc2.eval(&vars);
	  std::cout << "0. " << calc.eval(&vars) << " "
		    << calc2.eval(&vars) << std::endl;

	  fp_t exp1=log(factor_1/factor)/log(2.0);
	  fp_t exp2=log(factor_2/factor)/log(2.0);
	  fp_t exp3=log(factor_3/factor)/log(2.0);

	  fp_t hbar_fac=-(5*exp1+exp2-2*exp3)/6;
	  fp_t c_fac=-(-5*exp1+2*exp2+2*exp3)/3;
	  fp_t boltz_fac=-(exp1-exp2-exp3)/3;

	  std::cout << "1. " << exp1 << " " << exp2 << " " << exp3
		    << std::endl;
	  std::cout << "2. " << hbar_fac << " " << c_fac << " "
		    << boltz_fac << std::endl;

	  addl=pow(o2scl_const::hbar_f<fp_t>(),hbar_fac);
	  addl*=pow(o2scl_const::speed_of_light_f<fp_t>(),c_fac);
	  std::cout << "3. " << factor << " " << factor_1 << " "
		    << factor_2 << " " << factor_3 << std::endl;
	  addl*=pow(o2scl_mks::boltzmann,boltz_fac);
	  std::cout << "4. " << addl << " " << factor << std::endl;
	  
	  // Then set factor_s, factor_kg, factor_K equal to factor_m
	  // so the test below succeeds
	  factor_m=factor;
	  factor_s=factor;
	  factor_kg=factor;
	  factor_K=factor;
	  
	} else {
	
	  // Compute factor_K
	  set_vars(1.0,1.0,1.0,2.0,1.0,1.0,1.0,vars);
	  //vars["c"]=1.0;
	  //vars["hbar"]=1.0;
	  factor_K=calc.eval(&vars)/calc2.eval(&vars);
	
	  // Scale m, kg, and s at the same time
	  set_vars(2.0,0.5,2.0,1.0,1.0,1.0,1.0,vars);
	  //vars["c"]=1.0;
	  //vars["hbar"]=1.0;
	  fp_t factor_1=calc.eval(&vars)/calc2.eval(&vars);

	  if (false) {
	    std::cout << std::endl;
	    std::cout << "1. " << from << " " << to << std::endl;
	    std::cout << "2. " << vars["kg"] << " "
		      << vars["MeV"] << " " << vars["c"] << std::endl;
	    std::cout << "3. " << calc.eval(&vars) << " "
		      << calc2.eval(&vars) << std::endl;
	  }
	  
	  // Separately compute factor_s and factor_kg in order to
	  // determine how many factors of hbar and c we need
	  set_vars(1.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
	  //vars["c"]=0.5;
	  //vars["hbar"]=0.5;
	  fp_t factor_2=calc.eval(&vars)/calc2.eval(&vars);
	  
	  set_vars(1.0,2.0,1.0,1.0,1.0,1.0,1.0,vars);
	  fp_t factor_3=calc.eval(&vars)/calc2.eval(&vars);
	  
	  double exp1=log(factor_1/factor_2)/log(2.0);
	  double exp2=log(factor_3*factor_2/factor/factor)/log(2.0);

	  if (false) {
	    std::cout << "4. " << factor_1 << " " << factor_2 << " "
		      << factor_3 << std::endl;
	    std::cout << "5. " << exp1 << " " << exp2 << " "
		      << -exp1-exp2 << std::endl;
		      
	  }

	  addl=pow(o2scl_const::speed_of_light_f<fp_t>(),exp2);
	  addl*=pow(o2scl_const::hbar_f<fp_t>(),-exp1-exp2);
	  
	  //std::cout << "6. " << addl << " " << factor << std::endl;
	  
	  //std::cout << "\nexps: " << exp1 << " " << exp2 << " "
	  //<< exp2 << " " << -exp1-exp2 << std::endl;
	  
	  // Then set factor_s and factor_3 equal to factor_m so the
	  // test below succeeds
	  factor_m=factor;
	  factor_s=factor;
	  factor_kg=factor;

	}
      
      } else {

	// Compute factor_kg
	set_vars(1.0,2.0,1.0,1.0,1.0,1.0,1.0,vars);
	factor_kg=calc.eval(&vars)/calc2.eval(&vars);
    
	// Compute factor_K
	set_vars(1.0,1.0,1.0,2.0,1.0,1.0,1.0,vars);
	factor_K=calc.eval(&vars)/calc2.eval(&vars);
	
	// Scale m and s at the same time
	set_vars(2.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
	factor_m=calc.eval(&vars)/calc2.eval(&vars);
	
	// Separately compute factor_s in order to determine
	// how many factors of c we need
	set_vars(1.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
	factor_s=calc.eval(&vars)/calc2.eval(&vars);
	
	addl=pow(o2scl_mks::speed_of_light,
		 log(factor_s/factor_m)/log(2.0));
	
	// Then set factor_s equal to factor_m so the test below
	// succeeds
	factor_s=factor_m;

      }
      
    } else {
      
      // Compute factor_kg
      set_vars(1.0,2.0,1.0,1.0,1.0,1.0,1.0,vars);
      factor_kg=calc.eval(&vars)/calc2.eval(&vars);
      
      set_vars(2.0,1.0,1.0,1.0,1.0,1.0,1.0,vars);
      factor_m=calc.eval(&vars)/calc2.eval(&vars);
      
      set_vars(1.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
      factor_s=calc.eval(&vars)/calc2.eval(&vars);

      // Compute factor_K
      set_vars(1.0,1.0,1.0,2.0,1.0,1.0,1.0,vars);
      factor_K=calc.eval(&vars)/calc2.eval(&vars);
      
    }

    set_vars(1.0,1.0,1.0,1.0,2.0,1.0,1.0,vars);
    fp_t factor_A=calc.eval(&vars)/calc2.eval(&vars);

    set_vars(1.0,1.0,1.0,1.0,1.0,2.0,1.0,vars);
    fp_t factor_mol=calc.eval(&vars)/calc2.eval(&vars);

    set_vars(1.0,1.0,1.0,1.0,1.0,1.0,2.0,vars);
    fp_t factor_cd=calc.eval(&vars)/calc2.eval(&vars);

    if (verbose>1) {
      std::cout << "from: " << from << " to: " << to
		<< " factor: " << factor << " factor_m: " << factor_m
		<< "\n\tfactor_kg: " << factor_kg << " factor_s: " << factor_s 
		<< " factor_K: " << factor_K << "\n\tfactor_A: " << factor_A 
		<< " factor_mol: " << factor_mol
		<< " factor_cd: " << factor_cd << std::endl;
    }    
    
    if (fabs(factor/factor_m-1.0)<1.0e-14 &&
	fabs(factor/factor_kg-1.0)<1.0e-14 &&
	fabs(factor/factor_s-1.0)<1.0e-14 &&
	fabs(factor/factor_K-1.0)<1.0e-14 &&
	fabs(factor/factor_A-1.0)<1.0e-14 &&
	fabs(factor/factor_mol-1.0)<1.0e-14 &&
	fabs(factor/factor_cd-1.0)<1.0e-14) {
      factor*=addl;
      converted=factor*val;
      return 0;
    }
    
    return 5;
  }
  
  /// \name Basic usage
  //@{
  /** \brief Return the value \c val after converting using units \c
      from and \c to
  */
  virtual fp_t convert(std::string from, std::string to, fp_t val) {
    fp_t converted;
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
      std::string str=((std::string)"Conversion between ")+from+" and "+to+
      " not found in convert_units::convert().";
      O2SCL_ERR(str.c_str(),exc_enotfound);
    }
    return converted;
  }    

  /** \brief Return the value \c val after converting using units \c
      from and \c to (const version)
  */
  virtual fp_t convert_const(std::string from, std::string to,
			     fp_t val) const {
    fp_t converted;
    int ret=convert_ret_const(from,to,val,converted);
    if (ret==exc_enotfound) {
      std::string str=((std::string)"Conversion between ")+from+" and "+to+
	" not found in convert_units::convert().";
      O2SCL_ERR(str.c_str(),exc_enotfound);
    }
    return converted;
  }      
  //@}

  /// \name User settings
  //@{
  /// Verbosity (default 0)
  int verbose;

  /** \brief If true, use a system call to units to derive new
      conversions (default true)

      This also requires <tt>popen()</tt>.
  */
  bool use_gnu_units;

  /// If true, throw an exception when a conversion fails (default true)
  bool err_on_fail;

  /// If true, allow combinations of two conversions (default true)
  bool combine_two_conv;

  /// Command string to call units (default "units")
  std::string units_cmd_string;
  //@}

  /// \name Conversions which don't throw exceptions
  //@{
  /** \brief Return the value \c val after converting using units \c
      from and \c to, returning a non-zero value on failure
  */
  virtual int convert_ret(std::string from, std::string to, fp_t val,
			  fp_t &converted) {
    fp_t factor;
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

  /** \brief Return the value \c val after converting using units \c
      from and \c to, returning a non-zero value on failure
      (const version)
  */
  virtual int convert_ret_const(std::string from, std::string to,
				fp_t val, fp_t &converted) const {
    fp_t factor;
    bool new_conv;
  
    return convert_internal(from,to,val,converted,factor,new_conv);
  }
      
  //@}

  /// \name Manipulate cache and create units.dat files
  //@{
  /// Manually insert a unit conversion into the cache
  void insert_cache(std::string from, std::string to, fp_t conv) {

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

  /// Print the present unit cache to std::cout
  void print_cache() const {
    mciter m;
    if (mcache.size()==0) {
      std::cout << "No units in cache." << std::endl;
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

  /** \brief Desc
   */
  void test_cache_calc() const {
    mciter m;
    for (m=mcache.begin();m!=mcache.end();m++) {
      std::cout.setf(std::ios::left);
      std::cout.width(10);
      std::cout << m->second.f << " ";
      std::cout.width(10);
      std::cout << m->second.t << " ";
      std::cout.unsetf(std::ios::left);
      std::cout.precision(6);
      std::cout << m->second.c << " ";
      fp_t c, f;
      int ix=convert_calc(m->second.f,m->second.t,1.0,c,f);
      std::cout << ix << " " << f << std::endl;
    }
    return;
  }
  
  /// Manually remove a unit conversion into the cache
  void remove_cache(std::string from, std::string to) {

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
      O2SCL_ERR((((std::string)"Conversion ")+from+" -> "+to+
		 " not found in convert_units::remove_cache().").c_str(),
		exc_enotfound);
    }

    return;
  }

  /** \brief Add default conversions

      Where possible, this uses templates from constants.h to define
      the conversions exactly in the user-specified floating-point
      type.
  */
  void default_conversions() {
    
    // Default conversions are given here. Obviously GNU units is better
    // at handling these things, but it's nice to have some of the easy
    // conversions in by default rather than worrying about opening a
    // pipe, etc.
    
    fp_t sol_mks=o2scl_const::speed_of_light_f<fp_t>(o2scl_const::o2scl_mks);
    fp_t sol_cgs=o2scl_const::speed_of_light_f<fp_t>(o2scl_const::o2scl_cgs);
    fp_t hc=o2scl_const::hc_mev_fm_f<fp_t>();
    fp_t elem_charge=o2scl_const::elem_charge_f<fp_t>();
    
    // hbar=c=1 conversion from mass to inverse length
    insert_cache("kg","1/fm",1.0e-15/
		 o2scl_const::hbar_f<fp_t>(o2scl_const::o2scl_mks)*
		 sol_mks);
  
    // Simple mass/energy conversions with c^2=1
    
    insert_cache("kg","MeV",pow(sol_mks,2.0)/elem_charge*1.0e-6);
    insert_cache("kg","g",1.0e3);
    insert_cache("PeV","eV",1.0e15);
    insert_cache("TeV","eV",1.0e12);
    insert_cache("GeV","eV",1.0e9);
    insert_cache("MeV","eV",1.0e6);
    insert_cache("keV","eV",1.0e3);
    insert_cache("meV","eV",1.0e-3);
    insert_cache("Msun","kg",o2scl_mks::solar_mass);
    insert_cache("erg","kg",o2scl_mks::erg/pow(sol_mks,2.0));

    // Joules and Kelvin
  
    insert_cache("eV","J",elem_charge);
    insert_cache("K","J",o2scl_mks::boltzmann);
    insert_cache("K","kg",o2scl_mks::boltzmann/pow(sol_mks,2.0));

    // Energy density and pressure conversions

    insert_cache("atm","bar",o2scl_mks::std_atmosphere/o2scl_mks::bar);
    insert_cache("atm","Pa",o2scl_mks::std_atmosphere);
    insert_cache("kPa","Pa",1.0e3);
    insert_cache("Pa","kg/m^3",1.0/sol_mks/sol_mks);
    insert_cache("Pa","g/cm^3",10.0/sol_cgs/sol_cgs);
    insert_cache("Pa","MeV/fm^3",1.0e-44/o2scl_cgs::electron_volt);
    insert_cache("Pa","erg/cm^3",10.0);
    insert_cache("g/cm^3","Msun/km^3",1.0e12/o2scl_mks::solar_mass);
    insert_cache("erg/cm^3","Msun/km^3",1.0e12/sol_cgs/
		 sol_cgs/o2scl_mks::solar_mass);
    insert_cache("dyne/cm^2","Msun/km^3",1.0e12/sol_cgs/
		 sol_cgs/o2scl_mks::solar_mass);
    insert_cache("MeV/fm^3","Msun/km^3",
		 o2scl_cgs::electron_volt/sol_cgs/
		 sol_cgs/o2scl_mks::solar_mass*1.0e57);
    insert_cache("1/fm^4","Msun/km^3",hc*
		 o2scl_cgs::electron_volt/sol_cgs/
		 sol_cgs/o2scl_mks::solar_mass*1.0e57);

    // 1/fm^4 conversions using hbar*c
    
    insert_cache("1/fm^4","MeV/fm^3",hc);
    insert_cache("MeV/fm^3","MeV^2/fm^2",hc);
    insert_cache("MeV^2/fm^2","MeV^3/fm",hc);
    insert_cache("MeV^3/fm","MeV^4",hc);

    insert_cache("1/fm^4","MeV^2/fm^2",hc*hc);
    insert_cache("MeV/fm^3","MeV^3/fm",hc*hc);
    insert_cache("MeV^2/fm^2","MeV^4",hc*hc);
  
    insert_cache("1/fm^4","MeV^3/fm",hc*hc*hc);
    insert_cache("MeV/fm^3","MeV^4",hc*hc*hc);
  
    insert_cache("1/fm^4","MeV^4",pow(hc,4.0));
  
    // 1/fm^3 conversions using hbar*c

    insert_cache("1/fm^3","MeV/fm^2",hc);
    insert_cache("MeV/fm^2","MeV^2/fm",hc);
    insert_cache("MeV^2/fm","MeV^3",hc);

    insert_cache("1/fm^3","MeV^2/fm",hc*hc);
    insert_cache("MeV/fm^2","MeV^3",hc*hc);

    insert_cache("1/fm^3","MeV^3",pow(hc,3.0));

    // MeV*fm conversions

    insert_cache("MeV*fm","MeV*cm",1.0e-13);
    insert_cache("MeV*fm","MeV*m",1.0e-15);
    insert_cache("MeV*fm","eV*nm",1.0);
    insert_cache("MeV*fm","MeV*s",1.0e-15/o2scl_mks::speed_of_light);
    insert_cache("MeV*fm","erg*cm",elem_charge);
    insert_cache("MeV*fm","erg*s",elem_charge/o2scl_cgs::speed_of_light);
    insert_cache("MeV*fm","J*m",elem_charge/1.0e9);
    insert_cache("MeV*fm","J*s",elem_charge/1.0e9/o2scl_mks::speed_of_light);

    // Simple time conversions
    
    insert_cache("yr","s",31556926);
    insert_cache("wk","s",o2scl_mks::week);
    insert_cache("d","s",o2scl_mks::day);
    insert_cache("hr","s",o2scl_mks::hour);
    insert_cache("min","s",o2scl_mks::minute);

    // Simple powers of length conversions 

    insert_cache("AU","m",o2scl_mks::astronomical_unit);
    insert_cache("pc","m",o2scl_mks::parsec);
    insert_cache("kpc","m",o2scl_mks::parsec*1.0e3);
    insert_cache("km","m",1.0e3);
    insert_cache("cm","m",1.0e-2);
    insert_cache("fm","m",1.0e-15);
    insert_cache("mm","m",1.0e-3);
    insert_cache("nm","m",1.0e-9);
    insert_cache("lyr","m",o2scl_mks::light_year);
  
    insert_cache("AU^2","m^2",pow(o2scl_mks::astronomical_unit,2.0));
    insert_cache("pc^2","m^2",pow(o2scl_mks::parsec,2.0));
    insert_cache("kpc^2","m^2",pow(o2scl_mks::parsec*1.0e3,2.0));
    insert_cache("km^2","m^2",1.0e6);
    insert_cache("cm^2","m^2",1.0e-4);
    insert_cache("fm^2","m^2",1.0e-30);
    insert_cache("mm^2","m^2",1.0e-6);
    insert_cache("nm^2","m^2",1.0e-18);
    insert_cache("lyr^2","m^2",pow(o2scl_mks::light_year,2.0));
  
    insert_cache("AU^3","m^3",pow(o2scl_mks::astronomical_unit,3.0));
    insert_cache("pc^3","m^3",pow(o2scl_mks::parsec,3.0));
    insert_cache("kpc^3","m^3",pow(o2scl_mks::parsec*1.0e3,3.0));
    insert_cache("km^3","m^3",1.0e9);
    insert_cache("cm^3","m^3",1.0e-6);
    insert_cache("fm^3","m^3",1.0e-45);
    insert_cache("mm^3","m^3",1.0e-9);
    insert_cache("nm^3","m^3",1.0e-27);
    insert_cache("lyr^3","m^3",pow(o2scl_mks::light_year,3.0));

    // Simple inverse powers of length conversions 

    insert_cache("1/m","1/AU",o2scl_mks::astronomical_unit);
    insert_cache("1/m","1/pc",o2scl_mks::parsec);
    insert_cache("1/m","1/kpc",o2scl_mks::parsec*1.0e3);
    insert_cache("1/m","1/km",1.0e3);
    insert_cache("1/m","1/cm",1.0e-2);
    insert_cache("1/m","1/fm",1.0e-15);
    insert_cache("1/m","1/mm",1.0e-3);
    insert_cache("1/m","1/nm",1.0e-9);
    insert_cache("1/m","1/lyr",o2scl_mks::light_year);
    
    insert_cache("1/m^2","1/AU^2",pow(o2scl_mks::astronomical_unit,2.0));
    insert_cache("1/m^2","1/pc^2",pow(o2scl_mks::parsec,2.0));
    insert_cache("1/m^2","1/kpc^2",pow(o2scl_mks::parsec*1.0e3,2.0));
    insert_cache("1/m^2","1/km^2",1.0e6);
    insert_cache("1/m^2","1/cm^2",1.0e-4);
    insert_cache("1/m^2","1/fm^2",1.0e-30);
    insert_cache("1/m^2","1/mm^2",1.0e-6);
    insert_cache("1/m^2","1/nm^2",1.0e-18);
    insert_cache("1/m^2","1/lyr^2",pow(o2scl_mks::light_year,2.0));
    
    insert_cache("1/m^3","1/AU^3",pow(o2scl_mks::astronomical_unit,3.0));
    insert_cache("1/m^3","1/pc^3",pow(o2scl_mks::parsec,3.0));
    insert_cache("1/m^3","1/kpc^3",pow(o2scl_mks::parsec*1.0e3,3.0));
    insert_cache("1/m^3","1/km^3",1.0e9);
    insert_cache("1/m^3","1/cm^3",1.0e-6);
    insert_cache("1/m^3","1/fm^3",1.0e-45);
    insert_cache("1/m^3","1/mm^3",1.0e-9);
    insert_cache("1/m^3","1/nm^3",1.0e-27);
    insert_cache("1/m^3","1/lyr^3",pow(o2scl_mks::light_year,3.0));
    
    return;
  }
  
  /** \brief Make a GNU \c units.dat file from the GSL constants

      If \c c_1 is true, then the second is defined in terms of
      meters so that the speed of light is unitless. If \c hbar_1 is
      true, then the kilogram is defined in terms of <tt>s/m^2</tt>
      so that \f$ \hbar \f$ is unitless.

      \note While convert() generally works with the OSX version
      of 'units', the OSX version can't read units.dat files 
      created by this function.

      \note Not all of the GSL constants or the canonical GNU units 
      conversions are given here.
  */
  void make_units_dat(std::string fname, bool c_1=false, 
		      bool hbar_1=false, bool K_1=false) const {
  
    std::ofstream fout(fname.c_str());
    fout.precision(14);

    fp_t sol_mks=o2scl_const::speed_of_light_f<fp_t>(o2scl_const::o2scl_mks);
    fp_t elem_charge=o2scl_const::elem_charge_f<fp_t>();
    
#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ Fundamental units" << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# Fundamental units" << std::endl;
#endif
    fout << "m\t!" << std::endl;
    fout << "meter\tm" << std::endl;
    if (c_1==false) {
      fout << "s\t!" << std::endl;
      fout << "second\ts" << std::endl;
    } else {
      fout << "s\t" << sol_mks << " m" << std::endl;
      fout << "second\ts" << std::endl;
    }
    if (hbar_1==false) {
      fout << "kg\t!" << std::endl;
      fout << "kilogram\tkg" << std::endl;
    } else {
      fout << "kg\t" << 1.0/o2scl_const::hbar_f<fp_t>(o2scl_const::o2scl_mks) 
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
      
#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# SI and common prefixes" << std::endl;
#endif
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
    fout << "fp_t-\t\t2" << std::endl;
    fout << "triple-\t\t3" << std::endl;
    fout << "treble-\t\t3" << std::endl;
    fout << std::endl;

#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# SI prefix abbreviations" << std::endl;
#endif
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

#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# Basic numbers" << std::endl;
#endif
    fout << "one                     1" << std::endl;
    fout << "two                     2" << std::endl;
    fout << "fp_t                  2" << std::endl;
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

#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# Basic SI units" << std::endl;
#endif
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

#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# Dimensional analysis units" << std::endl;
#endif
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
      
#ifdef O2SCL_OSX
    fout << "/ ----------------------------------------------" 
    << "--------------" << std::endl;
    fout << "/ " << std::endl;
#else
    fout << "################################################" 
    << "##############" << std::endl;
    fout << "# GSL constants" << std::endl;
#endif
    fout << "schwarzchild_radius\t\t" << o2scl_mks::schwarzchild_radius
    << " m" << std::endl;
    fout << "Rschwarz\t\tschwarzchild_radius" << std::endl;
    fout << "speed_of_light\t\t" << sol_mks
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
    << o2scl_const::hbar_f<fp_t>(o2scl_const::o2scl_mks)
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
    fout << "electron_volt\t\t" << elem_charge
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
    
  /** \brief Exhaustive test the cache against GNU units
   */
  int test_cache() {
    err_on_fail=false;
    mciter m, m2;
    std::cout << "units_cmd_string: " << units_cmd_string << std::endl;
    for (m=mcache.begin();m!=mcache.end();m++) {
      for (m2=m;m2!=mcache.end();m2++) {
	std::string from=m->second.f;
	std::string to=m2->second.t;
	if (from!=to) {
	  fp_t v=1.0, c, f1, f2=0.0;
	  int cret=convert_cache(from,to,v,c,f1);
	  if (cret==0) {
	    bool new_conv;
	    int gret=convert_gnu_units(from,to,v,c,f2);
	    if (gret==0) {
	      if (fabs(f1-f2)/f1>1.0e-6) {
		std::cout << "* ";
	      } else {
		std::cout << "  ";
	      }
	      std::cout.width(10);
	      std::cout << from << " ";
	      std::cout.width(10);
	      std::cout << to << " " << f1 << " " << f2 << " "
			<< fabs(f1-f2)/f1 << std::endl;
	    }
	  }
	}
	to=m2->second.f;
	if (from!=to) {
	  fp_t v=1.0, c, f1, f2=0.0;
	  int cret=convert_cache(from,to,v,c,f1);
	  if (cret==0) {
	    bool new_conv;
	    int gret=convert_gnu_units(from,to,v,c,f2);
	    if (gret==0) {
	      if (fabs(f1-f2)/f1>1.0e-6) {
		std::cout << "* ";
	      } else {
		std::cout << "  ";
	      }
	      std::cout.width(10);
	      std::cout << from << " ";
	      std::cout.width(10);
	      std::cout << to << " " << f1 << " " << f2 << " "
			<< fabs(f1-f2)/f1 << std::endl;
	    }
	  }
	}
      }
    }
    return 0;
  }
  //@}
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
