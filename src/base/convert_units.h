/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2021, Andrew W. Steiner
  
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
#include <vector>
#include <map>

#include <o2scl/misc.h>
#include <o2scl/constants.h>
#include <o2scl/string_conv.h>
#include <o2scl/shunting_yard.h>
#include <o2scl/calc_utf8.h>
#include <o2scl/find_constants.h>
#include <o2scl/vector.h>

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
      for the requested conversion first. If this is successful, then the
      conversion factor is returned and the conversion is added to the
      cache.

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

      \note Only the const functions, \ref convert_const and
      \ref convert_ret_const are guaranteed to be thread-safe,
      since they are not allowed to update the unit cache.

      \note This class is designed to allow for higher-precision
      conversions, but of course not all of the unit conversions are
      known to high precision.

      \future Add G=1. 

      \future An in_cache() function to test
      to see if a conversion is currently in the cache. 
  */
  template<class fp_t=double> class convert_units {

  public:
    
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
      /// Long name
      std::string name;
    } der_unit;

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
    std::map<std::string,unit_t,std::less<std::string> > mcache;
    
    /// The iterator type
    typedef typename std::map<std::string,unit_t,
                              std::less<std::string> >::iterator miter;
      
    /// The const iterator type
    typedef typename std::map<std::string,unit_t,
                              std::less<std::string> >::const_iterator
    mciter;

    /// Number of SI prefixes
    static const size_t n_prefixes=21;

    /// SI prefix labels
    std::vector<std::string> prefixes;

    /// SI prefix labels
    std::vector<std::string> prefix_names;

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
    
      // (kilo)grams
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
        fp_t val=pow(m,SI[i].m)*pow(k,SI[i].k)*pow(s,SI[i].s)*SI[i].val;
        vars.insert(std::make_pair(SI[i].label,val));
        for(size_t j=0;j<n_prefixes;j++) {
          if (test_vars && vars.find(prefixes[j]+SI[i].label)!=vars.end()) {
            O2SCL_ERR((((std::string)"SI unit ")+prefixes[j]+SI[i].label+
                       " already found in list.").c_str(),o2scl::exc_esanity);
          }
          val=pow(m,SI[i].m)*pow(k,SI[i].k)*pow(s,SI[i].s)*SI[i].val;
          if (false && SI[i].label=="eV") {
            std::cout << "X1: " << SI[i].val << " "
                      << prefixes[j]+SI[i].label << " "
                      << pow(m,SI[i].m)*pow(k,SI[i].k)*pow(s,SI[i].s) << " "
                      << val << " " << prefix_facts[j] << " "
                      << std::endl;
          }
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

      if (verbose>=3) {
        std::cout << "vars: " << std::endl;
        for (typename std::map<std::string,fp_t>::iterator it=vars.begin();
             it!=vars.end();it++) {
          std::cout << it->first << " " << it->second << std::endl;
        }
        std::cout << std::endl;
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
        if (verbose>=2) {
          std::cout << "Function convert_units::convert_internal(): "
                    << "found conversion in cache." << std::endl;
        }
        new_conv=false;
        return 0;
      } else if (verbose>=2) {
        std::cout << "Function convert_units::convert_internal(): "
                  << "did not find conversion in cache." << std::endl;
      }

      // Compute conversion from convert_calc()
      int ret=convert_calc_hck(from,to,val,converted,factor);
      if (ret==0) {
        if (verbose>=2) {
          std::cout << "Function convert_units::convert_internal(): "
                    << "calculated conversion." << std::endl;
        }
        new_conv=true;
        return 0;
      } else if (verbose>=3) {
        std::cout << "Function convert_units::convert_internal(): "
                  << "failed to calculate conversion." << std::endl;
      }
      
      if (err_on_fail) {
        std::string str=((std::string)"Conversion between ")+from+" and "+to+
          " not found in convert_units::convert_internal().";
        O2SCL_ERR(str.c_str(),exc_enotfound);
      }
  
      return exc_enotfound;
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
        if (verbose>=2) {
          std::cout << "Function convert_units::convert_cache(): "
                    << "found forward conversion." << std::endl;
        }
        return 0;
      }

      // Look in cache for reverse conversion
      std::string both2=to+","+from;
      m3=mcache.find(both2);
      if (m3!=mcache.end()) {
        factor=1.0/m3->second.c;
        converted=val*factor;
        if (verbose>=2) {
          std::cout << "Function convert_units::convert_cache(): "
                    << "found inverse conversion." << std::endl;
        }
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
      units_cmd_string="units";
      err_on_fail=true;
      combine_two_conv=true;

      //prefixes={"Q","R","Y","Z","E","P","T","G","M",
      //"k","h","da","d","c","m","mu","μ","n",
      //"p","f","a","z","y","r","q"};
      prefixes={"Y","Z","E","P","T",
                "G","M","k","h","da",
                "d","c","m","mu","μ",
                "n","p","f","a","z",
                "y"};
    
      prefix_facts={1.0e24,1.0e21,1.0e18,1.0e15,1.0e12,
                    1.0e9,1.0e6,1.0e3,1.0e2,10.0,
                    0.1,1.0e-2,1.0e-3,1.0e-6,1.0e-6,
                    1.0e-9,1.0e-12,1.0e-15,1.0e-18,1.0e-21,
                    1.0e-24};

      prefix_names={"yotta","zetta","exa","peta","tera",
                    "giga","mega","kilo","hecto","deka",
                    "deci","centi","milli","micro","micro",
                    "nano","pico","femto","atto","zepto",
                    "yocto"};
    
      // SI derived units, in order m, kg, s, K, A, mol, and cd. Note
      // that, according to the SI page on Wikipedia, "newton" is left
      // lowercase even though it is named after a person.
      std::vector<der_unit> SI_=
        {{"rad",0,0,0,0,0,0,0,1.0,"radian"},
         {"sr",0,0,0,0,0,0,0,1.0,"steradian"},
         {"J",2,1,-2,0,0,0,0,1.0,"joule"},
         {"N",1,1,-2,0,0,0,0,1.0,"newton"},
         {"Pa",-1,1,-2,0,0,0,0,1.0,"pascal"},
         {"W",2,1,-3,0,0,0,0,1.0,"watt"},
         {"C",0,1,0,0,1,0,0,1.0,"coulomb"},
         {"V",2,1,-3,0,-1,0,0,1.0,"volt"},
         {"ohm",2,1,-3,0,-2,0,0,1.0,"ohm"},
         {"Ω",2,1,-3,0,-2,0,0,1.0,"ohm"},
         {"S",-2,-1,3,0,2,0,0,1.0,"siemens"},
         {"F",-2,-1,4,0,2,0,0,1.0,"farad"},
         {"Wb",2,1,-2,0,-1,0,0,1.0,"weber"},
         {"H",2,1,-2,0,-2,0,0,1.0,"henry"},
         {"T",0,1,-2,0,-1,0,0,1.0,"tesla"},
         {"Hz",0,0,-1,0,0,0,0,1.0,"hertz"},
         {"lm",0,0,0,0,0,0,1,1.0,"lumen"},
         {"lx",-2,0,0,0,0,0,1,1.0,"lux"},
         {"Bq",0,0,-1,0,0,0,0,1.0,"becquerel"},
         {"Gy",2,0,-2,0,0,0,0,1.0,"gray"},
         {"Sv",2,0,-2,0,0,0,0,1.0,"sievert"},
         {"kat",0,0,-1,0,0,1,0,1.0,"katal"},
         // liter: "l" and "L"
         {"l",3,0,0,0,0,0,0,1.0e-3,"liter"},
         {"L",3,0,0,0,0,0,0,1.0e-3,"liter"},
         // Daltons (atomic mass units)
         {"Da",0,1,0,0,0,0,0,o2scl_mks::unified_atomic_mass,"dalton"},
         // Electron volts
         {"eV",2,1,-2,0,0,0,0,o2scl_mks::electron_volt,"electron volt"}};
      
      SI=SI_;
    
      // Other units, in order m kg s K A mol cd
      std::vector<der_unit> other_=
        {

         // Length units
         {"ft",1,0,0,0,0,0,0,o2scl_mks::foot,"foot"},
         {"foot",1,0,0,0,0,0,0,o2scl_mks::foot,"foot"},
         {"in",1,0,0,0,0,0,0,o2scl_mks::inch,"inch"},
         {"yd",1,0,0,0,0,0,0,o2scl_mks::yard,"yard"},
         {"mi",1,0,0,0,0,0,0,o2scl_mks::mile,"mile"},
         {"nmi",1,0,0,0,0,0,0,o2scl_mks::nautical_mile,"nautical mile"},
         {"fathom",1,0,0,0,0,0,0,o2scl_mks::fathom,"fathom"},
         {"angstrom",1,0,0,0,0,0,0,o2scl_mks::angstrom,"angstrom"},
         {"mil",1,0,0,0,0,0,0,o2scl_mks::mil,"mil"},
         {"point",1,0,0,0,0,0,0,o2scl_mks::point,"point"},
         {"texpoint",1,0,0,0,0,0,0,o2scl_mks::texpoint,"texpoint"},
         {"micron",1,0,0,0,0,0,0,o2scl_mks::micron,"micron"},
         // AU's: "au" and "AU"
         {"AU",1,0,0,0,0,0,0,o2scl_mks::astronomical_unit,"astronomical unit"},
         {"au",1,0,0,0,0,0,0,o2scl_mks::astronomical_unit,"astronomical unit"},
         // light years: "ly" and "lyr"
         {"ly",1,0,0,0,0,0,0,o2scl_mks::light_year,"light year"},
         {"lyr",1,0,0,0,0,0,0,o2scl_mks::light_year,"light year"},
         // We add common SI-like prefixes for non-SI units
         {"Gpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e9,"gigaparsec"},
         {"Mpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e6,"megaparsec"},
         {"kpc",1,0,0,0,0,0,0,o2scl_mks::parsec*1.0e3,"kiloparsec"},
         {"pc",1,0,0,0,0,0,0,o2scl_mks::parsec,"parsec"},
         {"fermi",1,0,0,0,0,0,0,1.0e-15,"fermi"},

         // Area units
         // hectares, "ha" and "hectare"
         {"hectare",2,0,0,0,0,0,0,o2scl_mks::hectare,"hectare"},
         {"ha",2,0,0,0,0,0,0,1.0e4,"hectare"},
         // acre
         {"acre",2,0,0,0,0,0,0,o2scl_mks::acre,"acre"},
         // barn
         {"barn",2,0,0,0,0,0,0,o2scl_mks::barn,"barn"},
     
         // Volume units
         {"us_gallon",3,0,0,0,0,0,0,o2scl_mks::us_gallon,"gallon"},
         {"quart",3,0,0,0,0,0,0,o2scl_mks::quart,"quart"},
         {"pint",3,0,0,0,0,0,0,o2scl_mks::pint,"pint"},
         {"cup",3,0,0,0,0,0,0,o2scl_mks::cup,"cup"},
         {"tbsp",3,0,0,0,0,0,0,o2scl_mks::tablespoon,"tablespoon"},
         {"tsp",3,0,0,0,0,0,0,o2scl_mks::teaspoon,"teaspoon"},
         {"ca_gallon",3,0,0,0,0,0,0,o2scl_mks::canadian_gallon,
          "canadian gallon"},
         {"uk_gallon",3,0,0,0,0,0,0,o2scl_mks::uk_gallon,"uk gallon"},

         // Mass units
         // Solar masses, "Msun", and "Msolar"
         {"Msun",0,1,0,0,0,0,0,o2scl_mks::solar_mass,"solar mass"},
         {"Msolar",0,1,0,0,0,0,0,o2scl_mks::solar_mass,"solar mass"},
         {"pound",0,1,0,0,0,0,0,o2scl_mks::pound_mass,"pound"},
         {"ounce",0,1,0,0,0,0,0,o2scl_mks::ounce_mass,"ounce"},

         // The offical abbreviation for tonne is "t", but then SI
         // prefixes cause confusion between between "foot" and
         // "femtotonne". For now, we use "tonne" instead of "t".
         {"tonne",0,1,0,0,0,0,0,1.0e3,"(metric) tonne"},
         {"uk_ton",0,1,0,0,0,0,0,o2scl_mks::uk_ton,"uk ton"},
         {"troy_ounce",0,1,0,0,0,0,0,o2scl_mks::troy_ounce,"troy ounce"},
         {"carat",0,1,0,0,0,0,0,o2scl_mks::carat,"carat"},
     
         // Velocity units
         {"knot",1,0,-1,0,0,0,0,o2scl_mks::knot,"knot"},
         {"c",1,0,-1,0,0,0,0,o2scl_const::speed_of_light_f<fp_t>(),
          "speed of light"},
     
         // Energy units
         {"cal",2,1,-2,0,0,0,0,o2scl_mks::calorie,"calorie"},
         {"btu",2,1,-2,0,0,0,0,o2scl_mks::btu,"btu"},
         {"erg",2,1,-2,0,0,0,0,o2scl_mks::erg,"erg"},

         // Power units
         {"therm",2,1,-3,0,0,0,0,o2scl_mks::therm,"therm"},
         {"horsepower",2,1,-3,0,0,0,0,o2scl_mks::horsepower,"horsepower"},
         {"hp",2,1,-3,0,0,0,0,o2scl_mks::horsepower,"horsepower"},

         // Pressure units
         {"atm",-1,2,-2,0,0,0,0,o2scl_mks::std_atmosphere,"atmosphere"},
         {"bar",-1,1,-2,0,0,0,0,o2scl_mks::bar,"bar"},
         {"torr",-1,1,-2,0,0,0,0,o2scl_mks::torr,"torr"},
         {"psi",-1,1,-2,0,0,0,0,o2scl_mks::psi,"psi"},

         // Time units. Years are obtained from
         // https://pdg.lbl.gov/2021/reviews/contents_sports.html,
         // which references the 2020 Astronomical almanac.
         {"yr",0,0,1,0,0,0,0,31556925.1,"year (tropical)"},
         {"wk",0,0,1,0,0,0,0,o2scl_mks::week,"week"},
         {"d",0,0,1,0,0,0,0,o2scl_mks::day,"day"},

         // Angular units
         {"deg",0,0,0,0,0,0,0,o2scl_const::pi/180.0,"degree"},
         {"°",0,0,0,0,0,0,0,o2scl_const::pi/180.0,"degree"},
         {"′",0,0,0,0,0,0,0,o2scl_const::pi/10800.0,
          "arcminute (fraction of a degree)"},
         {"″",0,0,0,0,0,0,0,o2scl_const::pi/648000.0,
          "arcsecond (fraction of a degree)"},

         // Hours, "hr", to avoid confusion with Planck's constant
         {"hr",0,0,1,0,0,0,0,o2scl_mks::hour,"hour"},
         {"min",0,0,1,0,0,0,0,o2scl_mks::minute,"minute"},

         // Inverse time units
         {"curie",0,0,-1,0,0,0,0,o2scl_mks::curie,"curie"},
     
         // Force units
         {"dyne",1,1,-2,0,0,0,0,o2scl_mks::dyne,"dyne"},
     
         // Viscosity units
         {"poise",-1,1,-1,0,0,0,0,o2scl_mks::poise,"poise"},

         // Units of heat capacity or entropy
         // Boltzmann's constant. Could be confused with kilobytes?
         {"kB",2,1,-2,-1,0,0,0,o2scl_const::boltzmann_f<fp_t>(),
          "Boltzmann's constant"},
     
         {"hbar",2,1,-1,0,0,0,0,o2scl_const::hbar_f<fp_t>(),
          "reduced Planck constant"},
         // "Planck" instead of "h", to avoid confusing with hours
         {"Planck",2,1,-1,0,0,0,0,o2scl_const::planck_f<fp_t>(),
          "Planck constant"},
     
         // Gravitational constant. We cannot use "G" because of
         // confusion with "Gauss"
         {"GNewton",3,-1,-2,0,0,0,0,o2scl_mks::gravitational_constant,
          "gravitational constant"},
         
         // Gauss, and note the possible confusion with the gravitational
         // constant
         {"G",0,1,-2,0,-1,0,0,o2scl_mks::gauss,"gauss"},
         
         {"NA",0,0,0,0,0,-1,0,o2scl_const::avogadro,"Avogadro's number"}
     
        };
      
      other=other_;
      
      c_is_1=false;
      hbar_is_1=false;
      kb_is_1=false;
      G_is_1=false;
    }
    
    virtual ~convert_units() {}

    /// The associated find constants object
    find_constants fc;
    
    /** \brief Add a user-defined unit
     */
    void add_unit(const der_unit &d) {
      other.push_back(d);
      return;
    }
    
    /** \brief Remove a non-SI unit
     */
    void del_unit(std::string name) {
      size_t n_matches=0, i_match;
      for(size_t i=0;i<other.size();i++) {
        if (other[i].name==name) {
          n_matches++;
          i_match=i;
        }
      }
      if (n_matches==1) {
        typename std::vector<convert_units<fp_t>::der_unit>::iterator
          it=other.begin();
        it+=i_match;
        other.erase(it);
      }
      O2SCL_ERR2("Zero or more than one match found in ",
                "convert_units::del_unit().",o2scl::exc_efailed);
      return;
    }
    
    /** \brief Set natural units
     */
    void set_natural_units(bool c_is_one=true, bool hbar_is_one=true,
                           bool kb_is_one=true) {
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

        std::u32string c32;
        utf8_to_char32(SI[i].label,c32);
        int ns=c32.length();

        for(int j=0;j<15-ns;j++) out << ' ';
        
        out << SI[i].label << ' ';
        out.width(2);
        out << SI[i].m << ' ';
        out.width(2);
        out << SI[i].k << ' ';
        out.width(2);
        out << SI[i].s << ' ';
        out.width(2);
        out << SI[i].K << ' ';
        out.width(2);
        out << SI[i].A << ' ';
        out.width(3);
        out << SI[i].mol << ' ';
        out.width(2);
        out << SI[i].cd << ' ';
        out << SI[i].val << std::endl;
      }
      out << std::endl;

      out << "SI prefixes: prefix value" << std::endl;
      out << "------------------- ------------" << std::endl;
      for(size_t i=0;i<prefixes.size();i++) {
        
        std::u32string c32;
        utf8_to_char32(prefixes[i],c32);
        int ns=c32.length();

        for(int j=0;j<19-ns;j++) out << ' ';
        out << prefixes[i] << ' ';
        out << prefix_facts[i] << std::endl;
        
      }
      out << std::endl;
    
      out << "Other: label  m kg  s  K  A mol cd value        "
          << "full name" << std::endl;
      out << "------------ -- -- -- -- -- --- -- ------------ "
          << "---------------------" << std::endl;
          
    
      for(size_t i=0;i<other.size();i++) {
        
        std::u32string c32;
        utf8_to_char32(other[i].label,c32);
        int ns=c32.length();
        
        out << other[i].label;
        for(int j=0;j<13-ns;j++) out << ' ';
        
        out.width(2);
        out << other[i].m << ' ';
        out.width(2);
        out << other[i].k << ' ';
        out.width(2);
        out << other[i].s << ' ';
        out.width(2);
        out << other[i].K << ' ';
        out.width(2);
        out << other[i].A << ' ';
        out.width(3);
        out << other[i].mol << ' ';
        out.width(2);
        out << other[i].cd << ' ';
        out << other[i].val << ' ' << other[i].name << std::endl;
      }
    
      return;
    }

    /** \brief Print units to std::cout
     */
    void print_units_cout() {
      print_units(std::cout);
      return;
    }
    
    //void set_vars(fp_t m, fp_t k, fp_t s, fp_t K, fp_t A, fp_t mol,
    //fp_t cd, std::map<std::string, fp_t> &vars,
    //bool test_vars=false) const {
    
    void get_curr_unit_list(std::vector<std::string> &vs) const {

      vs.clear();
      
      // SI base units
      vs.push_back("m");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"m");
      }
      vs.push_back("m");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"m");
      }
      vs.push_back("g");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"g");
      }
      vs.push_back("s");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"s");
      }
      vs.push_back("K");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"K");
      }
      vs.push_back("A");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"A");
      }
      vs.push_back("mol");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"mol");
      }
      vs.push_back("cd");
      for(size_t i=0;i<n_prefixes;i++) {
        vs.push_back(prefixes[i]+"cd");
      }
      for(size_t i=0;i<SI.size();i++) {
        vs.push_back(SI[i].label);
        for(size_t j=0;j<n_prefixes;j++) {
          vs.push_back(prefixes[j]+SI[i].label);
        }
      }
      for(size_t i=0;i<other.size();i++) {
        vs.push_back(other[i].label);
      }
      return;
    }
    
    /** \brief Convert units, possibly using constants from
        the internal \ref o2scl::find_constants object
     */
    int convert_calc(std::string from, std::string to,
                      fp_t val, fp_t &converted,
                      fp_t &factor) {

      o2scl::calc_utf8 calc;
      o2scl::calc_utf8 calc2;

      int cret1=calc.compile_nothrow(from.c_str());
      if (cret1!=0) {
        if (verbose>0) {
          std::cout << "Compile from expression " << from << " failed."
                    << std::endl;
        }
        return 1;
      }
      int cret2=calc2.compile_nothrow(to.c_str());
      if (cret2!=0) {
        if (verbose>0) {
          std::cout << "Compile to expression " << to << " failed."
                    << std::endl;
        }
        return 2;
      }
      
      std::vector<std::u32string> vars=calc.get_var_list();
      std::vector<std::string> vars_str;
      for(size_t ij=0;ij<vars.size();ij++) {
        std::string tmp;
        char32_to_utf8(vars[ij],tmp);
        vars_str.push_back(tmp);
      }
      if (verbose>=2) {
        std::cout << "Unit list of from: ";
        o2scl::vector_out(std::cout,vars_str,true);
      }
      
      std::vector<std::u32string> vars2=calc2.get_var_list();
      std::vector<std::string> vars2_str;
      for(size_t ij=0;ij<vars2.size();ij++) {
        std::string tmp;
        char32_to_utf8(vars2[ij],tmp);
        vars2_str.push_back(tmp);
      }
      if (verbose>=2) {
        std::cout << "Unit list of to: ";
        o2scl::vector_out(std::cout,vars2_str,true);
      }

      std::vector<std::string> curr, new_units;
      get_curr_unit_list(curr);

      for(size_t i=0;i<vars_str.size();i++) {
        if (std::find(curr.begin(),curr.end(),vars_str[i])==curr.end() &&
            std::find(new_units.begin(),new_units.end(),
                      vars_str[i])==new_units.end()) {
          new_units.push_back(vars_str[i]);
        }
      }
      for(size_t i=0;i<vars2_str.size();i++) {
        if (std::find(curr.begin(),curr.end(),vars2_str[i])==curr.end() &&
            std::find(new_units.begin(),new_units.end(),
                      vars2_str[i])==new_units.end()) {
          new_units.push_back(vars2_str[i]);
        }
      }
      
      if (verbose>=2) {
        std::cout << "New units not found: ";
        o2scl::vector_out(std::cout,new_units,true);
      }
      
      std::vector<find_constants::const_entry> matches;
      for(size_t i=0;i<new_units.size();i++) {
        int fret=fc.find_nothrow(new_units[i],"mks",matches);
        if (fret==find_constants::one_exact_match_unit_match ||
            fret==find_constants::one_pattern_match_unit_match) {
          der_unit du;
          find_constants::const_entry &fcl=matches[0];
          du.label=new_units[i];
          du.m=fcl.m;
          du.k=fcl.k;
          du.s=fcl.s;
          du.K=fcl.K;
          du.A=fcl.A;
          du.mol=fcl.mol;
          du.cd=fcl.cd;
          du.val=fcl.val;
          du.name=fcl.names[0];
          if (verbose>=2) {
            std::cout << "Found constant " << new_units[i]
                      << " with value " << fcl.val << std::endl;
          }
          other.push_back(du);
        } else {
          if (verbose>=2) {
            std::cout << "New unit " << new_units[i]
                      << " not uniquely specified in constant list ("
                      << fret << ")." << std::endl;
          }
          return 1;
        }
      }

      return convert_calc_hck(from,to,val,converted,factor);
    }
    
    /** \brief Convert units, taking into account conversions
        which are allowed by setting hbar, c, or kB to 1
        
        This function is also useful as a const version of 
        convert_calc() .
    */
    int convert_calc_hck(std::string from, std::string to,
                         fp_t val, fp_t &converted,
                         fp_t &factor) const {
      
      if (verbose>=2) {
        std::cout << "Function convert_units::convert_calc(), "
                  << "kb_is_1,hbar_is_1,c_is_1: " << kb_is_1 << ' '
                  << hbar_is_1 << ' ' << c_is_1 << std::endl;
      }
      
      // These calculator objects have to be inside this
      // function to make the function const
      o2scl::calc_utf8 calc;
      o2scl::calc_utf8 calc2;

      std::map<std::string, fp_t> vars;

      set_vars(1.0,1.0,1.0,1.0,1.0,1.0,1.0,vars);

      if (verbose>=3) {
        std::cout << "Compile: " << from << std::endl;
      }
      int cret1=calc.compile_nothrow(from.c_str());
      if (verbose>=3) {
        std::cout << "Result: " << cret1 << std::endl;
        std::cout << "Compile: " << to << std::endl;
      }
      if (cret1!=0) return 1;
      int cret2=calc2.compile_nothrow(to.c_str());
      if (verbose>=3) {
        std::cout << "Result: " << cret2 << std::endl;
      }
      if (cret2!=0) return 2;
      fp_t before, after;
      int cret3=calc.eval_nothrow(&vars,before);
      if (verbose>=3) {
        std::cout << "Result: " << cret3 << std::endl;
        std::cout << "eval: " << before << std::endl;
      }
      if (cret3!=0) return 3;
      int cret4=calc2.eval_nothrow(&vars,after);
      if (verbose>=3) {
        std::cout << "Result: " << cret4 << std::endl;
        std::cout << "eval: " << after << std::endl;
      }
      if (cret4!=0) {
        return 4;
      }
      
      factor=before/after;
      converted=val*factor;
        
      // Now, having verified that a conversion is possible, we see
      // if any additional factors of k_B, hbar, or c are required
      // to perform the conversion.

      fp_t kb_power=0, factor_kb=1, hbar_power=0, factor_hbar=1;
      fp_t c_power=0, factor_c=1, two=2;
    
      if (kb_is_1) {

        // Determine how many powers of kB are required to match
        // powers of K
        set_vars(1.0,1.0,1.0,2.0,1.0,1.0,1.0,vars);
        kb_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
        if (verbose>1) {
          std::cout << "kb: " << before << ' ' << after << ' '
                    << calc.eval(&vars) << "\n\t"
                    << calc2.eval(&vars) << ' ' << kb_power << ' ';
        }
        kb_power=log(kb_power)/log(two);
        if (verbose>1) {
          std::cout << kb_power << ' ';
        }
        factor_kb=pow(o2scl_mks::boltzmann,kb_power);
        if (verbose>1) {
          std::cout << factor_kb << std::endl;
        }
        
      }

      if (hbar_is_1) {
        
        // Determine how many powers of hbar are required to match
        // powers of kg
        set_vars(1.0,2.0,1.0,1.0,1.0,1.0,1.0,vars);
        hbar_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
        if (verbose>1) {
          std::cout << "hbar: " << before << ' ' << after << ' '
                    << calc.eval(&vars) << "\n\t"
                    << calc2.eval(&vars) << ' ' << hbar_power << ' ';
        }
        hbar_power=-log(hbar_power)/log(two)-kb_power;
        if (verbose>1) {
          std::cout << hbar_power << ' ';
        }
        factor_hbar=pow(o2scl_mks::plancks_constant_hbar,hbar_power);
        if (verbose>1) {
          std::cout << factor_hbar << std::endl;
        }
        
      }
        
      if (c_is_1) {
        
        // Determine how many powers of c are required to match
        // powers of s
        set_vars(1.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
        c_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
        if (verbose>1) {
          std::cout << "c: " << before << ' ' << after << ' '
                    << calc.eval(&vars) << "\n\t"
                    << calc2.eval(&vars) << ' ' << c_power << ' ';
        }
        c_power=log(c_power)/log(two)-2*kb_power-hbar_power;
        if (verbose>1) {
          std::cout << c_power << ' ';
        }
        factor_c=pow(o2scl_mks::speed_of_light,c_power);
        if (verbose>1) {
          std::cout << factor_c << std::endl;
        }
        
      }

      // Determine the number of powers of length remaining
      // and make sure it's consistent
      set_vars(2.0,1.0,1.0,1.0,1.0,1.0,1.0,vars);
      fp_t m_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
      m_power=log(m_power)/log(two);
      if (verbose>1) {
        std::cout << "m_power: " << m_power << ' '
                  << 2*kb_power+2*hbar_power+c_power << ' '
                  << fabs(m_power+2*kb_power+2*hbar_power+c_power)
                  << std::endl;
      }

      // If this is non-zero, the units are inconsistent
      if (fabs(m_power+2*kb_power+2*hbar_power+c_power)>1.0e-12) {
        return 5;
      }

      // Determine the number of powers of mass remaining
      // and make sure it's consistent
      set_vars(1.0,2.0,1.0,1.0,1.0,1.0,1.0,vars);
      fp_t kg_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
      kg_power=log(kg_power)/log(two);
      if (verbose>1) {
        std::cout << "kg_power: " << kg_power << ' '
                  << kb_power+hbar_power << std::endl;
      }
          
      // If this is non-zero, the units are inconsistent
      if (fabs(kg_power+kb_power+hbar_power)>1.0e-12) {
        return 6;
      }

      // Determine the number of powers of time remaining
      // and make sure it's consistent
      set_vars(1.0,1.0,2.0,1.0,1.0,1.0,1.0,vars);
      fp_t s_power=calc.eval(&vars)/calc2.eval(&vars)/before*after;
      s_power=log(s_power)/log(two);
      if (verbose>1) {
        std::cout << "s_power: " << s_power << ' '
                  << -2*kb_power-hbar_power-c_power << std::endl;
      }
      
      // If this is non-zero, the units are inconsistent
      if (fabs(s_power-2*kb_power-hbar_power-c_power)>1.0e-12) {
        return 7;
      }
      
      factor*=factor_kb*factor_hbar*factor_c;
      converted*=factor_kb*factor_hbar*factor_c;

      /*
        set_vars(1.0,1.0,1.0,1.0,2.0,1.0,1.0,vars);
        fp_t factor_A=calc.eval(&vars)/calc2.eval(&vars);
        
        set_vars(1.0,1.0,1.0,1.0,1.0,2.0,1.0,vars);
        fp_t factor_mol=calc.eval(&vars)/calc2.eval(&vars);
        
        set_vars(1.0,1.0,1.0,1.0,1.0,1.0,2.0,vars);
        fp_t factor_cd=calc.eval(&vars)/calc2.eval(&vars);
      */

      if (verbose>1) {
        std::cout << "from: " << from << " to: " << to
                  << " before: " << before << " after: " << after 
                  << " factor: " << factor << std::endl;
        std::cout << "  factor_kb: " << factor_kb
                  << " factor_hbar: " << factor_hbar
                  << " factor_c: " << factor_c << std::endl; 
        //<< " factor_K: " << factor_K << "\n\tfactor_A: " << factor_A 
        //<< " factor_mol: " << factor_mol
        //<< " factor_cd: " << factor_cd << std::endl;
      }    

      return 0;
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
        std::cout << "----------------------------------"
                  << "----------------------------------" << std::endl;
        for (m=mcache.begin();m!=mcache.end();m++) {
          std::cout.setf(std::ios::left);
          std::cout.width(25);
          std::cout << m->second.f << ' ';
          std::cout.width(25);
          std::cout << m->second.t << ' ';
          std::cout.unsetf(std::ios::left);
          std::cout.precision(10);
          std::cout << m->second.c << std::endl;
        }
      }
      return;
    }

    /** \brief Test the cache and the new convert_calc() function

        \note We make this a template to avoid problems with
        circular headers.
     */
    template<class test_mgr_t> void test_cache_calc(test_mgr_t &t) const {
      mciter m;
      for (m=mcache.begin();m!=mcache.end();m++) {
        std::cout.setf(std::ios::left);
        std::cout.width(10);
        std::cout << m->second.f << ' ';
        std::cout.width(10);
        std::cout << m->second.t << ' ';
        std::cout.unsetf(std::ios::left);
        std::cout.precision(6);
        std::cout << m->second.c << ' ';
        fp_t c, f;
        int ix=convert_calc(m->second.f,m->second.t,1.0,c,f);
        std::cout << ix << ' ' << f << ' ';
        t.test_rel(f,m->second.c,1.0e-13,"test_cache_calc");
        std::cout << (fabs(f-m->second.c)/fabs(f)<1.0e-12) << std::endl;
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

      set_natural_units(true,true,true);
      
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

      // (tropical year)
      insert_cache("yr","s",o2scl_mks::tropical_year);
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
      // This quantity is constructed from the Julian year
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
    //@}
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
