/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2013, Andrew W. Steiner
  
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
    \brief File for \ref o2scl::convert_units
*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>

#include <o2scl/err_hnd.h>
#include <o2scl/misc.h>
#include <o2scl/constants.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Convert units

      Allow the user to convert between two different units after
      specifying a conversion factor. This class will also
      automatically combine two conversion factors to create a new
      unit conversion (but it cannot combine more than two).

      Conversions are performed by the \ref convert() function and the
      conversion factors must be specified beforehand using the \ref
      insert_cache() function.

      If the GNU units command is not in the local path, the user may
      modify \ref units_cmd_string to specify the full pathname. One
      can also modify \ref units_cmd_string to specify a different
      <tt>units.dat</tt> file.

      \future A remove_cache() and in_cache() function to test
      to see if a conversion is currently in the cache. 
      
      Example:
      \code
      convert_units cu;
      cu.insert_cache("in","cm",2.54);
      cout << "12 in is " << cu.convert("in","cm",12.0) << " cm. " << endl;
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

      \future Ideally, a real C++ API for the GNU units command
      would be better.
  */
  class convert_units {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The type for caching unit conversions
    typedef struct {
      /// The input unit
      std::string f;
      /// The output unit
      std::string t;
      /// The conversion factor
      double c;
    } unit_t;

    /// The cache where unit conversions are stored
    std::map<std::string,unit_t,string_comp> mcache;
    
    /// The iterator type
    typedef std::map<std::string,unit_t,string_comp>::iterator miter;
      
#endif

  public:

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

    convert_units();

    virtual ~convert_units() {}

    /** \brief Return the value \c val after converting using units \c
	from and \c to
    */
    virtual double convert(std::string from, std::string to, double val);

    /// Manually insert a unit conversion into the cache
    int insert_cache(std::string from, std::string to, double conv);

    /// Manually remove a unit conversion into the cache
    int remove_cache(std::string from, std::string to);
    
    /// Print the present unit cache to std::cout
    int print_cache();

    /** \brief Make a GNU \c units.dat file from the GSL constants

	If \c c_1 is true, then the second is defined in terms of
	meters so that the speed of light is unitless. If \c hbar_1 is
	true, then the kilogram is defined in terms of <tt>s/m^2</tt>
	so that \f$ \hbar \f$ is unitless.

	\note Not all of the GSL constants or the canonical GNU units 
	conversions are given here.
    */
    int make_units_dat(std::string fname, bool c_1=false, 
		       bool hbar_1=false, bool K_1=false);
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
