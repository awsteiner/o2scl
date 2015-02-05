/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_SHARED_PTR_H
#define O2SCL_SHARED_PTR_H

/** \file shared_ptr.h
    \brief File defining \ref o2scl::o2_shared_ptr
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#ifdef DOXYGEN

  /** \brief A struct to provide the shared_ptr type

      This object exists in order to provide the shared_ptr template
      type used in \o2. The full specification of a shared pointer
      in \o2 for an object of type \c T is thus 
      \verbatim
      o2scl::o2_shared_ptr<T>::type
      \endverbatim
      In a default \o2 installation, \ref type (as given below)
      is a typedef of
      \verbatim
      std::tr1::shared_ptr<T>
      \endverbatim
      Ifs <tt>O2SCL_HAVE_BOOST</tt>
      is defined, then it is a typedef of 
      \verbatim
      boost::shared_ptr<T>
      \endverbatim

      See also the discussion at http://www.gotw.ca/gotw/079.htm . This
      struct won't be necessary when C++ allows template typedef's as
      part of the C++11 standard http://en.wikipedia.org/wiki/C%2B%2B11
      , but very few compilers have implemented this standard yet.
  */
  template<class T> struct o2_shared_ptr {
    /// The actual shared_ptr type
    typedef std::tr1::shared_ptr<T> type;
  };

#endif

#ifndef DOXYGEN_NO_O2NS
}
#endif

// -------------------------------------------------------------------
// Define the o2_shared_ptr struct according to the installation settings

// AWS - 11/29/11: I can't remember if the #include statements have
// to be outside the o2scl namespace, but I make sure they're
// outside just in case it matters.

#ifndef O2SCL_NO_TR1_MEMORY

#include <tr1/memory>
namespace o2scl {
  template<class T> struct o2_shared_ptr {
    typedef std::tr1::shared_ptr<T> type;
  };
}

#else

#include <boost/shared_ptr.hpp>
namespace o2scl {
  template<class T> struct o2_shared_ptr {
    typedef boost::shared_ptr<T> type;
  };
}

// end of if O2SCL_NO_TR1_MEMORY
#endif

// -------------------------------------------------------------------

// end of ifdef O2SCL_SHARED_PTR_H
#endif 
