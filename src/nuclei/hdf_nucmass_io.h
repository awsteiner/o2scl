/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2014, Andrew W. Steiner

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
#ifndef O2SCL_HDF_NUCMASS_IO_H
#define O2SCL_HDF_NUCMASS_IO_H

/** \file hdf_nucmass_io.h
    \brief File for HDF input of the \o2 \ref o2scl::nucmass_ame and 
    \ref o2scl::mnmsk_mass data files
*/

#include <hdf5.h>
#include <hdf5_hl.h>

#include <o2scl/constants.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_frdm.h>

#ifndef DOXYGEN_NO_O2NS

// While ideally we would leave this o2scl_hdf namespace in the
// documentation, Doxygen seems to get confused with the o2scl_hdf
// documentation already in the parent O2scl library. 
// For now, we leave the o2scl_hdf namespace in

namespace o2scl_hdf {
#endif

  /** \brief Read data for \ref o2scl::nucmass_ame from an HDF table

      There are five data sets, selected by the specification of the
      \c version string
      - "95rmd" - "Recommended" data from \ref Audi95 (ame95rmd.o2)
      - "95exp" - "Experimental" data from \ref Audi95 (ame95exp.o2)
      - "03round" - "Rounded" data from \ref Audi03 (ame03round.o2)
      - "03" - Data from \ref Audi03 (ame03.o2)
      - "12" - Data from \ref Audi12 (default) (ame12.o2)
      
      If any string other than these five is used, the default data is
      loaded. If the constructor cannot find the data file (e.g.
      because of a broken installation), then nucmass_ame::is_loaded()
      returns false.
      
      \note This function is in the \ref o2scl_hdf namespace.

  */
  void ame_load(o2scl::nucmass_ame &ame, std::string file_name, 
		std::string table_name);

  /** \brief Read an AME mass table from an external file
   */
  void ame_load(o2scl::nucmass_ame &ame, std::string name="");
  
  /** \brief Read data for \ref o2scl::mnmsk_mass from an HDF table

      \note This function is in the \ref o2scl_hdf namespace.

      \todo Switch to the same 'external' file mechanism used in 
      hdf_eos_io.h. 
   */
  void mnmsk_load(o2scl::mnmsk_mass &mnmsk, std::string dir="");
  
  /** \brief Read data for \ref o2scl::nucmass_hfb from an HDF table
      
      Valid values of \c model at present are 2, 8, and 14, corresponding
      to the HFB2 (\ref Goriely02), HFB8 (\ref Samyn04), and HFB14 
      (\ref Goriely07). If a number other than these three is given,
      the error handler is called. 

      \note This function is in the \ref o2scl_hdf namespace.

      \todo Switch to the same 'external' file mechanism used in 
      hdf_eos_io.h. 
  */
  void hfb_load(o2scl::nucmass_hfb &hfb, size_t model=14, std::string dir="");

  /** \brief Read data for \ref o2scl::nucmass_hfb from an HDF table
      
      Valid values of \c model at present are 17 and 21, corresponding
      to the HFB17 (\ref Goriely02) and HFB21 (\ref Samyn04). If a
      number other than these two is given, the error handler is
      called.

      \note This function is in the \ref o2scl_hdf namespace.

      \todo Switch to the same 'external' file mechanism used in 
      hdf_eos_io.h. 
  */
  void hfb_sp_load(o2scl::hfb_sp_mass &hfb, size_t model=21, 
		   std::string dir="");

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
