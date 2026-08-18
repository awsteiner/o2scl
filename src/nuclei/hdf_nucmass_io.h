/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2026, Andrew W. Steiner

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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_HDF_NUCMASS_IO_H
#define O2SCL_HDF_NUCMASS_IO_H

/** \file hdf_nucmass_io.h
    \brief File for HDF input of the \ref o2scl::nucmass_ame and 
    \ref o2scl::nucmass_mnmsk data files
*/

#ifdef O2SCL_PLAIN_HDF5_HEADER
#include <hdf5.h>
#include <hdf5_hl.h>
#else
#ifdef O2SCL_LINUX
#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/hdf5_hl.h>
#else
#include <hdf5.h>
#include <hdf5_hl.h>
#endif
#endif

#include <o2scl/constants.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_frdm.h>

namespace o2scl_hdf {

#ifdef O2SCL_NEVER_DEFINED
  
  /** \brief Read data for \ref o2scl::nucmass_ame from an HDF table
      specified in a file
      
      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .
  */
  void ame_load_ext(o2scl::nucmass_ame &ame, std::string file_name, 
    std::string table_name, bool exp_only=false);

  /** \brief Read data for \ref o2scl::nucmass_ame from an HDF table
      specified in a file
      
      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .
  */
  void ame_load(o2scl::nucmass_ame &ame, std::string name="20",
                bool exp_only=false);

#endif
  
  /** \brief Read data for \ref o2scl::nucmass_mnmsk from an HDF table

      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .

      If the filename is unspecified, then the default O2scl
      data file is loaded.
   */
  void mnmsk_load(o2scl::nucmass_mnmsk &mnmsk, std::string model="",
                  std::string filename="");
  
  /** \brief Read data for \ref o2scl::nucmass_hfb from an HDF table
      
      Valid values of \c model at present are 2, 8, and 14, corresponding
      to the HFB2 (Goriely02), HFB8 (Samyn04), and HFB14 
      (Goriely07). If a number other than these three is given,
      the error handler is called. 

      \verbatim embed:rst
      See also [Goriely02]_, [Samyn04]_, and [Goriely07]_.
      \endverbatim

      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .
  */
  void hfb_load(o2scl::nucmass_hfb &hfb, size_t model=14,
                std::string filename="");

  /** \brief Read data for \ref o2scl::nucmass_hfb from an HDF table
      
      Valid values of \c model at present are 17, and 21 through 32.
      The first two correspond to the HFB17 (Goriely02) and HFB21
      (Samyn04). If a number outside this range is given, then 32
      is assumed.

      \verbatim embed:rst
      See also [Goriely02]_, [Samyn04]_, and [Goriely07]_.

      .. todo::

         In hfb_sp_load(): Document models 22 through 32.

      \endverbatim

      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .
  */
  void hfb_sp_load(o2scl::nucmass_hfb_sp &hfb, size_t model=27,
                   std::string filename="");

  /** \brief Read data for \ref o2scl::nucmass_hfb_sp from an HDF
      table of Brussels-Skyrme-on-a-Grid (BSkG) masses

      Valid values of \c model are 1, 2, 3, and 4, corresponding to
      BSkG1 [Scamps21]_, BSkG2 [Ryssens22]_, BSkG3 [Grams23]_, and
      BSkG4 [Grams24]_, respectively. If a number outside this
      range is given, the error handler is called.

      The BSkG tables fill in various extra fields in \ref
      o2scl::nucmass_hfb_sp::entry which are left unused by \ref
      hfb_sp_load() . BSkG3 fills in the extra deformation,
      separation-energy, odd-even staggering, and moment of
      inertia fields (\ref o2scl::nucmass_hfb_sp::entry::gamma
      through \ref o2scl::nucmass_hfb_sp::entry::I3). BSkG1, BSkG2,
      and BSkG4 instead fill in the rotational-correction,
      pairing-gap, experimental-radius, moment-of-inertia, and
      subsystem-parity fields (\ref
      o2scl::nucmass_hfb_sp::entry::Erot through \ref
      o2scl::nucmass_hfb_sp::entry::par_n), and BSkG1 and BSkG2 do
      not provide \ref o2scl::nucmass_hfb_sp::entry::beta30 or \ref
      o2scl::nucmass_hfb_sp::entry::beta32 . None of the BSkG
      tables provide a value for \ref
      o2scl::nucmass_hfb_sp::entry::def_wig, which is thus left at
      zero.

      \verbatim embed:rst
      See also [Scamps21]_, [Ryssens22]_, [Grams23]_, and
      [Grams24]_.
      \endverbatim

      \note This function is in the o2scl_hdf namespace,
      see \ref hdf_nucmass_io.h .
  */
  void bskg_load(o2scl::nucmass_hfb_sp &hfb, size_t model=3,
                 std::string filename="");

}

#endif
