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
#ifndef O2SCL_HDF_EOS_IO_H
#define O2SCL_HDF_EOS_IO_H

/** \file hdf_eos_io.h
    \brief HDF input of the \o2 EOS data files
*/

#include <hdf5.h>

#include <o2scl/constants.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>
#include <o2scl/skyrme_eos.h>
#include <o2scl/rmf_eos.h>
#include <o2scl/gogny_eos.h>

/** \brief Additional functions to read and write EOS data to HDF5 files
 */
namespace o2scl_hdf {

  /** \brief Read the Gogny EOS from a data file
   */
  void gogny_load(o2scl::gogny_eos &ge, std::string model, 
		  bool external=false);
  
  /** \brief Input a \ref o2scl::rmf_eos object from an HDF file

      If \ref external is <tt>false</tt> (the default), then the model
      is loaded from the \o2 data directory <tt>rmfdata</tt> with the
      suffix <tt>.o2</tt>. Otherwise, the parameter \c model is 
      taken to be the full pathname of the HDF5 file containing 
      the model to be loaded.
  */
  void rmf_load(o2scl::rmf_eos &rmf, std::string model, bool external=false);
  
  /** \brief Input a \ref o2scl::skyrme_eos object from an HDF file

      If \ref external is <tt>false</tt> (the default), then the model
      is loaded from the \o2 data directory <tt>skdata</tt> with the
      suffix <tt>.o2</tt>. Otherwise, the parameter \c model is 
      taken to be the full pathname of the HDF5 file containing 
      the model to be loaded.
  */
  void skyrme_load(o2scl::skyrme_eos &sk, std::string model, 
		   bool external=false);
  
  /** \brief Write a \ref o2scl::skyrme_eos object to an HDF file
   */
  void skyrme_write(hdf_file &hf, o2scl::skyrme_eos &sk, std::string name);
  
  /** \brief Write a \ref o2scl::skyrme_eos object to an HDF file in the \o2
      data directory
  */
  void skyrme_write(o2scl::skyrme_eos &sk, std::string model);
  
}

#endif


