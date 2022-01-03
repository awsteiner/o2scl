/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2022, Andrew W. Steiner

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

#ifdef O2SCL_PLAIN_HDF5_HEADER
#include <hdf5.h>
#else
#ifdef O2SCL_LINUX
#include <hdf5/serial/hdf5.h>
#else
#include <hdf5.h>
#endif
#endif

#include <o2scl/constants.h>
#include <o2scl/hdf_file.h>
#include <o2scl/lib_settings.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/eos_had_gogny.h>

/** \brief Additional functions to read and write EOS data to HDF5 files
 */
namespace o2scl_hdf {

  /** \brief Read the Gogny EOS from a data file

      If \c external is <tt>false</tt> (the default), then the model
      (either <tt>"d1n"</tt> or <tt>"d1s"</tt> is loaded from the \o2
      data directory in file <tt>gogny.o2</tt>. Otherwise, the
      parameter \c model is taken to be the full pathname of the HDF5
      file containing the EOS model data to be loaded.
  */
  void gogny_load(o2scl::eos_had_gogny &ge, std::string model,
		  std::string filename="");
  
  /** \brief Input a \ref o2scl::eos_had_rmf object from an HDF file

      If \c external is <tt>false</tt> (the default), then the model
      is loaded from the \o2 data directory <tt>rmfdata</tt> with the
      suffix <tt>.o2</tt>. Otherwise, the parameter \c model is 
      taken to be the full pathname of the HDF5 file containing 
      the EOS model data to be loaded.
  */
  void rmf_load(o2scl::eos_had_rmf &rmf, std::string model, 
		bool external=false);
  
  /** \brief Input a \ref o2scl::eos_had_skyrme object from an HDF file

      If \c external is <tt>false</tt> (the default), then the model
      is loaded from the \o2 data directory <tt>skdata</tt> with the
      suffix <tt>.o2</tt>. Otherwise, the parameter \c model is 
      taken to be the full pathname of the HDF5 file containing 
      the EOS model data to be loaded.

      The parameters <tt>b4</tt>, <tt>b4p</tt> and the reference
      are directly read from the file. 
      
      If the file does not contain an integer named <tt>dpfix</tt>,
      or the value of <tt>dpfix</tt> is false, then 
      - the parameters named <tt>x0</tt>, <tt>x1</tt>, <tt>x2</tt>, 
      <tt>x3</tt>, <tt>a</tt>, <tt>b</tt>, and <tt>alpha</tt>
      are directly read from the file,
      - and the parameters named <tt>t0</tt>, <tt>t1</tt>, <tt>t2</tt>, 
      and <tt>t3</tt> are presumed to be stored in the file with
      an extra factor of \f$ \hbar c \f$ and stored in fields named
      <tt>t0hc</tt>, <tt>t1hc</tt>, <tt>t2hc</tt>, and
      <tt>t3hc</tt> . 
      
      Alternatively if <tt>dpfix</tt> is present and greater than
      zero, then the values \f$ t_1=-t_2/3(5+4 x_2) \f$, \f$ x_1 =
      -(4+5 x_2)/ (5+4 x_2) \f$, \f$ \alpha=1/3 \f$, \f$ a=1 \f$ and
      \f$ b=0 \f$ are assumed. The values <tt>x0</tt>, <tt>x2</tt>,
      and <tt>x3</tt> are directly read and the values <tt>t0</tt>,
      <tt>t2</tt>, and <tt>t3</tt> are computed from fields named
      <tt>t0hc</tt>, <tt>t2hc</tt>, and <tt>t3hc</tt> .

      If the file contains an integer named <tt>pdmode</tt>
      and that integer is greater than zero, then 
      the parameter named <tt>W0</tt> is taken from
      the numbers named <tt>pairfn</tt> and <tt>pairfp</tt>
      using the relation
      \f[
      W_0 = \frac{(\mathrm{pairfn}+\mathrm{pairfp})}{4 \hbar c}
      \f]
      Otherwise, it is assumed that the file contains a 
      field named <tt>W0hc</tt> which stores the value of
      <tt>W0</tt> times \f$ \hbar c \f$ .
  */
  void skyrme_load(o2scl::eos_had_skyrme &sk, std::string model, 
		   bool external=false, int verbose=0);
  
  /** \brief Write a \ref o2scl::eos_had_skyrme object to an HDF file
   */
  void skyrme_write(hdf_file &hf, o2scl::eos_had_skyrme &sk,
		    std::string name);
  
  /** \brief Write a \ref o2scl::eos_had_skyrme object to an HDF file
      in the \o2 data directory
  */
  void skyrme_write(o2scl::eos_had_skyrme &sk, std::string model);

  /** \brief Return a pointer to an eos_had_base object 
      from two strings specifying type and name
  */
  o2scl::eos_had_base *eos_had_strings(std::string type,
				       std::string name="");

  /** \brief List EOSs understood by \ref eos_had_strings() .
  */
  void eos_had_strings_list();

}

#endif
