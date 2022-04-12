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
/** \file eos_sn.h
    \brief File defining \ref o2scl::eos_sn_base
*/
#ifndef O2SCL_EOS_SN_H
#define O2SCL_EOS_SN_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <o2scl/constants.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/table.h>
#include <o2scl/boson_eff.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/test_mgr.h>
#include <o2scl/convert_units.h>
#include <o2scl/interp2_direct.h>
#include <o2scl/eos_base.h>

namespace o2scl {

  /** \brief A base class for the supernova EOSs [abstract]

      This class is experimental.

      \verbatim embed:rst
      See also the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.
      \endverbatim

      \comment 
      \todo Allow logarithmic grids for any of nb, Ye, or T. 
      12/10/13: The member variables are in the parent class, but no 
      code is written to use them yet. What really are these for?
      1/3/14: In fact, the linear vs. logarithmic distinction isn't
      necessarily useful, because some of the grids (e.g. T for sht)
      aren't either purely linear or purely logarithmic.
      \endcomment

      \future Add option to rescale energies and chemical 
      potentials to different masses.
      \future Create a \ref o2scl::table object, possibly using 
      tensor_grid::vector_slice. 
      \future Show how matrix_slice and vector_slice can be used
      with this object (Is this done in ex_eos_sn or is this
      made obsolete by the tensor rearrange function?)
      \future Add option to load and store a separate lepton/photon
      EOS

      \comment
      \future Could this be a child of eos_had_temp_base and
      then directly used in nstar_cold()? Actually no, this doesn't
      work because of the nuclei. 
      \endcomment

  */
  class eos_sn_base {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    eos_sn_base();
    
    virtual ~eos_sn_base();

    /// \name Grid and data sizes
    //@{
    /// Size of baryon density grid
    size_t n_nB;
    /// Size of electron fraction grid
    size_t n_Ye;
    /// Size of temperature grid
    size_t n_T;
    /// Baryon density grid (in \f$ \mathrm{fm}^{-3} \f$)
    std::vector<double> nB_grid;
    /// Electron fraction grid
    std::vector<double> Ye_grid;
    /// Temperature grid (in \f$ \mathrm{MeV} \f$)
    std::vector<double> T_grid;
    /// Number of additional data sets
    size_t n_oth;
    /// Number of base data sets
    static const size_t n_base=16;
    //@}

    /// EOS of leptons and photons
    eos_leptons elep;
    
    /// \name Data
    //@{
    /** \brief Total free energy per baryon in MeV (without 
        baryon rest masses but including electron rest mass)

        By default, this energy is relative to 
        \f[
        m_n (1-Y_e) + m_p Y_e
        \f]
        where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
        is stored in \ref m_prot .
    */
    tensor_grid3<> F;
    /** \brief Free energy per baryon without lepton and photon 
        contributions in MeV

        By default, this energy is relative to 
        \f[
        m_n (1-Y_e) + m_p Y_e
        \f]
        where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
        is stored in \ref m_prot .
    */
    tensor_grid3<> Fint;
    /** \brief Total internal energy per baryon in MeV (without 
        baryon rest masses but including electron rest mass)

        By default, this energy is relative to 
        \f[
        m_n (1-Y_e) + m_p Y_e
        \f]
        where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
        is stored in \ref m_prot .
    */
    tensor_grid3<> E;
    /** \brief Internal energy per baryon without lepton and photon 
        contributions in MeV

        By default, this energy is relative to 
        \f[
        m_n (1-Y_e) + m_p Y_e
        \f]
        where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
        is stored in \ref m_prot .
    */
    tensor_grid3<> Eint;
    /// Total pressure in \f$ \mathrm{MeV}/\mathrm{fm}^3 \f$
    tensor_grid3<> P;
    /** \brief Pressure without lepton and photon contributions
        in \f$ \mathrm{MeV}/\mathrm{fm}^3 \f$
    */
    tensor_grid3<> Pint;
    /// Total entropy per baryon
    tensor_grid3<> S;
    /// Entry per baryon without lepton and photon contributions
    tensor_grid3<> Sint;
    /** \brief Neutron chemical potential in MeV

        By default this is relative to the neutron mass in
        \ref m_neut .
    */
    tensor_grid3<> mun;
    /** \brief Proton chemical potential in MeV

        By default this is relative to the proton mass in
        \ref m_prot .
    */
    tensor_grid3<> mup;
    /// Proton number
    tensor_grid3<> Z;
    /// Mass number
    tensor_grid3<> A;
    /// Neutron baryon fraction
    tensor_grid3<> Xn;
    /// Proton baryon fraction
    tensor_grid3<> Xp;
    /// Alpha particle baryon fraction
    tensor_grid3<> Xalpha;
    /// Heavy nuclei baryon fraction
    tensor_grid3<> Xnuclei;
    /// Other data sets
    tensor_grid3<> other[30];
    /// List of pointers to data
    tensor_grid3<> *arr[n_base+30];
    //@}

    /** \brief Check the table composition entries
     */
    void check_composition(double &max1, double &max2);

    /// \name Interpolation
    //@{
    /** \brief Set the interpolation type of all the 
        \ref o2scl::tensor_grid3 objects to type \c interp_type .

        \note This is used by the constructor to set all tensors
        to linear interpolation.
    */
    void set_interp_type(size_t interp_type);
    //@}

    /// \name Nucleon masses
    //@{
    /** \brief Neutron mass in \f$ \mathrm{MeV} \f$ 
        (defaults to o2scl_mks::mass_neutron times o2scl_const::hc_mev_fm)
    */
    double m_neut;
    
    /** \brief Proton mass in \f$ \mathrm{MeV} \f$ 
        (defaults to o2scl_mks::mass_proton times o2scl_const::hc_mev_fm)
    */
    double m_prot;
    //@}

    /// \name Electron and photon contribution
    //@{
    /// Photon
    boson photon;
    /// Electron
    fermion electron;
    /// Muon
    fermion muon;
    /// If true, include muons
    bool include_muons;
    /// Relativistic fermion thermodynamics
    fermion_rel relf;
    /** \brief Compute the electron and photon contribution for the full
        grid

        If \ref baryons_only is true, this function computes
        the data for <tt>E, P, S,</tt> and <tt>F</tt> by adding
        electrons and photons to the baryon contributions stored in
        <tt>Eint, Pint, Sint,</tt> and <tt>Fint</tt>. Otherwise,
        this function computes <tt>Eint, Pint, Sint,</tt> and 
        <tt>Fint</tt> by subtracting electron and photon 
        contributions from <tt>E, P, S,</tt> and <tt>F</tt>. 

        The electron contribution to the internal energy and free
        energy computed by this function includes the electron rest
        mass.
    */
    virtual void compute_eg();

    /** \brief Compute lepton contribution at one point

        The temperature is to be specified in \c MeV. An initial guess
        for the electron chemical potential, \c mue, should be stored
        in units of \f$ \mathrm{fm}^{-1} \f$. This value should also
        include the electron rest mass. After completion, the final
        value of the electron chemical potential will be stored in \c
        mue (including the rest mass).

        This function will fail if a sufficiently poor initial
        guess is given.
    */
    virtual void compute_eg_point(double nB, double Ye, double TMeV,
                                  thermo &th, double &mue);

    /** \brief Check electrons and photons
        
        This checks that the electron and photon thermodynamics
        generated by \o2 is consistent with the data in 
        \c E, \c Eint, \c F, \c Fint, \c P, \c Pint, \c S,
        and \c Sint.
    */
    virtual double check_eg();
    //@}

    /** \brief Test the free energy and store results in \c tm

        This checks that the data in \c Fint is consistent with that
        in \c Eint and \c Sint (if \ref baryons_only is true)
        and that \c F is consistent with that in \c E and \c S (if
        \ref with_leptons is true).
    */
    void check_free_energy(double &avg);
    
    /** \brief Verbosity parameter (default 1)
     */
    int verbose;

    /** \brief Compute properties of matter in beta equilibrium
        at fixed entropy per baryon

        This function just does a simple hard-coded linear
        interpolation.

        The temperature is returned in units of MeV.
    */
    virtual void beta_eq_sfixed
    (double nB, double entr, double &Ye, double &T);

    /** \brief Compute the electron fraction for beta-equilibrium
        at fixed density and temperature 
        
        This function just uses linear interpolation to 
        interpolate in baryon density and temperature and 
        the uses a quadratic to determine the minimum of the
        free energy.

        If \ref data_with_leptons() is <tt>false</tt>, then 
        \ref compute_eg() is used to compute the leptons. 
    */
    virtual void beta_eq_Tfixed(double nB, double T, double &Ye);

    /// Return true if data has been loaded
    bool is_loaded() {
      return loaded;
    }

    /// Free allocated memory
    void free();

    /// Return true if data with lepton information has been loaded
    bool data_with_leptons() {
      return with_leptons;
    }
    
    /// Return true if data with only baryon information has been loaded
    bool data_baryons_only() {
      return baryons_only;
    }

    /* \brief Load EOS from file named \c file_name
    */
    virtual void load(std::string fname, size_t mode);

    /* \brief Output EOS to file named \c file_name
    */
    virtual void output(std::string fname);

    /// Labels for the extra data sets included in current EOS
    std::vector<std::string> oth_names;

    /// Units for the extra data sets included in current EOS
    std::vector<std::string> oth_units;

    /** \brief A slice of data from \ref eos_sn_base for one index fixed
        
        This class allows one to easily construct a \ref
        o2scl::interp2_direct object automatically by fixing one index
        from one of the \ref o2scl::tensor_grid3 objects in a child of
        \ref o2scl::eos_sn_base .
    */
    class slice {
      
    public:
      
      /// Typedef for the matrix type
      typedef std::function<double &(size_t,size_t)> data_t;

      /// Data object in the form of a matrix
      data_t data;

      /// \name Grid vectors
      //@{
      ubvector grid_x, grid_y;
      //@}

      /** \brief The interpolation object
       */
      interp2_direct<ubvector,data_t,matrix_row_gen<data_t>,
                     matrix_column_gen<data_t> > it;
      
      /** \brief Set the slice to correspond to a matrix 
          in the form \f$ (n_B,T) \f$
      */
      void set_nB_T(tensor_grid3<> &tg3, size_t iYe) {
        /*
          data=std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
          (&tensor_grid3<>::get),tg3,std::placeholders::_1,iYe,
          std::placeholders::_2);
          size_t nx=tg3.get_size(0);
          grid_x.resize(nx);
          for(size_t i=0;i<nx;i++) grid_x[i]=tg3.get_grid(0,i);
          size_t ny=tg3.get_size(2);
          grid_y.resize(ny);
          for(size_t i=0;i<ny;i++) grid_y[i]=tg3.get_grid(2,i);
          it.set_data(nx,ny,grid_x,grid_y,data);
        */
        return;
      }
      
      /** \brief Set the slice to correspond to a matrix 
          in the form \f$ (n_B,Y_e) \f$
      */
      void set_nB_Ye(tensor_grid3<> &tg3, size_t iT) {
        /*
          data=std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
          (&tensor_grid3<>::get),tg3,std::placeholders::_1,
          std::placeholders::_2,iT);
          size_t nx=tg3.get_size(0);
          grid_x.resize(nx);
          for(size_t i=0;i<nx;i++) grid_x[i]=tg3.get_grid(0,i);
          size_t ny=tg3.get_size(1);
          grid_y.resize(ny);
          for(size_t i=0;i<ny;i++) grid_y[i]=tg3.get_grid(1,i);
          it.set_data(nx,ny,grid_x,grid_y,data);
        */
        return;
      }
      
      /** \brief Set the slice to correspond to a matrix 
          in the form \f$ (T,Y_e) \f$
      */
      void set_T_Ye(tensor_grid3<> &tg3, size_t inB) {
        /*
          data=std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
          (&tensor_grid3<>::get),tg3,inB,std::placeholders::_2,
          std::placeholders::_1);
          size_t nx=tg3.get_size(2);
          grid_x.resize(nx);
          for(size_t i=0;i<nx;i++) grid_x[i]=tg3.get_grid(2,i);
          size_t ny=tg3.get_size(1);
          grid_y.resize(ny);
          for(size_t i=0;i<ny;i++) grid_y[i]=tg3.get_grid(1,i);
          it.set_data(nx,ny,grid_x,grid_y,data);
        */
        return;
      }
      
    };

  protected:

    /** \brief Unit conversion object (set automatically in constructor)
     */
    convert_units<double> &cu;
    /// If true, a EOS table was successfully loaded (default false)
    bool loaded;
    /// True if thermodynamics with leptons has been loaded
    bool with_leptons;
    /// True if baryon-only thermodynamics has been loaded
    bool baryons_only;

    /// \name Memory allocation
    //@{
    /// Allocate memory
    void alloc();
    //@}

  };

  /** \brief The Lattimer-Swesty supernova EOS 

      This class is experimental.

      \note \o2 Does not contain the Lattimer-Swesty EOS, only
      provides some code to manipulate it. This class is designed
      to be used with the files <tt>ls.dat, sk1.dat, ska.dat</tt>
      and <tt>skm.dat</tt> as provided on Jim Lattimer's website,
      http://www.astro.sunysb.edu/lattimer/EOS/main.html .

      Note that the tables on this website are different
      than what is generated from the LS Fortran code. See
      \ref eos_sn_oo to read O'Connor and Ott's tables
      generated from the LS Fortran code.

      The four models are
      - LS (K=370, Sv=31)
      - SKI' (K=371, Sv=30.4)
      - SKa (K=263, Sv=34.5)
      - SKM* (K=217, Sv=31.4) 

      \note In the original table, the full internal energy per baryon
      (data section 4 of 26) is apparently based on a rest mass of
      \f$ Y_e m_p + (1-Y_e) m_n \f$, while the baryon part of the
      internal energy per baryon (data section 13 of 26) is based
      on a rest mass of \f$ m_n \f$. This means that 
      \f[
      E - E_{\mathrm{int}} = E_{\mathrm{eg}} - Y_e (m_n - m_p)
      \f]
      where \f$ E_{\mathrm{eg}} \f$ is the energy per baryon of
      electrons and photons. In order to keep things consistent with
      the other EOS tables, when the EOS table is loaded, \ref
      eos_sn_base::Eint is rescaled to a rest mass of \f$ Y_e m_p +
      (1-Y_e) m_n \f$ .

      \verbatim embed:rst
      See also the documentation at \ref eos_sn_base and
      the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.

      See also [Lattimer91]_ and [Lattimer85]_.

      .. todo:: 

         - In class eos_sn_ls: There are still a few points for which
           the electron/photon EOS seems to be off, but this may be
           the result of small inaccuracies from finite-differencing
           the LS table.

      \endverbatim
  */
  class eos_sn_ls : public eos_sn_base {

  public:
    
    /// \name Additional data included in this EOS
    //@{
    /// Filling factor for nuclei
    tensor_grid3<> &fill;
    /// Baryon number density inside nuclei in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3<> &nb_in;
    /// Derivative of pressure with respect to baryon density
    tensor_grid3<> &dPdn;
    /// Derivative of pressure with respect to temperature
    tensor_grid3<> &dPdT;
    /// Derivative of pressure with respect to electron fraction
    tensor_grid3<> &dPdY;
    /// Derivative of entropy with respect to temperature
    tensor_grid3<> &dsdT;
    /// Derivative of entropy with respect to electron fraction
    tensor_grid3<> &dsdY;
    /// Number of neutrons in skin
    tensor_grid3<> &Nskin;
    /// Baryon density outside nuclei in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3<> &nb_out;
    /// Proton fraction outside nuclei
    tensor_grid3<> &x_out;
    /** \brief Out of whackness parameter, 
        \f$ \mu_n-\mu_p-\mu_e+1.293~\mathrm{MeV} \f$, in MeV
    */
    tensor_grid3<> &mu;
    //@}

    eos_sn_ls() :
      fill(other[0]),
      nb_in(other[1]),
      dPdn(other[2]),
      dPdT(other[3]),
      dPdY(other[4]),
      dsdT(other[5]),
      dsdY(other[6]),
      Nskin(other[7]),
      nb_out(other[8]),
      x_out(other[9]),
      mu(other[10]) {
    }

    /// Load table from filename \c fname
    virtual void load(std::string fname, size_t mode);

    /** \brief Check electrons and photons

        This checks that the electron and photon thermodynamics
        generated by \o2 is consistent with the data in 
        \c E, \c Eint, \c F, \c Fint, \c P, \c Pint, \c S,
        and \c Sint.
    */
    virtual double check_eg();
  };

  /** \brief The EOS tables from O'Connor and Ott

      \verbatim embed:rst
      This class reads the HDF5 EOS tables generated by E. O'Connor
      and C. Ott in [OConnor10]_. The tables are available from 
      \endverbatim

      http://stellarcollapse.org/equationofstate

      and are available under a creative commons
      attribution-noncommercial-share alike license. This \o2 code to
      read those tables is licensed (along with all \o2 code) under
      the GPLv3 license (with permission from Evan O'Connor).

      The original README file from O'Connor and Ott's EOSdriver
      code is available in the \o2 
      distribution in <tt>doc/o2scl/eos/sphinx/static/scollapse_README</tt>
      and is reproduced below

      \verbatim embed:rst
      .. literalinclude:: ../static/scollapse_README
      \endverbatim

      \verbatim embed:rst
      See also the documentation at \ref eos_sn_base and
      the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.
      \endverbatim

      \future Loading an EOS currently requires loading the HDF5 file
      and then copying it. This wouldn't be necessary if the \o2
      tensor had the same ordering as the indices in the original
      HDF5 file.
  */
  class eos_sn_oo : public eos_sn_base {

  public:
    
    eos_sn_oo() :
      cs2(other[0]),
      dedt(other[1]),
      dpderho(other[2]),
      dpdrhoe(other[3]),
      gamma(other[4]),
      mu_e(other[5]),
      muhat(other[6]),
      munu(other[7]),
      XHe3(other[8]),
      XLi4(other[9]),
      Xt(other[10]), 
      Xd(other[11]) {
    }
    
    /// \name Additional data included in this EOS
    //@{
    /// Speed of sound in cm^2/s^2
    tensor_grid3<> &cs2;
    /// C_V in erg/g/K
    tensor_grid3<> &dedt;
    /// dpderho in dyn*g/cm^2/erg
    tensor_grid3<> &dpderho;
    /// dpdrhoe in dyn cm^3/cm^2/g
    tensor_grid3<> &dpdrhoe;
    /// Gamma
    tensor_grid3<> &gamma;
    /// Electron chemical potential per baryon including rest mass
    tensor_grid3<> &mu_e;
    /// mun - mup
    tensor_grid3<> &muhat;
    /// mue - mun + mup
    tensor_grid3<> &munu;
    /// Helion fraction
    tensor_grid3<> &XHe3;
    /// Lithium-4 fraction
    tensor_grid3<> &XLi4;
    /// Triton fraction
    tensor_grid3<> &Xt;
    /// Deuteron fraction
    tensor_grid3<> &Xd;
    /// The original mass density grid from the table in g/cm^3
    std::vector<double> rho;
    /// Energy shift for table storage in erg/g
    double energy_shift;
    //@}

    /// \name Table modes
    //@{
    /// Use the J. Lattimer et al. method for handling the chemical potentials
    static const size_t ls_mode=0;
    /// Use the H. Shen et al. method for handling the chemical potentials
    static const size_t stos_mode=1;
    /// Set for a Hempel et al. table with light nuclei 
    static const size_t hfsl_mode=2;
    /// Set for a G. Shen et al. table
    static const size_t sht_mode=3;
    //@}
    
    /// Load table from filename \c fname with mode \c mode
    virtual void load(std::string fname, size_t mode);

    /// Desc
    virtual void load_auto(std::string model, std::string directory);
    
  };

  /** \brief The H. Shen et al. supernova EOS
      
      This class is experimental.

      \note \o2 Does not contain the EOS, only provides some code to
      manipulate it. This class is designed to be used with the file
      which was originally called <tt>eos.tab</tt> and now referred to
      as <tt>eos1.tab</tt> and stored e.g. at
      http://user.numazu-ct.ac.jp/~sumi/eos/.
      
      In order to force the EOS to a uniform grid, linear
      interpolation is used to recast the variation in baryon density,
      choosing the grid in baryon density to be the same as the
      section in the table with T=0.1 MeV and \f$ Y_p = 0.1 \f$ for
      all temperature and proton fraction points.

      The data for \ref eos_sn_base::E, \ref eos_sn_base::F, \ref
      eos_sn_base::S, and \ref eos_sn_base::P is not stored in the table
      but can be computed with \ref eos_sn_base::compute_eg().

      \verbatim embed:rst
      See also the documentation at \ref eos_sn_base and
      the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.

      See also [Shen98]_ and [Shen98b]_.
      \endverbatim

      \note Thanks to Matthias Hempel for providing the correct
      temperature grid.

      \verbatim embed:rst
      .. todo:: 

         - In class eos_sn_stos: Add the T=0 and Ye=0 data to this
           class. Separate tables for these cases have been released,
           but I don't think this class can read them yet.

      \endverbatim
  */
  class eos_sn_stos : public eos_sn_base {

  public:

    /// \name Additional data included in this EOS
    //@{
    /** \brief Logarithm of baryon number density in 
        \f$ \mathrm{g}/\mathrm{cm}^3 \f$
    */
    tensor_grid3<> &log_rho;
    /// Baryon number density in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3<> &nB;
    /// Logarithm of proton fraction
    tensor_grid3<> &log_Y;
    /// Proton fraction
    tensor_grid3<> &Yp;
    /// Nucleon effective mass in MeV
    tensor_grid3<> &M_star;
    /// Fraction of quark matter
    tensor_grid3<> &quark_frac;
    //@}

    eos_sn_stos() :
      log_rho(other[0]),
      nB(other[1]),
      log_Y(other[2]),
      Yp(other[3]),
      M_star(other[4]),
      quark_frac(other[5]) {
      check_grid=true;
      m_neut=938.0;
      m_prot=938.0;
    }

    static const size_t orig_mode=0;
    static const size_t quark_mode=1;
    static const size_t fyss_mode=2;

    /// If true, check the grid after load() (default true)
    bool check_grid;

    /// Load table from filename \c fname with mode \c mode
    virtual void load(std::string fname, size_t mode);

    /// Desc
    virtual void load_auto(std::string model, std::string directory);

  };

  /** \brief A class to manipulate the G. Shen et al. EOS

      This class is experimental.

      \note \o2 Does not contain the EOS, only provides some code to
      manipulate it. This class was designed to be used with the FSU
      models given at
      http://cecelia.physics.indiana.edu/gang_shen_eos/FSU/fsu.html .
      The full list of files and the associated modes for the
      \ref load() function are:
      - <tt>"FSU1.7eos1.01.dat"</tt> (\ref mode_17)
      - <tt>"FSU2.1eos1.01.dat"</tt> (\ref mode_21)
      - <tt>"FSU1.7eosb1.01.dat"</tt> (\ref mode_17b)
      - <tt>"FSU2.1eosb1.01.dat"</tt> (\ref mode_21b)
      - <tt>"NL3eos1.03.dat"</tt> (\ref mode_NL3)
      - <tt>"NL3eosb1.03.dat"</tt> (\ref mode_NL3b)

      \verbatim embed:rst
      See also the documentation at \ref eos_sn_base and
      the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.
      \endverbatim

      The free energy per baryon neutron and proton chemical
      potentials are relative to a nucleon mass of 939 MeV. The values
      of \ref o2scl::eos_sn_base::m_neut and \ref
      o2scl::eos_sn_base::m_prot are set to 939 MeV accordingly. The
      electron chemical potential still includes its rest mass. All
      quantites are stored as in the original table, except that
      the values in \ref o2scl::eos_sn_base::E or \ref
      o2scl::eos_sn_base::Eint are computed directly from the
      thermodynamic identity.

      \verbatim embed:rst
      See also [Shen11]_.
      \endverbatim

      \warning The NL3 model is probably ruled out by nuclear mass
      data, neutron matter calculations, and neutron star mass and
      radius observations.
  */
  class eos_sn_sht : public eos_sn_base {

  public:

    /// \name Table modes
    //@{
    /// 1.7 solar masses with leptons and photons
    static const size_t mode_17=0;
    /// 2.1 solar masses with leptons and photons
    static const size_t mode_21=1;
    /// 1.7 solar masses without leptons and photons
    static const size_t mode_17b=2;
    /// 2.1 solar masses without leptons and photons
    static const size_t mode_21b=3;
    /// NL3 model with leptons and photons
    static const size_t mode_NL3=4;
    /// NL3 model with leptons and photons
    static const size_t mode_NL3b=5;
    //@}

    /// \name Additional data included in this EOS
    //@{
    /// Temperature in MeV
    tensor_grid3<> &T;
    /// Proton fraction
    tensor_grid3<> &Yp;
    /// Baryon number density in \f$ 1/\mathrm{fm}^3 \f$
    tensor_grid3<> &nB;
    /// Electron chemical potential in MeV
    tensor_grid3<> &mue;
    /// Nucleon effective mass (Dirac) in MeV
    tensor_grid3<> &M_star;
    //@}
    
    eos_sn_sht() :
      T(other[0]),
      Yp(other[1]),
      nB(other[2]),
      mue(other[3]),
      M_star(other[4]) {
      check_grid=true;
      m_neut=939.0;
      m_prot=939.0;
    }
      
    /// If true, check the grid after load() (default true)
    bool check_grid;

    /// Load table from filename \c fname with mode \c mode
    virtual void load(std::string fname, size_t mode);

  };

  /** \brief The Hempel et al. supernova EOSs
      
      This class is experimental.

      \note \o2 Does not contain the EOS, only provides some code to
      manipulate it. This class was designed to be used with the files
      <tt>dd2_frdm_eos_shen98format_v1.02.tab</tt>,
      <tt>fsg_roca_eos_shen98format_v1.0.tab</tt>, and
      <tt>nl3_lala_eos_shen98format_v1.0.tab</tt> as obtained from
      http://phys-merger.physik.unibas.ch/~hempel/eos.html.

      The free energy is stored with respect to the proton mass
      of 938 MeV, so \ref eos_sn_base::Fint is shifted by
      \f[
      938~\mathrm{MeV}-Y_e m_p-(1-Y_e) m_n
      \f]
      and the internal energy is stored with respect to an
      atomic mass unit so \ref eos_sn_base::Eint is shifted
      by 
      \f[
      931~\mathrm{MeV}-Y_e m_p-(1-Y_e) m_n
      \f]
      the rest of the file data is copied over directly from
      the file. 

      \verbatim embed:rst
      See also the documentation at \ref eos_sn_base and
      the general description in the
      :ref:`Finite-temperature Equation of State Tables`
      section of the User's guide.

      See also [Hempel10]_ and [Hempel12]_.
      \endverbatim
  */
  class eos_sn_hfsl : public eos_sn_base {

  public:

    /// The atomic mass unit
    double m_amu;

    /// \name Additional data included in this EOS
    //@{
    /** \brief Logarithm of baryon number density in 
        \f$ \mathrm{g}/\mathrm{cm}^3 \f$
    */
    tensor_grid3<> &log_rho;
    /// Baryon number density in \f$ 1/\mathrm{fm}^3 \f$
    tensor_grid3<> &nB;
    /// Logarithm of proton fraction
    tensor_grid3<> &log_Y;
    /// Proton fraction
    tensor_grid3<> &Yp;
    /// Nucleon effective mass in MeV
    tensor_grid3<> &M_star;
    /// Mass number of light fragments
    tensor_grid3<> &A_light;
    /// Proton number of light fragments
    tensor_grid3<> &Z_light;
    //@}

    /// If true, check the grid after load() (default true)
    bool check_grid;

    eos_sn_hfsl() :
      log_rho(other[0]),
      nB(other[1]),
      log_Y(other[2]),
      Yp(other[3]),
      M_star(other[4]),
      A_light(other[5]),
      Z_light(other[6]) {
      check_grid=true;
      m_neut=939.565346;
      m_prot=938.272013;
      m_amu=931.49432;
    }
    
    /// Load table from filename \c fname
    virtual void load(std::string fname, size_t mode=0);
    
  };

  /** \brief Dewsc
   */
  class eos_sn_compose : public eos_sn_base {
  
  public:
  
    eos_sn_compose() {
    }
  
    /// Load table from filename \c fname with mode \c mode
    virtual void load() {
    
      //wordexp_single_file(fname);

      std::vector<double> grid;
    
      std::ifstream fin("eos.nb");
      // the first entry is ignored
      fin >> n_nB >> n_nB;
      nB_grid.resize(n_nB);
      for(size_t j=0;j<n_nB;j++) {
        fin >> nB_grid[j];
        grid.push_back(nB_grid[j]);
      }
      fin.close();
    
      std::ifstream fin2("eos.yq");
      // the first entry is ignored
      fin2 >> n_Ye >> n_Ye;
      Ye_grid.resize(n_Ye);
      for(size_t j=0;j<n_Ye;j++) {
        fin2 >> Ye_grid[j];
        grid.push_back(Ye_grid[j]);
      }
      fin2.close();
    
      std::ifstream fin3("eos.t");
      // the first entry is ignored
      fin3 >> n_T >> n_T;
      T_grid.resize(n_T);
      for(size_t j=0;j<n_T;j++) {
        fin3 >> T_grid[j];
        grid.push_back(T_grid[j]);
      }
      fin3.close();

      alloc();
    
      for(size_t i=0;i<n_base+n_oth;i++) {
        arr[i]->set_grid_packed(grid);
      }

      std::ifstream fin4("eos.thermo");
      fin4 >> m_neut;
      fin4 >> m_prot;
    
      double dtemp, dtemp2;
      for(size_t m=0;m<n_T;m++) {
        for(size_t k=0;k<n_Ye;k++) {
          for(size_t j=0;j<n_nB;j++) {
          
            // Skip the grid points since we know them already
            fin4 >> dtemp;
            fin4 >> dtemp;
            fin4 >> dtemp;
          
            fin4 >> dtemp;
            P.set(j,k,m,dtemp*nB_grid[j]);
            fin4 >> dtemp;
            S.set(j,k,m,dtemp);
          
            fin4 >> dtemp;
            mun.set(j,k,m,(dtemp+1.0)*m_neut);
            fin4 >> dtemp2;
            mup.set(j,k,m,dtemp2*m_neut+(dtemp+1.0)*m_neut);
            // Skip the lepton chemical potential
            fin4 >> dtemp;
          
            fin4 >> dtemp;
            F.set(j,k,m,(dtemp+1.0)*nB_grid[j]*m_neut);
            fin4 >> dtemp;
            E.set(j,k,m,(dtemp+1.0)*nB_grid[j]*m_neut);

            // Skip the last column
            fin4 >> dtemp;
          }
        }
      }
      fin4.close();

      std::ifstream fin5("eos.compo");
    
      for(size_t m=0;m<n_T;m++) {
        for(size_t k=0;k<n_Ye;k++) {
          for(size_t j=0;j<n_nB;j++) {
          
            // Skip the grid points since we know them already
            fin5 >> dtemp;
            fin5 >> dtemp;
            fin5 >> dtemp;

            fin5 >> dtemp;
            fin5 >> dtemp;
            fin5 >> dtemp;
            fin5 >> dtemp;
            fin5 >> dtemp;
            fin5 >> dtemp;
          
            // This isn't right yet
            fin5 >> dtemp;
            A.set(j,k,m,dtemp);
            fin5 >> dtemp;
            Z.set(j,k,m,dtemp);
          
            // Skip the last column
            fin5 >> dtemp;
          }
        }
      }
      fin5.close();

      // Loaded must be set to true before calling set_interp()
      n_oth=0;
      loaded=true;
      with_leptons=true;
      baryons_only=false;
    
      if (n_oth!=oth_names.size()) {
        O2SCL_ERR2("Number of names does not match number of data sets ",
                   "in eos_sn_oo::load().",exc_efailed);
      }
    
      // It is important that 'loaded' is set to true before the call to
      // set_interp_type().
      set_interp_type(itp_linear);
    
      if (verbose>0) {
        std::cout << "Done in eos_sn_compose::load()." << std::endl;
      }
    
      return;
    }
  
  };

}

#endif
