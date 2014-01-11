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
#ifndef GEN_SN_EOS_H
#define GEN_SN_EOS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <o2scl/constants.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/table.h>
#include <o2scl/eff_boson.h>
#include <o2scl/rel_fermion.h>
#include <o2scl/eff_fermion.h>
#include <o2scl/test_mgr.h>
#include <o2scl/convert_units.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A base class for the supernova EOSs [abstract]

      This class is experimental.

      See also the general description in the \ref sneos_section 
      section of the User's guide.

      The EOSs are stored in a set of \ref tensor_grid3 objects on
      grids with baryon density in \f$ \mathrm{fm}^{-3} \f$, electron
      fraction (unitless) and temperature in \f$ \mathrm{MeV} \f$.

      Not all tabulated EOSs contain all columns, in which case the
      associated tensor_grid3 object may be empty. For example, EOSs
      which do not contain the leptonic contributions do not provide
      \ref gen_sn_eos::E, \ref gen_sn_eos::F, \ref gen_sn_eos::S, and
      \ref gen_sn_eos::P. In these case, the grid is set for these
      objects but the data is set to zero. To compute these from the
      data after loading the EOS table, use \ref
      gen_sn_eos::compute_eg().

      The functions named <tt>load()</tt> in the children classes load
      the entire EOS into memory. Memory allocation is automatically
      performed, but not deallocated until free() or the destructor is
      called.

      After loading, you can interpolate the EOS by using 
      \ref tensor_grid3::interp_linear() directly. For example,
      the following returns the mass number at an arbitrary
      baryon density, electron fraction, and temperature
      assuming the table is stored in <tt>skm.dat</tt>:
      \verbatim
      ls_eos ls;
      ls.load("skm.dat");
      double nb=0.01, Ye=0.2, T=10.0;
      cout << ls.A.interp_linear(nb,Ye,T) << endl;
      \endverbatim
      Interpolation for all EOSs is linear by default, however, some
      of the grids are logarithmic, so linear interpolation on a
      logarithmic grid leads to power-laws in between grid points.

      \todo Ensure all energies and chemical potentials are based on
      the same rest masses, and document the shifts accordingly.

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
      with this object.
      \future Add option to load and store a separate lepton/photon
      EOS
      \future Add pions?
      \future Create a standard output format? Output to
      stellarcollapse.org HDF5 format?

      \comment
      \future Could this be a child of hadronic_eos_temp and
      then directly used in cold_nstar()? Actually no, this doesn't
      work because of the nuclei. 
      \endcomment

  */
  class gen_sn_eos {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    gen_sn_eos();
    
    virtual ~gen_sn_eos();

    /// \name Grid and data sizes
    //@{
    /// Size of baryon density grid
    size_t n_nB;
    /// Size of electron fraction grid
    size_t n_Ye;
    /// Size of temperature grid
    size_t n_T;
    /// Number of additional data sets
    size_t n_oth;
    /// Number of base data sets
    static const size_t n_base=16;
    //@}

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
    tensor_grid3 F;
    /** \brief Free energy per baryon without lepton and photon 
	contributions in MeV

	By default, this energy is relative to 
	\f[
	m_n (1-Y_e) + m_p Y_e
	\f]
	where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
	is stored in \ref m_prot .
    */
    tensor_grid3 Fint;
    /** \brief Total internal energy per baryon in MeV (without 
	baryon rest masses but including electron rest mass)

	By default, this energy is relative to 
	\f[
	m_n (1-Y_e) + m_p Y_e
	\f]
	where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
	is stored in \ref m_prot .
    */
    tensor_grid3 E;
    /** \brief Internal energy per baryon without lepton and photon 
	contributions in MeV

	By default, this energy is relative to 
	\f[
	m_n (1-Y_e) + m_p Y_e
	\f]
	where \f$ m_n \f$ is stored in \ref m_neut and \f$ m_p \f$
	is stored in \ref m_prot .
    */
    tensor_grid3 Eint;
    /// Total pressure in \f$ \mathrm{MeV}/\mathrm{fm}^3 \f$
    tensor_grid3 P;
    /** \brief Pressure without lepton and photon contributions
	in \f$ \mathrm{MeV}/\mathrm{fm}^3 \f$
    */
    tensor_grid3 Pint;
    /// Total entropy per baryon
    tensor_grid3 S;
    /// Entry per baryon without lepton and photon contributions
    tensor_grid3 Sint;
    /** \brief Neutron chemical potential in MeV

	By default this is relative to the neutron mass in
	\ref m_neut .
     */
    tensor_grid3 mun;
    /** \brief Proton chemical potential in MeV

	By default this is relative to the proton mass in
	\ref m_prot .
    */
    tensor_grid3 mup;
    /// Proton number
    tensor_grid3 Z;
    /// Mass number
    tensor_grid3 A;
    /// Neutron fraction
    tensor_grid3 Xn;
    /// Proton fraction
    tensor_grid3 Xp;
    /// Alpha particle fraction
    tensor_grid3 Xalpha;
    /// Fraction of heavy nuclei
    tensor_grid3 Xnuclei;
    /// Other data sets
    tensor_grid3 other[20];
    /// List of pointers to data
    tensor_grid3 *arr[n_base+20];
    //@}

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
    rel_fermion relf;
    /** \brief Compute the electron and photon contribution for the full
	grid

	If \ref baryons_only_loaded is true, this function computes
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
    int compute_eg();
    //@}

    /** \brief Test the free energy and store results in \c tm

	This checks that the data in \c Fint is consistent with that
	in \c Eint and \c Sint (if \ref baryons_only_loaded is true)
	and that \c F is consistent with that in \c E and \c S (if
	\ref with_leptons_loaded is true).
    */
    void check_free_energy(double &avg);
    
    /** \brief Verbosity parameter (default 1)
     */
    int verbose;

    /** \brief Compute energy per baryon of matter in beta equilibrium
	at zero temperature at a fixed grid point [abstract]

	Given an index \c i for the baryon grid, between 0 and \ref
	n_nB (inclusive), this computes the properties of matter in
	beta equilibrium at zero temperature by finding the electron
	fraction grid point which minimizes \ref E. The baryon density
	is returned in \c nb, the energy per baryon in \c E_beta, the
	pressure in \c P_beta, the electron fraction in \c Ye_beta,
	the proton number in \c Z_beta and the mass number of the
	nucleus in \c A_beta.

	Some tables explicitly contain zero-temperature data which is
	used when it is available. Otherwise, linear interpolation is
	used to extrapolate down to zero from the lowest temperature 
	grid points.
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta)=0;

    /** \brief Compute properties of matter in beta equilibrium
	at s=4

	This function just does a simple hard-coded linear
	interpolation.
    */
    virtual void beta_eq_sfixed
      (size_t i, double entr, double &nb, double &E_beta, 
       double &P_beta, double &Ye_beta, double &Z_beta, double &A_beta,
       double &T_beta);

    /// Return true if data has been loaded
    bool is_loaded() {
      return loaded;
    }

    /// Free allocated memory
    void free();

    /// Return true if data with lepton information has been loaded
    bool data_with_leptons() {
      return with_leptons_loaded;
    }
    
    /// Return true if data with only baryon information has been loaded
    bool data_baryons_only() {
      return baryons_only_loaded;
    }

#ifdef O2SCL_NEVER_DEFINED

    /*
      This is an interesting idea, but e.g. interp2_direct doesn't yet
      handle generic vector and matrix types
    */
    template<class interp2_t> class slice {
      
    public:

      typedef std::function<double(size_t,size_t)> data_t;

      data_t data;

      ubvector grid_x, grid_y;

      interp2_t it;

      void fixed_Ye(tensor_grid3 &tg3, size_t iYe) {
	data=std::bind(std::mem_fn<double &(size_t,size_t,size_t)>
		       (&tensor_grid3::get),tg3,std::placeholders::_1,iYe,
		       std::placeholders::_2);
	size_t nx=tg3.get_size(0);
	grid_x.resize(nx);
	for(size_t i=0;i<nx;i++) grid_x[i]=tg3.get_grid(0,i);
	size_t ny=tg3.get_size(2);
	grid_y.resize(ny);
	for(size_t i=0;i<ny;i++) grid_y[i]=tg3.get_grid(2,i);
	it.set_data(nx,ny,grid_x,grid_y,data);
	return;
      }
      
    };
#endif

  protected:

    /** \brief Unit conversion object (set automatically in constructor)
     */
    convert_units &cu;
    /// If true, a EOS table was successfully loaded (default false)
    bool loaded;
    /// True if thermodynamics with leptons has been loaded
    bool with_leptons_loaded;
    /// True if baryon-only thermodynamics has been loaded
    bool baryons_only_loaded;

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
      \ref oo_eos to read O'Connor and Ott's tables
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
      gen_sn_eos::Eint is rescaled to a rest mass of \f$ Y_e m_p +
      (1-Y_e) m_n \f$ .

      See also the documentation at \ref gen_sn_eos and the
      \ref sneos_section section of the User's guide.

      See \ref Lattimer91 and \ref Lattimer85.

      \todo There are still a few points for which the electron/photon
      EOS seems to be off, but this may be the result of small
      inaccuracies from finite-differencing the LS table.
  */
  class ls_eos : public gen_sn_eos {

  public:
    
    /// \name Additional data included in this EOS
    //@{
    /// Filling factor for nuclei
    tensor_grid3 &fill;
    /// Baryon number density inside nuclei in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3 &nb_in;
    /// Derivative of pressure with respect to baryon density
    tensor_grid3 &dPdn;
    /// Derivative of pressure with respect to temperature
    tensor_grid3 &dPdT;
    /// Derivative of pressure with respect to electron fraction
    tensor_grid3 &dPdY;
    /// Derivative of entropy with respect to temperature
    tensor_grid3 &dsdT;
    /// Derivative of entropy with respect to electron fraction
    tensor_grid3 &dsdY;
    /// Number of neutrons in skin
    tensor_grid3 &Nskin;
    /// Baryon density outside nuclei in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3 &nb_out;
    /// Proton fraction outside nuclei
    tensor_grid3 &x_out;
    /** \brief Out of whackness parameter, 
	\f$ \mu_n-\mu_p-\mu_e+1.293~\mathrm{MeV} \f$, in MeV
    */
    tensor_grid3 &mu;
    //@}

  ls_eos() :
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
    virtual void load(std::string fname);

    /** \brief Check electrons and photons

	This checks that the electron and photon thermodynamics
	generated by \o2 is consistent with the data in 
	\c E, \c Eint, \c F, \c Fint, \c P, \c Pint, \c S,
	and \c Sint.
    */
    int check_eg(test_mgr &tm);

    /** \brief Compute properties of matter in beta equilibrium
	at zero temperature at a baryon density grid point
	
	This EOS table doesn't have T=0 results, so we extrapolate
	from the two low-temperature grid points.
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta) {

      if (loaded==false || with_leptons_loaded==false) {
	O2SCL_ERR2("No data loaded in ",
		   "ls_eos::beta_eq_T0().",exc_einval);
      }
      if (i>=n_nB) {
	O2SCL_ERR2("Too high for baryon grid in ",
		   "ls_eos::beta_eq_T0().",exc_einval);
      }
      // Get baryon density from grid
      nb=E.get_grid(0,i);
      // Find the minima for the two low-temperature grid points
      double ET1=E.get(i,0,0), PT1=P.get(i,0,0);
      double ET2=E.get(i,0,1), PT2=P.get(i,0,1);
      double ZT1=Z.get(i,0,0), ZT2=Z.get(i,0,1);
      double AT1=A.get(i,0,0), AT2=A.get(i,0,1);
      Ye_beta=E.get_grid(1,0);
      for(size_t j=1;j<n_Ye;j++) {
	double Enew=E.get(i,j,0);
	if (Enew<ET1) {
	  ET1=Enew;
	  ET2=E.get(i,j,1);
	  PT1=P.get(i,j,0);
	  PT2=P.get(i,j,1);
	  ZT1=Z.get(i,j,0);
	  ZT2=Z.get(i,j,1);
	  AT1=A.get(i,j,0);
	  AT2=A.get(i,j,1);
	  Ye_beta=E.get_grid(1,j);
	}
      }
      // Now extrapolate to T=0
      E_beta=ET1-(ET2-ET1)*E.get_grid(2,0)/
	(E.get_grid(2,1)-E.get_grid(2,0));
      P_beta=PT1-(PT2-PT1)*P.get_grid(2,0)/
	(P.get_grid(2,1)-P.get_grid(2,0));
      Z_beta=ZT1-(ZT2-ZT1)*Z.get_grid(2,0)/
	(Z.get_grid(2,1)-Z.get_grid(2,0));
      A_beta=AT1-(AT2-AT1)*A.get_grid(2,0)/
	(A.get_grid(2,1)-A.get_grid(2,0));
	
      return;
    }

  };

  /** \brief The EOS tables from O'Connor and Ott

      This class reads the HDF5 EOS tables generated by E. O'Connor
      and C. Ott in \ref OConnor10. The tables are available from 

      http://stellarcollapse.org/equationofstate

      and are available under a creative commons
      attribution-noncommercial-share alike license. This \o2 code to
      read those tables is licensed (along with all \o2 code) under
      the GPLv3 license (with permission from Evan O'Connor).

      The original README file from O'Connor and Ott's EOSdriver
      code is available in the \o2 
      distribution in <tt>doc/o2scl/eos/extras/scollapse_README</tt>
      and is reproduced below

      \verbinclude scollapse_README

      See also the documentation at \ref gen_sn_eos and the
      \ref sneos_section section of the User's guide.

      \future Loading an EOS currently requires loading the HDF5 file
      and then copying it. This wouldn't be necessary if the \o2
      tensor had the same ordering as the indices in the original
      HDF5 file.
  */
  class oo_eos : public gen_sn_eos {

  public:
    
  oo_eos() :
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
    tensor_grid3 &cs2;
    /// C_V in erg/g/K
    tensor_grid3 &dedt;
    /// dpderho in dyn*g/cm^2/erg
    tensor_grid3 &dpderho;
    /// dpdrhoe in dyn cm^3/cm^2/g
    tensor_grid3 &dpdrhoe;
    /// Gamma
    tensor_grid3 &gamma;
    /// Electron chemical potential per baryon including rest mass
    tensor_grid3 &mu_e;
    /// mun - mup
    tensor_grid3 &muhat;
    /// mue - mun + mup
    tensor_grid3 &munu;
    /// Helion fraction
    tensor_grid3 &XHe3;
    /// Lithium-4 fraction
    tensor_grid3 &XLi4;
    /// Triton fraction
    tensor_grid3 &Xt;
    /// Deuteron fraction
    tensor_grid3 &Xd;
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

    /** \brief Compute properties of matter in beta equilibrium
	at zero temperature at a baryon density grid point
	
	This EOS table doesn't have T=0 results, so we extrapolate
	from the two low-temperature grid points.
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta) {

      if (loaded==false || with_leptons_loaded==false) {
	O2SCL_ERR2("No data loaded in ",
		   "oo_eos::beta_eq_T0().",exc_einval);
      }
      if (i>=n_nB) {
	O2SCL_ERR2("Too high for baryon grid in ",
		   "oo_eos::beta_eq_T0().",exc_einval);
      }
      // Get baryon density from grid
      nb=E.get_grid(0,i);
      // Find the minima for the two low-temperature grid points
      double ET1=E.get(i,0,0), PT1=P.get(i,0,0);
      double ET2=E.get(i,0,1), PT2=P.get(i,0,1);
      double ZT1=Z.get(i,0,0), ZT2=Z.get(i,0,1);
      double AT1=A.get(i,0,0), AT2=A.get(i,0,1);
      Ye_beta=E.get_grid(1,0);
      for(size_t j=1;j<n_Ye;j++) {
	double Enew=E.get(i,j,0);
	if (Enew<ET1) {
	  ET1=Enew;
	  ET2=E.get(i,j,1);
	  PT1=P.get(i,j,0);
	  PT2=P.get(i,j,1);
	  ZT1=Z.get(i,j,0);
	  ZT2=Z.get(i,j,1);
	  AT1=A.get(i,j,0);
	  AT2=A.get(i,j,1);
	  Ye_beta=E.get_grid(1,j);
	}
      }

      // Now extrapolate to T=0
      E_beta=ET1-(ET2-ET1)*E.get_grid(2,0)/
	(E.get_grid(2,1)-E.get_grid(2,0));
      P_beta=PT1-(PT2-PT1)*P.get_grid(2,0)/
	(P.get_grid(2,1)-P.get_grid(2,0));
      Z_beta=ZT1-(ZT2-ZT1)*Z.get_grid(2,0)/
	(Z.get_grid(2,1)-Z.get_grid(2,0));
      A_beta=AT1-(AT2-AT1)*A.get_grid(2,0)/
	(A.get_grid(2,1)-A.get_grid(2,0));
	
      return;
    }

  };
  
  /** \brief The Shen et al. supernova EOS
      
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

      Also, the original EOS is tabulated for constant proton
      fraction, and this \o2 interface assumes that the electron
      fraction is equal to the proton fraction. Currently, this is a
      problem only at higher densities where muons might appear.
      
      The data for \ref gen_sn_eos::E, \ref gen_sn_eos::F, \ref
      gen_sn_eos::S, and \ref gen_sn_eos::P is not stored in the table
      but can be computed with \ref gen_sn_eos::compute_eg().

      See also the documentation at \ref gen_sn_eos and the
      \ref sneos_section section of the User's guide.

      See \ref Shen98 and \ref Shen98b .

      \note Thanks to Matthias Hempel for providing the correct
      temperature grid.

      \todo Add the T=0 and Ye=0 data to this class. Separate
      tables for these cases have been released, but I don't think
      this class can read them yet. 
  */
  class stos_eos : public gen_sn_eos {

  public:

    /// \name Additional data included in this EOS
    //@{
    /** \brief Logarithm of baryon number density in 
	\f$ \mathrm{g}/\mathrm{cm}^3 \f$
    */
    tensor_grid3 &log_rho;
    /// Baryon number density in \f$ \mathrm{fm}^{-3} \f$
    tensor_grid3 &nB;
    /// Logarithm of proton fraction
    tensor_grid3 &log_Y;
    /// Proton fraction
    tensor_grid3 &Yp;
    /// Nucleon effective mass in MeV
    tensor_grid3 &M_star;
    /// Fraction of quark matter
    tensor_grid3 &quark_frac;
    //@}

  stos_eos() :
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

    /// If true, check the grid after load() (default true)
    bool check_grid;

    /// Load table from filename \c fname with mode \c mode
    virtual void load(std::string fname, size_t mode);

    /** \brief Compute properties of matter in beta equilibrium
	at zero temperature at a baryon density grid point

	This EOS table doesn't have T=0 results, so we extrapolate
	from the two low-temperature grid points.
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta) {
      if (loaded==false || baryons_only_loaded==false) {
	O2SCL_ERR2("No data loaded in ",
		   "stos_eos::beta_eq_T0().",exc_einval);
      }
      if (i>=n_nB) {
	O2SCL_ERR2("Too high for baryon grid in ",
		   "stos_eos::beta_eq_T0().",exc_einval);
      }
      if (with_leptons_loaded==false) {
	compute_eg();
      }
      // Get baryon density from grid
      nb=E.get_grid(0,i);
      // Find the minima for the two low-temperature grid points
      double ET1=E.get(i,0,0), PT1=P.get(i,0,0);
      double ET2=E.get(i,0,1), PT2=P.get(i,0,1);
      double ZT1=Z.get(i,0,0), ZT2=Z.get(i,0,1);
      double AT1=A.get(i,0,0), AT2=A.get(i,0,1);
      Ye_beta=E.get_grid(1,0);
      for(size_t j=1;j<n_Ye;j++) {
	double Enew=E.get(i,j,0);
	if (Enew<ET1) {
	  ET1=Enew;
	  ET2=E.get(i,j,1);
	  PT1=P.get(i,j,0);
	  PT2=P.get(i,j,1);
	  ZT1=Z.get(i,j,0);
	  ZT2=Z.get(i,j,1);
	  AT1=A.get(i,j,0);
	  AT2=A.get(i,j,1);
	  Ye_beta=E.get_grid(1,j);
	}
      }
      // Now extrapolate to T=0
      E_beta=ET1-(ET2-ET1)*E.get_grid(2,0)/
	(E.get_grid(2,1)-E.get_grid(2,0));
      P_beta=PT1-(PT2-PT1)*P.get_grid(2,0)/
	(P.get_grid(2,1)-P.get_grid(2,0));
      Z_beta=ZT1-(ZT2-ZT1)*Z.get_grid(2,0)/
	(Z.get_grid(2,1)-Z.get_grid(2,0));
      A_beta=AT1-(AT2-AT1)*A.get_grid(2,0)/
	(A.get_grid(2,1)-A.get_grid(2,0));
	
      return;
    }
    
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

      See also the documentation at \ref gen_sn_eos and the
      \ref sneos_section section of the User's guide.

      The proton fraction is assumed to be equal to the electron
      fraction. The free energy per baryon neutron and proton chemical
      potentials are relative to a nucleon mass of 939 MeV. The values
      of \ref o2scl::gen_sn_eos::m_neut and \ref
      o2scl::gen_sn_eos::m_prot are set to 939 MeV accordingly. The
      electron chemical potential still includes its rest mass. All
      quantites are stored as in the original table, except that
      the values in \ref o2scl::gen_sn_eos::E or \ref
      o2scl::gen_sn_eos::Eint are computed directly from the
      thermodynamic identity.

      See \ref Shen11.

      \warning The NL3 model is probably ruled out by nuclear mass
      data, neutron matter calculations, and neutron star mass and
      radius observations.
  */
  class sht_eos : public gen_sn_eos {

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
    tensor_grid3 &T;
    /// Proton fraction
    tensor_grid3 &Yp;
    /// Baryon number density in \f$ 1/\mathrm{fm}^3 \f$
    tensor_grid3 &nB;
    /// Electron chemical potential in MeV
    tensor_grid3 &mue;
    /// Nucleon effective mass in MeV
    tensor_grid3 &M_star;
    //@}
    
  sht_eos() :
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

    /** \brief Compute properties of matter in beta equilibrium
	at zero temperature at a baryon density grid point
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta) {
      if (loaded==false || (with_leptons_loaded==false && 
			    baryons_only_loaded==false)) {
	O2SCL_ERR2("No data loaded in ",
		   "sht_eos::beta_eq_T0().",exc_einval);
      }
      if (i>=n_nB) {
	O2SCL_ERR2("Too high for baryon grid in ",
		   "sht_eos::beta_eq_T0().",exc_einval);
      }
      if (with_leptons_loaded==false) {
	compute_eg();
      }
      // Look up baryon density from grid
      nb=E.get_grid(0,i);
      // Minimize over all electron fractions
      E_beta=E.get(i,0,0);
      P_beta=P.get(i,0,0);
      Ye_beta=E.get_grid(1,0);
      Z_beta=Z.get(i,0,0);
      A_beta=A.get(i,0,0);
      for(size_t j=1;j<n_Ye;j++) {
	double Enew=E.get(i,j,0);
	if (Enew<E_beta) {
	  E_beta=Enew;
	  P_beta=P.get(i,j,0);
	  Ye_beta=E.get_grid(1,j);
	  Z_beta=Z.get(i,j,0);
	  A_beta=A.get(i,j,0);
	}
      }
      return;
    }
    
  };
  
  /** \brief The Hempel et al. supernova EOSs
      
      This class is experimental.

      \note \o2 Does not contain the EOS, only provides some code to
      manipulate it. This class was designed to be used with the files
      <tt>dd2_frdm_eos_shen98format_v1.02.tab</tt>,
      <tt>fsg_roca_eos_shen98format_v1.0.tab</tt>, and
      <tt>nl3_lala_eos_shen98format_v1.0.tab</tt> as obtained from
      http://phys-merger.physik.unibas.ch/~hempel/eos.html.

      See also the documentation at \ref gen_sn_eos and the
      \ref sneos_section section of the User's guide.

      See \ref Hempel10 and \ref Hempel11.
  */
  class hfsl_eos : public gen_sn_eos {

  public:

    /// The atomic mass unit
    double m_amu;

    /// \name Additional data included in this EOS
    //@{
    /** \brief Logarithm of baryon number density in 
	\f$ \mathrm{g}/\mathrm{cm}^3 \f$
    */
    tensor_grid3 &log_rho;
    /// Baryon number density in \f$ 1/\mathrm{fm}^3 \f$
    tensor_grid3 &nB;
    /// Logarithm of proton fraction
    tensor_grid3 &log_Y;
    /// Proton fraction
    tensor_grid3 &Yp;
    /// Nucleon effective mass in MeV
    tensor_grid3 &M_star;
    /// Mass number of light fragments
    tensor_grid3 &A_light;
    /// Proton number of light fragments
    tensor_grid3 &Z_light;
    //@}

    /// If true, check the grid after load() (default true)
    bool check_grid;

  hfsl_eos() :
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
    virtual void load(std::string fname);

    /** \brief Compute properties of matter in beta equilibrium
	at zero temperature at a baryon density grid point
    */
    virtual void beta_eq_T0(size_t i, double &nb, double &E_beta, 
			    double &P_beta, double &Ye_beta,
			    double &Z_beta, double &A_beta) {
      if (loaded==false || baryons_only_loaded==false) {
	O2SCL_ERR2("No data loaded in ",
		   "hfsl_eos::beta_eq_T0().",exc_einval);
      }
      if (i>=n_nB) {
	O2SCL_ERR2("Too high for baryon grid in ",
		   "hfsl_eos::beta_eq_T0().",exc_einval);
      }
      if (with_leptons_loaded==false) {
	compute_eg();
      }
      // Get baryon density from grid
      nb=E.get_grid(0,i);
      // Find the minima for the two low-temperature grid points
      double ET1=E.get(i,0,0), PT1=P.get(i,0,0);
      double ET2=E.get(i,0,1), PT2=P.get(i,0,1);
      double ZT1=Z.get(i,0,0), ZT2=Z.get(i,0,1);
      double AT1=A.get(i,0,0), AT2=A.get(i,0,1);
      Ye_beta=E.get_grid(1,0);
      for(size_t j=1;j<n_Ye;j++) {
	double Enew=E.get(i,j,0);
	if (Enew<ET1) {
	  ET1=Enew;
	  ET2=E.get(i,j,1);
	  PT1=P.get(i,j,0);
	  PT2=P.get(i,j,1);
	  ZT1=Z.get(i,j,0);
	  ZT2=Z.get(i,j,1);
	  AT1=A.get(i,j,0);
	  AT2=A.get(i,j,1);
	  Ye_beta=E.get_grid(1,j);
	}
      }
      // Now extrapolate to T=0
      E_beta=ET1-(ET2-ET1)*E.get_grid(2,0)/
	(E.get_grid(2,1)-E.get_grid(2,0));
      P_beta=PT1-(PT2-PT1)*P.get_grid(2,0)/
	(P.get_grid(2,1)-P.get_grid(2,0));
      Z_beta=ZT1-(ZT2-ZT1)*Z.get_grid(2,0)/
	(Z.get_grid(2,1)-Z.get_grid(2,0));
      A_beta=AT1-(AT2-AT1)*A.get_grid(2,0)/
	(A.get_grid(2,1)-A.get_grid(2,0));
	
      return;
    }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
