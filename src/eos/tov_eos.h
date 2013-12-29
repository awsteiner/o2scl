/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_TOV_EOS_H
#define O2SCL_TOV_EOS_H

#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/constants.h>
#include <o2scl/lib_settings.h>
#include <o2scl/interp.h>
#include <o2scl/table_units.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A EOS base class for the TOV solver
  */
  class tov_eos {
    
  public:
    
    tov_eos() {
      baryon_column=false;
      verbose=1;
    }

    virtual ~tov_eos() {}
    
    /// Control for output (default 1)
    int verbose;

    /** \brief Given the pressure, produce the energy and number densities

	The arguments \c P and \c e should always be in \f$
	M_{\odot}/\mathrm{km}^3 \f$ . The argument for \c nb should be
	in \f$ \mathrm{fm}^{-3} \f$ .
	
	If \ref baryon_column is false, then \c nb is unmodified.
    */
    virtual void get_eden(double P, double &e, double &nb)=0;
    
    /** \brief Given the pressure, produce all the remaining quantities 
	
	The argument \c P should always be in \f$
	M_{\odot}/\mathrm{km}^3 \f$ .
    */
    virtual void get_aux(double P, size_t &np, std::vector<double> &auxp) {
      np=0;
      return;
    }
    
    /** \brief Fill a list with strings for the names of the remaining 
	quanities
    */
    virtual void get_names_units(size_t &np, 
				std::vector<std::string> &pnames,
				std::vector<std::string> &punits) {
      np=0;
      return;
    }

    /** \brief Set to true if the baryon density is provided in the
	EOS (default false)
    */
    bool baryon_column;

  };

  /** \brief An EOS for the TOV solver using simple linear
      interpolation and an optional crust EOS

      The simplest usage of this class is simply to use \ref
      read_table() to read a tabulated EOS stored in a \ref
      table_units object and optionally specify a separate crust EOS.

      \note This stores a pointer to the user-specified table,
      so if that pointer becomes invalid, the interpolation will
      fail.

      Alternatively, the user can simply specify objects
      of type <tt>std::vector<double></tt> which store the energy
      density, pressure, and baryon density. This avoids the 
      overhead of the table and allows faster TOV integrations.

      There are two methods to handle the crust-core interface. The
      default, <tt>smooth_trans</tt> uses the crust below pressure \f$
      P_1 \f$ (equal to the value of \ref trans_pres divided by \ref
      trans_width) and the core above pressure \f$ P_2 \f$ (the value
      of \ref trans_pres times \ref trans_width) and then in between
      uses
      \f[
      \varepsilon(P) = [1-\chi(P)] \varepsilon_{\mathrm{crust}} + 
      \chi(P) \varepsilon_{\mathrm{core}}
      \f]
      where the value \f$ \chi(P) \f$ is determined by
      \f[
      \chi(P) = (P-P_1)/(P_2-P_1) \, .
      \f]
      
      Alternatively the <tt>match_line</tt> method uses
      \f$ \varepsilon_1=\varepsilon(P_1) \f$ and 
      \f$ \varepsilon_2=\varepsilon(P_2) \f$ and
      \f[
      \varepsilon(P) = (\varepsilon_2 - \varepsilon_1) \chi 
      + \varepsilon_1 \, .
      \f]

      Internally, energy and pressure are stored in units of \f$
      \mathrm{M}_{\odot}/\mathrm{km}^3 \f$ and baryon density is
      stored in units of \f$ \mathrm{fm}^{-3} \f$ . The user-specified
      EOS table is left as is, and unit conversion is performed as
      needed in get_eden() and other functions from the units
      specified in the input \ref table_units object.

      \todo It might be useful to exit more gracefully when non-finite
      values are obtained in interpolation, analogous to the
      err_nonconv mechanism elsewhere.
     
  */
  class tov_interp_eos : public tov_eos {
    
  public:
    
    tov_interp_eos();

    virtual ~tov_interp_eos();

    /// \name Mode of transitioning between crust and core EOS
    //@{
    int transition_mode;
    static const int smooth_trans=0;
    static const int match_line=1;
    //@}

    /// \name Basic usage
    //@{
    /// Specify the EOS through a table
    void read_table(table_units<> &eosat, std::string s_cole, 
		    std::string s_colp, std::string s_colnb="");

    /** \brief Read the EOS from a set of equal length
	vectors for energy density, pressure, and baryon density
     */
    void read_vectors(size_t n_core, std::vector<double> &core_ed, 
		      std::vector<double> &core_pr, 
		      std::vector<double> &core_nb);
    
    /** \brief Read the EOS from a pair of equal length
	vectors for energy density and pressure
    */
    void read_vectors(size_t n_core, std::vector<double> &core_ed, 
		      std::vector<double> &core_pr);
    //@}
    
    /// \name Crust EOS functions
    //@{

#ifdef O2SCL_NEVER_DEFINED
    /** \brief Set the crust EOS
     */
    void set_low_density_eos(std::string s_nvpath,
			    int s_nvcole=0, int s_nvcolp=1, int s_nvcolnb=2);
#endif

    /// Default crust EOS from \ref Negele73 and \ref Baym71
    void default_low_dens_eos();

    /// Crust EOS from \ref Shen11b
    void sho11_low_dens_eos();

    /** \brief Crust EOS from \ref Steiner12

	Current acceptable values for \c model are <tt>APR</tt>,
	<tt>Gs</tt>, <tt>Rs</tt> and <tt>SLy4</tt>.
    */
    void s12_low_dens_eos(std::string model="SLy4",
			      bool external=false);

    /** \brief Crust EOS from Goriely, Chamel, and Pearson
	
	From \ref Goriely10, \ref Pearson11, and \ref Pearson12 .
     */
    void gcp10_low_dens_eos(std::string model="BSk20",
			  bool external=false);

    /** \brief Crust EOS from \ref Newton13 given L in MeV

	Current acceptable values for \c model are <tt>PNM</tt>
	and <tt>J35</tt>. 
     */
    void ngl13_low_dens_eos(double L, std::string model="PNM",
			     bool external=false);
    
    /** \brief Crust EOS from \ref Newton13 given S and L in MeV
	and a transition density
     */
    int ngl13_low_dens_eos2(double S, double L, double nt,
			      std::string fname="");
    
    /// Compute with no crust EOS
    void no_low_dens_eos() {
      use_crust=false;
      return;
    }
    //@}

    /// \name Functions used by the tov_solve class
    //@{
    /** \brief Given the pressure, produce the energy and number densities

	The arguments \c P and \c e should always be in \f$
	M_{\odot}/\mathrm{km}^3 \f$ . The argument for \c nb should be
	in \f$ \mathrm{fm}^{-3} \f$ .
	
	If the baryon density is not specified, it should be set to
	zero or \ref baryon_column should be set to false
    */
    virtual void get_eden(double pres, double &ed, double &nb);

    /** \brief Given the pressure, produce all the remaining quantities 
	
	The argument \c P should always be in
	\f$ M_{\odot}/\mathrm{km}^3 \f$ .
    */
    virtual void get_aux(double P, size_t &nv, std::vector<double> &auxp);
    
    /** \brief Fill a list with strings for the names of the remaining 
	quanities
    */
    virtual void get_names_units(size_t &np, 
				 std::vector<std::string> &pnames,
				 std::vector<std::string> &punits);
    //@}

    /// \name Other functions
    //@{
    /** \brief Get the energy and baryon density at the specified pressure
	(after matching crust and core)

	If the core EOS came from a \ref o2scl::table_units object,
	then all arguments \c pres, \c ed, and \c nb are in the same
	units as those in the table.
    */
    virtual int get_eden_full(double pres, double &ed, double &nb);

    /** \brief Given the pressure, produce the energy and number
	densities from the user-specified EOS

	If the core EOS came from a \ref o2scl::table_units object,
	then all arguments \c pres, \c ed, and \c nb are in the same
	units as those in the table.
    */
    virtual void get_eden_high(double pres, double &ed, double &nb);

    /** \brief Given the pressure, produce the energy and number
	densities from the crust EOS
	
	If the core EOS came from a \ref o2scl::table_units object,
	then all arguments \c pres, \c ed, and \c nb are in the same
	units as those in the table.
    */
    virtual void get_eden_low(double pres, double &ed, double &nb);
    
    /** \brief Return limiting and transition pressures

	Returns, in order:
	- the highest pressure in the crust EOS
	- the transition pressure
	- the lowest pressure in the core EOS
	in the same pressure units as the user-specified EOS.

	If the user has not yet specified an EOS, then phi will
	be zero and units of \f$ \mathrm{M}_{\odot}/\mathrm{km}^3 \f$
	will be used for \c plow and \c ptrans.
    */
    void get_transition(double &plow, double &ptrans, double &phi);
    
    /** \brief Set the transition pressure and "width"

	Sets the transition pressure and the width (specified as a
	number greater than unity in \c pw) of the transition between
	the two EOSs. The transition should be in the same units of
	the user-specified EOS. The transition is done smoothly using
	linear interpolation between \f$ P=\mathrm{ptrans}/pmathrm{pw}
	\f$ and \f$ P=\mathrm{ptrans} \times pmathrm{pw} \f$.
     */
    void set_transition(double ptrans, double pw);
    //@}

    /// \name Phase designation
    //@{
    static const int icrust=-1;
    static const int itrans=0;
    static const int icore=1;
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Check that the EOS is valid
    void check_eos();
    
    /// \name Crust EOS variables
    //@{
    /// Set to true if we are using a crust EOS (default false)
    bool use_crust;

    /// Energy densities
    std::vector<double> crust_vece;
    /// Pressures
    std::vector<double> crust_vecp;
    /// Baryon densities
    std::vector<double> crust_vecnb;
    /// Number of EOS entries
    int crust_nlines;
    
    /** \brief Highest pressure in crust EOS 
	(in \f$ M_{\odot}/\mathrm{km}^3 \f$)
    */
    double crust_high_pres;
    //@}
    
    /// \name User EOS
    //@{
    /// Energy densities from core EOS
    std::vector<double> core_vece;
    /// Pressures from core EOS
    std::vector<double> core_vecp;
    /// Baryon densities from core EOS
    std::vector<double> core_vecnb;
    /// Full user EOS table
    table_units<> *core_table;
    /// Number of lines in EOS file
    int core_nlines;
    /// Column for energy density in EOS file
    int cole;
    /// Column for pressure in EOS file
    int colp;
    /// Column for baryon density in EOS file
    int colnb;
    /// True if an EOS has been specified
    bool eos_read;
    /// Number of additional columns in the core EOS
    int core_auxp;
    //@}

    /// \name Unit conversion factors for core EOS
    //@{
    /// Unit conversion factor for energy density (default 1.0)
    double efactor;
    /// Unit conversion factor for pressure (default 1.0)
    double pfactor;
    /// Unit conversion factor for baryon density (default 1.0)
    double nfactor;
    //@}

    /// \name Properties of transition
    //@{
    /** \brief Transition pressure (in \f$ M_{\odot}/\mathrm{km}^3 \f$)
     */
    double trans_pres;
    /// Transition width (unitless)
    double trans_width;
    //@}

    /// Internal linear interpolation function
    void interp(const std::vector<double> &x, const std::vector<double> &y, 
		double xx, double &yy, int n1, int n2);

#endif

  };

  /** \brief The Buchdahl EOS for the TOV solver

      Given the EOS
      \f[
      \varepsilon = 12 \sqrt{p_{*} P}- 5 P
      \f]
      the TOV equation has an analytical solution
      \f[
      R=(1-\beta) \sqrt{\frac{\pi}{288 p_{*} G (1-2 \beta)}}
      \f]
      where \f$ \beta = GM/R \f$.

      The central pressure and energy density are
      \f[
      P_c = 36 p_{*} \beta^2 
      \f]
      \f[
      {\varepsilon}_c = 72 p_{*} \beta (1-5 \beta/2) 
      \f]

      Physical solutions are obtained only for \f$ P< 25 p_{*}/144 \f$
      (ensuring that the argument to the square root is positive)
      and \f$ \beta<1/6 \f$ (ensuring that the EOS is not acausal). 

      Based on \ref Lattimer01 .

      \future Figure out what to do with the buchfun() function
  */
  class tov_buchdahl_eos : public tov_eos {
    
  public:

    tov_buchdahl_eos() {
      Pstar=3.2e-5;
    }

    virtual ~tov_buchdahl_eos() {}

    /** \brief The parameter with units of pressure in units of solar
	masses per km cubed (default value \f$ 3.2 \times 10^{-5} \f$
	)
    */
    double Pstar;
    
    /** \brief Given the pressure, produce the energy and number densities
	
	If the baryon density is not specified, it should be set to
	zero or \ref baryon_column should be set to false
    */
    virtual void get_eden(double P, double &e, double &nb) {
      e=12.0*sqrt(Pstar*P)-5.0*P;
      nb=0.0;
      return;
    }
    
    /// Given the pressure, produce all the remaining quantities 
    virtual void get_aux(double P, size_t &np, std::vector<double> &auxp) {
      np=0;
      return;
    }
    
    /** \brief Fill a list with strings for the names of the remaining 
	quanities
    */
    virtual void get_names(size_t &np, std::vector<std::string> &pnames) {
      np=0;
      return;
    }
    
#ifdef O2SCL_NEVER_DEFINED
    
    int buchfun(size_t bv, const std::vector<double> &bx, 
		std::vector<double> &by, int &pa) {
      double u, rp;
      u=bx[0];
      rp=bx[1];
      by[0]=rp*(1.0-beta+u)/(1.0-2.0*beta)-buchrad;
      by[1]=beta/biga/rp*sin(biga*rp);
      return gsl_success;
    }

#endif

  };

  /** \brief Standard polytropic EOS \f$ P = K \varepsilon^{1+1/n} \f$

      The quantity \f$ K \f$ must be in units of 
      \f$ \left(M_{\odot}/km^3\right)^{-1/n} \f$ .
  */
  class tov_polytrope_eos : public tov_eos {
    
  public:

    tov_polytrope_eos() {
      K=1.0;
      n=3.0;
    }

    virtual ~tov_polytrope_eos() {}

    /** \brief Coefficient (default 1.0)
    */
    double K;

    /// Index (default 3.0)
    double n;

    /** \brief Given the pressure, produce the energy and number densities
	
	If the baryon density is not specified, it should be set to
	zero or \ref baryon_column should be set to false
    */
    virtual void get_eden(double P, double &e, double &nb) {
      e=pow(P/K,1.0/(1.0+1.0/n));
      nb=0.0;
      return;
    }
    
    /// Given the pressure, produce all the remaining quantities 
    virtual void get_aux(double P, size_t &np, std::vector<double> &auxp) {
      np=0;
      return;
    }
    
    /** \brief Fill a list with strings for the names of the remaining 
	quanities
    */
    virtual void get_names(size_t &np, std::vector<std::string> &pnames) {
      np=0;
      return;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif


