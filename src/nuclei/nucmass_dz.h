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
#ifndef O2SCL_DZ_MASS_H
#define O2SCL_DZ_MASS_H

/** \file nucmass_dz.h
    \brief File defining \ref o2scl::nucmass_dz_table and other classes
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/nucmass.h>
#include <o2scl/tensor.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Duflo-Zuker mass formula from tables

      The mass formulas from \ref Duflo95 as given in the 
      files <tt>du_zu_28.feb95</tt> and <tt>du_zu_10.feb96</tt> 
      as obtained from http://amdc.in2p3.fr/web/dz.html . These
      data files have been reformatted for \o2 into HDF files
      with names <tt>du_zu_95.o2</tt> and <tt>du_zu_96.o2</tt>.
  */
  class nucmass_dz_table : public nucmass_table {
    
  public:
    
    /** \brief Create a new mass formula object

	The string \c model is either <tt>"95"</tt>
	or <tt>"96"</tt>. 
    */
    nucmass_dz_table(std::string model="96", bool external=false);

    virtual ~nucmass_dz_table();

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);
    
    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
    /// Verify that the constructor properly loaded the table
    bool is_loaded() { return (n>0); }
    
    /// Return the type, \c "nucmass_dz_table".
    virtual const char *type() { return "nucmass_dz_table"; }

    /// Return number of entries
    int get_nentries() { return n; }
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The reference for the original data
    std::string reference;
    
    /// Table containing the data
    table<> data;

    /// The last table index for caching
    int last;

    /// The total number of entries
    int n;
    
#endif
    
  };

  /** \brief The 10-parameter Duflo-Zuker mass formula

      This class is designed to provide essentially identical results
      to the original 10-parameter Duflo-Zuker code (see \ref Duflo95)
      at
      
      http://amdc.in2p3.fr/theory/du_zu_10.feb96fort

      The default values of \ref nucmass::m_neut and \ref
      nucmass::m_prot are adjusted to make sure that the mass
      excesses match the values given in the original.
      
      \todo This appears to be limited for large nuclei because 'i'
      becomes larger than imax and then statements like
      noc[i][j]=moc-ju and noc[i+1][j]=ju become invalid. This needs
      to be more carefully understood and documented. For now,
      is_included() just arbitrarily chooses 240 as a maximum for N
      and Z. 
      \comment
      Are there any bound nuclei for which the arrays aren't 
      sufficient? Maybe not, in which case there isn't really a
      problem. 
      \endcomment
      \todo Document each field.
  */
  class nucmass_dz_fit : public nucmass_densmat {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix<int> ubmatrix_int;

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /** \name Workspace vectors used internally

	These are allocated in the class constructor.
    */
    //@{
    /** \brief Desc
	
	Note that the first index is already 0 indexed in the DZ version
    */
    tensor3<> onp;

    /// Desc
    ubvector y;
    
    /// Desc
    ubvector pp;
    
    /// Desc
    ubvector oei;
    
    /// Desc
    ubvector dei;
    
    /// Desc
    ubvector qx;
    
    /// Desc
    ubvector dx;
    
    /// Desc
    ubvector op;
    
    /// Desc
    ubvector os;

    /// Desc
    ubvector dyda;

    /// Desc
    ubvector_int n2;

    /// Desc
    ubmatrix_int noc;
    //@}
    
#endif

  public:
    
    nucmass_dz_fit();

    virtual ~nucmass_dz_fit();

    /// Coefficients
    ubvector b;
    
    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Return the type, \c "nucmass_dz_fit".
    virtual const char *type() { return "nucmass_dz_fit"; }
    
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

    /** \brief Return the binding energy in MeV
	
	This function reproduces the original function called
	<tt>mass10()</tt>, except that, to be consistent
	with the other \o2 nuclear mass classes, it returns 
	the binding energy with the opposite sign from the
	original.
    */
    virtual double binding_energy(int Z, int N);    

    /** \brief Return the binding energy in MeV
     */
    virtual double binding_energy_d(double Z, double N);

    /** \brief Given \c Z and \c N, return the mass excess in MeV
     */
    virtual double mass_excess(int Z, int N);
    
    /** \brief Given \c Z and \c N, return the mass excess in MeV
     */
    virtual double mass_excess_d(double Z, double N);

  };

  /** \brief The 33-parameter Duflo-Zuker mass formula

      This class is designed to provide essentially identical results
      to the original Duflo-Zuker code (see \ref Duflo95) at
      
      http://amdc.in2p3.fr/theory/dz31.f

      (Two parameters were added by Duflo and Zuker after the fact to
      the original 31-parameter code, and still referred to as
      <tt>dz31.f</tt>.)

      The default values of \ref nucmass::m_neut and \ref
      nucmass::m_prot are adjusted to make sure that the mass
      excesses match the values given in the original.
      
      \todo Document each field.

      Some explanations about the individual terms come from 
      \ref MendozaTemis10 and the work by 
      G. Bertsch at 
      http://www.phys.washington.edu/users/bertsch/duflo2.ps

      - <tt>a[0]</tt>: "Full master term". Density sqaured divided by
      cube root of A. This is the master term which includes the bulk
      energy from the liquid droplet model and the harmonic oscillator
      effects
      - <tt>a[2]</tt>: "Full spin-orbit term +"
      - <tt>a[4]</tt>: "Full spin-orbit term -"
      - <tt>a[6]</tt>: "Full cross term"
      - <tt>a[8]</tt>: "Partial master term"
      - <tt>a[10]</tt>: "Partial spin-orbit term +"
      - <tt>a[12]</tt>: "Partial spin-orbit term -"
      - <tt>a[14]</tt>: "S3", polarizability of the valence spin-orbit shell
      - <tt>a[16]</tt>: "SQ", "QQM", a neutron-proton interaction
      - <tt>a[18], a[19]</tt>: "D3" balance of monopole effects
      - <tt>a[20], a[21], a[24], a[25]</tt>: "QQ+"/"QQ-", Quadrupole
      terms, corresponding to filling equidistant Nilsson orbits
      - <tt>a[22], a[23]</tt>: "D0", Loss of monopole and gain of quadrupole
      energy for intruders
      - <tt>a[26]</tt>: "TT", Symmetry energy
      - <tt>a[28]</tt>: "SS"
      - <tt>a[30]</tt>: "C", Coulomb energy
      - <tt>a[31]</tt>: "P0", First pairing energy term
      - <tt>a[32]</tt>: "P1", Second pairing energy term

      For odd parameters up to <tt>a[29]</tt>, the odd parameter is just
      the preceeding even term divided by the cube root of A. 

      Note that the original code states that, <tt>"for i even
      a(i,program) =a(i-1,paper)*a(i,paper)"</tt>.
  */
  class nucmass_dz_fit_33 : public nucmass_densmat {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix<int> ubmatrix_int;

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /** \name Workspace vectors used internally

	These are allocated in the class constructor.
    */
    //@{
    ubvector dyda, fyda, fyd0, onps, oei, dei, op2, ym, op1;
    ubvector shell, sshell;
    tensor3<> op, onp, ot;
    ubvector_int n4, nn, jup, jud, n2;
    ubmatrix_int noc;
    //@}
    
#endif

  public:
    
    nucmass_dz_fit_33();

    virtual ~nucmass_dz_fit_33();

    /// Coefficients
    ubvector a;
    
    /// Return the type, \c "nucmass_dz_fit_33".
    virtual const char *type() { return "nucmass_dz_fit_33"; }
    
    /// Fix parameters from an array for fitting
    virtual int fit_fun(size_t nv, const ubvector &x);
    
    /// Fill array with guess from present values for fitting
    virtual int guess_fun(size_t nv, ubvector &x);

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /** \brief Return the binding energy in MeV
	
	This function reproduces the original function called
	<tt>EMASSDZ()</tt>, except that, to be consistent
	with the other \o2 nuclear mass classes, it returns 
	the binding energy with the opposite sign from the
	original.
    */
    virtual double binding_energy(int Z, int N);
    
    /** \brief Return the binding energy in MeV
     */
    virtual double binding_energy_d(double Z, double N);

    /** \brief Given \c Z and \c N, return the mass excess in MeV
     */
    virtual double mass_excess(int Z, int N);
    
    /** \brief Given \c Z and \c N, return the mass excess in MeV
     */
    virtual double mass_excess_d(double Z, double N);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
