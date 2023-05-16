/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_AME_MASS_H
#define O2SCL_AME_MASS_H

/** \file nucmass_ame.h
    \brief File defining \ref o2scl::nucmass_ame
*/

#include <cmath>
#include <string>
#include <map>
#include <o2scl/nucleus.h>
#include <o2scl/constants.h>
#include <o2scl/table.h>
#include <o2scl/nucmass.h>

// Forward definition of the nucmass_ame class for HDF I/O
namespace o2scl {
  class nucmass_ame;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void ame_load_ext(o2scl::nucmass_ame &ame, std::string file_name, 
		    std::string table_name, bool exp_only);
  void ame_load(o2scl::nucmass_ame &ame, std::string name,
		bool exp_only);
}
  
//#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
  //#endif

  /** \brief Masses from the Atomic Mass Evaluation 
      
      \note This class requires data stored in an HDF file and
      thus requires HDF support for normal usage.

      \verbatim embed:rst
      This class provides an interface to the atomic mass table using
      data from [Audi95]_, [Audi03]_, [Audi12]_, [Wang12]_,
      [Huang17]_, [Wang17]_, [Huang21ta]_, and [Wang21ta]_.
      \endverbatim

      To load data from the \o2 HDF5 data files, use
      <tt>o2scl_hdf::ame_load()</tt> .

      The 1995 data provided the binding energy (stored in
      nucmass_ame::entry::be and nucmass_ame::entry::dbe), while the 2003
      data provided the binding energy divided by the mass number
      (stored in nucmass_ame::entry::beoa and nucmass_ame::entry::dbeoa).
      When the 1995 data is used, nucmass_ame::entry::beoa and
      nucmass_ame::entry::dbeoa are calculated automatically, and when
      the 2003 data is used nucmass_ame::entry::be and
      nucmass_ame::entry::dbe are calculated automatically. To indicate
      that \o2 has automatically calculated a value in this way, the
      associated accuracy field is set to \ref
      o2scl::nucmass_ame::intl_computed.
      
      Note that all uncertainties are 1 sigma uncertainties.

      The functions \ref mass_excess() and \ref
      o2scl::nucmass::mass_excess_d() directly return the value from the
      data. For consistency, the functions \ref
      o2scl::nucmass::binding_energy(), \ref
      o2scl::nucmass::binding_energy_d(), \ref
      o2scl::nucmass::total_mass(), and \ref
      o2scl::nucmass::total_mass_d() return values which are
      automatically computed from the mass excess with the neutron and
      proton mass in \ref m_neut and \ref m_prot. In order to obtain
      the value of the binding energy as reported in the original data
      instead of the value computed from the mass excess, you can use
      the function \ref get_ZN(), and access the corresponding entry
      from the data directly.

      In cases where the decimal point in the original table was
      replaced with a <tt>#</tt>, the associated accuracy field is set
      to \ref o2scl::nucmass_ame::estimated. In cases where the original
      table contained a asterisk to indicate a value was not
      calculable, the accuracy field is set to \ref
      o2scl::nucmass_ame::not_calculable and the value is set to zero. If
      \o2 internally computed the value because it was not present in
      the original table, the accuracy field is set to \ref
      o2scl::nucmass_ame::intl_computed. In cases where either \ref
      o2scl::nucmass_ame::entry::orig or \ref
      o2scl::nucmass_ame::entry::bdmode in the original table was blank,
      the string is set to <tt>"blank"</tt>.

      In the original table, binding energies are defined with a
      positive sign, so that lead has a binding energy of +8 MeV and
      this is what is stored in \ref o2scl::nucmass_ame::entry::be.
      However, for consistency with the other mass formulas, \ref
      o2scl::nucmass_ame::binding_energy() gives, e.g., about \f$ -8
      \f$ MeV for lead. See also the documentation for the class
      structure for each table entry in \ref
      o2scl::nucmass_ame::entry.
      
      \future Create a caching and more intelligent search system for
      the table. The table is sorted by A and then N, so we could
      probably just copy the search routine from mnmsk_mass, which is
      sorted by Z and then N (some code written for this, but 
      it doesn't work yet). 
      \future Should m_neut and m_prot be set to the neutron and
      proton masses from the table by default?
  */
  class nucmass_ame : public nucmass_table {
    
  public:

    friend void o2scl_hdf::ame_load_ext(nucmass_ame &ame,
					std::string file_name, 
					std::string table_name,
					bool exp_only);
    
    friend void o2scl_hdf::ame_load(nucmass_ame &ame, std::string name,
				    bool exp_only);

    /// Create an AME mass object
    nucmass_ame();

    ~nucmass_ame();
    
    /// \name Accuracy modes
    //@{
    /// Measured value from source data
    static const int measured=0;
    /// Value estimated in source data
    static const int estimated=1;
    /// Value listed in data as not calculable
    static const int not_calculable=2;
    /// Value computed by \o2
    static const int intl_computed=3;
    /// Value computed by \o2
    static const int unc_less_than_half_eV=4;
    //@}

    /** \brief Atomic mass entry structure

	Atomic mass entry data object for \ref o2scl::nucmass_ame.

	This has to be a struct, not a class, so that it can
	be processed by the HDF5 make table functions.
    */
    struct entry {

    public:

      /// N-Z
      int NMZ;
    
      /// Neutron number
      int N;
    
      /// Proton number
      int Z;
    
      /// Mass number
      int A;
    
      /// Element name
      char el[4];
    
      /// Data origin
      char orig[5];
    
      /// Mass excess (in keV)
      double mass;
    
      /// Mass excess uncertainty (in keV)
      double dmass;

      /// Mass accuracy flag 
      int mass_acc;
    
      /// Binding energy (in keV, given in the '95 data)
      double be;
    
      /// Binding energy uncertainty (in keV, given in the '95 data)
      double dbe;

      /// Binding energy accuracy flag
      int be_acc;
    
      /// Binding energy / A (in keV, given in the '03 data)
      double beoa;
    
      /// Binding energy / A uncertainty (in keV, given in the '03 data)
      double dbeoa;
    
      /// Binding energy / A accuracy flag
      int beoa_acc;

      /// Beta decay mode
      char bdmode[3];
    
      /// Beta-decay energy (in keV)
      double bde;

      /// Beta-decay energy uncertainty (in keV)
      double dbde;
    
      /// Beta-decay energy accuracy flag
      int bde_acc;

      /// Mass number (reported twice in original table)
      int A2;
    
      /// Atomic mass (in keV)
      double amass;
    
      /// Atomic mass uncertainty (in keV)
      double damass;
    
      /// Atomic mass accuracy flag
      int amass_acc;

    };
    
    /// Return the type, \c "nucmass_ame".
    virtual const char *type() { return "nucmass_ame"; }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
    /// Get element with Z=l_Z and N=l_N (e.g. 82,126).
    entry get_ZN(int l_Z, int l_N);
    
    /// Get element with Z=l_Z and A=l_A (e.g. 82,208).
    entry get_ZA(int l_Z, int l_A);
    
    /// Get element with name l_el and A=l_A (e.g. "Pb",208).
    entry get_elA(std::string l_el, int l_A);
    
    /// Get element with string (e.g. "Pb208")
    entry get(std::string nucleus);
    
    /// Returns true if data has been loaded
    bool is_loaded() { return (n>0); }

    /// Return number of entries
    size_t get_nentries() { return n; }

    /// Return the reference
    std::string get_reference() { return reference; }
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The array containing the mass data of length ame::n
	
	\comment
	Ideally I'd prefer to store a vector<entry> rather than 
	a pointer, but the pointer is required to read the
	HDF5 table.
	\endcomment
     */
    entry *mass;
    
    /// The last table index for caching
    int last;
    
#endif

  };

  /** \brief 
      
      \note This version avoids storing a raw pointer, supports
      the NUBASE data, and provides more field information
   */
  class nucmass_ame2 : public nucmass_table {
    
  public:

    /// Create an AME mass object
    nucmass_ame2();

    virtual ~nucmass_ame2();
    
    /// \name Accuracy modes
    //@{
    /// Measured value from source data
    static const int measured=0;
    /// Value estimated in source data
    static const int estimated=1;
    /// Value listed in data as not calculable
    static const int not_calculable=2;
    /// Value computed by \o2
    static const int intl_computed=3;
    /// The uncertainty is less than half of an eV
    static const int unc_less_than_half_eV=4;
    /// The value is blank
    static const int blank=5;
    /// The uncertainty is blank
    static const int blank_unc=6;
    /// The nucleus is unstable
    static const int unstable=7;
    /// 
    static const int part_unstable=8;
    /// 
    static const int lower_limit=9;
    /// 
    static const int upper_limit=10;
    /// 
    static const int does_not_exist=11;
    /// The stored value is only approximate
    static const int approximate=12;
    //@}

    /** \brief Data structure for a row in the NUBASE file

        \note All character arrays have space for an extra terminating
        character
     */
    class entry_nubase_20 {

    public:
      
      /// Note
      char Znote;
      /// Element name
      char A_el[6];
      /// Type field (isomer, resonance, etc.)
      char isomer;
      /// Mass excess (in keV)
      double mass;
      /// Mass excess uncertainty (in keV)
      double dmass;
      /// Desc
      int mass_acc;
      /// Excitation energy
      double exc_energy;
      /// Excitation energy uncertainty
      double dexc_energy;
      /// Desc
      int exc_energy_acc;
      /// Excitation energy origin
      char origin[3];
      /// Isomer uncertainty
      char isomer_unc;
      /// Isomer inversion
      char isomer_inv;
      /// Half-life
      double hlife;
      /// Desc
      int hlife_acc;
      /// Half-life unit
      char hl_unit[3];
      /// Half-life uncertainty
      char dhlife[8];
      /// Spin and parity
      char spinp[15];
      /// Isospin multiplet
      int ENSDF_year;
      /// Year of discovery
      int discovery;
      /// Decay mode and intensity
      char decay_intensity[91];
      
    };

    /** \brief Atomic mass entry structure

	Atomic mass entry data object for \ref o2scl::nucmass_ame.

	This has to be a struct, not a class, so that it can
	be processed by the HDF5 make table functions.
    */
    class entry {

    public:

      /// N-Z
      int NMZ;
    
      /// Neutron number
      int N;
    
      /// Proton number
      int Z;
    
      /// Mass number
      int A;
    
      /// Element name
      char el[4];
    
      /// Data origin
      char orig[5];
    
      /// Mass excess (in keV)
      double mass;
    
      /// Mass excess uncertainty (in keV)
      double dmass;

      /// Mass accuracy flag 
      int mass_acc;
    
      /// Binding energy (in keV, given in the '95 data)
      double be;
    
      /// Binding energy uncertainty (in keV, given in the '95 data)
      double dbe;

      /// Binding energy accuracy flag
      int be_acc;
    
      /// Binding energy / A (in keV, given in the '03 data)
      double beoa;
    
      /// Binding energy / A uncertainty (in keV, given in the '03 data)
      double dbeoa;
    
      /// Binding energy / A accuracy flag
      int beoa_acc;

      /// Beta decay mode
      char bdmode[3];
    
      /// Beta-decay energy (in keV)
      double bde;

      /// Beta-decay energy uncertainty (in keV)
      double dbde;
    
      /// Beta-decay energy accuracy flag
      int bde_acc;

      /// Mass number (reported twice in original table)
      int A2;
    
      /// Atomic mass (in keV)
      double amass;
    
      /// Atomic mass uncertainty (in keV)
      double damass;
    
      /// Atomic mass accuracy flag
      int amass_acc;

      /// NUBASE properties
      std::vector<entry_nubase_20> props;
      
    };
    
    /** \brief Read data for \ref o2scl::nucmass_ame from an HDF table
        specified in a file

        \todo Why was this moved here from hdf_nucmass.h ?
    */
    void load(std::string name="20", bool exp_only=false) {
      std::string file_name, nubase_name;
      file_name=o2scl::o2scl_settings.get_data_dir()+"/nucmass";
      if (name=="95exp") {
        file_name+="/ame95exp.o2";
      } else if (name=="95rmd") {
        file_name+="/ame95rmd.o2";
      } else if (name=="03round") {
        file_name+="/ame03round.o2";
      } else if (name=="03") {
        file_name+="/ame03.o2";
      } else if (name==((std::string)"12")) { 
        file_name+="/ame12.o2";
      } else if (name==((std::string)"16")) { 
        file_name+="/ame16.o2";
      } else if (name==((std::string)"16round")) { 
        file_name+="/ame16round.o2";
      } else if (name==((std::string)"20")) { 
        nubase_name=file_name;
        file_name+="/mass.mas20.txt";
        nubase_name+="/nubase_1.mas20.txt";
      } else if (name==((std::string)"20round")) { 
        file_name+="/massround.mas20.txt";
      } else {
        std::string s=((std::string)"Invalid name '")+name+
          "' in o2scl_hdf::ame_load().";
        O2SCL_ERR(s.c_str(),exc_einval);
      }
      
      load_ext(name,file_name,nubase_name,exp_only);
      return;
    }
    
    /** \brief Read data for \ref o2scl::nucmass_ame from an HDF table
        specified in a file
        
        \todo Why was this moved here from hdf_nucmass.h ?
    */
    void load_ext(std::string name, std::string filename,
                  std::string nubase_file, bool exp_only=false);
    
    /** \brief Parse strings \c s1 and \c s2 from the AME into a value,
        \c d1, an uncertainty, \c d2, and an accuracy flag, \c acc
        
        - If string \c s1 has an asterisk, then \c d1 and \c d2 are
        set to zero and \c acc is set to \ref nucmass_ame::not_calculable.
        - If string \c s2 contains the letter 'a', then \c d2 is set to
        zero and \c ass is set to \ref nucmass_ame::unc_less_than_half_eV.
        The value of d1 is computed from <tt>stod_nothrow()</tt>.
        - Otherwise, if string \c s1 has a pound sign, then \c acc is set
        to \ref nucmass_ame::estimated, otherwise, \c acc is set to \ref
        nucmass_ame::measured. The values of \c d1 and \c d2 are computed
        from \c s1 and \c s2 using <tt>stod_nothrow()</tt>.
        - If either of the stod_nothrow() calls returns a non-zero value,
        then the error handler is called.
    */
    int parse(std::string s1, std::string s2, double &d1, double &d2,
              int &acc) {
      //std::cout << "Hx: '" << s1 << "' '" << s2 << "'." << std::endl;
      if (count_words(s1)==0) {
        d1=0.0;
        d2=0.0;
        acc=nucmass_ame2::blank;
        return 0;
      }
      if (s1.find('*')!=std::string::npos) {
        d1=0.0;
        d2=0.0;
        acc=nucmass_ame2::not_calculable;
        return 0;
      } 
      if (s1.find('#')!=std::string::npos) {
        acc=nucmass_ame2::estimated;
      } else {
        acc=nucmass_ame2::measured;
      }
      //std::cout << "Hz: '" << s1 << "'" << std::endl;
      int ret1=o2scl::stod_nothrow(s1,d1);
      if (ret1!=0) {
        std::cerr << "Failed to convert: '" << s1 << "'." << std::endl;
        O2SCL_ERR("Failed to convert first string in parse().",
                  o2scl::exc_einval);
      }
      if (count_words(s2)==0) {
        d2=0.0;
        acc=nucmass_ame2::blank_unc;
      } else if (s2.find('a')!=std::string::npos) {
        d2=0.0;
        acc=nucmass_ame2::unc_less_than_half_eV;
      } else {
        //std::cout << "Hy: '" << s2 << "'" << std::endl;
        int ret2=o2scl::stod_nothrow(s2,d2);
        //std::cout << "ret2: " << ret2 << std::endl;
        if (ret2!=0) {
          std::cerr << "Failed to convert: '" << s2 << "'." << std::endl;
          O2SCL_ERR("Failed to convert second string in parse().",
                    o2scl::exc_einval);
        }
      }
      return 0;
    }
    
    /// Return the type, \c "nucmass_ame2".
    virtual const char *type() { return "nucmass_ame2"; }

    /** \brief Return false if the mass formula does not include 
	specified nucleus
    */
    virtual bool is_included(int Z, int N);

    /// Given \c Z and \c N, return the mass excess in MeV
    virtual double mass_excess(int Z, int N);
    
    /// Get element with Z=l_Z and N=l_N (e.g. 82,126).
    entry get_ZN(int l_Z, int l_N);
    
    /// Get element with Z=l_Z and A=l_A (e.g. 82,208).
    entry get_ZA(int l_Z, int l_A);
    
    /// Get element with name l_el and A=l_A (e.g. "Pb",208).
    entry get_elA(std::string l_el, int l_A);
    
    /// Get element with string (e.g. "Pb208")
    entry get(std::string nucleus);
    
    /// Returns true if data has been loaded
    bool is_loaded() { return (mass.size()>0); }

    /// Return number of entries
    size_t get_nentries() { return mass.size(); }

    /// Return the reference
    std::string get_reference() { return reference; }
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The array containing the mass data of length ame::n
     */
    std::vector<entry> mass;
    
    /// The last table index for caching
    int last;
    
#endif

  };
  
  //#ifndef DOXYGEN_NO_O2NS
}
//#endif

#endif
