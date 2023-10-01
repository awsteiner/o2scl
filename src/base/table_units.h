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
#ifndef O2SCL_TABLE_UNITS_H
#define O2SCL_TABLE_UNITS_H

/** \file table_units.h
    \brief File defining \ref o2scl::table_units 
*/

#include <o2scl/table.h>
#include <o2scl/lib_settings.h>

// For the MPI table send and receive functions below
#ifdef O2SCL_MPI
#include <mpi.h>
#endif

// Forward definition of the table_units class for HDF I/O
namespace o2scl {
  template<class vec_t> class table_units;
}

// Forward definition of HDF I/O to extend friendship in table_units
namespace o2scl_hdf { 

  class hdf_file; 
  
  template<class vec_t>
  void hdf_input(hdf_file &hf, o2scl::table_units<vec_t> &t, 
                 std::string name="");
  
  void hdf_output
    (hdf_file &hf, 
     o2scl::table_units<std::vector<double> > &t, 
     std::string name);

  template<class vec_t>
    void hdf_input_data(hdf_file &hf, o2scl::table_units<vec_t> &t);

  void hdf_output_data
    (hdf_file &hf, 
     o2scl::table_units<std::vector<double> > &t);

}

namespace o2scl {

  /** \brief Data \table class with units

      \future The unit conversion object is now always a pointer to
      the global conversion object. This could be modified, so that
      each table can have it's own, but this might require HDF output
      of unit conversion objects.

      \future Make table methods virtual? (not necessary yet since
      delete_column() isn't referred to internally)
  */
  template<class vec_t=std::vector<double> > 
  class table_units : public table<vec_t,double> {
    
  public:
  
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /** \brief Create a new table_units with space for nlines<=cmaxlines.
     */
  table_units(int cmaxlines=0) : table<vec_t,double>(cmaxlines) {
      cup=&o2scl_settings.get_convert_units();
    }
    
    virtual ~table_units() {
      utree.clear();
    }

    /// \name Copy constructors
    //@{
    /// Copy with constructor from \ref table_units
  table_units(const table_units &t) : table<vec_t,double>(t.get_nlines()) {
  
      // Copy constants 
      this->constants=t.constants;
  
      // Copy interpolation type
      this->itype=t.get_interp_type();
    
      // Copy the columns and data
      this->nlines=t.get_nlines();
      this->maxlines=this->nlines;
      for(size_t i=0;i<t.get_ncolumns();i++) {

	// Column name
	std::string cname=t.get_column_name(i);

	// Insert column into tree
	typename table<vec_t,double>::col s;
	s.dat.resize(this->nlines);
	s.index=this->atree.size();
	this->atree.insert(make_pair(cname,s));

	// Insert in iterator index
	typename table<vec_t,double>::aiter it=this->atree.find(cname);
	this->alist.push_back(it);
    
	// Insert in unit list
	utree.insert(make_pair(cname,t.get_unit(cname)));
    
	// Fill the data
	for(size_t j=0;j<t.get_nlines();j++) {
	  it->second.dat[j]=t.get(cname,j);
	}
    
      }

      this->intp_set=false;
      cup=&o2scl_settings.get_convert_units();

      this->is_valid();

      return;
    }

    /// Copy with constructor from \ref table
  table_units(const table<vec_t,double> &t) :
    table<vec_t,double>(t.get_nlines()) {
  
      // Copy constants 
      this->constants=t.constants;
  
      // Copy interpolation type
      this->itype=t.get_interp_type();
    
      // Copy the columns and data
      this->nlines=t.get_nlines();
      this->maxlines=this->nlines;
      for(size_t i=0;i<t.get_ncolumns();i++) {
	
	// Column name
	std::string cname=t.get_column_name(i);

	// Insert column into tree
	typename table<vec_t,double>::col s;
	s.dat.resize(this->nlines);
	s.index=this->atree.size();
	this->atree.insert(make_pair(cname,s));

	// Insert in iterator index
	typename table<vec_t,double>::aiter it=this->atree.find(cname);
	this->alist.push_back(it);
    
	// Fill the data
	for(size_t j=0;j<t.get_nlines();j++) {
	  it->second.dat[j]=t.get(cname,j);
	}
    
      }

      this->intp_set=false;
      cup=&o2scl_settings.get_convert_units();

      this->is_valid();

      return;
    
    }

    /// Copy with <tt>operator=</tt> from \ref table_units
    table_units &operator=(const table_units &t) {

      if (this!=&t) {

	this->clear();

	// Copy constants 
	this->constants=t.constants;
  
	// Copy interpolation type
	this->itype=t.get_interp_type();
	
	// Copy the columns and data
	this->nlines=t.get_nlines();
	this->maxlines=this->nlines;
	for(size_t i=0;i<t.get_ncolumns();i++) {

	  // Column name
	  std::string cname=t.get_column_name(i);

	  // Insert column into tree
	  typename table<vec_t,double>::col s;
	  s.dat.resize(this->nlines);
	  s.index=this->atree.size();
	  this->atree.insert(make_pair(cname,s));

	  // Insert in iterator index
	  typename table<vec_t,double>::aiter it=this->atree.find(cname);
	  this->alist.push_back(it);
    
	  // Insert in unit list
	  utree.insert(make_pair(cname,t.get_unit(cname)));
    
	  // Fill the data
	  for(size_t j=0;j<t.get_nlines();j++) {
	    it->second.dat[j]=t.get(cname,j);
	  }
    
	}

	if (this->intp_set) {
	  this->intp_set=false;
	  delete this->si;
	}

	cup=&o2scl_settings.get_convert_units();

      }

      this->is_valid();

      return *this;
    }

    /// Copy with <tt>operator=</tt> from \ref table
    table_units &operator=(const table<vec_t,double> &t) {

      if (this!=&t) {
  
	this->clear();

	// Copy constants 
	this->constants=t.constants;
  
	// Copy interpolation type
	this->itype=t.get_interp_type();
    
	// Copy the columns and data
	this->nlines=t.get_nlines();
	this->maxlines=this->nlines;
	for(size_t i=0;i<t.get_ncolumns();i++) {

	  // Column name
	  std::string cname=t.get_column_name(i);

	  // Insert column into tree
	  typename table<vec_t,double>::col s;
	  s.dat=new vec_t(this->nlines);
	  s.index=this->atree.size();
	  this->atree.insert(make_pair(cname,s));

	  // Insert in iterator index
	  typename table<vec_t,double>::aiter it=this->atree.find(cname);
	  this->alist.push_back(it);
    
	  // Fill the data
	  for(size_t j=0;j<t.get_nlines();j++) {
	    (*it->second.dat)[j]=t.get(cname,j);
	  }
    
	}

	if (this->intp_set) {
	  this->intp_set=false;
	  delete this->si;
	}

	cup=&o2scl_settings.get_convert_units();

      }

      this->is_valid();

      return *this;
    }
    //@}

    /** \brief Copy all rows matching a particular condition to
	a new table

	This function begins by ensuring that all columns in the
	current table are present in \c dest, creating new columns
	(and copying their units) in \c dest if necessary. It then
	copies all rows where \c func evaluates to a number greater
	than 0.5 to table \c dest by adding rows at the end of the
	table.
    */
    template<class vec2_t>
      void copy_rows(std::string func, table_units<vec2_t> &dest) {

      // Set up columns
      for(size_t i=0;i<this->get_ncolumns();i++) {
	std::string cname=this->get_column_name(i);
	if (dest.is_column(cname)==false) {
	  dest.new_column(cname);
	}
	dest.set_unit(cname,get_unit(cname));
      }
    
      size_t new_lines=dest.get_nlines();
      for(size_t i=0;i<this->get_nlines();i++) {
	double val=this->row_function(func,i);
	if (val>0.5) {
	  this->set_nlines_auto(new_lines+1);
	  for(size_t j=0;j<this->get_ncolumns();j++) {
	    std::string cname=this->get_column_name(j);
	    dest.set(cname,new_lines,this->get(cname,i));
	  }
	  new_lines++;
	}
      }

      return;
    }

    /// \name Unit manipulation
    //@{
    /// Get the unit for column \c scol
    std::string get_unit(std::string scol) const {
  
      uciter it=utree.find(scol);
      if (it==utree.end()) {
	// Not found in unit entry, look for column of data
	typename table<vec_t,double>::aciter at=this->atree.find(scol);
	if (at==this->atree.end()) {
	  O2SCL_ERR((((std::string)"Column '")+scol+
		     "' not found in table_units::get_unit().").c_str(),
		    exc_enotfound);
	} else {
	  return "";
	}
      }

      return it->second;
    }

    /** \brief Specify the units as a string separated by spaces
     */
    void line_of_units(std::string unit_line) {
      std::string unitval;
      
      std::istringstream is(unit_line);
      int icol=0;
      while(is >> unitval) {
	if (unitval!=std::string(".")) {
	  this->set_unit(this->get_column_name(icol),unitval);
	}
	icol++;
      } 
      return;
    }

    /** \brief Copy all names and units from \c src to
	the current table

	If a column in the source table is already present
	in the current table, then its unit is unchanged.
	If all columns in the source are already present
	in the current table, then this function silently
	does nothing.
    */
    void copy_names_and_units(table_units &src) {
      for(size_t i=0;i<src.get_ncolumns();i++) {
	std::string cname=src.get_column_name(i);
	if (this->is_column(cname)==false) {
	  this->new_column(cname);
	  std::string cunit=get_unit(i);
	  this->set_unit(cname,cunit);
	}
      }
      return;
    }

    /** \brief Get the unit for column with index i
	
	\future Is there a way to make this function have O(1) time
	rather than searching?
    */
    std::string get_unit(size_t i) const {
      return get_unit(this->get_column_name(i));
    }

    /// Remove the unit for column \c scol
    void remove_unit(std::string scol) {
  
      uiter it=utree.find(scol);
      if (it==utree.end()) {
	O2SCL_ERR((((std::string)"Column '")+scol+
		   "' not found in table_units::get_unit().").c_str(),
		  exc_enotfound);
      }
  
      if (utree.size()==0) {
	O2SCL_ERR("No units specified in table_units::remove_unit().",
		  exc_efailed);
      }
  
      utree.erase(it);
      return;
    }
    
    /// Set the unit for column \c scol to \c unit
    void set_unit(std::string scol, std::string unit) {
  
      uiter it=utree.find(scol);
      if (it==utree.end()) {
	typename table<vec_t,double>::aiter at=this->atree.find(scol);
	if (at==this->atree.end()) {
	  O2SCL_ERR((((std::string)"Column '")+scol+
		     "' not found in table_units::set_unit().").c_str(),
		    exc_enotfound);
	}
	utree.insert(make_pair(scol,unit));
      } else {
	it->second=unit;
      }

      return;
    }

    /// Convert the units of column \c scol to \c unit
    int convert_to_unit(std::string scol, std::string unit,
			bool err_on_fail=true) {

      // Find unit entry
      uiter it=utree.find(scol);
      if (it==utree.end()) {
	if (err_on_fail) {
	  O2SCL_ERR((((std::string)"Column '")+scol+"' not found in "+
		     "table_units::convert_to_unit().").c_str(),
		    exc_enotfound);
	} else {
	  return exc_enotfound;
	}
      };
      
      // If the units are equal, just do nothing
      if (it->second==unit) return success;

      // Find column of data
      typename table<vec_t,double>::aiter at=this->atree.find(scol);
      if (at==this->atree.end()) {
	if (err_on_fail) {
	  O2SCL_ERR((((std::string)"Column '")+scol+"' not found in "+
		     "table_units::convert_to_unit().").c_str(),
		    exc_enotfound);
	} else {
	  return exc_enotfound;
	}
      };
      
      // Perform conversion
      vec_t &vec=at->second.dat;
      double conv=cup->convert(it->second,unit,1.0);

      for(size_t i=0;i<this->get_nlines();i++) {
	vec[i]*=conv;
      }

      // Set new unit entry
      it->second=unit;
      
      return success;
    }
    
    /// Return the number of columns with units
    size_t get_nunits() {
      return utree.size();
    }

    /*
    // Get the conversion factor from \c old_unit to \c new_unit
    double get_conv(std::string old_unit, std::string new_unit);
   
    // Set the convert units object
    void set_convert(convert_units &c) {
    cup=&c;
    return;
    }

    // Show the unit cache as given by \ref convert_units::print_cache()
    void show_units() {
    cup->print_cache();
    return;
    }

    // The default object for unit conversions
    convert_units def_cu;
    */
    //@}

    /// \name Virtual functions from \ref table
    //@{

    /** \brief Clear the table and the column names and units 
	(but leave constants)
    */
    virtual void clear_table() {
      this->atree.clear();
      utree.clear();
      this->alist.clear();
      this->nlines=0;
      if (this->intp_set==true) {
	delete this->si;
	this->intp_set=false;
      }
      return;
    }

    /// Delete column named \c scol
    virtual void delete_column(std::string scol) {

      // Find the tree iterator for the element we want to erase
      typename table<vec_t,double>::aiter it=this->atree.find(scol);
      if (it==this->atree.end()) {
	O2SCL_ERR((((std::string)"Column '")+scol+
		   " not found in table_units::delete_column().").c_str(),
		  exc_enotfound);
	return;
      }

      // Remove the unit for the specified column
      if (get_unit(scol).length()>0) remove_unit(scol);

      if (true) {

        int ix_match=it->second.index;
        
        for(typename table<vec_t,double>::aiter it2=this->atree.begin();
            it2!=this->atree.end();it2++) {
          if (it2->second.index>ix_match) {
            it2->second.index=it2->second.index-1;
          }
        }
        
        // Erase the elements from the list and the tree
        this->atree.erase(it);

        
        // Resize the list to the correct size
        this->alist.resize(this->atree.size());
        
      } else {
        
        // Find the corresponding list iterator
        typename table<vec_t,double>::aviter vit=this->alist.begin();
        vit+=it->second.index;
        
        // Find the last element in the list and it's corresponding table
        // entry. Change it's index to the index of the column to be
        // deleted.
        this->alist[this->alist.size()-1]->second.index=it->second.index;
        
        // Erase the elements from the list and the tree
        this->atree.erase(it);
        this->alist.erase(vit);

      }
      
      // Reset the list to reflect the proper iterators
      this->reset_list();

      if ((this->intp_colx==scol || this->intp_coly==scol) &&
          this->intp_set==true) {
        delete this->si;
        this->intp_set=false;
      }

      this->is_valid();
      
      return;
    }

    /** \brief Rename column named \c src to \c dest
	\f$ {\cal O}(C) \f$
    */
    virtual void rename_column(std::string src, std::string dest) {
      std::string unit=get_unit(src);
      table<vec_t,double>::rename_column(src,dest);
      set_unit(dest,unit);
      return;
    }

    /// Output a summary of the information stored
    virtual void summary(std::ostream *out, size_t ncol=79) const {

      if (this->constants.size()==1) {
	(*out) << "1 constant:" << std::endl;
      } else {
	(*out) << this->constants.size() << " constants:" << std::endl;
      }
      std::map<std::string,double>::const_iterator mit;
      for(mit=this->constants.begin();mit!=this->constants.end();mit++) {
	(*out) << mit->first << " " << mit->second << std::endl;
      }

      // Output number of columns and preprend column numbers
      size_t nh=this->get_ncolumns(), nh2;

      if (nh==0) {

	(*out) << "No columns." << std::endl;

      } else {

	if (nh==1) {
	  (*out) << "1 column: " << std::endl;
	} else {
	  (*out) << nh << " columns: " << std::endl;
	}
	std::string *h=new std::string[nh];
	for(size_t i=0;i<nh;i++) {
	  h[i]=szttos(i)+". "+this->get_column_name(i)+" ["+
	    get_unit(this->get_column_name(i))+"]";
	}
	
	std::vector<std::string> h2;
	// Convert to string with width 'ncol'
	screenify(nh,h,h2,ncol);
	nh2=h2.size();
  
	// Output column names
	for(size_t i=0;i<nh2;i++) {
	  (*out) << h2[i] << std::endl;
	}

	delete[] h;
  
      }
  
      if (this->get_nlines()==0) {
	(*out) << "No lines of data." << std::endl;
      } else if (this->get_nlines()==1) {
	(*out) << "One line of data." << std::endl;
      } else {
	(*out) << this->get_nlines() << " lines of data." << std::endl;
      }
  
      return;
    }
    //@}

    /// Return the type, \c "table_units".
    virtual const char *type() { return "table_units"; }

    /** \brief Copy data from column named \c src to column named \c
	dest, creating a new column if necessary \f$ {\cal O}(R
	\log(C)) \f$

	This function also sets the units of column \c dest to be the
	same as that in \c src, even if the column named \c dest
	already exists and previously had different units.
    */
    virtual void copy_column(std::string src, std::string dest) {
      if (!this->is_column(dest)) this->new_column(dest);

      typedef typename std::map<std::string,
	typename table<vec_t,double>::col,
                                std::greater<std::string> >::iterator aiter2;

      aiter2 its=this->atree.find(src);
      if (its==this->atree.end()) {
	O2SCL_ERR((((std::string)"Column '")+src+
		   " not found in table_units::copy_column().").c_str(),
		  exc_enotfound);
	return;
      }
      aiter2 itd=this->atree.find(dest);
      if (itd==this->atree.end()) {
	O2SCL_ERR((((std::string)"Destination column '")+dest+
		   " not found in table_units::copy_column().").c_str(),
		  exc_esanity);
	return;
      }
      this->set_unit(dest,this->get_unit(src));
      for(size_t i=0;i<this->nlines;i++) {
	itd->second.dat[i]=its->second.dat[i];
      }
      return;
    }
    
    /// Clear the current table and read from a generic data file
    virtual int read_generic(std::istream &fin, int verbose=0) {
	
      this->clear();
      
      double data;
      std::string line;
      
      std::string stemp;

      // Read the first line
      getline(fin,line);
      
      // Determine if there are constants
      std::vector<std::string> vsc;
      split_string_delim(line,vsc,' ');
      if (vsc.size()>1 && (vsc[1]=="constants." ||
                           vsc[1]=="constant.")) {
        size_t n_const=o2scl::stoszt(vsc[0]);
        std::string name;
        double val;
        for(size_t i=0;i<n_const;i++) {
          fin >> name >> val;
          this->add_constant(name,val);
        }
        // Read the remaining carriage return at the end of the
        // constant list
        getline(fin,line);
        // Read the next full line
        getline(fin,line);
      }

      // Determine if the interpolation type was specified
      std::vector<std::string> vsi;
      split_string_delim(line,vsi,' ');
      if (vsi.size()>1 && vsi[0]=="Interpolation:") {
        this->set_interp_type(o2scl::stoszt(vsi[1]));
        // Read the next full line
        getline(fin,line);
      }
      
      // Read first line and into list
      std::vector<std::string> onames, nnames;
      std::istringstream is2(line);
      while (is2 >> stemp) {
	onames.push_back(stemp);
	if (verbose>2) {
	  std::cout << "Read possible column name: " << stemp << std::endl;
	}
      }

      // Count number of likely numbers in the first row
      size_t n_nums=0;
      for(size_t i=0;i<onames.size();i++) {
	if (is_number(onames[i])) n_nums++;
      }

      int irow=0;

      if (n_nums==onames.size()) {

	if (verbose>0) {
	  std::cout << "First row looks like it contains numerical values." 
		    << std::endl;
	  std::cout << "Creating generic column names: ";
	}

	for(size_t i=0;i<onames.size();i++) {
	  nnames.push_back(((std::string)"c")+szttos(i+1));
	  if (verbose>0) std::cout << nnames[i] << " ";
      
	}
	if (verbose>0) std::cout << std::endl;

	// Make columns
	for(size_t i=0;i<nnames.size();i++) {
	  this->new_column(nnames[i]);
	}

	// Add first row of data
	this->set_nlines_auto(irow+1);
	for(size_t i=0;i<onames.size();i++) {
	  this->set(i,irow,o2scl::stod(onames[i]));
	}
	irow++;

      } else {

	// Ensure good column names
	for(size_t i=0;i<onames.size();i++) {
	  std::string temps=onames[i];
	  this->make_fp_varname(temps);
	  this->make_unique_name(temps,nnames);
	  nnames.push_back(temps);
	  if (temps!=onames[i] && verbose>0) {
	    std::cout << "Converted column named '" << onames[i] << "' to '" 
		      << temps << "'." << std::endl;
	  }
	}

	// Make columns
	for(size_t i=0;i<nnames.size();i++) {
	  this->new_column(nnames[i]);
	}

	// Read another line, and see if it looks like units
	std::vector<std::string> units;
	getline(fin,line);
        std::istringstream is(line);
	int num_units=0;
	while (is >> stemp) {
	  units.push_back(stemp);
	  if (stemp[0]=='[') num_units++;
	  if (verbose>2) {
	    std::cout << "Read word in second row: " << stemp << std::endl;
	  }
	}

	if (units.size()!=nnames.size()) {
	  std::cout << "Second row appears not to have same number of "
		    << "entries as the first." << std::endl;
	  std::cout << "Aborting." << std::endl;
	  return -1;
	}

	if (num_units==((int)units.size()) || num_units>2) {
	  if (verbose>2) {
	    std::cout << "Looks like second row contains units." << std::endl;
	  }
	  for(size_t i=0;i<units.size();i++) {
	    // Remove brackets
	    stemp=units[i];
	    if (stemp[0]=='[') stemp=stemp.substr(1,stemp.length()-1);
	    if (stemp[stemp.length()-1]==']') {
	      stemp=stemp.substr(0,stemp.length()-1);
	    }
	    // Set the units
	    set_unit(nnames[i],stemp);
	    if (verbose>2) {
	      std::cout << "Name,unit: " << nnames[i] << " [" << stemp << "]" 
			<< std::endl;
	    }
	  }
      
	} else {

	  // Otherwise, assume this is a row of data
	  this->set_nlines_auto(1);
	  for(size_t i=0;i<units.size();i++) {
	    this->set(i,0,o2scl::stod(units[i]));
	  }
	  irow++;

	}

      }

      // Read remaining rows
      while ((fin) >> data) {
	this->set_nlines_auto(irow+1);
	this->set(0,irow,data);
	for(size_t i=1;i<this->get_ncolumns();i++) {
	  (fin) >> data;
	  this->set(i,irow,data);
	}
	irow++;
      }

      return 0;
    }
    
    /** \brief Insert columns from a source table into the new
	table by interpolation (or extrapolation)

      This takes all of the columns in \c source, and adds them into
      the current table using interpolation, using the columns \c
      src_index and \c dest_index as the independent variable. The
      column named \c src_index is the column of the independent
      variable in \c source and the column named \c dest_index
      is the column of the independent variable in the current table.
      If \c dest_index is empty (the default) then the names in the
      two tables are taken to be the same.

      If necessary, columns are created in the current table for the
      dependent variable columns in \c source. Columns in the current
      table which do not correspond to dependent variable columns in
      \c source are left unchanged.

      If \c allow_extrap is false, then extrapolation is not allowed,
      and rows in the current table which have values of the independent
      variable which are outside the source table are unmodified. 

      If a column for a dependent variable in \c source has the
      same name as \c dest_index, then it is ignored and not inserted
      into the current table.

      If the column named \c src_index cannot be found in 
      \c source or the column names \c dest_index cannot be found
      in the current table, then the error handler is called.

      If the \c allow_extrap is false and either the minimum or
      maximum values of the column named \c src_index in the \c source
      table are not finite, then the error handler is called.
    */
    template<class vec2_t>
      void insert_table(table_units<vec2_t> &source, std::string src_index,
			bool allow_extrap=true, std::string dest_index="") {
      
      if (dest_index=="") dest_index=src_index;
      
      if (!source.is_column(src_index)) {
	O2SCL_ERR("Source indep. var. column not found.",o2scl::exc_einval);
      }
      if (!this->is_column(dest_index)) {
	O2SCL_ERR("Dest. indep. var. column not found.",o2scl::exc_einval);
      }

      // Find limits to avoid extrapolation if necessary
      double min=source.min(src_index);
      double max=source.max(src_index);
      if (allow_extrap==false) {
	if (!std::isfinite(min) || !std::isfinite(max)) {
	  O2SCL_ERR2("Minimum or maximum of source index not finite ",
		     "in table_units::insert_table().",exc_einval);
	}
      }
      
      // Create list of columns to interpolate
      std::vector<std::string> col_list;
      for(size_t i=0;i<source.get_ncolumns();i++) {
	std::string colx=source.get_column_name(i);
	if (colx!=src_index && colx!=dest_index) {
	  col_list.push_back(colx);
	}
      }
      
      // Create new columns and perform interpolation
      for(size_t i=0;i<col_list.size();i++) {
	if (!this->is_column(col_list[i])) this->new_column(col_list[i]);
	set_unit(col_list[i],source.get_unit(col_list[i]));
	for(size_t j=0;j<this->get_nlines();j++) {
	  double val=this->get(dest_index,j);
	  if (allow_extrap || (val>=min && val<=max)) {
	    this->set(col_list[i],j,source.interp(src_index,val,col_list[i]));
	  }
	}
      }
      
      return;
    }
    
    // ---------
    // Allow HDF I/O functions to access table_units data
    friend void o2scl_hdf::hdf_output
      (o2scl_hdf::hdf_file &hf, table_units<> &t, std::string name);
    
    template<class vecf_t> friend void o2scl_hdf::hdf_input
      (o2scl_hdf::hdf_file &hf, table_units<vecf_t> &t, std::string name);

    friend void o2scl_hdf::hdf_output_data
      (o2scl_hdf::hdf_file &hf, table_units<> &t);
    
    template<class vecf_t> friend void o2scl_hdf::hdf_input_data
      (o2scl_hdf::hdf_file &hf, table_units<vecf_t> &t);

    // ---------

#ifndef DOXYGEN_INTERNAL

  protected:
    

    /// The pointer to the convert units object
    convert_units<double> *cup;

    /// \name Unit map iterator types
    //@{
    typedef std::map<std::string,std::string,
      std::greater<std::string> >::iterator uiter;
    typedef std::map<std::string,std::string,
      std::greater<std::string> >::const_iterator uciter;
    //@}
  
    /// Unit map
    std::map<std::string,std::string,std::greater<std::string> > utree;
  
#endif

  };

#ifdef O2SCL_MPI

  /** \brief Send a \ref o2scl::table_units object to
      MPI rank \c dest_rank
  */
  template<class vec_t>
    void o2scl_table_mpi_send(o2scl::table_units<vec_t> &t, size_t dest_rank) {

    int tag;
    int ibuffer;
    std::vector<double> dbuffer;
    std::string cbuffer;

    // --------------------------------------------------------------
    // Constant names and values
  
    for(size_t i=0;i<t.get_nconsts();i++) {
      std::string name;
      double val;
      t.get_constant(i,name,val);
      cbuffer+=name;
      if (i<t.get_nconsts()-1) {
	cbuffer+=' ';
      }
      dbuffer.push_back(val);
    }

    ibuffer=cbuffer.length();
  
    tag=0;
    MPI_Send(&ibuffer,1,MPI_INT,dest_rank,
	     tag,MPI_COMM_WORLD);
    tag=1;
    MPI_Send(&(cbuffer[0]),cbuffer.length(),MPI_CHAR,dest_rank,
	     tag,MPI_COMM_WORLD);
    tag=2;
    MPI_Send(&(dbuffer[0]),dbuffer.size(),MPI_DOUBLE,dest_rank,
	     tag,MPI_COMM_WORLD);
    cbuffer.clear();
    dbuffer.clear();

    // --------------------------------------------------------------
    // Column names
  
    for(size_t i=0;i<t.get_ncolumns();i++) {
      std::string name=t.get_column_name(i);
      cbuffer+=name;
      if (i<t.get_ncolumns()-1) {
	cbuffer+=' ';
      }
    }

    ibuffer=cbuffer.length();
    tag=3;
    MPI_Send(&ibuffer,1,MPI_INT,dest_rank,
	     tag,MPI_COMM_WORLD);
    tag=4;
    MPI_Send(&(cbuffer[0]),cbuffer.length(),MPI_CHAR,dest_rank,
	     tag,MPI_COMM_WORLD);
    cbuffer.clear();
    
    // --------------------------------------------------------------
    // Column units
  
    for(size_t i=0;i<t.get_ncolumns();i++) {
      std::string unit=t.get_unit(t.get_column_name(i));
      cbuffer+=unit;
      if (i<t.get_ncolumns()-1) {
	cbuffer+=' ';
      }
    }

    ibuffer=cbuffer.length();
    tag=5;
    MPI_Send(&ibuffer,1,MPI_INT,dest_rank,
	     tag,MPI_COMM_WORLD);
    tag=6;
    MPI_Send(&(cbuffer[0]),cbuffer.length(),MPI_CHAR,dest_rank,
	     tag,MPI_COMM_WORLD);
    cbuffer.clear();

    // --------------------------------------------------------------
    // Interpolation type

    tag=7;
    ibuffer=t.get_interp_type();
    MPI_Send(&ibuffer,1,MPI_INT,dest_rank,
	     tag,MPI_COMM_WORLD);
    
    // --------------------------------------------------------------
    // Column data

    ibuffer=t.get_nlines();
    tag=8;
    MPI_Send(&ibuffer,1,MPI_INT,dest_rank,
	     tag,MPI_COMM_WORLD);
  
    for(size_t i=0;i<t.get_ncolumns();i++) {
      tag++;
      MPI_Send(&(t[i][0]),t.get_nlines(),MPI_DOUBLE,dest_rank,
	       tag,MPI_COMM_WORLD);
    }

    return;
  }

  /** \brief Receive a \ref o2scl::table_units object from
      MPI rank \c src_rank
  */
  template<class vec_t>
    void o2scl_table_mpi_recv(size_t src_rank,
			      o2scl::table_units<vec_t> &t) {
    int tag;
    int ibuffer;
    // Not std::string because we need to guarantee contiguous storage?!
    std::vector<char> cbuffer;
    std::vector<double> dbuffer;

    std::vector<std::string> names;
    std::string stemp;

    // --------------------------------------------------------------
    // Constants

    // Names
    tag=0;
    MPI_Recv(&ibuffer,1,MPI_INT,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    cbuffer.resize(ibuffer);
    tag=1;
    MPI_Recv(&(cbuffer[0]),ibuffer,MPI_CHAR,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    // Parse into std::vector<string>
    for(size_t i=0;i<cbuffer.size();i++) {
      if (cbuffer[i]!=' ') {
	stemp+=cbuffer[i];
      } else {
	if (stemp.size()>0) {
	  names.push_back(stemp);
	  stemp="";
	}
      }
    }
    if (stemp.size()>0) {
      names.push_back(stemp);
    }
    stemp="";

    // Load values
    dbuffer.resize(names.size());
  
    tag=2;
    MPI_Recv(&(dbuffer[0]),dbuffer.size(),MPI_DOUBLE,
	     src_rank,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Set constants
    for(size_t i=0;i<names.size();i++) {
      t.add_constant(names[i],dbuffer[i]);
    }
    names.clear();
    dbuffer.clear();
  
    // --------------------------------------------------------------
    // Column names
  
    tag=3;
    MPI_Recv(&ibuffer,1,MPI_INT,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    cbuffer.resize(ibuffer);
    tag=4;
    MPI_Recv(&(cbuffer[0]),ibuffer,MPI_CHAR,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    // Parse into std::vector<string>
    for(size_t i=0;i<cbuffer.size();i++) {
      if (cbuffer[i]!=' ') {
	stemp+=cbuffer[i];
      } else {
	if (stemp.size()>0) {
	  names.push_back(stemp);
	  stemp="";
	}
      }
    }
    if (stemp.size()>0) {
      names.push_back(stemp);
    }
    stemp="";

    for(size_t i=0;i<names.size();i++) {
      t.new_column(names[i]);
    }
    names.clear();

    // --------------------------------------------------------------
    // Column units
  
    tag=5;
    MPI_Recv(&ibuffer,1,MPI_INT,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    cbuffer.resize(ibuffer);
    tag=6;
    MPI_Recv(&(cbuffer[0]),ibuffer,MPI_CHAR,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    // Parse into std::vector<string>
    for(size_t i=0;i<cbuffer.size();i++) {
      if (cbuffer[i]!=' ') {
	stemp+=cbuffer[i];
      } else {
	if (stemp.size()>0) {
	  names.push_back(stemp);
	  stemp="";
	}
      }
    }
    if (stemp.size()>0) {
      names.push_back(stemp);
    }
    stemp="";

    for(size_t i=0;i<names.size();i++) {
      t.set_unit(t.get_column_name(i),names[i]);
    }
    names.clear();

    // --------------------------------------------------------------
    // Column data

    tag=7;
    MPI_Recv(&ibuffer,1,MPI_INT,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    t.set_interp_type(ibuffer);

    // --------------------------------------------------------------
    // Column data

    tag=8;
    MPI_Recv(&ibuffer,1,MPI_INT,src_rank,tag,MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    t.set_nlines(ibuffer);

    for(size_t i=0;i<t.get_ncolumns();i++) {
      tag++;
      std::vector<double> v(t.get_maxlines());
      MPI_Recv(&(v[0]),ibuffer,MPI_DOUBLE,src_rank,
	       tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      t.swap_column_data(t.get_column_name(i),v);
    }

    return;
  }
  
#endif  
  
}

#endif
