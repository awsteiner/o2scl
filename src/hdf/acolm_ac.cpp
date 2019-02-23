/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int acol_manager::comm_assign(std::vector<std::string> &sv, bool itive_com) {
  if (type!="table" && type!="table3d") {
    cerr << "No table/table3d object to add a constant to." << endl;
    return exc_efailed;
  }

  if (sv.size()==2) {
    if (verbose>0) {
      cout << "Removing constant named '" << sv[1] << "'." << endl;
    }
    if (type=="table3d") {
      table3d_obj.remove_constant(sv[1]);
    } else {
      table_obj.remove_constant(sv[1]);
    }
    return 0;
  }

  vector<string> pr, in;
  pr.push_back("Name of constant");
  pr.push_back("Value");
  int ret=get_input(sv,pr,in,"assign",itive_com);
  if (ret!=0) return ret;
  
  if (type=="table3d") {
    table3d_obj.add_constant(sv[1],function_to_double(sv[2]));
  } else {
    table_obj.add_constant(sv[1],function_to_double(sv[2]));
  }

  return ret;
}

int acol_manager::comm_autocorr(std::vector<std::string> &sv,
				bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table with columns to compute "
	   << "autocorrelations with." << endl;
      return exc_efailed;
    }

    vector<string> in, pr;
    if (sv.size()<4) {
      pr.push_back("Enter output column for autocorrelations");
      pr.push_back("Enter output column for 5*tau/m");
      pr.push_back("Enter vector specification for data");
      int ret=get_input(sv,pr,in,"autocorr",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t i=1;i<sv.size();i++) in.push_back(sv[i]);
    }

    vector<vector<double> > v_all, ac_all, ftom_all;
    size_t max_ftom_size=0;
    
    for(size_t ix=2;ix<in.size();ix++) {
      
      vector<double> v, ac, ftom;
      
      if (in[ix].find(':')==std::string::npos &&
	  table_obj.is_column(in[ix])==false) {
	cerr << "Could not find column named '" << in[ix] << "'." << endl;
	return exc_efailed;
      }
      
      if (!table_obj.is_column(in[0])) {
	table_obj.new_column(in[0]);
      }
      if (!table_obj.is_column(in[1])) {
	table_obj.new_column(in[1]);
      }
      
      if (in[ix].find(':')==std::string::npos) {
	v.resize(table_obj.get_nlines());
	for(size_t i=0;i<table_obj.get_nlines();i++) {
	  v[i]=table_obj.get(in[ix],i);
	}
      } else {
	int vs_ret=o2scl_hdf::vector_spec(in[ix],v,verbose,false);
	if (vs_ret!=0) {
	  cout << "Vector specification failed." << endl;
	  return 1;
	}
      }
      
      // Compute autocorrelation length and sample size
      vector_autocorr_vector(v,ac);
      size_t len=vector_autocorr_tau(ac,ftom);
      if (len>0) {
	cout << "Autocorrelation length: " << len << " sample size: "
	     << table_obj.get_nlines()/len << endl;
      } else {
	cout << "Autocorrelation length determination failed." << endl;
      }

      if (ftom.size()>max_ftom_size) {
	max_ftom_size=ftom.size();
      }
      v_all.push_back(v);
      ac_all.push_back(ac);
      ftom_all.push_back(ftom);
    }

    vector<double> ac_avg(max_ftom_size), ftom_avg(max_ftom_size);
    for(size_t i=0;i<max_ftom_size;i++) {
      size_t n=0;
      ac_avg[i]=0.0;
      ftom_avg[i]=0.0;
      for(size_t j=0;j<ac_all.size();j++) {
	if (i<ac_all[j].size() && i<ftom_all[j].size()) {
	  n++;
	  ac_avg[i]+=ac_all[j][i];
	  ftom_avg[i]+=ftom_all[j][i];
	}
      }
      if (n==0) {
	cerr << "Failed to find any data in 'autocorr'." << endl;
	return 1;
      }
      ac_avg[i]/=((double)n);
      ftom_avg[i]/=((double)n);
    }
    
    // Add autocorrelation and ftom data to table, replacing the
    // values with zero when we reach the end of the vectors given by
    // vector_autocorr_tau() .
    for(size_t i=0;i<table_obj.get_nlines();i++) {
      if (i<ac_avg.size()) {
	table_obj.set(in[0],i,ac_avg[i]);
      } else {
	table_obj.set(in[0],i,0.0);
      }
      if (i<ftom_avg.size()) {
	table_obj.set(in[1],i,ftom_avg[i]);
      } else {
	table_obj.set(in[1],i,0.0);
      }
    }

  } else if (type=="double[]") {

    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(ac_vec,ftom);
    if (len>0) {
      cout << "Autocorrelation length: " << len << " sample size: "
	   << doublev_obj.size()/len << endl;
    } else {
      cout << "Autocorrelation length determination failed." << endl;
    }

    doublev_obj=ac_vec;

  } else if (type=="int[]") {

    vector_copy(intv_obj,doublev_obj);
    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(ac_vec,ftom);
    if (len>0) {
      cout << "Autocorrelation length: " << len << " sample size: "
	   << doublev_obj.size()/len << endl;
    } else {
      cout << "Autocorrelation length determination failed." << endl;
    }

    command_del();
    clear_obj();
    doublev_obj=ac_vec;
    command_add("double[]");
    type="double[]";
    
  } else if (type=="size_t[]") {
    
    vector_copy(size_tv_obj,doublev_obj);
    vector<double> ac_vec, ftom;
    vector_autocorr_vector(doublev_obj,ac_vec);
    size_t len=vector_autocorr_tau(ac_vec,ftom);
    if (len>0) {
      cout << "Autocorrelation length: " << len << " sample size: "
	   << doublev_obj.size()/len << endl;
    } else {
      cout << "Autocorrelation length determination failed." << endl;
    }

    command_del();
    clear_obj();
    doublev_obj=ac_vec;
    command_add("double[]");
    type="double[]";
    
  } else {
    
    vector<string> in, pr;
    if (sv.size()<2) {
      pr.push_back("Enter vector specification for data");
      int ret=get_input(sv,pr,in,"autocorr",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t i=1;i<sv.size();i++) in.push_back(sv[i]);
    }

    vector<vector<double> > v_all, ac_all, ftom_all;
    size_t max_ftom_size=0;
    
    for(size_t ix=0;ix<in.size();ix++) {
      
      vector<double> v, ac, ftom;
      
      if (in[0].find(':')==std::string::npos) {
	v.resize(table_obj.get_nlines());
	for(size_t i=0;i<table_obj.get_nlines();i++) {
	  v[i]=table_obj.get(in[0],i);
	}
      } else {
	int vs_ret=o2scl_hdf::vector_spec(in[0],v,verbose,false);
	if (vs_ret!=0) {
	  cout << "Vector specification failed." << endl;
	  return 1;
	}
      }
      
      // Compute autocorrelation length and sample size
      vector_autocorr_vector(v,ac);
      size_t len=vector_autocorr_tau(ac,ftom);
      if (len>0) {
	cout << "Autocorrelation length: " << len << " sample size: "
	     << table_obj.get_nlines()/len << endl;
      } else {
	cout << "Autocorrelation length determination failed." << endl;
      }

      if (ftom.size()>max_ftom_size) {
	max_ftom_size=ftom.size();
      }
      v_all.push_back(v);
      ac_all.push_back(ac);
      ftom_all.push_back(ftom);
    }
    
    command_del();
    clear_obj();
    
    doublev_obj.resize(max_ftom_size);
    for(size_t i=0;i<max_ftom_size;i++) {
      size_t n=0;
      doublev_obj[i]=0.0;
      for(size_t j=0;j<ac_all.size();j++) {
	if (i<ac_all[j].size() && i<ftom_all[j].size()) {
	  n++;
	  doublev_obj[i]+=ac_all[j][i];
	}
      }
      if (n==0) {
	cerr << "Failed to find any data in 'autocorr'." << endl;
	return 2;
      }
      doublev_obj[i]/=((double)n);
    }
    
    command_add("double[]");
    type="double[]";
  }    
  
  return 0;
}

#ifdef O2SCL_NEVER_DEFINED
int acol_manager::comm_comment(std::vector<std::string> &sv, 
			       bool itive_com) {

  if (sv.size()==2) {
    hdf_file hf;
    hf.open(sv[1]);
    std::string def, s;
    hf.gets_def("comment",def,s);
    if (s==def) {
      cout << "No comment in file " << sv[1] << endl;
    } else {
      cout << "Comment in file " << sv[1] << " :" << endl;
      cout << s << endl;
    }
    hf.close();
    return 0;
  }
  
  hdf_file hf;
  // Make sure to open with write access
  hf.open(sv[1],1);
  
  // If it's already present as a fixed length string,
  // then we need to double check
  std::string def, s;
  int iret=hf.gets_def_fixed("comment",def,s);
  if (s!=def) {
    if (iret==1) {
      size_t len=s.length();
      if (sv[2].length()>len) {
	cerr << "Size of new comment (" << sv[2].length()
	     << ") longer than size of current "
	     << "fixed length string " << len << "." << endl;
	hf.close();
	return 1;
      } else {
	while (sv[2].length()<len) sv[2]+=' ';
      }
      hf.sets_fixed("comment",sv[2]);
    } else {
      hf.sets("comment",sv[2]);
    }
  } else {
    // String is not present so just set
    hf.sets("comment",sv[2]);
  }
  cout << "Set comment in file " << sv[1] << " to " << endl;
  cout << sv[2] << endl;
  hf.close();
  return 0;
}
#endif

int acol_manager::comm_calc(std::vector<std::string> &sv, bool itive_com) {

  std::string i1;
  if (sv.size()>1) {
    i1=sv[1];
  } else if (itive_com) {
    i1=cl->cli_gets("Enter expression to compute (or blank to stop): ");
    if (i1.length()==0) {
      if (verbose>0) cout << "Command 'calc' cancelled." << endl;
      return 0;
    }
  } else {
    cerr << "No expression to compute in 'calc'." << endl;
    return exc_efailed;
  }
  double d=o2scl::function_to_double(i1);
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  if (verbose>0) cout << "Result: ";
  cout << d << endl;
  return 0;
}

int acol_manager::comm_cat(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cerr << "Not enough arguments to cat." << endl;
    return exc_efailed;
  }
  string file2=sv[1];
  
  if (type=="table3d") {

    if (type!="table3d") {
      cerr << "No table3d to add to in command 'cat'." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Read the new table3d
    
    table3d tab2;

    hdf_file hf;
    std::string name2;
    if (sv.size()>=3) name2=sv[2];

    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    // ---------------------------------------------------------------------
    // Copy constants from the new table3d

    for(size_t i=0;i<tab2.get_nconsts();i++) {
      string tnam;
      double tval;
      tab2.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      table3d_obj.add_constant(tnam,tval);
    }

    // ---------------------------------------------------------------------
    // Copy slices over if not already present in the current table3d

    const ubvector &xg=tab2.get_x_data();
    const ubvector &yg=tab2.get_y_data();

    for(size_t k=0;k<tab2.get_nslices();k++) {
      std::string sl_name=tab2.get_slice_name(k);
      size_t slix;
      if (!table3d_obj.is_slice(sl_name,slix)) {
	table3d_obj.new_slice(sl_name);
	for(size_t i=0;i<tab2.get_nx();i++) {
	  for(size_t j=0;j<tab2.get_ny();j++) {
	    double x=xg[i];
	    double y=yg[j];
	    table3d_obj.set_val(x,y,sl_name,tab2.get(i,j,sl_name));
	  }
	}
      }
    }

  } else if (type=="table") {

    if (table_obj.get_nlines()==0) {
      cerr << "No table to add to in command 'cat'." << endl;
      return exc_efailed;
    }

    // ---------------------------------------------------------------------
    // Read the new table 

    table_units<> tab2;

    hdf_file hf;
    std::string name2;
    if (sv.size()>=3) name2=sv[2];

    int hfret=hf.open(file2,false,false);
    if (hfret!=0) {
      cerr << "Failed to read file named " << file2 << endl;
      return exc_efailed;
    }
    hdf_input(hf,tab2,name2);
    hf.close();

    // ---------------------------------------------------------------------
    // Copy constants from the new table

    for(size_t i=0;i<tab2.get_nconsts();i++) {
      string tnam;
      double tval;
      tab2.get_constant(i,tnam,tval);
      if (verbose>2) {
	cout << "Adding constant " << tnam << " = " << tval << endl;
      }
      table_obj.add_constant(tnam,tval);
    }

    // ---------------------------------------------------------------------

    size_t n1=table_obj.get_nlines();
    size_t n2=tab2.get_nlines();
    table_obj.set_nlines(n1+n2);
    for(size_t j=0;j<tab2.get_ncolumns();j++) {
      std::string col_name=tab2.get_column_name(j);
      if (!table_obj.is_column(col_name)) {
	table_obj.new_column(col_name);
	for(size_t i=0;i<n1+n2;i++) table_obj.set(col_name,i,0.0);
      }
      for(size_t i=0;i<n2;i++) {
	table_obj.set(col_name,i+n1,tab2.get(col_name,i));
      }
    }

    if (verbose>0) {
      cout << "Table with " << n1 << " lines now has "
	   << n1+n2 << " lines." << endl;
    }
    
  } else {

    cerr << "Cannot 'cat' with object of type " << type << endl;
    return exc_efailed;
    
  }
  
  return 0;
}

int acol_manager::comm_commands(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()==2) {
    cout << "Commands argument: " << sv[1] << endl;
    string temp_type=sv[1];
    string cur_type=type;

    command_del();
    command_add(temp_type);
    
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
    int ret=cl->comm_option_commands(sv,itive_com);

    command_del();
    command_add(cur_type);
    return ret;
  }
  return cl->comm_option_commands(sv,itive_com);
}

int acol_manager::comm_contours(std::vector<std::string> &sv, bool itive_com) {

  if (type!="table3d" && type!="hist_2d") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  bool frac_mode=false;
  
  if (sv.size()>=2 && sv[1]=="frac") {
    cout << "Fraction mode is true." << endl;
    frac_mode=true;
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
  } else {
    cout << "Fraction mode is false." << endl;
  }
  if (sv.size()<2 && itive_com) {
    string temp=((string)"Enter \"frac\" for fractions of total sum and ")
      +"\"abs\" for absolute scale", i1;
    int ret=get_input_one(sv,temp,i1,"contours",itive_com);
    if (ret!=0) return ret;
    if (i1=="frac") frac_mode=true;
  }
    
  std::string svalue, file, name="contours";

  if (type=="table3d") {

    std::string slice;

    if (sv.size()<3) {
      svalue=cl->cli_gets("Contour value (or blank to cancel): ");
      slice=cl->cli_gets("Slice (or blank to cancel): ");
      if (svalue.length()==0) return 1;
      file=cl->cli_gets("Filename (or blank to keep): ");
      if (file.length()>0) {
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==3) {
      svalue=sv[1];
      slice=sv[2];
    } else if (sv.size()==4) {
      svalue=sv[1];
      slice=sv[2];
      file=sv[3];
    } else {
      svalue=sv[1];
      slice=sv[2];
      file=sv[3];
      name=sv[4];
    }

    ubvector levs(1);
    levs[0]=o2scl::function_to_double(svalue);
    size_t nlev=1;
    
    if (frac_mode) {
      cout << "Fraction mode not implemented with table3d objects." << endl;
    } else {
      if (file.length()>0) {
	std::vector<contour_line> clines;
	table3d_obj.slice_contours(slice,1,levs,clines);
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,clines,name);
	hf.close();
      } else {
	table3d_obj.slice_contours(slice,1,levs,cont_obj);
	command_del();
	clear_obj();
	command_add("vector<contour_line>");
	type="vector<contour_line>";
      }
    }
    
  } else if (type=="hist_2d") {

    if (sv.size()<2) {
      svalue=cl->cli_gets("Contour value (or blank to cancel): ");
      if (svalue.length()==0) return 1;
      file=cl->cli_gets("Filename (or blank to keep): ");
      if (file.length()>0) {
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==2) {
      svalue=sv[1];
    } else if (sv.size()==3) {
      svalue=sv[1];
      file=sv[2];
    } else {
      svalue=sv[1];
      file=sv[2];
      name=sv[3];
    }

    
    ubvector levs(1);
    levs[0]=o2scl::function_to_double(svalue);
    size_t nlev=1;

    if (frac_mode) {

      // Get references to the histogram data
      size_t nx=hist_2d_obj.size_x();
      size_t ny=hist_2d_obj.size_y();
      const ubmatrix &m=hist_2d_obj.get_wgts();
      const ubvector &xbins=hist_2d_obj.get_x_bins();
      const ubvector &ybins=hist_2d_obj.get_y_bins();

      // Compute the total integral and the target fraction
      double min, max;
      o2scl::matrix_minmax(m,min,max);
      double sum=hist_2d_obj.integ_wgts();
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  sum-=min*(xbins[i+1]-xbins[i])*(ybins[j+1]-ybins[j]);
	}
      }
      double target=levs[0]*sum;
      if (verbose>1) {
	cout << "sum,target: " << sum << " " << target << endl;
      }

      // Setup the vectors to interpolate the target integral
      uniform_grid_end<double> ug(min,max,100);
      ubvector integx, integy;
      ug.vector(integx);
      size_t N=integx.size();
      if (verbose>1) {
	cout << "N integx[0] integx[1]: " << N << " "
	     << integx[0] << " " << integx[1] << endl;
      }
      integy.resize(N);

      // Fill the interpolation vectors
      for(size_t k=0;k<N;k++) {
	integy[k]=0.0;
	for(size_t i=0;i<nx;i++) {
	  for(size_t j=0;j<ny;j++) {
	    if (m(i,j)>integx[k]) {
	      integy[k]+=(m(i,j)-min)*(xbins[i+1]-xbins[i])*
		(ybins[j+1]-ybins[j]);
	    }
	  }
	}
	if (verbose>1) {
	  cout << k << " " << integx[k] << " " << integy[k] << endl;
	}
      }

      // Perform the interpolation
      bool found=false;
      double level=0.0;
      for(size_t k=0;k<N-1;k++) {
	if (integy[k]>target && integy[k+1]<target) {
	  found=true;
	  level=integx[k]+(integx[k+1]-integx[k])*(target-integy[k])/
	    (integy[k+1]-integy[k]);
	}
      }
      
      // Return if the interpolation failed
      if (found==false) {
	cerr << "Failed to find a level matching requested fraction."
	     << endl;
	return 2;
      }
      
      if (verbose>1) {
	cout << "Found: " << level << endl;
      }
      // Set level from interpolated value
      levs[0]=level;
      
    }

    contour co;
    co.set_levels(nlev,levs);
    
    ubvector xreps(hist_2d_obj.size_x());
    for (size_t i=0;i<hist_2d_obj.size_x();i++) {
      xreps[i]=hist_2d_obj.get_x_rep_i(i);
    }
    ubvector yreps(hist_2d_obj.size_y());
    for (size_t i=0;i<hist_2d_obj.size_y();i++) {
      yreps[i]=hist_2d_obj.get_y_rep_i(i);
    }
    co.set_data(hist_2d_obj.size_x(),hist_2d_obj.size_y(),xreps,yreps,
		hist_2d_obj.get_wgts());

    if (file.length()>0) {
      std::vector<contour_line> clines;
      co.calc_contours(clines);
      
      hdf_file hf;
      hf.open_or_create(file);
      hdf_output(hf,clines,name);
      hf.close();
    } else {
      command_del();
      clear_obj();
      co.calc_contours(cont_obj);
      command_add("vector<contour_line>");
      type="vector<contour_line>";
    }
    
  }
  
  return 0;
}

int acol_manager::comm_convert_unit
(std::vector<std::string> &sv, bool itive_com) {
  
  if (type!="table") {
    cerr << "Not implemented for " << type << " objects." << endl;
    return exc_efailed;
  }
  
  if (table_obj.get_nlines()==0) {
    cerr << "No table to convert units in." << endl;
    return exc_efailed;
  }

  vector<string> in, pr;
  pr.push_back("Column in which to convert units");
  pr.push_back("New unit");
  int ret=get_input(sv,pr,in,"convert-unit",itive_com);
  if (ret!=0) return ret;
  
  if (table_obj.is_column(in[0])==false) {
    cerr << "Could not find column named '" << in[0] << "'." << endl;
    return exc_efailed;
  }

  if (unit_fname.length()>0) {
    cng.units_cmd_string=((string)"units -f ")+unit_fname;
  }
  ret=table_obj.convert_to_unit(in[0],in[1],false);
  if (ret!=0) {
    cerr << "Could not find column or column does not have unit." << endl;
  }

  return 0;
}

int acol_manager::comm_create(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype, tval;

  // Delete previous object
  command_del();
  clear_obj();
  
  int ret=get_input_one(sv,"Enter type of object to create",ctype,"create",
			itive_com);
  if (ret!=0) return ret;

  vector<string> sv2=sv;
  vector<string>::iterator it=sv2.begin();
  sv2.erase(it+1);
  
  if (ctype=="int") {

    int ret=get_input_one(sv2,"Enter integer",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    int_obj=o2scl::stoi(tval);
    type="int";
    command_add("int");
    obj_name="int";
    
  } else if (ctype=="size_t") {

    int ret=get_input_one(sv2,"Enter size_t",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    size_t_obj=o2scl::stoszt(tval);
    type="size_t";
    command_add("size_t");
    obj_name="size_t";
    
  } else if (ctype=="char") {

    int ret=get_input_one(sv2,"Enter char",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    char_obj=tval[0];
    type="char";
    command_add("char");
    obj_name="char";
    
  } else if (ctype=="double") {

    int ret=get_input_one(sv2,"Enter double",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    double_obj=o2scl::function_to_double(tval);
    type="double";
    command_add("double");
    obj_name="double";

  } else if (ctype=="string") {

    int ret=get_input_one(sv2,"Enter string",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    string_obj=tval;
    type="string";
    command_add("string");
    obj_name="string";
    
  } else if (ctype=="double[]") {

    std::string in1;
    int ret1=get_input_one(sv2,"Size of vector or full vector specification",
			   in1,"create",itive_com);
    if (ret1!=0) return ret1;

    if (in1.find(':')==std::string::npos) {

      vector<string> sv3=sv2;
      vector<string>::iterator it=sv3.begin();
      sv3.erase(it+1);
      
      std::string in2;

      int ret2=get_input_one(sv3,"Function of i (starting with zero)",
			     in2,"create",itive_com);
      if (ret2!=0) return ret2;
      
      calculator calc;
      std::map<std::string,double> vars;
      std::map<std::string,double>::const_iterator mit;
      size_t nn=o2scl::stoszt(in1);
      doublev_obj.clear();
      calc.compile(in2.c_str(),&vars);
      for(size_t i=0;i<nn;i++) {
	vars["i"]=((double)i);
	doublev_obj.push_back(calc.eval(&vars));
      }
      
    } else {

      int ret2=vector_spec(in1,doublev_obj,2,false);
      if (ret2!=0) {
	cerr << "Function vector_spec() failed." << endl;
	return ret2;
      }
      
    }
    
    command_add("double[]");
    type="double[]";
    
  } else if (ctype=="int[]") {
    
    vector<string> in, pr;
    pr.push_back("Size");
    pr.push_back("Function of i (starting with zero)");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    calculator calc;
    std::map<std::string,double> vars;
    std::map<std::string,double>::const_iterator mit;
    size_t nn=o2scl::stoszt(in[0]);
    intv_obj.clear();
    calc.compile(in[1].c_str(),&vars);
    for(size_t i=0;i<nn;i++) {
      vars["i"]=((double)i);
      intv_obj.push_back(((int)(calc.eval(&vars))));
    }
    command_add("int[]");
    type="int[]";
    
  } else if (ctype=="size_t[]") {
    
    vector<string> in, pr;
    pr.push_back("Size");
    pr.push_back("Function of i (starting with zero)");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    calculator calc;
    std::map<std::string,double> vars;
    std::map<std::string,double>::const_iterator mit;
    size_t nn=o2scl::stoszt(in[0]);
    size_tv_obj.clear();
    calc.compile(in[1].c_str(),&vars);
    for(size_t i=0;i<nn;i++) {
      vars["i"]=((double)i);
      size_tv_obj.push_back(((size_t)(calc.eval(&vars))));
    }
    command_add("size_t[]");
    type="size_t[]";
    
  } else if (ctype=="table") {
    
    vector<string> in, pr;
    pr.push_back("Name of new column");
    pr.push_back("Value for first row");
    pr.push_back("Maximum value");
    pr.push_back("Increment");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;
    
    double d2=function_to_double(in[1]);
    double d3=function_to_double(in[2]);
    double d4=function_to_double(in[3]);
    d3+=d4/1.0e4;
    int cnl=((int)((d3-d2)/d4))+1;
    
    table_obj.clear();
    table_obj.line_of_names(in[0]);
    table_obj.set_nlines(cnl);
    
    for(int li=0;li<cnl;li++) {
      table_obj.set(in[0],li,d2+((double)li)*d4);
    }
    command_add("table");
    type="table";
    
  } else if (ctype=="tensor_grid") {

    std::string i1;
    int ret=get_input_one(sv2,"Enter rank",i1,"create",itive_com);
    if (ret!=0) return ret;
    size_t rank=o2scl::stoszt(sv2[1]);

    if (sv2.size()<2+rank) {
      vector<string> pr, in;
      for(size_t k=0;k<rank;k++) {
	pr.push_back(((std::string)"Enter size for rank ")+
		     o2scl::szttos(rank));
      }
      int ret=get_input(sv2,pr,in,"create",itive_com);
      if (ret!=0) return ret;
    }
    
    vector<size_t> sarr(rank);
    for(size_t k=0;k<rank;k++) {
      sarr[k]=o2scl::stoszt(sv2[2+k]);
    }

    tensor_grid_obj.resize(rank,sarr);
    command_add("tensor_grid");
    type="tensor_grid";

  } else if (ctype=="table3d") {
    
    vector<string> in;
    vector<string> pr=
      {"Enter x-axis name (or blank to stop): ",
       "Enter x-axis lower limit (or blank to stop): ",
       "Enter x-axis upper limit (or blank to stop): ",
       "Enter x-axis step size (or blank to stop): ",
       "Enter y-axis name (or blank to stop): ",
       "Enter y-axis lower limit (or blank to stop): ",
       "Enter y-axis upper limit (or blank to stop): ",
       "Enter y-axis step size (or blank to stop): ",
       "Enter slice name (or blank to stop): ",
       "Enter slice function (or blank to stop): "};
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;
    
    std::string xname=in[0];
    double x0=function_to_double(in[1]);
    double x1=function_to_double(in[2]);
    double dx=function_to_double(in[3]);
    uniform_grid_end<double> ugx(x0,x1,((size_t)(x1-x0)/(dx*(1.0-1.0e-14))));
    
    std::string yname=in[4];
    double y0=function_to_double(in[5]);
    double y1=function_to_double(in[6]);
    double dy=function_to_double(in[7]);
    uniform_grid_end<double> ugy(y0,y1,((size_t)(y1-y0)/(dy*(1.0-1.0e-14))));
    
    std::string zname=in[8];
    std::string zfunc=in[9];
    
    table3d_obj.set_xy(xname,ugx,yname,ugy);
    
    table3d_obj.function_slice(zfunc,zname);
    
    command_add("table3d");
    type="table3d";

  } else {

    cerr << "Cannot create object of type " << ctype << endl;
    return 1;
      
  }

  return 0;
}

