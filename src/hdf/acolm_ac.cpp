/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

  double d;
  int retx=function_to_double_nothrow(sv[2],d);
  if (retx!=0) {
    cerr << "Converting " << sv[2] << " to value failed." << endl;
    return 1;
  }
  if (type=="table3d") {
    table3d_obj.add_constant(sv[1],d);
  } else {
    table_obj.add_constant(sv[1],d);
  }

  return ret;
}

int acol_manager::comm_autocorr(std::vector<std::string> &sv,
				bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "Table has no lines of data to compute "
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

    vector<vector<double> > ac_all;
    size_t max_ac_size=0;
    
    for(size_t ix=2;ix<in.size();ix++) {
      
      vector<double> v, ac, ftom;

      // If the argument is not a vector specification, then look
      // for the column in the table
      if (in[ix].find(':')==std::string::npos &&
	  table_obj.is_column(in[ix])==false) {
	cerr << "Could not find column named '" << in[ix] << "'." << endl;
	return exc_efailed;
      }

      // Determine vector from table column (requires copy) or
      // vector specification
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

      if (ix==2 || ac.size()>max_ac_size) {
	max_ac_size=ac.size();
      }
      ac_all.push_back(ac);
    }

    if (max_ac_size==0) {
      cerr << "Failed to find any data in 'autocorr'." << endl;
      return 1;
    }
    vector<double> ac_avg(max_ac_size);
    for(size_t i=0;i<max_ac_size;i++) {
      size_t n=0;
      ac_avg[i]=0.0;
      for(size_t j=0;j<ac_all.size();j++) {
	if (i<ac_all[j].size()) {
	  n++;
	  ac_avg[i]+=ac_all[j][i];
	}
      }
      ac_avg[i]/=((double)n);
    }
    
    // Now report autocorrelation length and sample size from
    // averaged result
    vector<double> ftom2;
    
    // Compute autocorrelation length
    size_t len=vector_autocorr_tau(ac_avg,ftom2);
    cout << "Averaged data, ";
    if (len>0) {
      cout << "autocorrelation length: " << len << endl;
    } else {
      cout << "autocorrelation length determination failed." << endl;
    }

    // Create new columns if necessary
    if (!table_obj.is_column(in[0])) {
      table_obj.new_column(in[0]);
    }
    if (!table_obj.is_column(in[1])) {
      table_obj.new_column(in[1]);
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
      if (i<ftom2.size()) {
	table_obj.set(in[1],i,ftom2[i]);
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

    command_del(type);
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

    command_del(type);
    clear_obj();
    doublev_obj=ac_vec;
    command_add("double[]");
    type="double[]";
    
  } else {
    
    vector<string> in, pr;
    if (sv.size()<2) {
      pr.push_back("Enter multiple vector specification for data");
      int ret=get_input(sv,pr,in,"autocorr",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t i=1;i<sv.size();i++) in.push_back(sv[i]);
    }

    vector<vector<double> > ac_all;
    size_t max_ac_size=0;
    
    for(size_t ix=0;ix<in.size();ix++) {
      
      vector<vector<double> > v;

      int vs_ret=o2scl_hdf::mult_vector_spec(in[ix],v,verbose,false);
      if (vs_ret!=0) {
	cout << "Vector specification failed (returned " << vs_ret << ")."
	     << endl;
	return 1;
      }

      for(size_t j=0;j<v.size();j++) {
	
	vector<double> ac, ftom;
	
	// Compute autocorrelation length and sample size
	vector_autocorr_vector(v[j],ac);
	size_t len=vector_autocorr_tau(ac,ftom);
	if (len>0) {
	  cout << "Autocorrelation length: " << len << " sample size: "
	       << v[j].size()/len << endl;
	} else {
	  cout << "Autocorrelation length determination failed." << endl;
	}
	
	if (j==0 || max_ac_size<ac.size()) {
	  max_ac_size=ac.size();
	}
	
	ac_all.push_back(ac);
      }
    }
    
    if (verbose>0) {
      cout << "Storing autocorrelation coefficient averaged "
	   << "over all data sets in\n  double[] object." << endl;
    }
    
    if (max_ac_size==0) {
      cerr << "Failed to find any data in 'autocorr'." << endl;
      return 2;
    }
    
    // Average over all of the autocorrelation coefficients from the
    // specified data sets and store the result in doublev_obj .
    doublev_obj.resize(max_ac_size);
    for(size_t i=0;i<max_ac_size;i++) {
      size_t n=0;
      doublev_obj[i]=0.0;
      for(size_t j=0;j<ac_all.size();j++) {
	if (i<ac_all[j].size()) {
	  n++;
	  doublev_obj[i]+=ac_all[j][i];
	}
      }
      doublev_obj[i]/=((double)n);
    }

    // Now report autocorrelation length and sample size from
    // averaged result
    vector<double> ftom2;
    
    // Compute autocorrelation length
    size_t len=vector_autocorr_tau(doublev_obj,ftom2);
    cout << "Averaged data, ";
    if (len>0) {
      cout << "autocorrelation length: " << len << endl;
    } else {
      cout << "autocorrelation length determination failed." << endl;
    }

    // If no current object, then set current object as double[]
    if (type=="") {
      command_del(type);
      clear_obj();
      command_add("double[]");
      type="double[]";
    }
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
  double d;
  int retx=o2scl::function_to_double_nothrow(i1,d);
  if (retx!=0) {
    cerr << "Converting " << i1 << " to value failed." << endl;
    return 1;
  }
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(prec);
  if (verbose>0) cout << "Result: ";
  cout << d << endl;
  return 0;
}

int acol_manager::comm_clear(std::vector<std::string> &sv, bool itive_com) {

  command_del(type);

  // The clear_obj() function sets type to an empty string.
  clear_obj();
  
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

    command_del(cur_type);
    command_add(temp_type);
    
    std::vector<std::string>::iterator it=sv.begin();
    it++;
    sv.erase(it);
    int ret=cl->comm_option_commands(sv,itive_com);

    command_del(temp_type);
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

  /*
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
  */
    
  std::string svalue, file, name="contours";

  if (type=="table3d") {

    std::string slice;

    if (sv.size()<3) {
      // If not enough arguments were given, then prompt for them
      vector<string> pr, in;
      pr.push_back("Contour value or \"frac\" and contour value");
      pr.push_back("Slice name");
      pr.push_back("Filename (or \"none\")");
      int ret=get_input(sv,pr,in,"contours",itive_com);
      if (ret!=0) return ret;
      
      if (in[0].find("frac ")==0) {
	in[0]=in[0].substr(5,in[0].length()-5);
	frac_mode=true;
      }

      if (in[2]!="none") {
	file=in[2];
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==3) {
      if (sv[1].find("frac ")==0) {
	sv[1]=sv[1].substr(5,sv[1].length()-5);
	frac_mode=true;
      }
      svalue=sv[1];
      slice=sv[2];
    } else {
      if (sv[1]=="frac") {
	svalue=sv[2];
	slice=sv[3];
	if (sv.size()>4) file=sv[4];
	if (sv.size()>5) name=sv[5];
      } else {
	if (sv[1].find("frac ")==0) {
	  sv[1]=sv[1].substr(5,sv[1].length()-5);
	  frac_mode=true;
	}
	svalue=sv[1];
	slice=sv[2];
	file=sv[3];
	if (sv.size()>4) name=sv[4];
      }
    }
    
    ubvector levs(1);
    int retx=o2scl::function_to_double_nothrow(svalue,levs[0]);
    if (retx!=0) {
      cerr << "Failed to convert " << svalue << " to value." << endl;
      return 1;
    }
    size_t nlev=1;
    
    if (frac_mode) {
      cout << "Fraction mode not implemented with table3d objects." << endl;
      // Get references to the histogram data
      size_t nx=table3d_obj.get_nx();
      size_t ny=table3d_obj.get_ny();
      const ubmatrix &m=table3d_obj.get_slice(slice);
      const ubvector &xd=table3d_obj.get_x_data();
      const ubvector &yd=table3d_obj.get_y_data();

      // Construct bin vectors
      ubvector xbins(nx+1);
      if (xd[1]>xd[0]) {
	xbins[0]=xd[0]-(xd[1]-xd[0])/2.0;
	xbins[nx]=xd[nx-1]+(xd[nx-1]-xd[nx-2])/2.0;
      } else {
	xbins[0]=xd[0]+(xd[0]-xd[1])/2.0;
	xbins[nx]=xd[nx-1]-(xd[nx-2]-xd[nx-1])/2.0;
      }
      for(size_t i=1;i<nx-1;i++) {
	xbins[i]=(xd[i-1]+xd[i])/2.0;
      }
      ubvector ybins(ny+1);
      if (yd[1]>yd[0]) {
	ybins[0]=yd[0]-(yd[1]-yd[0])/2.0;
	ybins[ny]=yd[ny-1]+(yd[ny-1]-yd[ny-2])/2.0;
      } else {
	ybins[0]=yd[0]+(yd[0]-yd[1])/2.0;
	ybins[ny]=yd[ny-1]-(yd[ny-2]-yd[ny-1])/2.0;
      }
      for(size_t i=1;i<ny-1;i++) {
	ybins[i]=(yd[i-1]+yd[i])/2.0;
      }

      // Compute the total integral and the target fraction
      double min, max;
      o2scl::matrix_minmax(m,min,max);
      double sum=matrix_sum<ubmatrix,double>(nx,ny,m);
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
    
    if (file.length()>0) {
      std::vector<contour_line> clines;
      table3d_obj.slice_contours(slice,1,levs,clines);
      if (clines.size()>0) {
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,clines,name);
	hf.close();
      }
    } else {
      table3d_obj.slice_contours(slice,1,levs,cont_obj);
      if (cont_obj.size()>0) {
	command_del(type);
	clear_obj();
	command_add("vector<contour_line>");
	type="vector<contour_line>";
      } else {
	cout << "No contours found. Leaving table3d object unmodified."
	     << endl;
	return 1;
      }
    }
    
  } else if (type=="hist_2d") {

    if (sv.size()<2) {
      // If not enough arguments were given, then prompt for them
      vector<string> pr, in;
      pr.push_back("Contour value or \"frac\" and contour value");
      pr.push_back("Filename (or \"none\")");
      int ret=get_input(sv,pr,in,"contours",itive_com);
      if (ret!=0) return ret;
      
      if (in[0].find("frac ")==0) {
	in[0]=in[0].substr(5,in[0].length()-5);
	frac_mode=true;
      }
      
      if (in[1]!="none") {
	file=in[1];
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==2) {
      if (sv[1].find("frac ")==0) {
	sv[1]=sv[1].substr(5,sv[1].length()-5);
	frac_mode=true;
      }
      svalue=sv[1];
    } else {
      if (sv[1]=="frac") {
	svalue=sv[2];
	if (sv.size()>3) file=sv[3];
	if (sv.size()>4) name=sv[4];
      } else {
	if (sv[1].find("frac ")==0) {
	  sv[1]=sv[1].substr(5,sv[1].length()-5);
	  frac_mode=true;
	}
	svalue=sv[1];
	file=sv[2];
	if (sv.size()>3) name=sv[3];
      }
    }

    ubvector levs(1);
    int retx=o2scl::function_to_double_nothrow(svalue,levs[0]);
    if (retx!=0) {
      cerr << "Failed to convert " << svalue << " to value." << endl;
      return 1;
    }
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
      if (clines.size()>0) {
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,clines,name);
	hf.close();
      }
    } else {
      co.calc_contours(cont_obj);
      if (cont_obj.size()>0) {
	command_del(type);
	clear_obj();
	command_add("vector<contour_line>");
	type="vector<contour_line>";
      } else {
	cout << "No contours found. Leaving hist_2d object unmodified."
	     << endl;
	return 1;
      }
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
  command_del(type);
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
    int vsret=value_spec(tval,double_obj,verbose,false);
    if (vsret!=0) {
      cerr << "Function value_spec() failed." << endl;
      return 1;
    }
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

      int ret2=vector_spec(in1,doublev_obj,verbose,false);
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
    pr.push_back("Vector specification for column");
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    std::vector<double> d;
    int vs_ret=vector_spec(in[1],d,verbose,false);
    if (vs_ret!=0) {
      cerr << "Vector specification " << in[1] << " failed." << endl;
    }
    
    table_obj.clear();
    table_obj.line_of_names(in[0]);
    table_obj.set_nlines(d.size());
    
    for(size_t li=0;li<d.size();li++) {
      table_obj.set(in[0],li,d[li]);
    }
    command_add("table");
    type="table";
    
  } else if (ctype=="tensor") {

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

    tensor_obj.resize(rank,sarr);
    tensor_obj.set_all(0.0);
    command_add("tensor");
    type="tensor";

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
    vector<double> grid;
    for(size_t k=0;k<rank;k++) {
      for(size_t j=0;j<sarr[k];j++) {
	grid.push_back((double)j);
      }
    }
    tensor_grid_obj.set_grid_packed(grid);
    
    tensor_grid_obj.set_all(0.0);
    command_add("tensor_grid");
    type="tensor_grid";

  } else if (ctype=="table3d") {

    vector<string> in;
    vector<string> pr=
      {"Enter x-axis name",
       "Enter x-axis vector specification",
       "Enter y-axis name",
       "Enter y-axis vector specification",
       "Enter slice name",
       "Enter slice function"};
    int ret=get_input(sv2,pr,in,"create",itive_com);
    if (ret!=0) return ret;

    std::vector<double> dx, dy;
    
    std::string xname=in[0];
    int vs_ret=vector_spec(in[1],dx,0,false);
    if (vs_ret!=0) {
      cerr << "Vector specification " << in[1] << " failed." << endl;
    }

    std::string yname=in[2];
    vs_ret=vector_spec(in[3],dy,0,false);
    if (vs_ret!=0) {
      cerr << "Vector specification " << in[3] << " failed." << endl;
    }
    
    std::string zname=in[4];
    std::string zfunc=in[5];
    
    table3d_obj.set_xy(xname,dx.size(),dx,yname,dy.size(),dy);
    
    table3d_obj.function_slice(zfunc,zname);
    
    command_add("table3d");
    type="table3d";

  } else {

    cerr << "Cannot create object of type " << ctype << endl;
    return 1;
      
  }

  return 0;
}

