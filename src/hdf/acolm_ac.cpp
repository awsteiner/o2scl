/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/acolm.h>
#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>
#include <o2scl/set_fftw.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/*
  int acol_manager::comm_ac_len()
  int acol_manager::comm_add_vec()
  int acol_manager::comm_assign()
  int acol_manager::comm_autocorr()
  int acol_manager::comm_average_rows()
  int acol_manager::comm_binary()
  int acol_manager::comm_calc()
  int acol_manager::comm_cat()
  int acol_manager::comm_clear()
  int acol_manager::comm_commands()
  int acol_manager::comm_constant()
  int acol_manager::comm_contours()
  int acol_manager::comm_convert()
  int acol_manager::comm_correl()
  int acol_manager::comm_convert_unit()
  int acol_manager::comm_create()
*/

int acol_manager::comm_add_vec(std::vector<std::string> &sv, bool itive_com) {

  if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "No table to add a column to." << endl;
      return exc_efailed;
    }
    
    vector<string> pr, in;
    pr.push_back("Enter vector specification for new column");
    pr.push_back("Enter name for new column");
    int ret=get_input(sv,pr,in,"function",itive_com);
    if (ret!=0) return ret;
    
    // Remove single or double quotes just in case
    if (in[0].size()>=3 && ((in[0][0]=='\'' && in[0][in[0].size()-1]=='\'') ||
			    (in[0][0]=='\"' && in[0][in[0].size()-1]=='\"'))) {
      in[0]=in[0].substr(1,in[0].size()-2);
    }

    std::vector<double> vec;
    int vret=vector_spec(in[0],vec,false,verbose,false);
    if (vret!=0) {
      cerr << "Vector spec failed in add-vec." << endl;
      return exc_efailed;
    }
    if (vec.size()<table_obj.get_nlines()) {
      cerr << "Not enough entries in vector." << endl;
      cout << vec.size() << " " << table_obj.get_nlines() << endl;
      return exc_efailed;
    }

    table_obj.new_column(in[1]);
    for(size_t i=0;i<table_obj.get_nlines();i++) {
      table_obj.set(in[1],i,vec[i]);
    }
    
  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  return 0;
}

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
  int retx=function_to_double_nothrow(sv[2],d,0,&rng);
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

  string options="";
  if (sv.size()>=2) {
    options=sv[1];
  }
  vector<string> voptions;
  split_string_delim(options,voptions,',');
  
  size_t ix;

  bool store=false;
  if (vector_search(voptions,"store",ix)==true) {
    store=true;
    if (verbose>0) {
      cout << "autocorr storing autocorrelation vector(s)." << endl;
    }
  }

  // Either 'def', 'acor', or 'fft'
  string alg;
  if (vector_search(voptions,"acor",ix)==true) {
    alg="acor";
    if (verbose>0) {
      cout << "autocorr using acor algorithm." << endl;
    }
  } else if (vector_search(voptions,"fft",ix)==true) {
    alg="fft";
    if (verbose>0) {
      cout << "autocorr using FFT algorithm." << endl;
    }
  } else {
    alg="def";
    if (verbose>0) {
      cout << "autocorr using default algorithm." << endl;
    }
  }
  
  // Either 'max' or 'avg'
  string combine;
  if (vector_search(voptions,"avg",ix)==true) {
    combine="avg";
    if (verbose>0) {
      cout << "autocorr will average autocorrelation data." << endl;
    }
  } else {
    combine="max";
    if (verbose>0) {
      cout << "autocorr will report maximum autocorrelation length." << endl;
    }
  }

  // First, collect all the data into the object name vvd
  vector<vector<double>> vvd;

  if (type=="int[]") {

    vvd.resize(1);
    vector_copy(intv_obj,vvd[0]);
    
  } else if (type=="size_t[]") {
    
    vvd.resize(1);
    vector_copy(size_tv_obj,vvd[0]);
    
  } else if (type=="double[]") {
    
    vvd.resize(1);
    vector_copy(doublev_obj,vvd[0]);
    
  } else if (type=="table") {
    
    if (table_obj.get_nlines()==0) {
      cerr << "Table has no lines of data to compute "
	   << "autocorrelations with." << endl;
      return exc_efailed;
    }

    vector<string> in;
    for(size_t i=2;i<sv.size();i++) in.push_back(sv[i]);
    cout << "in: ";
    vector_out(cout,in,true);
    
    if (in.size()==0) {
      cerr << "No columns or vectors specified." << endl;
      return 2;
    }

    for(size_t jx=0;jx<in.size();jx++) {

      // Determine vector from table column (requires copy) or
      // multiple vector specification
      
      if (in[jx].find(':')==std::string::npos) {
        cout << "column: " << in[jx] << endl;

        if (table_obj.is_column(in[jx])==false) {
          cerr << "Could not find column named '" << in[jx] << "'." << endl;
          return exc_efailed;
        }

        vector<double> v;
	v.resize(table_obj.get_nlines());
	for(size_t i=0;i<table_obj.get_nlines();i++) {
	  v[i]=table_obj.get(in[jx],i);
	}
        vvd.push_back(v);
        
      } else {

        vector<vector<double> > vvd2;
        
        int vs_ret=o2scl_hdf::mult_vector_spec(in[jx],vvd2,verbose,false);
	if (vs_ret!=0) {
	  cerr << "Multiple vector specification failed." << endl;
	  return 1;
	}

        for(size_t kk=0;kk<vvd2.size();kk++) {
          vvd.push_back(vvd[kk]);
        }
      }
      
    }
    
  } else {
    
    vector<string> in, pr;
    if (sv.size()<3) {
      pr.push_back("Enter options");
      pr.push_back("Enter multiple vector specification for data");
      int ret=get_input(sv,pr,in,"autocorr",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t i=1;i<sv.size();i++) in.push_back(sv[i]);
    }

    vector<vector<double> > vvd2;
    
    int vs_ret=o2scl_hdf::mult_vector_spec(in[1],vvd2,false,verbose,false);
    if (vs_ret!=0) {
      cerr << "Multiple vector specification failed." << endl;
      return 1;
    }
    
    for(size_t kk=0;kk<vvd2.size();kk++) {
      vvd.push_back(vvd2[kk]);
    }
    
  }

  // In this section, compute the autocorrelation vector for
  // each vector inside "vvd" and store it in "ac_vec".
  // Also, compute "max_ac_size" which is the maximum size of
  // an object in "ac_vec".
  
  vector<vector<double>> ac_vec;
  
  size_t max_ac_size=0;
  ac_vec.resize(vvd.size());
  for(size_t jj=0;jj<vvd.size();jj++) {
    if (alg=="acor") {
      double mean, sigma, tau;
      vector_acor<vector<double>>(vvd[jj].size(),vvd[jj],
                                  mean,sigma,tau);
    } else if (alg=="fft") {
#ifdef O2SCL_SET_FFTW
      double mean=vector_mean(vvd[jj]);
      double stddev=vector_stddev(vvd[jj]);
      vector_autocorr_vector_fftw(vvd[jj],ac_vec[jj],mean,stddev);
#else
      cerr << "FFTW support not included." << endl;
      return 1;
#endif
    } else {
      vector_autocorr_vector(vvd[jj].size(),vvd[jj],ac_vec[jj]);
    }
    if (ac_vec[jj].size()>max_ac_size) {
      max_ac_size=ac_vec[jj].size();
    }
  }

  // Temporary vectors for autocorrelation length computation
  vector<vector<double>> ftom;
  
  if (combine=="max") {

    // In "max" mode, compute all the autocorrelation coefficients and
    // just report the largest autocorrelation length
    size_t len_max=0;
  
    ftom.resize(vvd.size());
    for(size_t jj=0;jj<vvd.size();jj++) {
      size_t len=vector_autocorr_tau(ac_vec[jj],ftom[jj]);
      if (len>0) {
        cout << "Autocorrelation length: " << len << " sample size: "
             << vvd[jj].size()/len << endl;
      } else {
        cout << "Autocorrelation length determination failed." << endl;
      }
      if (len>len_max) len=len_max;
    }
    
  } else {

    // In "avg" mode, compute all of the autocorrelation coefficients
    // and then compute the average.
    
    vector<double> ac_avg(max_ac_size);
    for(size_t i=0;i<max_ac_size;i++) {
      size_t n=0;
      ac_avg[i]=0.0;
      for(size_t j=0;j<ac_vec.size();j++) {
	if (i<ac_vec[j].size()) {
	  n++;
	  ac_avg[i]+=ac_vec[j][i];
          if (verbose>2 && i<10) {
            cout << "i,j,ac: " << i << " " << j << " " << ac_vec[j][i]
                 << endl;
          }
	}
      }
      ac_avg[i]/=((double)n);
      if (verbose>2 && i<10) {
        cout << "avg: " << ac_avg[i] << endl;
      }
    }

    // Use the average to compute the autocorrelation length
    
    ftom.resize(1);
    size_t len=vector_autocorr_tau(ac_avg,ftom[0]);
    if (len>0) {
      cout << "Autocorrelation length: " << len << "." << endl;
    } else {
      cout << "Autocorrelation length determination failed." << endl;
    }

    // Add the average to the list of autocorrelation coefficient
    // vectors so it can be accessed by the user
    ac_vec.push_back(ac_avg);
  }

  if (store) {
    std::swap(ac_vec,vvdouble_obj);
    command_del(type);
    clear_obj();
    command_add("vec_vec_double");
    type="vec_vec_double";
  }
    
  return 0;
}

int acol_manager::comm_average_rows(std::vector<std::string> &sv,
				    bool itive_com) {

  if (type!="table") {
    cerr << "No table to average_rows." << endl;
    return exc_efailed;
  }

  if (sv.size()==3) {
    
    if (sv[1]=="*" && sv.size()>=4 && o2scl::stob(sv[3])) {
      table_obj.average_rows(o2scl::stoszt(sv[2]));
    } else {
      table_obj.average_col_roll(sv[1],o2scl::stoszt(sv[2]));
    }

  } else {
    
    vector<string> pr, in;
    pr.push_back("Column name, or '*' for all");
    pr.push_back("Window size (greater than 1)");
    pr.push_back("True for block average, false for rolling average");
    int ret=get_input(sv,pr,in,"average-rows",itive_com);
    if (ret!=0) return ret;
    
    if (o2scl::stob(in[2])==true) {
      if (in[0]=="*") {
	table_obj.average_rows(o2scl::stoszt(in[1]),true);
      } else {
	cerr << "Cannot specify column name with block averages." << endl;
	return 2;
      }
    } else {
      table_obj.average_col_roll(in[0],o2scl::stoszt(in[1]));
    }

  }
  
  return 0;
}

int acol_manager::comm_binary(std::vector<std::string> &sv, bool itive_com) {

  if (type=="tensor_grid") {

    std::string function, fname, oname;

    vector<string> pr, in;
    pr.push_back("Filename of second object");
    pr.push_back("Name of second object.");
    pr.push_back("Enter function of i0,i1,... and x0,x1,...");
    int ret=get_input(sv,pr,in,"binary",itive_com);

    fname=in[0];
    oname=in[1];
    function=in[2];

    hdf_file hf;
    wordexp_single_file(fname);
    hf.open(fname);
    tensor_grid<> tg;
    hdf_input(hf,tg,oname);
    hf.close();

    if (tg.get_rank()!=tensor_grid_obj.get_rank()) {
      cerr << "Ranks do not match." << endl;
      return 2;
    }

    if (tg.total_size()!=tensor_grid_obj.total_size()) {
      cerr << "Sizes do not match." << endl;
      return 3;
    }

#ifdef O2SCL_SET_OPENMP
#pragma omp parallel
    {
#endif
      
      // Parse function(s)
      calc_utf8<> calc;
      calc.set_rng(rng);
      std::map<std::string,double> vars;
      calc.compile(function.c_str(),&vars);
      
      // Set
      size_t rk=tensor_grid_obj.get_rank();
      vector<size_t> ix(rk);
      
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
      for(size_t i=0;i<tensor_grid_obj.total_size();i++) {
        tensor_grid_obj.unpack_index(i,ix);
        vector<double> xa;
        for(size_t j=0;j<rk;j++) {
          vars[((string)"i")+szttos(j)]=ix[j];
          vars[((string)"x")+szttos(j)]=tensor_grid_obj.get_grid(j,ix[j]);
          xa.push_back(tensor_grid_obj.get_grid(j,ix[j]));
        }
        vars["v"]=tensor_grid_obj.get(ix);
        vars["w"]=tg.interp_linear(xa);
        tensor_grid_obj.set(ix,calc.eval(&vars));
      }
      
#ifdef O2SCL_SET_OPENMP
      // End of parallel region
    }
#endif
    
  } else {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }

  return 0;
}

int acol_manager::comm_calc(std::vector<std::string> &sv, bool itive_com) {

  // The format_float class can't handle higher precisions
  if (precision<=14) {
    ff.set_sig_figs(precision+1);
  }

#ifdef O2SCL_SET_MULTIP
  
  if (sv.size()>2 && o2scl::stob(sv[2])==true) {
    
    std::string i1=sv[1];
    
    funct_multip_string fms;
    fms.verbose=verbose;
    fms.set_function(i1,"x");
    funct_multip_string *fmsp=&fms;
    
    funct_multip fm2;
    fm2.verbose=verbose;
    
    // Note the funct_multip_string object uses a tolerance of
    // pow(10.0,-std::numeric_limits<fp_t>::digits10+1), a factor of
    // 10 different from digits10, and then when cout.precision is 6,
    // actually 7 significant figures are output, so we need a two
    // digit buffer thus, e.g., anything over 33 digit precision
    // requires 35-digit floats.
    
    if (precision>48) {
      
      cerr << "Requested precision too large for the calc "
           << "command (maximum is 48)." << endl;
      return 2;

    } else if (precision>33) {
      
      cpp_dec_float_50 d=0, err;
      int retx=fm2.eval_tol_err([fmsp](auto &&t) mutable
      { return (*fmsp)(t); },d,d,err);
        
      if (retx!=0) {
        cerr << "Converting " << i1 << " to value failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_50): ";
      cout << dtos(d,precision) << endl;
      return 0;
      
    } else if (precision>23) {
      
      cpp_dec_float_35 d=0, err;
      int retx=fm2.eval_tol_err([fmsp](auto &&t) mutable
      { return (*fmsp)(t); },d,d,err);
        
      if (retx!=0) {
        cerr << "Converting " << i1 << " to value failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_35): ";
      cout << dtos(d,precision) << endl;
      return 0;
      
    } else if (precision>16) {
      
      cpp_dec_float_25 d=0, err;
      int retx=fm2.eval_tol_err([fmsp](auto &&t) mutable
      { return (*fmsp)(t); },d,d,err);
        
      if (retx!=0) {
        cerr << "Converting " << i1 << " to value failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (cpp_dec_float_25): ";
      cout << dtos(d,precision) << endl;
      
      return 0;
      
    } else if (precision>13) {
      
      long double d=0, err;
      int retx=fm2.eval_tol_err([fmsp](auto &&t) mutable
      { return (*fmsp)(t); },d,d,err);
        
      if (retx!=0) {
        cerr << "Converting " << i1 << " to value failed." << endl;
        return 1;
      }
      if (verbose>0) cout << "Result (long double): ";
      cout << dtos(d,precision) << " " << endl;
      
      return 0;
    }
    
    double d=0, err;
    int retx=fm2.eval_tol_err([fmsp](auto &&t) mutable
    { return (*fmsp)(t); },d,d,err);
      
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (scientific) cout.setf(ios::scientific);
    else cout.unsetf(ios::scientific);
    cout.precision(precision);
    if (verbose>0) cout << "Result: ";
    cout << d << " (" << ff.convert(d) << ")" << endl;
    return 0;
    
  }

#else

  cout << "No multiprecision support was included."
       << " was defined." << endl;
  
#endif
  
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

  std::cout << "Precision: " << precision << std::endl;
  
  if (precision>100) {
    cerr << "Requested precision too large for the calc "
         << "command." << endl;
    return 2;
    
#ifdef O2SCL_SET_MULTIP
    
  } else if (precision>50) {
    cpp_dec_float_100 d;
    convert_units<cpp_dec_float_100> cu100;
    int retx=o2scl::function_to_fp_nothrow<cpp_dec_float_100>
      (i1,d,cu100,verbose);
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (verbose>0) cout << "Result (cpp_dec_float_100): ";
    cout << dtos(d,precision) << endl;
    return 0;
  } else if (precision>35) {
    cpp_dec_float_50 d;
    convert_units<cpp_dec_float_50> cu50;
    int retx=o2scl::function_to_fp_nothrow<cpp_dec_float_50>
      (i1,d,cu50,verbose);
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (verbose>0) cout << "Result (cpp_dec_float_50): ";
    cout << dtos(d,precision) << endl;
    return 0;
  } else if (precision>25) {
    cpp_dec_float_35 d;
    convert_units<cpp_dec_float_35> cu35;
    int retx=o2scl::function_to_fp_nothrow<cpp_dec_float_35>
      (i1,d,cu35,verbose);
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (verbose>0) cout << "Result (cpp_dec_float_35): ";
    cout << dtos(d,precision) << endl;
    return 0;
  } else if (precision>18) {
    cpp_dec_float_25 d;
    convert_units<cpp_dec_float_25> cu25;
    int retx=o2scl::function_to_fp_nothrow<cpp_dec_float_25>
      (i1,d,cu25,verbose);
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (verbose>0) cout << "Result (cpp_dec_float_25): ";
    cout << dtos(d,precision) << endl;
    return 0;
  } else if (precision>15) {
    long double d;
    convert_units<long double> culd;
    int retx=o2scl::function_to_fp_nothrow<long double>
      (i1,d,culd,verbose,&rng);
    if (retx!=0) {
      cerr << "Converting " << i1 << " to value failed." << endl;
      return 1;
    }
    if (verbose>0) cout << "Result (long double): ";
    cout << dtos(d,precision) << endl;
    return 0;
#endif
  }
  
  double d;
  int retx=o2scl::function_to_double_nothrow(i1,d,verbose,&rng);
  if (retx!=0) {
    cerr << "Converting " << i1 << " to value failed." << endl;
    return 1;
  }
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(precision);
  if (verbose>0) cout << "Result: ";
  cout << d;
  if (precision<=14) {
    cout << " (" << ff.convert(d) << ")";
  }
  cout << endl;
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

    const ubvector &xg1=table3d_obj.get_x_data();
    const ubvector &yg1=table3d_obj.get_y_data();
    const ubvector &xg2=tab2.get_x_data();
    const ubvector &yg2=tab2.get_y_data();
    bool grids_match=(vectors_equal(xg1,xg2) &&
                      vectors_equal(yg1,yg2));
    
    for(size_t k=0;k<tab2.get_nslices();k++) {
      std::string sl_name=tab2.get_slice_name(k);
      size_t slix;
      if (!table3d_obj.is_slice(sl_name,slix)) {
	table3d_obj.new_slice(sl_name);
	for(size_t i=0;i<table3d_obj.get_nx();i++) {
	  for(size_t j=0;j<table3d_obj.get_ny();j++) {
            if (grids_match) {
              table3d_obj.set_val(i,j,sl_name,tab2.get(i,j,sl_name));
            } else {
              double x=xg1[i];
              double y=yg1[j];
              table3d_obj.set_val(i,j,sl_name,tab2.interp(x,y,sl_name));
            }
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

    // Make space in the new table, copying to a new memory segment if
    // necessary.
    table_obj.set_nlines(n1+n2);

    // Loop through all the column names in the second table
    for(size_t j=0;j<tab2.get_ncolumns();j++) {
      
      std::string col_name=tab2.get_column_name(j);
      
      // If the current table does not have a column with that name
      // then create the new column and set the first n1 entries to
      // zero. 
      if (!table_obj.is_column(col_name)) {
        if (verbose>0) {
          cout << "Original table does not contain column " << col_name
               << ", so this new column will be created and initialized to "
               << "zero before the concatenation is performed." << endl;
        }
	table_obj.new_column(col_name);
	for(size_t i=0;i<n1;i++) table_obj.set(col_name,i,0.0);
      }

      // Copy the information from the second table to the current table
      for(size_t i=0;i<n2;i++) {
	table_obj.set(col_name,i+n1,tab2.get(col_name,i));
      }
    }
    
    // Look for columns in the new table which did not have data in
    // the second table and set their new entries to zero.
    for(size_t j=0;j<table_obj.get_ncolumns();j++) {
      std::string col_name=table_obj.get_column_name(j);
      if (!tab2.is_column(col_name)) {
        if (verbose>0) {
          cout << "The table to be added does not contain column "
               << col_name << ", so the additional rows in this column "
               << "will be set to zero." << endl;
        }
	for(size_t i=n1;i<n1+n2;i++) table_obj.set(col_name,i,0.0);
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

int acol_manager::comm_clear(std::vector<std::string> &sv, bool itive_com) {

  command_del(type);

  // The clear_obj() function sets type to an empty string.
  clear_obj();

  return 0;
}

int acol_manager::comm_commands(std::vector<std::string> &sv,
				bool itive_com) {

  // FIXME: there is some code duplication in this function
  // which could be removed.
  
  terminal ter;

  if (type.length()>0) {
    cout << "Current object named \"" << obj_name << "\" has type "
         << type_color << type << default_color
         << ".\n" << endl;
  } else {
    cout << "No current object.\n" << endl;
  }
    
  if (sv.size()==2) {

    if (sv[1]=="all") {
      
      cout << "Commands which do not require a current object:\n" << endl;

      // Delete the current type from the cli object temporarily
      // so we can get a list of commands which don't require an
      // object.
      
      std::string curr_type=type;
      command_del(curr_type);
      std::vector<std::string> comm_list=cl->get_option_list();
      command_add(curr_type);

      // Sort, screenify, and output
      
      std::vector<std::string> comm_out;
      vector_sort<std::vector<std::string>,std::string>
        (comm_list.size(),comm_list);
      for(size_t j=0;j<comm_list.size();j++) {
	comm_list[j]=command_color+comm_list[j]+default_color;
      }
      screenify(comm_list.size(),comm_list,comm_out);
      for(size_t j=0;j<comm_out.size();j++) {
	cout << comm_out[j] << endl;
      }
      cout << endl;

      // Proceed through all of the types, go through type_comm_list
      // to output all of the type-specific commands (this code
      // presumes that object is already sorted).
      
      std::map<std::string,std::vector<std::string> >::iterator it;
      for(it=type_comm_list.begin();it!=type_comm_list.end();it++) {
	cout << "Commands for an object of type " << type_color
	     << it->first << default_color
	     << ":\n" << endl;
	std::vector<std::string> clist=it->second;
	for(size_t j=0;j<clist.size();j++) {
	  clist[j]=command_color+clist[j]+default_color;
	}
	comm_out.clear();
	screenify(clist.size(),clist,comm_out);
	for(size_t j=0;j<comm_out.size();j++) {
	  cout << comm_out[j] << endl;
	}
	cout << endl;
      }
      
      return 0;
      
    } else {
      
      string temp_type=sv[1];
      string curr_type=type;

      if (std::find(type_list.begin(),type_list.end(),
                    temp_type)==type_list.end()) {
        cerr << "Type " << sv[1] << " is not a valid "
             << cl->cmd_name << " type." << endl;
        return 1;
      }
      
      cout << "Commands which do not require a current object:\n" << endl;

      // Delete the current type from the cli object temporarily
      // so we can get a list of commands which don't require an
      // object
      
      command_del(curr_type);
      std::vector<std::string> comm_list=cl->get_option_list();
      command_add(curr_type);
      
      // Sort, decorate, screenify, and output

      std::vector<std::string> comm_out;
      vector_sort<std::vector<std::string>,std::string>
        (comm_list.size(),comm_list);
      std::vector<std::string> comm_list_decor(comm_list.size());
      for(size_t j=0;j<comm_list.size();j++) {
	comm_list_decor[j]=command_color+
          comm_list[j]+default_color;
      }
      screenify(comm_list_decor.size(),comm_list_decor,comm_out);
      for(size_t j=0;j<comm_out.size();j++) {
	cout << comm_out[j] << endl;
      }
      cout << endl;
      
      cout << "Commands for an object of type " << type_color
           << temp_type << default_color
           << ":\n" << endl;
      
      command_del(curr_type);
      command_add(temp_type);

      // The list of commands for an object of type 'temp_type'
      std::vector<std::string> comm_list_type=cl->get_option_list();

      // Subtract out the generic list of commands
      std::vector<std::string> comm_list_type_only;
      for(size_t i=0;i<comm_list_type.size();i++) {
        if (std::find(comm_list.begin(),comm_list.end(),
                      comm_list_type[i])==comm_list.end()) {
          comm_list_type_only.push_back(comm_list_type[i]);
        }
      }
      vector_sort<vector<string>,string>(comm_list_type_only.size(),
                                         comm_list_type_only);

      // Decorate, screenify, and output
      comm_out.clear();
      for(size_t j=0;j<comm_list_type_only.size();j++) {
	comm_list_type_only[j]=command_color+
          comm_list_type_only[j]+default_color;
      }
      screenify(comm_list_type_only.size(),comm_list_type_only,comm_out);
      for(size_t j=0;j<comm_out.size();j++) {
	cout << comm_out[j] << endl;
      }
      cout << endl;

      // Restore the commands for the original type
      command_del(temp_type);
      command_add(curr_type);
      
      return 0;
    }
  }

  // This section applies if there are no arguments to the 'commands'
  // option.
  
  if (type!="") {
    cout << "Commands which do not require a current object or "
	 << "which apply to\n  objects of type " << type << ".\n" << endl;
  } else {
    cout << "Commands which do not require a current object:\n" << endl;
  }
  int ret=cl->comm_option_commands(sv,itive_com);
  cout << "Use '-commands all' for a list of all commands "
       << "for the various types." << endl;
  
  return ret;
}

int acol_manager::comm_constant(std::vector<std::string> &sv,
				bool itive_com) {
  
  convert_units<double> &cu=o2scl_settings.get_convert_units();

  if (sv.size()<2) {
    vector<string> pr, in;
    pr.push_back("Name or search pattern");
    pr.push_back("Unit (or 'none' for any)");
    int ret=get_input(sv,pr,in,"constant",itive_com);
    if (ret!=0) return ret;
    if (verbose>=1) {
      cout << "constant: Looking up constant " << in[0]
           << " with unit " << in[1] << endl;
    }
    cu.verbose=verbose;
    if (precision>50) {
      cerr << "Requested precision too large for the constant "
           << "command (the maximum is 50)." << endl;
#ifdef O2SCL_SET_MULTIP
    } else if (precision>35) {
      convert_units<cpp_dec_float_50> cu50;
      cu50.find_print(in[0],in[1],precision,false);
    } else if (precision>25) {
      convert_units<cpp_dec_float_35> cu35;
      cu35.find_print(in[0],in[1],precision,false);
    } else if (precision>18) {
      convert_units<cpp_dec_float_25> cu25;
      cu25.find_print(in[0],in[1],precision,false);
    } else if (precision>15) {
      convert_units<long double> culd;
      culd.find_print(in[0],in[1],precision,false);
#endif
    } else {
      cu.find_print(in[0],in[1],precision,false);
    }
    return 0;
  }
  
  if (sv.size()>=3 && sv[2]!="none" && sv[2]!="None" && sv[1]!="add" &&
      sv[1]!="del") {
    // If a unit was specified, then search for the constant with
    // the specified unit
    if (verbose>=1) {
      cout << "constant: Looking up constant " << sv[1]
           << " with unit " << sv[2] << endl;
    }
    cu.verbose=verbose;
    cu.find_print(sv[1],sv[2],precision,false);
  } else if (sv[1]=="list") {
    cout.precision(precision);
    cu.fc.output_list(cout);
  } else if (sv[1]=="list-full") {
    cout.precision(precision);
    cu.fc.output_list_full(cout);
  } else if (sv[1]=="add") {
    if (sv.size()<4) {
      cerr << "Argument 'add' given to command 'constant' implies add "
           << "a constant but not enough arguments were given." << endl;
      return 1;
    }
    find_constants<>::const_entry f;
    f.names.push_back(sv[2]);
    f.val=o2scl::stod(sv[3]);
    if (sv.size()>=5) {
      f.unit=sv[4];
    }
    if (sv.size()>=6) {
      f.unit_flag=stoszt(sv[5]);
    } else {
      f.unit_flag=0;
    }
    if (sv.size()>=7) {
      f.source=sv[6];
    }
    if (sv.size()>=8) {
      f.m=o2scl::stoi(sv[7]);
    }
    if (sv.size()>=9) {
      f.k=o2scl::stoi(sv[8]);
    }
    if (sv.size()>=10) {
      f.s=o2scl::stoi(sv[9]);
    }
    if (sv.size()>=11) {
      f.K=o2scl::stoi(sv[10]);
    }
    if (sv.size()>=12) {
      f.A=o2scl::stoi(sv[11]);
    }
    if (sv.size()>=13) {
      f.mol=o2scl::stoi(sv[12]);
    }
    if (sv.size()>=14) {
      f.cd=o2scl::stoi(sv[13]);
    }
    if (sv.size()>=15) {
      for(size_t j=14;j<sv.size();j++) {
        f.names.push_back(sv[j]);
      }
    }
    if (verbose>=1) {
      cout << "constant: Adding constant named " << sv[2]
           << " with value " << sv[3] << endl;
    }
    cu.fc.add_constant(f,verbose);
  } else if (sv[1]=="del") {
    if (sv.size()==2) {
      cerr << "Argument 'del' given to command 'constant' implies delete "
           << "a constant but no name was given." << endl;
      return 1;
    }
    if (verbose>=1) {
      cout << "constant: Deleting constant named " << sv[2] << endl;
    }
    cu.fc.del_constant(sv[2],verbose);
  } else {
    if (verbose>=1) {
      cout << "constant: Printing constant named " << sv[1]
           << " (unit unspecified)" << endl;
    }
    cu.verbose=verbose;
    cu.find_print(sv[1],"",precision,false);
  }

  return 0;
}

int acol_manager::comm_contours(std::vector<std::string> &sv, bool itive_com) {

  // The implementation for the prob_dens_mdim_gaussian type is
  // unfinished.
  
  if (type!="table3d" && type!="hist_2d" &&
      type!="prob_dens_mdim_gaussian") {
    cerr << "Not implemented for type " << type << " ." << endl;
    return exc_efailed;
  }
  
  // Note: Input parsing is handled separately for each type below
  
  std::string svalue, file, name="contours";
  
  if (type=="table3d") {

    /*
      The calculation of fractional integrals is more difficult for a
      table3d object than for a hist_2d object because one can imagine
      integrals of complicated regions inside of contour lines. In the
      histogram case, this is simpler because we can just count bins.
      In order to avoid these complications, fractional integrals are
      now only allowed for hist_2d types. The user may convert a
      table3d slice to a hist_2d object and then compute the contours
      from a fractional integral if they need to do so.
    */
    
    std::string slice;

    if (sv.size()<3) {
      // If not enough arguments were given, then prompt for them
      vector<string> pr, in;
      pr.push_back("Contour value");
      pr.push_back("Slice name");
      pr.push_back("Filename (or \"none\")");
      int ret=get_input(sv,pr,in,"contours",itive_com);
      if (ret!=0) return ret;
      
      if (in[2]!="none") {
	file=in[2];
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==3) {
      svalue=sv[1];
      slice=sv[2];
    } else {
      svalue=sv[1];
      slice=sv[2];
      file=sv[3];
      if (sv.size()>4) name=sv[4];
    }

    // Convert the specified contour level string into a double
    ubvector levs(1);
    int retx=o2scl::function_to_double_nothrow(svalue,levs[0],0,&rng);
    if (retx!=0) {
      cerr << "Failed to convert " << svalue << " to value." << endl;
      return 1;
    }
    size_t nlev=1;

    // Compute the contours
    table3d_obj.slice_contours(slice,1,levs,cont_obj);

    if (file.length()>0) {
      // Write to a file
      if (cont_obj.size()>0) {
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,cont_obj,name);
	hf.close();
      } else {
        cout << "File specified, but no contours found, so no file was "
             << "written." << endl;
        return 0;
      }
    } else {
      // Store as a new object
      if (cont_obj.size()>0) {
	command_del(type);
	clear_obj();
	command_add("vector<contour_line>");
	type="vector<contour_line>";
      } else {
	cout << "No contours found. Leaving table3d object unmodified."
	     << endl;
	return 0;
      }
    }
    
  } else if (type=="hist_2d") {
    
    /*
      Note that when "frac" or "frac2" is present, the contour level
      is not computed with full double precision accuracy, The
      integrals for one hundred levels between the min and max are
      computed by directly summing the histogram weights. This data is
      then linearly interpolated to obtain the appropriate absolute
      contour level. However, the contour calculation is also an
      approximate method based on linear interpolation so no more
      accuracy is likely warranted here.
    */
    int frac_mode=0;

    if (sv.size()<2) {
      // If not enough arguments were given, then prompt for them
      vector<string> pr, in;
      pr.push_back("Contour value or \"frac\" and contour value");
      pr.push_back("Filename (or \"none\")");
      int ret=get_input(sv,pr,in,"contours",itive_com);
      if (ret!=0) return ret;
      
      if (in[0].find("frac ")==0) {
	in[0]=in[0].substr(5,in[0].length()-5);
	frac_mode=1;
      } else if (in[0].find("frac2 ")==0) {
	in[0]=in[0].substr(6,in[0].length()-6);
	frac_mode=2;
      } else if (in[0].find("fracx ")==0) {
	in[0]=in[0].substr(6,in[0].length()-6);
	frac_mode=3;
      } else if (in[0].find("fracy ")==0) {
	in[0]=in[0].substr(6,in[0].length()-6);
	frac_mode=4;
      }
      
      if (in[1]!="none") {
	file=in[1];
	name=cl->cli_gets("Object name (or blank for \"contours\"): ");
	if (name.length()==0) name="contours";
      }
    } else if (sv.size()==2) {
      if (sv[1].find("frac ")==0) {
	sv[1]=sv[1].substr(5,sv[1].length()-5);
	frac_mode=1;
      } else if (sv[1].find("frac2 ")==0) {
	sv[1]=sv[1].substr(6,sv[1].length()-6);
	frac_mode=2;
      } else if (sv[1].find("fracx ")==0) {
	sv[1]=sv[1].substr(6,sv[1].length()-6);
	frac_mode=3;
      } else if (sv[1].find("fracy ")==0) {
	sv[1]=sv[1].substr(6,sv[1].length()-6);
	frac_mode=4;
      }
      svalue=sv[1];
    } else {
      if (sv[1]=="frac") {
	svalue=sv[2];
	if (sv.size()>3) file=sv[3];
	if (sv.size()>4) name=sv[4];
	frac_mode=1;
      } else if (sv[1]=="frac2") {
	svalue=sv[2];
	if (sv.size()>3) file=sv[3];
	if (sv.size()>4) name=sv[4];
	frac_mode=2;
      } else if (sv[1]=="fracx") {
	svalue=sv[2];
	if (sv.size()>3) file=sv[3];
	if (sv.size()>4) name=sv[4];
	frac_mode=3;
      } else if (sv[1]=="fracy") {
	svalue=sv[2];
	if (sv.size()>3) file=sv[3];
	if (sv.size()>4) name=sv[4];
	frac_mode=4;
      } else {
	if (sv[1].find("frac ")==0) {
	  sv[1]=sv[1].substr(5,sv[1].length()-5);
	  frac_mode=1;
	} else if (sv[1].find("frac2 ")==0) {
	  sv[1]=sv[1].substr(6,sv[1].length()-6);
	  frac_mode=2;
	} else if (sv[1].find("fracx ")==0) {
	  sv[1]=sv[1].substr(6,sv[1].length()-6);
	  frac_mode=3;
	} else if (sv[1].find("fracy ")==0) {
	  sv[1]=sv[1].substr(6,sv[1].length()-6);
	  frac_mode=4;
	}
	svalue=sv[1];
	file=sv[2];
	if (sv.size()>3) name=sv[3];
      }
    }

    // Convert the specified contour level string into a double
    ubvector levs(1);
    int retx=o2scl::function_to_double_nothrow(svalue,levs[0],0,&rng);
    if (retx!=0) {
      cerr << "Failed to convert " << svalue << " to value." << endl;
      return 1;
    }

    if (frac_mode==3) {
      
      // Get references to the histogram data
      size_t nx=hist_2d_obj.size_x();
      size_t ny=hist_2d_obj.size_y();

      vector<double> clow, chigh, cy;
      for(size_t j=0;j<ny;j++) {
        vector<double> row_x, row_y;
        bool nonzero=false;
        for(size_t i=0;i<nx;i++) {
          row_x.push_back(hist_2d_obj.get_x_rep_i(i));
          row_y.push_back(hist_2d_obj.get_wgt_i(i,j));
          if (hist_2d_obj.get_wgt_i(i,j)>0.0) {
            nonzero=true;
          }
          if (i<10 && verbose>1) {
            cout << i << " " << hist_2d_obj.get_wgt_i(i,j) << endl;
          }
        }
        if (nonzero) {
          double clt, cht;
          int vbf_ret=vector_bound_fracint(nx,row_x,row_y,
                                           levs[0],clt,cht,0,verbose,false);
          if (vbf_ret==0) {
            clow.push_back(clt);
            chigh.push_back(cht);
            cy.push_back(hist_2d_obj.get_y_rep_i(j));
          }
        } else {
          cout << "Skipping row " << j << " because it's zero."
               << endl;
        }
      }

      if (clow.size()==0) {
        cout << "No regions found to compute the contours for." << endl;
        return 2;
      }

      cont_obj.clear();
      
      contour_line cline;
      cline.level=levs[0];
      cline.x=clow;
      cline.y=cy;
      cont_obj.push_back(cline);
      cline.level=levs[0];
      cline.x=chigh;
      cline.y=cy;
      cont_obj.push_back(cline);
      
    } else if (frac_mode==4) {
      
      // Get references to the histogram data
      size_t nx=hist_2d_obj.size_x();
      size_t ny=hist_2d_obj.size_y();

      vector<double> clow, chigh, cx;
      for(size_t i=0;i<nx;i++) {
        bool nonzero=false;
        vector<double> col_x, col_y;
        for(size_t j=0;j<ny;j++) {
          col_x.push_back(hist_2d_obj.get_wgt_i(i,j));
          col_y.push_back(hist_2d_obj.get_y_rep_i(i));
          if (hist_2d_obj.get_wgt_i(i,j)>0.0) {
            nonzero=true;
          }
          if (i<10 && verbose>1) {
            cout << i << " " << hist_2d_obj.get_wgt_i(i,j) << endl;
          }
        }
        if (nonzero) {
          double clt, cht;
          int vbf_ret=vector_bound_fracint(nx,col_x,col_y,
                                           levs[0],clt,cht,0,verbose,false);
          if (vbf_ret==0) {
            clow.push_back(clt);
            chigh.push_back(cht);
            cx.push_back(hist_2d_obj.get_x_rep_i(i));
          }
        } else {
          cout << "Skipping column " << i << " because it's zero."
               << endl;
        }
      }

      if (clow.size()==0) {
        cout << "No regions found to compute the contours for." << endl;
        return 2;
      }

      cont_obj.clear();
      
      contour_line cline;
      cline.level=levs[0];
      cline.x=cx;
      cline.y=clow;
      cont_obj.push_back(cline);
      cline.level=levs[0];
      cline.x=cx;
      cline.y=chigh;
      cont_obj.push_back(cline);
      
    } else {

      if (frac_mode>=1) {
        
        // Compute the contour level from the fractional integral
        
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
        if (frac_mode==1) {
          for(size_t i=0;i<nx;i++) {
            for(size_t j=0;j<ny;j++) {
              sum-=min*(xbins[i+1]-xbins[i])*(ybins[j+1]-ybins[j]);
            }
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
                if (frac_mode==1) {
                  integy[k]+=(m(i,j)-min)*(xbins[i+1]-xbins[i])*
                    (ybins[j+1]-ybins[j]);
                } else {
                  integy[k]+=m(i,j)*(xbins[i+1]-xbins[i])*
                    (ybins[j+1]-ybins[j]);
                }
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
          if (verbose>0) {
            cout << "target: " << target << endl;
            cout << "k x y: " << endl;
            for(size_t k=0;k<N;k++) {
              cout << k << " " << integx[k] << " " << integy[k] << endl;
            }
          }
          return 2;
        }
        
        if (verbose>1) {
          cout << "Found: " << level << endl;
        }
        // Set level from interpolated value
        levs[0]=level;
        
        // End of the fractional integral calculation
      }
      
      // Compute the contour levels
      contour co;
      co.set_levels(levs.size(),levs);
      
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
      
      co.calc_contours(cont_obj);
      
    }

    // Write to a file or store the contours
    if (file.length()>0) {
      if (cont_obj.size()>0) {
	hdf_file hf;
	hf.open_or_create(file);
	hdf_output(hf,cont_obj,name);
	hf.close();
      }
    } else {
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
    
  } else if (type=="prob_dens_mdim_gaussian") {

    if (pdmg_obj.dim()!=2) {
      cerr << "Command contours only works with dim=2." << endl;
      return 1;
    }

    prob_dens_mdim_biv_gaussian<> biv=pdmg_obj.make_biv();
    
  }
  
  return 0;
}

int acol_manager::comm_convert
(std::vector<std::string> &sv, bool itive_com) {

  if (verbose>=2) {
    cng.verbose=verbose;
  } else {
    cng.verbose=0;
  }
  
  if (sv.size()==2 && sv[1]=="list") {
    cng.print_units(std::cout);
    cout << endl;
    cng.print_cache();
    return 0;
  }

  vector<string> in, pr;
  if (sv.size()>=3) {
    for(size_t j=1;j<sv.size();j++) {
      in.push_back(sv[j]);
    }
  } else {
    std::vector<std::string> sv2;
    std::string in2;
    int ret=get_input_one(sv2,"Old unit (or \"add\", \"del\", or \"nat\")",
                          in2,"convert",itive_com);
    if (ret!=0) return ret;
    in.push_back(in2);
    if (in2=="add") {
      ret=get_input_one(sv2,"New unit",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of length",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of mass",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of time",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of temperature",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of current",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of moles",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Power of luminous intensity",in2,
                        "convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Value",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Long name",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
    } else if (in2=="del") {
      ret=get_input_one(sv2,"Unit to delete",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
    } else if (in2=="nat") {
      ret=get_input_one(sv2,"Treat c as 1 (true or false)",
                        in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Treat ħ as 1 (true or false)",
                        in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Treat kb as 1 (true or false)",
                        in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
    } else {
      ret=get_input_one(sv2,"New unit",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
      ret=get_input_one(sv2,"Value",in2,"convert",itive_com);
      if (ret!=0) return ret;
      in.push_back(in2);
    }
  }

  // Set the proper output precision and mode
  if (scientific) cout.setf(ios::scientific);
  else cout.unsetf(ios::scientific);
  cout.precision(precision);
  
  if (in[0]=="add") {

    cout << "Add unit." << endl;

    double val=1.0;
    if (in.size()>=3) {
      int ret2=function_to_double_nothrow(in[9],val,0,&rng);
      if (ret2!=0) {
        cerr << "Converting " << in[9] << " to value failed." << endl;
        return 2;
      }
    }

    convert_units<double>::der_unit d;
    d.label=in[1];
    d.m=o2scl::stod(in[2]);
    d.k=o2scl::stod(in[3]);
    d.s=o2scl::stod(in[4]);
    d.K=o2scl::stod(in[5]);
    d.A=o2scl::stod(in[6]);
    d.mol=o2scl::stod(in[7]);
    d.cd=o2scl::stod(in[8]);
    d.val=val;
    d.name=in[10];

    cng.add_unit(d);

    return 0;
    
  } else if (in[0]=="del") {

    cout << "Delete unit." << endl;

    cng.del_unit(in[2]);

    return 0;
    
  } else if (in[0]=="nat") {

    if (o2scl::stob(in[1])) {
      if (o2scl::stob(in[2])) {
        if (o2scl::stob(in[3])) {
          cout << "Setting c, ħ, and kB to 1." << endl;
        } else {
          cout << "Setting c and ħ to 1." << endl;
        }
      } else {
        if (o2scl::stob(in[3])) {
          cout << "Setting c and kB to 1." << endl;
        } else {
          cout << "Setting c to 1." << endl;
        }
      }
    } else {
      if (o2scl::stob(in[2])) {
        if (o2scl::stob(in[3])) {
          cout << "Setting ħ and kB to 1." << endl;
        } else {
          cout << "Setting ħ to 1." << endl;
        }
      } else {
        if (o2scl::stob(in[3])) {
          cout << "Setting kB to 1." << endl;
        } else {
          cout << "No natural units will be used." << endl;
        }
      }
    }
    
    cng.set_natural_units(o2scl::stob(in[1]),
                          o2scl::stob(in[2]),
                          o2scl::stob(in[3]));
    
    return 0;
    
  }
  
  double val=1.0;
  
  if (in[0]=="2") {
    
    cout << "Using convert()." << endl;
    
    if (in.size()>=4) {
      int ret2=function_to_double_nothrow(in[3],val,0,&rng);
      if (ret2!=0) {
        cerr << "Converting " << in[3] << " to value failed." << endl;
        return 2;
      }
    }
    
    double val_out;
    int cret=cng.convert_ret(in[1],in[2],val,val_out);
    if (cret!=0) {
      cerr << "Conversion failed." << endl;
      return 1;
    }

    cout << val << " " << in[1] << " = " << val_out << " " << in[2] << endl;
    return 0;
    
  } else {
    
    if (in.size()>=3) {
      int ret2=function_to_double_nothrow(in[2],val,0,&rng);
      if (ret2!=0) {
        cerr << "Converting " << in[2] << " to value failed." << endl;
        return 2;
      }
    }

    // If cng.verbose is non-zero, then cng.convert may output
    // verbose information to cout
    double val_out;
    int cret=cng.convert_ret(in[0],in[1],val,val_out);
    if (cret!=0) {
      cerr << "Conversion failed." << endl;
      return 1;
    }
    cout << val << " " << in[0] << " = " << val_out << " " << in[1] << endl;
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

  if (verbose>=3) cng.verbose=2;
  else cng.verbose=0;
  
  ret=table_obj.convert_to_unit(in[0],in[1],false);
  if (ret!=0) {
    cerr << "Could not find column or column does not have unit." << endl;
  }

  return 0;
}

int acol_manager::comm_correl(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()>=3) {
    
    double c=vector_correlation(table_obj.get_nlines(),table_obj[sv[1]],
				table_obj[sv[2]]);
    cout << "Correlation coefficient: " << c << endl;
    
  } else if (sv.size()>=2 && sv[1]==((string)"table3d")) {
      
    table3d_obj.clear();
    size_t n=table_obj.get_ncolumns();
    uniform_grid_end_width<double> ug(0,n-1,1.0);
    table3d_obj.set_xy("x",ug,"y",ug);
    table3d_obj.new_slice("correl");
    
    for(size_t i=0;i<n;i++) {
      for(size_t j=i;j<n;j++) {
        if (i==j) {
          table3d_obj.set(i,j,"correl",1.0);
        } else {
          double c=vector_correlation(table_obj.get_nlines(),
                                      table_obj[i],table_obj[j]);
          if (!std::isfinite(c)) c=0.0;
          table3d_obj.set(i,j,"correl",c);
          table3d_obj.set(j,i,"correl",c);
        }
      }
    }
    
    command_del(type);
    clear_obj();
    command_add("table3d");
    type="table3d";
    
  } else {
      
    vector<string> labels;
    vector<double> coeffs, abs_coeffs;
    
    size_t n=table_obj.get_ncolumns();
    
    for(size_t i=0;i<n;i++) {
      cout << "Computing correlations for column "
           << i+1 << " of " << n << endl;
      for(size_t j=i+1;j<n;j++) {
        labels.push_back(table_obj.get_column_name(i)+", "+
                         table_obj.get_column_name(j));
        double c=vector_correlation(table_obj.get_nlines(),
                                    table_obj[i],table_obj[j]);
        if (!std::isfinite(c)) c=0.0;
        coeffs.push_back(c);
        abs_coeffs.push_back(fabs(c));
        /*
          cout << labels[labels.size()-1] << " " << c << endl;
          char ch;
          cin >> ch;
        */
      }
    }
    cout << endl;
    
    vector<size_t> indexes(coeffs.size());
    vector_sort_index(coeffs.size(),abs_coeffs,indexes);
    
    for(size_t j=0;j<coeffs.size();j++) {
      size_t k=indexes[coeffs.size()-1-j];
      cout << j << " ";
      cout << labels[k] << " "
           << coeffs[k] << " " << abs_coeffs[k] << endl;;
    }
      
  }
  
  return 0;
}
  
int acol_manager::comm_create(std::vector<std::string> &sv, bool itive_com) {
  std::string ctype, tval;

  // Delete previous object
  command_del(type);
  clear_obj();
  
  int retx=get_input_one(sv,"Enter type of object to create",ctype,"create",
			itive_com);
  if (retx!=0) return retx;

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
    
  } else if (ctype=="string[]") {

    int ret=get_input_one(sv2,"Enter string spec.",tval,"create",
			  itive_com);
    if (ret!=0) return ret;
    
    int ssret=strings_spec(tval,stringv_obj,verbose,false);
    if (ssret!=0) {
      cerr << "Strings specification failed." << endl;
    } else {
      type="string[]";
      command_add("string[]");
      obj_name="string[]";
    }
    
  } else if (ctype=="vec_vec_double") {

    int ret=get_input_one(sv2,"Enter mult-vector spec.",tval,"create",
			  itive_com);
    if (ret!=0) return ret;

    int mvs_ret=mult_vector_spec(tval,vvdouble_obj,false,verbose,false);
    if (mvs_ret!=0) {
      cerr << "Multiple vector specification failed." << endl;
    } else {
      type="vec_vec_double";
      command_add("vec_vec_double");
      obj_name="vec_vec_double";

      for(size_t j=2;j<sv2.size();j++) {
        vector<vector<double>> vvd2;
        mvs_ret=mult_vector_spec(sv2[j],vvd2,false,verbose,false);
        if (mvs_ret==0) {
          for(size_t k=0;k<vvd2.size();k++) {
            vvdouble_obj.push_back(vvd2[k]);
          }
        }
      }
    }
    
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
      vector<string>::iterator itx=sv3.begin();
      sv3.erase(itx+1);
      
      std::string in2;

      int ret2=get_input_one(sv3,"Function of i (starting with zero)",
			     in2,"create",itive_com);
      if (ret2!=0) return ret2;

      calc_utf8<> calc;
      calc.set_rng(rng);
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

      int ret2=vector_spec(in1,doublev_obj,false,verbose,false);
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

    calc_utf8<> calc;
    calc.set_rng(rng);
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

    calc_utf8<> calc;
    calc.set_rng(rng);
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
    int vs_ret=vector_spec(in[1],d,false,verbose,false);
    if (vs_ret!=0) {
      cerr << "Vector specification " << in[1] << " failed." << endl;
      return 1;
    }
    
    table_obj.clear();
    table_obj.line_of_names(in[0]);
    table_obj.set_nlines(d.size());
    
    for(size_t li=0;li<d.size();li++) {
      table_obj.set(in[0],li,d[li]);
    }
    command_add("table");
    type="table";
    
  } else if (ctype=="table-mv") {
    
    vector<string> in, pr;
    if (sv2.size()<2) {
      pr.push_back("Specification of column names");
      pr.push_back("Multiple vector specification for column data");
      int ret=get_input(sv2,pr,in,"create",itive_com);
      if (ret!=0) return ret;
    } else {
      for(size_t j=1;j<sv2.size();j++) {
	in.push_back(sv2[j]);
      }
    }

    table_obj.clear();
    
    size_t nk=sv2.size()/2;
    
    if (verbose>2) {
      std::cout << "nk: " << nk << std::endl;
    }
    
    vector<string> c;
    vector<std::vector<double> > d;

    for(size_t ik=0;ik<nk;ik++) {

      int ss_ret=strings_spec(in[ik*2],c,verbose,false);
      if (ss_ret!=0) {
	cerr << "String specification " << in[ik*2] << " failed." << endl;
	return 1;
      }

      int vs_ret=mult_vector_spec(in[ik*2+1],d,false,verbose,false);
      if (vs_ret!=0) {
	cerr << "Multiple vector specification "
	     << in[ik*2+1] << " failed." << endl;
	return 2;
      }
    }

    if (c.size()!=d.size()) {
      cerr << "Mismatch between number of column names and "
	   << "number of data vectors." << endl;
      return 3;
    }

    for(size_t i=0;i<c.size();i++) {
      table_obj.new_column(c[i]);
    }

    size_t max_size=0;
    for(size_t i=0;i<d.size();i++) {
      if (d[i].size()>max_size) max_size=d[i].size();
    }

    if (max_size==0) {
      cerr << "Data vectors all have size 0." << endl;
      return 4;
    }
    
    table_obj.set_nlines(max_size);

    for(size_t i=0;i<d.size();i++) {
      
      if (d[i].size()==max_size) {
	// If the vector has enough data, then copy it into
	// the new table
	for(size_t li=0;li<d[i].size();li++) {
	  table_obj.set(i,li,d[i][li]);
	}
      } else {
	// If the vector size doesn't match, use the internal
	// interpolation type to expand the vector to fit in the table
	vector_index_vector<double> viv;
	interp_vec<vector_index_vector<double>,vector<double> > iv
	  (d[i].size(),viv,d[i],interp_type);
	for(size_t li=0;li<max_size;li++) {
	  double dval=((double)li)/((double)(max_size-1))*
	    ((double)(d[i].size()-1));
	  table_obj.set(i,li,iv.eval(dval));
	}
      }
    }
    
    command_add("table");
    type="table";
    
  } else if (ctype=="tensor") {

    std::string i1;
    int rety=get_input_one(sv2,"Enter rank",i1,"create",itive_com);
    if (rety!=0) return rety;
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
    int retz=get_input_one(sv2,"Enter rank",i1,"create",itive_com);
    if (retz!=0) return retz;
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
    int vs_ret=vector_spec(in[1],dx,false,verbose,false);
    if (vs_ret!=0) {
      cerr << "Vector specification " << in[1] << " failed." << endl;
    }

    std::string yname=in[2];
    vs_ret=vector_spec(in[3],dy,false,verbose,false);
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
