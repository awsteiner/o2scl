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
#ifndef O2SCL_TEST_MGR_H
#define O2SCL_TEST_MGR_H

/** \file test_mgr.h
    \brief File defining \ref o2scl::test_mgr
*/

#include <string>

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/number.hpp>
#endif

#include <o2scl/string_conv.h>
#include <o2scl/misc.h>
#include <o2scl/table_units.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_matrix.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief A class to manage testing and record success and failure
   */
  class test_mgr {

  protected:

#ifndef DOXYGEN_INTERNAL

    /// The number of tests performed
    int ntests;

    /// The output level
    int output_level;
  
    /// A helper function for processing tests
    void process_test(bool ret, std::string d2, std::string description);

    /// True if all tests have passed
    bool success;

    /// The description of the last failed test
    std::string last_fail;

#endif
  
  public:

    /// Create a \ref test_mgr object
    test_mgr(bool success_l=true, std::string last_fail_l="", 
	     int ntests_l=0, int output_level_l=1) {
      success=success_l;
      last_fail=last_fail_l;
      ntests=ntests_l;
      output_level=output_level_l;
    }

    /** \brief Provide a report of all tests so far.

	This function reports on whether or not all tests have passed
	according to the current output level. It returns true if all
	tests have passed and false if at least one test failed.
    */
    bool report() const;

    /// \name Individual get and set methods
    //@{
    /// Return true if all tests have succeeded
    bool get_success() const {
      return success;
    }

    /// Return the last failure description
    std::string get_last_fail() const {
      return last_fail;
    }

    /// Return the output level
    int get_output_level() const {
      return output_level;
    }

    /// Returns the description of the last test that failed.
    std::string get_last_fail() { return last_fail; };

    /** \brief Set the output level

	Possible values:
	- 0 = No output
	- 1 = Output only tests that fail
	- 2 = Output all tests
    */
    void set_output_level(int l) { output_level=l; };

    /// Return the number of tests performed so far
    int get_ntests() const { return ntests; };
    //@}

    /// \name Main testing methods
    //@{
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	<\mathrm{abs\_error}\f$
    */
    template<class data_t=double>
      bool test_abs(data_t result, data_t expected, data_t abs_error,
		    std::string description) {
      bool ret;
      if (std::isnan(expected)) {
	ret=(std::isnan(expected)==std::isnan(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	  "\n "+description;
      } else if (std::isinf(expected)) {
	ret=(std::isinf(expected)==std::isinf(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	"\n "+description;
      } else {
	ret=(std::abs(expected-result)<abs_error);
	if (ret) {
	  description=dtos<data_t>(result)+" vs. "+
	    dtos<data_t>(expected)+" : "
	    +dtos<data_t>(std::abs(expected-result))+" < "+
	    dtos<data_t>(abs_error)+
	    "\n "+description;
	} else {
	  description=dtos<data_t>(result)+" vs. "+
	  dtos<data_t>(expected)+" : "
	  +dtos<data_t>(std::abs(expected-result))+" > "+
	  dtos<data_t>(abs_error)+
	  "\n "+description;
	}
      }
  
      process_test(ret,"absolute",description);

      return ret;
    }

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$
    */
    template<class data_t=double>
      bool test_rel(data_t result, data_t expected, data_t rel_error,
		    std::string description) {
      bool ret;
      if (std::isnan(expected)) {
	ret=(std::isnan(expected)==std::isnan(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	  "\n "+description;
      } else if (std::isinf(expected)) {
	ret=(std::isinf(expected)==std::isinf(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	"\n "+description;
      } else if (expected==0.0) {
	ret=test_abs<data_t>(result,expected,rel_error,description);
	return ret;
      } else {
	double obs_err=std::abs(expected-result)/std::abs(expected);
	ret=(obs_err<rel_error);
	if (ret) {
	  description=dtos<data_t>(result)+" vs. "+dtos<data_t>(expected)+
	    " : "+dtos<data_t>(obs_err)+
	    " < "+dtos<data_t>(rel_error)+"\n "+description;
	} else {
	  description=dtos<data_t>(result)+" vs. "+dtos<data_t>(expected)+
	  " : "+dtos<data_t>(obs_err)+
	  " > "+dtos<data_t>(rel_error)+"\n "+description;
	}
      }
      
      process_test(ret,"relative",description);
      return ret;
    }
    
#if defined(O2SCL_LD_TYPES) || defined(DOXYGEN)
    
    /** \brief Testing functions for \c boost::multiprecision numbers

	This function is similar to \ref test_abs(), but replaces 
	\c isnan and related functions with boost versions from
	\c boost::multiprecision::number .
     */
    template<class data_t=double>
      bool test_abs_boost(data_t result, data_t expected, data_t abs_error,
			  std::string description) {
      bool ret;
      if (boost::multiprecision::isnan(expected)) {
	ret=(boost::multiprecision::isnan(expected)==
	     boost::multiprecision::isnan(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	  "\n "+description;
      } else if (boost::multiprecision::isinf(expected)) {
	ret=(boost::multiprecision::isinf(expected)==
	     boost::multiprecision::isinf(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	"\n "+description;
      } else {
	ret=(boost::multiprecision::abs(expected-result)<abs_error);
	if (ret) {
	  description=dtos<data_t>(result)+" vs. "+
	    dtos<data_t>(expected)+" : "
	    +dtos<data_t>(boost::multiprecision::abs(expected-result))+" < "+
	    dtos<data_t>(abs_error)+
	    "\n "+description;
	} else {
	  description=dtos<data_t>(result)+" vs. "+
	  dtos<data_t>(expected)+" : "
	  +dtos<data_t>(boost::multiprecision::abs(expected-result))+" > "+
	  dtos<data_t>(abs_error)+
	  "\n "+description;
	}
      }
  
      process_test(ret,"absolute",description);

      return ret;
    }

    /** \brief Testing functions for \c boost::multiprecision numbers

	This function is similar to \ref test_abs(), but replaces 
	\c isnan and related functions with boost versions from
	\c boost::multiprecision::number .
    */
    template<class data_t=double>
      bool test_rel_boost(data_t result, data_t expected, data_t rel_error,
			  std::string description) {
      bool ret;
      if (boost::multiprecision::isnan(expected)) {
	ret=(boost::multiprecision::isnan(expected)==
	     boost::multiprecision::isnan(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	  "\n "+description;
      } else if (boost::multiprecision::isinf(expected)) {
	ret=(boost::multiprecision::isinf(expected)==
	     boost::multiprecision::isinf(result));
	description=dtos<data_t>(result)+" vs. "+ dtos<data_t>(expected)+
	"\n "+description;
      } else if (expected==0.0) {
	ret=test_abs_boost<data_t>(result,expected,rel_error,description);
	return ret;
      } else {
	double obs_err=boost::multiprecision::abs(expected-result)/boost::multiprecision::abs(expected);
	ret=(obs_err<rel_error);
	if (ret) {
	  description=dtos<data_t>(result)+" vs. "+dtos<data_t>(expected)+
	    " : "+dtos<data_t>(obs_err)+
	    " < "+dtos<data_t>(rel_error)+"\n "+description;
	} else {
	  description=dtos<data_t>(result)+" vs. "+dtos<data_t>(expected)+
	  " : "+dtos<data_t>(obs_err)+
	  " > "+dtos<data_t>(rel_error)+"\n "+description;
	}
      }
      
      process_test(ret,"relative",description);
      return ret;
    }
    
#endif

    /** \brief  Test for \f$1/\mathrm{factor} < \mathrm{result/expected} 
	< \mathrm{factor}\f$
    */
    template<class data_t>
      bool test_fact(data_t result, data_t expected, data_t factor,
		     std::string description) {
      bool ret;
      double ratio;
      if (std::isnan(expected)) {
	ret=(std::isnan(expected)==std::isnan(result));
      } else if (std::isinf(expected)) {
	ret=(std::isinf(expected)==std::isinf(result));
      } else {
	ratio=expected/result;
	ret=(ratio<factor && ratio>1.0/factor);
      }

      description= dtos(result)+" vs. "+ dtos(expected)+"\n "+
	description;
      process_test(ret,"factor",description);

      return ret;
    }

    /// Test for \f$\mathrm{result}=\mathrm{expected}\f$
    bool test_str(std::string result, std::string expected, 
		  std::string description);

    /// Test for \f$\mathrm{result}=\mathrm{expected}\f$
    bool test_gen(bool value, std::string description);
    //@}

    /// \name Vector testing methods
    //@{
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class vec_t, class vec2_t, class data_t>
      bool test_rel_vec(int nv, const vec_t &result, const vec2_t &expected, 
			data_t rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i;
  
      for(i=0;i<nv;i++) {
	if (std::isnan(expected[i])) {
	  ret=(ret && (std::isnan(expected[i])==std::isnan(result[i])));
	} else if (std::isinf(expected[i])) {
	  ret=(ret && (std::isinf(expected[i])==std::isinf(result[i])));
	} else if (expected[i]==0.0) {
	  ret=(ret && test_abs(result[i],expected[i],rel_error,description));
	  if (std::abs(result[i]-expected[i])>max) {
	    max=std::abs(result[i]-expected[i]);
	  }
	} else {
	  ret=(ret && ((std::abs(expected[i]-result[i]))/
		       std::abs(expected[i])<rel_error));
	  if (std::abs(expected[i]-result[i])/std::abs(expected[i])>max) {
	    max=std::abs(expected[i]-result[i])/std::abs(expected[i]);
	  }
	}
      }
      
      description=((std::string)"max=")+o2scl::dtos(max)+
	"\n "+description;
      process_test(ret,"relative array",description);

      return ret;
    }

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	<\mathrm{abs\_error}\f$ over each element
	of an array
    */
    template<class vec_t, class vec2_t, class data_t>
      bool test_abs_vec(int nv, const vec_t &result, const vec2_t &expected, 
			data_t abs_error, std::string description) {
      bool ret=true;
      int i;
  
      for(i=0;i<nv;i++) {
	if (std::isnan(expected[i])) {
	  ret=(ret && (std::isnan(expected[i])==std::isnan(result[i])));
	} else if (std::isinf(expected[i])) {
	  ret=(ret && (std::isinf(expected[i])==std::isinf(result[i])));
	} else {
	  ret=(ret && (std::abs(expected[i]-result[i])<abs_error));
	}
      }
  
      description="\n "+description;
      process_test(ret,"absolute array",description);
  
      return ret;
    }

    /** \brief Test for \f$ 1/factor < result/expected < factor \f$ 
	over each element of an array
    */
    template<class vec_t, class vec2_t, class data_t>
      bool test_fact_vec(int nv, const vec_t &result, const vec2_t &expected, 
			 data_t factor, std::string description) {
      bool ret=true;
      int i;
      double ratio;
  
      for(i=0;i<nv;i++) {
	if (std::isnan(expected[i])) {
	  ret=(ret && (std::isnan(expected[i])==std::isnan(result[i])));
	} else if (std::isinf(expected[i])) {
	  ret=(ret && (std::isinf(expected[i])==std::isinf(result[i])));
	} else {
	  ratio=expected[i]/result[i];
	  ret=(ret && (ratio<factor && ratio>1.0/factor));
	}
      }
  
      description="\n "+description;
      process_test(ret,"factor array",description);
  
      return ret;
    }
    
    /// Test for equality of a generic array
    template<class vec_t>
      bool test_gen_vec(int nv, const vec_t &result, const vec_t &expected, 
			std::string description) {
      bool ret=true;
      int i;
      
      for(i=0;i<nv;i++) {
	ret=(ret && (result[i]==expected[i]));
      }
      
      description="\n "+description;
      process_test(ret,"generic array",description);
      
      return ret;
    }
    //@}

    /// \name Matrix testing methods
    //@{
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	in a matrix
    */
    template<class mat_t, class mat2_t, class data_t>
      bool test_rel_mat(int nr, int nc, const mat_t &result, 
			const mat2_t &expected, 
			data_t rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
  
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (std::isnan(expected(i,j))) {
	    ret=(ret && (std::isnan(expected(i,j))==
			 std::isnan(result(i,j))));
	  } else if (std::isinf(expected(i,j))) {
	    ret=(ret && (std::isinf(expected(i,j))==
			 std::isinf(result(i,j))));
	  } else if (expected(i,j)==0.0) {
	    ret=(ret && test_abs(result(i,j),expected(i,j),rel_error,
				 description));
	    if (std::abs(result(i,j)-expected(i,j))>max) {
	      max=std::abs(result(i,j)-expected(i,j));
	    }
	  } else {
	    ret=(ret && ((std::abs(expected(i,j)-result(i,j)))/
			 std::abs(expected(i,j))<rel_error));
	    if (std::abs(expected(i,j)-result(i,j))/
		std::abs(expected(i,j))>max) {
	      max=std::abs(expected(i,j)-result(i,j))/
		std::abs(expected(i,j));
	    }
	  }
	}
      }
      
      description=((std::string)"max=")+o2scl::dtos(max)+
	"\n "+description;
      process_test(ret,"relative matrix",description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	in a matrix larger than a specified tolerance
    */
    template<class mat_t, class mat2_t, class data_t>
      bool test_rel_nonzero_mat(int nr, int nc, const mat_t &result, 
				const mat2_t &expected, 
				data_t error, data_t zero_tol,
				std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
  
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (std::isnan(expected(i,j))) {
	    ret=(ret && (std::isnan(expected(i,j))==
			 std::isnan(result(i,j))));
	  } else if (std::isinf(expected(i,j))) {
	    ret=(ret && (std::isinf(expected(i,j))==
			 std::isinf(result(i,j))));
	  } else if (expected(i,j)<zero_tol) {
	    ret=(ret && test_abs(result(i,j),expected(i,j),error,
				 description));
	    if (std::abs(result(i,j)-expected(i,j))>max) {
	      max=std::abs(result(i,j)-expected(i,j));
	    }
	  } else {
	    ret=(ret && ((std::abs(expected(i,j)-result(i,j)))/
			 std::abs(expected(i,j))<error));
	    if (std::abs(expected(i,j)-result(i,j))/
		std::abs(expected(i,j))>max) {
	      max=std::abs(expected(i,j)-result(i,j))/
		std::abs(expected(i,j));
	    }
	  }
	}
      }
      
      description=((std::string)"max=")+o2scl::dtos(max)+
	"\n "+description;
      process_test(ret,"rel_nonzero matrix",description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}| <
	\mathrm{abs\_error} \f$ over each element in a matrix
    */
    template<class mat_t, class mat2_t, class data_t>
      bool test_abs_mat(int nr, int nc, const mat_t &result, 
			const mat2_t &expected, data_t abs_error, 
			std::string description) {
			
      bool ret=true;
      double max=0.0;
      int i, j;
      
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (std::isnan(expected(i,j))) {
	    ret=(ret && (std::isnan(expected(i,j))==
			 std::isnan(result(i,j))));
	  } else if (std::isinf(expected(i,j))) {
	    ret=(ret && (std::isinf(expected(i,j))==
			 std::isinf(result(i,j))));
	  } else if (expected(i,j)==0.0) {
	    ret=(ret && test_abs(result(i,j),expected(i,j),abs_error,
				 description));
	    if (std::abs(result(i,j)-expected(i,j))>max) {
	      max=std::abs(result(i,j)-expected(i,j));
	    }
	  } else {
	    ret=(ret && ((std::abs(expected(i,j)-result(i,j)))<abs_error));
	    if (std::abs(expected(i,j)-result(i,j))>max) {
	      max=std::abs(expected(i,j)-result(i,j));
	    }
	  }
	}
      }
      
      description=((std::string)"max=")+o2scl::dtos(max)+
	"\n "+description;
      process_test(ret,"absolute matrix",description);
      
      return ret;
      
    }
    //@}

    /** \brief Compare entries in \c expected to see if they match
	those in table \c result.

	If the numbers in the \c expected table have an absolute value
	less than \c zero_tol, then the absolute value of the
	difference is used for the comparison. Otherwise, the absolute
	value of the relative difference is used to make the
	comparison.
    */
    template<class vec_t, class data_t>
      bool test_rel_nonzero_table(const table_units<vec_t> &result,
				  const table_units<vec_t> &expected,
				  data_t error, data_t zero_tol,
				  std::string description) {
      bool ret=true;
      int i, j;
      int nr=result.get_nlines();
      int nc=result.get_ncolumns();
      std::vector<double> max(nc);
      
      for(i=0;i<nc;i++) {
	std::string col_name=expected.get_column_name(i);
	for(j=0;j<nr;j++) {
	  std::string desc1=description+" col: "+col_name+" row: "+
	    o2scl::itos(j);
	  if (std::isnan(expected.get(i,j))) {
	    bool ret1=test_gen(std::isnan(expected.get(i,j))==
			       std::isnan(result.get(i,j)),desc1);
	    ret=(ret && ret1);
	  } else if (std::isinf(expected.get(i,j))) {
	    bool ret1=test_gen(std::isinf(expected.get(i,j))==
			       std::isinf(result.get(i,j)),desc1);
	    ret=(ret && ret1);
	  } else if (expected.get(i,j)<zero_tol) {
	    bool ret1=test_abs(result.get(i,j),expected.get(i,j),error,
			       desc1);
	    if (std::abs(result.get(i,j)-expected.get(i,j))>max[i]) {
	      max[i]=std::abs(result.get(i,j)-expected.get(i,j));
	    }
	  } else {
	    bool ret1=test_rel(result.get(i,j),expected.get(i,j),error,
			       desc1);
	    ret=(ret && ret1);
	    if (std::abs(expected.get(i,j)-result.get(i,j))/
		std::abs(expected.get(i,j))>max[i]) {
	      max[i]=std::abs(expected.get(i,j)-result.get(i,j))/
		std::abs(expected.get(i,j));
	    }
	  }
	}
      }

      if ((output_level>=1 && ret==false) || output_level>=2) {
	for(i=0;i<nc;i++) {
	  std::cout << "Max for " << expected.get_column_name(i) << " is "
		    << max[i] << std::endl;
	}
      }
      process_test(ret,"rel_nonzero table",description);
      
      return ret;
      
    }

    /** \brief Add two test_mgr objects (if either failed, the sum fails)

	The output level is set to the maximum value of left and
	right operand and the number of tests is set equal to
	the sum. The last failure descriptions of both operands
	are appended with a <tt>operator+()</tt> prefix, or blank
	if there were no failures from either. 
    */
    friend const test_mgr operator+(const test_mgr& left, 
				    const test_mgr& right);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
