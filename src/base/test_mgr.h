/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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

#include <o2scl/string_conv.h>
#include <o2scl/misc.h>

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
	\mathrm{expected}<\mathrm{rel\_error}\f$
    */
    template<class data_t>
      bool test_rel(data_t result, data_t expected, data_t rel_error,
		    std::string description) {
      bool ret;
      if (std::isnan(expected)) {
	ret=(std::isnan(expected)==std::isnan(result));
	description=dtos(result)+" vs. "+ dtos(expected)+
	  "\n "+description;
      } else if (std::isinf(expected)) {
	ret=(std::isinf(expected)==std::isinf(result));
	description=dtos(result)+" vs. "+ dtos(expected)+
	  "\n "+description;
      } else if (expected==0.0) {
	ret=test_abs(result,expected,rel_error,description);
	return ret;
      } else {
	ret=((fabs(expected-result))/fabs(expected)<rel_error);	
	description=dtos(result)+" vs. "+dtos(expected)+
          " is "+dtos(fabs(expected-result)/fabs(expected))+
    	  " > "+dtos(rel_error)+"\n "+description;
      }
      
      process_test(ret,"relative",description);
      return ret;
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	<\mathrm{abs\_error}\f$
    */
    template<class data_t>
      bool test_abs(data_t result, data_t expected, data_t abs_error,
		    std::string description) {
      bool ret;
      if (std::isnan(expected)) {
	ret=(std::isnan(expected)==std::isnan(result));
	description=dtos(result)+" vs. "+ dtos(expected)+
	  "\n "+description;
      } else if (std::isinf(expected)) {
	ret=(std::isinf(expected)==std::isinf(result));
	description=dtos(result)+" vs. "+ dtos(expected)+
	  "\n "+description;
      } else {
	ret=(fabs(expected-result)<abs_error);
	description=dtos(result)+" vs. "+ dtos(expected)+" is "
	  +dtos(fabs(expected-result))+" > "+dtos(abs_error)+
	  "\n "+description;
      }
  
      process_test(ret,"absolute",description);

      return ret;
    }

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
	if (o2scl::is_nan(expected[i])) {
	  ret=(ret && (o2scl::is_nan(expected[i])==o2scl::is_nan(result[i])));
	} else if (o2scl::is_inf(expected[i])) {
	  ret=(ret && (o2scl::is_inf(expected[i])==o2scl::is_inf(result[i])));
	} else if (expected[i]==0.0) {
	  ret=(ret && test_abs(result[i],expected[i],rel_error,description));
	  if (fabs(result[i]-expected[i])>max) {
	    max=fabs(result[i]-expected[i]);
	  }
	} else {
	  ret=(ret && ((fabs(expected[i]-result[i]))/
		       fabs(expected[i])<rel_error));
	  if (fabs(expected[i]-result[i])/fabs(expected[i])>max) {
	    max=fabs(expected[i]-result[i])/fabs(expected[i]);
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
	if (o2scl::is_nan(expected[i])) {
	  ret=(ret && (o2scl::is_nan(expected[i])==o2scl::is_nan(result[i])));
	} else if (o2scl::is_inf(expected[i])) {
	  ret=(ret && (o2scl::is_inf(expected[i])==o2scl::is_inf(result[i])));
	} else {
	  ret=(ret && (fabs(expected[i]-result[i])<abs_error));
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
	if (o2scl::is_nan(expected[i])) {
	  ret=(ret && (o2scl::is_nan(expected[i])==o2scl::is_nan(result[i])));
	} else if (o2scl::is_inf(expected[i])) {
	  ret=(ret && (o2scl::is_inf(expected[i])==o2scl::is_inf(result[i])));
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
	of an array
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
	  if (o2scl::is_nan(expected(i,j))) {
	    ret=(ret && (o2scl::is_nan(expected(i,j))==
			 o2scl::is_nan(result(i,j))));
	  } else if (o2scl::is_inf(expected(i,j))) {
	    ret=(ret && (o2scl::is_inf(expected(i,j))==
			 o2scl::is_inf(result(i,j))));
	  } else if (expected(i,j)==0.0) {
	    ret=(ret && test_abs(result(i,j),expected(i,j),rel_error,
				 description));
	    if (fabs(result(i,j)-expected(i,j))>max) {
	      max=fabs(result(i,j)-expected(i,j));
	    }
	  } else {
	    ret=(ret && ((fabs(expected(i,j)-result(i,j)))/
			 fabs(expected(i,j))<rel_error));
	    if (fabs(expected(i,j)-result(i,j))/fabs(expected(i,j))>max) {
	      max=fabs(expected(i,j)-result(i,j))/fabs(expected(i,j));
	    }
	  }
	}
      }
      
      description=((std::string)"max=")+o2scl::dtos(max)+
	"\n "+description;
      process_test(ret,"relative matrix",description);
      
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
	  if (o2scl::is_nan(expected(i,j))) {
	    ret=(ret && (o2scl::is_nan(expected(i,j))==
			 o2scl::is_nan(result(i,j))));
	  } else if (o2scl::is_inf(expected(i,j))) {
	    ret=(ret && (o2scl::is_inf(expected(i,j))==
			 o2scl::is_inf(result(i,j))));
	  } else if (expected(i,j)==0.0) {
	    ret=(ret && test_abs(result(i,j),expected(i,j),abs_error,
				 description));
	    if (fabs(result(i,j)-expected(i,j))>max) {
	      max=fabs(result(i,j)-expected(i,j));
	    }
	  } else {
	    ret=(ret && ((fabs(expected(i,j)-result(i,j)))<abs_error));
	    if (fabs(expected(i,j)-result(i,j))>max) {
	      max=fabs(expected(i,j)-result(i,j));
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
