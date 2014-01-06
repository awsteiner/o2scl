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

      \future test_mgr::success and test_mgr::last_fail should be protected,
      but that breaks the operator+() function. Can this be fixed?
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

#endif
  
  public:

    test_mgr() {
      success=true; 
      ntests=0; 
      output_level=1;
    }

    /** \brief Provide a report of all tests so far.

	Returns true if all tests have passed and false if at least
	one test failed.
    */
    bool report();

    /// Returns the description of the last test that failed.
    std::string get_last_fail() {return last_fail;};

    /** \brief Set the output level

	Possible values:
	- 0 = No output
	- 1 = Output only tests that fail
	- 2 = Output all tests
    */
    void set_output_level(int l) { output_level=l; };

    /// Return the number of tests performed so far
    int get_ntests() {return ntests;};

    /// \name The testing methods
    ///@{

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$
    */
    bool test_rel(double result, double expected, double rel_error,
		  std::string description);

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	<\mathrm{abs\_error}\f$
    */
    bool test_abs(double result, double expected, double abs_error,
		  std::string description);

    /** \brief  Test for \f$1/\mathrm{factor} < \mathrm{result/expected} 
	< \mathrm{factor}\f$
    */
    bool test_fact(double result, double expected, double factor,
		   std::string description);

    /// Test for \f$\mathrm{result}=\mathrm{expected}\f$
    bool test_str(std::string result, std::string expected, 
		  std::string description);

    /// Test for \f$\mathrm{result}=\mathrm{expected}\f$
    bool test_gen(bool value, std::string description);

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class vec_t, class vec2_t>
      bool test_rel_arr(int nv, const vec_t &result, const vec2_t &expected, 
			double rel_error, std::string description) {
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
      
      process_test(ret,((std::string)"relative array, max=")+o2scl::dtos(max),
		   description);
      
      return ret;
      
    }

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    bool test_rel_arrgslgsl(int nv, gsl_vector *result, gsl_vector *expected, 
			 double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i;
      
      for(i=0;i<nv;i++) {
	if (o2scl::is_nan(gsl_vector_get(expected,i))) {
	  ret=(ret && (o2scl::is_nan(gsl_vector_get(expected,i))==
		       o2scl::is_nan(gsl_vector_get(result,i))));
	} else if (o2scl::is_inf(gsl_vector_get(expected,i))) {
	  ret=(ret && (o2scl::is_inf(gsl_vector_get(expected,i))==
		       o2scl::is_inf(gsl_vector_get(result,i))));
	} else if (gsl_vector_get(expected,i)==0.0) {
	  ret=(ret && test_abs(gsl_vector_get(result,i),
			       gsl_vector_get(expected,i),
			       rel_error,description));
	  if (fabs(gsl_vector_get(result,i)-
		   gsl_vector_get(expected,i))>max) {
	    max=fabs(gsl_vector_get(result,i)-gsl_vector_get(expected,i));
	  }
	} else {
	  ret=(ret && ((fabs(gsl_vector_get(expected,i)-
			     gsl_vector_get(result,i)))/
		       fabs(gsl_vector_get(expected,i))<rel_error));
	  if (fabs(gsl_vector_get(expected,i)-
		   gsl_vector_get(result,i))/
	      fabs(gsl_vector_get(expected,i))>max) {
	    max=fabs(gsl_vector_get(expected,i)-
		     gsl_vector_get(result,i))/
	      fabs(gsl_vector_get(expected,i));
	  }
	}
      }
      
      process_test(ret,((std::string)"relative array, max=")+dtos(max),
		   description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class vec_t>
      bool test_rel_arrgsl(int nv, const vec_t &result, gsl_vector *expected, 
			   double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i;
      
      for(i=0;i<nv;i++) {
	if (o2scl::is_nan(gsl_vector_get(expected,i))) {
	  ret=(ret && (o2scl::is_nan(gsl_vector_get(expected,i))==
		       o2scl::is_nan(result[i])));
	} else if (o2scl::is_inf(gsl_vector_get(expected,i))) {
	  ret=(ret && (o2scl::is_inf(gsl_vector_get(expected,i))==
		       o2scl::is_inf(result[i])));
	} else if (gsl_vector_get(expected,i)==0.0) {
	  ret=(ret && test_abs(result[i],gsl_vector_get(expected,i),
			       rel_error,description));
	  if (fabs(result[i]-gsl_vector_get(expected,i))>max) {
	    max=fabs(result[i]-gsl_vector_get(expected,i));
	  }
	} else {
	  ret=(ret && ((fabs(gsl_vector_get(expected,i)-result[i]))/
		       fabs(gsl_vector_get(expected,i))<rel_error));
	  if (fabs(gsl_vector_get(expected,i)-result[i])/
	      fabs(gsl_vector_get(expected,i))>max) {
	    max=fabs(gsl_vector_get(expected,i)-result[i])/
	      fabs(gsl_vector_get(expected,i));
	  }
	}
      }
      
      process_test(ret,((std::string)"relative array, max=")+dtos(max),
		   description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    bool test_rel_matgslgsl(int nr, int nc, gsl_matrix *result, 
			      gsl_matrix *expected, 
			      double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
      
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (o2scl::is_nan(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_nan(gsl_matrix_get(expected,i,j))==
			 o2scl::is_nan(gsl_matrix_get(result,i,j))));
	  } else if (o2scl::is_inf(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_inf(gsl_matrix_get(expected,i,j))==
			 o2scl::is_inf(gsl_matrix_get(result,i,j))));
	  } else if (gsl_matrix_get(expected,i,j)==0.0) {
	    ret=(ret && test_abs(gsl_matrix_get(result,i,j),
				 gsl_matrix_get(expected,i,j),rel_error,
				 description));
	    if (fabs(gsl_matrix_get(result,i,j)-
		     gsl_matrix_get(expected,i,j))>max) {
	      max=fabs(gsl_matrix_get(result,i,j)-
		       gsl_matrix_get(expected,i,j));
	    }
	  } else {
	    ret=(ret && ((fabs(gsl_matrix_get(expected,i,j)-
			       gsl_matrix_get(result,i,j)))/
			 fabs(gsl_matrix_get(expected,i,j))<rel_error));
	    if (fabs(gsl_matrix_get(expected,i,j)-
		     gsl_matrix_get(result,i,j))/
		fabs(gsl_matrix_get(expected,i,j))>max) {
	      max=fabs(gsl_matrix_get(expected,i,j)-
		       gsl_matrix_get(result,i,j))/
		fabs(gsl_matrix_get(expected,i,j));
	    }
	  }
	}
      }
      
      process_test(ret,((std::string)"relative matrix, max=")+dtos(max),
		   description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class mat_t>
      bool test_rel_matgsl(int nr, int nc, const mat_t &result, 
			   gsl_matrix *expected, 
			   double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
  
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (o2scl::is_nan(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_nan(gsl_matrix_get(expected,i,j))==
			 o2scl::is_nan(result(i,j))));
	  } else if (o2scl::is_inf(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_inf(gsl_matrix_get(expected,i,j))==
			 o2scl::is_inf(result(i,j))));
	  } else if (gsl_matrix_get(expected,i,j)==0.0) {
	    ret=(ret && test_abs(result(i,j),
				 gsl_matrix_get(expected,i,j),rel_error,
				 description));
	    if (fabs(result(i,j)-gsl_matrix_get(expected,i,j))>max) {
	      max=fabs(result(i,j)-gsl_matrix_get(expected,i,j));
	    }
	  } else {
	    ret=(ret && ((fabs(gsl_matrix_get(expected,i,j)-result(i,j)))/
			 fabs(gsl_matrix_get(expected,i,j))<rel_error));
	    if (fabs(gsl_matrix_get(expected,i,j)-
		     result(i,j))/fabs(gsl_matrix_get(expected,i,j))>max) {
	      max=fabs(gsl_matrix_get(expected,i,j)-
		       result(i,j))/fabs(gsl_matrix_get(expected,i,j));
	    }
	  }
	}
      }
      
      process_test(ret,((std::string)"relative matrix, max=")+o2scl::dtos(max),
		   description);
      
      return ret;
      
    }

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|
	<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class mat_t>
      bool test_abs_matgsl(int nr, int nc, const mat_t &result, 
			   gsl_matrix *expected, 
			   double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
  
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (o2scl::is_nan(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_nan(gsl_matrix_get(expected,i,j))==
			 o2scl::is_nan(result(i,j))));
	  } else if (o2scl::is_inf(gsl_matrix_get(expected,i,j))) {
	    ret=(ret && (o2scl::is_inf(gsl_matrix_get(expected,i,j))==
			 o2scl::is_inf(result(i,j))));
	  } else {
	    ret=(ret && ((fabs(gsl_matrix_get(expected,i,j)-result(i,j)))
			 <rel_error));
	    if (fabs(gsl_matrix_get(expected,i,j)-result(i,j))>max) {
	      max=fabs(gsl_matrix_get(expected,i,j)-result(i,j));
	    }
	  }
	}
      }
      
      process_test(ret,((std::string)"relative matrix, max=")+o2scl::dtos(max),
		   description);
      
      return ret;
      
    }

    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	\mathrm{expected}<\mathrm{rel\_error}\f$ over each element
	of an array
    */
    template<class mat_t, class mat2_t>
      bool test_rel_mat(int nr, int nc, const mat_t &result, 
			const mat2_t &expected, 
			double rel_error, std::string description) {
      bool ret=true;
      double max=0.0;
      int i, j;
  
      for(i=0;i<nr;i++) {
	for(j=0;j<nc;j++) {
	  if (o2scl::is_nan(expected(i,j))) {
	    ret=(ret && (o2scl::is_nan(expected(i,j))==o2scl::is_nan(result(i,j))));
	  } else if (o2scl::is_inf(expected(i,j))) {
	    ret=(ret && (o2scl::is_inf(expected(i,j))==o2scl::is_inf(result(i,j))));
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
      
      process_test(ret,((std::string)"relative matrix, max=")+o2scl::dtos(max),
		   description);
      
      return ret;
      
    }
    
    /** \brief Test for \f$|\mathrm{result}-\mathrm{expected}|/
	<\mathrm{abs\_error}\f$ over each element
	of an array
    */
    template<class vec_t, class vec2_t>
      bool test_abs_arr(int nv, const vec_t &result, const vec2_t &expected, 
			double rel_error, std::string description) {
      bool ret=true;
      int i;
  
      for(i=0;i<nv;i++) {
	if (o2scl::is_nan(expected[i])) {
	  ret=(ret && (o2scl::is_nan(expected[i])==o2scl::is_nan(result[i])));
	} else if (o2scl::is_inf(expected[i])) {
	  ret=(ret && (o2scl::is_inf(expected[i])==o2scl::is_inf(result[i])));
	} else {
	  ret=(ret && (fabs(expected[i]-result[i])<rel_error));
	}
      }
  
      process_test(ret,"absolute array",description);
  
      return ret;
    }

    /** \brief Test for \f$ 1/factor < result/expected < factor \f$ 
	over each element of an array
    */
    template<class vec_t, class vec2_t>
      bool test_fact_arr(int nv, const vec_t &result, const vec2_t &expected, 
			 double factor, std::string description) {
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
  
      process_test(ret,"factor array",description);
  
      return ret;
    }
    
    /// Test for equality of a generic array
    template<class vec_t>
      bool test_gen_arr(int nv, const vec_t &result, const vec_t &expected, 
			std::string description) {
      bool ret=true;
      int i;
      
      for(i=0;i<nv;i++) {
	ret=(ret && (result[i]==expected[i]));
      }
      
      process_test(ret,"generic array",description);
      
      return ret;
    }
    ///@}

    /// Add two test_mgr objects (if either failed, the sum fails)
    friend const test_mgr operator+(const test_mgr& left, 
				    const test_mgr& right);

    /// True if all tests have passed
    bool success;

    /// The description of the last failed test
    std::string last_fail;

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
