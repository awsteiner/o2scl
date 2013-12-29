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

#include <o2scl/vec_arith.h>
#include <o2scl/cx_arith.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_NEVER_DEFINED
{
#endif

  O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector,o2scl::ovector_const_view,
			  o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector,o2scl::ovector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector,o2scl::uvector_const_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::uvector,o2scl::uvector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::uvector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_ADD(o2scl::uvector_cx,o2scl::uvector_cx_view,
			    o2scl::uvector_cx_view)
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector,o2scl::ovector_const_view,
			  o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector,o2scl::ovector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector,o2scl::uvector_const_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::uvector,o2scl::uvector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::uvector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_SUB(o2scl::uvector_cx,o2scl::uvector_cx_view,
			    o2scl::uvector_cx_view)
      
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_MAT_VEC_MULT(o2scl::ovector,o2scl::ovector_const_view,
			   o2scl::omatrix_const_view)
    O2SCL_OPSRC_MAT_VEC_MULT(o2scl::ovector,o2scl::ovector_const_view,
			     o2scl::umatrix_const_view)
    O2SCL_OPSRC_MAT_VEC_MULT(o2scl::ovector,o2scl::uvector_const_view,
			     o2scl::omatrix_const_view)
    O2SCL_OPSRC_MAT_VEC_MULT(o2scl::uvector,o2scl::uvector_const_view,
			     o2scl::umatrix_const_view)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
			     o2scl::omatrix_cx_view)
    O2SCL_OPSRC_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
			       o2scl::umatrix_cx_view)
    O2SCL_OPSRC_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::uvector_cx_view,
			       o2scl::omatrix_cx_view)
    O2SCL_OPSRC_CMAT_CVEC_MULT(o2scl::uvector_cx,o2scl::uvector_cx_view,
			       o2scl::umatrix_cx_view)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector,o2scl::ovector_const_view,
			   o2scl::omatrix_const_view)
    O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector,o2scl::ovector_const_view,
			     o2scl::umatrix_const_view)
    O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector,o2scl::uvector_const_view,
			     o2scl::omatrix_const_view)
    O2SCL_OPSRC_VEC_MAT_MULT(o2scl::uvector,o2scl::uvector_const_view,
			     o2scl::umatrix_const_view)
    //O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
    //o2scl::omatrix_cx_view)
    //O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
    //o2scl::umatrix_cx_view)
    //O2SCL_OPSRC_VEC_MAT_MULT(o2scl::ovector_cx,o2scl::uvector_cx_view,
    //o2scl::omatrix_cx_view)
    //O2SCL_OPSRC_VEC_MAT_MULT(o2scl::uvector_cx,o2scl::uvector_cx_view,
    //o2scl::umatrix_cx_view)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_TRANS_MULT(o2scl::ovector,o2scl::ovector_const_view,
			 o2scl::omatrix_const_view)
    O2SCL_OPSRC_TRANS_MULT(o2scl::ovector,o2scl::ovector_const_view,
			   o2scl::umatrix_const_view)
    O2SCL_OPSRC_TRANS_MULT(o2scl::ovector,o2scl::uvector_const_view,
			   o2scl::omatrix_const_view)
    O2SCL_OPSRC_TRANS_MULT(o2scl::uvector,o2scl::uvector_const_view,
			   o2scl::umatrix_const_view)
    //O2SCL_OPSRC_TRANS_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
    //o2scl::omatrix_cx_view)
    //O2SCL_OPSRC_TRANS_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
    //o2scl::umatrix_cx_view)
    //O2SCL_OPSRC_TRANS_MULT(o2scl::ovector_cx,o2scl::uvector_cx_view,
    //o2scl::omatrix_cx_view)
    //O2SCL_OPSRC_TRANS_MULT(o2scl::uvector_cx,o2scl::uvector_cx_view,
    //o2scl::umatrix_cx_view)
	    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_DOT_PROD(double,o2scl::ovector_const_view,
		       o2scl::ovector_const_view)
    O2SCL_OPSRC_DOT_PROD(double,o2scl::ovector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OPSRC_DOT_PROD(double,o2scl::uvector_const_view,
			 o2scl::ovector_const_view)
    O2SCL_OPSRC_DOT_PROD(double,o2scl::uvector_const_view,
			 o2scl::uvector_const_view)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_CX_DOT_PROD(gsl_complex,o2scl::ovector_cx_view,
			  o2scl::ovector_cx_view)
    O2SCL_OPSRC_CX_DOT_PROD(gsl_complex,o2scl::ovector_cx_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_CX_DOT_PROD(gsl_complex,o2scl::uvector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_CX_DOT_PROD(gsl_complex,o2scl::uvector_cx_view,
			    o2scl::uvector_cx_view)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_SCA_VEC_MULT(double,o2scl::ovector_const_view,o2scl::ovector)
    O2SCL_OPSRC_SCA_VEC_MULT(double,o2scl::uvector_const_view,o2scl::uvector)
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  O2SCL_OPSRC_VEC_SCA_MULT(double,o2scl::ovector_const_view,o2scl::ovector)
    O2SCL_OPSRC_VEC_SCA_MULT(double,o2scl::uvector_const_view,o2scl::uvector)

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
  
  O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector,o2scl::ovector_const_view,
			  o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector,o2scl::ovector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector,o2scl::uvector_const_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::uvector,o2scl::uvector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_cx_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_cx_view,
			    o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_const_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_const_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_const_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_const_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::uvector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OPSRC_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_cx_view,
			    o2scl::uvector_cx_view)
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
  
  O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			    o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			      o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			      o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			      o2scl::uvector_const_view)
    /*
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
      o2scl::ovector_const_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
      o2scl::uvector_const_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
      o2scl::ovector_const_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
      o2scl::uvector_const_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_const_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_const_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_const_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_const_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
      o2scl::uvector_cx_view)
    */
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
  
  O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			     o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			       o2scl::uvector_const_view)
    O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			       o2scl::ovector_const_view)
    O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			       o2scl::uvector_const_view)
    /*
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
      o2scl::ovector_const_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
      o2scl::uvector_const_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
      o2scl::ovector_const_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
      o2scl::uvector_const_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
      o2scl::uvector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
      o2scl::ovector_cx_view)
      O2SCL_OPSRC_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
      o2scl::uvector_cx_view)
    */
    
#ifdef O2SCL_NEVER_DEFINED
    }
#endif

