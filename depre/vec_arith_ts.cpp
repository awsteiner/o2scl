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
#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <o2scl/cx_arith.h>
#include <o2scl/vec_arith.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/ovector_cx_tlate.h>
#include <o2scl/omatrix_tlate.h>
#include <o2scl/omatrix_cx_tlate.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  // -----------------------------------------------
  // Test vector arithmetic
  // -----------------------------------------------
  
#define TSUM(t1,t2,t3)					\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    v1[0]=3.0;						\
    v1[1]=1.0;						\
    v1[2]=4.0;						\
    v1[3]=1.0;						\
    v2[0]=5.0;						\
    v2[1]=9.0;						\
    v2[2]=2.0;						\
    v2[3]=6.0;						\
    t1 v3=v1+v2;					\
    t.test_rel(v3[0],8.0,1.0e-12,"sum1");		\
    t.test_rel(v3[1],10.0,1.0e-12,"sum2");		\
    t.test_rel(v3[2],6.0,1.0e-12,"sum3");		\
    t.test_rel(v3[3],7.0,1.0e-12,"sum4");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Sum:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TSUM(ovector,ovector,ovector);
  TSUM(ovector,ovector,uvector);
  TSUM(ovector,uvector,ovector);
  TSUM(uvector,uvector,uvector);
  
#define TSUM_CX(t1,t2,t3)				\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    gsl_complex gc5={{5.0,3.0}};			\
    gsl_complex gc6={{5.0,8.0}};			\
    gsl_complex gc7={{9.0,7.0}};			\
    gsl_complex gc8={{9.0,3.0}};			\
    v1[0]=gc1;						\
    v1[1]=gc2;						\
    v1[2]=gc3;						\
    v1[3]=gc4;						\
    v2[0]=gc5;						\
    v2[1]=gc6;						\
    v2[2]=gc7;						\
    v2[3]=gc8;						\
    t1 v3=v1+v2;					\
    t.test_rel(v3[0].dat[0],8.0,1.0e-12,"sumcx1");	\
    t.test_rel(v3[0].dat[1],4.0,1.0e-12,"sumcx2");	\
    t.test_rel(v3[1].dat[0],9.0,1.0e-12,"sumcx3");	\
    t.test_rel(v3[1].dat[1],9.0,1.0e-12,"sumcx4");	\
    t.test_rel(v3[2].dat[0],14.0,1.0e-12,"sumcx5");	\
    t.test_rel(v3[2].dat[1],16.0,1.0e-12,"sumcx6");	\
    t.test_rel(v3[3].dat[0],11.0,1.0e-12,"sumcx7");	\
    t.test_rel(v3[3].dat[1],9.0,1.0e-12,"sumcx8");	\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Sum_cx: " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TSUM_CX(ovector_cx,ovector_cx,ovector_cx);
  TSUM_CX(ovector_cx,ovector_cx,uvector_cx);
  TSUM_CX(ovector_cx,uvector_cx,ovector_cx);
  TSUM_CX(uvector_cx,uvector_cx,uvector_cx);

#define TDIF(t1,t2,t3)					\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    v1[0]=3.0;						\
    v1[1]=1.0;						\
    v1[2]=4.0;						\
    v1[3]=1.0;						\
    v2[0]=5.0;						\
    v2[1]=9.0;						\
    v2[2]=2.0;						\
    v2[3]=6.0;						\
    t1 v3=v1-v2;					\
    t.test_rel(v3[0],-2.0,1.0e-12,"dif1");		\
    t.test_rel(v3[1],-8.0,1.0e-12,"dif2");		\
    t.test_rel(v3[2],2.0,1.0e-12,"dif3");		\
    t.test_rel(v3[3],-5.0,1.0e-12,"dif4");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Dif:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TDIF(ovector,ovector,ovector);
  TDIF(ovector,ovector,uvector);
  TDIF(ovector,uvector,ovector);
  TDIF(uvector,uvector,uvector);
  
#define TDIF_CX(t1,t2,t3)				\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    gsl_complex gc5={{5.0,3.0}};			\
    gsl_complex gc6={{5.0,8.0}};			\
    gsl_complex gc7={{9.0,7.0}};			\
    gsl_complex gc8={{9.0,3.0}};			\
    v1[0]=gc1;						\
    v1[1]=gc2;						\
    v1[2]=gc3;						\
    v1[3]=gc4;						\
    v2[0]=gc5;						\
    v2[1]=gc6;						\
    v2[2]=gc7;						\
    v2[3]=gc8;						\
    t1 v3=v1-v2;					\
    t.test_rel(v3[0].dat[0],-2.0,1.0e-12,"difcx1");	\
    t.test_rel(v3[0].dat[1],-2.0,1.0e-12,"difcx2");	\
    t.test_rel(v3[1].dat[0],-1.0,1.0e-12,"difcx3");	\
    t.test_rel(v3[1].dat[1],-7.0,1.0e-12,"difcx4");	\
    t.test_rel(v3[2].dat[0],-4.0,1.0e-12,"difcx5");	\
    t.test_rel(v3[2].dat[1],2.0,1.0e-12,"difcx6");	\
    t.test_rel(v3[3].dat[0],-7.0,1.0e-12,"difcx7");	\
    t.test_rel(v3[3].dat[1],3.0,1.0e-12,"difcx8");	\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Dif_cx: " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TDIF_CX(ovector_cx,ovector_cx,ovector_cx);
  TDIF_CX(ovector_cx,ovector_cx,uvector_cx);
  TDIF_CX(ovector_cx,uvector_cx,ovector_cx);
  TDIF_CX(uvector_cx,uvector_cx,uvector_cx);

#define TRMULT(t1,tm2,t3)				\
  {							\
    tm2 v1(4,4);					\
    t3 v2(4);						\
    v1[0][0]=3.0;					\
    v1[0][1]=1.0;					\
    v1[0][2]=4.0;					\
    v1[0][3]=1.0;					\
    v1[1][0]=5.0;					\
    v1[1][1]=9.0;					\
    v1[1][2]=2.0;					\
    v1[1][3]=6.0;					\
    v1[2][0]=5.0;					\
    v1[2][1]=3.0;					\
    v1[2][2]=5.0;					\
    v1[2][3]=8.0;					\
    v1[3][0]=9.0;					\
    v1[3][1]=7.0;					\
    v1[3][2]=9.0;					\
    v1[3][3]=3.0;					\
    v2[0]=2.0;						\
    v2[1]=3.0;						\
    v2[2]=8.0;						\
    v2[3]=4.0;						\
    t1 v3=v1*v2;					\
    const char *ch1=#t1;				\
    const char *ch2=#tm2;				\
    const char *ch3=#t3;				\
    t.test_rel(v3[0],45.0,1.0e-12,"rmult1");            \
    t.test_rel(v3[1],77.0,1.0e-12,"rmult2");            \
    t.test_rel(v3[2],91.0,1.0e-12,"rmult3");            \
    t.test_rel(v3[3],123.0,1.0e-12,"rmult4");           \
    cout << "Rmult:  " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TRMULT(ovector,omatrix,ovector);
  TRMULT(ovector,umatrix,ovector);
  TRMULT(ovector,omatrix,uvector);
  TRMULT(uvector,umatrix,uvector);

#define TRMULT_CX(t1,tm2,t3)				\
  {							\
    tm2 v1(4,4);					\
    t3 v2(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    v1[0][0]=gc1;					\
    v1[0][1]=gc2;					\
    v1[0][2]=gc3;					\
    v1[0][3]=gc4;					\
    v1[1][0]=gc2;					\
    v1[1][1]=gc3;					\
    v1[1][2]=gc4;					\
    v1[1][3]=gc1;					\
    v1[2][0]=gc3;					\
    v1[2][1]=gc4;					\
    v1[2][2]=gc1;					\
    v1[2][3]=gc2;					\
    v1[3][0]=gc4;					\
    v1[3][1]=gc1;					\
    v1[3][2]=gc2;					\
    v1[3][3]=gc3;					\
    v2[0]=gc1;						\
    v2[1]=gc2;						\
    v2[2]=gc3;						\
    v2[3]=gc4;						\
    t1 v3=v1*v2;					\
    const char *ch1=#t1;				\
    const char *ch2=#tm2;				\
    const char *ch3=#t3;				\
    t.test_rel(v3[0].dat[0],-65.0,1.0e-12,"rmultcx1");  \
    t.test_rel(v3[0].dat[1],128.0,1.0e-12,"rmultcx2");  \
    t.test_rel(v3[1].dat[0],-22.0,1.0e-12,"rmultcx3");  \
    t.test_rel(v3[1].dat[1],116.0,1.0e-12,"rmultcx4");  \
    t.test_rel(v3[2].dat[0],16.0,1.0e-12,"rmultcx5");	\
    t.test_rel(v3[2].dat[1],116.0,1.0e-12,"rmultcx6");	\
    t.test_rel(v3[3].dat[0],-22.0,1.0e-12,"rmultcx7");  \
    t.test_rel(v3[3].dat[1],116.0,1.0e-12,"rmultcx8");  \
    cout << "Rmultcx:  " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TRMULT_CX(ovector_cx,omatrix_cx,ovector_cx);
  TRMULT_CX(ovector_cx,omatrix_cx,uvector_cx);
  TRMULT_CX(ovector_cx,umatrix_cx,ovector_cx);
  TRMULT_CX(uvector_cx,umatrix_cx,uvector_cx);
  
#define TLMULT(t1,tm2,t3)				\
  {							\
    tm2 v1(4,4);					\
    t3 v2(4);						\
    v1[0][0]=3.0;					\
    v1[0][1]=1.0;					\
    v1[0][2]=4.0;					\
    v1[0][3]=1.0;					\
    v1[1][0]=5.0;					\
    v1[1][1]=9.0;					\
    v1[1][2]=2.0;					\
    v1[1][3]=6.0;					\
    v1[2][0]=5.0;					\
    v1[2][1]=3.0;					\
    v1[2][2]=5.0;					\
    v1[2][3]=8.0;					\
    v1[3][0]=9.0;					\
    v1[3][1]=7.0;					\
    v1[3][2]=9.0;					\
    v1[3][3]=3.0;					\
    v2[0]=2.0;						\
    v2[1]=3.0;						\
    v2[2]=8.0;						\
    v2[3]=4.0;						\
    t1 v3=v2*v1;					\
    const char *ch1=#t1;				\
    const char *ch2=#tm2;				\
    const char *ch3=#t3;				\
    t.test_rel(v3[0],97.0,1.0e-12,"lmult1");            \
    t.test_rel(v3[1],81.0,1.0e-12,"lmult2");            \
    t.test_rel(v3[2],90.0,1.0e-12,"lmult3");            \
    t.test_rel(v3[3],96.0,1.0e-12,"lmult4");		\
    cout << "Lmult:  " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TLMULT(ovector,omatrix,ovector);
  TLMULT(ovector,umatrix,ovector);
  TLMULT(ovector,omatrix,uvector);
  TLMULT(uvector,umatrix,uvector);
  
#define TLMULT2(t1,tm2,t3)				\
  {							\
    tm2 v1(4,4);					\
    t3 v2(4);						\
    v1[0][0]=3.0;					\
    v1[0][1]=1.0;					\
    v1[0][2]=4.0;					\
    v1[0][3]=1.0;					\
    v1[1][0]=5.0;					\
    v1[1][1]=9.0;					\
    v1[1][2]=2.0;					\
    v1[1][3]=6.0;					\
    v1[2][0]=5.0;					\
    v1[2][1]=3.0;					\
    v1[2][2]=5.0;					\
    v1[2][3]=8.0;					\
    v1[3][0]=9.0;					\
    v1[3][1]=7.0;					\
    v1[3][2]=9.0;					\
    v1[3][3]=3.0;					\
    v2[0]=2.0;						\
    v2[1]=3.0;						\
    v2[2]=8.0;						\
    v2[3]=4.0;						\
    t1 v3=trans_mult(v2,v1);				\
    const char *ch1=#t1;				\
    const char *ch2=#tm2;				\
    const char *ch3=#t3;				\
    t.test_rel(v3[0],97.0,1.0e-12,"lmult2.1");          \
    t.test_rel(v3[1],81.0,1.0e-12,"lmult2.2");          \
    t.test_rel(v3[2],90.0,1.0e-12,"lmult2.3");          \
    t.test_rel(v3[3],96.0,1.0e-12,"lmult2.4");		\
    cout << "Lmult2:  " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TLMULT2(ovector,omatrix,ovector);
  TLMULT2(ovector,umatrix,ovector);
  TLMULT2(ovector,omatrix,uvector);
  TLMULT2(uvector,umatrix,uvector);
  
#define TDOT(t1,t2,t3)					\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    v1[0]=3.0;						\
    v1[1]=1.0;						\
    v1[2]=4.0;						\
    v1[3]=1.0;						\
    v2[0]=5.0;						\
    v2[1]=9.0;						\
    v2[2]=2.0;						\
    v2[3]=6.0;						\
    t1 v3=dot(v1,v2);					\
    t.test_rel(v3,38,1.0e-12,"dot1");			\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Dot:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TDOT(double,ovector,ovector);
  TDOT(double,ovector,uvector);
  TDOT(double,uvector,ovector);
  TDOT(double,uvector,uvector);

#define TDOTC(t1,t2,t3)					\
  {							\
    t2 v1(4);						\
    t3 v2(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    gsl_complex gc5={{5.0,3.0}};			\
    gsl_complex gc6={{5.0,8.0}};			\
    gsl_complex gc7={{9.0,7.0}};			\
    gsl_complex gc8={{9.0,3.0}};			\
    v1[0]=gc1;						\
    v1[1]=gc2;						\
    v1[2]=gc3;						\
    v1[3]=gc4;						\
    v2[0]=gc5;						\
    v2[1]=gc6;						\
    v2[2]=gc7;						\
    v2[3]=gc8;						\
    t1 v3=dot(v1,v2);					\
    t.test_rel(v3.dat[0],6,1.0e-12,"dot1");		\
    t.test_rel(v3.dat[1],227,1.0e-12,"dot1");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Dot cx:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }
  
  TDOTC(gsl_complex,ovector_cx,ovector_cx);
  TDOTC(gsl_complex,ovector_cx,uvector_cx);
  TDOTC(gsl_complex,uvector_cx,ovector_cx);
  TDOTC(gsl_complex,uvector_cx,uvector_cx);

#define TLMULTSV(t1,t2,t3)				\
  {							\
    t1 v1=3.0;						\
    t2 v2(4);						\
    v2[0]=1.0;						\
    v2[1]=4.0;						\
    v2[2]=1.0;						\
    v2[3]=5.0;						\
    t3 v3=v1*v2;					\
    t.test_rel(v3[0],3.0,1.0e-12,"lmultsv1");		\
    t.test_rel(v3[1],12.0,1.0e-12,"lmultsv2");		\
    t.test_rel(v3[2],3.0,1.0e-12,"lmultsv3");		\
    t.test_rel(v3[3],15.0,1.0e-12,"lmultsv4");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Lmultsv:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TLMULTSV(double,ovector,ovector);
  TLMULTSV(double,uvector,uvector);

#define TRMULTSV(t1,t2,t3)				\
  {							\
    t1 v1=3.0;						\
    t2 v2(4);						\
    v2[0]=1.0;						\
    v2[1]=4.0;						\
    v2[2]=1.0;						\
    v2[3]=5.0;						\
    t3 v3=v2*v1;					\
    t.test_rel(v3[0],3.0,1.0e-12,"rmultsv1");		\
    t.test_rel(v3[1],12.0,1.0e-12,"rmultsv2");		\
    t.test_rel(v3[2],3.0,1.0e-12,"rmultsv3");		\
    t.test_rel(v3[3],15.0,1.0e-12,"rmultsv4");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Rmultsv:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TRMULTSV(double,ovector,ovector);
  TRMULTSV(double,uvector,uvector);

#define TPPRR(t1,t2,t3)					\
  {							\
    t2 v2(4);						\
    t3 v3(4);						\
    v2[0]=3.0;						\
    v2[1]=1.0;						\
    v2[2]=4.0;						\
    v2[3]=1.0;						\
    v3[0]=5.0;						\
    v3[1]=9.0;						\
    v3[2]=2.0;						\
    v3[3]=6.0;						\
    t1 v1=pair_prod(v2,v3);				\
    t.test_rel(v1[0],15.0,1.0e-12,"tpprr1");		\
    t.test_rel(v1[1],9.0,1.0e-12,"tpprr2");		\
    t.test_rel(v1[2],8.0,1.0e-12,"tpprr3");		\
    t.test_rel(v1[3],6.0,1.0e-12,"tpprr4");		\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Tpprr:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TPPRR(ovector,ovector,ovector);
  TPPRR(ovector,uvector,ovector);
  TPPRR(ovector,ovector,uvector);
  TPPRR(uvector,uvector,uvector);

#define TPPRC(t1,t2,t3)					\
  {							\
    t2 v2(4);						\
    t3 v3(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    v2[0]=5.0;						\
    v2[1]=4.0;						\
    v2[2]=5.0;						\
    v2[3]=6.0;						\
    v3[0]=gc1;						\
    v3[1]=gc2;						\
    v3[2]=gc3;						\
    v3[3]=gc4;						\
    t1 v1=pair_prod(v2,v3);				\
    t.test_rel(v1[0].dat[0],15.0,1.0e-12,"tpprc1");	\
    t.test_rel(v1[0].dat[1],5.0,1.0e-12,"tpprc2");	\
    t.test_rel(v1[1].dat[0],16.0,1.0e-12,"tpprc3");	\
    t.test_rel(v1[1].dat[1],4.0,1.0e-12,"tpprc4");	\
    t.test_rel(v1[2].dat[0],25.0,1.0e-12,"tpprc5");	\
    t.test_rel(v1[2].dat[1],45.0,1.0e-12,"tpprc6");	\
    t.test_rel(v1[3].dat[0],12.0,1.0e-12,"tpprc7");	\
    t.test_rel(v1[3].dat[1],36.0,1.0e-12,"tpprc8");	\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Tpprc:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TPPRC(ovector_cx,ovector,ovector_cx);
  TPPRC(ovector_cx,uvector,ovector_cx);
  TPPRC(ovector_cx,ovector,uvector_cx);
  TPPRC(uvector_cx,uvector,uvector_cx);

#define TPPCR(t1,t2,t3)					\
  {							\
    t2 v2(4);						\
    t3 v3(4);						\
    gsl_complex gc1={{3.0,1.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    v3[0]=5.0;						\
    v3[1]=4.0;						\
    v3[2]=5.0;						\
    v3[3]=6.0;						\
    v2[0]=gc1;						\
    v2[1]=gc2;						\
    v2[2]=gc3;						\
    v2[3]=gc4;						\
    t1 v1=pair_prod(v2,v3);				\
    t.test_rel(v1[0].dat[0],15.0,1.0e-12,"tppcr1");	\
    t.test_rel(v1[0].dat[1],5.0,1.0e-12,"tppcr2");	\
    t.test_rel(v1[1].dat[0],16.0,1.0e-12,"tppcr3");	\
    t.test_rel(v1[1].dat[1],4.0,1.0e-12,"tppcr4");	\
    t.test_rel(v1[2].dat[0],25.0,1.0e-12,"tppcr5");	\
    t.test_rel(v1[2].dat[1],45.0,1.0e-12,"tppcr6");	\
    t.test_rel(v1[3].dat[0],12.0,1.0e-12,"tppcr7");	\
    t.test_rel(v1[3].dat[1],36.0,1.0e-12,"tppcr8");	\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Tppcr:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TPPCR(ovector_cx,ovector_cx,ovector);
  TPPCR(ovector_cx,uvector_cx,ovector);
  TPPCR(ovector_cx,ovector_cx,uvector);
  TPPCR(uvector_cx,uvector_cx,uvector);

#define TPPCC(t1,t2,t3)					\
  {							\
    t2 v2(4);						\
    t3 v3(4);						\
    gsl_complex gc1={{3.0,2.0}};			\
    gsl_complex gc2={{4.0,1.0}};			\
    gsl_complex gc3={{5.0,9.0}};			\
    gsl_complex gc4={{2.0,6.0}};			\
    v3[0]=gc4;						\
    v3[1]=gc3;						\
    v3[2]=gc2;						\
    v3[3]=gc1;						\
    v2[0]=gc1;						\
    v2[1]=gc2;						\
    v2[2]=gc3;						\
    v2[3]=gc4;						\
    t1 v1=pair_prod(v2,v3);				\
    t.test_rel(v1[0].dat[0],-6.0,1.0e-12,"tppcc1");	\
    t.test_rel(v1[0].dat[1],22.0,1.0e-12,"tppcc2");	\
    t.test_rel(v1[1].dat[0],11.0,1.0e-12,"tppcc3");	\
    t.test_rel(v1[1].dat[1],41.0,1.0e-12,"tppcc4");	\
    t.test_rel(v1[2].dat[0],11.0,1.0e-12,"tppcc5");	\
    t.test_rel(v1[2].dat[1],41.0,1.0e-12,"tppcc6");	\
    t.test_rel(v1[3].dat[0],-6.0,1.0e-12,"tppcc7");	\
    t.test_rel(v1[3].dat[1],22.0,1.0e-12,"tppcc8");	\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    const char *ch3=#t3;				\
    cout << "Tppcc:    " << ch1 << " " << ch2 << " "	\
	 << ch3 << endl;				\
  }

  TPPCC(ovector_cx,ovector_cx,ovector_cx);
  TPPCC(ovector_cx,uvector_cx,ovector_cx);
  TPPCC(ovector_cx,ovector_cx,uvector_cx);
  TPPCC(uvector_cx,uvector_cx,uvector_cx);

#define TEQRR(t1,t2)					\
  {							\
    t1 v1(4);						\
    t2 v2(4);						\
    t2 v3(4);						\
    v1[0]=3.0;						\
    v1[1]=1.0;						\
    v1[2]=4.0;						\
    v1[3]=1.0;						\
    v2[0]=3.0;						\
    v2[1]=1.0;						\
    v2[2]=4.0;						\
    v2[3]=1.0;						\
    v3[0]=5.0;						\
    v3[1]=9.0;						\
    v3[2]=2.0;						\
    v3[3]=6.0;						\
    bool b1=(v1==v2);					\
    bool b2=(v1==v3);					\
    t.test_gen(b1==true,"teqrr3");			\
    t.test_gen(b2==false,"teqrr4");			\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    cout << "Teqrr:    " << ch1 << " " << ch2 << endl;	\
  }

  TEQRR(ovector,ovector);
  TEQRR(uvector,ovector);
  TEQRR(ovector,uvector);
  TEQRR(uvector,uvector);

#define TNEQRR(t1,t2)					\
  {							\
    t1 v1(4);						\
    t2 v2(4);						\
    t2 v3(4);						\
    v1[0]=3.0;						\
    v1[1]=1.0;						\
    v1[2]=4.0;						\
    v1[3]=1.0;						\
    v2[0]=3.0;						\
    v2[1]=1.0;						\
    v2[2]=4.0;						\
    v2[3]=1.0;						\
    v3[0]=5.0;						\
    v3[1]=9.0;						\
    v3[2]=2.0;						\
    v3[3]=6.0;						\
    bool b1=(v1!=v2);					\
    bool b2=(v1!=v3);					\
    t.test_gen(b1==false,"tneqrr3");			\
    t.test_gen(b2==true,"tneqrr4");			\
    const char *ch1=#t1;				\
    const char *ch2=#t2;				\
    cout << "Tneqrr:    " << ch1 << " " << ch2 << endl;	\
  }

  TNEQRR(ovector,ovector);
  TNEQRR(uvector,ovector);
  TNEQRR(ovector,uvector);
  TNEQRR(uvector,uvector);

  t.report();
  
  return 0;
}

