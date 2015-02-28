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
#include <o2scl/pinside.h>

using namespace std;
using namespace o2scl;

int pinside::intersect(line P, line Q) { 
  
  double a, b, c, d; 
  
  a= (Q.p1.y - P.p1.y)*(P.p2.x - P.p1.x) 
    -(P.p2.y - P.p1.y)*(Q.p1.x - P.p1.x);
  b= (Q.p2.y - P.p1.y)*(P.p2.x - P.p1.x) 
    -(P.p2.y - P.p1.y)*(Q.p2.x - P.p1.x);
  c= (P.p1.y - Q.p1.y)*(Q.p2.x - Q.p1.x) 
    -(Q.p2.y - Q.p1.y)*(P.p1.x - Q.p1.x);
  d= (P.p2.y - Q.p1.y)*(Q.p2.x - Q.p1.x) 
    -(Q.p2.y - Q.p1.y)*(P.p2.x - Q.p1.x);


  if (a == 0 && b == 0) { 
    if ((((P.p1.x<=Q.p1.x && Q.p1.x<=P.p2.x) 
	  || (P.p2.x<=Q.p1.x && Q.p1.x<=P.p1.x)) && 
	 ((P.p1.y<=Q.p1.y && Q.p1.y<=P.p2.y) 
	  || (P.p2.y<=Q.p1.y && Q.p1.y<=P.p1.y))) || 
	(((P.p1.x<=Q.p2.x && Q.p2.x<=P.p2.x) 
	  || (P.p2.x<=Q.p2.x && Q.p2.x<=P.p1.x)) && 

	 ((P.p1.y<=Q.p2.y && Q.p2.y<=P.p2.y) 
	  || (P.p2.y<=Q.p2.y && Q.p2.y<=P.p1.y)))) {
      return 1; 
    }
      
    if ((((Q.p1.x<=P.p1.x && P.p1.x<=Q.p2.x) 
	  || (Q.p2.x<=P.p1.x && P.p1.x<=Q.p1.x)) && 
	 ((Q.p1.y<=P.p1.y && P.p1.y<=Q.p2.y) 
	  || (Q.p2.y<=P.p1.y && P.p1.y<=Q.p1.y))) || 
	(((Q.p1.x<=P.p2.x && P.p2.x<=Q.p2.x) 
	  || (Q.p2.x<=P.p2.x && P.p2.x<=Q.p1.x)) && 

	 ((Q.p1.y<=P.p2.y && P.p2.y<=Q.p2.y) 
	  || (Q.p2.y<=P.p2.y && P.p2.y<=Q.p1.y)))) {
      return 1; 
    }

    return success;
  }

  if (a*b <= 0 && c*d <= 0) {
    return 1;
  }

  return success;
} 

int pinside::inside(point t, point p[], int N) { 
  int i, count = 0, j = 0;
  double max = p[1].x; 
  line lt, lp; 

  p[0] = p[N];
  lt.p1 = t; lt.p2 = t; 

  for (i = 2; i <= N; i++) {
    if (max < p[i].x) {
      max = p[i].x; 
    }
  }

  lt.p2.x = max + 1; 

  for (i = 1; i <= N; i++) {

    lp.p1 = p[i]; lp.p2 = p[i];

    if (!intersect(lp, lt)) {
      
      lp.p2 = p[j]; 
      
      if (intersect(lp, lt) || 
	  (j < i-1 && (lp.p1.y-t.y)*(lp.p2.y-t.y) < 0)) {
	count++; 
      }
      
      j = i;
    }
    
  }
  
  if (j < i-1 && (p[j].y - t.y) * (p[1].y - t.y) > 0) {
    // This was printed as "count-" in Lewis' original article.
    // It has been corrected here to "count--"
    count--;
  }
  
  return (count%2!=0); 
} 

int pinside::inside(double x, double y, const ubvector &xa,
		    const ubvector &ya) {

  int N=xa.size(), ix;
  point t, *p=new point[N+1];
  t.x=x;
  t.y=y;

  // We have to copy the vectors so we can rearrange them because
  // they are const
  ubvector xb(N), yb(N);
  vector_copy(N,xa,xb);
  vector_copy(N,ya,yb);

  // Ensure that (yb[0],ya[0]) is the point with the smallest x
  // coordinate among all the points with the smallest y coordinate
  double xmin=xb[0];
  double ymin=yb[0];
  ix=0;
  for(int i=0;i<N;i++) {
    if (yb[i]<ymin) {
      ymin=yb[i];
      ix=i;
    }
  }
  for(int i=0;i<N;i++) {
    if (yb[i]==ymin && xb[i]<xmin) {
      xmin=xb[i];
      ix=i;
    }
  }
  vector_rotate<ubvector,double>(N,xb,ix);
  vector_rotate<ubvector,double>(N,yb,ix);

  // Copy to p[]
  for(int i=0;i<N;i++) {
    p[i+1].x=xb[i];
    p[i+1].y=yb[i];
  }
    
  int ret=inside(t,p,N);
  delete[] p;

  return ret;
}

int pinside::test(test_mgr &t) {
  line tx1={{0,0},{2,2}};
  line ty1={{0,1},{1,0}};
  t.test_gen(intersect(tx1,ty1)==1,"i1");
  line tx2={{0,0},{0.49,0.49}};
  line ty2={{0,1},{1,0}};
  t.test_gen(intersect(tx2,ty2)==0,"i2");
  line tx3={{0,0},{-1,-1}};
  line ty3={{0,1},{1,0}};
  t.test_gen(intersect(tx3,ty3)==0,"i3");

  ubvector x(12), y(12);
  x[0]=1; y[0]=0;
  x[1]=2; y[1]=0;
  x[2]=2; y[2]=1;
  x[3]=3; y[3]=1;
  x[4]=3; y[4]=2;
  x[5]=2; y[5]=2;
  x[6]=2; y[6]=3;
  x[7]=1; y[7]=3;
  x[8]=1; y[8]=2;
  x[9]=0; y[9]=2;
  x[10]=0; y[10]=1;
  x[11]=1; y[11]=1;
  
  // Test with ubvectors

  t.test_gen(inside(1.5,1.5,x,y)==1,"in1");
  t.test_gen(inside(2.01,1.99,x,y)==1,"in2");
  t.test_gen(inside(1.99,2.01,x,y)==1,"in3");
  t.test_gen(inside(1.99,1.99,x,y)==1,"in4");

  t.test_gen(inside(2.01,2.01,x,y)==0,"out1");
  t.test_gen(inside(0.5,0.5,x,y)==0,"out2");
  t.test_gen(inside(2.5,2.5,x,y)==0,"out3");

  t.test_gen(inside(0.5,1.5,x,y)==1,"c1");
  t.test_gen(inside(0.5,2.5,x,y)==0,"c2");
  t.test_gen(inside(1.5,0.5,x,y)==1,"c3");
  t.test_gen(inside(1.5,2.5,x,y)==1,"c4");
  t.test_gen(inside(2.5,0.5,x,y)==0,"c5");
  t.test_gen(inside(2.5,1.5,x,y)==1,"c6");

  // Test template version

  t.test_gen(inside(1.5,1.5,x.size(),x,y)==1,"in1");
  t.test_gen(inside(2.01,1.99,x.size(),x,y)==1,"in2");
  t.test_gen(inside(1.99,2.01,x.size(),x,y)==1,"in3");
  t.test_gen(inside(1.99,1.99,x.size(),x,y)==1,"in4");

  t.test_gen(inside(2.01,2.01,x.size(),x,y)==0,"out1");
  t.test_gen(inside(0.5,0.5,x.size(),x,y)==0,"out2");
  t.test_gen(inside(2.5,2.5,x.size(),x,y)==0,"out3");

  t.test_gen(inside(0.5,1.5,x.size(),x,y)==1,"c1");
  t.test_gen(inside(0.5,2.5,x.size(),x,y)==0,"c2");
  t.test_gen(inside(1.5,0.5,x.size(),x,y)==1,"c3");
  t.test_gen(inside(1.5,2.5,x.size(),x,y)==1,"c4");
  t.test_gen(inside(2.5,0.5,x.size(),x,y)==0,"c5");
  t.test_gen(inside(2.5,1.5,x.size(),x,y)==1,"c6");

  return success;
}
