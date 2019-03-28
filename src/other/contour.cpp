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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/contour.h>

using namespace std;
using namespace o2scl;

contour::contour() {
  nx=0;
  lev_adjust=1.0e-8;
  verbose=0;
  debug_next_point=false;
}

contour::~contour() {
}

void contour::check_data() {
  bool increasing=false;
  if (xfun[0]<xfun[1]) increasing=true;
  for(int j=1;j<nx;j++) {
    if (increasing && xfun[j-1]>=xfun[j]) {
      O2SCL_ERR("The 'x' array is not monotonic",exc_einval);
    } else if (!increasing && xfun[j-1]<=xfun[j]) {
      O2SCL_ERR("The 'x' array is not monotonic",exc_einval);
    }
  }
  increasing=false;
  if (yfun[0]<yfun[1]) increasing=true;
  for(int k=1;k<ny;k++) {
    if (increasing && yfun[k-1]>=yfun[k]) {
      O2SCL_ERR("The 'y' array is not monotonic",exc_einval);
    } else if (!increasing && yfun[k-1]<=yfun[k]) {
      O2SCL_ERR("The 'y' array is not monotonic",exc_einval);
    }
  }
  return;
}

void contour::regrid_data(size_t xfact, size_t yfact, size_t interp_type) {
  
  if (nx==0) {
    O2SCL_ERR("Data not set in regrid_data().",exc_einval);
  }
  if ((xfact<=1 && yfact<=1) || xfact==0 || yfact==0) {
    O2SCL_ERR("Parameters xfact and/or yfact too small in regrid_data().",
	      exc_einval);
  }

  int newx=(nx-1)*xfact+1, newy=(ny-1)*yfact+1;

  ubvector xfun_old=xfun, yfun_old=yfun;
  ubmatrix data_old=data;

  interp_vec<ubvector> si;

  // Create new xfun
  ubvector xidx(nx);
  xfun.resize(newx);
  for(int i=0;i<nx;i++) xidx[i]=((double)i);
  si.set(nx,xidx,xfun_old,interp_type);
  for(int i=0;i<newx;i++) {
    xfun[i]=si.eval(((double)i)/(newx-1+xfact)*nx);
  }

  // Create new yfun
  ubvector yidx(ny);
  yfun.resize(newy);
  for(int i=0;i<ny;i++) yidx[i]=((double)i);
  si.set(ny,yidx,yfun_old,interp_type);
  for(int i=0;i<newy;i++) {
    yfun[i]=si.eval(((double)i)/(newy-1+yfact)*ny);
  }
  
  // Perform the first round, expanding the number of rows from 
  // nx -> newx 
  ubmatrix bc(newx,ny);
  for(int i=0;i<newy;i+=yfact) {
    // Get the column of data_old with index i/yfact
    ubvector at(nx);
    for(int ij=0;ij<nx;ij++) at[ij]=data_old(ij,i/yfact);
    // Use it to interpolate a column into bc
    si.set(nx,xfun_old,at,interp_type);
    for(int j=0;j<newx;j++) {
      bc(j,i/yfact)=si.eval(xfun[j]);
    }
  }
  
  // Perform the second round, expanding the number of columns
  // from ny -> newy
  data.resize(newx,newy);
  for(int i=0;i<newx;i++) {
    ubvector at2(ny);
    for(int ij=0;ij<ny;ij++) at2[ij]=bc(i,ij);
    si.set(ny,yfun_old,at2,interp_type);
    for(int j=0;j<newy;j++) {
      data(i,j)=si.eval(yfun[j]);
    }
  }

  // Modify size
  nx=newx;
  ny=newy;
  
  return;
}

int contour::find_next_point_y_direct(int j, int k, int &jnext, int &knext, 
				      int &dir_next, int nsw,
				      edge_crossings &xedges,
				      edge_crossings &yedges) {
  
  double closest=0.0, next;
  bool found=false;

  if (true) {

    size_t xlo_count=0, xhi_count=0;
    if (j>=1) {
      if (yedges.status(j-1,k)!=empty) {
	xlo_count++;
      }
      if (xedges.status(j-1,k)!=empty) {
	xlo_count++;
      }
      if (xedges.status(j-1,k+1)!=empty) {
	xlo_count++;
      }
    }
    if (j<nx-1) {
      if (xedges.status(j,k)!=empty) {
	xhi_count++;
      }
      if (xedges.status(j,k+1)!=empty) {
	xhi_count++;
      }
      if (yedges.status(j+1,k)!=empty) {
	xhi_count++;
      }
    }

    if (j>0 && j<nx-1 && xhi_count+xlo_count%2==1) {
      O2SCL_ERR("Malformed edge in find_next_point_y_direct().",exc_esanity);
    }

  }
  
  if (j<nx-1) {
    
    double xscale=fabs(xfun[j+1]-xfun[j]);
    double yscale=fabs(yfun[k+1]-yfun[k]);
    
    // down and left
    if (xedges.status(j,k)==nsw) {
      closest=sqrt(pow(yedges.values(j,k)-yfun[k],2.0)/yscale/yscale+
		   pow(xfun[j]-xedges.values(j,k),2.0)/xscale/xscale);
      found=true;
      jnext=j;
      knext=k;
      dir_next=dxdir;
    }
    
    // down
    if (yedges.status(j+1,k)==nsw) {
      next=sqrt(pow(yedges.values(j,k)-yedges.values(j+1,k),2.0)/yscale/
		yscale+pow(xfun[j]-xfun[j+1],2.0)/xscale/xscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j+1;
	knext=k;
	dir_next=dydir;
	closest=next;
      }
    }
    
    // down and right
    if (xedges.status(j,k+1)==nsw) {
      next=sqrt(pow(yedges.values(j,k)-yfun[k+1],2.0)/yscale/yscale+
		pow(xfun[j]-xedges.values(j,k+1),2.0)/xscale/xscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j;
	knext=k+1;
	dir_next=dxdir;
	closest=next;
      }
    }
    
  }

  if (j>0) {
  
    double xscale=fabs(xfun[j]-xfun[j-1]);
    double yscale=fabs(yfun[k+1]-yfun[k]);
    
    // up
    if (yedges.status(j-1,k)==nsw) { 
      next=sqrt(pow(yedges.values(j,k)-yedges.values(j-1,k),2.0)/yscale/
		yscale+pow(xfun[j]-xfun[j-1],2.0)/xscale/xscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j-1;
	knext=k;
	dir_next=dydir;
	closest=next;
      }
    }
    
    // up and left
    if (xedges.status(j-1,k)==nsw) {
      next=sqrt(pow(yedges.values(j,k)-yfun[k],2.0)/yscale/yscale+
		pow(xfun[j]-xedges.values(j-1,k),2.0)/xscale/xscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j-1;
	knext=k;
	dir_next=dxdir;
	closest=next;
      }
    }
    
    // up and right
    if (xedges.status(j-1,k+1)==nsw) {
      next=sqrt(pow(yedges.values(j,k)-yfun[k+1],2.0)/yscale/yscale+
		pow(xfun[j]-xedges.values(j-1,k+1),2.0)/xscale/xscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j-1;
	knext=k+1;
	dir_next=dxdir;
	closest=next;
      }
    }
  }

  if (debug_next_point) {

    cout.precision(4);
    cout << "y edges: " << endl;
    cout << "j,k,nx,ny: " << j << " " << k << " " << nx << " " << ny << endl;
    if (j>=1) {

      // Top points and edge above
      cout << "(" << xfun[j-1] << "," << yfun[k] << ")" << " ";
      if (yedges.status(j-1,k)!=empty) {
	cout << "[" << xfun[j-1] << "," 
	     << yedges.values(j-1,k) << "]" << " ";
      } else {
	cout << "[" << xfun[j-1] << "," << "          " << "]" << " ";
      }
      cout << "(" << xfun[j-1] << "," << yfun[k+1] << ")" << endl;
      
      // Bottom edges above
      if (xedges.status(j-1,k)!=empty) {
	cout << "[" << xedges.values(j-1,k) << "," << yfun[k] << "]";
      } else {
	cout << "[" << "          " << "," << yfun[k] << "]";
      }
      cout << "                         ";
      if (xedges.status(j-1,k+1)!=empty) {
	cout << "[" << xedges.values(j-1,k+1) << "," << yfun[k+1] 
	     << "]" << endl;
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << endl;
      }
    }
    
    // Right edge
    cout << "(" << xfun[j] << "," << yfun[k] << ")" << " ";
    cout << "[" << xfun[j] << "," << yedges.values(j,k) << "]" << " ";
    cout << "(" << xfun[j] << "," << yfun[k+1] << ")" << endl;
    
    if (j<nx-1) {
      
      // Bottom edges below
      if (xedges.status(j,k)!=empty) {
	cout << "[" << xedges.values(j,k) << "," << yfun[k] << "]";
      } else {
	cout << "[" << "          " << "," << yfun[k] << "]";
      }
      cout << "                         ";
      if (xedges.status(j,k+1)!=empty) {
	cout << "[" << xedges.values(j,k+1) << "," 
	     << yfun[k+1] << "]" << endl;
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << endl;
      }
      
      // Right points and edge below
      cout << "(" << xfun[j+1] << "," << yfun[k] << ")" << " ";
      if (yedges.status(j+1,k)!=empty) {
	cout << "[" << xfun[j+1] << "," 
	     << yedges.values(j+1,k) << "]" << " ";
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k+1] << ")" << endl;

    }
    
    cout << "found: " << found << endl;
    cout << "(" << xfun[j] << "," << yedges.values(j,k) << ")";
    if (found) {
      cout << " -> (";
      if (dir_next==dxdir) {
	cout << xedges.values(jnext,knext) << ","
	     << yfun[knext] << ")";
      } else {
	cout << xfun[jnext] << ","
	     << yedges.values(jnext,knext) << ")";
      }
    }
    cout << endl;
    
    cout << endl;
    cout.precision(6);
  }

  if (found==true) return efound;
  return enot_found;
}

int contour::find_next_point_x_direct(int j, int k, int &jnext, int &knext, 
				      int &dir_next, int nsw,
				      edge_crossings &xedges,
				      edge_crossings &yedges) {

  double closest=0.0, next;
  bool found=false;

  if (true) {

    size_t ylo_count=0, yhi_count=0;

    if (k>0) {
      if (xedges.values(j,k-1)!=empty) {
	ylo_count++;
      }
      if (yedges.status(j,k-1)!=empty) {
	ylo_count++;
      }
      if (yedges.status(j+1,k-1)!=empty) {
	ylo_count++;
      }
    }
    if (k<ny-1) {
      if (yedges.status(j,k)!=empty) {
	yhi_count++;
      }
      if (yedges.status(j+1,k)!=empty) {
	yhi_count++;
      }
      if (xedges.values(j,k+1)!=empty) {
	yhi_count++;
      }
    }
    
    if (k>0 && k<ny-1 && ylo_count+yhi_count%2==1) {
      O2SCL_ERR("Malformed edge in find_next_point_x_direct().",
		exc_esanity);
    }
  }
  
  if (k<ny-1) {
    
    double xscale=fabs(xfun[j+1]-xfun[j]);
    double yscale=fabs(yfun[k+1]-yfun[k]);
    
    // right and up
    if (yedges.status(j,k)==nsw) {
      closest=sqrt(pow(xedges.values(j,k)-xfun[j],2.0)/xscale/xscale+
		   pow(yfun[k]-yedges.values(j,k),2.0)/yscale/yscale);
      found=true;
      jnext=j;
      knext=k;
      dir_next=dydir;
    }
    
    // right and down
    if (yedges.status(j+1,k)==nsw) {
      next=sqrt(pow(xedges.values(j,k)-xfun[j+1],2.0)/xscale/xscale+
		pow(yfun[k]-yedges.values(j+1,k),2.0)/yscale/yscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j+1;
	knext=k;
	dir_next=dydir;
	closest=next;
      }
    }
    
    // right
    if (xedges.status(j,k+1)==nsw) {
      next=sqrt(pow(xedges.values(j,k)-xedges.values(j,k+1),2.0)/xscale/
		xscale+pow(yfun[k]-yfun[k+1],2.0)/yscale/yscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j;
	knext=k+1;
	dir_next=dxdir;
	closest=next;
      }
    }
  }
    
  if (k>0) {

    double xscale=fabs(xfun[j+1]-xfun[j]);
    double yscale=fabs(yfun[k]-yfun[k-1]);
    
    // left and up
    if (yedges.status(j,k-1)==nsw) {
      next=sqrt(pow(xedges.values(j,k)-xfun[j],2.0)/xscale/xscale+
		pow(yfun[k]-yedges.values(j,k-1),2.0)/yscale/yscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j;
	knext=k-1;
	dir_next=dydir;
	closest=next;
      }
    }
    
    // left and down
    if (yedges.status(j+1,k-1)==nsw) {
      next=sqrt(pow(xedges.values(j,k)-xfun[j+1],2.0)/xscale/xscale+
		pow(yfun[k]-yedges.values(j+1,k-1),2.0)/yscale/yscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j+1;
	knext=k-1;
	dir_next=dydir;
	closest=next;
      }
    }
    
    // left
    if (xedges.status(j,k-1)==nsw) {
      next=sqrt(pow(xedges.values(j,k)-xedges.values(j,k-1),2.0)/xscale/
		xscale+pow(yfun[k]-yfun[k-1],2.0)/yscale/yscale);
      if ((found==true && next<closest) || found==false) {
	found=true;
	jnext=j;
	knext=k-1;
	dir_next=dxdir;
	closest=next;
      }
    }
    
  }
  
  if (debug_next_point) {
    
    cout.precision(4);
    cout << "x edges: " << endl;
    cout << "j,k,nx,ny: " << j << " " << k << " " << nx << " " << ny << endl;
    if (k>0) {

      // Left points and edge 
      cout << "(" << xfun[j] << "," << yfun[k-1] << ")" << " ";
      if (xedges.values(j,k-1)!=empty) {
	cout << "[" << xedges.values(j,k-1) << "," 
	     << yfun[k-1] << "]" << " ";
      } else {
	cout << "[" << "          " << "," << yfun[k-1] << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k-1] << ")" << endl;

      // Left edges
      if (yedges.status(j,k-1)!=empty) {
	cout << "[" << xfun[j] << " " << yedges.values(j,k-1) << "]";
      } else {
	cout << "[" << xfun[j] << "," << "          " << "]";
      }
      cout << "                         ";
      if (yedges.status(j+1,k-1)!=empty) {
	cout << "[" << xfun[j+1] << "," << yedges.values(j+1,k-1) 
	     << "]" << endl;
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << endl;
      }
    }
    
    cout << "(" << xfun[j] << "," << yfun[k] << ")" << " ";
    cout << "[" << xedges.values(j,k) << "," << yfun[k] << "]" << " ";
    cout << "(" << xfun[j+1] << "," << yfun[k] << ")" << endl;
    
    if (k<ny-1) {

      // Right edges
      if (yedges.status(j,k)!=empty) {
	cout << "[" << xfun[j] << " " << yedges.values(j,k) << "]";
      } else {
	cout << "[" << xfun[j] << "," << "          " << "]";
      }
      cout << "                         ";
      if (yedges.status(j+1,k)!=empty) {
	cout << "[" << xfun[j+1] << "," << yedges.values(j+1,k) 
	     << "]" << endl;
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << endl;
      }

      // Right points and edge
      cout << "(" << xfun[j] << "," << yfun[k+1] << ")" << " ";
      if (xedges.values(j,k+1)!=empty) {
	cout << "[" << xedges.values(j,k+1) << "," 
	     << yfun[k+1] << "]" << " ";
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k+1] << ")" << endl;

    }
    
    cout << "found: " << found << endl;
    cout << "(" << xedges.values(j,k) << "," << yfun[k] << ")";
    if (found) {
      cout << " -> (";
      if (dir_next==dxdir) {
	cout << xedges.values(jnext,knext) << ","
	     << yfun[knext] << ")";
      } else {
	cout << xfun[jnext] << ","
	     << yedges.values(jnext,knext) << ")";
      }
    }
    cout << endl;
    
    cout << endl;
    cout.precision(6);
  }
  
  if (found==true) return efound;
  return enot_found;
}

void contour::find_intersections(size_t ilev, double &level,
				 edge_crossings &xedges, 
				 edge_crossings &yedges) {
  
  // Adjust the specified contour level to ensure none of the data
  // points is exactly on a contour
  bool level_corner;
  do {
    
    // Look for a match
    level_corner=false;
    for(int k=0;k<ny;k++) {
      for(int j=0;j<nx;j++) {
	if (data(j,k)==level) {
	  level_corner=true;
	}
      }
    }

    // If we found a match, adjust the level
    if (level_corner==true) {
      if (verbose>0) {
	cout << "Found intersection of contour level and corner." << endl;
	cout << "Adjusting level " << level << " to ";
      }
      if (nlev==1) {
	// If there's only one contour level, then make a
	// simple adjustment
	level*=1.0+lev_adjust;
      } else {
	// Otherwise, compute the appropriate adjustment by 
	// finding the closest contour level not equal to the
	// current contour level
	double diff;
	int iclosest=0;
	if (ilev==0) {
	  diff=fabs(level-levels[1]);
	  iclosest=1;
	} else {
	  diff=fabs(level-levels[0]);
	  iclosest=0;
	}
	for(int ik=0;ik<nlev;ik++) {
	  if (ik!=((int)ilev)) {
	    if (fabs(level-levels[ik])<diff) {
	      iclosest=ik;
	      diff=fabs(level-levels[ik]);
	    }
	  }
	}
	level+=fabs(level-levels[iclosest])*lev_adjust;
      }
      if (verbose>0) cout << level << endl;
    }

  } while (level_corner==true);
  
  // Find all level crossings
  for(int k=0;k<ny;k++) {
    for(int j=0;j<nx;j++) {
      if (j<nx-1) {
	xedges.status(j,k)=empty;
	if ((data(j,k)-level)*(data(j+1,k)-level)<0.0) {
	  xedges.status(j,k)=edge;
	  if (verbose>1) {
	    cout << "Vertical edge for level   " << level << " between (" 
		 << k << "," << j << ") and (" << k << "," 
		 << j+1 << ")" << endl;
	  }
	} 
      }
      if (k<ny-1) {
	yedges.status(j,k)=empty;
	if ((data(j,k)-level)*(data(j,k+1)-level)<0.0) {
	  if (verbose>1) {
	    cout << "Horizontal edge for level " << level << " between (" 
		 << k << "," << j << ") and (" << k+1 << "," 
		 << j << ")" << endl;
	  }
	  yedges.status(j,k)=edge;
	}
      }
    }
  }

  return;
}

void contour::edges_in_y_direct(double level, interp<ubvector> &si,
				edge_crossings &yedges) {

  for(int k=0;k<ny-1;k++) {
    for(int j=0;j<nx;j++) {
      
      // For each iright edge
      if (yedges.status(j,k)==edge) {
		
	// Find two points for linear interpolation
	int ileft=k;
	if (ileft==ny-1) ileft--;
	int nint=2;
	
	ubvector xi(nint), yi(nint);
	
	xi[0]=data(j,k);
	xi[1]=data(j,k+1);
	yi[0]=yfun[ileft];
	yi[1]=yfun[ileft+1];
	
	// Do the interpolation and set the edge appropriately
	yedges.values(j,k)=si.eval(level,nint,xi,yi);
	
	if (verbose>1) {
	  cout << "Horizontal edge: (" << k << "," << j << ") -> ("
	       << k+1 << "," << j << ")" << endl;
	  cout << " coords: " << yfun[ileft] << " "
	       << yedges.values(j,k) << " " << yfun[ileft+1] << endl;
	  cout << "   data: " << data(j,k) << " "
	       << level << " " << data(j,k+1) << endl;
	}
      } else {
	yedges.values(j,k)=0.0;
      }
    }
  }

  return;
}

void contour::edges_in_x_direct(double level, interp<ubvector> &si,
				edge_crossings &xedges) {

  for(int k=0;k<ny;k++) {
    for(int j=0;j<nx-1;j++) {
      
      // For each bottom edge
      if (xedges.status(j,k)==edge) {
	
	// Find two points for linear interpolation
	int ileft=j;
	if (ileft==nx-1) ileft--;
	int nint=2;
	
	ubvector xi(nint), yi(nint);
	
	xi[0]=data(j,k);
	xi[1]=data(j+1,k);
	yi[0]=xfun[ileft];
	yi[1]=xfun[ileft+1];

	// Do the interpolation and set the edge appropriately
	xedges.values(j,k)=si.eval(level,nint,xi,yi);
	
	if (verbose>1) {
	  cout << "Vertical edge:   (" << k << "," << j << ") -> ("
	       << k+1 << "," << j << ")" << endl;
	  cout << " coords: " << xfun[ileft] << " "
	       << xedges.values(j,k) << " " << xfun[ileft+1] << endl;
	  cout << "   data: " << data(j,k) << " "
	       << level << " " << data(j+1,k) << endl;
	}
      } else {
	xedges.values(j,k)=0.0;
      }
    }
  }

  return;
}

void contour::process_line(int j, int k, int dir, std::vector<double> &x, 
			   std::vector<double> &y, bool first, 
			   edge_crossings &xedges, edge_crossings &yedges) {
  
  std::vector<double> xt, yt;

  int fp, jnext, knext, dir_next;
  bool zero_points;

  // If we've found a new intersection
  if (dir==dydir) {
    fp=find_next_point_y_direct(j,k,jnext,knext,dir_next,edge,xedges,yedges);
  } else {
    fp=find_next_point_x_direct(j,k,jnext,knext,dir_next,edge,xedges,yedges);
  }
  if (fp==enot_found) zero_points=true;
  else zero_points=false;
  while (fp==efound) {

    j=jnext;
    k=knext;
    dir=dir_next;
    if (dir==dydir) {
      if (verbose>0) {
	cout << "(" << xfun[j] << ", " << yedges.values(j,k) << ")" << endl;
      }

      if (first) {
	x.push_back(xfun[j]);
	y.push_back(yedges.values(j,k));
      } else {
	xt.push_back(xfun[j]);
	yt.push_back(yedges.values(j,k));
      }

      yedges.status(j,k)=contourp;
      fp=find_next_point_y_direct(j,k,jnext,knext,dir_next,edge,xedges,yedges);
    } else {
      if (verbose>0) {
	cout << "(" << xedges.values(j,k) << ", " << yfun[k] << ")" << endl;
      }

      if (first) {
	x.push_back(xedges.values(j,k));
	y.push_back(yfun[k]);
      } else {
	xt.push_back(xedges.values(j,k));
	yt.push_back(yfun[k]);
      }

      xedges.status(j,k)=contourp;
      fp=find_next_point_x_direct(j,k,jnext,knext,dir_next,edge,xedges,yedges);
    }
  }
  
  // Set the last point to an endpoint
  if (dir==dydir) {
    yedges.status(j,k)=endpoint;
  } else {
    xedges.status(j,k)=endpoint;
  }
  
  if (first==false) {

    std::vector<double> xt2, yt2;

    for(int i=xt.size()-1;i>=0;i--) {
      xt2.push_back(xt[i]);
      yt2.push_back(yt[i]);
    }
    for(size_t i=0;i<x.size();i++) {
      xt2.push_back(x[i]);
      yt2.push_back(y[i]);
    }
    x=xt2;
    y=yt2;

  }

  // If this is the second half of the line, see if we
  // need to connect two endpoints
  if (first==false && x.size()>0) {

    // Look for an endpoint next to the last endpoint
    if (dir==dydir) {
      fp=find_next_point_y_direct(j,k,jnext,knext,dir_next,endpoint,
				  xedges,yedges);
    } else {
      fp=find_next_point_x_direct(j,k,jnext,knext,dir_next,endpoint,
				  xedges,yedges);
    }
    
    // If we found two adjacent end-points, close the contour. We 
    // also ensure the size is greater than two in order to avoid
    // erroneously closing a line with only one segment.
    if (fp==1 && x.size()>2) {

      // We found a connection, so change the edge status
      // of the two endpoints to internal contour points
      if (dir==dydir) {
	yedges.status(j,k)=contourp;
      } else {
	xedges.status(j,k)=contourp;
      }
      if (dir_next==dydir) {
	yedges.status(jnext,knext)=contourp;
      } else {
	xedges.status(jnext,knext)=contourp;
      }
      
      // Add the first point to the end to close the contour
      x.push_back(x[0]);
      y.push_back(y[0]);

    }

  }
  
  if (verbose>1) {
    cout << "Endpoint reached." << endl;
    if (verbose>2) {
      char ch;
      cin >> ch;
    }
  }
  return;
}

void contour::calc_contours(std::vector<contour_line> &clines) {

  // Check that we're ready
  if (nx==0) {
    O2SCL_ERR("Data not set in calc_contours().",exc_einval);
    return;
  }
  if (levels_set==false) {
    O2SCL_ERR("Contour levels not set in calc_contours().",exc_einval);
    return;
  }

  // Clear edge storage
  yed.clear();
  xed.clear();

  // Clear contour lines object
  clines.clear();
  
  // The interpolation object (only works with linear interpolation
  // at the moment)
  interp<ubvector> oi(itp_linear);
  
  // For each level
  for(int i=0;i<nlev;i++) {

    // Make space for the edges
    edge_crossings xedges, yedges;
    xedges.status.resize(nx-1,ny);
    xedges.values.resize(nx-1,ny);
    yedges.status.resize(nx,ny-1);
    yedges.values.resize(nx,ny-1);

    if (verbose>1) {
      std::cout << "\nLooking for edges for level: " 
		<< levels[i] << std::endl;
    }

    // Examine the each of the rows for an intersection
    find_intersections(i,levels[i],xedges,yedges);
	
    if (verbose>1) {
      std::cout << "\nInterpolating edge intersections for level: " 
		<< levels[i] << std::endl;
    }

    // Process the edges in the x direction
    edges_in_x_direct(levels[i],oi,xedges);

    // Process the edges in the y direction
    edges_in_y_direct(levels[i],oi,yedges);
    
    if (verbose>1) {
      std::cout << "\nPiecing together contour lines for level: " 
		<< levels[i] << std::endl;
    }

    // Now go through and one side of the line
    bool foundline=true;
    while(foundline==true) {
      for(int j=0;j<nx;j++) {
	for(int k=0;k<ny;k++) {
	  foundline=false;

	  contour_line c;
	  c.level=levels[i];
	      
	  // A line beginning with a right edge
	  if (k<ny-1 && yedges.status(j,k)==edge) {
	    if (verbose>0) {
	      std::cout << "Starting contour line for level "
			<< levels[i] << ":" << std::endl;
	      std::cout << "(" << xfun[j] << ", " << yedges.values(j,k) 
			<< ")" << std::endl;
	    }
	    c.x.push_back(xfun[j]);
	    c.y.push_back(yedges.values(j,k));
	    yedges.status(j,k)++;
	    
	    // Go through both sides
	    process_line(j,k,dydir,c.x,c.y,true,xedges,yedges);
	    if (verbose>0) {
	      std::cout << "Computing other side of line." << std::endl;
	    }
	    process_line(j,k,dydir,c.x,c.y,false,xedges,yedges);
	    foundline=true;
	  }

	  // A line beginning with a bottom edge
	  if (j<nx-1 && foundline==false && xedges.status(j,k)==edge) {
	    if (verbose>0) {
	      std::cout << "Starting contour line for level "
			<< levels[i] << ":" << std::endl;
	      std::cout << "(" << xedges.values(j,k) << ", " << yfun[k] 
			<< ")" << std::endl;
	    }
	    c.x.push_back(xedges.values(j,k));
	    c.y.push_back(yfun[k]);
	    xedges.status(j,k)++;
	    
	    // Go through both sides
	    process_line(j,k,dxdir,c.x,c.y,true,xedges,yedges);
	    if (verbose>0) {
	      std::cout << "Computing other side of line." << std::endl;
	    }
	    process_line(j,k,dxdir,c.x,c.y,false,xedges,yedges);
	    foundline=true;
	  }

	  // Add line to list
	  if (foundline==true) {
	    clines.push_back(c);
	  }
	}
      }

    }

    // Store edge information from this level
    xed.push_back(xedges);
    yed.push_back(yedges);

    if (verbose>0) {
      std::cout << "Processing next level." << std::endl;
    }
  }
  
  return;
}

void contour::print_edges_yhoriz(edge_crossings &xedges,
				 edge_crossings &yedges) {
  
  size_t ny2=xedges.status.size2();
  size_t nx2=yedges.status.size1();
  cout << " ";
  for(size_t j=0;j<ny2;j++) {
    cout << j%10 << " ";
  }
  cout << endl;
  for(int i=nx2-1;i>=0;i--) {
    if (i<((int)(nx2))-1) {
      cout << " ";
      for(size_t j=0;j<ny2;j++) {
	if (xedges.status(i,j)==empty) cout << "|";
	else if (xedges.status(i,j)==edge) cout << "+";
	else if (xedges.status(i,j)==contourp) cout << "*";
	else cout << "E";
	cout << " ";
      }
      cout << endl;
    }
    cout << i%10;
    for(size_t j=0;j<ny2;j++) {
      if (j<ny2-1) {
	cout << " ";
	if (yedges.status(i,j)==empty) cout << "-";
	else if (yedges.status(i,j)==edge) cout << "+";
	else if (yedges.status(i,j)==contourp) cout << "*";
	else cout << "E";
      }
    }
    cout << " " << i%10 << endl;
  }
  cout << " ";
  for(size_t j=0;j<ny2;j++) {
    cout << j%10 << " ";
  }
  cout << endl;

  return;
}

void contour::print_edges_xhoriz(edge_crossings &xedges,
				 edge_crossings &yedges) {
  
  size_t ny2=xedges.status.size2();
  size_t nx2=yedges.status.size1();
  cout << " ";
  for(size_t i=0;i<nx2;i++) {
    cout << i%10 << " ";
  }
  cout << endl;
  for(int j=ny2-1;j>=0;j--) {
    if (j<((int)ny2)-1) {
      cout << " ";
      for(size_t i=0;i<nx2;i++) {
	if (yedges.status(i,j)==empty) cout << "|";
	else if (yedges.status(i,j)==edge) cout << "+";
	else if (yedges.status(i,j)==contourp) cout << "*";
	else cout << "E";
	cout << " ";
      }
      cout << endl;
    }
    cout << j%10;
    for(size_t i=0;i<nx2;i++) {
      if (i<nx2-1) {
	cout << " ";
	if (xedges.status(i,j)==empty) cout << "-";
	else if (xedges.status(i,j)==edge) cout << "+";
	else if (xedges.status(i,j)==contourp) cout << "*";
	else cout << "E";
      }
    }
    cout << " " << j%10 << endl;
  }
  cout << " ";
  for(size_t i=0;i<nx2;i++) {
    cout << i%10 << " ";
  }
  cout << endl;
  
  return;
}

