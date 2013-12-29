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

  interp_vec<ubvector> *si;

  // Create new xfun
  ubvector xidx(nx);
  xfun.resize(newx);
  for(int i=0;i<nx;i++) xidx[i]=((double)i);
  si=new interp_vec<ubvector>(nx,xidx,xfun_old,interp_type);
  for(int i=0;i<newx;i++) {
    xfun[i]=si->eval(((double)i)/(newx-1+xfact)*nx);
  }
  delete si;

  // Create new yfun
  ubvector yidx(ny);
  yfun.resize(newy);
  for(int i=0;i<ny;i++) yidx[i]=((double)i);
  si=new interp_vec<ubvector>(ny,yidx,yfun_old,interp_type);
  for(int i=0;i<newy;i++) {
    yfun[i]=si->eval(((double)i)/(newy-1+yfact)*ny);
  }
  delete si;
  
  // Perform the first round, expanding the number of rows from 
  // nx -> newx 
  ubmatrix bc(newx,ny);
  for(int i=0;i<newy;i+=yfact) {
    // Get the column of data_old with index i/yfact
    ubvector at(nx);
    for(int ij=0;ij<nx;ij++) at[ij]=data_old(ij,i/yfact);
    // Use it to interpolate a column into bc
    si=new interp_vec<ubvector>(nx,xfun_old,at,interp_type);
    for(int j=0;j<newx;j++) {
      bc(j,i/yfact)=si->eval(xfun[j]);
    }
    delete si;
  }
  
  // Perform the second round, expanding the number of columns
  // from ny -> newy
  data.resize(newx,newy);
  for(int i=0;i<newx;i++) {
    ubvector at2(ny);
    for(int ij=0;ij<ny;ij++) at2[ij]=bc(i,ij);
    si=new interp_vec<ubvector>(ny,yfun_old,at2,interp_type);
    for(int j=0;j<newy;j++) {
      data(i,j)=si->eval(yfun[j]);
    }
    delete si;
  }

  // Modify size
  nx=newx;
  ny=newy;
  
  return;
}

int contour::find_next_point_right(int j, int k, int &jnext, int &knext, 
				   int &dir_next, int nsw,
				   edge_crossings &right,
				   edge_crossings &bottom) {

  double closest=0.0, next;
  bool found=false;

  if (true) {

    size_t top_count=0, bottom_count=0;
    if (j>=1) {
      if (right.status(k,j-1)!=empty) {
	top_count++;
      }
      if (bottom.status(k,j-1)!=empty) {
	top_count++;
      }
      if (bottom.status(k+1,j-1)!=empty) {
	top_count++;
      }
    }
    if (j<nx-1) {
      if (bottom.status(k,j)!=empty) {
	bottom_count++;
      }
      if (bottom.status(k+1,j)!=empty) {
	bottom_count++;
      }
      if (right.status(k,j+1)!=empty) {
	bottom_count++;
      }
    }

    if (j>0 && j<nx-1 && bottom_count+top_count%2==1) {
      O2SCL_ERR("Malformed edge in find_next_point_right().",exc_efailed);
    }

  }

  if (false) {

    cout.precision(4);
    cout << "right: " << endl;
    cout << "j,k,nx,ny: " << j << " " << k << " " << nx << " " << ny << endl;
    if (j>=1) {

      // Top points and edge above
      cout << "(" << xfun[j-1] << "," << yfun[k] << ")" << " ";
      if (right.status(k,j-1)!=empty) {
	cout << "[" << xfun[j-1] << "," 
	     << right.values(k,j-1) << "]" << " ";
      } else {
	cout << "[" << xfun[j-1] << "," << "          " << "]" << " ";
      }
      cout << "(" << xfun[j-1] << "," << yfun[k+1] << ")" << endl;
      
      // Bottom edges above
      if (bottom.status(k,j-1)!=empty) {
	cout << "[" << bottom.values(k,j-1) << "," << yfun[k] << "]";
      } else {
	cout << "[" << "          " << "," << yfun[k] << "]";
      }
      cout << "                         ";
      if (bottom.status(k+1,j-1)!=empty) {
	cout << "[" << bottom.values(k+1,j-1) << "," << yfun[k+1] 
	     << "]" << endl;
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << endl;
      }
    }
    
    // Right edge
    cout << "(" << xfun[j] << "," << yfun[k] << ")" << " ";
    cout << "[" << xfun[j] << "," << right.values(k,j) << "]" << " ";
    cout << "(" << xfun[j] << "," << yfun[k+1] << ")" << endl;
    
    if (j<nx-1) {
      
      // Bottom edges below
      if (bottom.status(k,j)!=empty) {
	cout << "[" << bottom.values(k,j) << "," << yfun[k] << "]";
      } else {
	cout << "[" << "          " << "," << yfun[k] << "]";
      }
      cout << "                         ";
      if (bottom.status(k+1,j)!=empty) {
	cout << "[" << bottom.values(k+1,j) << "," 
	     << yfun[k+1] << "]" << endl;
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << endl;
      }
      
      // Right points and edge below
      cout << "(" << xfun[j+1] << "," << yfun[k] << ")" << " ";
      if (right.status(k,j+1)!=empty) {
	cout << "[" << xfun[j+1] << "," 
	     << right.values(k,j+1) << "]" << " ";
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k+1] << ")" << endl;

    }
    
    cout << endl;
    cout.precision(6);
  }

  // down and left
  if (j<nx-1 && bottom.status(k,j)==nsw) {
    closest=sqrt(pow(right.values(k,j)-yfun[k],2.0)+
		 pow(xfun[j]-bottom.values(k,j),2.0));
    found=true;
    jnext=j;
    knext=k;
    dir_next=dbottom;
  }
  // down
  if (j<nx-1 && right.status(k,j+1)==nsw) {
    next=sqrt(pow(right.values(k,j)-right.values(k,j+1),2.0)+
	      pow(xfun[j]-xfun[j+1],2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j+1;
      knext=k;
      dir_next=dright;
      closest=next;
    }
  }
  // down and right
  if (j<nx-1 && bottom.status(k+1,j)==nsw) {
    next=sqrt(pow(right.values(k,j)-yfun[k+1],2.0)+
	      pow(xfun[j]-bottom.values(k+1,j),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j;
      knext=k+1;
      dir_next=dbottom;
      closest=next;
    }
  }
  // up
  if (j>0 && right.status(k,j-1)==nsw) { 
    next=sqrt(pow(right.values(k,j)-right.values(k,j-1),2.0)+
	      pow(xfun[j]-xfun[j-1],2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j-1;
      knext=k;
      dir_next=dright;
      closest=next;
    }
  }
  // up and left
  if (j>0 && bottom.status(k,j-1)==nsw) {
    next=sqrt(pow(right.values(k,j)-yfun[k],2.0)+
	      pow(xfun[j]-bottom.values(k,j-1),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j-1;
      knext=k;
      dir_next=dbottom;
      closest=next;
    }
  }
  // up and right
  if (j>0 && bottom.status(k+1,j-1)==nsw) {
    next=sqrt(pow(right.values(k,j)-yfun[k+1],2.0)+
	      pow(xfun[j]-bottom.values(k+1,j-1),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j-1;
      knext=k+1;
      dir_next=dbottom;
      closest=next;
    }
  }
  if (found==true) return efound;
  return enot_found;
}

int contour::find_next_point_bottom(int j, int k, int &jnext, int &knext, 
				    int &dir_next, int nsw,
				    edge_crossings &right,
				    edge_crossings &bottom) {

  double closest=0.0, next;
  bool found=false;

  if (true) {

    size_t left_count=0, right_count=0;

    if (k>0) {
      if (bottom.values(k-1,j)!=empty) {
	left_count++;
      }
      if (right.status(k-1,j)!=empty) {
	left_count++;
      }
      if (right.status(k-1,j+1)!=empty) {
	left_count++;
      }
    }
    if (k<ny-1) {
      if (right.status(k,j)!=empty) {
	right_count++;
      }
      if (right.status(k,j+1)!=empty) {
	right_count++;
      }
      if (bottom.values(k+1,j)!=empty) {
	right_count++;
      }
    }
    
    if (k>0 && k<ny-1 && left_count+right_count%2==1) {
      O2SCL_ERR("Malformed edge in find_next_point_bottom().",exc_efailed);
    }
  }

  if (false) {

    cout.precision(4);
    cout << "bottom: " << endl;
    cout << "j,k,nx,ny: " << j << " " << k << " " << nx << " " << ny << endl;
    if (k>0) {

      // Left points and edge 
      cout << "(" << xfun[j] << "," << yfun[k-1] << ")" << " ";
      if (bottom.values(k-1,j)!=empty) {
	cout << "[" << bottom.values(k-1,j) << "," 
	     << yfun[k-1] << "]" << " ";
      } else {
	cout << "[" << "          " << "," << yfun[k-1] << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k-1] << ")" << endl;

      // Left edges
      if (right.status(k-1,j)!=empty) {
	cout << "[" << xfun[j] << " " << right.values(k-1,j) << "]";
      } else {
	cout << "[" << xfun[j] << "," << "          " << "]";
      }
      cout << "                         ";
      if (right.status(k-1,j+1)!=empty) {
	cout << "[" << xfun[j+1] << "," << right.values(k-1,j+1) 
	     << "]" << endl;
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << endl;
      }
    }
    
    cout << "(" << xfun[j] << "," << yfun[k] << ")" << " ";
    cout << "[" << bottom.values(k,j) << "," << yfun[k] << "]" << " ";
    cout << "(" << xfun[j+1] << "," << yfun[k] << ")" << endl;
    
    if (k<ny-1) {

      // Right edges
      if (right.status(k,j)!=empty) {
	cout << "[" << xfun[j] << " " << right.values(k,j) << "]";
      } else {
	cout << "[" << xfun[j] << "," << "          " << "]";
      }
      cout << "                         ";
      if (right.status(k,j+1)!=empty) {
	cout << "[" << xfun[j+1] << "," << right.values(k,j+1) 
	     << "]" << endl;
      } else {
	cout << "[" << xfun[j+1] << "," << "          " << "]" << endl;
      }

      // Right points and edge
      cout << "(" << xfun[j] << "," << yfun[k+1] << ")" << " ";
      if (bottom.values(k+1,j)!=empty) {
	cout << "[" << bottom.values(k+1,j) << "," 
	     << yfun[k+1] << "]" << " ";
      } else {
	cout << "[" << "          " << "," << yfun[k+1] << "]" << " ";
      }
      cout << "(" << xfun[j+1] << "," << yfun[k+1] << ")" << endl;

    }
    
    cout << endl;
    cout.precision(6);
  }

  // right and up
  if (k<ny-1 && right.status(k,j)==nsw) {
    closest=sqrt(pow(bottom.values(k,j)-xfun[j],2.0)+
		 pow(yfun[k]-right.values(k,j),2.0));
    found=true;
    jnext=j;
    knext=k;
    dir_next=dright;
  }
  // right and down
  if (k<ny-1 && right.status(k,j+1)==nsw) {
    next=sqrt(pow(bottom.values(k,j)-xfun[j+1],2.0)+
	      pow(yfun[k]-right.values(k,j+1),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j+1;
      knext=k;
      dir_next=dright;
      closest=next;
    }
  }
  // right
  if (k<ny-1 && bottom.status(k+1,j)==nsw) {
    next=sqrt(pow(bottom.values(k,j)-bottom.values(k+1,j),2.0)+
	      pow(yfun[k]-yfun[k+1],2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j;
      knext=k+1;
      dir_next=dbottom;
      closest=next;
    }
  }
  // left and up
  if (k>0 && right.status(k-1,j)==nsw) {
    next=sqrt(pow(bottom.values(k,j)-xfun[j],2.0)+
	      pow(yfun[k]-right.values(k-1,j),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j;
      knext=k-1;
      dir_next=dright;
      closest=next;
    }
  }
  // left and down
  if (k>0 && right.status(k-1,j+1)==nsw) {
    next=sqrt(pow(bottom.values(k,j)-xfun[j+1],2.0)+
	      pow(yfun[k]-right.values(k-1,j+1),2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j+1;
      knext=k-1;
      dir_next=dright;
      closest=next;
    }
  }
  // left
  if (k>0 && bottom.status(k-1,j)==nsw) {
    next=sqrt(pow(bottom.values(k,j)-bottom.values(k-1,j),2.0)+
	      pow(yfun[k]-yfun[k-1],2.0));
    if ((found==true && next<closest) || found==false) {
      found=true;
      jnext=j;
      knext=k-1;
      dir_next=dbottom;
      closest=next;
    }
  }
  
  if (found==true) return efound;
  return enot_found;
}

void contour::find_intersections(size_t ilev, double &level,
				edge_crossings &right, 
				edge_crossings &bottom) {

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
	bottom.status(k,j)=empty;
	if ((data(j,k)-level)*(data(j+1,k)-level)<0.0) {
	  bottom.status(k,j)=edge;
	  if (verbose>1) {
	    cout << "Vertical edge for level   " << level << " between (" 
		 << k << "," << j << ") and (" << k << "," 
		 << j+1 << ")" << endl;
	  }
	} 
      }
      if (k<ny-1) {
	right.status(k,j)=empty;
	if ((data(j,k)-level)*(data(j,k+1)-level)<0.0) {
	  if (verbose>1) {
	    cout << "Horizontal edge for level " << level << " between (" 
		 << k << "," << j << ") and (" << k+1 << "," 
		 << j << ")" << endl;
	  }
	  right.status(k,j)=edge;
	}
      }
    }
  }

  //print_edges(right,bottom);

  return;
}

void contour::right_edges(double level, 
			 interp<ubvector> &si,
			 edge_crossings &right) {

  for(int k=0;k<ny-1;k++) {
    for(int j=0;j<nx;j++) {
      
      // For each iright edge
      if (right.status(k,j)==edge) {
		
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
	right.values(k,j)=si.eval(level,nint,xi,yi);
	
	if (verbose>1) {
	  cout << "Horizontal edge: (" << k << "," << j << ") -> ("
	       << k+1 << "," << j << ")" << endl;
	  cout << " coords: " << yfun[ileft] << " "
	       << right.values(k,j) << " " << yfun[ileft+1] << endl;
	  cout << "   data: " << data(j,k) << " "
	       << level << " " << data(j,k+1) << endl;
	}
      } else {
	right.values(k,j)=0.0;
      }
    }
  }

  return;
}

void contour::bottom_edges(double level, 
			  interp<ubvector> &si,
			  edge_crossings &bottom) {

  for(int k=0;k<ny;k++) {
    for(int j=0;j<nx-1;j++) {
      
      // For each bottom edge
      if (bottom.status(k,j)==edge) {
	
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
	bottom.values(k,j)=si.eval(level,nint,xi,yi);
	
	if (verbose>1) {
	  cout << "Vertical edge:   (" << k << "," << j << ") -> ("
	       << k+1 << "," << j << ")" << endl;
	  cout << " coords: " << xfun[ileft] << " "
	       << bottom.values(k,j) << " " << xfun[ileft+1] << endl;
	  cout << "   data: " << data(j,k) << " "
	       << level << " " << data(j+1,k) << endl;
	}
      } else {
	bottom.values(k,j)=0.0;
      }
    }
  }

  return;
}

void contour::process_line(int j, int k, int dir, std::vector<double> &x, 
			   std::vector<double> &y, bool first, 
			   edge_crossings &right,
			   edge_crossings &bottom) {
  std::vector<double> xt, yt;

  int fp, jnext, knext, dir_next;
  bool zero_points;

  // If we've found a new intersection
  if (dir==dright) {
    fp=find_next_point_right(j,k,jnext,knext,dir_next,edge,right,bottom);
  } else {
    fp=find_next_point_bottom(j,k,jnext,knext,dir_next,edge,right,bottom);
  }
  if (fp==enot_found) zero_points=true;
  else zero_points=false;
  while (fp==efound) {

    //print_edges(right,bottom);

    j=jnext;
    k=knext;
    dir=dir_next;
    if (dir==dright) {
      if (verbose>0) {
	cout << "(" << xfun[j] << ", " << right.values(k,j) << ")" << endl;
      }

      if (first) {
	x.push_back(xfun[j]);
	y.push_back(right.values(k,j));
      } else {
	xt.push_back(xfun[j]);
	yt.push_back(right.values(k,j));
      }

      right.status(k,j)=contourp;
      fp=find_next_point_right(j,k,jnext,knext,dir_next,edge,right,bottom);
    } else {
      if (verbose>0) {
	cout << "(" << bottom.values(k,j) << ", " << yfun[k] << ")" << endl;
      }

      if (first) {
	x.push_back(bottom.values(k,j));
	y.push_back(yfun[k]);
      } else {
	xt.push_back(bottom.values(k,j));
	yt.push_back(yfun[k]);
      }

      bottom.status(k,j)=contourp;
      fp=find_next_point_bottom(j,k,jnext,knext,dir_next,edge,right,bottom);
    }
  }
  
  // Set the last point to an endpoint
  if (dir==dright) {
    right.status(k,j)=endpoint;
  } else {
    bottom.status(k,j)=endpoint;
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
    if (dir==dright) {
      fp=find_next_point_right(j,k,jnext,knext,dir_next,endpoint,
			       right,bottom);
    } else {
      fp=find_next_point_bottom(j,k,jnext,knext,dir_next,endpoint,
				right,bottom);
    }
    
    // If we found two adjacent end-points, close the contour. We 
    // also ensure the size is greater than two in order to avoid
    // erroneously closing a line with only one segment.
    if (fp==1 && x.size()>2) {

      // We found a connection, so change the edge status
      // of the two endpoints to internal contour points
      if (dir==dright) {
	right.status(k,j)=contourp;
      } else {
	bottom.status(k,j)=contourp;
      }
      if (dir_next==dright) {
	right.status(knext,jnext)=contourp;
      } else {
	bottom.status(knext,jnext)=contourp;
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
  red.clear();
  bed.clear();

  // The interpolation object (only works with linear interpolation
  // at the moment)
  interp<ubvector> oi(itp_linear);
  
  // For each level
  for(int i=0;i<nlev;i++) {

    // Make space for the edges
    edge_crossings right, bottom;
    right.status.resize(ny-1,nx);
    right.values.resize(ny-1,nx);
    bottom.status.resize(ny,nx-1);
    bottom.values.resize(ny,nx-1);

    if (verbose>1) {
      std::cout << "\nLooking for edges for level: " 
		<< levels[i] << std::endl;
    }

    // Examine the each of the rows for an intersection
    find_intersections(i,levels[i],right,bottom);
	
    if (verbose>1) {
      std::cout << "\nInterpolating edge intersections for level: " 
		<< levels[i] << std::endl;
    }

    // Process the right edges
    right_edges(levels[i],oi,right);
    
    // Process the bottom edges
    bottom_edges(levels[i],oi,bottom);

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
	  if (k<ny-1 && right.status(k,j)==edge) {
	    if (verbose>0) {
	      std::cout << "Starting contour line for level "
			<< levels[i] << ":" << std::endl;
	      std::cout << "(" << xfun[j] << ", " << right.values(k,j) 
			<< ")" << std::endl;
	    }
	    c.x.push_back(xfun[j]);
	    c.y.push_back(right.values(k,j));
	    right.status(k,j)++;
	    
	    // Go through both sides
	    process_line(j,k,dright,c.x,c.y,true,right,bottom);
	    if (verbose>0) {
	      std::cout << "Computing other side of line." << std::endl;
	    }
	    process_line(j,k,dright,c.x,c.y,false,right,bottom);
	    foundline=true;
	  }

	  // A line beginning with a bottom edge
	  if (j<nx-1 && foundline==false && bottom.status(k,j)==edge) {
	    if (verbose>0) {
	      std::cout << "Starting contour line for level "
			<< levels[i] << ":" << std::endl;
	      std::cout << "(" << bottom.values(k,j) << ", " << yfun[k] 
			<< ")" << std::endl;
	    }
	    c.x.push_back(bottom.values(k,j));
	    c.y.push_back(yfun[k]);
	    bottom.status(k,j)++;
	    
	    // Go through both sides
	    process_line(j,k,dbottom,c.x,c.y,true,right,bottom);
	    if (verbose>0) {
	      std::cout << "Computing other side of line." << std::endl;
	    }
	    process_line(j,k,dbottom,c.x,c.y,false,right,bottom);
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
    red.push_back(right);
    bed.push_back(bottom);

    if (verbose>0) {
      std::cout << "Processing next level." << std::endl;
    }
  }
  
  return;
}

void contour::print_edges(edge_crossings &right,
			  edge_crossings &bottom) {
  
  size_t ny2=bottom.status.size1();
  size_t nx2=right.status.size2();
  for(size_t j=0;j<nx2;j++) {
    for(size_t i=0;i<ny2;i++) {
      if (i<ny2-1) {
	cout << " ";
	if (right.status(i,j)==empty) cout << "_";
	else if (right.status(i,j)==edge) cout << "+";
	else if (right.status(i,j)==contourp) cout << "*";
	else cout << "E";
      }
    }
    cout << endl;
    for(size_t i=0;i<ny2;i++) {
      if (j<nx2-1) {
	if (bottom.status(i,j)==empty) cout << "|";
	else if (bottom.status(i,j)==edge) cout << "+";
	else if (bottom.status(i,j)==contourp) cout << "*";
	else cout << "E";
	cout << " ";
      }
    }
    cout << endl;
  }

  return;
}

