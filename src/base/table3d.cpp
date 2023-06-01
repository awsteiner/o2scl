/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/table3d.h>
#include <o2scl/hist_2d.h>
#include <o2scl/vec_stats.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

hist_2d table3d::to_hist_2d(std::string slice, int verbose) {
  hist_2d h;
  
  ubvector bin_edges_x, bin_edges_y;

  vector_to_bins(xval,bin_edges_x,verbose);
  vector_to_bins(yval,bin_edges_y,verbose);

  h.set_bin_edges(numx+1,bin_edges_x,numy+1,bin_edges_y);

  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      h.set_wgt_i(i,j,get(i,j,slice));
    }
  }
  
  return h;
}

hist table3d::to_hist(std::string slice, size_t n_bins, int verbose) {
  
  hist h;
  
  const ubmatrix &sl=this->get_slice(slice);
  
  double min=1.0, max=0.0;
  for(size_t i=0;i<this->get_nx();i++) {
    for(size_t j=0;j<this->get_ny();j++) {
      if (i==0 && j==0) {
        min=sl(i,j);
        max=sl(i,j);
      } else if (sl(i,j)<min) {
        min=sl(i,j);
      } else if (sl(i,j)>max) {
        max=sl(i,j);
      }
    }
  }
  
  uniform_grid<double> ug=uniform_grid_end<double>(min,max,n_bins);
  h.set_bin_edges(ug);
  
  for(size_t i=0;i<this->get_nx();i++) {
    for(size_t j=0;j<this->get_ny();j++) {
      h.update(sl(i,j));
    }
  }
  
  return h;
}

table3d::table3d() {
  xy_set=false;
  size_set=false;
  has_slice=false;
  itype=itp_cspline;
}

table3d::~table3d() {
  tree.clear();
  for(size_t i=0;i<list.size();i++) {
    list[i].clear();
  }
  list.clear();
  
  if (xy_set) {
    xval.clear();
    yval.clear();
  }
}

table3d::table3d(table_units<> &t, std::string colx, std::string coly) {

  // Default constructor
  xy_set=false;
  size_set=false;
  has_slice=false;
  itype=itp_cspline;
  
  // Create grid vectors

  std::vector<double>::iterator it;

  std::vector<double> xgrid=t.get_column(colx);
  std::sort(xgrid.begin(),xgrid.end());
  it=std::unique(xgrid.begin(),xgrid.end());
  xgrid.resize(std::distance(xgrid.begin(),it));

  std::vector<double> ygrid=t.get_column(coly);
  std::sort(ygrid.begin(),ygrid.end());
  it=std::unique(ygrid.begin(),ygrid.end());
  ygrid.resize(std::distance(ygrid.begin(),it));

  // Check sizing
  if (t.get_nlines()!=xgrid.size()*ygrid.size()) {
    O2SCL_ERR("Sizes incommensurate.",exc_efailed);
  }

  // Set up grid
  set_xy(colx,xgrid.size(),xgrid,coly,ygrid.size(),ygrid);

  // Find minimum grid spacings
  double min_x=0.0, min_y=0.0;
  for(size_t i=0;i<xgrid.size()-1;i++) {
    if (fabs(xgrid[i+1]-xgrid[i])<min_x) {
      min_x=fabs(xgrid[i+1]-xgrid[i]);
    }
  }
  for(size_t i=0;i<ygrid.size()-1;i++) {
    if (fabs(ygrid[i+1]-ygrid[i])<min_y) {
      min_y=fabs(ygrid[i+1]-ygrid[i]);
    }
  }

  // Create new slices
  for(size_t j=0;j<t.get_ncolumns();j++) {
    if (t.get_column_name(j)!=colx &&
	t.get_column_name(j)!=coly) {
      new_slice(t.get_column_name(j));
    }
  }

  // Loop through lines
  for(size_t i=0;i<t.get_nlines();i++) {
    // Check that we're on the a grid point

    std::vector<double>::iterator itx, ity;
    itx=std::find(xgrid.begin(),xgrid.end(),t.get(colx,i));
    ity=std::find(xgrid.begin(),xgrid.end(),t.get(coly,i));

    double xtarget=*itx;
    double ytarget=*ity;
    if (fabs(t.get(colx,i)-xtarget)>min_x/10.0 ||
	fabs(t.get(coly,i)-ytarget)>min_y/10.0) {
      O2SCL_ERR("Row not on grid point.",exc_efailed);
    }
    // If so, loop through all columns to set slice data
    for(size_t j=0;j<t.get_ncolumns();j++) {
      string col_name=t.get_column_name(j);
      if (col_name!=colx && col_name!=coly) {
	set_val(xtarget,ytarget,col_name,t.get(col_name,i));
      }
    }
  }

  // Copy constants over
  for(size_t i=0;i<t.get_nconsts();i++) {
    string cname;
    double cval;
    t.get_constant(i,cname,cval);
    set_constant(cname,cval);
  }

}

table3d::table3d(const table3d &t) {
      
  // Copy constants
  constants=t.constants;
      
  // Copy interpolation type
  itype=t.itype;
      
  // Copy grid
  numx=t.numx;
  numy=t.numy;
  xname=t.xname;
  yname=t.yname;
  xval=t.xval;
  yval=t.yval;
  xy_set=t.xy_set;
  size_set=t.size_set;
  has_slice=t.has_slice;
      
  for(size_t i=0;i<t.get_nslices();i++) {
	
    // Slice name
    std::string sl_name=t.get_slice_name(i);
	
    new_slice(sl_name);
	
    // Fill the data
    for(size_t j=0;j<t.get_nx();j++) {
      for(size_t k=0;k<t.get_ny();k++) {
	set(j,k,sl_name,t.get(j,k,sl_name));
      }
    }
	
  }
      
  return;
}
    
table3d &table3d::operator=(const table3d &t) {
      
  if (this!=&t) {
	
    clear();
	
    // Copy constants
    constants=t.constants;
	
    // Copy interpolation type
    itype=t.itype;
	
    // Copy grid
    numx=t.numx;
    numy=t.numy;
    xname=t.xname;
    yname=t.yname;
    xval=t.xval;
    yval=t.yval;
    xy_set=t.xy_set;
    size_set=t.size_set;
    has_slice=t.has_slice;
	
    for(size_t i=0;i<t.get_nslices();i++) {
	  
      // Slice name
      std::string sl_name=t.get_slice_name(i);
	  
      new_slice(sl_name);
	  
      // Fill the data
      for(size_t j=0;j<t.get_nx();j++) {
	for(size_t k=0;k<t.get_ny();k++) {
	  set(j,k,sl_name,t.get(j,k,sl_name));
	}
      }
	  
    }
	
  }
      
  return *this;
}

void table3d::set_xy(std::string x_name, uniform_grid<double> gx, 
		     std::string y_name, uniform_grid<double> gy) {

  if (has_slice && (size_set || xy_set) && 
      (gx.get_npoints()!=numx || gy.get_npoints()!=numy)) {
    O2SCL_ERR("Size cannot be reset in table3d::set_xy().",
	      o2scl::exc_einval);
    return;
  }

  if (xy_set) {
    xval.clear();
    yval.clear();
  }
  numx=gx.get_npoints();
  numy=gy.get_npoints();
  xname=x_name;
  yname=y_name;
  xval.resize(numx);
  yval.resize(numy);
  gx.vector(xval);
  gy.vector(yval);
  size_set=true;
  xy_set=true;
}

int table3d::read_gen3_list(std::istream &fin, int verbose, double eps) {
      
  double data;
  std::string line;
  std::string cname, xname_loc="x", yname_loc="y";

  // Read first line and into object called 'onames'
  std::vector<std::string> onames, nnames;
  getline(fin,line);
  std::istringstream is(line);
  while (is >> cname) {
    onames.push_back(cname);
    if (verbose>1) {
      std::cout << "Read possible name: " << cname << std::endl;
    }
  }

  // Return error if there aren't enough columns
  if (onames.size()<3) {
    std::cout << "Not enough columns of data." << std::endl;
    return o2scl::exc_efailed;
  }
      
  // We store a full copy of the entire data set. This is required
  // because we don't know the X-Y grid until the full data file is
  // read. Here, we make space for all of the columns so we can
  // fill them later

  std::vector<std::vector<double> > odata;
  std::vector<double> empty;
  for(size_t i=0;i<onames.size();i++) {
    odata.push_back(empty);
  }
	
  // Count number of likely numbers in the first row
  size_t n_nums=0;
  for(size_t i=0;i<onames.size();i++) {
    if (is_number(onames[i])) n_nums++;
  }

  // Count rows of data
  int irow=0;

  // All of the entries appear to be numbers
  if (n_nums==onames.size()) {
	
    if (verbose>0) {
      std::cout << "First row looks like it contains numerical values." 
		<< std::endl;
      std::cout << "Creating generic slice names: ";
    }

    for(size_t i=2;i<onames.size();i++) {
      nnames.push_back(((std::string)"s")+szttos(i-1));
      if (verbose>0) std::cout << nnames[i-2] << " ";
	  
    }
    if (verbose>0) std::cout << std::endl;
	
    // Add first row of data
    for(size_t i=0;i<onames.size();i++) {
      std::cout << "Adding: " << o2scl::stod(onames[i]) << std::endl;
      odata[i].push_back(o2scl::stod(onames[i]));
    }
    
    irow++;

  } else {

    // Presume the first row contains column names

    // Grid names
    xname_loc=onames[0];
    yname_loc=onames[1];
    if (verbose>0) {
      std::cout << "X grid name: " << onames[0] << endl;
      std::cout << "Y grid name: " << onames[1] << endl;
    }
	
    // Make slices
    for(size_t i=2;i<onames.size();i++) {
      nnames.push_back(onames[i]);
      if (verbose>0) {
	std::cout << "Slice " << i-2 << " name: " << onames[i]
		  << std::endl;
      }
    }
    
  }

  // Read remaining rows
  while ((fin) >> data) {
    if (verbose>2) {
      std::cout << "data: " << 0 << " " << data << std::endl;
    }
    odata[0].push_back(data);
    for(size_t i=1;i<onames.size();i++) {
      (fin) >> data;
      if (verbose>2) {
	std::cout << "data: " << i << " " << data << std::endl;
      }
      odata[i].push_back(data);
    }
    irow++;
  }

  // Setup x and y grid vectors from data
  std::vector<double> xgrid, ygrid;
  for(size_t i=0;i<odata[0].size();i++) {
    bool found=false;
    for(size_t j=0;j<xgrid.size();j++) {
      if (xgrid[j]<eps) {
	if (fabs(odata[0][i]-xgrid[j])<eps) {
	  found=true;
	}
      } else {
	if (fabs(odata[0][i]-xgrid[j])/fabs(xgrid[j])<eps) {
	  found=true;
	}
      }
    }
    if (found==false) {
      if (verbose>1) {
	cout << "Adding " << odata[0][i] << " to xgrid." << endl;
      }
      xgrid.push_back(odata[0][i]);
    }
    found=false;
    for(size_t j=0;j<ygrid.size();j++) {
      if (ygrid[j]<eps) {
	if (fabs(odata[1][i]-ygrid[j])<eps) {
	  found=true;
	}
      } else {
	if (fabs(odata[1][i]-ygrid[j])/fabs(ygrid[j])<eps) {
	  found=true;
	}
      }
    }
    if (found==false) {
      if (verbose>1) {
	cout << "Adding " << odata[1][i] << " to ygrid." << endl;
      }
      ygrid.push_back(odata[1][i]);
    }
  }

  // Sort grids
  vector_sort_double(xgrid.size(),xgrid);
  vector_sort_double(ygrid.size(),ygrid);
  
  if (verbose>1) {
    cout << "x grid (size " << xgrid.size() << "):" << endl;
    for(size_t k=0;k<xgrid.size();k++) {
      std::cout << k << " " << xgrid[k] << std::endl;
    }
    cout << "y grid (size " << ygrid.size() << "):" << endl;
    for(size_t k=0;k<ygrid.size();k++) {
      std::cout << k << " " << ygrid[k] << std::endl;
    }
  }
      
  // Set grid from x and y vectors
  set_xy(xname_loc,xgrid.size(),xgrid,yname_loc,ygrid.size(),ygrid);

  // Create new slices
  for(size_t i=0;i<nnames.size();i++) {
    if (verbose>0) {
      std::cout << "New slice: " << nnames[i] << std::endl;
    }
    new_slice(nnames[i]);
    set_slice_all(nnames[i],0.0);
  }

  // Set the data
  for(size_t j=2;j<odata.size();j++) {
    for(size_t i=0;i<odata[j].size();i++) {
      if (verbose>2) {
	std::cout << "Set value: " << odata[j][i] << std::endl;
      }
      set_val(odata[0][i],odata[1][i],nnames[j-2],odata[j][i]);
    }
  }

  return 0;
}

std::string table3d::get_x_name() const {
  return xname;
}

std::string table3d::get_y_name() const {
  return yname;
}

void table3d::set_x_name(std::string name) {
  xname=name;
  return;
}
    
void table3d::set_y_name(std::string name) {
  yname=name;
  return;
}

const ubvector &table3d::get_x_data() const {
  return xval;
}

const ubvector &table3d::get_y_data() const {
  return yval;
}

size_t table3d::get_nx() const {
  return numx;
}
    
size_t table3d::get_ny() const {
  return numy;
}

bool table3d::is_size_set() const {
  return size_set;
}

bool table3d::is_xy_set() const {
  return xy_set;
}

void table3d::add_slice_from_table(table3d &source, std::string slice,
				   std::string dest_slice, int verbose) {
      
  if (dest_slice.length()==0) dest_slice=slice;

  if (xy_set==false) {
    set_xy(source.get_x_name(),source.get_nx(),source.get_x_data(),
	   source.get_y_name(),source.get_ny(),source.get_y_data());
    new_slice(dest_slice);
    for(size_t i=0;i<numx;i++) {
      for(size_t j=0;j<numx;j++) {
	set(i,j,dest_slice,source.get(i,j,slice));
      }
    }
    return;
  }

  size_t szt_tmp;
  
  if (verbose>2) {
    std::cout << "Creating new slice: " << dest_slice << endl;
  }
  
  if (!is_slice(dest_slice,szt_tmp)) new_slice(dest_slice);
  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      double val=source.interp(get_grid_x(i),get_grid_y(j),slice);
      if (verbose>2) {
        //std::cout << i << " " << j << " " << val << std::endl;
      }
      set(i,j,dest_slice,val);
    }
  }
  return;
}

void table3d::set_size(size_t nx, size_t ny) {
  if ((has_slice && size_set) || xy_set) {
    O2SCL_ERR("Size cannot be reset in table3d::set_xy().",
	      exc_einval);
    return;
  }
  if (xy_set) {
    xval.clear();
    yval.clear();
  }
  numx=nx;
  numy=ny;
  size_set=true;
  return;
}

void table3d::set_slice_all(std::string name, double val) {
  size_t z=lookup_slice(name);
  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      (list[z])(i,j)=val;
    }
  }
  return;
}

void table3d::set(size_t ix, size_t iy, std::string name, double val) {
  size_t z=lookup_slice(name);
  (list[z])(ix,iy)=val;
  return;
}

void table3d::set_val(double x, double y, std::string name, double val) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  
  size_t z=lookup_slice(name);
  (list[z])(ix,iy)=val;
  return;
}
    
void table3d::set_val_ret(double &x, double &y, std::string name, double val) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];

  size_t z=lookup_slice(name);
  (list[z])(ix,iy)=val;
  return;
}
    
void table3d::set(size_t ix, size_t iy, size_t z, double val) {
  (list[z])(ix,iy)=val;
  return;
}

void table3d::set_val(double x, double y, size_t z, double val) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];
  
  (list[z])(ix,iy)=val;
  return;
}
    
void table3d::set_val_ret(double &x, double &y, size_t z, double val) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);

  (list[z])(ix,iy)=val;
  return;
}
    
double &table3d::get(size_t ix, size_t iy, std::string name) {
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}

const double &table3d::get(size_t ix, size_t iy, std::string name) const {
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}

double &table3d::get_val_ret(double &x, double &y, std::string name) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}

const double &table3d::get_val_ret(double &x, double &y, 
				   std::string name) const {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}

double &table3d::get_val(double x, double y, std::string name) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}
    
const double &table3d::get_val(double x, double y, std::string name) const {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  size_t z=lookup_slice(name);
  return (list[z])(ix,iy);
}
    
double &table3d::get(size_t ix, size_t iy, size_t z) {
  return (list[z])(ix,iy);
}

const double &table3d::get(size_t ix, size_t iy, size_t z) const {
  return (list[z])(ix,iy);
}

double &table3d::get_val_ret(double &x, double &y, size_t z) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];
  return (list[z])(ix,iy);
}

const double &table3d::get_val_ret(double &x, double &y, size_t z) const {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  x=xval[ix];
  y=yval[iy];
  return (list[z])(ix,iy);
}

double &table3d::get_val(double x, double y, size_t z) {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  return (list[z])(ix,iy);
}

const double &table3d::get_val(double x, double y, size_t z) const {
  size_t ix=0, iy=0;
  lookup_x(x,ix);
  lookup_y(y,iy);
  return (list[z])(ix,iy);
}

void table3d::set_grid_x(size_t ix, double val) {
  if (ix<numx) {
    (xval)[ix]=val;
    return;
  }
  O2SCL_ERR((((string)"Index '")+itos(ix)+"' out of range ('"+itos(numx)+
	     "') in table3d::set_grid_x().").c_str(),exc_einval);
  return;
}
    
void table3d::set_grid_y(size_t iy, double val) {
  if (iy<numy) {
    (yval)[iy]=val;
    return;
  }
  O2SCL_ERR((((string)"Index '")+itos(iy)+"' out of range ('"+itos(numy)+
	     "') in table3d::set_grid_y().").c_str(),exc_einval);
  return;
}
    
double table3d::get_grid_x(size_t ix) const {
  if (ix<numx) {
    return (xval)[ix];
  }
  O2SCL_ERR((((string)"Index '")+itos(ix)+"' out of range ('"+itos(numx)+
	     "') in table3d::get_grid_x().").c_str(),exc_einval);
  return 0.0;
}
    
double table3d::get_grid_y(size_t iy) const {
  if (iy<numy) {
    return (yval)[iy];
  }
  O2SCL_ERR((((string)"Index '")+itos(iy)+"' out of range ('"+itos(numy)+
	     "') in table3d::get_grid_y().").c_str(),exc_einval);
  return 0.0;
}

void table3d::get_size(size_t &nx, size_t &ny) const {
  if (!size_set) {
    nx=0;
    ny=0;
    return;
  }
  nx=numx;
  ny=numy;
  return;
}

size_t table3d::get_nslices() const {
  return tree.size();
}

void table3d::line_of_names(std::string names) {
  std::string head;
      
  std::istringstream is(names);
      
  while(is >> head) {
    new_slice(head);
  }
      
  return;
}

std::string table3d::get_slice_name(size_t col) const {
  for(map_const_iter mit=tree.begin();mit!=tree.end();mit++) {
    if (mit->second==col) return mit->first;
  }
  O2SCL_ERR("Not found in table3d::get_slice_name().",exc_einval);
  return "";
}

void table3d::new_slice(std::string name) {
  if (size_set==false) {
    O2SCL_ERR2("Size not set before new_slice() ",
	       "in table3d::new_slice().",exc_einval);
  }
  size_t z;
  bool found=is_slice(name,z);
  if (found==true) {
    O2SCL_ERR("Already a slice with the same name in table3d::new_slice().",
	      exc_einval);
  }
  ubmatrix mp(numx,numy);
  list.push_back(mp);
  tree.insert(make_pair(name,list.size()-1));
  has_slice=true;
  return;
}

size_t table3d::lookup_slice(std::string name) const {
  for(map_const_iter mit=const_begin();mit!=const_end();mit++) {
    if (mit->first==name) {
      return mit->second;
    }
  }
  O2SCL_ERR((((string)"Failed to find slice named '")+name+
	     "' in table3d::lookup_slice().").c_str(),exc_enotfound);
  return 0;
}

bool table3d::is_slice(std::string name, size_t &ix) const {
  for(map_const_iter mit=const_begin();mit!=const_end();mit++) {
    if (mit->first==name) {
      ix=mit->second;
      return true;
    }
  }
  return false;
}

void table3d::delete_slice(std::string sl) {
  size_t ix;
  for(map_const_iter mit=const_begin();mit!=const_end();mit++) {
    if (mit->first==sl) {
      ix=mit->second;
      list.erase(list.begin()+ix);
      for(map_iter mit2=tree.begin();mit2!=tree.end();mit2++) {
        if (mit2->second>ix) mit2->second--;
      }
      tree.erase(mit);
      cout << "Here: " << list.size() << " " << tree.size() << endl;
      return;
    }
  }
  O2SCL_ERR((((string)"Failed to find slice named '")+sl+
             "' in table3d::delete_slice().").c_str(),exc_efailed);
  return;
}

void table3d::rename_slice(std::string olds, std::string news) {

  if (news==olds) return;
  
  size_t oi=0;
  if (!is_slice(olds,oi)) {
    O2SCL_ERR((((string)"Failed to find slice named '")+olds+
	       "' in table3d::rename_slice().").c_str(),exc_efailed);
  }

  new_slice(news);
  size_t ni=lookup_slice(news);
  oi=lookup_slice(olds);
  
  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      list[ni](i,j)=list[oi](i,j);
    }
  }

  delete_slice(olds);
  
  return;
}

void table3d::copy_slice(std::string src, std::string dest) {
  size_t sl1=lookup_slice(src);
  
  new_slice(dest);
  size_t sl2=lookup_slice(dest);

  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      (list[sl2])(i,j)=(list[sl1])(i,j);
    }
  }

  return;
}
  
void table3d::init_slice(std::string scol, double val) {
  if (!std::isfinite(val)) {
    O2SCL_ERR("Value not finite in table3d::init_slice()",
	      exc_einval);
  }
  size_t sl1=lookup_slice(scol);
  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      (list[sl1])(i,j)=val;
    }
  }
  return;
}
  
void table3d::lookup_x(double val, size_t &ix) const {
  if (numx>0) {
    double dist=fabs(val-(xval)[0]);
    ix=0;
    for(size_t i=1;i<numx;i++) {
      if (fabs(val-(xval)[i])<dist) {
	ix=i;
	dist=fabs(val-(xval)[i]);
      }
    }
    return;
  }
  O2SCL_ERR("No grid specified in table3d::lookup_x().",exc_einval);
  return;
}
  
void table3d::lookup_y(double val, size_t &iy) const {
  if (numy>0) {
    double dist=fabs(val-(yval)[0]);
    iy=0;
    for(size_t i=1;i<numy;i++) {
      if (fabs(val-(yval)[i])<dist) {
	iy=i;
	dist=fabs(val-(yval)[i]);
      }
    }
    return;
  }
  O2SCL_ERR("No grid specified in table3d::lookup_y().",exc_einval);
  return;
}

void table3d::lookup(double val, std::string slice, size_t &ix, 
		     size_t &iy) const {
  if (numx>0 && numy>0) {
    double dist=fabs(get(0,0,slice)-val);
    ix=0;
    iy=0;
    for(size_t i=0;i<numx;i++) {
      for(size_t j=0;j<numy;j++) {
	if (fabs(get(i,j,slice)-val)<dist) {
	  dist=fabs(get(i,j,slice)-val);
	  ix=i;
	  iy=j;
	}
      }
    }
    return;
  }
  O2SCL_ERR("No grid specified in table3d::lookup().",exc_einval);
  return;
}

void table3d::set_interp_type(size_t interp_type) {
  itype=interp_type;
  return;
}

size_t table3d::get_interp_type() const {
  return itype;
}

double table3d::interp(double x, double y, std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_column> itp;
  
  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    // If we don't need to interpolate, just perform the lookup
    if (numx==1 && x==xval[0]) {
      icol[i]=col[0];
    } else {
      itp.set(numx,xval,col,itype);
      icol[i]=itp.eval(x);
    }
  }

  if (numy==1 && y==yval[0]) {
    result=icol[0];
  } else {
    interp_vec<ubvector> siy(numy,yval,icol,itype);
    result=siy.eval(y);
  }

  return result;
}

void table3d::deriv_y(std::string fname, std::string fpname) {

  size_t z=lookup_slice(fname);
  size_t zp;
  if (!is_slice(fpname,zp)) {
    new_slice(fpname);
    zp=lookup_slice(fpname);
  }
  
  interp_vec<ubvector,ubmatrix_row> itp;

  for(size_t i=0;i<numx;i++) {
    ubmatrix_row row(list[z],i);
    itp.set(numy,yval,row,itype);
    for(size_t j=0;j<numy;j++) {
      set(i,j,zp,itp.deriv(xval[i]));
    }
  }
  
  return;
}

void table3d::deriv_x(std::string fname, std::string fpname) {

  size_t z=lookup_slice(fname);
  size_t zp;
  if (!is_slice(fpname,zp)) {
    new_slice(fpname);
    zp=lookup_slice(fpname);
  }
  
  interp_vec<ubvector,ubmatrix_column> itp;

  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    for(size_t j=0;j<numx;j++) {
      set(j,i,zp,itp.deriv(xval[j]));
    }
  }
  
  return;
}

double table3d::deriv_x(double x, double y, std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_column> itp;

  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    icol[i]=itp.deriv(x);
  }
      
  interp_vec<ubvector> siy(numy,yval,icol,itype);
  result=siy.eval(y);

  return result;
}

double table3d::deriv_y(double x, double y, std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_column> itp;

  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    icol[i]=itp.eval(x);
  }
      
  interp_vec<ubvector> siy(numy,yval,icol,itype);
  result=siy.deriv(y);

  return result;
}

double table3d::integ_x(double x1, double x2, double y,
			std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_row> itp;

  ubvector icol(numx);
  for(size_t i=0;i<numx;i++) {
    ubmatrix_row row(list[z],i);
    itp.set(numy,yval,row,itype);
    icol[i]=itp.eval(y);
  }
      
  interp_vec<ubvector> siy(numx,xval,icol,itype);
  result=siy.integ(x1,x2);

  return result;
}

double table3d::integ_y(double x, double y1, double y2,
			std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_column> itp;

  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    icol[i]=itp.eval(x);
  }
      
  interp_vec<ubvector> siy(numy,yval,icol,itype);
  result=siy.integ(y1,y2);

  return result;
}

double table3d::deriv_xy(double x, double y, std::string name) const {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector,ubmatrix_column> itp;

  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    icol[i]=itp.deriv(x);
  }
      
  interp_vec<ubvector> siy(numy,yval,icol,itype);
  result=siy.deriv(y);

  return result;
}

void table3d::zero_table() {
  for(size_t i=0;i<list.size();i++) {
    for(size_t j=0;j<numx;j++) {
      for(size_t k=0;k<numy;k++) {
	list[i](j,k)=0.0;
      }
    }
  }
  return;
}

void table3d::clear() {
  clear_data();
  xval.clear();
  yval.clear();
  numx=0;
  numy=0;
  has_slice=false;
  xy_set=false;
  size_set=false;
  return;
}

void table3d::clear_data() {
  tree.clear();
  for(size_t i=0;i<list.size();i++) {
    list[i].clear();
  }
  list.clear();
      
  has_slice=false;
  return;
}

void table3d::summary(std::ostream *out, int ncol) const {

  if (out==0) out=&std::cout;
  
  if (constants.size()==1) {
    (*out) << "1 constant:" << std::endl;
  } else {
    (*out) << constants.size() << " constants:" << std::endl;
  }
  std::map<std::string,double>::const_iterator mit;
  for(mit=constants.begin();mit!=constants.end();mit++) {
    (*out) << mit->first << " " << mit->second << std::endl;
  }
      
  // Output number of slices
  int nh=get_nslices();
  if (nh==1) {
    (*out) << "1 slice: " << std::endl;
  } else {
    (*out) << nh << " slices: " << std::endl;
  }
  std::vector<std::string> h(nh);
  for(int i=0;i<nh;i++) {
    h[i]=itos(i)+". "+get_slice_name(i);
  }
      
  // Convert to string with width 'ncol'
  vector<std::string> h2;
  o2scl::screenify(nh,h,h2);
  size_t nh2=h2.size();
      
  // Output slice names
  for(size_t i=0;i<nh2;i++) {
    (*out) << h2[i] << std::endl;
  }
      
  if (numx>0) {
    (*out) << "X-grid " << xname << ", size: " 
	   << numx << " First: " << (xval)[0] 
	   << " Last: " << (xval)[numx-1] << std::endl;
  } else {
    (*out) << "No x-grid." << std::endl;
  }

  if (numy>0) {
    (*out) << "Y-grid " << yname << ", size: " 
	   << numy << " First: " << (yval)[0] 
	   << " Last: " << (yval)[numy-1] << std::endl;
  } else {
    (*out) << "No y-grid." << std::endl;
  }

  return;
}

bool table3d::is_constant(std::string name) const {
  if (constants.find(name)==constants.end()) {
    return false;
  }
  return true;
}


void table3d::get_constant(size_t ix, std::string &name, 
			   double &val) const {
  if (ix<constants.size()) {
    std::map<std::string,double>::const_iterator cit=constants.begin();
    for(size_t i=0;i<ix;i++) cit++;
    name=cit->first;
    val=cit->second;
    return;
  }
  O2SCL_ERR("Index too large in table3d::get_constant().",exc_eindex);
}

void table3d::add_constant(std::string name, double val) {
  if (constants.find(name)!=constants.end()) {
    constants.find(name)->second=val;
    return;
  }
  constants.insert(make_pair(name,val));
  return;
}

void table3d::remove_constant(std::string name) {
  constants.erase(name);
  return;
}

int table3d::set_constant(std::string name, double val,
			  bool err_on_notfound) {
  if (constants.find(name)!=constants.end()) {
    constants.find(name)->second=val;
    return 0;
  }
  if (err_on_notfound) {
    std::string err=((std::string)"No constant with name '")+name+
      "' in table3d::set_constant().";
    O2SCL_ERR(err.c_str(),exc_enotfound);
  }
  return exc_enotfound;
}

double table3d::get_constant(std::string name) {
  if (constants.find(name)==constants.end()) {
    std::string err=((std::string)"No constant with name '")+name+
      "' in table3d::get_constant().";
    O2SCL_ERR(err.c_str(),exc_einval);
  }
  return constants.find(name)->second;
}

void table3d::extract_x(double x, table<> &t) {

  if (get_ny()==0) {
    O2SCL_ERR2("Cannot extract slice from table with zero y-grid size ",
	       "in table3d::extract_x().",exc_einval);
  }
  if (get_nx()<2) {
    O2SCL_ERR2("Cannot extract slice from table with x-grid size<2 ",
	       "in table3d::extract_x().",exc_einval);
  }
  
  t.clear();
  
  string s=yname+" ";
  for(size_t i=0;i<tree.size();i++) {
    s+=get_slice_name(i)+" ";
  }
  t.line_of_names(s);
  t.set_nlines(numy);

  for(size_t i=0;i<numy;i++) {
    t.set(yname,i,yval[i]);
    for(size_t j=0;j<tree.size();j++) {
      string tmp=get_slice_name(j);
      double val=interp(x,yval[i],tmp);
      t.set(tmp,i,val);
    }
  }
  
  return;
}

void table3d::extract_y(double y, table<> &t) {
  
  if (get_nx()==0) {
    O2SCL_ERR2("Cannot extract slice from table with zero x-grid size ",
	       "in table3d::extract_x().",exc_einval);
  }
  if (get_ny()<2) {
    O2SCL_ERR2("Cannot extract slice from table with y-grid size<2 ",
	       "in table3d::extract_x().",exc_einval);
  }
  
  t.clear_table();
  
  string s=xname+" ";
  for(size_t i=0;i<tree.size();i++) {
    s+=get_slice_name(i)+" ";
  }
  t.line_of_names(s);
  t.set_nlines(numx);

  for(size_t i=0;i<numx;i++) {
    t.set(xname,i,xval[i]);
    for(size_t j=0;j<tree.size();j++) {
      string tmp=get_slice_name(j);
      t.set(tmp,i,interp(xval[i],y,tmp));
    }
  }
  
  return;
}
   
const boost::numeric::ublas::matrix<double> &table3d::get_slice
(std::string name) const {
  size_t z=lookup_slice(name);
  return list[z];
}

const boost::numeric::ublas::matrix<double> &table3d::get_slice
(size_t iz) const {
  return list[iz];
}

boost::numeric::ublas::matrix<double> &table3d::get_slice
(std::string name) {
  size_t z=lookup_slice(name);
  return list[z];
}

boost::numeric::ublas::matrix<double> &table3d::get_slice(size_t iz) {
  return list[iz];
}

const vector<boost::numeric::ublas::matrix<double> > &table3d::get_data() {
  return list;
}

void table3d::function_slice(string function, string scol) {

  size_t ic;
  if (is_slice(scol,ic)==false) {
    new_slice(scol);
    ic=lookup_slice(scol);
  }

  function_matrix(function,list[ic]);

  return;
}

table3d table3d::slice_to_uniform_grid(std::string slice, size_t xpts,
				       bool log_x, size_t ypts,
				       bool log_y) {
  uniform_grid<double> ugx, ugy;
  if (log_x) {
    ugx=uniform_grid_log_end<double>(xval[0],xval[numx-1],xpts-1);
  } else {
    ugx=uniform_grid_end<double>(xval[0],xval[numx-1],xpts-1);
  }
  if (log_y) {
    ugy=uniform_grid_log_end<double>(yval[0],yval[numy-1],ypts-1);
  } else {
    ugy=uniform_grid_end<double>(yval[0],yval[numy-1],ypts-1);
  }
  table3d t3d;
  t3d.set_xy(xname,ugx,yname,ugy);
  t3d.new_slice(slice);
  for(size_t i=0;i<t3d.get_nx();i++) {
    for(size_t j=0;j<t3d.get_ny();j++) {
      t3d.set(i,j,slice,this->interp(t3d.get_grid_x(i),
				     t3d.get_grid_y(j),slice));
    }
  }
  return t3d;
}

table3d table3d::table_to_uniform_grid(size_t xpts, bool log_x, 
				       size_t ypts, bool log_y) {

  uniform_grid<double> ugx, ugy;
  if (log_x) {
    ugx=uniform_grid_log_end<double>(xval[0],xval[numx-1],xpts-1);
  } else {
    ugx=uniform_grid_end<double>(xval[0],xval[numx-1],xpts-1);
  }
  if (log_y) {
    ugy=uniform_grid_log_end<double>(yval[0],yval[numy-1],ypts-1);
  } else {
    ugy=uniform_grid_end<double>(yval[0],yval[numy-1],ypts-1);
  }

  table3d t3d;
  t3d.set_xy(xname,ugx,yname,ugy);
  for(size_t k=0;k<this->get_nslices();k++) {
    std::string sl_name=this->get_slice_name(k);
    t3d.new_slice(sl_name);
    for(size_t i=0;i<t3d.get_nx();i++) {
      for(size_t j=0;j<t3d.get_ny();j++) {
	t3d.set(i,j,sl_name,this->interp(t3d.get_grid_x(i),
					 t3d.get_grid_y(j),sl_name));
      }
    }
  }
  return t3d;
}
