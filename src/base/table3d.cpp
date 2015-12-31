/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

#include <o2scl/table3d.h>

using namespace std;
using namespace o2scl;

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

int table3d::read_gen3_list(std::istream &fin, int verbose) {
      
  double data;
  std::string line;
  std::string cname, xname="x", yname="y";
      
  std::vector<std::vector<double> > odata;
      
  // Read first line and into list
  std::vector<std::string> onames, nnames;
  getline(fin,line);
  std::istringstream is(line);
  while (is >> cname) {
    onames.push_back(cname);
    if (verbose>2) {
      std::cout << "Read possible name: " << cname << std::endl;
    }
  }

  if (onames.size()<3) {
    std::cout << "Not enough columns of data." << std::endl;
    return o2scl::exc_efailed;
  }
      
  // Create odata vectors
  std::vector<double> empty;
  for(size_t i=0;i<onames.size();i++) {
    odata.push_back(empty);
  }
	
  // Count number of likely numbers in the first row
  size_t n_nums=0;
  for(size_t i=0;i<onames.size();i++) {
    if (is_number(onames[i])) n_nums++;
  }
      
  int irow=0;
      
  if (n_nums==onames.size()) {
	
    if (verbose>0) {
      std::cout << "First row looks like it contains numerical values." 
		<< std::endl;
      std::cout << "Creating generic slice names: ";
    }

    for(size_t i=2;i<onames.size();i++) {
      std::cout << "Here: " << onames.size() << std::endl;
      std::cout << "Here2: " << nnames.size() << std::endl;
      nnames.push_back(((std::string)"s")+szttos(i-1));
      std::cout << "Here: " << onames.size() << std::endl;
      std::cout << "Here2: " << nnames.size() << std::endl;
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

    // Ensure good names
    for(size_t i=0;i<onames.size();i++) {
      std::string temps=onames[i];
      //make_fp_varname(onames[i]);
      //make_unique_name(onames[i],onames);
      if (temps!=onames[i] && verbose>0) {
	std::cout << "Converted slice named '" << onames[i] << "' to '" 
		  << temps << "'." << std::endl;
      }
    }

    // Grid names
    xname=onames[0];
    yname=onames[1];
	
    // Make slices
    for(size_t i=2;i<onames.size();i++) {
      nnames.push_back(onames[i]);
    }
	
  }
      
  // Read remaining rows
  while ((fin) >> data) {
    if (verbose>1) {
      std::cout << "data: " << 0 << " " << data << std::endl;
    }
    odata[0].push_back(data);
    for(size_t i=1;i<onames.size();i++) {
      (fin) >> data;
      if (verbose>1) {
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
      if (fabs(odata[0][i]-xgrid[j])/fabs(xgrid[j])<1.0e-12) {
	found=true;
      }
    }
    if (found==false) {
      xgrid.push_back(odata[0][i]);
    }
    found=false;
    for(size_t j=0;j<ygrid.size();j++) {
      if (fabs(odata[1][i]-ygrid[j])/fabs(ygrid[j])<1.0e-12) {
	found=true;
      }
    }
    if (found==false) {
      ygrid.push_back(odata[1][i]);
    }
  }

  if (verbose>1) {
    for(size_t k=0;k<xgrid.size();k++) {
      std::cout << k << " " << xgrid[k] << std::endl;
    }
    for(size_t k=0;k<ygrid.size();k++) {
      std::cout << k << " " << ygrid[k] << std::endl;
    }
  }
      
  // Set grid from x and y vectors
  set_xy(xname,xgrid.size(),xgrid,yname,ygrid.size(),ygrid);

  // Create new slices
  for(size_t i=0;i<nnames.size();i++) {
    if (verbose>0) {
      std::cout << "New slice: " << nnames[i] << std::endl;
    }
    new_slice(nnames[i]);
  }

  // Set the data
  for(size_t j=2;j<odata.size();j++) {
    for(size_t i=0;i<odata[j].size();i++) {
      if (verbose>1) {
	std::cout << "Set value: " << odata[j][i] << std::endl;
      }
      set_val(odata[0][i],odata[1][i],nnames[j-2],odata[j][i]);
    }
  }

  return 0;
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
    
double table3d::get_grid_x(size_t ix) {
  if (ix<numx) {
    return (xval)[ix];
  }
  O2SCL_ERR((((string)"Index '")+itos(ix)+"' out of range ('"+itos(numx)+
	     "') in table3d::get_grid_x().").c_str(),exc_einval);
  return 0.0;
}
    
double table3d::get_grid_y(size_t iy) {
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

void table3d::rename_slice(std::string olds, std::string news) {
  size_t oi=0;
  if (!is_slice(olds,oi)) {
    O2SCL_ERR((((string)"Failed to find slice named '")+olds+
	       "' in table3d::rename_slice().").c_str(),exc_efailed);
  }

  new_slice(news);
  size_t ni=lookup_slice(news);
  
  for(size_t i=0;i<numx;i++) {
    for(size_t j=0;j<numy;j++) {
      list[ni](i,j)=list[oi](i,j);
    }
  }
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

double table3d::interp(double x, double y, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

  ubvector icol(numy);
  for(size_t i=0;i<numy;i++) {
    ubmatrix_column col(list[z],i);
    itp.set(numx,xval,col,itype);
    icol[i]=itp.eval(x);
  }
      
  interp_vec<ubvector> siy(numy,yval,icol,itype);
  result=siy.eval(y);

  return result;
}

double table3d::deriv_x(double x, double y, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

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

double table3d::deriv_y(double x, double y, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

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

double table3d::integ_x(double x1, double x2, double y, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

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

double table3d::integ_y(double x, double y1, double y2, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

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

double table3d::deriv_xy(double x, double y, std::string name) {
  double result;
  
  size_t z=lookup_slice(name);
  
  interp_vec<ubvector> itp;

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

void table3d::clear_table() {
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

void table3d::set_constant(std::string name, double val) {
  if (constants.find(name)!=constants.end()) {
    constants.find(name)->second=val;
    return;
  }
  O2SCL_ERR("No constant with specified name in table3d::set_constant().",
	    exc_einval);
}

double table3d::get_constant(std::string name) {
  if (constants.find(name)!=constants.end()) {
    return constants.find(name)->second;
  }
  O2SCL_ERR("No constant with specified name in table3d::get_constant().",
	    exc_einval);
  return 0.0;
}

void table3d::extract_x(double x, table<> &t) {
  t.clear_table();
  
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
      t.set(tmp,i,interp(x,yval[i],tmp));
    }
  }
  
  return;
}

void table3d::extract_y(double y, table<> &t) {
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
      t.set(tmp,i,interp(y,xval[i],tmp));
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

void table3d::function_slice(string function, string scol) {
  
  new_slice(scol);
  size_t ic=lookup_slice(scol);

  function_matrix(function,list[ic]);

  return;
}
