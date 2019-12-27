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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/interp2_seq.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

using namespace std;
using namespace o2scl;

interp2_seq::interp2_seq() {
}

interp2_seq::~interp2_seq() {
}

void interp2_seq::set_data(size_t n_x, size_t n_y, ubvector &x_grid,
			 ubvector &y_grid, ubmatrix &data, 
			 size_t interp_type, bool x_first) {
  // Set new data
  itype=interp_type;
  xfirst=x_first;
  nx=n_x;
  ny=n_y;
  xfun=&x_grid;
  yfun=&y_grid;
  datap=&data;

  // Set interpolation objects

  for(size_t i=0;i<itps.size();i++) {
    delete itps[i];
    delete vecs[i];
  }
  itps.clear();

  //if (x_first==true) {
  // If we interpolate along the x-axis first, then we want to fix the
  // first index, to get nx rows of size ny
  vecs.resize(nx);
  itps.resize(nx);
  for(size_t i=0;i<nx;i++) {
    vecs[i]=new ubmatrix_row
      (o2scl::matrix_row<ubmatrix,ubmatrix_row>(*datap,i));
    itps[i]=new interp_vec<ubvector,ubmatrix_row>(ny,*yfun,*vecs[i],itype);
  }
  /*
    } else {
    // If we interpolate along the y-axis first, then we want to fix the
    // first index, to get ny columns of size nx
    vecs.resize(ny);
    itps.resize(ny);
    for(size_t i=0;i<ny;i++) {
    vecs[i]=ubmatrix_column(*datap,i);
    itps[i]=new interp_vec<ubvector>(nx,*xfun,vecs[i],itype);
    }
    }
  */

  data_set=true;

  return;
}

void interp2_seq::reset_interp() {

  if (data_set) {

    for(size_t i=0;i<itps.size();i++) {
      delete itps[i];
      delete vecs[i];
    }
    itps.clear();

    // Set interpolation objects
    //if (xfirst==true) {
    vecs.resize(nx);
    itps.resize(nx);
    for(size_t i=0;i<nx;i++) {
      vecs[i]=new ubmatrix_row
	(o2scl::matrix_row<ubmatrix,ubmatrix_row>(*datap,i));
      itps[i]=new interp_vec<ubvector,ubmatrix_row>(ny,*yfun,*vecs[i],itype);
    }
    /*
      } else {
      vecs.resize(ny);
      itps.resize(ny);
      for(size_t i=0;i<ny;i++) {
      vecs[i]=ubmatrix_column(*datap,i);
      itps[i]=new interp_vec<ubvector>(nx,*xfun,vecs[i],itype);
      }
      }
    */
    
  } else {
    O2SCL_ERR("Data not set in interp2_seq::reset_interp().",exc_einval);
  }

  return;
}

double interp2_seq::operator()(double x, double y) const {
  return eval(x,y);
}

double interp2_seq::eval(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->eval(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.eval(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->eval(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.eval(y);
  }
  return result;
}

double interp2_seq::deriv_x(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->eval(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.deriv(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->deriv(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.eval(y);
  }
  return result;
}

double interp2_seq::deriv_xy(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::eval().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->deriv(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.deriv(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->deriv(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.deriv(y);
  }
  return result;
}

double interp2_seq::eval_gen(int ix, int iy, double x0, double x1,
			       double y0, double y1) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::interp_gen().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      if (iy==-1) {
	icol[i]=itps[i]->integ(y0,y1);
      } else if (iy==0) {
	icol[i]=itps[i]->eval(y0);
      } else if (iy==1) {
	icol[i]=itps[i]->deriv(y0);
      } else if (iy==2) {
	icol[i]=itps[i]->deriv2(y0);
      } else {
	O2SCL_ERR2("Invalid value of 'iy' for interp2_seq::",
		   "interp_gen(). (xfirst=true)",exc_einval);
      }
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    if (ix==-1) {
      result=six.integ(x0,x1);
    } else if (ix==0) {
      result=six.eval(x0);
    } else if (ix==1) {
      result=six.deriv(x0);
    } else if (ix==2) {
      result=six.deriv2(x0);
    } else {
      O2SCL_ERR2("Invalid value of 'ix' for interp2_seq::",
		 "interp_gen(). (xfirst=true)",exc_einval);
    }
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      if (ix==-1) {
	icol[i]=itps[i]->integ(x0,x1);
      } else if (ix==0) {
	icol[i]=itps[i]->eval(x0);
      } else if (ix==1) {
	icol[i]=itps[i]->deriv(x0);
      } else if (ix==2) {
	icol[i]=itps[i]->deriv2(x0);
      } else {
	O2SCL_ERR2("Invalid value of 'ix' for interp2_seq::",
		   "interp_gen(). (xfirst=false)",exc_einval);
      }
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    if (iy==-1) {
      result=siy.integ(y0,y1);
    } else if (iy==0) {
      result=siy.eval(y0);
    } else if (iy==1) {
      result=siy.deriv(y0);
    } else if (iy==2) {
      result=siy.deriv2(y0);
    } else {
      O2SCL_ERR2("Invalid value of 'iy' for interp2_seq::",
		 "interp_gen(). (xfirst=false)",exc_einval);
    }
  }
  return result;
}

double interp2_seq::deriv_y(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::deriv_y().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->deriv(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.eval(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->eval(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.deriv(y);
  }
  return result;
}

double interp2_seq::deriv_xx(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::deriv_xx().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->eval(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.deriv2(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->deriv2(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.eval(y);
  }
  return result;
}

double interp2_seq::deriv_yy(double x, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::deriv_yy().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->deriv2(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.eval(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->eval(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.deriv2(y);
  }
  return result;
}

double interp2_seq::integ_x(double x0, double x1, double y) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::integ_x().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->eval(y);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.integ(x0,x1);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->integ(x0,x1);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.eval(y);
  }
  return result;
}

double interp2_seq::integ_y(double x, double y0, double y1) const {
  if (data_set==false) {
    O2SCL_ERR("Data not set in interp2_seq::integ_y().",exc_efailed);
    return 0.0;
  }
  double result;
  if (xfirst==true) {
    ubvector icol(nx);
    for(size_t i=0;i<nx;i++) {
      icol[i]=itps[i]->integ(y0,y1);
    }
    interp_vec<ubvector> six(nx,*xfun,icol,itype);
    result=six.eval(x);
  } else {
    ubvector icol(ny);
    for(size_t i=0;i<ny;i++) {
      icol[i]=itps[i]->eval(x);
    }
    interp_vec<ubvector> siy(ny,*yfun,icol,itype);
    result=siy.integ(y0,y1);
  }
  return result;
}
