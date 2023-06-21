 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2023, Julien Garaud and Andrew W. Steiner
  
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
#include <o2scl/ode_bv_multishoot.h>
#include <o2scl/ode_funct.h>

#include <cmath>

using namespace std;
using namespace o2scl;

/* Example: ex_hydrogen.cpp
   -------------------------------------------------------------------
   Compute the first bound state of the hydrogen atom. 

   * This example is still in progress *
*/

class Boundary_cl {

public:

  ovector_int index;
  
  static const int Left   = 0;
  static const int Right  = 1;
  static const int Derivs = 2;
  static const int Extra  = 3;

  int n_left;
  int n_right;
  int n_derivs;
  int n_extra;

  int info(const int &where, const ovector_base &y_in, ovector_base &v);

};

Boundary_cl Boundary;

int Boundary_cl::info(const int &where, const ovector_base &y_in,
		      ovector_base &v) {
  /*
    This function sends the right parameters at the right routines
  */
  switch(where){
  case 0:
    if(int(v.size()) != n_left) {
      O2SCL_ERR("Incorrect size of v in Boundary_cl::info()",gsl_einval);
      return gsl_einval;}
    
    break;
  case 1:
    if(int(v.size()) != n_right) {
      O2SCL_ERR("Incorrect size of v in Boundary_cl::info()",gsl_einval);
      return gsl_einval;}
    v[0]=y_in[0]; //beta
    v[1]=y_in[1]; // A
    break;
  case 2:
    if(int(v.size()) != n_derivs) {
      O2SCL_ERR("Incorrect size of v in Boundary_cl::info()",gsl_einval);
      return gsl_einval;}
    v[0]=y_in[0]; //beta
    break;
  case 3:
    if(int(v.size()) != n_extra) {
      O2SCL_ERR("Incorrect size of v in Boundary_cl::info()",gsl_einval);
      return gsl_einval;}
    v[0]=y_in[2];
    v[1]=y_in[3];
    break;
  default : 
    O2SCL_ERR("Incorrect location in Boundary_cl::info()",gsl_einval);
    return gsl_einval;      
  }
  return 0;
}

/*  
    We want that derivs alway know the values of parameters in order
    to pass the eigenvalue beta. So we define a constant pointer to
    the vector of variables. this vector is passed to the function
    info wich returns the vector v and then we have acces to the
    eigenvalue... (or to some other global parameters)

    The radial equation is
    

*/
struct Deriv_class {

  const ovector *y_ptr;

  Deriv_class(const ovector *yy_ptr):y_ptr(yy_ptr) {}

  int derivs(double x, size_t nv, const ovector_base &y, ovector_base &dydx) {
    
    const int where = Boundary.Derivs;
    ovector v(Boundary.n_derivs);
    Boundary.info(where,y_ptr[0],v);
    
    int l=0;//int(param[1]);
    double beta=v[0];
    dydx[0]=y[1];
    dydx[1]= (-1/x+l*(l+1)/(x*x)+beta*beta)*y[0];
     
    return 0;
  }
};

int factorial(int number) {
  int temp;
  
  if(number <= 1) return 1;
  
  temp = number * factorial(number - 1);
  return temp;
}

int Left_boundary(double x, size_t nv, const ovector_base &y, 
		  ovector_base &y_out) {

  const int where = Boundary.Left;
  ovector v(Boundary.n_left);
  Boundary.info(where,y,v);


  y_out[0] = x;
  // FIXME
  y_out[1] = 0.0;

  /* Here the wave function is not normalized */
  return 0;
}

int Right_boundary(double x, size_t nv, const ovector_base &y, 
		   ovector_base &y_out) {

  const int where = Boundary.Right;
  ovector v(Boundary.n_right);
  Boundary.info(where,y,v);
  double beta=v[0];
  double A=v[1];

  /* Exponentially decreasing massive modes */
  y_out[0] = A*exp(-beta*x);
  y_out[1] = -beta*y_out[0];
  return 0;
}

int Extra_boundary(double x, size_t nv, const ovector_base &y, 
		   ovector_base &y_out) {

  const int where = Boundary.Extra;
  ovector v(Boundary.n_extra);
  Boundary.info(where,y,v);

  y_out[0] = v[0];
  y_out[1] = v[1]; 

  /*
    These parameters are useless and set to zero. This conditions
    represents the fact that the link between functions and
    derivatives is already done in the functions Left and Right.
  */
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);
  
  /* Mesh */
  const int N_mesh = 10;
  double X_min     = 0.00001;
  double X_max     = 30.0;
  double Step      = (X_max-X_min)/(N_mesh-1);
  ovector Mesh(N_mesh);
  for(int i=0;i<N_mesh;i++) {
    Mesh[i] = X_min+i*Step;
  }

  /* Starting vector of parameters */
  int N_equations = 2;
  const int N_var = N_equations*(N_mesh);
  ovector Y_start(N_var);

  /*
    The first N_equations points stored in Y_start
    are the shooting parameters of the problem and the next
    N_equations ones are dedicated to additionnal constraints.

    In this implementation, The value of the functions in term of the
    shooting parameters are not explicitly computed so the extra
    boundary conditions are unused. The extra boundary conditions are
    useful for problems with additional constraints

    Another strategy would be not to implement Left and Right and to
    keep all the information in Extra. That is to say give all
    boundaries in the same way in ode_bv_multishoot and they would be
    constrained by the N_equations constraints in the function Extra.

    The latter strategy is more natural in the point of view of the
    implementation. It will certainly be more efficient to deal with
    easy problems.

    Anyway, when considering more involved situations I think the
    first strategy is better. Because of the great flexibility you can
    implement any type of bvp with and without eigenvalues (for
    example). The only thing you have to know for this is that the
    first 2*N_equations parameters are different from the others.
  */

  /* initial values of parameters*/
  Y_start[0] = 0.5; 
  Y_start[1] = 0.1; 
  Y_start[2] = 0.0; 
  Y_start[3] = 0.0; 
 
  /*
    In order to converge, you have to provide an initial guess for the
    functions.
  */
  for(int i=2*N_equations;i<N_var;i+=2) {
    Y_start[i] = 0.01*i; 
    Y_start[i+1] = 0.01*i; 
  }

  /* Free Parameters */
  const int N_param = 2;
  ovector Param(N_param);
  // Here is a global constant parameter (n)
  Param[0] = 1.0; 
  // Here is a global constant parameter (l)
  Param[1] = 0.0; 
  
  /* Indexing */
  // number of shooting parameters living on the leftmost point
  const int n_left   = 0; 
  // number of shooting parameters living on the rightmost point
  const int n_right  = 2; 
  // number of shooting parameters living in the whole interval
  const int n_derivs = 1; 
  // number of shooting parameters living in global constraints
  const int n_extra  = 2; 

  Boundary.n_left   = n_left;
  Boundary.n_right  = n_right;
  Boundary.n_derivs = n_derivs;
  Boundary.n_extra  = n_extra;

  /* Now create the objects for multishooting */
 
  Deriv_class deriv_class_ptr(&Y_start);
  ode_funct_mfptr<Deriv_class> Derivs(&deriv_class_ptr,&Deriv_class::derivs);

  //  ode_funct_fptr<> Derivs(derivs);
  ode_funct_fptr<> Left(Left_boundary);
  ode_funct_fptr<> Right(Right_boundary);
  ode_funct_fptr<> Extra(Extra_boundary);

  /* vector and matrix to save the global profile at the end of integration */
  int n_grid=80;
  ovector x_save(n_grid);
  omatrix y_final(n_grid,N_equations);

  /* Multishooter */
  ode_bv_multishoot<> Multi;
  Multi.solve(Mesh,N_equations,Y_start,Left,Right,Extra,
	      Derivs,x_save,y_final);


  /* Output */
  for(int i=0;i<n_grid;i++) {
    if(x_save[i]!=0.0 || i==0) {
      cout << x_save[i] <<  "\t"
	   << y_final[i][0] << "\t"
	   << y_final[i][1] << endl; 
    }
  }
  cout << "beta = " << Y_start[0] << endl;
  cout << "A = " << Y_start[1] << endl;
  
  return 0;
}

