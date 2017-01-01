/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2017, Andrew W. Steiner
  
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
/* With gcc: add an underscore to the end of the function name and
   convert all arguments into pointers.  The Fortran routines use the
   transpose of matrices as they are stored in C, so you can either
   transpose the matrix before input or remember that the output is
   transposed.  I use the liblapack.a library from the Debian lapack-dev
   package, so I did not have to compile lapack myself, although I could
   have done that with g77.  
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

extern "C" {

  void zgesvd_(char *jobu, char *jobvt, int *m, int *n, 
		     double *a, int *lda, double *s, 
		     double *u, int *ldu, double *vt, 
		     int *ldvt, double *work, 
		     int *lwork, double *rwork, int *info);
}

int main () {

  size_t M = 5, N = 5;
  size_t L = GSL_MIN (M, N);

  gsl_matrix_complex *m = gsl_matrix_complex_alloc (M, N);
  gsl_vector *s = gsl_vector_alloc (L);
  gsl_matrix_complex *u = gsl_matrix_complex_alloc (M, N);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (M, N);

  int lwork = 2 * L + GSL_MAX (M, N);
  gsl_vector_complex *work = gsl_vector_complex_alloc (lwork);
  gsl_vector_complex *rwork = gsl_vector_complex_alloc (5 * L);

  int status;
  size_t i, j;

  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++)
      {
	gsl_complex z = gsl_complex_rect (1.0 / (i + j + 1.0), 0.0);
	gsl_matrix_complex_set (m, i, j, z);
      }

  zgesvd_ ("A", "A",		
	  (int *) &(m->size1), (int *) &(m->size2),
	  m->data, (int *) &(m->tda),
	  s->data,
	  u->data, (int *) &(u->tda),
	  v->data, (int *) &(v->tda),
	  work->data, &lwork, rwork->data, &status);

  gsl_matrix_complex_fprintf (stdout, m, "%g");
  printf ("\n");
  gsl_matrix_complex_fprintf (stdout, u, "%g");
  printf ("\n");
  gsl_matrix_complex_fprintf (stdout, v, "%g");
  printf ("\n");

  gsl_vector_fprintf (stdout, s, "%g");
  printf ("\n");

  return 0;
}
