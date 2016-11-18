/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016, Andrew W. Steiner
  
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
#ifndef O2SCL_PYTHON_H
#define O2SCL_PYTHON_H
/** \file python.h
    \brief Desc
*/

#ifdef O2SCL_PYTHON
  
#include <Python.h>

/** \brief Function to be called from Python
 */
extern "C" PyObject* py_myFunction(PyObject* self, PyObject* args);

/** \brief Another function to be called from Python
 */
static PyObject* py_myOtherFunction(PyObject* self, PyObject* args);
  
/*
  Bind Python function names to our C functions
*/
/*
  static PyMethodDef myModule_methods[] = {
  {"myFunction", py_myFunction, METH_VARARGS},
  {"myOtherFunction", py_myOtherFunction, METH_VARARGS},
  {NULL, NULL}
  };
*/

/** \brief Python calls this to let us initialize our module
 */
void initmyModule();

#endif
  
#endif

