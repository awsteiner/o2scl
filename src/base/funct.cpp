/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

#include <o2scl/set_python.h>

#ifdef O2SCL_SET_PYTHON
#include <Python.h>
#endif

#include <o2scl/funct.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_SET_PYTHON

funct_python::funct_python(std::string module, std::string func,
                           int v) {

  verbose=v;
  
  if (o2scl_settings.py_initialized==false) {
    if (verbose>0) {
      cout << "Running py_init()." << endl;
    }
    o2scl_settings.py_init();
  }
  set_function(module,func);
}

funct_python::~funct_python() {
  if (verbose>0) {
    cout << "Decref func." << endl;
  }
  Py_DECREF(pFunc);
  if (verbose>0) {
    cout << "Decref args." << endl;
  }
  Py_DECREF(pArgs);
  if (verbose>0) {
    cout << "Decref module." << endl;
  }
  Py_DECREF(pModule);
  if (verbose>0) {
    cout << "Done in funct_python destructor." << endl;
  }
}

int funct_python::set_function(std::string module, std::string func) {

  pModule=o2scl_settings.py_import_module(module,verbose);
  
  // Setup the arguments to the python function
  if (verbose>0) {
    cout << "Getting arguments for python function." << endl;
  }
  pArgs=PyTuple_New(1);
  if (pArgs==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "funct_python::set_function().",o2scl::exc_efailed);
  }

  // Load the python function
  if (verbose>0) {
    cout << "Loading python function." << endl;
  }
  pFunc=PyObject_GetAttrString(pModule,func.c_str());
  if (pFunc==0) {
    O2SCL_ERR2("Get function failed in ",
               "funct_python::set_function().",o2scl::exc_efailed);
  }
  
  return 0;
}

double funct_python::operator()(double x) const {

  // Create a python object from the value
  if (verbose>0) {
    cout << "Creating python object from double." << endl;
  }
  PyObject *pValue=PyFloat_FromDouble(x);
  if (pValue==0) {
    O2SCL_ERR2("Value creation failed in ",
               "funct_python::operator().",o2scl::exc_efailed);
  }

  // Set the python function arguments
  if (verbose>0) {
    cout << "Setting item in tuple." << endl;
  }
  int ret=PyTuple_SetItem(pArgs,0,pValue);
  if (ret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "funct_python::operator().",o2scl::exc_efailed);
  }

  // Call the python function
  if (verbose>0) {
    cout << "Call python function." << endl;
  }
  PyObject *result=PyObject_CallObject(pFunc,pArgs);
  if (result==0) {
    O2SCL_ERR2("Function call failed in ",
               "funct_python::operator().",o2scl::exc_efailed);
  }

  double y=PyFloat_AsDouble(result);

  if (verbose>0) {
    cout << "Decref value and result." << endl;
  }
  Py_DECREF(pValue);
  Py_DECREF(result);
  
  if (verbose>0) {
    cout << "Done in funct_python::operator()." << endl;
  }
    
  return y;
}

funct_python_method::funct_python_method(std::string module,
                                         std::string class_name,
                                         std::string func, int v) {

  verbose=v;
  
  if (o2scl_settings.py_initialized==false) {
    if (verbose>0) {
      cout << "Running py_init()." << endl;
    }
    o2scl_settings.py_init();
  }
  set_function(module,class_name,func);
}

funct_python_method::~funct_python_method() {
  if (verbose>0) {
    cout << "Decref func." << endl;
  }
  Py_DECREF(pFunc);
  if (verbose>0) {
    cout << "Decref instance." << endl;
  }
  Py_DECREF(pInstance);
  if (verbose>0) {
    cout << "Decref class." << endl;
  }
  Py_DECREF(pClass);
  if (verbose>0) {
    cout << "Decref args." << endl;
  }
  Py_DECREF(pArgs);
  if (verbose>0) {
    cout << "Decref module." << endl;
  }
  Py_DECREF(pModule);
  if (verbose>0) {
    cout << "Done in funct_python_method destructor." << endl;
  }
}

int funct_python_method::set_function(std::string module,
                                      std::string class_name,
                                      std::string func) {
    
  pModule=o2scl_settings.py_import_module(module,verbose);
  
  // Setup the arguments to the python function
  if (verbose>0) {
    cout << "Getting arguments for python function." << endl;
  }
  pArgs=PyTuple_New(1);
  if (pArgs==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "funct_python_method::set_function().",o2scl::exc_efailed);
  }

  // Load the python class
  if (verbose>0) {
    cout << "Loading python class." << endl;
  }
  pClass=PyObject_GetAttrString(pModule,class_name.c_str());
  if (pClass==0) {
    O2SCL_ERR2("Get class failed in ",
               "funct_python_method::set_function().",o2scl::exc_efailed);
  }

  // Create an instance of the class
  if (verbose>0) {
    cout << "Loading python class instance." << endl;
  }
  if (PyCallable_Check(pClass)==false) {
    O2SCL_ERR2("Check class callable failed in ",
               "funct_python_method::set_function().",o2scl::exc_efailed);
  }
  
  pInstance=PyObject_CallObject(pClass,0);
  if (pInstance==0) {
    O2SCL_ERR2("Instantiate class failed in ",
               "funct_python_method::set_function().",o2scl::exc_efailed);
  }
  
  // Load the python function
  if (verbose>0) {
    cout << "Loading python function." << endl;
  }
  pFunc=PyObject_GetAttrString(pInstance,func.c_str());
  if (pFunc==0) {
    O2SCL_ERR2("Get function failed in ",
               "funct_python_method::set_function().",o2scl::exc_efailed);
  }
  
  return 0;
}

double funct_python_method::operator()(double x) const {

  // Create a python object from the value
  if (verbose>0) {
    cout << "Creating python object from double." << endl;
  }
  PyObject *pValue=PyFloat_FromDouble(x);
  if (pValue==0) {
    O2SCL_ERR2("Value creation failed in ",
               "funct_python_method::operator().",o2scl::exc_efailed);
  }

  // Set the python function arguments
  if (verbose>0) {
    cout << "Setting item in tuple." << endl;
  }
  int ret=PyTuple_SetItem(pArgs,0,pValue);
  if (ret!=0) {
    O2SCL_ERR2("Tuple set failed in ",
               "funct_python_method::operator().",o2scl::exc_efailed);
  }

  // Call the python function
  if (verbose>0) {
    cout << "Call python function." << endl;
  }
  PyObject *result=PyObject_CallObject(pFunc,pArgs);
  if (result==0) {
    O2SCL_ERR2("Function call failed in ",
               "funct_python_method::operator().",o2scl::exc_efailed);
  }

  double y=PyFloat_AsDouble(result);

  if (verbose>0) {
    cout << "Decref value and result." << endl;
  }
  Py_DECREF(pValue);
  Py_DECREF(result);
  
  if (verbose>0) {
    cout << "Done in funct_python_method::operator()." << endl;
  }
    
  return y;
}

#endif
