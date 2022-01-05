/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#ifdef O2SCL_PYTHON
#include <Python.h>
#endif

#include <o2scl/funct.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_PYTHON

funct_python::funct_python(std::string module, std::string func,
                           std::string path) {
  if (o2scl_settings.py_initialized==false) {
    o2scl_settings.py_init();
  }
  set_function(module,func,path);
}

funct_python::~funct_python() {
  Py_DECREF(pFunc);
  Py_DECREF(pArgs);
  Py_DECREF(pModule);
  Py_DECREF(pName);
}

int funct_python::set_function(std::string module, std::string func,
                               std::string path) {

  // Import the system module so that we can modify the search path
  // to import the module where the function is located
  PyObject *sys_mod=PyImport_ImportModule("sys");
  if (sys_mod==0) {
    O2SCL_ERR2("Import system module failed in",
              "funct_python::set_function().",
              o2scl::exc_efailed);
  }

  // Obtain the path and the number of elements 
  PyObject *sys_path=PyObject_GetAttrString(sys_mod,"path");
  if (sys_path==0) {
    O2SCL_ERR2("Obtain sys.path failed in",
              "funct_python::set_function().",
              o2scl::exc_efailed);
  }
  Py_ssize_t path_size=PySequence_Size(sys_path);
  if (path_size==-1) {
    O2SCL_ERR2("Getting sys.path sequence size failed in",
              "funct_python::set_function().",
              o2scl::exc_efailed);
  }

  // Iterate through each element in the path, and compare with
  // the user-specified path
  
  bool path_found=false;
  
  for(int j=0;j<path_size;j++) {
    PyObject *item=PySequence_GetItem(sys_path,j);
    if (item==0) {
      O2SCL_ERR2("Getting sequence item failed in",
                "funct_python::set_function().",
                o2scl::exc_efailed);
    }
    PyObject *item2=PyObject_Repr(item);
    if (item2==0) {
      O2SCL_ERR2("Getting sequence item unicode failed in",
                "funct_python::set_function().",
                o2scl::exc_efailed);
    }
    PyObject *str=PyUnicode_AsEncodedString(item2,"utf-8","Error");
    if (str==0) {
      O2SCL_ERR2("Getting encoded sequence item failed in",
                "funct_python::set_function().",
                o2scl::exc_efailed);
    }
    const char *cstr=PyBytes_AS_STRING(str);
    if (cstr==0) {
      O2SCL_ERR2("Getting C string from sequence item failed in",
                "funct_python::set_function().",
                o2scl::exc_efailed);
    }
    string cppstr=cstr;
    if (cppstr==path) path_found=true;
    
    Py_DECREF(str);
    Py_DECREF(item2);
    Py_DECREF(item);
  }

  // Is this necessary?
  //PyRun_SimpleString("import sys");

  // If necessary, add the user-specified path to sys.path
  if (path_found==false) {
    std::string pys=((std::string)"sys.path.append('")+path+"')";
    PyRun_SimpleString(pys.c_str());
  }

  // Free the memory associated with the system path and module
  Py_DECREF(sys_path);
  Py_DECREF(sys_mod);
  
  // Get the Unicode name of the user-specified module
  pName=PyUnicode_FromString(module.c_str());
  if (pName==0) {
    O2SCL_ERR2("Create module name failed in ",
              "funct_python::set_function().",o2scl::exc_efailed);
  }

  // Import the user-specified module
  pModule=PyImport_Import(pName);
  if (pModule==0) {
    O2SCL_ERR2("Load module failed in ",
              "funct_python::set_function().",o2scl::exc_efailed);
  }

  // Setup the arguments to the python function
  pArgs=PyTuple_New(1);
  if (pArgs==0) {
    O2SCL_ERR2("Create arg tuple failed in ",
               "funct_python::set_function().",o2scl::exc_efailed);
  }

  // Load the python function
  PyObject *pFunc=PyObject_GetAttrString(pModule,func.c_str());
  if (pFunc==0) {
    O2SCL_ERR2("Get function failed in ",
               "funct_python::set_function().",o2scl::exc_efailed);
  }
  
  return 0;
}

double funct_python::operator()(double x) const {

  // Create a python object from the value
  PyObject *pValue=PyFloat_FromDouble(x);
  if (pValue==0) {
    cout << "Value creation failed." << endl;
  }

  // Set the python function arguments
  int ret=PyTuple_SetItem(pArgs,0,pValue);
  if (ret!=0) {
    cout << "Tuple set failed." << endl;
  }

  // Call the python function
  PyObject *result=PyObject_CallObject(pFunc,pArgs);
  cout << "result: " << result << endl;

  double y=PyFloat_AsDouble(result);

  Py_DECREF(pValue);
  Py_DECREF(result);
    
  return y;
}

#endif
