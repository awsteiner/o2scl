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

int funct_python::set_function(std::string module, std::string func,
                               std::string path) {
  
  PyRun_SimpleString("import sys");

  // FIXME, we want to make sure this path is not added twice
  std::string pys=((std::string)"sys.path.append('")+path+"')";
  PyRun_SimpleString(pys.c_str());
  
  pName=PyUnicode_FromString(module.c_str());
  if (pName==0) {
    O2SCL_ERR("Create module name failed.",o2scl::exc_efailed);
  }
  pModule=PyImport_Import(pName);
  if (pModule==0) {
    O2SCL_ERR("Load module failed.",o2scl::exc_efailed);
  }
  pArgs=PyTuple_New(1);
  if (pArgs==0) {
    O2SCL_ERR("Create arg tuple failed.",o2scl::exc_efailed);
  }

  func_name=func;
  
  return 0;
}

double funct_python::operator()(double x) const {

  PyObject *pValue=PyFloat_FromDouble(x);
  if (pValue==0) {
    cout << "Value creation failed." << endl;
  }
  int ret=PyTuple_SetItem(pArgs,0,pValue);
  if (ret!=0) {
    cout << "Tuple set failed." << endl;
  }
  
  PyObject *pFunc=PyObject_GetAttrString(pModule,func_name.c_str());
  if (pFunc==0) {
    cout << "Get function failed." << endl;
  }
  PyObject *result=PyObject_CallObject(pFunc,pArgs);
  cout << "result: " << result << endl;
  
  return x;
}

#endif
