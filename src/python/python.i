%define DOCSTRING
"O2scl docstring."
%enddef
%module(docstring=DOCSTRING) o2scl_base

namespace o2scl {
  double fermi_function(double E, double mu, double T, double limit=40.0);
}

  

