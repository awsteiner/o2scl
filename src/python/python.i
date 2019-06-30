%define DOCSTRING
"O2scl docstring."
%enddef
%module(docstring=DOCSTRING) o2scl

%{
#define SWIG_FILE_WITH_INIT
#include "../base/misc.h"
%}

%include "../base/misc.h"

