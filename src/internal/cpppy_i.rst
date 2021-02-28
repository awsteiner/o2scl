# Interface guide
# 
# [] arguments are optional
# <> arguments are required and do not have any whitespace
# {} arguments can have whitespace but no carriage returns

# Types, e.g. {type}, {parameter type}, {return type} are specified in
# the following way:
# 
# [static] [const] [std::shared_ptr] <type name> [*] [&] [**]

# Non-alphabetic characters in class names are always converted to
# underscores.

# Comments (all lines beginning with '#' are comments). Comments # may
appear anywhere, including inside a class or function definition.

# Header items:

# Namespace specification
namespace <name>

# Template for class documentation, using %name% to refer to the
# class name. Specified either as

py_class_doc {}

# or as

py_class_doc |
| {}
| {} 
...

# Name of dll to load

dll_name <name>

# Header for .rst files

rst_header {}

or 

rst_header |
| {}
| {} 
...

# Include statements for C++ header file

h_include <file (including quotes or angle brackets)>

# Include statement for C++ source code

cpp_include <file (including quotes or angle brackets)>

# Namespaces to use in C++ source code

cpp_using <namespace>

# Additional python header lines

py_header {}

# Class definitions

class <class name> ["abstract"]

# Python name of class (optional)

- py_name <name>

# Optional line which should be specified if the class defines both 
# class (const class &) and class &operator=(const class &) . This
# allows one to define the python __deepcopy__ method.
  
- std_cc                             

# Optional line which should be specified if the class has no
# default constructor
  
- no_def_cons

# Parent class

- parent <parent class name>

# Python documentation for this class (overrides template
# specification above)

- py_class_doc {}

# or

- py_class_doc |
| {}
| {} 
...

# Class member data

- {type} <name>
  
# Class member function definitions are of the following form.
# The return type and parameter specifications must begin with
# two spaces.

- function <function name>
  - {return type}
  - {parameter type} <parameter name>
  ...
    
# Extra python code for the class

- extra_py {}

# or

- extra_py |
| {}
| {}
...

# Class constructor with parameters. The parameter specifications must
# begin with two spaces.

- cons <python constructor name>
  - {parameter type} <parameter name>
  ...

# Specification of a shared pointer

- shared_ptr <class name>

# Python name of class for the shared pointer

- py_name <name>


