Yet ANother Interface between C++ and python
============================================

.. warning:: Very experimental. I probably shouldn't have written my
             own code to interface python and C++, but I did anyway,
             and here it is.
             
``yanic: <interface file> <c++ output prefix> <python prefix> <rst
prefix>``
             
Guide for interface files
-------------------------
 
- ``[]`` arguments are optional
- ``<>`` arguments are required and do not have any whitespace
- ``{}`` arguments can have whitespace but no carriage returns

Types, e.g. {type}, {parameter type}, {return type} are specified in
the following way::

  [static] [const] [std::shared_ptr] <type name> [*] [&] [**]

Non-alphabetic characters in class names are always converted to
underscores.

Comments (all lines beginning with ``#`` are comments). Comments may
appear anywhere, including inside a class or function definition.

Header items
------------

- Namespace specification::

    namespace <name>

  The namespace is used to name the automatically name the ``extern
  C`` functions which access the class.
    
- Template for class documentation, using ``%name%`` to refer to the
  class name. Specified either as::

    py_class_doc {}

  or as::

    py_class_doc |
    | {}
    | {} 
    ...

- Name of dll to load::

    dll_name <name>

- Header for .rst files::

    rst_header {}

  or::

    rst_header |
    | {}
    | {} 
    ...

- Include statements for C++ header file::
    
    h_include <file (including quotes or angle brackets)>

- Include statement for C++ source code::

    cpp_include <file (including quotes or angle brackets)>

- Namespaces to use in C++ source code::

    cpp_using <namespace>

- Additional python header lines::

    py_header {}

Functions
---------

- Function definitions::

    function <function name>
    - {return type}
    - py_name <python name>      
    - {parameter type} <parameter name> [default value]
    ...

  For each function, an ``extern C`` wrapper is created with the
  suffix ``_wrapper`` and then a python function is created to call
  that wrapper. The python name is optional, but when present must
  be after the return type and before the variable list.

Classes
-------
    
- Class definitions::

    class <class name> ["abstract"]

  If the ``abstract`` label is appended, then ``__init__()`` is tagged
  as an ``@abstractmethod`` in python. If the class is not abstract and the
  ``no_def_cons`` tag is not given (see below), then a ``create``
  function is created to create an object. A ``free`` function is
  always created to destroy an object. If ``std_cc`` is specified,
  then a C function named ``copy`` is created along with an analogous
  ``__deepcopy__`` python method.
    
- Python name of class (optional)::

    - py_name <name>

- Optional line which should be specified if the class defines both
  ``class (const class &)`` and ``class &operator=(const class &)``.
  This allows one to define the python ``__deepcopy__`` method::
  
  - std_cc                             

- Optional line which should be specified if the class has no default
  constructor::
  
  - no_def_cons

- Parent class (multiple parents not currently supported)::

    - parent <parent class name>

- Python documentation for this class (overrides template
  specification above)::

    - py_class_doc {}

  or::

    - py_class_doc |
    | {}
    | {} 
    ...

- Class member data::

  - {type} <name>

  Get and set methods for class member data are generated. For
  standard C types, ``std::string`` objects, and shared pointer
  objects, the get and set methods pass by value (i.e. they imply
  a copy). For classes from the interface, the get methods return 
  references and the set methods pass by value.
  
- Class member function definitions are of the following form.
  The return type and parameter specifications must begin with
  two spaces::

    - function <function name>
      - {return type}
      - {parameter type} <parameter name> [default value]
      ...
    
- Extra python code for the class::

    - extra_py {}

  or::

    - extra_py |
    | {}
    | {}
    ...

  The extra python code is prepended by four spaces to conform
  with the indentation style used by yanic.

- Class constructor with parameters. The parameter specifications must
  begin with two spaces::

    - cons <python constructor name>
      - py_name <python name>      
      - {parameter type} <parameter name>
      ...

Other objects
-------------
      
- Specification of a shared pointer::

    - shared_ptr <class name>

  Shared pointers imply the creation of a ``create`` function to
  create a shared pointer to a default object, a ``free`` function to
  free the memory associated with the shared pointer (which may or may
  not free the underlying object), and a pointer function which gets a
  raw pointer to the underlying object. Using shared pointers for
  objects which do not have a default constructor is not yet
  supported.

  * Python name of class for the shared pointer (must begin with
    two spaces)::

      - py_name <name>

Constraints
-----------

- Global functions and member functions may be overloaded, but
  only if they are given different python names.

Todos
-----

.. todo:: 

   In yanic:

   - Need to fix function names in case where there is no namespace.
   - Simplify code duplication in parsing: reading global and member
     functions should be the same
   - Allow use of numpy.arange for uniform_grid arguments
   - Document .i format
   - Make sure data members named 'del' are properly renamed without
     hacking, e.g. with a py_name argument
   - Make sure negative size_t arguments are rejected in python
     
Details
-------

Handling of function arguments:

- C type (bool, char, double, float, int, size_t): Convert from
  a Python object to ``ctypes.c_<type>`` in the Python wrapper.
- reference to C-type: The handling of these references depends
  on whether the reference is labeled as ``io`` (input and output)
  or ``out`` (output only). These references are included in
  return values in the python wrapper function, and if
  the label ``io`` is given, then they are also input parameters.
- pointer to C-type: not yet implemented
- ``std::string``: Use ``char *`` in C the wrapper. Convert python
  string to bytes object and then to char * in python code.
- reference to ``std::string``: Use ``void *&`` in the C wrapper ...
- std_vector & - table::line_of_data: convert to an
  std::vector<double>
- std::vector<double> & - uniform_grid::vector

Return values:

- C type bool: return Python ``True`` or ``False``
- C type char: return a one-character bytes object
- C type double or float: return a Python float
- C type int or size_t: return a Python int  
- reference to C type: this is supported currently only for operator[]
  and operator()
- ``std::string``: Return a Python bytes object
- ``std::string &``: Return a std_string object

Special funcions:

- Array-indexing, ``operator[]`` functions are translated to
  ``__getitem__`` functions on the python side. If the ``operator[]``
  function returns a non-const reference rather than a value,
  then a ``__setitem__`` function is also created.

