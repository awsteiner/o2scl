String Manipulation
===================

:ref:`O2scl <o2scl>`

There are a couple classes and functions to help manipulate strings of
text found in ``src/base/string_conv.h``. 

The :ref:`columnify <columnify>` class converts a set of strings into
nicely formatted columns by padding with the necessary amount of
spaces. This class operates on string objects of type ``std::string``,
and also works will for formatting columns of floating-point numbers.
This class is used to provide output for matrices in the functions
:cpp:func:`o2scl::matrix_out()`.

The :ref:`format_float <format_float>` class will reformat double
precision numbers into a form appropriate for HTML or LaTeX documents.

A related function, :cpp:func:`o2scl::screenify()`, reformats a column
of strings into many columns stored row-by-row in a new string array.
It operates very similar to the way the classic Unix command ``ls``
organizes files and directories in multiple columns in order to save
screen space.
    
The function :cpp:func:`o2scl::function_to_double()` converts strings
like ``"pi/3.0"`` and ``"exp(cos(-1.0e-2))"`` to double-precision
floating point numbers using :ref:`find_constants <find_constants>`
and :ref:`calc_utf8 <calc_utf8>`. An alternate version which won't
call the error handler is
:cpp:func:`o2scl::function_to_double_nothrow()`. This latter function
is the one used by ``acol -calc``.

There are also a set of string conversion functions. They partially
duplicate the functionality of C++ Standard Library functions like
``std::stod()`` and ``std::to_string``, but often provide alternate
forms and different options. These functions include:

- :cpp:func:`o2scl::btos()` - boolean value to string
- :cpp:func:`o2scl::dtos()` - double to string 
- :cpp:func:`o2scl::itos()` - integer to string
- :cpp:func:`o2scl::ptos()` - pointer to string
- :cpp:func:`o2scl::szttos()` - ``size_t`` to string
- :cpp:func:`o2scl::stob()` - string to boolean value
- :cpp:func:`o2scl::stod()` - string to double 
  (uses ``std::stod()``)
- :cpp:func:`o2scl::stoi()` - string to integer 
  (uses ``std::stoi()``)
- :cpp:func:`o2scl::stoszt()` - string to ``size_t``
- :cpp:func:`o2scl::s32tod_nothrow()` - char32 string to double which never
  throws an exception
- :cpp:func:`o2scl::stod_nothrow()` - string to double which never
  throws an exception
- :cpp:func:`o2scl::stoi_nothrow()` - string to integer which never
  throws an exception
- :cpp:func:`o2scl::stoszt_nothrow()` - ``size_t`` to string which never
  throws an exception

Other functions in ``src/base/string_conv.h`` are:

- :cpp:func:`o2scl::size_of_exponent()`
- :cpp:func:`o2scl::count_words()`
- :cpp:func:`o2scl::has_minus_sign()`
- :cpp:func:`o2scl::is_number()`
- :cpp:func:`o2scl::find_constant()`
- :cpp:func:`o2scl::split_string()`
- :cpp:func:`o2scl::split_string_delim()`
- :cpp:func:`o2scl::rewrap()`
- :cpp:func:`o2scl::rewrap_ignore_vt100()`
- :cpp:func:`o2scl::utf8_to_char32()`
- :cpp:func:`o2scl::char32_to_utf8()`
- :cpp:func:`o2scl::rewrap_keep_endlines()`
- :cpp:func:`o2scl::string_to_uint_list()`
- :cpp:func:`o2scl::parse_fortran_format()`
- :cpp:func:`o2scl::string_to_char_array()`

String example
--------------

.. literalinclude:: ../../../examples/ex_string.cpp
   :language: c++		    
   :start-after: sphinx-example-start


