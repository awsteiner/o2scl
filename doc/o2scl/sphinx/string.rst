String Manipulation
===================

:ref:`O2scl <o2scl>`

.. contents:: 

There are a couple classes and functions to help manipulate
strings of text. Conversion routines for \c std::string 
objects are given in ``src/base/string_conv.h`` and include

- :cpp:func:`o2scl::btos()` - boolean value to string
- :cpp:func:`o2scl::dtos()` - double to string
- :cpp:func:`o2scl::itos()` - integer to string
- :cpp:func:`o2scl::ptos()` - pointer to string
- :cpp:func:`o2scl::stob()` - string to boolean value
- :cpp:func:`o2scl::stod()` - string to double 
  (uses C++11 ``std::stod()``)
- :cpp:func:`o2scl::stoi()` - string to integer 
  (uses C++11 ``std::stoi()``)
- :cpp:func:`o2scl::stoszt()` - string to ``size_t``
- :cpp:func:`o2scl::szttos()` - ``size_t`` to string

(While :cpp:func:`o2scl::dtos()` and similar functions have been
implemented in ``std::to_string``, the O\ :sub:`2`\ scl versions have
been written to allow a bit more flexibility.)

There are also a set of conversion functions which return
and integer error code instead of throwing an exception

- :cpp:func:`o2scl::stoi_nothrow()` - string to integer 
- :cpp:func:`o2scl::stoszt_nothrow()` - ``size_t`` to string

The :ref:`columnify <columnify>` class converts a set of
strings into nicely formatted columns by padding with the
necessary amount of spaces. This class operates on string objects
of type ``std::string``, and also works will for formatting columns
of floating-point numbers.  This class is used to provide output
for matrices in the functions :cpp:func:`o2scl::matrix_out()`.

The :ref:`format_float <format_float>` class will reformat double
precision numbers into a form appropriate for HTML or LaTeX documents.

A related function, :cpp:func:`o2scl::screenify()`, reformats a column
of strings into many columns stored row-by-row in a new string array.
It operates very similar to the way the classic Unix command ``ls``
organizes files and directories in multiple columns in order to save
screen space.
    
The function :cpp:func:`o2scl::function_to_double()` converts strings
like ``"1.0/3.0"`` and ``"exp(cos(-1.0e-2))"`` to double-precision
floating point numbers using :ref:`calculator <calculator>`.

The function :cpp:func:`o2scl::size_of_exponent()` returns 2 or 3,
depending on the number of characters in the exponent when a floating
point number is output to the screen.

Finally, the function :cpp:func:`o2scl::count_words()` counts the
number of "words" in a string, which are delimited by whitespace.
