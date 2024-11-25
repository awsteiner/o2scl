# Interface file for std::string class
#
namespace o2scl
py_class_doc |
| Python interface for C++ class ``%name%``.
dll_name o2scl
rst_header |
| .. _string:
|
| String class
| ============
|
| :ref:`O2sclpy <o2sclpy>`
|
| This python interface is not intended to provide the full 
| functionality of the corresponding C++ class.
# 
# Include statements for C++ header file
# 
h_include <o2scl/lib_settings.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/string_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
#
# Additional python headers
#
#
# The interpolation types
py_header itp_linear=1
py_header itp_cspline=2
py_header itp_cspline_peri=3
py_header itp_akima=4
py_header itp_akima_peri=5
py_header itp_monotonic=6
py_header itp_steffen=7
py_header itp_nearest_neigh=8
#
# ------------------------------------------------------
#
# Class string
# 
class std::string
- py_class_doc |
| Note that std_string objects are not "immutable" like Python
| strings.
- py_name std_string
- std_cc                             
- function length
  - size_t
- function operator[]
  - char &
  - size_t n
- function resize
  - void
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: an int
|     """
|     return self.length()
| 
| def init_bytes(self,s):
|     """
|     Initialize the string from a Python bytes object
|
|     | Parameters:
|     | *s*: a Python bytes string
|     """
|
|     f=self._link.o2scl.o2scl_char_p_to_string
|     # AWS, 11/23/24: Note that, in ctypes, POINTER(char) is not 
|     # the same as char_p, which is a bit confusing. 
|     f.argtypes=[ctypes.c_int,ctypes.POINTER(ctypes.c_char),
|                 ctypes.c_void_p]
|     f(ctypes.c_int(len(s)),ctypes.c_char_p(s),self._ptr)
| 
|     return
|
| def to_bytes(self):
|     """
|     Copy the string to a Python bytes object
|
|     Returns: a Python bytes string
|     """
|
|     n=ctypes.c_int(self.length())
|                          
|     # AWS, 11/23/24: The function create_string_buffer()
|     # is the special ctypes function to create a POINTER(char) 
|     # object for use in ctypes. Note that, in ctypes,
|     # POINTER(char) is not the same as char_p, which is a bit
|     # confusing.
| 
|     b=ctypes.create_string_buffer(self.length())
|     f=self._link.o2scl.o2scl_string_to_char_p
|     f.argtypes=[ctypes.c_void_p,
|                 ctypes.POINTER(ctypes.c_int),
|                 ctypes.POINTER(ctypes.c_char)]
|     f(self._ptr,ctypes.byref(n),b)
| 
|     return b.value
