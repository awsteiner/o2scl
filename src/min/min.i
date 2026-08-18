# Interface file for o2scl minimizer classes
#
namespace o2scl
py_class_doc |
| Python interface for O2scl class ``%name%``.
| See
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _min:
|
| Minimizer classes
| ==================
|
| :ref:`O2sclpy <o2sclpy>`
#
# Include statements for C++ header file
#
h_include <o2scl/mmin.h>
h_include <o2scl/mmin_simp2.h>
h_include <o2scl/min.h>
h_include <o2scl/min_brent_gsl.h>
#
# Include statement for C++ source code
#
cpp_include <o2scl/min_python.h>
#
# Namespace to use in C++ source code
#
cpp_using std
cpp_using o2scl
#
# ------------------------------------------------------
#
# Class mmin_base
#
class mmin_base<> abstract
- py_class_doc |
| Python interface for O2scl class ``mmin_base``,
| see
| https://awsteiner.org/code/o2scl/html/class/mmin_base.html .
- py_name mmin_base
- int verbose
- int ntrial
- double tol_rel
- double tol_abs
- int last_ntrial
- bool err_nonconv
#
# ------------------------------------------------------
#
# Class mmin_simp2
#
class mmin_simp2<>
- py_class_doc |
| Python interface for O2scl class ``mmin_simp2``,
| see
| https://awsteiner.org/code/o2scl/html/class/mmin_simp2.html .
- py_name mmin_simp2
- parent mmin_base<>
- double size
- double fval
- int print_simplex
#
# ------------------------------------------------------
#
# Class min_base
#
class min_base<> abstract
- py_class_doc |
| Python interface for O2scl class ``min_base``,
| see
| https://awsteiner.org/code/o2scl/html/class/min_base.html .
- py_name min_base
- int verbose
- int ntrial
- double tol_rel
- double tol_abs
- int last_ntrial
- bool err_nonconv
#
# ------------------------------------------------------
#
# Class min_brent_gsl
#
class min_brent_gsl<>
- py_class_doc |
| Python interface for O2scl class ``min_brent_gsl``,
| see
| https://awsteiner.org/code/o2scl/html/class/min_brent_gsl.html .
- py_name min_brent_gsl
- parent min_base<>
- double x_minimum
- double x_lower
- double x_upper
- double f_minimum
- double f_lower
- double f_upper
