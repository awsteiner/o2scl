Python Integration
==================

:ref:`O2scl <o2scl>`

O₂scl has an associated experimental python package, O₂sclpy. O₂sclpy
assists in plotting HDF5 files generated by O₂scl using matplotlib and
yt. Many O₂scl classes also have a Python interface created by
:ref:`Yet ANother Interface between C++ and python` in O₂sclpy.
Finally, O₂sclpy also has a few experimental Python machine learning
classes based on ``sklearn``, ``tensorflow``, and ``torch`` which have C++
interface in O₂scl. O₂sclpy is separately documented at
https://awsteiner.org/code/o2sclpy . O₂sclpy's version numbers
which loosely follow O₂scl releases. The release version of the python
package can be installed with, e.g. ``pip3 install o2sclpy`` or
obtained from https://pypi.python.org/pypi/o2sclpy. The development
version of O₂sclpy can be obtained from
https://github.com/awsteiner/o2sclpy .

All of the plots in this documentation are created by ``o2graph``
which is part of O₂sclpy.

There are several global functions with extern "C" linkage used to
communicate between python and the
:cpp:class:`o2scl_acol::acol_manager` class defined in acolm.h.

If python support is enabled when installing O₂scl, then the classes
:ref:`interpm_python <interpm_python>`, :ref:`gmm_python
<gmm_python>`, :ref:`kde_python <kde_python>`, and
:ref:`nflows_python <nflows_python>` provide a C++
interface to similarly named Python classes in O₂sclpy. See
:ref:`Python support` for more details regarding installing
O₂scl with Python support.

All of the plots in this documentation are created by ``o2graph``
which is part of O₂sclpy.
