:ref:`O2scl <o2scl>`

Simulated Annealing
===================

Minimization by simulated annealing is performed by descendants of
sim_anneal (see \ref o2scl::anneal_gsl). 

\comment
Because simulated
annealing is particularly well-suited to parallelization, a
multi-threaded minimizer analogous to \ref o2scl::anneal_gsl is
given in \ref o2scl::anneal_mt. This header-only class uses the
Boost libraries and requires it for use.
\endcomment

\section ex_anneal_sect Simulated annealing example

This example minimizes the function
\f[
f(x,y) = J_0(x-2) J_0(y+3)
\f]
over \f$ (x,y) \f$ where \f$ J_0(x) \f$ is the Bessel function
given in \c gsl_sf_bessel_J0. The initial guess at \f$ (9,9) \f$
is far away from the global minimum. 

The plot below plots the function above, the initial guess, and
the minimum obtained by the example program.

\image html ex_anneal_plot.png "Results of simulated annealing minimization"
\comment
\image latex ex_anneal_plot.eps "Results of simulated annealing minimization" width=9cm
\endcomment

\dontinclude ex_anneal.cpp
\skip Example:
\until End of example
