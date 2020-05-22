:ref:`O2scl <o2scl>`

Probability Distributions and Markov Chain Monte Carlo
======================================================

\section probdist_subsect Probability distributions

Probability distributions are provided by the C++ standard
library, but multidimensional distributions are not provided. For
the time being, some experimental probability distributions are
being included in \o2.

The one-dimensional probability distributions are children of \ref
o2scl::prob_dens_func and are modeled similar to the C++
standard distributions. The following distributions are 
currently included
- \ref o2scl::prob_dens_gaussian
- \ref o2scl::prob_dens_lognormal
- \ref o2scl::prob_dens_uniform
- \ref o2scl::prob_dens_hist

Multi-dimensional distributions are children of \ref
o2scl::prob_dens_mdim , including
- \ref o2scl::prob_dens_mdim_factor
- \ref o2scl::prob_dens_mdim_gaussian
- \ref o2scl::prob_dens_mdim_biv_gaussian
- \ref o2scl::prob_dens_mdim_amr

Conditional probability distributions are children of 
\ref o2scl::prob_cond_mdim , including 
- \ref o2scl::prob_cond_mdim_gaussian
- \ref o2scl::prob_cond_mdim_fixed_step
- \ref o2scl::prob_cond_mdim_indep

\section mcmc_subsect Markov chain Monte Carlo

The class \ref o2scl::mcmc_para_base performs generic 
MCMC simulations and the child class \ref o2scl::mcmc_para_table 
performs MCMC simulations and stores the results in
a \ref o2scl::table_units object. These classes contain
support for OpenMP and MPI (or both). The class 
\ref o2scl::mcmc_para_cli is a specialized class which 
automatically handles the command-line interface when
using MCMC.

\section ex_mcmc_sect MCMC example

\dontinclude ex_mcmc.cpp
\skip Example:
\until End of example
