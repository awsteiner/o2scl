Nuclei and nuclear masses
=========================

:ref:`O2scl_part <o2sclp>`

A basic data structure for nuclei, class :ref:`nucleus <nucleus>`, is
implemented as a child of :ref:`part_tl <part_tl>`.

Nuclear masses are given as children of :ref:`nucmass <nucmass>` and
are generally of two types: tables of masses (children of
:ref:`nucmass_table <nucmass_table>`) or formulas which generate
masses from a set of parameters and can be fit to data (
:ref:`nucmass_fit <nucmass_fit>`).

There are eleven mass table types currently included.

- :ref:`nucmass_ame <nucmass_ame>`: data from 
  \ref [Audi95]_, [Audi03]_, [Audi12]_ and [Wang12]_, or
  [Huang17]_ and [Wang17]_.
- :ref:`nucmass_mnmsk <nucmass_mnmsk>` and :ref:`nucmass_mnmsk_exp
  <nucmass_mnmsk_exp>`: masses from [Moller95]_ and [Moller97]_.
- :ref:`nucmass_hfb <nucmass_hfb>` and :ref:`nucmass_hfb_sp
  <nucmass_hfb_sp>`: masses from [Goriely02]_, [Samyn04]_, or
  [Goriely07]_.
- :ref:`nucmass_dz_table : masses from [Duflo95]_
- :ref:`nucmass_ktuy : masses from [Koura00]_ and [Koura05]_
- :ref:`nucmass_wlw : masses from \ref Wang10, \ref Wang10b,
\ref Liu11, \ref Wang11, or \ref Wang14
- :ref:`nucmass_sdnp : masses from \ref Stoitsov03 or \ref
Dobaczewski04.
- :ref:`nucmass_dglg : masses from \ref Delaroche10 

The mass formulas which can be fit to data are
- :ref:`nucmass_semi_empirical : simple 5 parameter 
semi-empirical method
- :ref:`nucmass_frdm : macroscopic part of FRDM from \ref Moller95
- :ref:`nucmass_dz_fit and :ref:`nucmass_dz_fit_33 : 
10- and 33-parameter mass formulas from \ref Duflo95.
- :ref:`nucmass_dvi : 10-parameter formula from \ref Dieperink09
with :ref:`nucmass_ibm_shell for shell effects
    
In order to create a set of nuclei stored in a <tt>std::vector</tt>
object, one can use \ref o2scl_part::nucdist_set().

\section ex_nucmass_fit_sect Nuclear mass fit example
    
\dontinclude ex_nucmass_fit.cpp
\skip Example:
\until End of example

\section ex_nucmass_sect Nuclear mass example

\dontinclude ex_nucmass.cpp
\skip Example:
\until End of example

\image html ex_nucmass_se.png ""
\image html ex_nucmass_mnmsk.png ""
\image html ex_nucmass_dz96.png ""
\image html ex_nucmass_ame03.png ""
\image html ex_nucmass_hfb14.png ""
\image html ex_nucmass_hfb21.png ""
\image html ex_nucmass_hfb27.png ""
\image html ex_nucmass_ktuy05.png ""
\image html ex_nucmass_dvi.png ""
\image html ex_nucmass_ws32.png ""
\image html ex_nucmass_ws36.png ""
