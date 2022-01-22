Cold Neutron Star Structure
===========================

:ref:`O2scl_eos <o2scle>`
     
The class :ref:`nstar_cold <nstar_cold>` computes the structure of
non-rotating zero-temperature spherically-symmetric neutron stars,
given a core hadronic equation of state (of type :ref:`eos_had_base
<eos_had_base>`) It automatically tabulates the core EOS, adds a crust
EOS (if necessary) and then uses :ref:`tov_solve <tov_solve>` to
compute the mass-radius curve. It also computes the adiabatic index,
the speed of sound, and determines the possibility of the direct Urca
process as a function of density or radius.

Cold neutron star example
-------------------------
    
.. literalinclude:: ../../../../examples/ex_nstar_cold.cpp
   :language: c++		    
   :start-after: sphinx-example-start

.. literalinclude:: ../../../../examples/ex_nstar_cold.scr      
   :language: none

