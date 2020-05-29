Neutron Star Outer Crust
========================

The class :ref:`eos_crust <eos_crust>` computes the composition and
EOS of the outer crust of a neutron star from
:math:`10^{6}~\mathrm{g}/\mathrm{cm}^3` up to the neutron drip density
of :math:`4 \times 10^{11}~\mathrm{g}/\mathrm{cm}^3`. It uses a simple
expression for the lattice energy and works with any nuclear mass
formula specified as an object of type :ref:`nucmass <nucmass>`.

Outer crust example
-------------------

This example computes the outer crust EOS from the 
:ref:`eos_crust <eos_crust>` class and then compares it with the numerical
values from the original paper in [Baym71tg]_ as stored in the
file ``examples/ex_bps_input.txt``. Then it computes several
crusts from [Shen11b]_, [Steiner12]_, and [Newton13]_.
    
.. literalinclude:: ../../../examples/ex_eos_crust.cpp
   :language: c++		    
   :start-after: sphinx-example-start

