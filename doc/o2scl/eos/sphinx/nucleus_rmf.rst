Nuclear structure in the Hartree approximation
==============================================

See class :ref:`nucleus_rmf <nucleus_rmf>`.

Nucleus in the Hartree approximation example
--------------------------------------------

This example uses the NL3 EOS as implemented in :ref:`eos_had_rmf
<eos_had_rmf>` and computes the structure of :math:`^{208}\mathrm{Pb}`
using :ref:`nucleus_rmf <nucleus_rmf>`. The results are stored in
:ref:`table_units <o2scl:table_units>` objects and output to HDF
files. It also computes the energy of three unoccupied neutron and
proton levels.

.. literalinclude:: ../../../examples/ex_nucleus_rmf.cpp
   :language: c++		    
   :start-after: sphinx-example-start

.. image:: ../../../examples/plot/ex_nuc_prof.png
   :width: 60%	   
   :alt: alt text

Typical output:

.. literalinclude:: ../../../../examples/ex_nucleus_rmf.scr      
   :language: none
