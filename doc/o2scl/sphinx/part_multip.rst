Multiprecision Support for Particles and EOSs
=============================================

:ref:`O2scl <o2scl>`

Multiprecision support for :ref:`O2scl <o2scl>` classes is in
progress. In the future, classes with a ``_tl`` suffix are compatible
with a range of floating-point types. Support for types beyond
``double`` is still experimental.

List of classes which currently use multiprecision:

- Fermionic statistics :ref:`fermion_rel_ld <fermion_rel_ld>` and
  :ref:`fermion_rel_cdf25 <fermion_rel_cdf25>` and 
  :ref:`fermion_deriv_rel_tl <fermion_deriv_rel_tl>`.
- Lepton equation of state: :ref:`eos_leptons <eos_leptons>`

.. I have to reference fermion_rel_deriv_tl rather than the typedef
   because sphinx doesn't resolve the references to typedefs.
