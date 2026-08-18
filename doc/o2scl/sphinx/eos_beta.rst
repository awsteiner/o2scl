Equations of State in Beta Equilibrium
======================================

:ref:`O2scl <o2scl>`
     
There are a few ways to compute equations of state (EOSs) in beta
equilibrium.

The most user-friendly is the :ref:`nstar_cold <nstar_cold>` EOS,
which will automatically compute a :ref:`table <table>` object which
holds the beta-equilibrium EOS for any child of the :ref:`eos_had_base
<eos_had_base>` class (but this currently only supports EOSs which
describe neutrons and protons only).

The second method is :cpp:func:`eos_base::beta_eq_T0()` which
is designed to handle beta equilibrium for EOSs involving hadrons,
quarks, or hyperons. This infrastructure is still experimental.

