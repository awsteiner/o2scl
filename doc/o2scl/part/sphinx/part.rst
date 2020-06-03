Particles
=========
    
:ref:`O2scl_part <o2sclp>`

These classes in the library O\ :sub:`2`\ scl_part calculate the
thermodynamic properties of interacting and non-interacting
quantum and \classical particles.

Particle data classes
---------------------

The class :cpp:class:`o2scl::part_tl` is the basic structure for
a particle:

- :cpp:var:`o2scl::part_tl::m` - mass (i.e. rest mass), :math:`m`
- :cpp:var:`o2scl::part_tl::g` - degeneracy factor (e.g. :math:`g=2j+1`)
- :cpp:var:`o2scl::part_tl::n` - number density, :math:`n`
- :cpp:var:`o2scl::part_tl::ed` - energy density, :math:`\varepsilon`
- :cpp:var:`o2scl::part_tl::pr` - pressure, :math:`P`
- :cpp:var:`o2scl::part_tl::en` - entropy density, :math:`s`
- :cpp:var:`o2scl::part_tl::ms` - effective mass, :math:`m^{*}`
- :cpp:var:`o2scl::part_tl::mu` - chemical potential, :math:`\mu`
- :cpp:var:`o2scl::part_tl::nu` - effective chemical potential, :math:`\nu`
- :cpp:var:`o2scl::part_tl::inc_rest_mass` - True if the rest mass is included
  (default true)
- :cpp:var:`o2scl::part_tl::non_interacting` - False if the particle 
  includes interactions (default true)

The data members :cpp:var:`o2scl::part_tl::ms` and
:cpp:var:`o2scl::part_tl::nu` allow one to specify modifications to
the mass and the chemical potential due to interactions. This allows
one to calculate the properties of particle due to interactions so
long as the basic form of the free-particle dispersion relation is
unchanged, i.e.

.. math::

   \sqrt{k^2+m^2} - \mu \rightarrow \sqrt{k^2+m^{* 2}} - \nu 

If the particle is non-interacting, then :cpp:var:`o2scl::part_tl::nu` and
:cpp:var:`o2scl::part_tl::ms` are sometimes used by O\ :sub:`2`\ scl_part
functions for temporary storage.

If :cpp:var:`o2scl::part_tl::inc_rest_mass` is \c true (as is the
default in all of the classes except :cpp:class:`o2scl::nucleus`),
then all functions include the rest mass (stored in
:cpp:var:`o2scl::part_tl::m`) energy density in the energy density,
the "mu" functions expect that the rest mass is included in
:cpp:var:`o2scl::part_tl::mu` or 
:cpp:var:`o2scl::part_tl::nu` as input and the "density" functions
output  :cpp:var:`o2scl::part_tl::mu` or
:cpp:var:`o2scl::part_tl::nu` including the rest mass. Note that it
is assumed that :cpp:var:`o2scl::part_tl::m` contains the rest mass
even if the particle is interacting and an effective mass is stored in
:cpp:var:`o2scl::part_tl::ms`.
    
When :cpp:var:`o2scl::part_tl::inc_rest_mass` is true, antiparticles are
implemented by choosing the antiparticle chemical potential to be
:math:`- \mu`. When :cpp:var:`o2scl::part_tl::inc_rest_mass` is false,
there is an ambiguity in the relative definitions of the rest mass
contribution for the antiparticles and the combination of both
particles and antiparticles. Define energy density for particles
including the rest mass contribution as :math:`\varepsilon_+`, and
the energy density without the rest mass contribution as 
:math:`\tilde{\varepsilon}_{+} = \varepsilon_{+} - n_{+} m` .
Similarly, for antiparticles, we have :math:`\tilde{\varepsilon}_- =
\varepsilon_- - n_- m`. The total energy density including the
rest mass contribution is then :math:`\varepsilon = \varepsilon_{+} +
\varepsilon_-` and without the rest mass contribution 
:math:`\tilde{\varepsilon} \equiv \varepsilon - (n_{+}-n_-) m`. Then,

.. math::

   \begin{eqnarray}
   \tilde{\varepsilon} & = & 
   \varepsilon_+ - n_{+} m + \varepsilon_- + n_- m \nonumber \\
   & = & \varepsilon_+ - n_{+} m + \varepsilon_- - n_- m + 
   2 n_- m \nonumber \\
   & = & \tilde{\varepsilon}_+ + \tilde{\varepsilon}_- + 2 n_- m 
   \nonumber
   \end{eqnarray}
   
Similarly, for the chemical potentials, we have 

.. math::

   \tilde{\mu}_+ \equiv \frac{\partial \tilde{\varepsilon}_+}{n_+} = 
   \mu_+ - m \quad \mathrm{and} \quad 
   \tilde{\mu}_- \equiv \frac{\partial \tilde{\varepsilon}_-}{n_-} = 
   \mu_- - m 

thus :math:`\tilde{\mu}_- = - \tilde{\mu}_+ - 2 m` . This bookkeeping
is handled by :cpp:func:`o2scl::part_tl::anti()`, the
:cpp:func:`o2scl::fermion_thermo_tl::pair_mu()`, and the
:cpp:func:`o2scl::fermion_thermo_tl::pair_density()`, functions.

The thermodynamic identity used to compute the pressure for
interacting particles is

.. math::
   
   P = -\varepsilon + s T + \nu n

where :cpp:var:`o2scl::part_tl::nu` is used. This way, the particle class
doesn't need to know about the structure of the interactions to
ensure that the thermodynamic identity is satisfied. Note that in
the \o2e library, where in the equations of state the normal
thermodynamic identity is used

.. math::
   
   P = -\varepsilon + s T + \mu n

Frequently, the interactions which create an effective chemical
potential which is different than :cpp:var:`o2scl::part_tl::mu` thus create
extra terms in the pressure and the energy density for the given
equation of state.

The :cpp:class:`o2scl::fermion_tl` class is a child of
:cpp:class:`o2scl::part_tl` which contains data members for the Fermi
momentum and energy gap. The :cpp:class:`o2scl::boson` class contains
an extra data member for the condensate. The :cpp:class:`o2scl::quark`
class is a descendant of the :cpp:class:`o2scl::fermion_tl` class
which contains extra data members for the quark condensate and the
contribution to the bag constant. Nuclei are represented by the
:cpp:class:`o2scl::nucleus` class and documented in nuclei_section.

Units
-----

Factors of :math:`\hbar, c` and :math:`k_B` have been removed
everywhere, so that mass, energy, and temperature all have the
same units. Number and entropy densities have units of mass cubed
(or energy cubed). The particle classes can be used with any
system of units which is based on powers of one unit, i.e. 
:math:`[n] = [T]^3 = [m]^3 = [P]^{3/4} = [\varepsilon]^{3/4}`, etc.

Classes for particle thermodynamics
-----------------------------------

At zero temperature, the thermodynamic properties of fermions can
be computed using :cpp:class:`o2scl::fermion_zerot_tl`. The class 
:cpp:class:`o2scl::classical_thermo_tl` computes the properties of particles
in the classical limit.

At finite temperature, there are different classes corresponding to
different approaches to computing the integrals over the distribution
functions. The approximation scheme from Johns96 is used in
:cpp:class:`o2scl::boson_eff` and :cpp:class:`o2scl::fermion_eff`. An
exact method employing direct integration of the distribution
functions is used in :cpp:class:`o2scl::boson_rel` and
:cpp:class:`o2scl::fermion_rel_tl`, but these are necessarily quite a
bit slower. All of these classes use expansions to give ensure
comparably accurate results in the degenerate and non-degenerate
limits.

The class :cpp:class:`o2scl::fermion_eff` usually works to within about 1
part in :math:`10^4`, but can be as bad as 1 part in :math:`10^2`
in some more extreme cases. The default settings for 
:cpp:class:`o2scl::fermion_rel_tl` give an accuracy of at least 1 part in
:math:`10^6` (and frequently better than this). For 
:cpp:class:`o2scl::fermion_rel_tl`, the accuracy can be improved to 1 part in
:math:`10^{10}` by decreasing the integration tolerances.

The class :cpp:class:`o2scl::fermion_nonrel_tl` assumes a non-relativistic
dispersion relation for fermions. It uses an exact method for both
zero and finite temperatures. The non-relativistic integrands are
much simpler and :cpp:class:`o2scl::fermion_nonrel_tl` uses the appropriate
GSL functions (which are nearly exact) to compute them.

Thermodynamics with derivatives
-------------------------------

Sometimes it is useful to know derivatives like :math:`ds/dT` in
addition to the energy and pressure.
The class :cpp:class:`o2scl::part_deriv_press_tl` stores the three
derivatives which correspond to second derivatives
of the pressure

.. math::
   
   \left(\frac{\partial n}{\partial \mu}\right)_{T}, \quad
   \left(\frac{\partial n}{\partial T}\right)_{\mu}, \quad
   \mathrm{and} \quad \left(\frac{\partial s}{\partial
   T}\right)_{\mu} \quad . 

All other first derivatives of the thermodynamic functions can be
written in terms of these three.

The new data classes are 
:cpp:class:`o2scl::part_deriv_tl` and :cpp:class:`o2scl::fermion_deriv_tl` 
which store the basic particle thermodynamics described
above with these additional three derivatives.

There are three classes which compute these derivatives for
fermions and classical particles. The class 
:cpp:class:`o2scl::classical_deriv_thermo_tl` handles the nondegenerate limit,
:cpp:class:`o2scl::fermion_deriv_rel_tl` handles fermions and 
:cpp:class:`o2scl::fermion_deriv_nr_tl` handles nonrelativistic fermions.
The class :cpp:class:`o2scl::fermion_deriv_thermo_tl` is a base
class for :cpp:class:`o2scl::fermion_deriv_rel_tl` and uses
degenerate and nondegenerate expansions to evaluate
both the base thermodynamic quantities and the three 
derivatives from :cpp:class:`o2scl::part_deriv_press_tl` .

The function :cpp:func:`o2scl::part_deriv_tl::deriv_f()` computes
the derivatives which are second derivatives of the
free energy from the three computed above.

Other derivatives
-----------------
    
For the derivative of the entropy with respect
to the chemical potential, there is a Maxwell relation

.. math::
   
   \left(\frac{\partial s}{\partial \mu}\right)_{T,V} =
   \left(\frac{\partial^2 P}{\partial \mu \partial T}\right)_{V} =
   \left(\frac{\partial^2 P}{\partial T \partial \mu}\right)_{T,V} =
   \left(\frac{\partial n}{\partial T}\right)_{\mu,V}

The first derivatives of the energy density can be computed using
the thermodynamic identity:

.. math::

   \left(\frac{\partial \varepsilon}{\partial \mu}\right)_{T,V}=
   \mu \left(\frac{\partial n}{\partial \mu}\right)_{T,V}+
   T \left(\frac{\partial s}{\partial \mu}\right)_{T,V}

.. math::

   \left(\frac{\partial \varepsilon}{\partial T}\right)_{\mu,V}=
   \mu \left(\frac{\partial n}{\partial T}\right)_{\mu,V}+
   T \left(\frac{\partial s}{\partial T}\right)_{\mu,V}
    
Most of the other common derivatives which are used 
are those which can be obtained by second derivatives
of the Gibbs free energy, :math:`G = F + P V`.

.. math::
   
   \begin{eqnarray}
   \left(\frac{\partial^2 G}{\partial T^2}\right)_{P,\{N_i\}} &=&
   -\left( \frac{\partial S}{\partial T} \right)_{P,\{N_i\}}
   = - \frac{N c_P}{T} 
   \nonumber \\
   \left(\frac{\partial^2 G}{\partial T \partial P}\right)_{\{N_i\}} &=&
   \left( \frac{\partial V}{\partial T} \right)_{P,\{N_i\}}
   = V \alpha
   \nonumber \\
   \left(\frac{\partial^2 G}{\partial P^2}\right)_{T,\{N_i\}} &=&
   \left( \frac{\partial V}{\partial P} \right)_{T,\{N_i\}}
   = - V \kappa_T
   \nonumber
   \end{eqnarray}

Other common derivatives are the heat capacity per particle at
constant volume, :math:`c_V`, and the speed of sound, :math:`( d P / d
\varepsilon)_{\{N_i\},S}`. These derivatives are computed by functions
in :cpp:class:`o2scl::deriv_thermo_base_tl` from the three second
derivatives of the pressure stored in a
:cpp:class:`o2scl::part_deriv_tl` or
:cpp:class:`o2scl::fermion_deriv_tl` object.

- :cpp:func:`o2scl::deriv_thermo_base_tl::heat_cap_ppart_const_vol()`
- :cpp:func:`o2scl::deriv_thermo_base_tl::heat_cap_ppart_const_press()`
- :cpp:func:`o2scl::deriv_thermo_base_tl::compress_adiabatic()`
- :cpp:func:`o2scl::deriv_thermo_base_tl::compress_const_tptr()`
- :cpp:func:`o2scl::deriv_thermo_base_tl::coeff_thermal_exp()`
- :cpp:func:`o2scl::deriv_thermo_base_tl::squared_sound_speed()`

..
   (begin comment)
   I think the expression below only works for fermions? 
   I'm taking this section out until it's better commented
    
   In the case where the particle is interacting, the 
   derivative of the density with respect to the effective mass is
   \f[
   \left(\frac{dn}{dm^{*}}\right)_{\mu,T} = 
   \left(\frac{3 n}{m^{*}}\right) - 
   \frac{T}{m^{*}} \left(\frac{dn}{dT}\right)_{m^{*},\mu} -
   \frac{\nu}{m^{*}} \left(\frac{dn}{d\mu}\right)_{m^{*},T} 
   \f]
   This relation holds whether or not the mass is included in the
   chemical potential :math:`\nu`, as the rest mass is held
   constant even though the effective mass is varying. This
   relation also holds in the case where the particle is
   non-interacting, so long as one does not allow the rest mass in
   the chemical potential to vary. This derivative is useful, for
   example, in models of quark matter where the quark mass is
   dynamically generated.
   (end comment)

Particle example
----------------

.. literalinclude:: ../../../../examples/ex_part.cpp
   :language: c++		    
   :start-after: sphinx-example-start
