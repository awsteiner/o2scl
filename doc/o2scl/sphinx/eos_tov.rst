Equations of State for the TOV equations
========================================

:ref:`O2scl <o2scl>`

Background
----------

In the simplest models, neutron stars consist only of neutrons,
protons, and electrons at :math:`T=0`. The thermodynamic identity
is

.. math::

   P + \varepsilon = n_n \mu_n + n_p \mu_p + n_e \mu_e

where :math:`P` is the pressure, :math:`\varepsilon` is the energy
density, :math:`n_i` is the number density of particle :math:`i`, and
:math:`\mu_i` is the chemical potential of particle :math:`i`.
Presuming charge neutrality, :math:`n_p=n_e` and beta-equilibrium,
:math:`\mu_n=\mu_p+\mu_e`, this simplifies to
   
.. math::

   P + \varepsilon = n_B \mu_B

where the baryon density is :math:`n_B \equiv n_n+n_p` and the baryon
chemical potential is :math:`\mu_B \equiv \mu_n`. Thus, in this simple
model, the thermodynamic properties of neutron star matter behave as
if matter only has one component. Given the relationship between
energy density and baryon density, one can obtain the baryon chemical
potential and the pressure directly

.. math::

   \mu_B = \frac{\partial \varepsilon}{\partial n_B} \qquad \mathrm{and}
   \qquad P = n_B^2 \frac{\partial (\varepsilon/n_B)}{\partial n_B}

The speed of sound is

.. math::

   c_s^2 = c^2 \left( \frac{\partial P}
   {\partial \varepsilon}\right)_{\mathrm{comp}}

where :math:`\mathrm{comp}` indicates that the derivative is to be
performed at fixed composition (e.g. a fixed ratio between the
neutron, proton, and electron densities). 
     
   
Class infrastructure
--------------------
     
The TOV solver requires the EOS to be specified as an object of type
:ref:`eos_tov <eos_tov>`. The documentation of this parent class
contains more information. The class :ref:`eos_tov_interp
<eos_tov_interp>` is used most frequently. It uses linear
interpolation to interpolate a user-specified :ref:`table <table>`
object. A faster lower-level EOS interpolation is performed by
:ref:`eos_tov_vectors <eos_tov_vectors>`. The Buchdahl EOS is given in
:ref:`eos_tov_buchdahl <eos_tov_buchdahl>`, a single polytrope EOS is
given in :ref:`eos_tov_polytrope <eos_tov_polytrope>`, a linear EOS is
given in :ref:`eos_tov_linear <eos_tov_linear>`, and an EOS with a
polynomial form for the speed of sound is given in :ref:`eos_cs2_poly
<eos_cs2_poly>`

From pressure and energy density
--------------------------------
     
Given a relationship between pressure and energy density in
beta-equilibrium, one can obtain the baryon density and the baryon
chemical potential up to a constant. To see this, start with the
thermodynamic identity

.. math::

   P + \varepsilon = \mu_B n_B = n_B \left(\frac{\partial
   \varepsilon}{\partial n_B}\right)

Then, expressing the pressure in terms of the energy density

.. math::

   \frac{d \varepsilon}{P(\varepsilon)+\varepsilon} =
   \frac{d n_B}{n_B}

Now we integrate, to get

.. math::

   \int \frac{d \varepsilon}{P(\varepsilon)+\varepsilon} =
   \ln n_B + C

where :math:`C` is an undetermined constant. If we create a
new integration variable for the energy density and
presume that :math:`\varepsilon_1 \equiv \varepsilon(n_B = n_{B,1})`,
then we can fix :math:`C`

.. math::

   \int_{\varepsilon_1}^{\varepsilon}
   \frac{d \tilde{\varepsilon}}{P(\tilde{\varepsilon})+\tilde{\varepsilon}} =
   \ln n_B - \ln n_{B,1}

Going back to the original thermodynamic identity, we can also
write

.. math::

   P + \varepsilon = \mu_B \left( \frac{\partial
   P}{\partial \mu_B}\right)

to obtain (similar to the method above)

.. math::

   \int \frac{d P}{\varepsilon(P)+P} =
   \ln \mu_B + C

and thus if :math:`P_1=P(\mu_B=\mu_{B,1})`:

.. math::

   \int_{P_1}^{P}
   \frac{d \tilde{P}}{\varepsilon(\tilde{P})+\tilde{P}} =
   \ln \mu_B - \ln \mu_{B,1}

At low densities, if one assumes the low-density equation of state
made of an ideal gas of nuclei, then at :math:`P=0`, :math:`\mu_{B}` is
about 931 MeV (approximately the atomic mass unit, :math:`m_u`). Then
we get

.. math::

   \int_{0}^{P}
   \frac{d \tilde{P}}{\varepsilon(\tilde{P})+\tilde{P}} =
   \ln \mu_B - \ln m_u


   
   
