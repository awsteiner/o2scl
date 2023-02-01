Equations of State for the TOV equations
========================================

:ref:`O2scl <o2scl>`

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


   
   
