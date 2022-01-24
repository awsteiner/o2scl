Solution of the Tolman-Oppenheimer-Volkov Equations
===================================================

:ref:`O2scl <o2scl>`
     
The class :ref:`tov_solve <tov_solve>` provides a solution to the
Tolman-Oppenheimer-Volkov (TOV) equations given an equation of state
(EOS), provided as an object of type :ref:`eos_tov <eos_tov>`. These
classes are particularly useful for static neutron star structure:
given any equation of state one can calculate the mass vs. radius
curve and the properties of any star of a given mass.

The EOS is typically specified using :ref:`eos_tov_interp <eos_tov_interp>`
which uses linear interpolation to interpolate a user-specified
:ref:`table <o2scl:table>` object. The Buchdahl EOS is given in 
:ref:`eos_tov_buchdahl <eos_tov_buchdahl>`, a single polytrope EOS is given
in :ref:`eos_tov_polytrope <eos_tov_polytrope>`, and a linear EOS is given
in :ref:`eos_tov_linear <eos_tov_linear>`.

In units where :math:`c=1`, the most general static and 
spherically symmetric metric is of the form

.. math::

   ds^2 = - e^{2 \Phi(r)} d t^2 + e^{2 \Lambda(r)} d r^2 + 
   r^2 d \theta^2 + r^2 \sin^2 \theta~d \phi^2

where :math:`\theta` is the polar angle and :math:`\phi`
is the azimuthal angle. Often we will not write explicitly
the radial dependence for many of the quantities defined
below, i.e. :math:`\Phi \equiv \Phi(r)`.

This leads to the TOV equation (i.e. Einstein's
equations for a static and spherically symmetric star)

.. math::

  \frac{d P}{d r} = - \frac{G \varepsilon m}{r^2}
  \left( 1+\frac{P}{\varepsilon}\right)
  \left( 1+\frac{4 \pi P r^3}{m} \right)
  \left( 1-\frac{2 G m}{r}\right)^{-1}

where :math:`r` is the radial coordinate, :math:`m(r)` is the
gravitational mass enclosed within a radius :math:`r`, and
:math:`\varepsilon(r)` and :math:`P(r)` are the energy density and
pressure at :math:`r`, and :math:`G` is the gravitational constant.
The mass enclosed, :math:`m(r)`, is related to the energy density
through

.. math::
   
   \frac{d m}{d r} = 4 \pi r^2 \varepsilon

and these two differential equations can be solved simultaneously
given an equation of state, :math:`P(\varepsilon)`.
The total gravitational mass is given by

.. math::
   
   M = \int_0^R 4 \pi r^2 \varepsilon d r

The boundary conditions are :math:`m(r=0)=0` and the condition
:math:`P(r=R)=0` for some fixed radius :math:`R`. These boundary
conditions give a one-dimensional family solutions to the TOV
equations as a function of the central pressure. Each central
pressure implies a gravitational mass, :math:`M`, and radius,
:math:`R`, and thus defines a mass-radius curve.

The metric function :math:`\Lambda` is

.. math::

   e^{2 \Lambda} = \left( 1-\frac{2 G m}{r}\right)^{-1}

The other metric function, :math:`\Phi(r)` is sometimes referred
to as the gravitational potential. In vacuum above the star, it is

.. math::

   e^{2 \Phi} = \left( 1-\frac{2 G M}{r}\right)

and inside the star it is determined by

.. math::
   
   \frac{d \Phi}{d r} = - \frac{1}{\varepsilon}
   \frac{ d P}{d r} \left(1+\frac{P}{\varepsilon}\right)^{-1} =
   \frac{G m}{r^2} \left( 1+\frac{4 \pi P r^3}{m} \right)
   \left( 1-\frac{2 G m}{r}\right)^{-1}

Alternatively, that this can be rewritten as

.. math::

   -d \Phi = \frac{d P}{P+\varepsilon} \, .

In this form, :math:`\Phi` has no explicit dependence on :math:`r`
so it can be computed (up to a constant) directly from the 
EOS.

If the neutron star is at zero temperature and there is
only one conserved charge, (i.e. baryon number), then

.. math::
   
   -d \Phi = \frac{d P}{\mu n} = \frac{d \mu}{\mu}

and this implies that :math:`\mu e^{\Phi}` is everywhere
constant in the star. If one defines the
"enthalpy" by

.. math::
   
   d h = \frac{dP}{P + \varepsilon} 

then

.. math::
   
   -d \Phi = dh

and thus :math:`\mu \propto e^{h}` or :math:`h = \ln \mu + C`.
This is the enthalpy used by the :ref:`nstar_rot <nstar_rot>` class.

Keep in mind that this enthalpy is determined by integrating
the quantities in the stellar profile (which may be, for example,
in beta-equilibrium). Thus, 
this is not equal the usual thermodynamic enthalpy which is

.. math::
   
   H(P,S,N) = E + P V = T S + \sum_i \mu_i N_i

or in differential form

.. math::

   d H = T dS + V dP + \sum_i \mu_i d N_i \, .

The proper boundary condition for the gravitational potential
is

.. math::
   
   \Phi(r=R) = \frac{1}{2} \ln \left(1-\frac{2 G M}{R} \right)

which ensures that :math:`\Phi(r)` matches the metric
above in vacuum. Since the expression for :math:`d\Phi/dr`
is independent of :math:`\Phi`, the differential equation
can be solved for an arbitrary value of :math:`\Phi(r=0)`
and then shifted afterwards to obtain the correct
boundary condition at :math:`r=R` .

The surface gravity is defined to be

.. math::

   g = \frac{G m}{r^2}\left(1-\frac{2 G m}{r}\right)^{-1/2}

which is computed in units of inverse
kilometers, and the redshift is defined to be

.. math::

   z = \left(1-\frac{2 G m}{r}\right)^{-1/2} - 1

which is unitless.

The baryonic mass is typically defined by

.. math::
   
   M_B = \int_0^R 4 \pi r^2 n_B m_B 
   \left(1-\frac{2 G m}{r}\right)^{-1/2} d r

where :math:`n_B(r)` is the baryon number density at radius :math:`r`
and :math:`m_B` is the mass one baryon (taken to be the mass of the
proton by default and stored in
:cpp:var:`o2scl::tov_solve::baryon_mass`). If the EOS specifies the
baryon density (i.e. if :cpp:var:`o2scl::eos_tov::baryon_column` is
true), then :ref:`tov_solve <tov_solve>` will compute the associated
baryonic mass for you.
