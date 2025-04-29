Moment of Inertia in the Slowly-Rotating Approximation
======================================================

:ref:`O2scl <o2scl>`
     
The differential equations for slow rigid rotation are solved by
:ref:`tov_solve <tov_solve>` if :cpp:var:`o2scl::tov_solve::ang_vel`
is set to ``true``.

In the case of slow rigid rotation with angular velocity
:math:`\Omega`, the moment of inertia is
      
.. math::
   
   I = \frac{8 \pi}{3} \int_0^R dr~r^4\left(\varepsilon+P\right)
   \left(\frac{\bar{\omega}}{\Omega}\right)
   e^{\Lambda-\Phi}
   = \frac{8 \pi}{3} \int_0^R dr~r^4\left(\varepsilon+P\right)
   \left(\frac{\bar{\omega}}{\Omega}\right)
   \left(1-\frac{2 G m}{r}\right)^{-1/2} e^{-\Phi} 

where :math:`\omega(r)` is the rotation rate of the inertial
frame, :math:`\Omega` is the angular velocity in the fluid
frame, and :math:`\bar{\omega}(r) \equiv \Omega - \omega(r)` 
is the angular velocity of a fluid element at infinity.
The function :math:`\bar{\omega}(r)` is the solution of

.. math::
   
   \frac{d}{dr} \left( r^4 j \frac{d \bar{\omega}}{dr}\right)
   + 4 r^3 \frac{d j}{dr} \bar{\omega} = 0

where the function :math:`j(r)` is defined by

.. math::
   
   j = e^{-\Lambda-\Phi} =
   \left( 1-\frac{2 G m}{r} \right)^{1/2} e^{-\Phi} \, .

Note that :math:`j(r=R) = 1`. 
The boundary conditions for :math:`\bar{\omega}` are
:math:`d \bar{\omega}/dr = 0` at :math:`r=0` and
      
.. math::
   
   \bar{\omega}(R) = \Omega - \left(\frac{R}{3}\right)
   \left(\frac{d \bar{\omega}}{dr}\right)_{r=R} \, .

One can use the TOV equation to rewrite the moment of 
inertia as

.. math::
   
   I= \left(\frac{d \bar{\omega}}{dr}\right)_{r=R} 
   \frac{R^4}{6 G \Omega} \, .

The star's angular momentum is just :math:`J = I \Omega`.

