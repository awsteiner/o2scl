Moment of Inertia in the Slowly-Rotating Approximation
======================================================

The differential equations for slow rigid rotation are solved
by \ref o2scl::tov_solve if \ref o2scl::tov_solve::ang_vel is
set to <tt>true</tt>. 

In the case of slow rigid rotation with angular velocity
\f$ \Omega \f$, the moment of inertia is
\f[
I = \frac{8 \pi}{3} \int_0^R dr~r^4\left(\varepsilon+P\right)
\left(\frac{\bar{\omega}}{\Omega}\right)
e^{\Lambda-\Phi}
= \frac{8 \pi}{3} \int_0^R dr~r^4\left(\varepsilon+P\right)
\left(\frac{\bar{\omega}}{\Omega}\right)
\left(1-\frac{2 G m}{r}\right)^{-1/2} e^{-\Phi} 
\f]
where \f$ \omega(r) \f$ is the rotation rate of the inertial
frame, \f$ \Omega \f$ is the angular velocity in the fluid
frame, and \f$ \bar{\omega}(r) \equiv \Omega - \omega(r) \f$ 
is the angular velocity of a fluid element at infinity.
The function \f$ \bar{\omega}(r) \f$ is the solution of 
\f[
\frac{d}{dr} \left( r^4 j \frac{d \bar{\omega}}{dr}\right)
+ 4 r^3 \frac{d j}{dr} \bar{\omega} = 0
\f]
where the function \f$ j(r) \f$ is defined by 
\f[
j = e^{-\Lambda-\Phi} =
\left( 1-\frac{2 G m}{r} \right)^{1/2} e^{-\Phi} \, .
\f]
Note that \f$ j(r=R) = 1 \f$. 
The boundary conditions for \f$ \bar{\omega} \f$ are
\f$ d \bar{\omega}/dr = 0 \f$ at \f$ r=0 \f$ and 
\f[
\bar{\omega}(R) = \Omega - \left(\frac{R}{3}\right)
\left(\frac{d \bar{\omega}}{dr}\right)_{r=R} \, .
\f]
One can use the TOV equation to rewrite the moment of 
inertia as 
\f[
I= \left(\frac{d \bar{\omega}}{dr}\right)_{r=R} 
\frac{R^4}{6 G \Omega} \, .
\f]
The star's angular momentum is just \f$ J = I \Omega \f$.

