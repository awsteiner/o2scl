Fermions
========

:ref:`O2scl_part <o2sclp>`

Relativistic and Non-relativistic Fermions
------------------------------------------

There are a few distinctions between how relativistic and
nonrelativistic fermions are handled in :ref:`o2scl_part <o2sclp>`
which are worth noting. For an interacting relativistic fermion, the
effective mass, :math:`m^{*}`, and the effective chemical potential,
:math:`\nu` are defined so that the energy density is

.. math::

   {\varepsilon}_{\mathrm{R}} = \frac{g}{2 \pi^2} \int
   dk~\frac{k^2 \sqrt{k^2+m^{* 2}}}
   { 1+\exp\left[\left(\sqrt{k^2+m^{*2}}-
   \nu_{\mathrm{R}}\right)/T\right]}
   = \frac{g}{2 \pi^2} \int
   dk~\frac{k^2 \sqrt{k^2+m^{* 2}} }
   {1+\exp\left[\left(\sqrt{k^2+m^{*2}}-
   \bar{\nu}_{\mathrm{R}}-m\right)/T\right]}

for a relativistic fermion, where we define the chemical potential
without the rest mass with :math:`\bar{\nu}_{\mathrm{R}} \equiv
\nu_{\mathrm{R}}-m`. When :cpp:var:`o2scl::part_tl::inc_rest_mass` is
true, :cpp:var:`o2scl::part_tl::nu` is equal to :math:`\nu_{\mathrm{R}}`
and when :cpp:var:`o2scl::part_tl::inc_rest_mass` is false,
:cpp:var:`o2scl::part_tl::nu` is equal to
:math:`\bar{\nu}_{\mathrm{R}}` . If we define :math:`\psi_{\mathrm{R}}
= (\nu_{\mathrm{R}}-m^{*})/T`, :math:`\phi = m^{*}/T`, :math:`u \equiv
k/T` then the energy density is

.. math::

   {\varepsilon}_{\mathrm{R}} = \frac{g T^4}{2 \pi^2} \int
   du~\frac{u^2 \sqrt{u^2+\phi^2}}
   { 1+\exp\left[\sqrt{u^2+\phi^2} - \psi_{\mathrm{R}} - \phi \right]}

This expression is used for
:cpp:var:`o2scl::part_calibrate_class::part_calibrate()` because
:math:`\varepsilon_{\mathrm{R}}/(g T^4)` depends only on
:math:`\psi_{\mathrm{R}}` and :math:`\phi`. For a nonrelativistic
fermion,

.. math::

   \bar{\varepsilon}_{\mathrm{NR}} = 
   \frac{g}{2 \pi^2} \int dk~
   \frac{k^4}{2 m^{*}}
   \left\{ 1+\exp\left[\left(\frac{k^2}{2 m^{*}}-
   \bar{\nu}_{\mathrm{NR}}\right)/T\right] \right\}^{-1}
   = \frac{g}{2 \pi^2} \int dk~
   \frac{k^4}{2 m^{*}} 
   \left\{ 1+\exp\left[\left(\frac{k^2}{2 m^{*}}-
   \nu_{\mathrm{NR}}+m\right)/T\right] \right\}^{-1}

where :math:`\bar{\nu}_{\mathrm{NR}} = \nu_{\mathrm{NR}} - m` . Note
that the rest mass energy density is
:math:`\varepsilon_{\mathrm{rest}} = n m` (not :math:`n m^{*}`) in
both cases, but it is included in :math:`\varepsilon_{\mathrm{R}}`
while it is not included in :math:`\bar{\varepsilon}_{\mathrm{NR}}` .
Taking the nonrelativsitic limit of the relativistic energy density
shows that :math:`\nu_{\mathrm{R}} - m^{*} = \bar{\nu}_{\mathrm{NR}}`.
Thus the class :ref:`fermion_nonrel_tl <fermion_nonrel_tl>` uses the
value stored in :cpp:var:`o2scl::part_tl::nu` slightly differently
than does :ref:`fermion_rel_tl <fermion_rel_tl>` and :ref:`fermion_eff
<fermion_eff>` . The Fermi momentum is also handled slightly
differently, :math:`k_{F,\mathrm{R}} \equiv
\sqrt{\nu_{\mathrm{R}}^2-m^{* 2}}` and :math:`k_{F,\mathrm{NR}} \equiv
\sqrt{2 \bar{\nu}_{\mathrm{NR}} m^{*}}`.

Now if we define :math:`u_{\mathrm{NR}} \equiv k^2/(2 m^{*} T)` 
and :math:`\psi_{\mathrm{NR}} \equiv (\nu_{\mathrm{NR}}-m^{*})/T`
then the argument of the exponential is 

.. math::

   \frac{k^2}{2 m^{*} T } - \frac{\bar{\nu}_{\mathrm{NR}}}{T} = 
   u_{\mathrm{NR}} - \psi_{\mathrm{NR}} + \frac{m}{T}- \phi

which is inconvenient because then :math:`\varepsilon_{\mathrm{NR}}/(g
T^4)` is no longer a function of :math:`\psi_{\mathrm{NR}}` and
:math:`\phi` alone. Thus we define :math:`\psi_{\mathrm{NR}} \equiv
\bar{\nu}_{\mathrm{NR}}/T` and then the energy density is

.. math::

   \bar{\varepsilon}_{\mathrm{NR}} = \frac{g T^4}{2 \pi^2} \int
   du_{\mathrm{NR}}~\frac{\sqrt{2}~u_{\mathrm{NR}}^{3/2} \phi^{3/2}}
   { 1+\exp\left[u_{\mathrm{NR}} - \psi_{\mathrm{NR}} \right]}

which is now a function of :math:`\psi_{\mathrm{NR}}` and
:math:`\phi`alone. This is the form used to compute the energy density
in :ref:`fermion_nonrel_tl <fermion_nonrel_tl>` and the definition of
:math:`\psi_{\mathrm{NR}}` used for nonrelativistic fermions in \ref
:cpp:func:`o2scl::part_calibrate_class::part_calibrate()`.

Fermion integrations
--------------------

(Slowly moving the documentation from fermion_rel
and fermion_deriv_rel to this location.)

In many cases, the non-interacting expressions for fermion
thermodynamics can be used in interacting systems as long as one
replaces the mass with an effective mass, :math:`m^{*}` and the
chemical potential with an effective chemical potential, :math:`\nu` .
In the case where :math:`\nu` includes the rest mass (still denoted
:math:`m`), the fermionic distribution function is

.. math::

   f = \frac{1}{1+e^{(\sqrt{k^2+m^{* 2}}-\nu)/T}}
   \quad ; \quad
   f = \frac{1}{1+e^{(\sqrt{k^2+m^{* 2}}-\nu-m)/T}}

where the left expression is used when the chemical potential includes
the rest mass and the energy density includes the rest mass energy
density, (:cpp:var:`o2scl::part_tl::inc_rest_mass` is ``true``) and
the right expression is used when the rest mass is not included
(:cpp:var:`o2scl::part_tl::inc_rest_mass` is ``false``). For
convenience, we define :math:`E^{*} \equiv \sqrt{k^2+m^{* 2}}`.

Upper limits
------------

The fermionic integrands vanish when the argument of
the exponential becomes large compared to a positive
number :math:`\zeta`.
This condition is

.. math::

   \sqrt{k^2+m^{* 2}}-\nu \gg \zeta T \quad ; \quad
   \sqrt{k^2+m^{* 2}}-\nu-m \gg \zeta T

Thus solving
for the momentum, an upper limit, :math:`k_{\mathrm{ul}}` is

.. math::

   k_{\mathrm{ul}} = \sqrt{\left(\zeta T + \nu\right)^2-m^{* 2}}
   \quad ; \quad
   k_{\mathrm{ul}} = \sqrt{\left(\zeta T + m + \nu\right)^2-m^{* 2}}
    
The entropy is only significant at the Fermi surface, thus
in the degenerate case, the lower limit of the entropy
integral can be given be determined by the value of :math:`k` 
which solves

.. math::

   - \zeta = \frac{\sqrt{k^2+m^{* 2}}-\nu}{T} 
   \quad ; \quad
   - \zeta = \frac{\sqrt{k^2+m^{* 2}}-\nu-m}{T} 

The solution is 

.. math::

   k_{\mathrm{ll}} = \sqrt{(-\zeta T+{\nu})^2-m^{*,2}}
   \quad ; \quad
   k_{\mathrm{ll}} = \sqrt{(-\zeta T + m +\nu)^2-m^{*,2}}

which is a valid lower limit only if the argument under
the square root is positive.

Integrands
----------
    
The energy density is

.. math::

   \varepsilon = \frac{g}{2 \pi^2} \int_0^{\infty} 
   k^2~dk~\sqrt{k^2+m^{* 2}} f 
   \quad ; \quad
   \varepsilon = \frac{g}{2 \pi^2} \int_0^{\infty} 
   k^2~dk~\left(\sqrt{k^2+m^{* 2}}-m\right) f \, ,

the number density is

.. math::

   n = \frac{g}{2 \pi^2} \int_0^{\infty} 
   k^2~dk~f \, ,

and the entropy density is

.. math::

   s = \frac{g}{2 \pi^2} \int_0^{\infty} 
   dk~(-k^2 {\cal S})

where

.. math::

   {\cal S}\equiv f \ln f +(1-f) \ln (1-f)
   \quad ; \quad
   \frac{\partial {\cal S}}{\partial f} = \ln 
   \left(\frac{f}{1-f}\right) \, .

The derivative can also be written

.. math::

   \frac{\partial {\cal S}}{\partial f} = 
   \left(\frac{\nu-E^{*}}{T}\right)
   \quad ; \quad
   \frac{\partial {\cal S}}{\partial f} = 
   \left(\frac{\nu-E^{*}+m}{T}\right)
    
In the degenerate regime, :math:`{\cal S}`, can lose precision when
:math:`(E^{*} - \nu)/T` is negative and sufficiently large in absolute
magnitude. Thus when :math:`(E^{*} - \nu)/T < \xi` (for :math:`\xi
\rightarrow - \infty` ) an alternative expression

.. math::

   {\cal S} \approx 
   e^{(E^{*}-\nu)/T}
   \left( \frac{E^{*} -\nu-T}{T} \right)
   \quad ; \quad
   {\cal S} \approx 
   e^{(E^{*}-\nu-m)/T}
   \left( \frac{E^{*} -\nu-m-T}{T} \right)
   \, 

can be used.
    
Non-degenerate integrands
-------------------------
    
..
   (begin comment)
   It's not at all clear that this dimensionless form is more
   accurate than other potential alternatives. On the other hand,
   it seems that the uncertainties in the integrations are larger
   than the errors made by the integrand at present.
   (end comment)

The integrands in the non-degenerate regime are written in a
dimensionless form, by defining :math:`u=(E^{*}-m^{*})/T` (this choice
ensures :math:`k=0` corresponds to :math:`u=0`), :math:`y \equiv \nu/
T` (or :math:`y = (\nu+m)/T` if the chemical potential does not
include the mass), and :math:`\eta \equiv m^{*}/T`. Then :math:`k/T =
\sqrt{u^2+2 u \eta}`, :math:`(1/T) dk = E^{*}/k du =
(u+\eta)/\sqrt{u^2+2 u \eta}~du`, and :math:`f = 1/(1+e^{u+\eta-y})` .
The density is

.. math::

   n = \frac{g T^3}{2 \pi^2} \int_0^{\infty}~du~
   \sqrt{u^2+2 u \eta} (u+\eta)
   \left(1+e^{u+\eta-y}\right)^{-1}

the energy density is 

.. math::

   \varepsilon = \frac{g T^4}{2 \pi^2} \int_0^{\infty}~du~
   \sqrt{u^2+2 u \eta} (u+\eta)^2
   \left(1+e^{u+\eta-y}\right)^{-1}

and the entropy density is 

.. math::

   s = -\frac{g T^3}{2 \pi^2} \int_0^{\infty}~du~
   \sqrt{u^2+2 u \eta} (u+\eta) {\cal S}
    
Evaluation of the derivatives
-----------------------------
    
The relevant
derivatives of the distribution function are

.. math::

   \frac{\partial f}{\partial T}=
   f(1-f)\frac{E^{*}-\nu}{T^2}
   \quad ; \quad
   \frac{\partial f}{\partial T}=
   f(1-f)\frac{E^{*}-m-\nu}{T^2}

.. math::

   \frac{\partial f}{\partial \nu}=
   f(1-f)\frac{1}{T}

.. math::
   
   \frac{\partial f}{\partial k}=
   -f(1-f)\frac{k}{E^{*} T}
   
.. math::

   \frac{\partial f}{\partial m^{*}}=
   -f(1-f)\frac{m^{*}}{E^{*} T}
    
The derivatives can be integrated directly direct) or they may be
converted to integrals over the distribution function through an
integration by parts

.. math::

   \int_a^b f(k) \frac{d g(k)}{dk} dk = \left.f(k) g(k)\right|_{k=a}^{k=b}
   - \int_a^b g(k) \frac{d f(k)}{dk} dk 

using the distribution function for :math:`f(k)` and 0 and 
:math:`\infty` as the limits, we have

.. math::

   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{d g(k)}{dk} f dk =
   \frac{g}{2 \pi^2} \int_0^{\infty} g(k) f (1-f) \frac{k}{E^{*} T} dk 

as long as :math:`g(k)` vanishes at :math:`k=0` .
Rewriting using :math:`g(k) = h(k) E^{*} T/k` 

.. math::

   \frac{g}{2 \pi^2} \int_0^{\infty} h(k) f (1-f) dk =
   \frac{g}{2 \pi^2} \int_0^{\infty} f \frac{T}{k} 
   \left[ h^{\prime} E^{*}-\frac{h E^{*}}{k}+\frac{h k}{E^{*}} \right] dk

as long as :math:`h(k)/k` vanishes at :math:`k=0` .
    
Explicit forms
--------------
    
1) The derivative of the density wrt the chemical potential

.. math::

   \left(\frac{d n}{d \mu}\right)_T = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2}{T} f (1-f) dk

Using :math:`h(k)=k^2/T` we get

.. math::

   \left(\frac{d n}{d \mu}\right)_T = 
   \frac{g}{2 \pi^2} \int_0^{\infty} 
   \left(\frac{k^2+E^{*2}}{E^{*}}\right) f dk
    
2) The derivative of the density wrt the temperature

.. math::

   \left(\frac{d n}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-\nu)}{T^2} 
   f (1-f) dk
   \quad ; \quad
   \left(\frac{d n}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{k^2(E^{*}-m-\nu)}{T^2} 
   f (1-f) dk

Using :math:`h(k)=k^2(E^{*}-\nu)/T^2` we get

.. math::

   \left(\frac{d n}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
   \left[2 k^2+E^{*2}-E^{*} \nu -
   k^2 \left(\frac{\nu}{E^{*}}\right)\right] dk
   \quad ; \quad
   \left(\frac{d n}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f}{T} 
   \left[2 k^2+E^{*2}-E^{*}\left(\nu+m\right)-
   k^2 \left(\frac{\nu+m}{E^{*}}\right)\right] dk
    
3) The derivative of the entropy wrt the chemical potential

.. math::

   \left(\frac{d s}{d \mu}\right)_T = 
   \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
   \frac{(E^{*}-\nu)}{T^2} dk
   \quad ; \quad
   \left(\frac{d s}{d \mu}\right)_T = 
   \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
   \frac{(E^{*}-m-\nu)}{T^2} dk

This verifies the Maxwell relation

.. math::

   \left(\frac{d s}{d \mu}\right)_T =
   \left(\frac{d n}{d T}\right)_{\mu}
    
4) The derivative of the entropy wrt the temperature

.. math::

   \left(\frac{d s}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
   \frac{(E^{*}-\nu)^2}{T^3} dk
   \quad ; \quad
   \left(\frac{d s}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} k^2 f (1-f) 
   \frac{(E^{*}-m-\nu)^2}{T^3} dk

Using :math:`h(k)=k^2 (E^{*}-\nu)^2/T^3` 

.. math::

   \left(\frac{d s}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f(E^{*}-\nu)}{E^{*}T^2} 
   \left[E^{* 3}+3 E^{*} k^2- (E^{* 2}+k^2)\nu\right] d k
   \quad ; \quad
   \left(\frac{d s}{d T}\right)_{\mu} = 
   \frac{g}{2 \pi^2} \int_0^{\infty} \frac{f(E^{*}-m-\nu)}{E^{*}T^2} 
   \left[E^{* 3}+3 E^{*} k^2- (E^{* 2}+k^2)(\nu+m)\right] d k
    
5) The derivative of the density wrt the effective mass

.. math::

   \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
   -\frac{g}{2 \pi^2} \int_0^{\infty} 
   \frac{k^2 m^{*}}{E^{*} T} f (1-f) dk

Using :math:`h(k)=-(k^2 m^{*})/(E^{*} T)` we get

.. math::

   \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = 
   -\frac{g}{2 \pi^2} \int_0^{\infty} 
   m^{*} f dk

..
   (begin comment)
   This derivative may be written in terms of the 
   others
   \f[
   \left(\frac{d n}{d m^{*}}\right)_{T,\mu} = \frac{3 n}{m^{*}}
   - \frac{T}{m^{*}}\left[ \left(\frac{d n}{d T}\right)_{\mu}
   +\frac{\mu}{T} \left(\frac{d n}{d \mu}\right)_{T}
   \right] - \left(\frac{d n}{d \mu}\right)_{T}
   \f]
   (end comment)
    
Expansions for Fermions
-----------------------

Presuming the chemical potential includes the rest mass,
and :math:`E=\sqrt{k^2+m^2}`,
the pressure for non-interacting fermions with degeneracy :math:`g` is

.. math::

   P = \frac{g T}{2 \pi^2} \int_0^{\infty} 
   k^2~dk~\ln \left[ 1 + e^{-(E-\mu)/T}\right] = 
   \frac{g}{2 \pi^2} \int_0^{\infty} k^2\left(\frac{k^2}{3 E}\right)~dk~
   \frac{1}{1 + e^{(E-\mu)/T}} \, ,

where the second form is obtained with an integration by parts. We use
units where :math:`\hbar=c=1`. The variable substitutions from \ref
Johns96 are :math:`\ell = k/m`, :math:`\psi = (\mu-m)/T`, and
:math:`t=T/m`. (Presumably this choice of variables gives better
results for non-relativistic fermions because the mass is separated
from the chemical potential in the definition of :math:`\psi`, but I
haven't checked this.) These replacements give

.. math::

   P = \frac{g m^4}{2 \pi^2} 
   \int_0^{\infty} d\ell~\frac{\ell^4}{3 \sqrt{\ell^2+1}}
   \left( \frac{1}{1 + e^{z/t-\psi}} \right)

where :math:`z = \sqrt{\ell^2+1}-1` . 
Re-expressing in terms of :math:`z`, one obtains

.. math::

   \frac{\ell^4}{3 \sqrt{\ell^2+1}} = \frac{z^2(2+z)^2}
   {3 (1+z)} \quad\mathrm{and}\quad 
   \frac{d \ell}{d z} = \frac{1+z}{\sqrt{z(2+z)}} \, .

The pressure is

.. math::

   P = \frac{g m^4}{2 \pi^2} 
   \int_0^{\infty} dz~\frac{1}{3}[z(2+z)]^{3/2}
   \left[ \frac{1}{1 + e^{(z-x)/t}} \right] \, .

where :math:`x = \psi t = (\mu-m)/m`. 

Degenerate expansion
--------------------

The Sommerfeld expansion for :math:`t \rightarrow 0` is

.. math::

   \begin{eqnarray}
   \int_0^{\infty} dz~\frac{f(z)}{1 + e^{(z-x)/t}} &=&
   \int_0^{x} f(z) + \frac{\pi^2 t^2}{6} f^{\prime}(x) +
   \frac{7 \pi^4 t^4}{360} f^{(3)}(x) +
   \frac{31 \pi^6 t^6}{15120} f^{(5)}(x) + \ldots \nonumber \\
   &=& \int_0^{x} f(z) + \sum_{n=1}^{\infty}
   \pi^{2n}t^{2n} \left[f^{(2n -1)}(x) \right] 
   \left[ \frac{2 (-1)^{1+n}(2^{2n-1}-1)B_{2n}}{(2n)!} \right] \nonumber
   \end{eqnarray}
   
This is an asymptotic expansion, and must thus be used with care.
Define :math:`\tilde{P}(x,t) \equiv 2 \pi^2 P/(g m^4)`. The first term
in the Sommerfeld expansion for :math:`\tilde{P}` depends only on
:math:`x` alone:

.. math::

   P_0 \equiv \frac{1}{24} (1+x)\sqrt{x(2+x)} \left[ -3 + 2 x(2+x)\right]
   + \frac{1}{4} \log \left[ \frac{
   \sqrt{x}+\sqrt{2+x}}{\sqrt{2}} \right]

where :math:`x = \psi t` . This expression cannot be used when
:math:`x` is small, but a Taylor series expansion can be used
instead. A few terms are

.. math::

   \frac{2 \pi^2 P}{g m^4} = P_0 + \frac{\pi^2 t^2}{6} \sqrt{x(2+x)}(1 + x) +
   \frac{7 \pi^4 t^4}{360} \left\{\frac{(1+x)(2
   x^2+4x-1)}{[x(2+x)]^{3/2}} \right\}
   -\frac{31\pi^6 t^6}{1008} \frac{(1+x)\sqrt{x(2+x)}}{x^4 (2+x)^4} + 
   \ldots

The number density is

.. math::

   n = \frac{dP}{d \mu} = \frac{d P}{d x} \frac{d x}{d \mu} = 
   \frac{1}{m} \left(\frac{d P}{d x}\right)_t

Note that because the density is a derivative, it is possible
that the terms in the density fail before the terms in the 
pressure, thus we should use one less term for the density
when using the expansion. The entropy is

.. math::

   s = \frac{dP}{d T} = \frac{d P}{d t} \frac{d t}{d T} = 
   \frac{1}{m} \left(\frac{d P}{d t}\right)_x

The derivative of the number density with respect to the 
chemical potential is

.. math::

   \frac{d n}{d \mu} = \frac{d^2P}{d \mu^2} = \frac{d}{d \mu}
   \left(\frac{d P}{d x} \frac{d x}{d \mu}\right) = 
   \frac{d^2 P}{d x^2} \left(\frac{d x}{d \mu}\right)^2 +
   \frac{d P}{d x} \frac{d^2 x}{d \mu^2} = 
   \frac{1}{m^2} \left(\frac{d^2 P}{d x^2}\right)_t \, .

The derivative of the number density with respect to the
temperature is

.. math::

   \frac{d n}{d T} = \frac{d^2P}{d \mu dT} = 
   \frac{1}{m^2} \frac{d^2 P}{d x d t} \, ,

and the derivative of the entropy density with respect to 
the temperature is

.. math::

   \frac{d s}{d T} = \frac{d^2P}{d T^2} = 
   \frac{1}{m^2} \left(\frac{d^2 P}{d t^2}\right)_x \, .

Finally, the derivative of the number density with respect to the mass
is more involved because of the mass-dependent prefactor.

.. math::

   \begin{eqnarray}
   \frac{d n}{d m} &=& \frac{4 n}{m}+
   \left(\frac{g m^4}{2 \pi^2}\right) \frac{d}{d m}
   \left(\frac{1}{m}\frac{d \tilde{P}}{d x} \right) =
   \frac{4 n}{m} + 
   \left(\frac{g m^4}{2 \pi^2}\right)
   \left[\frac{1}{m}\left(\frac{d^2\tilde{P}}{dx^2}\frac{dx}{dm}+
   \frac{d^2\tilde{P}}{dt dx}\frac{dt}{dm}\right)-
   \frac{1}{m^2}\frac{d \tilde{P}}{d x}\right] \nonumber \\
   &=& \frac{4 n}{m} - \left(\frac{g m^2}{2 \pi^2}\right)
   \left( \frac{d\tilde{P}}{dx}
   +\frac{\mu}{m} \frac{d^2\tilde{P}}{dx^2}
   +\frac{T}{m} \frac{d^2\tilde{P}}{dt dx} \right) = 
   \frac{3n}{m} -\left[(x+1) \left(\frac{dn}{d\mu}\right) +
   t \left(\frac{dn}{dT}\right) \right] \nonumber
   \end{eqnarray}

These expansions are used in
:cpp:func:`o2scl::fermion_thermo_tl::calc_mu_deg()` and
:cpp:func:`o2scl::fermion_deriv_thermo_tl::calc_mu_deg()`.

Nondegenerate Expansion
-----------------------

There is a useful identity ([Chandrasekhar10]_ and [Tooper69]_)

.. math::

   \int_0^{\infty} \frac{x^4 \left(x^2+z^2\right)^{-1/2}~dx}
   {1+e^{\sqrt{x^2+z^2}-\phi}} = 
   3 z^2 \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^2} e^{n \phi} K_2(n z)

which works well when :math:`\phi-z < -1`. This result directly 
gives the sum in  Johns96

.. math::

   P = \frac{g m^4}{2 \pi^2} \sum_{k=1}^{\infty} P_k \equiv 
   \frac{g m^4}{2 \pi^2} \left[ \sum_{k=1}^{\infty}
   \frac{t^2 (-1)^{k+1}}{k^2} e^{k x/t} e^{k/t} K_2\left(\frac{k}{t}\right)
   \right]

The function :math:`e^{y} K_2(y)` is implemented in GSL as
``gsl_sf_bessel_Kn_scaled()``. In the case that one
wants to include antiparticles, the result is
similar

.. math::

   P = \frac{g m^4}{2 \pi^2} \sum_{k=1}^{\infty} \bar{P}_k \equiv 
   \frac{g m^4}{2 \pi^2} \left\{ \sum_{k=1}^{\infty}
   \frac{2 t^2 (-1)^{k+1}}{k^2} e^{-k/t} \mathrm{cosh}
   \left[k(x+1)/t\right] \left[ e^{k/t} 
   K_2\left(\frac{k}{t}\right) \right]
   \right\}

where the scaled Bessel function has been separated out.
Similarly defining

.. math::

   n = \frac{g m^3}{2 \pi^2} \sum_{k=1}^{\infty} n_k  \, ,

the terms in the expansion for the density (without and
with antiparticles) are

.. math::

   \begin{eqnarray}
   n_k &=& \frac{k}{t}{P_k}
   \nonumber \\
   \bar{n}_k &=& \frac{k}{t}{\bar{P}_k} 
   \mathrm{tanh} \left[k (x+1)/t\right]
   \end{eqnarray}

The entropy terms (with and without antiparticles) are
   
.. math::

   \begin{eqnarray}
   s_k &=& \left( \frac{4t-kx-k}{kt}\right) n_k +
   \frac{(-1)^{k+1}}{k} e^{k x/t} \left[ e^{k/t} K_1(k/t) \right]
   \nonumber \\
   \bar{s}_k &=& 
   -\frac{(1+x)\bar{n}_k}{t} +
   \frac{2(-1)^{k+1}}{k}  e^{-k/t} \mathrm{cosh}[k(x+1)/t] 
   \left[ e^{k/t} K_3(k/t) \right]
   \end{eqnarray}

included. For the derivatives, no additional Bessel functions are
required
   
.. math::

   \begin{eqnarray}
   \left(\frac{dn}{d\mu}\right)_k &=& 
   \frac{k}{t}{n_k} \\
   \left(\frac{d\bar{n}}{d\mu}\right)_k &=&
   \frac{k}{t}{\bar{n}_k} \\
   \left(\frac{dn}{dT}\right)_k &=& 
   \frac{k}{t} s_k - \frac{1}{t} n_k \\
   \left(\frac{d\bar{n}}{dT}\right)_k &=& 
   \frac{k}{t} \bar{s}_k \mathrm{tanh}\left[k(x+1)/t\right]
   - \left\{ t+2 k (1+x) \mathrm{csch}\left[k(x+1)/t\right]
   \right\} \frac{\bar{n}_k}{t^2} \\
   \left(\frac{ds}{dT}\right)_k &=& 
   \left[ \frac{3t -2k x -2 k}{t^2}\right] s_k
   + \left[ \frac{5 k t - 2 k^2 x +5 k t x - k^2 x^2}{k t^3}\right] n_k \\
   \left(\frac{d\bar{s}}{dT}\right)_k &=& 
   \left\{2 k (1+x) \mathrm{tanh}\left[ k(1+x)/t\right] - 3 t\right\}
   \frac{\bar{s}_k}{t^2} +
   \left\{2 k^2 (1+x)^2 \mathrm{tanh}\left[ k(1+x)/t\right] - 
   \right. \nonumber \\
   && \left.
   k^2 (2 + 2 x + x^2) \mathrm{coth}\left[ k(1+x)/t\right] -
   5 k(1+x) t \right\}
   \frac{\bar{n}_k}{k t^3}
   \end{eqnarray}

These expansions are used in
:cpp:func:`o2scl::fermion_thermo_tl::calc_mu_ndeg()`.
 
