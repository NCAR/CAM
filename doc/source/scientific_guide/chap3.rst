.. _dynamics:

*************
 Dynamics
*************

Finite Volume Dynamical Core
----------------------------

Overview
~~~~~~~~

This document describes the Finite-Volume (FV) dynamical core that was
initially developed and used at the NASA Data Assimilation Office (DAO)
for data assimilation, numerical weather predictions, and climate
simulations. The finite-volume discretization is local and entirely in
physical space. The horizontal discretization is based on a conservative
“*flux-form semi-Lagrangian*” scheme described by (hereafter LR96) and
(hereafter LR97). The vertical discretization can be best described as
*Lagrangian* with a conservative re-mapping, which essentially makes it
*quasi-Lagrangian*. The *quasi-Lagrangian* aspect of the vertical
coordinate is transparent to model users or physical parameterization
developers, and it functions exactly like the :math:`\eta -coordinate`
(a hybrid :math:`\sigma -p` coordinate) used by other dynamical cores within CAM.

In the current implementation for use in CAM, the FV dynamics and
physics are “time split” in the sense that all prognostic variables are
updated sequentially by the “dynamics” and then the “physics”. The time
integration within the FV dynamics is fully explicit, with sub-cycling
within the 2D Lagrangian dynamics to stabilize the fastest wave (see
section [FVvdisc]). The transport for tracers, however, can take a much
larger time step (*e.g.*, 30 minutes as for the physics).

The governing equations for the hydrostatic atmosphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For reference purposes, we present the continuous differential equations
for the hydrostatic 3D atmospheric flow on the sphere for a general
vertical coordinate :math:`\zeta` (*e.g*., ). Using standard
notations, the hydrostatic balance equation is given as follows:

.. math::

   \label{hydro}
   \frac{1}{\rho }\frac{\partial p}{\partial z}+g=0 ,

where :math:`\rho` is the density of the air, *p* the pressure, and
*g* the gravitational constant. Introducing the “*pseudo-density*”
:math:`\pi =\frac{\partial p}{\partial \zeta }` (*i.e.*, the vertical pressure gradient in the general coordinate),
from the hydrostatic balance equation the *pseudo-density* and the true
density are related as follows:

.. math::
   :label: hodrostatic-pi

   \pi = -\frac{\partial \Phi }{\partial \zeta }\rho ,

where :math:`\Phi =gz` is the geopotential. Note that :math:`\pi`
reduces to the “true density” if :math:`\zeta =-gz`, and the “surface
pressure” :math:`P_{s}` if :math:`\zeta =\sigma` (:math:`\sigma
=\frac{p}{P_{s}}`). The conservation of total air mass using
:math:`\pi` as the prognostic variable can be written as

.. math::
   :label: mass-pi

   \frac{\partial }{\partial t}\pi +\nabla \cdot (\overrightarrow{V}\pi ) =0 ,

where :math:`\overrightarrow{V}=(u,v,\frac{d\zeta }{dt})`. Similarly,
the mass conservation law for tracer species (or water vapor) can be
written as

.. math::

   \label{tracer-pi}
   \frac{\partial }{\partial t}(\pi q)+\nabla \cdot
   (\overrightarrow{V}\pi q) =0 ,

where *q* is the mass mixing ratio (or specific humidity) of the tracers
(or water vapor).

Choosing the (virtual) potential temperature :math:`\Theta` as the
thermodynamic variable, the first law of thermodynamics is written as

.. math::

   \label{thermo-pi}
   \frac{\partial }{\partial t}(\pi \Theta )+\nabla \cdot (
   \overrightarrow{V}\pi \Theta ) =0 .

Letting :math:`(\lambda ,\theta )` denote the (longitude, latitude)
coordinate, the momentum equations can be written in the
“vector-invariant form” as follows:

.. math::

   \label{u-pi}
   \frac{\partial }{\partial t}u=\Omega v-\frac{1}{Acos\theta }[
     \frac{\partial }{\partial \lambda }( \kappa +\Phi -\nu D)
     +\frac{1}{\rho }\frac{\partial }{\partial \lambda }p]
   -\frac{d\zeta }{dt}\frac{\partial u}{\partial \zeta } ,

.. math::

   \label{v-pi}
   \frac{\partial }{\partial t}v=-\Omega u-\frac{1}{A}[
     \frac{\partial }{\partial \theta }( \kappa +\Phi -\nu D)
     +\frac{1}{\rho }\frac{\partial }{\partial \theta }p]
   -\frac{d\zeta }{dt}\frac{\partial v}{\partial \zeta } ,

where *A* is the radius of the earth, :math:`\nu` is the coefficient
for the optional divergence damping, *D* is the horizontal divergence

.. math::

   D=\frac{1}{Acos\theta }[ \frac{\partial
   }{\partial \lambda }(u)+\frac{\partial }{\partial \theta }(v\,
   cos\theta )] ,

.. math:: \kappa =\frac{1}{2}( u^{2}+v^{2}) ,

and :math:`\Omega`, the vertical component of the absolute vorticity,
is defined as follows:

.. math::

   \Omega =2\omega \, sin\theta +\frac{1}{A\, cos\theta }[
   \frac{\partial }{\partial \lambda }v-\frac{\partial }{\partial \theta
   }(u\, cos\theta )] ,

where :math:`\omega` is the angular velocity of the earth. Note that
the last term in ([u-pi]) and ([v-pi]) vanishes if the vertical
coordinate :math:`\zeta` is a conservative quantity (*e.g*., entropy
under adiabatic conditions or an imaginary conservative tracer), and the
3D divergence operator becomes 2D along constant :math:`\zeta`
surfaces. The discretization of the 2D horizontal transport process is
described in section [FVhdisc]. The complete dynamical system using the
Lagrangian control-volume vertical discretization is described in
section [FVvdisc]. A mass, momentum, and total energy conservative
mapping algorithm is described in section [FVmap].

Horizontal discretization of the transport process on the sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the vertical transport term would vanish after the introduction of
the vertical Lagrangian control-volume discretization (see section
[FVvdisc]), we shall present here only the 2D (horizontal) forms of the
FFSL transport algorithm for the transport of density ([mass-pi]) and
mixing ratio-like quantities ([tracer-pi]) on the sphere. The governing
equation for the pseudo-density ([mass-pi]) becomes

.. math::

   \label{pi-2d}
   \frac{\partial }{\partial t}\pi +\frac{1}{Acos\theta }[
   \frac{\partial }{\partial \lambda }(u\pi )+\frac{\partial }{\partial
   \theta }(v\pi \, cos\theta )] =0 .

The finite-volume (*integral*) representation of the continuous
:math:`\pi` field is defined as follows:

Given the *exact* 2D wind field :math:`
\overrightarrow{V}(t;\lambda ,\theta )=(U,V)` the 2D integral
representation of the conservation law for :math:`\widetilde{\pi }`
can be obtained by integrating ([pi-2d]) in time and in space

The above 2D transport equation is still *exact* *for the finite-volume
under consideration*. To carry out the contour integral, certain
approximations must be made. LR96 essentially decomposed the flux
integral using two orthogonal 1D flux-form transport operators.
Introducing the following difference operator

.. math:: \delta _{x}q=q(x+\frac{\Delta x}{2})-q(x-\frac{\Delta x}{2}),

and assuming :math:`(u^{*},v^{*})` is the time-averaged (from time
:math:`t` to time :math:`t+\Delta t`) :math:`\overrightarrow{V}` on the
C-grid (*e.g*., Fig. 1 in LR96), the 1-D finite-volume flux-form
transport operator *F* in the :math:`\lambda` -direction is

.. math::

   \label{xtp}
   F(u^{*},\Delta t,\widetilde{\pi })=-\frac{1}{A\Delta \lambda cos\theta
   }\, \delta _{\lambda }[ \int _{t}^{t+\Delta t}\pi U\, dt]
   =-\frac{\Delta t}{A\Delta \lambda cos\theta }\, \delta _{\lambda
   }[ \chi (u^{*},\Delta t;\pi )] ,

where :math:`\chi` , the time-accumulated (from *t* to
*t*\ +\ *:math:`\Delta`\ t*) mass flux across the cell wall, is
defined as follows,

.. math::

   \label{xmass}
   \chi (u^{*},\Delta t;\pi )=\frac{1}{\Delta t}\int _{t}^{t+\Delta t}\pi
   U\, dt\equiv u^{*}\pi ^{*}(u^{*},\Delta t,\widetilde{\pi }) ,

and

.. math::

   \label{pi-*}
   \pi ^{*}(u^{*},\Delta t;\widetilde{\pi })\approx \frac{1}{\Delta
   t}\int _{t}^{t+\Delta t}\pi \, dt

can be interpreted as a time mean (from time :math:`t` to time :math:`
t+\Delta t`) pseudo-density value of all material that passed through
the cell edge from the upwind direction.

Note that the above *time integration* is to be carried out along the
*backward-in-time* trajectory of the cell edge position from
:math:`t=t+\Delta t` (the arrival point; (*e.g*., point B in Fig. 3 of
LR96) back to time :math:`t` (the departure point; *e.g*., point B’ in
Fig. 3 of LR96). The very essence of the 1D finite-volume algorithm is
to construct, based on the given initial cell-mean values of :math:`
\widetilde{\pi }`, an approximated subgrid distribution of the true
:math:`\pi` field, to enable an analytic integration of ([pi-\*]).
Assuming there is no error in obtaining the time-mean wind
:math:`(u^{*})`, the only error produced by the 1D transport scheme
would be solely due to the approximation to the continuous distribution
of :math:`\pi` within the subgrid under consideration. From this
perspective, it can be said that the 1D finite-volume transport
algorithm combines the time-space discretization in the approximation of
the time-mean cell-edge values :math:`\pi ^{*}`. The physically
correct way of approximating the integral ([pi-\*]) must be “upwind”, in
the sense that it is integrated along the backward trajectory of the
cell edges. For example, a center difference approximation to ([pi-\*])
would be physically incorrect, and consequently numerically unstable
unless artificial numerical diffusion is added.

Central to the accuracy and computational efficiency of the
finite-volume algorithms is the degrees of freedom that describe the
subgrid distribution. The first order upwind scheme, for example, has
zero degrees of freedom within the volume as it is assumed that the
subgrid distribution is piecewise constant having the same value as the
given volume-mean. The second order finite-volume scheme (*e.g*., )
assumes a piece-wise linear subgrid distribution, which allows one
degree of freedom for the specification of the “slope” of the linear
distribution to improve the accuracy of integrating ([pi-\*]). The
Piecewise Parabolic Method (PPM, ) has two degrees of freedom in the
construction of the second order polynomial within the volume, and as a
result, the accuracy is significantly enhanced. The PPM appears to
strike a good balance between computational efficiency and accuracy.
Therefore, the PPM is the basic 1D scheme we chose. (An extension of the
standard PPM by S.-J. Lin has also been documented in ). Note that the
subgrid PPM distributions are compact, and do not extend beyond the
volume under consideration. The accuracy is therefore significantly
better than the order of the chosen polynomials implies. While the PPM
scheme possesses all the desirable attributes (mass conserving,
monotonicity preserving, and high-order accuracy) in 1D, it is important
that a solution be found to avoid the directional splitting in the
multi-dimensional problem of modeling the dynamics and transport
processes of the Earth’s atmosphere.

The first step for reducing the splitting error is to apply the two
orthogonal 1D flux-form operators in a directionally symmetric way.
After symmetry is achieved, the “inner operators” are then replaced with
corresponding advective-form operators. A consistent advective-form
operator in the :math:`\lambda -`\ direction can be derived from its
flux-form counterpart (:math:`F`) as follows:

.. math::

   \label{xadv}
   f(u^{\*},\Delta t,\widetilde{\pi })=F(u^{\*},\Delta t,\widetilde{\pi
   })+\widetilde{\rho }\, F(u^{\*},\Delta t,\widetilde{\pi }\equiv
   1)=F(u^{\*},\Delta t,\widetilde{\pi })+\widetilde{\pi }\,
   C_{def}^{\lambda } ,

.. math::

   \label{xdef}
   C^{\lambda }_{def}=\frac{\Delta t\, \delta _{\lambda }u^{\*}}{A\Delta
   \lambda cos\theta } ,

where :math:`C_{def}^{\lambda }` is a dimensionless number indicating
the degree of the flow deformation in the :math:`\lambda`-direction.
The above derivation of :math:`f` is slightly different from LR96’s
approach, which adopted the traditional 1D advective-form
semi-Lagrangian scheme. The advantage of using ([xadv]) is that
computation of winds at cell centers (Eq. 2.25 in LR96) are avoided.

Analogously, the **1D flux-form transport operator** *G* **in the
latitudinal (:math:`\theta`) direction is derived as follows:**

.. math::

   \label{ytp}
   G(v^{*},\Delta t,\widetilde{\pi })=-\frac{1}{A\Delta \theta cos\theta
   }\, \delta _{\theta }[ \int _{t}^{t+\Delta t}\pi Vcos\theta \,
   dt] =-\frac{\Delta t}{A\Delta \theta cos\theta }\, \delta
   _{\theta }[ v^{*}cos\theta \, \pi ^{*}] ,

and likewise the advective-form operator,

.. math::

   \label{yadv}
   g(v^{*},\Delta t,\widetilde{\pi })=G(v^{*},\Delta t,\widetilde{\pi
   })+\widetilde{\pi }\, C_{def}^{\theta } ,

 where

.. math::

   \label{ydef}
   C^{\theta }_{def}=\frac{\Delta t\, \delta _{\theta }[
   v^{*}cos\theta ] }{A\Delta \theta cos\theta } .

To complete the construction of the 2D algorithm on the sphere, we
introduce the following short hand notations:

.. math::
   :label: def-1

   (\, )^{\theta }=(\, )^{n}+\frac{1}{2}g[ v^{*},\Delta t,\, (\,
   )^{n}] ,

.. math::
   :label: def-2

   (\, )^{\lambda }=(\, )^{n}+\frac{1}{2}f[ u^{*},\Delta t,\, (\,)^{n}] .

The 2D transport algorithm (*cf*, Eq. 2.24 in LR96) can then be written as

.. math::
   :label: den-gf

   \widetilde{\pi }^{n+1}=\widetilde{\pi }^{n}+F[ u^{\*},\Delta
   t,\widetilde{\pi }^{\theta }] +G[ v^{\*},\Delta t,\widetilde{\pi }^{\lambda }] .

Using explicitly the mass fluxes :math:`( \chi ,Y)`, ([den-gf]) is rewritten as

.. math::
   :label: air

   \widetilde{\pi }^{n+1}=\widetilde{\pi }^{n}-\frac{\Delta t}{Acos\theta
   }\{ \frac{1}{\Delta \lambda }\delta _{\lambda }[ \chi
   (u^{\*},\Delta t;\widetilde{\pi }^{\theta })] +\frac{1}{\Delta
   \theta }\delta _{\theta }[ cos\theta \, Y(v^{*},\Delta
   t;\widetilde{\pi }^{\lambda })] \} ,

where :math:`Y`, the mass flux in the meridional direction, is
defined in a similar fashion as :math:`\chi` ([xmass]). It can be
verified that in the special case of constant density flow
(:math:`\widetilde{\pi}=constant)` the above equation degenerates to the finite-difference
representation of the *incompressibility condition* of the “time mean”
wind field :math:`(u^{*},v^{*})`, *i.e*.,

.. math::
   :label: div=0

   \frac{1}{\Delta \lambda }\delta _{\lambda }u^{\*}+\frac{1}{\Delta
   \theta }\delta _{\theta }( v^{\*}cos\theta ) =0 .

The fulfillment of the above *incompressibility condition* for constant
density flows is crucial to the accuracy of the 2D flux-form
formulation. For transport of volume mean mixing ratio-like quantities
:math:`(\widetilde{q})` the mass fluxes :math:`(\chi ,Y)` as defined
previously should be used as follows

.. math::

   \label{tracer}
   \widetilde{q}^{n+1}=\frac{1}{\widetilde{\pi }^{n+1}}[
   \widetilde{\pi }^{n}\widetilde{q}^{n}+F(\chi ,\Delta
   t,\widetilde{q}^{\theta })+G(Y,\Delta t,\widetilde{q}^{\lambda
   })] .

Note that the above form of the tracer transport equation consistently
degenerates to ([den-gf]) if :math:`\widetilde{q}\equiv 1` (*i.e*.,
the tracer density equals to the background air density), which is
another important condition for a flux-form transport algorithm to be
able to avoid generation of noise (*e.g*., creation of artificial
gradients) and to maintain mass conservation.

A *vertically Lagrangian* and *horizontally Eulerian* control-volume discretization of the hydrodynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The very idea of using Lagrangian vertical coordinate for formulating
governing equations for the atmosphere is not entirely new. ) is likely
the first to have formulated, in the *continuous differential form*, the
governing equations using a Lagrangian coordinate. Starr did not make
use of the *discrete* Lagrangian control-volume concept for
discretization nor did he present a solution to the problem of computing
the pressure gradient forces. In the *finite-volume discretization* to
be described here, the Lagrangian surfaces are treated as the bounding
material surfaces of the Lagrangian control-volumes within which the
finite-volume algorithms developed in LR96, LR97, and L97 will be
directly applied.

To use a vertical Lagrangian coordinate system to reduce the 3D
governing equations to the 2D forms, one must first address the issue of
whether it is an inertial coordinate or not. For hydrostatic flows, it
is. This is because both the right-hand-side and the left-hand-side of
the vertical momentum equation vanish for purely hydrostatic flows.

Realizing that the earth’s surface, for all practical modeling purposes,
can be regarded as a non-penetrable material surface, it becomes
straightforward to construct a terrain-following Lagrangian
control-volume coordinate system. In fact, any commonly used
terrain-following coordinates can be used as the starting reference
(*i.e*., fixed, Eulerian coordinate) of the floating Lagrangian
coordinate system. To close the coordinate system, the model top (at a
prescribed constant pressure) is also assumed to be a Lagrangian
surface, which is the same assumption being used by practically all
global hydrostatic models.

The basic idea is to start the time marching from the chosen
terrain-following Eulerian coordinate (*e.g*., pure :math:`\sigma` or
hybrid :math:`\sigma`-p), *treating the initial coordinate surfaces as
material surfaces*, the finite-volumes bounded by two coordinate
surfaces, *i.e*., the Lagrangian control-volumes, are free vertically,
to float, compress, or expand with the flow as dictated by the
hydrostatic dynamics.

By choosing an imaginary conservative tracer :math:`\zeta` that is a
monotonic function of height and constant on the initial reference
coordinate surfaces (*e.g*., the value of “:math:`\eta`” in the hybrid
:math:`\sigma -p` coordinate used in CAM), the 3D governing equations
written for the general vertical coordinate in section 1.2 can be
reduced to 2D forms. After factoring out the constant :math:`\delta
\zeta`, ([mass-pi]), the conservation law for the pseudo-density
(:math:`\pi=\frac{\delta p}{\delta \zeta }`), becomes

.. math::

   \label{mass-lcv}
   \frac{\partial }{\partial t}\delta p+\frac{1}{Acos\theta }[
   \frac{\partial }{\partial \lambda }(u\delta p)+\frac{\partial
   }{\partial \theta }(v\delta p\, cos\theta )] =0 ,

where the symbol :math:`\delta` represents the vertical difference
between the two neighboring Lagrangian surfaces that bound the finite
control-volume. From ([hydro]), the pressure thickness
:math:`\delta p` of that control-volume is proportional to the total
mass, *i.e*., :math:`\delta p=-\rho g\delta z`. Therefore, it can be said that the
Lagrangian control-volume vertical discretization has the hydrostatic
balance built-in, and :math:`\delta p` can be regarded as the
“pseudo-density” for the discretized Lagrangian vertical coordinate system.

Similarly, ([tracer-pi]), the mass conservation law for all tracer
species, is

.. math::

   \label{tracer-lcv}
   \frac{\partial }{\partial t}(q\delta p)+\frac{1}{Acos\theta }[
   \frac{\partial }{\partial \lambda }(uq\delta p)+\frac{\partial
   }{\partial \theta }(vq\delta p\, cos\theta )] =0,

the thermodynamic equation, ([thermo-pi]), becomes

.. math::

   \label{Thermo-lcv}
   \frac{\partial }{\partial t}(\Theta \delta p)+\frac{1}{Acos\theta
   }[ \frac{\partial }{\partial \lambda }(u\Theta \delta
   p)+\frac{\partial }{\partial \theta }(v\Theta \delta p\, cos\theta
   )] =0,

and ([u-pi]) and ([v-pi]), the momentum equations, are reduced to

.. math::

   \label{u-lcv}
   \frac{\partial }{\partial t}u=\Omega v-\frac{1}{Acos\theta }[
   \frac{\partial }{\partial \lambda }( \kappa +\Phi -\nu D)
   +\frac{1}{\rho }\frac{\partial }{\partial \lambda }p] ,

.. math::

   \label{v-lcv}
   \frac{\partial }{\partial t}v=-\Omega u-\frac{1}{A}[
   \frac{\partial }{\partial \theta }( \kappa +\Phi -\nu D)
   +\frac{1}{\rho }\frac{\partial }{\partial \theta }p] .

Given the prescribed pressure at the model top :math:`P_{\infty }`,
the position of each Lagrangian surface :math:`P_{l}` (horizontal
subscripts omitted) is determined in terms of the hydrostatic pressure
as follows:

.. math::

   \label{L-coord}
   P_{l}=P_{\infty }+\sum ^{l}_{k=1}\delta P_{k},\, \, \, \, \, (for\,
   l=1,\, 2,\, 3,\, ...,\, N) ,

where the subscript *:math:`l`* is the vertical index ranging from 1
at the lower bounding Lagrangian surface of the first (the highest)
layer to *:math:`N`* at the Earth’s surface. There are :math:`N`\ +1
Lagrangian surfaces to define a total number of *:math:`N`* Lagrangian
layers. The surface pressure, which is the pressure at the lowest
Lagrangian surface, is easily computed as :math:`P_{N}` using
([L-coord]). The surface pressure is needed for the physical
parameterizations and to define the reference Eulerian coordinate for
the mapping procedure (to be described in section [FVmap]).

With the exception of the pressure-gradient terms and the addition of a
thermodynamic equation, the above 2D Lagrangian dynamical system is the
same as the shallow water system described in LR97. The conservation law
for the depth of fluid :math:`h` in the shallow water system of LR97
is replaced by ([mass-lcv]) for the pressure thickness :math:`
\delta p`. The ideal gas law, the mass conservation law for air mass,
the conservation law for the potential temperature ([Thermo-lcv]),
together with the modified momentum equations ([u-lcv]) and ([v-lcv])
close the 2D Lagrangian dynamical system, which are vertically coupled
only by the hydrostatic relation (see ([hydro-PT]), section [FVmap]).

The time marching procedure for the 2D Lagrangian dynamics follows
closely that of the shallow water dynamics fully described in LR97. For
computational efficiency, we shall take advantage of the stability of
the FFSL transport algorithm by using a much larger time step
(:math:`\Delta t)` for the transport of all tracer species (including
water vapor). As in the shallow water system, the Lagrangian dynamics
uses a relatively small time step, :math:`\Delta \tau
=\Delta t/m`, where :math:`m` is the number of the sub-cycling needed
to stabilize the fastest wave in the system. We shall describe here this
time-split procedure for the *prognostic variables* :math:`
[ \delta p,\Theta ,u,v;q]` on the D-grid. Discretization on
the C-grid for obtaining the *diagnostic variables*, the time-averaged
winds :math:`(u^{*},v^{*})`, is analogous to that of the D-grid (see
also LR97).

Introducing the following short hand notations (*cf*, ([def-1]) and ([def-2])):

.. math::

   (\, )_{i}^{\theta }=(\,
   )^{n+\frac{i-1}{m}}+\frac{1}{2}g[v_{i}^{*},\Delta \tau ,(\,
   )^{n+\frac{i-1}{m}}],

.. math::

   (\, )_{i}^{\lambda }=(\,
   )^{n+\frac{i-1}{m}}+\frac{1}{2}f[u_{i}^{*},\Delta \tau ,(\,
   )^{n+\frac{i-1}{m}}],

and applying directly ([air]), the update of “pressure thickness”
:math:`\delta p`, using the fractional time step
:math:`\Delta \tau =\Delta
t/m`, can be written as

.. math::
   :label: mass

   \delta p^{n+\frac{i}{m}}=\delta p^{n+\frac{i-1}{m}}-\frac{\Delta \tau
   }{Acos\theta }\{ \frac{1}{\Delta \lambda }\delta _{\lambda
   }[ x_{i}^{*}(u_{i}^{*},\Delta \tau ;\delta p_{i}^{\theta
   })] +\frac{1}{\Delta \theta }\delta _{\theta }[ cos\theta
   \, y_{i}^{*}(v_{i}^{*},\Delta \tau ;\delta p_{i}^{\lambda })] \}

.. math:: (for\, i=1,...,m),

where :math:`[ x_{i}^{*},y_{i}^{*}]` are the background
air mass fluxes, which are then used as input to Eq. 24 for transport of
the potential temperature :math:`\Theta`:

.. math::

   \label{pt}
   \Theta ^{n+\frac{i}{m}}=\frac{1}{\delta p^{n+\frac{i}{m}}}[
   \delta p^{n+\frac{i-1}{m}}\Theta ^{n+\frac{i-1}{m}}+F(x_{i}^{*},\Delta
   \tau ;\Theta _{i}^{\theta })+G(y_{i}^{*},\Delta \tau ,\Theta
   _{i}^{\lambda })] .

The discretized momentum equations for the shallow water system (*cf*,
Eq. 16 and Eq. 17 in LR97) are modified for the pressure gradient terms
as follows:

.. math::

   \label{u}
   u^{n+\frac{i}{m}}=u^{n+\frac{i-1}{m}}+\Delta \tau \, [
   y_{i}^{*}( v_{i}^{*},\Delta \tau ;\Omega ^{\lambda })
   -\frac{1}{A\Delta \lambda cos\theta }\delta _{\lambda }(\kappa
   ^{*}-\nu D^{*})+\widehat{P_{\lambda }}] ,

.. math::

   \label{v}
   v^{n+\frac{i}{m}}=v^{n+\frac{i-1}{m}}-\Delta \tau \, [
   x_{i}^{*}( u_{i}^{*},\Delta \tau ;\Omega ^{\theta })
   +\frac{1}{A\Delta \theta }\delta _{\theta }(\kappa ^{*}-\nu
   D^{*})-\widehat{P_{\theta }}] ,

where :math:`\kappa ^{*}` is the upwind-biased “kinetic energy” (as
defined by Eq. 18 in LR97), and :math:`D^{*}`, the horizontal
divergence on the D-grid, is discretized as follows:

.. math::

   D^{*}=\frac{1}{Acos\theta }[ \frac{1}{\Delta \lambda }\delta
   _{\lambda }u^{n+\frac{i-1}{m}}+\frac{1}{\Delta \theta }\delta _{\theta
   }( v^{n+\frac{i-1}{m}}cos\theta ) ] .

The finite-volume mean pressure-gradient terms in ([u]) and ([v]) are
computed as follows:

.. math::

   \label{px}
   \widehat{P_{\lambda }}=\frac{\oint _{\Pi leftharpoons \lambda
   }\phi d\Pi }{Acos\theta \, \oint _{\Pi leftharpoons \lambda }\Pi
   d\lambda } ,

.. math::

   \label{py}
   \widehat{P_{\theta }}=\frac{\oint _{\Pi leftharpoons \theta
   }\phi d\Pi }{A\, \oint _{\Pi leftharpoons \theta }\Pi d\theta } ,

where :math:`\Pi =p^{\kappa }\, (\kappa =R/C_{p})`, and the symbols
“:math:`\Pi leftharpoons \lambda`” and “:math:`\Pi` :math:`leftharpoons \theta`” 
indicate that the contour integrations are
to be carried out, using the finite-volume algorithm described in L97,
in the :math:`(\Pi ,\lambda )` and :math:`(\Pi ,\theta )` space, respectively.

To complete one time step, equations ([mass]-[v]), together with their
counterparts on the C-grid are cycled :math:`m` times using the
fractional time step :math:`\Delta \tau`, which are followed by the
tracer transport using ([tracer-lcv]) with the large-time-step
:math:`\Delta t`.

Mass fluxes :math:`(x^{*},y^{*})` and the winds
:math:`(u^{*},v^{*})` on the C-grid are accumulated for the
large-time-step transport of tracer species (including water vapor)
:math:`q` as

.. math::

   \label{tracers}
   q_{}^{n+1}=\frac{1}{\delta p^{n+1}}[ q_{}^{n}\delta
   p^{n}+F(X^{*},\Delta t,q_{}^{\theta })+G(Y^{*},\Delta t,q_{}^{\lambda
   })] ,

where the time-accumulated mass fluxes :math:`(X^{*},Y^{*})` are
computed as

.. math::

   \label{x-mass}
   X^{*}=\sum ^{m}_{i=1}x_{i}^{*}(u_{i}^{*},\, \Delta \tau ,\, \delta
   p_{i}^{\theta }) ,

.. math::

   \label{y-mass}
   Y^{*}=\sum _{i=1}^{m}y_{i}^{*}(v_{i}^{*},\, \Delta \tau ,\, \delta
   p_{i}^{\lambda }) .

The time-averaged winds :math:`(U^{*},V^{*})`, defined as follows, are
to be used as input for the computations of :math:`q^{\lambda }` and
:math:`q^{\theta }:`

.. math::

   \label{u-wind}
   U^{*}=\frac{1}{m}\sum ^{m}_{i=1}u_{i}^{*} ,

.. math::

   \label{v-wind}
   V^{*}=\frac{1}{m}\sum ^{m}_{i=1}v_{i}^{*} .

The use of the time accumulated mass fluxes and the time-averaged winds
for the large-time-step tracer transport in the manner described above
ensures the conservation of the tracer mass and maintains the highest
degree of consistency possible given the time split integration
procedure.

The algorithm described here can be readily applied to a regional model
if appropriate boundary conditions are supplied. There is formally no
Courant number related time step restriction associated with the
transport processes. There is, however, a stability condition imposed by
the gravity-wave processes. For application on the whole sphere, it is
computationally advantageous to apply a polar filter to allow a dramatic
increase of the size of the small time step :math:`\Delta
\tau`. The effect of the polar filter is to stabilize the
short-in-wavelength (and high-in-frequency) gravity waves that are being
unnecessarily and unidirectionally resolved at very high latitudes in
the zonal direction. To minimize the impact to meteorologically
significant larger scale waves, the polar filter is highly scale
selective and is applied only to the diagnostic variables on the
auxiliary C-grid and the tendency terms in the D-grid momentum
equations. No polar filter is applied directly to any of the prognostic
variables.

The design of the polar filter follows closely that of for the C-grid
Arakawa type dynamical core (*e.g*., ). For the the fast-fourier
transform component of the polar filtering has replaced the algebraic
form at all filtering latitudes. Because our prognostic variables are
computed on the D-grid and the fact that the FFSL transport scheme is
stable for Courant number greater than one, in realistic test cases the
maximum size of the time step is about two to three times larger than a
model based on Arakawa and Lamb’s C-grid differencing scheme. It is
possible to avoid the use of the polar filter if, for example, the
“Cubed grid” is chosen, instead of the current latitude-longitude grid.
However, this would require a significant rewrite of the rest of the
model codes including physics parameterizations, the land model, and
most of the post processing packages.

The size of the small time step for the Lagrangian dynamics is only a
function of the horizontal resolution. Applying the polar filter, for
the 2-degree horizontal resolution, a small-time-step size of 450
seconds can be used for the Lagrangian dynamics. From the
large-time-step transport perspective, the small-time-step integration
of the 2D Lagrangian dynamics can be regarded as a very accurate
iterative solver, with *m* iterations, for computing the time mean winds
and the mass fluxes, analogous in functionality to a semi-implicit
algorithm’s elliptic solver (*e.g*., ). Besides accuracy, the merit of
an “explicit” versus “semi-implicit” algorithm ultimately depends on the
computational efficiency of each approach. In light of the advantage of
the explicit algorithm in parallelization, we do not regard the explicit
algorithm for the Lagrangian dynamics as an impedance to computational
efficiency, particularly on modern parallel computing platforms.
Furthermore, it may be possible to further increase the size of the
small time step via vertical mode decomposition. This approach is one of
the algorithm design issues we plan to revisit.

A mass, momentum, and total energy conserving mapping algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Lagrangian surfaces that bound the finite-volume will eventually
deform, particularly in the presence of persistent diabatic
heating/cooling, in a time scale of a few hours to a day depending on
the strength of the heating and cooling, to a degree that it will
negatively impact the accuracy of the
horizontal-to-Lagrangian-coordinate transport and the computation of the
pressure gradient forces. Therefore, a key to the success of the
Lagrangian control-volume discretization is an accurate and conservative
algorithm for mapping the deformed Lagrangian coordinate back to a fixed
reference Eulerian coordinate.

There are some degrees of freedom in the design of the vertical mapping
algorithm. To ensure conservation, our current (and recommended) mapping
algorithm is based on the reconstruction of the “mass” (pressure
thickness :math:`\delta p`), zonal and meridional “winds”, “tracer
mixing ratios”, and “total energy” (volume integrated sum of the
internal, potential, and kinetic energy), using the monotonic Piecewise
Parabolic sub-grid distributions with the hydrostatic pressure (as
defined by ([L-coord])) as the mapping coordinate. We outline the
mapping procedure as follows.

- **Step 1**: Define a suitable Eulerian reference coordinate. The
  mass in each layer (:math:`\delta p`) is then distributed
  vertically according to the chosen Eulerian coordinate. The surface
  pressure typically plays an “anchoring” role in defining the terrain
  following Eulerian vertical coordinate. The hybrid
  :math:`\eta -coordinate` used in the NCAR CCM3 is adopted in the current model setup.

- **Step 2**: Construct the piece-wise continuous vertical subgrid
  profiles of tracer mixing ratios (:math:`q`), zonal and meridional
  winds (*u* and *v*), and total energy (:math:`\Gamma`) in the
  Lagrangian control-volume coordinate based on the Piece-wise
  Parabolic Method (PPM, ). The total energy :math:`\Gamma` is
  computed as the sum of the finite-volume integrated geopotential
  :math:`\phi`, internal energy :math:`(C_{v}T)`, and the kinetic energy (:math:`K`) as follows:

  .. math::
     
     \label{int-t}
     \Gamma =\frac{1}{\delta p}\int [ C_{v}T+\phi +\frac{1}{2}(u^{2}+v^{2}) ] dp .
     
  Applying integration by parts and the ideal gas law, the above
  integral can be rewritten as

  .. math::

     \label{TE_fv}
     \Gamma =C_{p}\overline{T}+\frac{1}{\delta p}\delta ( p\phi ) +K ,

  where :math:`\overline{T}` is the layer mean temperature,
  :math:`K` is the kinetic energy, :math:`p` is the pressure at layer edges,
  and :math:`C_{v}` and :math:`C_{p}` are the specific heat of the
  air at constant volume and at constant pressure, respectively. Layer
  mean values of :math:`q`, (*u*, *v*), and :math:`\Gamma` in the
  Eulerian coordinate system are obtained by integrating analytically
  the sub-grid distributions, in the vertical direction, from model
  top to the surface, layer by layer. Since the hydrostatic pressure
  is chosen as the mapping coordinate, tracer mass, momentum, and
  total energy are locally and globally conserved.

- **Step 3**: Compute kinetic energy in the Eulerian coordinate system
  for each layer. Substituting kinetic energy and the hydrostatic
  relationship into ([TE:sub:`f`\ v]), the layer mean temperature
  :math:`\overline{T}_{k}` for layer :math:`k` in the Eulerian coordinate
  is then retrieved from the reconstructed total energy (done in Step
  2) by a fully explicit integration procedure starting from the
  surface up to the model top as follows:

  .. math::

     \label{map-t}
     \overline{T}_{k}=\frac{\Gamma _{k}-K_{k}-\phi
     _{k+\frac{1}{2}}}{C_{p}[ 1-\kappa \, p_{k-\frac{1}{2}}\frac{ln\,
     p_{k+\frac{1}{2}}-ln\,
     p_{k-\frac{1}{2}}}{p_{k+\frac{1}{2}}-p_{k-\frac{1}{2}}}] } .

To convert the potential temperature :math:`\Theta` to the layer mean
temperature the conversion factor is obtained by equating the following
two equivalent forms of the hydrostatic relation for :math:`\Theta` and :math:`\overline{T}:`

.. math::

   \label{hydro-PT}
   \delta \phi =-C_{p}\Theta \, \delta \Pi ,

.. math::

   \label{hydro-T}
   \delta \phi =-R\overline{T}\, \delta ln\, \Pi ,

where :math:`\Pi =p^{\kappa }`. The conversion formula between layer
mean temperature and layer mean potential temperature is obtained as
follows:

.. math::

   \label{convt}
   \Theta =\kappa \frac{\delta ln\Pi }{\delta \Pi }\overline{T} .

The physical implication of retrieving the layer mean temperature from
the total energy as described in Step 3 is that the dissipated kinetic
energy, if any, is locally converted into internal energy via the
vertically sub-grid mixing (dissipation) processes. Due to the
monotonicity preserving nature of the sub-grid reconstruction the
column-integrated kinetic energy inevitably decreases (dissipates),
which leads to local frictional heating. The frictional heating is a
physical process that maintains the conservation of the total energy in
a closed system.

As viewed by an observer riding on the Lagrangian surfaces, the mapping
procedure essentially performs the physical function of the
relative-to-the-Eulerian-coordinate vertical transport, by vertically
redistributing (air and tracer) mass, momentum, and total energy from
the Lagrangian control-volume back to the Eulerian framework.

As described in section [FVvdisc], the model time integration cycle
consists of :math:`m` small time steps for the 2D Lagrangian dynamics
and one large time step for tracer transport. The mapping time step can
be much larger than that used for the large-time-step tracer transport.
In tests using the Held-Suarez forcing , a three-hour mapping time
interval is found to be adequate. In the full model integration, one may
choose the same time step used for the physical parameterizations so as
to ensure the input state variables to physical parameterizations are in
the usual “Eulerian” vertical coordinate.

Adjustment of specific humidity to conserve water
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The physics parameterizations operate on a model state provided by the
dynamics, and are allowed to update specific humidity. However, the
surface pressure remains fixed throughout the physics updates, and since
there is an explicit relationship between the surface pressure and the
air mass within each layer, the total air mass must remain fixed as
well. This implies a change of dry air mass at the end of the physics
updates. We impose a restriction that dry air mass and water mass be
conserved as follows:

The total pressure :math:`p` is

.. math:: p = d + e .

with dry pressure :math:`d`, water vapor pressure :math:`e`. The
specific humidity is

.. math:: q = \frac{e}{p} = \frac{e}{d+e}, \qquad d = (1-q) p .

We define a layer thickness as :math:`\delta^k p \equiv p^{k+1/2} -
p^{k-1/2}`, so

.. math:: \delta^k d = (1-q^k)\delta^k p .

We are concerned about 3 time levels: :math:`q_n` is input to physics,
:math:`q_{n*}` is output from physics, :math:`q_{n+1}` is the adjusted
value for dynamics.

Dry mass is the same at :math:`n` and :math:`n+1` but not at :math:`n*`.
To conserve dry mass, we require that

.. math:: \delta^k d_n  = \delta^k d_{n+1}

or

.. math::

   (1-q^k_n)\delta^k p_n = (1-q^k_{n+1})\delta^k p_{n+1}
     \label{eq:drydif} .

Water mass is the same at :math:`n*` and :math:`n+1`, but not at
:math:`n`. To conserve water mass, we require that

.. math:: q^k_{n*} \delta^k p_n = q^k_{n+1}\delta^k p_{n+1} \label{eq:wetdif} .

Substituting ([eq:wetdif]) into ([eq:drydif]),

.. math:: (1-q^k_n)\delta^k p_n  = \delta^k p_{n+1} - q^k_{n*} \delta^k p_n

.. math:: \delta^k p_{n+1} = (1 - q^k_n +  q^k_{n*})\delta^k p_n

 which yields a modified specific humidity for the dynamics:

.. math::

   q^k_{n+1} = q^k_n \frac{\delta^k p_n}{\delta^k p_{n+1}} =
    \frac{q^k_{n*}}{1 - q^k_n +  q^k_{n*}} .

Further discussion
~~~~~~~~~~~~~~~~~~

There are still aspects of the numerical formulation in the finite
volume dynamical core that can be further improved. For example, the
choice of the horizontal grid, the computational efficiency of the
split-explicit time marching scheme, the choice of the various
monotonicity constraints, and how the conservation of total energy is
achieved.

The impact of the non-linear diffusion associated with the monotonicity
constraint is difficult to assess. All discrete schemes must address the
problem of subgrid-scale mixing. The finite-volume algorithm contains a
non-linear diffusion that mixes strongly when monotonicity principles
are locally violated. However, the effect of nonlinear diffusion due to
the imposed monotonicity constraint diminishes quickly as the resolution
matches better to the spatial structure of the flow. In other numerical
schemes, however, an explicit (and tunable) linear diffusion is often
added to the equations to provide the subgrid-scale mixing as well as to
smooth and/or stabilize the time marching.

The finite-volume dynamical core as implemented in CAM and described
here conserves the dry air and all other tracer mass exactly without a
“mass fixer”. The vertical Lagrangian discretization and the associated
remapping conserves the total energy exactly. The only remaining issue
regarding conservation of the total energy is the horizontal
discretization and the use of the “diffusive” transport scheme with
monotonicity constraint. To compensate for the loss of total energy due
to horizontal discretization, we apply a global fixer to add the loss in
kinetic energy due to “diffusion” back to the thermodynamic equation so
that the total energy is conserved. However, it should be noted that
even without the “energy fixer” the loss in total energy (in flux unit)
is found to be less than 2 :math:`(W/m^{2}`) with the 2 degrees
resolution, and much smaller with higher resolution. In the future, we
may consider using the total energy as a transported prognostic variable
so that the total energy could be automatically conserved.

Spectral Element Dynamical Core
-------------------------------

The CAM includes an optional dynamical core from HOMME, NCAR’s
High-Order Method Modeling Environment . The stand-alone HOMME is used
for research in several different types of dynamical cores. The
dynamical core incorporated into CAM4 uses HOMME’s continuous Galerkin
spectral finite element method , here abbreviated to the spectral
element method (SEM). This method is designed for fully unstructured
quadrilateral meshes. The current configurations in the CAM are based on
the cubed-sphere grid. The main motivation for the inclusion of HOMME is
to improve the scalability of the CAM by introducing quasi-uniform grids
which require no polar filters . HOMME is also the first dynamical core
in the CAM which locally conserves energy in addition to mass and
two-dimensional potential vorticity .

HOMME represents a large change in the horizontal grid as compared to
the other dynamical cores in CAM. Almost all other aspects of HOMME are
based on a combination of well-tested approaches from the Eulerian and
FV dynamical cores. For tracer advection, HOMME is modeled as closely as
possible on the FV core. It uses the same conservation form of the
transport equation and the same vertically Lagrangian discretization .
The HOMME dynamics are modeled as closely as possible on Eulerian core.
They share the same vertical coordinate, vertical discretization,
hyper-viscosity based horizontal diffusion, top-of-model dissipation,
and solve the same moist hydrostatic equations. The main differences are
that HOMME advects the surface pressure instead of its logarithm (in
order to conserve mass and energy), and HOMME uses the vector-invariant
form of the momentum equation instead of the vorticity-divergence
formulation. Several dry dynamical cores including HOMME are evaluated
in using a grid-rotated version of the baroclinic instability test case.

The timestepping in HOMME is a form of dynamics/tracer/physics
subcycling, achieved through the use of multi-stage 2nd order accurate
Runge-Kutta methods. The tracers and dynamics use the same timestep
which is controlled by the maximum anticipated wind speed, but the
dynamics uses more stages than the tracers in order to maintain
stability in the presence of gravity waves. The forcing is applied using
a time-split approach. The optimal forcing strategy in HOMME has not yet
been determined, so HOMME supports several options. The first option is
modeled after the FV dynamical core and the forcing is applied as an
adjustment at each physics timestep. The second option is to convert all
forcings into tendencies which are applied at the end of each
dynamics/tracer timestep. If the physics timestep is larger than the
tracer timestep, then the tendencies are held fixed and only updated at
each physics timestep. Finally, a hybrid approach can be used where the
tracer tendencies are applied as in the first option and the dynamics
tendencies are applied as in the second option.

Continuum Formulation of the Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

HOMME uses a conventional vector-invariant form of the moist primitive
equations. For the vertical discretization it uses the hybrid
:math:`\eta` pressure vertical coordinate system modeled after
[eul:terrain] The formulation here differs only in that surface pressure
is used as a prognostic variable as opposed to its logarithm.

In the :math:`\eta`-coordinate system, the pressure is given by

.. math:: p(\eta) = A(\eta) p_0 + B(\eta) p_s.

The hydrostatic approximation
:math:`\partial p / \partial z  = - g \rho` is used to replace the
mass density :math:`\rho` by an :math:`\eta`-coordinate pseudo-density
:math:`\partial p / \partial \eta`. The material derivative in
:math:`\eta`-coordinates can be written (e.g. , Sec.3.3),

.. math:: \frac{ D X }{D t} = {{\frac{\partial {X}}{\partial t}}} + {{{\smash[t]{\vec{u}}}}}\cdot {\nabla}X + {{\dot\eta}}{\frac{\partial {X}}{\partial \eta}}

where the :math:`{\nabla}()` operator (as well as
:math:`{\nabla\cdot}()` and :math:`{{{\nabla}\times}}()` below) is the
two-dimensional gradient on constant :math:`\eta`-surfaces,
:math:`\partial / \partial \eta` is the vertical derivative,
:math:`{{\dot\eta}}= D\eta/Dt` is a vertical flow velocity and
:math:`{{{\smash[t]{\vec{u}}}}}` is the horizontal velocity component
(tangent to constant :math:`z`-surfaces, not :math:`\eta`-surfaces).

The :math:`\eta` -coordinate atmospheric primitive equations, neglecting
dissipation and forcing terms can then be written as

.. math::
   :label: E:PEmom-E:PEtemp-E:PEcont-E:PEq

   {{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial t}}} + ( {{\mathbf{\zeta}}}+ f ) {\hat{k}}{{\times}}{{{\smash[t]{\vec{u}}}}}+
   {\nabla}( \frac12 {{{\smash[t]{\vec{u}}}}}^2 + \Phi )  +
   {{\dot\eta}}{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial \eta}} + \frac{RT_v}{p} {\nabla}p & = 0  \\
   {{\frac{\partial {T}}{\partial t}}} + {{{\smash[t]{\vec{u}}}}}\cdot {\nabla}T + {{\dot\eta}}{\frac{\partial {T}}{\partial \eta}} - 
   \frac{RT_v}{c^*_p p} \omega  & = 0   \\
   {{\frac{\partial {}}{\partial t}}}({\frac{\partial {p}}{\partial \eta}}) + {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta}}{{{\smash[t]{\vec{u}}}}}) + 
   {\frac{\partial {}}{\partial \eta}} ( {{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}) & = 0  \\
   {{\frac{\partial {}}{\partial t}}} ( {\frac{\partial {p}}{\partial \eta}}q ) +  {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta}}q {{{\smash[t]{\vec{u}}}}}) + 
   {\frac{\partial {}}{\partial \eta}} ( {{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}q ) & = 0. 

These are prognostic equations for :math:`{{{\smash[t]{\vec{u}}}}}`, the
temperature :math:`T`, density
:math:`{\frac{\partial {p}}{\partial \eta}}`, and
:math:`{\frac{\partial {p}}{\partial \eta}}q` where :math:`q` is the
specific humidity. The prognostic variables are functions of time
:math:`t`, vertical coordinate :math:`\eta` and two coordinates
describing the surface of the sphere. The unit vector normal to the
surface of the sphere is denoted by :math:`{\hat{k}}`. This formulation
has already incorporated the hydrostatic equation and the ideal gas law,
:math:`p = \rho R T_v`. There is a no-flux (:math:`{{\dot\eta}}= 0`)
boundary condition at :math:`\eta=1` and :math:`\eta=\eta_\text{top}`.
The vorticity is denoted by
:math:`\zeta = {\hat{k}}\cdot {{{\nabla}\times}}{{{\smash[t]{\vec{u}}}}}`,
:math:`f` is a Coriolis term and :math:`\omega = Dp/Dt` is the pressure
vertical velocity. The virtual temperature :math:`T_v` and
variable-of-convenience :math:`c^*_p` are defined as in [eul:terrain].

The diagnostic equations for the geopotential height field :math:`\Phi` is

.. math::
   :label: E:hydrostatic

   \Phi = \Phi_s + \int_{\eta}^{1} \frac{R T_v }{p} {\frac{\partial {p}}{\partial \eta}}\, d\eta

where :math:`\Phi_s` is the prescribed surface geopotential height
(given at :math:`\eta=1`). To complete the system, we need diagnostic
equations for :math:`{{\dot\eta}}` and :math:`\omega`, which come from
integrating with respect to :math:`\eta`. In fact, can be replaced by a
diagnostic equation for
:math:`{{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}` and a
prognostic equation for surface pressure :math:`p_s`

.. math::
   :label: E:PEcont2a

   & {{\frac{\partial {}}{\partial t}}}p_s +  \int_{\eta_\text{top}}^{1} {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta}}{{{\smash[t]{\vec{u}}}}}) \, d\eta = 0     \\
   & {{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}= - {{\frac{\partial {p}}{\partial t}}} - \int_{\eta_\text{top}}^\eta {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta'}}{{{\smash[t]{\vec{u}}}}}) \, d\eta', 

where is evaluated at the model bottom (:math:`\eta=1`) after using that
:math:`\partial p / \partial t = B(\eta) \partial p_s / \partial t` and
:math:`{{\dot\eta}}(1)=0, B(1)=1`. Using Eq [E:PEcont2c], we can derive
a diagnostic equation for the pressure vertical velocity
:math:`\omega = Dp/Dt`,

.. math::

   \omega =   {{\frac{\partial {p}}{\partial t}}} +  {{{\smash[t]{\vec{u}}}}}\cdot {\nabla}p + {{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}=
   {{{\smash[t]{\vec{u}}}}}\cdot {\nabla}p  - \int_{\eta_\text{top}}^\eta {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta}}{{{\smash[t]{\vec{u}}}}}) \, d\eta'

Finally, we rewrite as

.. math::
   :label: E:PEcont2

   {{\dot\eta}}{\frac{\partial {p}}{\partial \eta}}= B(\eta) \int_{\eta_\text{top}}^{1} {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta}}{{{\smash[t]{\vec{u}}}}}) \, d\eta
   - \int_{\eta_\text{top}}^\eta {\nabla\cdot}( {\frac{\partial {p}}{\partial \eta'}}{{{\smash[t]{\vec{u}}}}}) \, d\eta',

Conserved Quantities
~~~~~~~~~~~~~~~~~~~~

The equations have infinitely many conserved quantities, including mass,
tracer mass, potential temperature defined by

.. math:: 

   M_X =  \iint {\frac{\partial {p}}{\partial \eta}}X \, d\eta {d\cal{A}}

with (:math:`X = 1, q` or :math:`(p/p_0)^{-\kappa} T`) and the total
moist energy :math:`E` defined by

.. math::

   E =  \iint {\frac{\partial {p}}{\partial \eta}}( \frac12 {{{\smash[t]{\vec{u}}}}}^2 + c_p^* T  ) \, d\eta {d\cal{A}}+  \int p_s \Phi_s \, {d\cal{A}} \label{E:E1}

where :math:`{d\cal{A}}` is the spherical area measure. To compute
these quantities in their traditional units they should be divided by
the constant of gravity :math:`g`. We have omitted this scaling since
:math:`g` has also been scaled out from –. We note that in this
formulation of the primitive equations, the pressure :math:`p` is a
moist pressure, representing the effects of both dry air and water
vapor. The unforced equations conserve both the moist air mass
(:math:`X=1` above) and the dry air mass (:math:`X=1-q` ). However, in
the presence of a forcing term in (representing sources and sinks of
water vapor as would be present in a full model) a corresponding forcing
term must be added to to ensure that dry air mass is conserved.

The energy is specific to the hydrostatic equations. We have omitted
terms from the physical total energy which are constant under the
evolution of the unforced hydrostatic equations . It can be converted
into a more universal form involving
:math:`\tfrac12 {{{\smash[t]{\vec{u}}}}}^2 + c^*_v T + \Phi`, with
:math:`c^*_v` defined similarly to :math:`c^*_p`, so that
:math:`c^*_v = c_v + (c_{vv}-c_v) q` where :math:`c_v` and
:math:`c_{vv}` are the specific heats of dry air and water vapor defined
at constant volume. We note that :math:`c_p = R + c_v` and
:math:`c_{pv} =  R_v + c_{vv}` so that
:math:`c_p^* T = c_v^* T + R T_v`. Expanding :math:`c_p^* T` with this
expression, integrating by parts with respect to :math:`\eta` and making
use of the fact that the model top is at a constant pressure

.. math::

   \int {\frac{\partial {p}}{\partial \eta}}R T_v   \, d\eta = -\int p \frac{\partial \Phi}{\partial \eta}   \, d\eta =  \int {\frac{\partial {p}}{\partial \eta}}\Phi    \, d\eta  - (  p \Phi  )   \Big| ^{\eta=1}_{\eta=\eta_\text{top}}

and thus

.. math::

   E =  \iint {\frac{\partial {p}}{\partial \eta}}( \frac12 {{{\smash[t]{\vec{u}}}}}^2 + c^*_v T  + \Phi ) \, d\eta {d\cal{A}}+   \int  p_\text{top} \Phi(\eta_\text{top})  \,  {d\cal{A}}.    \label{E:E2}

The model top boundary term in vanishes if :math:`p_\text{top}=0`.
Otherwise it must be included to be consistent with the hydrostatic
equations. It is present due to the fact that the hydrostatic momentum
equation neglects the vertical pressure gradient.

Horizontal Discretization: Functional Spaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the finite element method, instead of constructing discrete
approximations to derivative operators, one constructs a discrete
functional space, and then finds the function in this space which solves
the equations of interest in a minimum residual sense. As compared to
finite volume methods, there is less choice in how one constructs the
discrete derivative operators in this setting, since functions in the
discrete space are represented in terms of known basis functions whose
derivatives are known, often analytically.

Let :math:`x^\alpha` and :math:`{{\smash[t]{\vec{x}}}}=x^1{\smash[t]{\vec{e}}}_1 + x^2{\smash[t]{\vec{e}}}_2`
be the Cartesian coordinates and position vector of a point in the
reference square :math:`{{[-1,1]^2}}` and let :math:`r^\alpha` and
:math:`{{\smash[t]{\vec{r}}}}` be the coordinates and position vector of
a point on the surface of the sphere, denoted by :math:`\Omega`. We mesh
:math:`\Omega` using the cubed-sphere grid (Fig. [f:sphere4]) first used
in . Each cube face is mapped to the surface of the sphere with the
equal-angle gnomonic projection . The map from the reference element
:math:`[-1,1]^2` to the cube face is a translation and scaling. The
composition of these two maps defines a :math:`{\cal C^{1}}` map
from the spherical elements to the reference element
:math:`{{[-1,1]^2}}`. We denote this map and its inverse by

.. math::
   :label: e-map

   {{\smash[t]{\vec{r}}}}= {{\smash[t]{\vec{r}}}}({{\smash[t]{\vec{x}}}};m), \qquad  {{\smash[t]{\vec{x}}}}= {{\smash[t]{\vec{x}}}}({{\smash[t]{\vec{r}}}};m).

[f:sphere4] Tiling the surface of the sphere with
quadrilaterals. An inscribed cube is projected to the surface of the
sphere. The faces of the cubed sphere are further subdivided to form a
quadrilateral grid of the desired resolution. Coordinate lines from the
gnomonic equal-angle projection are shown. 

We now define the discrete space used by the SEM. First we denote the
space of polynomials up to degree :math:`d` in :math:`{{[-1,1]^2}}` by

.. math::

   {\cal P_d}=  {\mathop{\mathrm{span}}}_{i,j=0}^d (x^1)^i (x^2)^j  =  
   {\mathop{\mathrm{span}}}\limits_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}}\phi_{{\smash[t]{\vec{\imath}}}}({{\smash[t]{\vec{x}}}}),

where :math:`\mathbb{I} = \{0,\ldots,d\}^2` contains all the degrees and
:math:`\phi_{{\smash[t]{\vec{\imath}}}}({{\smash[t]{\vec{x}}}})= \varphi_{i^1}(x^1) \varphi_{i^2}(x^2)`,
:math:`i^\alpha=0,\dots,d`, are the cardinal functions, namely
polynomials that interpolate the tensor-product of degree-\ :math:`d`
Gauss-Lobatto-Legendre (GLL) nodes
:math:`{{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}} =  \xi_{i^1}{\smash[t]{\vec{e}}}_1 + \xi_{i^2}{\smash[t]{\vec{e}}}_2`.
The GLL nodes used within an element for :math:`d=3` are shown in
Fig. [f:GLLnodes]. The cardinal-function expansion coefficients of a
function :math:`g` are its GLL nodal values, so we have

.. math::
   :label: e:cfvec

   g({{\smash[t]{\vec{x}}}})= \sum_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}} g({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}})  \phi_{{\smash[t]{\vec{\imath}}}}({{\smash[t]{\vec{x}}}}).

We can now define the piecewise-polynomial SEM spaces
:math:`{\cal V^{0}_{}}` and :math:`{\cal V^{1}_{}}` as

.. math::
   :label: e:Hzero

   {\cal V^{0}_{}}& =  \{f \in{\cal L^2}(\Omega) :  f({{\smash[t]{\vec{r}}}}(\cdot;m)) \in {\cal P_d}, \forall m\}
   = {\mathop{\mathrm{span}}}_{m=1}^M\{\phi_{{\smash[t]{\vec{\imath}}}}({{\smash[t]{\vec{x}}}}(\cdot;m))\}_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}} \\
    \text{and}\qquad{\cal V^{1}_{}}& = {\cal C^{0}}(\Omega)\cap{\cal V^{0}_{}}.
    \notag

Functions in :math:`{\cal V^{0}_{}}` are polynomial within each
element but may be discontinuous at element boundaries and
:math:`{\cal V^{1}_{}}` is the subspace of continuous function in
:math:`{\cal V^{0}_{}}`. We take
:math:`M_d = \dim {\cal V^{0}_{}} = (d+1)^3 M`, and
:math:`L = \dim {\cal V^{1}_{}} <  M_d`. We then construct a set of
:math:`L` unique points by

.. math::

   \{{{\smash[t]{\vec{r}}}}_\ell\}_{\ell=1}^L = \bigcup_{m=1}^M{{\smash[t]{\vec{r}}}}(\{{{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}}\}_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}};m),
   \label{e:GlInt}

For every point :math:`{{\smash[t]{\vec{r}}}}_\ell`, there exists at
least one element :math:`\Omega_m` and at least one GLL node
:math:`{{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}}={{\smash[t]{\vec{x}}}}({{\smash[t]{\vec{r}}}}_\ell;m)`.
In 2D, if :math:`{{\smash[t]{\vec{r}}}}_\ell` belongs to exactly one
:math:`\Omega_m` it is an element-interior node. If it belongs to
exactly two :math:`\Omega_m`\ s, it is an element-edge interior node.
Otherwise it is a vertex node.

[f:GLLnodes] A :math:`4 \times 4` tensor product grid of GLL
nodes used within each element, for a degree :math:`d=3` discretization.
Nodes on the boundary are shared by neighboring elements. |

We also define similar spaces for 2D vectors. We introduce two families
of spaces, with a subscript of either *con* or *cov*, denoting if the
contravariant or covariant components of the vectors are piecewise
polynomial, respectively.

.. math::

   {{\cal V^{0}_{\rm con}}} & =  \{{{{\smash[t]{\vec{u}}}}}\in{\cal L^2}(\Omega)^2 :  u^\alpha \in {\cal V^{0}_{}},\;\alpha=1,2\} \\
   \text{and}\qquad{{\cal V^{1}_{\rm con}}} & = {\cal C^{0}}(\Omega)^2\cap{{\cal V^{0}_{\rm con}}},

where :math:`u^1, u^2` are the contravariant components of
:math:`{{{\smash[t]{\vec{u}}}}}` defined below. Vectors in
:math:`{{\cal V^{1}_{\rm con}}}` are globally continuous and their
contravariant components are polynomials in each element. Similarly,

.. math::

   {{{\cal V^{0}_{\rm cov}}}}& =  \{{{{\smash[t]{\vec{u}}}}}\in{\cal L^2}(\Omega)^2 :  u_\beta \in {\cal V^{0}_{}},\;\beta=1,2\}
   \\
   \text{and}\qquad{{{\cal V^{1}_{\rm cov}}}}& = {\cal C^{0}}(\Omega)^2\cap{{{\cal V^{0}_{\rm cov}}}}.

The SEM is a Galerkin method with respect to the
:math:`{\cal V^{1}_{}}` subspace and it can be formulated solely in
terms of functions in :math:`{\cal V^{1}_{}}`. In CAM-HOMME, the
typical configuration is to run with :math:`d=3` which achieves a 4th
order accurate horizontal discretization . All variables in the
CAM-HOMME initial condition and history files as well as variables
passed to the physics routines are represented by their grid point
values at the points
:math:`\{{{\smash[t]{\vec{r}}}}_\ell\}_{\ell=1}^L`. However, for some
intermediate quantities and internally in the dynamical core it is
useful to consider the larger :math:`{\cal V^{0}_{}}` space, where
variables are represented by their grid point values at the :math:`M_d`
mapped GLL nodes. This later representation can also be considered as
the cardinal-function expansion of a function :math:`f` local to each
element,

.. math::

   f({{\smash[t]{\vec{r}}}})  =
   \sum_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}} f({{\smash[t]{\vec{r}}}}({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}};m)) \phi_{{\smash[t]{\vec{\imath}}}}({{\smash[t]{\vec{x}}}}({{\smash[t]{\vec{r}}}};m))
   \label{e:localexpand}

since the expansion coefficients are the function values at the mapped
GLL nodes. Functions :math:`f` in :math:`{\cal V^{0}_{}}` can be
multiple-valued at GLL nodes that are *redundant* (i.e., shared by more
than one element), while for :math:`f \in {\cal V^{1}_{}}`, the
values at any redundant points must all be the same.

Horizontal Discretization: Differential Operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use the standard curvilinear coordinate formulas for vector operators
following . Given the :math:`2\times2` Jacobian of the the mapping from
:math:`{{[-1,1]^2}}` to :math:`\Omega_m`, we denote its
determinant-magnitude by

.. math::

   {J}=| {\frac{\partial {{{\smash[t]{\vec{r}}}}}}{\partial {{\smash[t]{\vec{x}}}}}} |.
   \label{e:Jd}

A vector :math:`{{{\smash[t]{\vec{v}}}}}` may be written in terms of
physical or covariant or contravariant components,
:math:`{\newcommand\arraystretch{.2}v[\begin{array}{@{}l@{}}\vphantom{\scriptstyle\alpha}\\\scriptstyle\gamma\\\vphantom{\scriptstyle\alpha}\end{array}]}`
or :math:`v_\beta` or :math:`v^\alpha`,

.. math::

   {{{\smash[t]{\vec{v}}}}}=\sum_{\gamma=1}^3{\newcommand\arraystretch{.2}v[\begin{array}{@{}l@{}}\vphantom{\scriptstyle\alpha}\\\scriptstyle\gamma\\\vphantom{\scriptstyle\alpha}\end{array}]}{\frac{\partial {{{\smash[t]{\vec{r}}}}}}{\partial r^\gamma}}=
   \sum_{\beta=1}^3v_\beta{{\smash[t]{\vec{g}}}}^\beta=
   \sum_{\alpha=1}^3v^\alpha{{\smash[t]{\vec{g}}}}_\alpha,
   \label{e:3comps}

that are related by
:math:`v_\beta={{{\smash[t]{\vec{v}}}}}{\mathbin{\mathbf{\cdot}}}{{\smash[t]{\vec{g}}}}_\beta`
and
:math:`v^\alpha={{{\smash[t]{\vec{v}}}}}{\mathbin{\mathbf{\cdot}}}{{\smash[t]{\vec{g}}}}^\alpha`,
where :math:`{{\smash[t]{\vec{g}}}}^\alpha={\nabla}x^\alpha` is a
contravariant basis vector and
:math:`{{\smash[t]{\vec{g}}}}_\beta={\frac{\partial {{{\smash[t]{\vec{r}}}}}}{\partial x^\beta}}`
is a covariant basis vector.

The dot product and contravariant components of the cross product are

.. math::

   {{{\smash[t]{\vec{u}}}}}{\mathbin{\mathbf{\cdot}}}{{{\smash[t]{\vec{v}}}}}=  \sum_{\alpha=1}^3 u_\alpha v^\alpha\qquad\text{and}\qquad
   ( {{{\smash[t]{\vec{u}}}}}{{\times}}{{{\smash[t]{\vec{v}}}}})^\alpha  = 
   \frac 1 {J}\sum_{\beta,\gamma=1}^3\epsilon^{\alpha \beta \gamma } u_\beta v_\gamma
   \label{e:dotandcross}

where :math:`\epsilon^{\alpha\beta\gamma}\in\{0,\pm1\}` is the
Levi-Civita symbol. The divergence, covariant coordinates of the
gradient and contravariant coordinates of the curl are

.. math::

   {\nabla\cdot}{{{\smash[t]{\vec{v}}}}}= \frac 1 {J}\sum_\alpha{\frac{\partial {}}{\partial x^\alpha}}({J}v^\alpha),
   \quad
   ( {\nabla}f )_\alpha = {\frac{\partial {f}}{\partial x^\alpha}}
   \quad\text{and}\quad
   ( {{{\nabla}\times}}{{{\smash[t]{\vec{v}}}}})^\alpha =  \frac 1 {J}\sum_{\beta,\gamma}
   \epsilon^{\alpha \beta \gamma} {\frac{\partial {v_\gamma}}{\partial x^\beta}}.
   \label{e:divgradcurl}

In the SEM, these operators are all computed in terms of the
derivatives with respect to :math:`{{\smash[t]{\vec{x}}}}` in the
reference element, computed exactly (to machine precision) by
differentiating the local element expansion . For the gradient, the
covariant coordinates of :math:`{\nabla}f, f \in
{\cal V^{0}_{}}` are thus computed exactly within each element. Note
that :math:`{\nabla}f \in {{{\cal V^{0}_{\rm cov}}}}`, but may not
be in :math:`{{{\cal V^{1}_{\rm cov}}}}` even for
:math:`f \in {\cal V^{1}_{}}` due to the fact that its components
will be multi-valued at element boundaries because :math:`{\nabla}f`
computed in adjacent elements will not necessarily agree along their
shared boundary. In the case where :math:`{J}` is constant within each
element, the SEM curl of
:math:`{{{\smash[t]{\vec{v}}}}}\in {{{\cal V^{0}_{\rm cov}}}}` and
the divergence of
:math:`{{{\smash[t]{\vec{u}}}}}\in {{\cal V^{0}_{\rm con}}}` will
also be exact, but as with the gradient, multiple-valued at element
boundaries.

For non-constant :math:`{J}`, these operators may not be computed
exactly by the SEM due to the Jacobian factors in the operators and the
Jacobian factors that appear when converting between covariant and
contravariant coordinates. We follow and evaluate these operators in the
form shown in . The quadratic terms that appear are first projected into
:math:`{\cal V^{0}_{}}` via interpolation at the GLL nodes and then
this interpolant is differentiated exactly using . For example, to
compute the divergence of
:math:`{{{\smash[t]{\vec{v}}}}}\in {{\cal V^{0}_{\rm con}}}`, we
first compute the interpolant
:math:`{\mathcal I}({J}v^\alpha)\in{\cal V^{0}_{}}` of
:math:`{J}v^\alpha`, where the GLL interpolant of a product :math:`fg`
derives simply from the product of the GLL nodal values of :math:`f` and
:math:`g`. This operation is just a reinterpretation of the nodal values
and is essentially free in the SEM. The derivatives of this interpolant
are then computed exactly from . The sum of partial derivatives are then
divided by :math:`{J}` at the GLL nodal values and thus the SEM
divergence operator :math:`{{\nabla_{\rm h} \cdot }}()` is given by

.. math::

   {\nabla\cdot}{{{\smash[t]{\vec{v}}}}}\approx {{\nabla_{\rm h} \cdot }}{{{\smash[t]{\vec{v}}}}}= {\mathcal I}( 
   \frac 1 {J}\sum_\alpha{\frac{\partial {{\mathcal I}( {J}v^\alpha)}}{\partial x^\alpha}} 
   )  \in {\cal V^{0}_{}}.
   \label{e:SEMdiv}

Similarly, the gradient and curl are approximated by

.. math::

    ( {\nabla}f )_\alpha  \approx
   ( {{\nabla_{\rm h}}}f )_\alpha & = 
   {\frac{\partial {f}}{\partial x^\alpha}} 
   \label{e:hgrad}\\
   \text{and}\qquad
   ( {{{\nabla}\times}}{{{\smash[t]{\vec{v}}}}})^\alpha \approx
   ( {{{\nabla_{\rm h}}}{{\times}}}{{{\smash[t]{\vec{v}}}}})^\alpha  &  = 
   \sum_{\beta,\gamma}
   \epsilon^{\alpha \beta \gamma} {\mathcal I}(\frac 1 {J}{\frac{\partial {v_\gamma}}{\partial x^\beta}})
   \label{e:hcurl}

with :math:`{{\nabla_{\rm h}}}f \in{{{\cal V^{0}_{\rm cov}}}}` and
:math:`{{{\nabla_{\rm h}}}{{\times}}}{{{\smash[t]{\vec{v}}}}}\in{{\cal V^{0}_{\rm con}}}`.
The SEM is well known for being quite efficient in computing these types
of operations. The SEM divergence, gradient and curl can all be
evaluated at the :math:`(d+1)^3` GLL nodes within each element in
:math:`\mathcal{O}(d)` operations per node using the tensor-product
property of these points .

Horizontal Discretization: Discrete Inner-Product
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of using exact integration of the basis functions as in a
traditional finite-element method, the SEM uses a GLL quadrature
approximation for the integral over :math:`\Omega`, that we denote by
:math:`{\langle}\cdot {\rangle}`. We can write this integral as a sum of
area-weighted integrals over the set of elements
:math:`\{ \Omega_m \}_{m=1}^M` used to decompose the domain,

.. math:: \int f g \,{d\cal{A}}= \sum_{m=1}^M\int_{\Omega_m} f g \, {d\cal{A}}.

The integral over a single element :math:`\Omega_m` is written as an
integral over :math:`{{[-1,1]^2}}` by

.. math::

   \int_{\Omega_m} f g \, {d\cal{A}}= \iint_{{{[-1,1]^2}}}f({{\smash[t]{\vec{r}}}}(\cdot;m))g({{\smash[t]{\vec{r}}}}(\cdot;m)){J}_m \,d x^1 \,d x^2 
   \approx{\langle}fg{\rangle}_{\Omega_m},
   \label{e:el2sq}

where we approximate the integral over :math:`{{[-1,1]^2}}` by GLL
quadrature,

.. math::

   {\langle}f g {\rangle}_{\Omega_m} =
   \sum_{{\smash[t]{\vec{\imath}}}\in\mathbb{I}}w_{i^1}w_{i^2}{J}_m({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}})
   f({{\smash[t]{\vec{r}}}}({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}};m))g({{\smash[t]{\vec{r}}}}({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}};m))
   \label{E:intomega}

The SEM approximation to the global integral is then naturally defined
as

.. math::

   \int f g \,{d\cal{A}}\approx  \sum_{m=1}^M {\langle}f g {\rangle}_{\Omega_m}
   ={\langle}fg{\rangle}\label{e:dip}

When applied to the product of functions
:math:`f,g \in {\cal V^{0}_{}}`, the quadrature approximation
:math:`{\langle}f g {\rangle}` defines a discrete inner-product in the
usual manner.

Horizontal Discretization: The Projection Operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`{{P}}: {\cal V^{0}_{}} arrow {\cal V^{1}_{}}`
be the unique orthogonal (self-adjoint) projection operator from
:math:`{\cal V^{0}_{}}` onto :math:`{\cal V^{1}_{}}` w.r.t. the
SEM discrete inner product . The operation :math:`{{P}}` is essentially
the same as the common procedure in the SEM described as *assembly* , or
*direct stiffness summation* . Thus the SEM assembly procedure is not an
ad-hoc way to remove the redundant degrees of freedom in
:math:`{\cal V^{0}_{}}`, but is in fact the natural projection
operator :math:`{{P}}`. Applying the projection operator in a finite
element method requires inverting the finite element mass matrix. A
remarkable fact about the SEM is that with the GLL based discrete inner
product and the careful choice of global basis functions, the mass
matrix is diagonal . The resulting projection operator then has a very
simple form: at element interior points, it leaves the nodal values
unchanged, while at element boundary points shared by multiple elements
it is a Jacobian-weighted average over all redundant values .

To apply the projection
:math:`{{P}}: {{{\cal V^{0}_{\rm cov}}}} arrow {{{\cal V^{1}_{\rm cov}}}}`
to vectors :math:`{{{\smash[t]{\vec{u}}}}}`, one cannot project the
covariant components since the corresponding basis vectors
:math:`{{\smash[t]{\vec{g}}}}_\beta` and
:math:`{{\smash[t]{\vec{g}}}}^\alpha` do not necessarily agree along
element faces. Instead we must define the projection as acting on the
components using a globally continuous basis such as the
latitude-longitude unit vectors :math:`\hat\theta` and
:math:`\hat\lambda`,

.. math::

   {{P}}({{{\smash[t]{\vec{u}}}}}) = 
   {{P}}( {{{\smash[t]{\vec{u}}}}}\cdot \hat\lambda) \hat\lambda 
   +
   {{P}}( {{{\smash[t]{\vec{u}}}}}\cdot \hat\theta) \hat\theta.

Horizontal Discretization: Galerkin Formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SEM solves a Galerkin formulation of the equations of interest.
Given the discrete differential operators described above, the primitive
equations can be written as an ODE for a generic prognostic variable
:math:`U` and right-hand-side (RHS) terms

.. math:: {{\frac{\partial {U}}{\partial t}}} = {\rm RHS}.

The SEM solves this equation in integral form with respect to the SEM
inner product. That is, for a :math:`\rm{RHS} \in {\cal V^{0}_{}}`,
the SEM finds the unique
:math:`{{\frac{\partial {U}}{\partial t}}} \in {\cal V^{1}_{}}` such
that

.. math:: {\langle}\phi  {{\frac{\partial {U}}{\partial t}}}  {\rangle}= {\langle}\phi \, {\rm RHS} {\rangle}\qquad \forall \phi \in {\cal V^{1}_{}}.

As the prognostic variable is assumed to belong to
:math:`{\cal V^{1}_{}}`, the RHS will in general belong to
:math:`{\cal V^{0}_{}}` since it contains derivatives of the
prognostic variables, resulting in the loss of continuity at the element
boundaries. If one picks a suitable basis for
:math:`{\cal V^{1}_{}}`, this discrete integral equation results in
a system of :math:`L` equations for the :math:`L` expansion coefficients
of :math:`{{\frac{\partial {U}}{\partial t}}}`. The SEM solves these
equations exactly, and the solution can be written in terms of the SEM
projection operator as

.. math:: {{\frac{\partial {U}}{\partial t}}} = {{P}}(  {\rm RHS} ).

The projection operator commutes with any time-stepping scheme, so the
equations can be solved in a two step process, illustrated here for
simplicity with the forward Euler method

-  Step 1:

   .. math:: U^* = U^t + \Delta t  \, {\rm RHS} \qquad U^* \in {\cal V^{0}_{}}

-  Step 2:

   .. math:: U^{t+1} = {{P}}( U^* )  \qquad U^{t+1} \in {\cal V^{1}_{}}

For compactness of notation, we will denote this two step procedure in
what follows by

.. math:: {{P}}^{-1} {{\frac{\partial {U}}{\partial t}}} =   {\rm RHS}.

Note that :math:`{{P}}` maps a :math:`M_d` dimensional space
:math:`{\cal V^{0}_{}}` into a :math:`L` dimensional space
:math:`{\cal V^{1}_{}}`, so here :math:`{{P}}^{-1}` denotes the left
inverse of :math:`{{P}}`. This inverse will never be computed, it is
only applied as in step 2 above.

This two step Galerkin solution process represents a natural separation
between computation and communication for the implementation of the SEM
on a parallel computer. The computations in step 1 are all local to the
data contained in a single element. Assuming an element-based
decomposition so that each processor contains at least one element, no
inter-processor communication is required in step 1. All inter-processor
communication in HOMME is isolated to the projection operator step, in
which element boundary data must be exchanged between adjacent elements.

Vertical Discretization
~~~~~~~~~~~~~~~~~~~~~~~

The vertical coordinate system uses a Lorenz staggering of the variables
as shown in [figure:1]. Let :math:`K` be the total number of layers,
with variables :math:`{{{\smash[t]{\vec{u}}}}}, T, q, \omega, \Phi` at
layer mid points denoted by :math:`k=1,2,\dots,K`. We denote layer
interfaces by :math:`k+\tfrac12,  k=0,1,\dots,K`, so that
:math:`\eta_{1/2}=\eta_\text{top}` and :math:`\eta_{K+1/2}=1`. The
:math:`\eta`-integrals will be replaced by sums. We will use
:math:`{ \mathop{\delta_\eta}}` to denote the discrete
:math:`\partial / \partial \eta` operator. The
:math:`{ \mathop{\delta_\eta}}` operator uses centered differences to
compute derivatives with respect to :math:`\eta` at layer mid point from
layer interface values,
:math:`{ \mathop{\delta_\eta}}(X)_k = (X_{k+1/2} - X_{k-1/2})/(\eta_{k+1/2}-\eta_{k-1/2})`.
We will use the over-bar notation for vertical averaging,
:math:`\overline q_{k+1/2} = (q_{k+1}+q_k)/2`. We also introduce the
symbol :math:`{\pi}` to denote the discrete pseudo-density
:math:`{\frac{\partial {p}}{\partial \eta}}` given by

.. math:: {\pi}_{k} = { \mathop{\delta_\eta}}( p )_k

.

We will use :math:`{ \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}` to
denote the discrete form of the
:math:`{{\dot\eta}}\partial / \partial\eta` operator. We use the
discretization given in [eul:econserve]. This operator acts on
quantities defined at layer mid-points and returns a result also at
layer mid-points,

.. math::

   { \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}(X)_k = 
   \frac1{2 {\pi}_k \Delta\eta_k}
   [ 
   ( {{\dot\eta}}{\pi})_{k+1/2}  (  X_{k+1} - X_k  )
   + 
   ( {{\dot\eta}}{\pi})_{k-1/2} (  X_{k} - X_{k-1} ) 
   ]
   \label{E:dnhprimealt}

where :math:`\Delta\eta_k = \eta_{k+1/2} - \eta_{k-1/2}`. We use the
over-bar notation since the formula can be seen as a
:math:`{\pi}`-weighted average of a layer interface centered difference
approximation to :math:`{{\dot\eta}}\partial/\partial \eta`. This
formulation was constructed in in order to ensure mass and energy
conservation. Here we will use an equivalent expression that can be
written in terms of :math:`{ \mathop{\delta_\eta}}`,

.. math::

   { \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}(X)_k = 
   \frac1{{\pi}_k} \Big[ { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi}\overline X )_k - X { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi})_k \Big]
   .
   \label{E:dnhprime}

Discrete formulation: Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We discretize the equations exactly in the form shown in , , and ,
obtaining

.. math::
   :label: E:PEmomD

   {{P}}^{-1}
   {{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial t}}} & = 
   -
   ( {{\mathbf{\zeta}}}+ f ) {\hat{k}}{{\times}}{{{\smash[t]{\vec{u}}}}}+ {\nabla_{\rm h}}( \frac12 {{{\smash[t]{\vec{u}}}}}^2 + \Phi )
     - { \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}({{{\smash[t]{\vec{u}}}}})  - \frac{RT_v}{p} {\nabla_{\rm h}}( p )  \\
   {{P}}^{-1}
   {{\frac{\partial {T}}{\partial t}}} & = 
   - {{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}( T  )  - { \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}(T)  + \frac{RT_v}{c^*_p p} \omega  \\
   {{P}}^{-1}
   {{\frac{\partial {p_s}}{\partial t}}} & = 
   -  \sum_{j = 1}^{K} {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_{j} \Delta\eta_{j}   \\
   ( {{\dot\eta}}{\pi})_{i+1/2} & =
   B(\eta_{i+1/2}) 
    \sum_{j = 1}^{K} {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_{j} \Delta\eta_{j}
   - \sum_{j=1}^i {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_{j} \Delta\eta_{j}.

We consider :math:`({{\dot\eta}}{\pi})` a single quantity given at layer
interfaces and defined by . The no-flux boundary condition is
:math:`({{\dot\eta}}{\pi})_{1/2} = ({{\dot\eta}}{\pi})_{K+1/2} = 0`. In
, we used a midpoint quadrature rule to evaluate the indefinite integral
from . In practice :math:`\Delta\eta` can be eliminated from the
discrete equations by scaling :math:`{\pi}`, but here we retain them so
as to have a direct correspondence with the continuum form of the
equations written in terms of
:math:`{\frac{\partial {p}}{\partial \eta}}`.

Finally we give the approximations for the diagnostic equations. We
first integrate to layer interface :math:`i-\frac12` using the same
mid-point rule as used to derive , and then add an additional term
representing the integral from :math:`i-\tfrac12` to :math:`i`:

.. math::
   :label: E:omegaD1

   \omega_i & =   ({{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}p)_i  - 
   \sum_{j=1}^{i-1} {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_j \Delta\eta_j
   +  {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_i \frac{\Delta\eta_i}{2}   \\
   & =   ({{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}p)_i  - \sum_{j=1}^{K} C_{ij} {\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}})_j

where

.. math::

   C_{ij} = 
   \begin{cases}
   \Delta\eta_j &  \quad i>j\\
   {\Delta\eta_j}/{2} &  \quad i = j \\
   0&  \quad i<j\\
   \end{cases}

and similar for :math:`\Phi`,

.. math::

   (\Phi - \Phi_s)_i & =   
   ( \frac{R T_v}{p} {\pi})_i  \frac{\Delta\eta_i}{2}
   + \sum_{j=i+1}^{K} ( \frac{R T_v}{p} {\pi})_j \Delta\eta_j
   \\
    & =    \sum_{j=1}^{K} H_{ij} ( \frac{R T_v}{p} {\pi})_j
   \label{E:Phi}

where

.. math::

   H_{ij} = 
   \begin{cases}
   \Delta\eta_j &  \quad i<j\\
   {\Delta\eta_j}/{2} &  \quad i = j \\
   0&  \quad i>j\\
   \end{cases}

Similar to [eul:econserve], we note that

.. math::

   \Delta\eta_i \, C_{ij}  = \Delta\eta_j \, H_{ji} 
   \label{E:quadconsistency}

which ensures energy conservation .

Consistency
~~~~~~~~~~~

It is important that the discrete equations be as consistent as
possible. In particular, we need a discrete version of , the
non-vertically averaged continuity equation. Equation implicitly implies
such an equation. To see this, apply :math:`{ \mathop{\delta_\eta}}` to
and using that
:math:`\partial p / \partial t = B(\eta) \partial p_s / \partial t`
then we can derive, at layer mid-points,

.. math::

   {{P}}^{-1} {{\frac{\partial {{\pi}}}{\partial t}}} = -{\nabla_{\rm h} \cdot }( {\pi}{{{\smash[t]{\vec{u}}}}}) - 
   { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi}).
   \label{E:PEcontD}

A second type of consistency that has been identified as important is
that , the discrete equation for :math:`\omega`, be consistent with ,
the discrete continuity equation . The two discrete equations should
imply a reasonable discretization of :math:`\omega = Dp/Dt`. To show
this, we take the average of at layers :math:`i-1/2` and :math:`i+1/2`
and combine this with (at layer mid-points :math:`i`) and assuming that
:math:`B(\eta_i) = B(\eta_{i-1/2}) +  B(\eta_{i+1/2})` we obtain

.. math::

   {{P}}^{-1} {{\frac{\partial {p}}{\partial t}}}  = \omega_i -  ({{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}p)_i
   - \frac12 ( ({{\dot\eta}}{ \mathop{\delta_\eta}})_{i-1/2} + ({{\dot\eta}}{ \mathop{\delta_\eta}})_{i+1/2} ).

which, since :math:`{{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}p` is
given at layer mid-points and :math:`{{\dot\eta}}{\pi}` at layer
interfaces, is the SEM discretization of
:math:`w  = \partial p / \partial t + {{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}p + {{\dot\eta}}{\pi}`.

Time Stepping
~~~~~~~~~~~~~

Applying the SEM discretization to - results in a system of ODEs. These
are solved with an :math:`N`-stage Runge-Kutta method. This method
allows for a gravity-wave based CFL number close to :math:`N-1`,
(normalized so that the largest stable timestep of the Robert filtered
Leapfrog method has a CFL number of 1.0). The value of :math:`N` is
chosen large enough so that the dynamics will be stable at the same
timestep used by the tracer advection scheme. To determine :math:`N`, we
first note that the tracer advection scheme uses a less efficient (in
terms of maximum CFL) strong stability preserving Runge-Kutta method
described below. It is stable at an advective CFL number of 1.4. Let
:math:`u_0` be a maximum wind speed and :math:`c_0` be the maximum
gravity wave speed. The gravity wave and advective CFL conditions are

.. math::

   \Delta t \le (N-1) \Delta x / c_0,
   \qquad
   \Delta t \le 1.4 \Delta x / u_0.

In the case where :math:`\Delta t` is chosen as the largest stable
timestep for advection, then we require :math:`N \ge 1 + 1.4 c_0/u_0`
for a stable dynamics timestep. Using a typical values :math:`u_0=120`
m/s and :math:`c_0 = 340`\ m/s gives :math:`N=5`. CAM places additional
restrictions on the timestep (such as that the physics timestep must be
an integer multiple of :math:`\Delta t`) which also influence the choice
of :math:`\Delta t` and :math:`N`.

Dissipation
~~~~~~~~~~~

A horizontal hyper-viscosity operator, modeled after [eul:hdiff] is
applied to the momentum and temperature equations. It is applied in a
time-split manor after each dynamics timestep. The hyper-viscosity step
for vectors can be written as

.. math:: {{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial t}}} = -\nu \Delta^2 {{{\smash[t]{\vec{u}}}}}.

An integral form of this equation suitable for the SEM is obtained using
a mixed finite element formulation (following ) which writes the
equation as a system of equations involving only first derivatives. We
start by introduced an auxiliary vector :math:`{{\smash[t]{\vec{f}}}}`
and using the identity

:math:
   \Delta {{{\smash[t]{\vec{u}}}}}= {\nabla}( {\nabla\cdot}{{{\smash[t]{\vec{u}}}}}) - {{{\nabla}\times}}( {{{\nabla}\times}}{{{\smash[t]{\vec{u}}}}}),

.. math::
   :label: E:HV

   {{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial t}}} & = -\nu ( {\nabla}( {\nabla\cdot}{{\smash[t]{\vec{f}}}}) - {{{\nabla}\times}}{\hat{k}}({{{\nabla}\times}}{{\smash[t]{\vec{f}}}}) )  \\
   {{\smash[t]{\vec{f}}}}& =   {\nabla}({\nabla\cdot}{{{\smash[t]{\vec{u}}}}}) - {{{\nabla}\times}}({{{\nabla}\times}}{{{\smash[t]{\vec{u}}}}}) {\hat{k}}.

Integrating the gradient and curl operators by parts gives

.. math::
   :label: E:weakHV1

   \iint {{\smash[t]{\vec{\phi}}}}\cdot {{\frac{\partial {{{{\smash[t]{\vec{u}}}}}}}{\partial t}}} \,{d\cal{A}}& = \nu \iint [ 
   ({\nabla\cdot}{{\smash[t]{\vec{\phi}}}}) ( {\nabla\cdot}{{\smash[t]{\vec{f}}}}) 
   + ({{{\nabla}\times}}{{\smash[t]{\vec{\phi}}}}) \cdot   {\hat{k}}( {{{\nabla}\times}}{{\smash[t]{\vec{f}}}}) ]
   \,{d\cal{A}} \\
   \iint {{\smash[t]{\vec{\phi}}}}\cdot {{\smash[t]{\vec{f}}}}\,{d\cal{A}} & =  
   - \iint [ ({\nabla\cdot}{{\smash[t]{\vec{\phi}}}}) ({\nabla\cdot}{{{\smash[t]{\vec{u}}}}}) 
   + ({{{\nabla}\times}}{{\smash[t]{\vec{\phi}}}})\cdot {\hat{k}}({{{\nabla}\times}}{{{\smash[t]{\vec{u}}}}}) ] \, {d\cal{A}}.

The SEM Galerkin solution of this integral equation is most naturally
written in terms of an inverse mass matrix instead of the projection
operator. It can be written in terms of the SEM projection operator by
first testing with the product of the element cardinal functions and the
contravariant basis vector
:math:`{{\smash[t]{\vec{\phi}}}}= \phi_{{\smash[t]{\vec{\imath}}}} {{\smash[t]{\vec{g}}}}_\alpha`.
With this type of test function, the RHS of can be defined as a weak
Laplacian operator
:math:`{{\smash[t]{\vec{f}}}}= D({{{\smash[t]{\vec{u}}}}}) \in {{{\cal V^{0}_{\rm cov}}}}`.
The covariant components of :math:`{{\smash[t]{\vec{f}}}}` given by
:math:`f_\alpha = {{\smash[t]{\vec{f}}}}\cdot {{\smash[t]{\vec{g}}}}_\alpha`
are then

.. math::

   f_\alpha ({{\smash[t]{\vec{r}}}}({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}} ; m) ) = \frac{-1}{w_{i^1}w_{i^2}{J}_m({{\smash[t]{\vec{\xi}}}}_{{\smash[t]{\vec{\imath}}}})} 
   {\langle}({{\nabla_{\rm h} \cdot }}\phi_{{\smash[t]{\vec{\imath}}}} {{\smash[t]{\vec{g}}}}_{\alpha}) ({{\nabla_{\rm h} \cdot }}{{{\smash[t]{\vec{u}}}}}) 
   + ({{{\nabla_{\rm h}}}{{\times}}}\phi_{{\smash[t]{\vec{\imath}}}} {{\smash[t]{\vec{g}}}}_\alpha )\cdot {\hat{k}}({{{\nabla_{\rm h}}}{{\times}}}{{{\smash[t]{\vec{u}}}}}).
   {\rangle}

 Then the SEM solution to and is given by

.. math::

   {{{\smash[t]{\vec{u}}}}}(t + \Delta t)   = {{{\smash[t]{\vec{u}}}}}(t) - \nu \Delta t {{P}}\Bigg(  D 
   \Big( {{P}}\big( D({{{\smash[t]{\vec{u}}}}})  \big) 
   \Big ) 
   \Bigg).

Because of the SEM tensor product decomposition, the expression for
:math:`D` can be evaluated in only :math:`O(d)` operations per grid
point, and in CAM-HOMME typically :math:`d=3`.

Following [eul:hdiff], a correction term is added so the hyper-viscosity
does not damp rigid rotation. The hyper-viscosity formulation used for
scalars such as :math:`T` is much simpler, since instead of the vector
Laplacian identity we use :math:`\Delta
T = {\nabla\cdot}{\nabla}T`. Otherwise the approach is identical to that
used above so we omit the details. The correction for terrain following
coordinates given in [eul:hdiff] is not yet implemented in CAM-HOMME.

Discrete formulation: Tracer Advection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All tracers, including specific humidity, are advected with a
discretized version of . HOMME uses the vertically Lagrangian approach
(see [FVvdisc]) from . At the beginning of each timestep, the tracers
are assumed to be given on the :math:`\eta`-coordinate layer mid points.
The tracers are advanced in time on a moving vertical coordinate system
:math:`\eta'` defined so that :math:`{{\dot\eta}}' = 0`. At the end of
the timestep, the tracers are remapped back to the
:math:`\eta`-coordinate layer mid points using the monotone remap
algorithm from .

The horizontal advection step consists of using the SEM to solve

.. math::

   {{\frac{\partial {}}{\partial t}}} ( {\pi}q ) = 
   - {\nabla_{\rm h} \cdot }( {\overline{({\pi}{{{\smash[t]{\vec{u}}}}})}}q  )  
   \label{E:PEqDlagrange}

on the surfaces defined by the :math:`\eta'` layer mid points. The
quantity :math:`{\overline{({\pi}{{{\smash[t]{\vec{u}}}}})}}`
is the mean flux computed during the dynamics update. The mean flux used
in , combined with a suitable mean vertical flux used in the remap stage
allows HOMME to preserve mass/tracer-mass consistency: The tracer
advection of :math:`{\pi}q` with :math:`q=1` will be identical to the
advection of :math:`{\pi}` implied from . The mass/tracer-mass
consistency capability is not in the version of HOMME included in CAM
4.0, but should be in all later versions.

The equation is discretized in time using the optimal 3 stage strong
stability preserving (SSP) second order Runge-Kutta method from . The
RK-SSP method is chosen because it will preserve the monotonicity
properties of the horizontal discretization. RK-SSP methods are convex
combinations of forward-Euler timesteps, so each stage :math:`s` of the
RK-SSP timestep looks like

.. math::

   ( {\pi}q )^{s+1} = 
   ( {\pi}q )^s 
   -  \Delta t {\nabla_{\rm h} \cdot }( {\overline{({\pi}{{{\smash[t]{\vec{u}}}}})}}q^s  )  
   \label{E:PEqDlagrange2}

Simply discretizing this equation with the SEM will result in locally
conservative, high-order accurate but oscillatory transport scheme. A
limiter is added to reduce or eliminate these oscillations . HOMME
supports both monotone and sign-preserving limiters, but the most
effective limiter for HOMME has not yet been determined. The default
configuration in CAM4 is to use the sign-preserving limiter to prevent
negative values of :math:`q` coupled with a sign-preserving
hyper-viscosity operator which dissipates :math:`q^2`.

Conservation and Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The SEM is compatible, meaning it has a discrete version of the
divergence theorem, Stokes theorem and curl/gradient annihilator
properties . The divergence theorem is the key property of the
horizontal discretization that is needed to show conservation. For an
arbitrary scalar :math:`h` and vector :math:`{{{\smash[t]{\vec{u}}}}}`
at layer mid-points, the divergence theorem (or the divergence/gradient
adjoint relation) can be written

.. math:: \int h {\nabla\cdot}{{{\smash[t]{\vec{u}}}}}\,{d\cal{A}}+ \int {{{\smash[t]{\vec{u}}}}}{\nabla}h \,{d\cal{A}}= 0.

The discrete version obeyed by the SEM discretization, using , is given
by

.. math::

   {\langle}h {\nabla_{\rm h} \cdot }{{{\smash[t]{\vec{u}}}}}{\rangle}+ {\langle}{{{\smash[t]{\vec{u}}}}}\cdot {\nabla_{\rm h}}h {\rangle}= 0.
   \label{E:IBPDA1}

The discrete divergence and Stokes theorem apply locally at the element
with the addition of an element boundary integral. The local form is
used to show local conservation of mass and that the horizontal
advection operator locally conserves the two-dimensional potential
vorticity .

In the vertical, showed that the :math:`{ \mathop{\delta_\eta}}` and
:math:`{ \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}` operators
needed to satisfy two integral identities to ensure conservation. For
any :math:`{{\dot\eta}}` layer interface velocity which satisfies
:math:`{{\dot\eta}}_{1/2}={{\dot\eta}}_{K+1/2}=0` and :math:`f,g`
arbitrary functions of layer mid points. The first identity is the
adjoint property (compatibility) for :math:`{ \mathop{\delta_\eta}}` and
:math:`{\pi}`,

.. math::
   :label: E:IBPDN1

   \sum_{i=1}^K \Delta \eta_i\,   {\pi}_i { \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}(f) + 
   \sum_{i=1}^K \Delta \eta_i\,  f_i { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi}) = 0

which follows directly from the definition of the
:math:`{ \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}` difference
operator given in . The second identity we write in terms of
:math:`{ \mathop{\delta_\eta}}`,

.. math::

   \sum_{i=1}^K \Delta\eta_i \,  
    f g { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi})  = 
   \sum_{i=1}^K \Delta\eta_i \,  
   f  { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi}\overline g) 
   +
   \sum_{i=1}^K \Delta\eta_i \,  
   g  { \mathop{\delta_\eta}}( {{\dot\eta}}{\pi}\overline f)
   \label{E:IBPDN2}

which is a discrete integrated-by-parts analog of :math:`\partial (fg) = f \partial g + g \partial f.` 
Construction of methods with both properties on a staggered unequally
spaced grid is the reason behind the complex definition for
:math:`{ \mathop{\overline{ {{\dot\eta}}\delta_\eta }}}` in .

The energy conservation properties of CAM-HOMME were studied in using
the aqua planet test case . CAM-HOMME uses

.. math::

   E =   
    {\langle}\sum_{i=1}^K \Delta \eta_i {\pi}_i  ( \frac12 {{{\smash[t]{\vec{u}}}}}^2 + c_p^* T  )_i
   {\rangle}+
   {\langle}p_s \Phi_s 
   {\rangle}

as the discretization of the total moist energy . The conservation of
:math:`E` is *semi-discrete*, meaning that the only error in
conservation is the time truncation error. In the adiabatic case (with
no hyper-viscosity and no limiters), running from a fully spun up
initial condition, the error in conservation decreases to machine
precision at a second-order rate with decreasing timestep. In the full
non-adiabatic case with a realistic timestep,
:math:`dE/dt \sim 0.013 \text{W/m}^2`.

The CAM physics conserve a dry energy :math:`E_{\text dry}` from which
is not conserved by the moist primitive equations. Although
:math:`E-E_{\text dry}` is small, adiabatic processes in the primitive
equations result in a net heating
:math:`dE_{\text dry}/dt \sim 0.5 \text{W/m}^2` . If it is desired that
the dynamical core conserve :math:`E_\text{dry}` instead of :math:`E`,
HOMME uses the energy fixer from [energyfixer].

Eulerian Dynamical Core
-----------------------

The hybrid vertical coordinate that has been implemented in is described
in this section. The hybrid coordinate was developed by in order to
provide a general framework for a vertical coordinate which is terrain
following at the Earth’s surface, but reduces to a pressure coordinate
at some point above the surface. The hybrid coordinate is more general
in concept than the modified :math:`\sigma` scheme of , which is used in
the GFDL SKYHI model. However, the hybrid coordinate is normally
specified in such a way that the two coordinates are identical.

The following description uses the same general development as , who
based their development on the generalized vertical coordinate of . A
specific form of the coordinate (the hybrid coordinate) is introduced at
the latest possible point. The description here differs from in allowing
for an upper boundary at finite height (nonzero pressure), as in the
original development by Kasahara. Such an upper boundary may be required
when the equations are solved using vertical finite differences.

Generalized terrain-following vertical coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Deriving the primitive equations in a generalized terrain-following
vertical coordinate requires only that certain basic properties of the
coordinate be specified. If the surface pressure is :math:`\pi`, then we
require the generalized coordinate :math:`\eta(p,\pi)` to satisfy:

#. :math:`\eta(p,\pi)` is a monotonic function of :math:`p`.

#. :math:`\eta(\pi,\pi)=1`

#. :math:`\eta(0,\pi)=0`

#. :math:`\eta(p_t,\pi)=\eta_t` where :math:`p_t` is the top of the
   model.

The latter requirement provides that the top of the model will be a
pressure surface, simplifying the specification of boundary conditions.
In the case that :math:`p_t=0`, the last two requirements are identical
and the system reduces to that described in . The boundary conditions
that are required to close the system are:

.. math::
   :label: 3.a.1-2

   \dot\eta(\pi,\pi) & = 0, \\ 
   \dot\eta(p_t,\pi) & = \omega(p_t) = 0. 

Given the above description of the coordinate, the continuous system of
equations can be written following and . The prognostic equations are:

.. math::
   :label: 3.a3-7

   \frac{\partial\zeta}{\partial t}  & = \mathbf{k}\cdot\nabla\times ({\mathbf{n}}/\cos\phi) +  F_{\zeta_H} \,  \\ 
   \frac{\partial\delta}{\partial t} & = \nabla\cdot ({\mathbf{n}}/\cos\phi) -\nabla^2(E+\Phi ) + F_{\delta_H} \,  \\
   \frac{\partial T}{\partial t}     & = \frac{-1}{a\cos^2\phi} [\frac{\partial}{\partial\lambda} (UT) + \cos\phi \frac{\partial}{\partial\phi} (VT) ] + T\delta - \dot\eta  \frac{\partial T}{\partial\eta} + \frac{R}{c_p^*}{T_v}  \frac{\omega}{p} \nonumber \\ 
                                     &   \phantom{=} + Q + F_{T_H} + F_{F_H} ,  \\ 
   \frac{\partial q}{\partial t}     & =  \frac{-1}{a\cos^2\phi} [ \frac{\partial}{\partial\lambda} (Uq) + \cos\phi \frac{\partial}{\partial\phi } (Vq) ] + q\delta - \dot\eta \frac{\partial q}{\partial\eta} + S  \\
   \frac{\partial \pi}{\partial t}   & = \int_1^{\eta_t} {\mathbf{\nabla}\cdot} ( \frac{\partial p}{\partial\eta} {\mathbf{V}} ) d\eta . 

The notation follows standard conventions, and the following terms have been introduced with :math:`{\mathbf{n}} = (n_U,n_V)`:

.. math::
   :label: 3.a.8-13

   n_U    & = + (\zeta+f)V - \dot\eta \frac{\partial U}{\partial\eta} R\frac{T_v}{p}\frac{1}{a} - \frac{\partial p}{\partial\lambda} + F_U  \\ 
   n_V    & = - (\zeta+f)U - \dot\eta \frac{\partial V}{\partial\eta} - R\frac{T_v}{p}\frac{\cos\phi}{a} \frac{\partial p}{\partial\phi} + F_V  \\ 
   E      & = \frac{U^2+V^2}{2\cos^2\phi}   \\ 
   (U,V)  & = (u,v)\cos\phi   \\ 
   T_v    & = [ 1 + ( \frac{R_v}{R} -1 ) q ] T   \\ 
   c_p^\* & = [ 1 + ( \frac{c_{p_v}}{c_p} -1 ) q ] c_p  . 

The terms :math:`F_U, F_V, Q`, and :math:`S` represent the sources and
sinks from the parameterizations for momentum (in terms of :math:`U` and
:math:`V`), temperature, and moisture, respectively. The terms
:math:`F_{\zeta_H}` and :math:`F_{\delta_H}` represent sources due to
horizontal diffusion of momentum, while :math:`F_{T_H}` and
:math:`F_{F_H}` represent sources attributable to horizontal diffusion
of temperature and a contribution from frictional heating (see sections
on horizontal diffusion and horizontal diffusion correction).

In addition to the prognostic equations, three diagnostic equations are
required:

.. math::
   :label: 3.a.14-16

   \Phi                                     & = \Phi_s + R\int_{p(\eta)}^{p(1)}{T_v} d\ln p  \\ 
   \dot\eta \frac{\partial p}{\partial\eta} & = -\frac{\partial p}{\partial t} - \int^\eta_{\eta_t} {\mathbf{\nabla}\cdot}(\frac{\partial p}{\partial\eta}{\mathbf{V}}) d\eta  \\ 
   \omega                                   & = {\mathbf{V} \cdot\nabla}p - \int^\eta_{\eta_t} {\mathbf{\nabla}\cdot} ( \frac{\partial p}{\partial\eta} {\mathbf{V}} ) d\eta 

Note that the bounds on the vertical integrals are specified as values
of :math:`\eta` (:math:`\eta_t`, 1) or as functions of :math:`p` (:math:`p` (1), which is the pressure at :math:`\eta = 1`).

Conversion to final form
~~~~~~~~~~~~~~~~~~~~~~~~

Equations ([3.a.1])-([3.a.16]) are the complete set which must be solved
by a GCM. However, in order to solve them, the function
:math:`\eta(p,\pi)` must be specified. In advance of actually specifying
:math:`\eta(p,\pi)`, the equations will be cast in a more convenient
form. Most of the changes to the equations involve simple applications
of the chain rule for derivatives, in order to obtain terms that will be
easy to evaluate using the predicted variables in the model. For
example, terms involving horizontal derivatives of :math:`p` must be
converted to terms involving only :math:`\partial p/\partial\pi` and
horizontal derivatives of :math:`\pi`. The former can be evaluated once
the function :math:`\eta(p,\pi)` is specified.

The vertical advection terms in ([3.a.5]), ([3.a.6]), ([3.a.8]), and
([3.a.9]) may be rewritten as:

.. math::

   \dot\eta \frac{\partial \psi}{\partial\eta} = \dot\eta \frac{\partial
    p}{\partial\eta}
   \frac{\partial\psi}{\partial p}\, , \label{3.a.17}

since :math:`\dot\eta {\partial p/\partial\eta}` is given by ([3.a.15]).
Similarly, the first term on the right-hand side of ([3.a.15]) can be
expanded as

.. math::

   \frac{\partial p}{\partial t} = \frac{\partial p}{\partial\pi}
    \frac{\partial\pi}{\partial t} ,\label{3.a.18}

and ([3.a.7]) invoked to specify :math:`\partial\pi/\partial t`.

The integrals which appear in ([3.a.7]), ([3.a.15]), and ([3.a.16]) can
be written more conveniently by expanding the kernel as

.. math::

   {\mathbf{\nabla}\cdot} ( \frac{\partial p}{\partial\eta}
    {\mathbf{V}} ) = {\mathbf{V}\cdot\nabla} (\frac{\partial
    p}{\partial\eta}) + \frac{\partial p}{\partial\eta}
    {\mathbf{\nabla}\cdot\mathbf{V}} \ .\label{3.a.19}

The second term in ([3.a.19]) is easily treated in vertical integrals,
since it reduces to an integral in pressure. The first term is expanded
to:

.. math::

   {\mathbf{V}\cdot\nabla} (\frac{\partial p}{\partial\eta})
   & = {\mathbf{V}\cdot}\frac{\partial}{\partial\eta}(\nabla
   p) \nonumber \\
   & = {\mathbf{V}\cdot}\frac{\partial}{\partial\eta}
         (\frac{\partial p}{\partial\pi}\nabla\pi) \nonumber
         \\
   & = {\mathbf{V}\cdot}\frac{\partial}{\partial\eta}
         (\frac{\partial p}{\partial\pi}) \nabla\pi
    + {\mathbf{V}\cdot}\frac{\partial p}{\partial\pi}
        \nabla(\frac{\partial\pi}{\partial\eta})\,
        . \label{3.a.20}

The second term in ([3.a.20]) vanishes because
:math:`\partial\pi/\partial\eta=0`, while the first term is easily
treated once :math:`\eta(p,\pi)` is specified. Substituting ([3.a.20])
into ([3.a.19]), one obtains:

.. math::
   :label: 3.a.21

   {\mathbf{\nabla}\cdot} ( \frac{\partial p}{\partial\eta}
    {\mathbf{V}} )
    = \frac{\partial}{\partial\eta} (\frac{\partial
         p}{\partial\pi}) {\mathbf{V}\cdot\nabla}\pi
    + \frac{\partial p}{\partial\eta} {\mathbf{\nabla}\cdot V} 

Using ([3.a.21]) as the kernel of the integral in ([3.a.7]), ([3.a.15]),
and ([3.a.16]), one obtains integrals of the form

.. math::
   :label: 3.a.22

   \int {\mathbf{\nabla}\cdot} ( \frac{\partial p}{\partial\eta}  {\mathbf{V}}) d\eta  
    & = \int [ \frac{\partial}{\partial\eta} (\frac{\partial p}{\partial\pi}) {\mathbf{V}\cdot\nabla}\pi + \frac{\partial p}{\partial\eta} {\mathbf{\nabla}\cdot V}] d\eta \nonumber \\
    & = \int {\mathbf{V}\cdot\nabla}\pi d(\frac{\partial p}{\partial\pi}) + \int \delta dp.

The original primitive equations ([3.a.3])-([3.a.7]), together with
([3.a.8]), ([3.a.9]), and ([3.a.14])-([3.a.16]) can now be rewritten
with the aid of ([3.a.17]), ([3.a.18]), and ([3.a.22]).

.. math::
   :label: 3.a.23-3.a.30

    \frac{\partial\zeta}{\partial t} & = {\mathbf{k}}\cdot\nabla\times
   ({\mathbf{n}}/\cos\phi) + F_{\zeta_H}  \\
   \frac{\partial\delta}{\partial t} & = {\mathbf{\nabla}\cdot
   (\mathbf{n}/\cos\phi)}
    -\nabla^2(E+\Phi ) + F_{\delta_H}  \\
   \frac{\partial T}{\partial t} & = \frac{-1}{a\cos^2\phi}
    [ \frac{\partial}{\partial\lambda} (UT) + \cos\phi
      \frac{\partial}{\partial\phi} (VT) ] + T\delta
     - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
     T}{\partial p} + \frac{R}{c_p^*}{T_v}\frac{\omega}{p} \nonumber \\
     & \phantom{=} + Q + F_{T_H} + F_{F_H}  \\ 
     \frac{\partial
     q}{\partial t}
   & = \frac{-1}{a\cos^2\phi}
        [ \frac{\partial}{\partial\lambda} (Uq) + \cos\phi
            \frac{\partial}{\partial\phi} (Vq) ]
      + q\delta - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      q}{\partial p} + S  \\
   \frac{\partial \pi}{\partial t}
   & =-\int_{(\eta_t)}^{(1)} {\mathbf{V}\cdot\nabla}\pi
           d(\frac{\partial p}{\partial\pi})
     -\int_{p(\eta_t)}^{p(1)} \delta dp  \\ 
     n_U & = + (\zeta+f)V 
      - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      - U}{\partial p}
      - R\frac{T_v}{a}\frac{1}{p} \frac{\partial p}{\partial\pi}
         \frac{\partial\pi}{\partial\lambda}
      + F_U  \\ 
	n_V & = - (\zeta+f)U
      - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      - V}{\partial p} R\frac{{T_v}\cos \phi}{a} \frac{1}{p}
       \frac{\partial p}{\partial\pi} \frac{\partial\pi}{\partial\phi} +
      F_V  \\
   \Phi & = \Phi_s + R\int_{p(\eta)}^{p(1)}{T_v} d\ln p  \\ 
   \dot\eta \frac{\partial p}{\partial\eta}
   & = \frac{\partial p}{\partial\pi}
     [ \int_{(\eta_t)}^{(1)} {\mathbf{V}}\cdot\nabla\pi
      d(\frac{\partial p}{\partial\pi})
     +\int_{p(\eta_t)}^{p(1)} \delta dp ] \\ 
     \nonumber
   & \phantom{=} -\int_{(\eta_t)}^{(\eta)} {\mathbf{V}}\cdot\nabla\pi
      d( \frac{\partial p}{\partial\pi})
      -\int_{p(\eta_t)}^{p(\eta)} \delta dp  \\
   \omega & = \frac{\partial p}{\partial\pi} {\mathbf{V} \cdot\nabla}\pi
      - \int_{(\eta_t)}^{(\eta)} {\mathbf{V}\cdot\nabla}\pi
             d(\frac{\partial p}{\partial\pi})
      - \int_{p(\eta_t)}^{p(\eta)} \delta dp 

Once :math:`\eta(p,\pi)` is specified, then
:math:`\partial p/\partial\pi` can be determined and
([3.a.23])-([3.a.32]) can be solved in a GCM.

In the actual definition of the hybrid coordinate, it is not necessary
to specify :math:`\eta(p,\pi)` explicitly, since ([3.a.23])-([3.a.32])
only requires that :math:`p` and :math:`\partial
p/\partial\pi` be determined. It is sufficient to specify
:math:`p(\eta,\pi)` and to let :math:`\eta` be defined implicitly. This
will be done in section [ssec:finitediffeqs]. In the case that
:math:`p(\eta,\pi)=\sigma\pi` and :math:`\eta_t=0`,
([3.a.23])-([3.a.32]) can be reduced to the set of equations solved by
CCM1.

Continuous equations using :math:`\partial\ln(\pi)/\partial t`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In practice, the solutions generated by solving the above equations are
excessively noisy. This problem appears to arise from aliasing problems
in the hydrostatic equation ([3.a.30]). The :math:`\ln p` integral
introduces a high order nonlinearity which enters directly into the
divergence equation ([3.a.24]). Large gravity waves are generated in the
vicinity of steep orography, such as in the Pacific Ocean west of the
Andes.

The noise problem is solved by converting the equations given above,
which use :math:`\pi` as a prognostic variable, to equations using
:math:`\Pi=\ln(\pi)`. This results in the hydrostatic equation becoming
only quadratically nonlinear except for moisture contributions to
virtual temperature. Since the spectral transform method will be used to
solve the equations, gradients will be obtained during the transform
from wave to grid space. Outside of the prognostic equation for
:math:`\Pi`, all terms involving :math:`\nabla\pi` will then appear as
:math:`\pi\nabla\Pi`.

Equations ([3.a.23])-([3.a.32]) become:

.. math::
   :label: 3.a.33-3.a.41

   \frac{\partial\zeta}{\partial t} & = {\mathbf{k}\cdot\nabla\times
   (\mathbf{n}/\cos\phi)} + F_{\zeta_H} \\
   \frac{\partial\delta}{\partial t} & = {\mathbf{\nabla}\cdot
   (\mathbf{n}/\cos\phi)}
      -\nabla^2(E+\Phi ) + F_{\delta_H}  \\
   \frac{\partial T}{\partial t}
   & = \frac{-1}{a\cos^2\phi} [ \frac{\partial}{\partial\lambda}
      (UT) + \cos\phi\frac{\partial}{\partial\phi} (VT) ] + T\delta
      - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      T}{\partial p} + \frac{R}{c_p^*}{T_v}\frac{\omega}{p} \\ 
      \nonumber
   & \phantom{=} + Q + F_{T_H} + F_{F_H}  \\
   \frac{\partial q}{\partial t}
   & = \frac{-1}{a\cos^2\phi} [ \frac{\partial}{\partial\lambda}
      (Uq) + \cos\phi \frac{\partial}{\partial\phi} (Vq) ]
     + q\delta - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
    q}{\partial p}  + S  \\ 
     \frac{\partial \Pi}{\partial t}
   & =-\int_{(\eta_t)}^{(1)} {\mathbf{V}\cdot\nabla}\Pi
      d(\frac{\partial p}{\partial\pi})
     -\frac{1}{\pi}\int_{p(\eta_t)}^{p(1)} \delta dp  \\
   n_U & = + (\zeta+f)V
      - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      - U}{\partial p} R \frac{T_v}{a} \frac{\pi}{p}
      \frac{\partial p}{\partial\pi} \frac{\partial\Pi}{\partial\lambda}
      + F_U  \\
   n_V & = - (\zeta+f)U
      - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
      - V}{\partial p} R\frac{{T_v}\cos\phi}{a} \frac{\pi}{p}
       \frac{\partial p}{\partial\pi} \frac{\partial\Pi}{\partial\phi} +
      F_V  \\
   \Phi & = \Phi_s + R\int_{p(\eta)}^{p(1)}{T_v} d\ln p  \\ 
   \dot\eta \frac{\partial p}{\partial\eta}
   & = \frac{\partial p}{\partial\pi} [ \int_{(\eta_t)}^{(1)} \pi
     {\mathbf{V}}\cdot\nabla\Pi d(\frac{\partial
     p}{\partial\pi})
    +\int_{p(\eta_t)}^{p(1)} \delta dp ] \\ 
    \nonumber  & \phantom{=} -\int_{(\eta_t)}^{(\eta)} \pi
      {\mathbf{V}}\cdot\nabla\Pi d(\frac{\partial
      p}{\partial\pi}) -\int_{p(\eta_t)}^{p(\eta)} \delta dp  \\
   \omega & = \frac{\partial p}{\partial\pi} \pi {\mathbf{V}
   \cdot\nabla}\Pi
      - \int_{(\eta_t)}^{(\eta)} \pi {\mathbf{V}\cdot\nabla}\Pi
       d(\frac{\partial p}{\partial\pi})
      - \int_{p(\eta_t)}^{p(\eta)} \delta dp 

The above equations reduce to the standard :math:`\sigma` equations used
in CCM1 if :math:`\eta=\sigma` and :math:`\eta_t=0`. (Note that in this
case :math:`\partial p / \partial\pi = p/\pi = \sigma`.)

Semi-implicit formulation
~~~~~~~~~~~~~~~~~~~~~~~~~

The model described by ([3.a.33])-([3.a.42]), without the horizontal
diffusion terms, together with boundary conditions ([3.a.1]) and
([3.a.2]), is integrated in time using the semi-implicit leapfrog scheme
described below. The semi-implicit form of the time differencing will be
applied to ([3.a.34]) and ([3.a.35]) without the horizontal diffusion
sources, and to ([3.a.37]). In order to derive the semi-implicit form,
one must linearize these equations about a reference state. Isolating
the terms that will have their linear parts treated implicitly, the
prognostic equations ([3.a.33]), ([3.a.34]), and ([3.a.37]) may be
rewritten as:

.. math::
   :label: 3.a.43-45

   \frac{\partial\delta}{\partial t} & = - {R{T_v}} \nabla^2 \ln p  -\nabla^2\Phi + X_1  \\ 
   \frac{\partial T}{\partial t}     & = + \frac{R}{c_p^*}{T_v}\frac{\omega}{p} - \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial T}{\partial  p} + Y_1  \\
   \frac{\partial\Pi}{\partial t}    & = - \frac{1}{\pi}\int_{p(\eta_t)}^{p(1)} \delta dp + Z_1 

where :math:`X_1, Y_1, Z_1` are the remaining nonlinear terms not
explicitly written in ([3.a.43])-([3.a.45]). The terms involving
:math:`\Phi` and :math:`\omega` may be expanded into vertical integrals
using ([3.a.40]) and ([3.a.42]), while the :math:`\nabla^2 \ln p` term
can be converted to :math:`\nabla^2\Pi`, giving:

.. math::
   :label: 3.a.46-48

   \frac{\partial\delta}{\partial t} & = -{RT} \frac{\pi}{p}\frac{\partial p}{\partial \pi} \nabla^2 \Pi  -R\nabla^2\int_{p(\eta)}^{p(1)}T d\ln p\ + X_2 \\
   \frac{\partial T}{\partial t}     & = - \frac{R}{c_p}\frac{T}{p} \int_{p(\eta_t)}^{p(\eta)}\delta dp - [ \frac{\partial p}{\partial\pi} \int_{p(\eta_t)}^{p(1)} \delta dp  - \int_{p(\eta_t)}^{p(\eta)} \delta dp ] \frac{\partial T}{\partial p}  + Y_2  \\ 
   \frac{\partial\Pi}{\partial t}    & = - \frac{1}{pi}\int_{p(\eta_t)}^{p(1)} \delta dp + Z_2 

Once again, only terms that will be linearized have been explicitly
represented in ([3.a.46])-([3.a.48]), and the remaining terms are
included in :math:`X_2`, :math:`Y_2`, and :math:`Z_2`. Anticipating the
linearization, :math:`T_v` and :math:`c_p^*` have been replaced by
:math:`T` and :math:`c_p` in ([3.a.46]) and ([3.a.47]). Furthermore, the
virtual temperature corrections are included with the other nonlinear
terms.

In order to linearize ([3.a.46])-([3.a.48]), one specifies a reference
state for temperature and pressure, then expands the equations about the
reference state:

.. math::
   :label: 3.a.49-51

   T   & = T^r + T^\prime  \\ 
   \pi & = \pi^r + \pi^\prime  \\ 
   p   & = p^r(\eta,\pi^r) + p^\prime 

In the special case that :math:`p(\eta,\pi)=\sigma\pi`,
([3.a.46])-([3.a.48]) can be converted into equations involving only
:math:`\Pi=\ln\pi` instead of :math:`p`, and ([3.a.50]) and ([3.a.51])
are not required. This is a major difference between the hybrid
coordinate scheme being developed here and the :math:`\sigma` coordinate
scheme in CCM1.

Expanding ([3.a.46])-([3.a.48]) about the reference state
([3.a.49])-([3.a.51]) and retaining only the linear terms explicitly,
one obtains:

.. math::
   label: 3.a.52-54

   \frac{\partial\delta}{\partial t}
   & = -R\nabla^2 [ T^{r} \frac{\pi^r}{p^r} (\frac{\partial
     p}{\partial\pi} )^r \Pi
    + \int_{p^r(\eta)}^{p^r(1)}T^\prime d\ln p^r
    + \int_{p^\prime(\eta)}^{p^\prime(1)}\frac{T^r}{p^r} dp^\prime
      ]
    + X_3  \\ 
      \frac{\partial T}{\partial t}
   & = - \frac{R}{c_p}\frac{T^r}{p^r}
       \int_{p^r(\eta_t)}^{p^r(\eta)}\delta dp^r
     - [ (\frac{\partial p}{\partial\pi})^r
       \int_{p^r(\eta_t)}^{p^r(1)} \delta dp^r
      - \int_{p^r(\eta_t)}^{p^r(\eta)} \delta dp^r ] \frac{\partial
       T^r}{\partial p^r}
     + Y_3 \\ 
       \frac{\partial \Pi}{\partial t} & = -
   \frac{1}{\pi^r}\int_{p^r(\eta_t)}^{p^r(1)} \delta dp^r + Z_3

The semi-implicit time differencing scheme treats the linear terms in
([3.a.52])-([3.a.54]) by averaging in time. The last integral in
([3.a.52]) is reduced to purely linear form by the relation

.. math::
   :label: 3.a.55

   dp^\prime = \pi^\prime d (\frac{\partial p}{\partial\pi})^r  + x \, 

In the hybrid coordinate described below, :math:`p` is a linear
function of :math:`\pi`, so :math:`x` above is zero.

We will assume that centered differences are to be used for the
nonlinear terms, and the linear terms are to be treated implicitly by
averaging the previous and next time steps. Finite differences are used
in the vertical, and are described in the following sections. At this
stage only some very general properties of the finite difference
representation must be specified. A layering structure is assumed in
which field values are predicted on :math:`K` layer midpoints denoted by
an integer index, :math:`\eta_k` (see Figure [figure:1]). The interface
between :math:`\eta_k` and :math:`\eta_{k+1}` is denoted by a
half-integer index, :math:`\eta_{k+1/2}`. The model top is at
:math:`\eta_{1/2}=\eta_t`, and the Earth’s surface is at
:math:`\eta_{K+1/2}=1`. It is further assumed that vertical integrals
may be written as a matrix (of order :math:`K`) times a column vector
representing the values of a field at the :math:`\eta_k` grid points in
the vertical. The column vectors representing a vertical column of grid
points will be denoted by underbars, the matrices will be denoted by
bold-faced capital letters, and superscript :math:`T` will denote the
vector transpose.

.. todo:: put in figure3-1 Vertical level structure of

The finite difference forms of ([3.a.52])-([3.a.54]) may then be written
down as:

.. math::
   label: 3.a.56-58

   \underline{\delta}^{n+1}
   & = \underline{\delta}^{n-1} + 2\Delta t \underline{X}^n \nonumber \\
   & \phantom{=} - 2\Delta t R \underline{b}^r \nabla^2 (
         \frac{\Pi^{n-1} + \Pi^{n+1}}{2} - \Pi^{n} ) \nonumber \\
   & \phantom{=} - 2\Delta t R{\mathbf{H}}^r \nabla^2 (
         \frac{(\underline{T}^\prime)^{n-1} +
         (\underline{T}^\prime)^{n+1}}{2}
              - (\underline{T}^\prime)^{n} ) \nonumber \\
   & \phantom{=} - 2\Delta t R \underline{h}^r \nabla^2
         ( \frac{\Pi^{n-1} + \Pi^{n+1}}{2} - \Pi^{n} )  \\ 
	 \underline{T}^{n+1}
   & =  \underline{T}^{n-1} + 2 \Delta t \underline{Y}^n
    - 2\Delta t {\mathbf{D}}^r ( \frac{\underline{\delta}^{n-1} +
          \underline{\delta}^{n+1}}{2}
               - \underline{\delta}^n )  \\ 
   \Pi^{n+1} & =  \Pi^{n-1} + 2\Delta t Z^n
    - 2\Delta t ( \frac{\underline{\delta}^{n-1} +
        \underline{\delta}^{n+1}}{2}
             - \underline{\delta}^n )^T \frac{1}{\Pi^r}
      \underline{\Delta p}^r 

where :math:`()^n` denotes a time varying value at time step :math:`n`.
The quantities :math:`\underline{X}^n, \underline{Y}^n,` and :math:`Z^n`
are defined so as to complete the right-hand sides of
([3.a.43])-([3.a.45]). The components of :math:`\underline{\Delta
p}^r` are given by
:math:`\Delta p^r_k = p^r_{k + \frac{1}{2}} - p^r_{k -
\frac{1}{2}}`. This definition of the vertical difference operator
:math:`\Delta` will be used in subsequent equations. The reference
matrices :math:`{\mathbf{H}}^r` and :math:`{\mathbf{D}}^r`,
and the reference column vectors :math:`\displaystyle\underline{b}^r`
and :math:`\displaystyle\underline{h}^r`, depend on the precise
specification of the vertical coordinate and will be defined later.

Energy conservation
~~~~~~~~~~~~~~~~~~~

We shall impose a requirement on the vertical finite differences of the
model that they conserve the global integral of total energy *in the
absence of sources and sinks*. We need to derive equations for kinetic
and internal energy in order to impose this constraint. The momentum
equations (more painfully, the vorticity and divergence equations)
without the :math:`F_U, F_V, F_{\zeta_H}` and :math:`F_{\delta_H}`
contributions, can be combined with the continuity equation

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta}
         )
    + \nabla\cdot ( \frac{\partial p}{\partial\eta} {\mathbf{V}}
         )
    + \frac{\partial}{\partial \eta} ( \frac{\partial
         p}{\partial\eta} \dot\eta )
    = 0 \label{3.a.59}

to give an equation for the rate of change of kinetic energy:

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta} E
         )
   & = -\nabla\cdot ( \frac{\partial p}{\partial\eta} E
         {\mathbf{V}} )
    - \frac{\partial}{\partial \eta} ( \frac{\partial
         p}{\partial\eta} E \dot\eta ) \nonumber \\
   & \phantom{=}- \frac{R{T_v}}{p} \frac{\partial p}{\partial\eta}
   {\mathbf{V}}\cdot\nabla p
    - \frac{\partial p}{\partial\eta} {\mathbf{V}}\cdot\nabla\Phi \,
    - . \label{3.a.60}

The first two terms on the right-hand side of ([3.a.60]) are transport
terms. The horizontal integral of the first (horizontal) transport term
should be zero, and it is relatively straightforward to construct
horizontal finite difference schemes that ensure this. For spectral
models, the integral of the horizontal transport term will not vanish in
general, but we shall ignore this problem.

The vertical integral of the second (vertical) transport term on the
right-hand side of ([3.a.60]) should vanish. Since this term is obtained
from the vertical advection terms for momentum, which will be finite
differenced, we can construct a finite difference operator that will
ensure that the vertical integral vanishes.

The vertical advection terms are the product of a vertical velocity
(:math:`\dot\eta \partial p/\partial\eta`) and the vertical derivative
of a field (:math:`\partial\psi/\partial p`). The vertical velocity is
defined in terms of vertical integrals of fields ([3.a.41]), which are
naturally taken to interfaces. The vertical derivatives are also
naturally taken to interfaces, so the product is formed there, and then
adjacent interface values of the products are averaged to give a
midpoint value. It is the definition of the average that must be correct
in order to conserve kinetic energy under vertical advection in
([3.a.60]). The derivation will be omitted here, the resulting vertical
advection terms are of the form:

.. math::
   :label: 3.a.61-2

   ( \dot\eta \frac{\partial p}{\partial\eta} \frac{\partial
         \psi}{\partial p} )_{k}
   & = \frac{1}{2\Delta p_k}
      [ ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k+1/2} ( \psi_{k+1} - \psi_k )
         + ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k-1/2} ( \psi_k - \psi_{k-1} )
      ] , \ \\ 
      \Delta p_k & = p_{k+1/2} - p_{k-1/2}
   .

The choice of definitions for the vertical velocity at interfaces is not
crucial to the energy conservation (although not completely arbitrary),
and we shall defer its definition until later. The vertical advection of
temperature is not required to use ([3.a.61]) in order to conserve mass
or energy. Other constraints can be imposed that result in different
forms for temperature advection, but we will simply use ([3.a.61]) in
the system described below.

The last two terms in ([3.a.60]) contain the conversion between kinetic
and internal (potential) energy and the form drag. Neglecting the
transport terms, under assumption that global integrals will be taken,
noting that :math:`\nabla p/p = \frac{\pi}{p} \frac{\partial
p}{\partial\pi} \nabla \Pi`, and substituting for the geopotential using
([3.a.40]), ([3.a.60]) can be written as:

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta} E
         )
   & =- {R{T_v}} \frac{\partial p}{\partial\eta} {\mathbf{V}} \cdot
      ( \frac{\pi}{p} \frac{\partial p}{\partial\pi} \nabla \Pi
      ) \\
   \nonumber & \phantom{=} - \frac{\partial p}{\partial\eta}
   {\mathbf{V}}\cdot\nabla\Phi_s
    - \frac{\partial p}{\partial\eta} {\mathbf{V}}\cdot\nabla
         \int_{p(\eta)}^{p(1)}R{T_v} d\ln p
    + \, \ldots \label{3.a.63}

The second term on the right-hand side of ([3.a.63]) is a source (form
drag) term that can be neglected as we are only interested in internal
conservation properties. The last term on the right-hand side of
([3.a.63]) can be rewritten as

.. math::

   \frac{\partial p}{\partial\eta} {\mathbf{V}}\cdot\nabla
         \int_{p(\eta)}^{p(1)}R{T_v} d\ln p
   = \nabla\cdot
        \{ \frac{\partial p}{\partial\eta} {\mathbf{V}}
           \int_{p(\eta)}^{p(1)}R{T_v} d\ln p
        \}
   - \nabla\cdot
        ( \frac{\partial p}{\partial\eta} {\mathbf{V}}
        ) \int_{p(\eta)}^{p(1)}R{T_v} d\ln p \, . \label{3.a.64}

The global integral of the first term on the right-hand side of
([3.a.64]) is obviously zero, so that ([3.a.63]) can now be written as:

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta} E
         )
   =- {R{T_v}} \frac{\partial p}{\partial\eta} {\mathbf{V}} \cdot (
      \frac{\pi}{p} \frac{\partial p}{\partial\pi} \nabla \Pi )
   + \nabla\cdot ( \frac{\partial p}{\partial\eta} {\mathbf{V}}
        ) \int_{p(\eta)}^{p(1)}R{T_v} d\ln p
    + ...  \label{3.a.65}

We now turn to the internal energy equation, obtained by combining the
thermodynamic equation ([3.a.35]), without the :math:`Q`,
:math:`F_{T_H}`, and :math:`F_{F_H}` terms, and the continuity equation
([3.a.59]):

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta}
         c_p^* T )
   = -\nabla\cdot ( \frac{\partial p}{\partial\eta} c_p^* T
         {\mathbf{V}} )
    - \frac{\partial}{\partial \eta} ( \frac{\partial
         p}{\partial\eta} c_p^* T \dot\eta )
   + {R{T_v}} \frac{\partial p}{\partial\eta} \frac{\omega}{p} \, .
   \label{3.a.66}

As in ([3.a.60]), the first two terms on the right-hand side are
advection terms that can be neglected under global integrals. Using
([3.a.16]), ([3.a.66]) can be written as:

.. math::

   \frac{\partial}{\partial t} ( \frac{\partial p}{\partial\eta}
         c_p^* T )
   = {R{T_v}} \frac{\partial p}{\partial\eta} {\mathbf{V}} \cdot (
     \frac{\pi}{p} \frac{\partial p}{\partial\pi} \nabla \Pi )
   - {R{T_v}} \frac{\partial p}{\partial\eta} \frac{1}{p}
        \int_{\eta_t}^{\eta} \nabla\cdot
        ( \frac{\partial p}{\partial\eta} {\mathbf{V}}
        ) d\eta + ... \label{3.a.67}

The rate of change of total energy due to internal processes is obtained
by adding ([3.a.65]) and ([3.a.67]) and must vanish. The first terms on
the right-hand side of ([3.a.65]) and ([3.a.67]) obviously cancel in the
continuous form. When the equations are discretized in the vertical, the
terms will still cancel, providing that the same definition is used for
:math:`(1/p\,\,\partial p/\partial\pi)_k` in the nonlinear terms of the
vorticity and divergence equations ([3.a.38]) and ([3.a.39]), and in the
:math:`\omega` term of ([3.a.35]) and ([3.a.42]).

The second terms on the right-hand side of ([3.a.65]) and ([3.a.67])
must also cancel in the global mean. This cancellation is enforced
locally in the horizontal on the column integrals of ([3.a.65]) and
([3.a.67]), so that we require:

.. math::

   \int^1_{\eta_t} \{ \nabla\cdot
        ( \frac{\partial p}{\partial\eta} {\mathbf{V}}
        ) \int_{p(\eta)}^{p(1)}R{T_v} d\ln p
   \} d\eta
   =\int^1_{\eta_t} \{ {R{T_v}} \frac{\partial p}{\partial\eta}
        \frac{1}{p} \int_{\eta_t}^{\eta} \nabla\cdot
        ( \frac{\partial p}{\partial\eta^\prime} {\mathbf{V}}
        ) d\eta^\prime \} d\eta . \label{3.a.68}

The inner integral on the left-hand side of ([3.a.68]) is derived from
the hydrostatic equation ([3.a.40]), which we shall approximate as

.. math::
   :label: 3.a.69-70

   \Phi_k           & = \Phi_s + R\sum_{\ell=k}^K H_{k\ell}{T_v}_\ell , \nonumber \\ 
                    & = \Phi_s + R\sum_{\ell=1}^K H_{k\ell}{T_v}_\ell ,  \\ 
   \underline{\Phi} & = \Phi_s \underline{1} + R {\mathbf{H}} \underline {T_v} , 

where :math:`H_{k\ell}=0` for :math:`\ell<k`. The quantity
:math:`\underline{1}` is defined to be the unit vector. The inner
integral on the right-hand side of ([3.a.68]) is derived from the
vertical velocity equation ([3.a.42]), which we shall approximate as

.. math::

   ( \frac{\omega}{p} )_k
   = ( \frac{\pi}{p} \frac{\partial p}{\partial\pi} )_k
            {\mathbf{V}}_k \cdot \nabla\Pi
   - \sum_{\ell=1}^K C_{k\ell}
         [ \delta_\ell \Delta p_\ell
            + \pi ({\mathbf{V}}_\ell \cdot \nabla \Pi )
              \Delta ( \frac{\partial p}{\partial\pi} )_\ell
         ] , \label{3.a.71}

where :math:`C_{k\ell}=0` for :math:`\ell>k`, and :math:`C_{k\ell}` is
included as an approximation to :math:`1/p_k` for :math:`\ell \leq k`
and the symbol :math:`\Delta` is similarly defined as in ([3.a.62]).
:math:`C_{k\ell}` will be determined so that :math:`\omega` is
consistent with the discrete continuity equation following . Using
([3.a.69]) and ([3.a.71]), the finite difference analog of ([3.a.68]) is

.. math::

   \nonumbereqn{\sum_{k=1}^K
     \{ \frac{1}{\Delta\eta_k}
         [ \delta_k \Delta p_k
            + \pi ({\mathbf{V}}_k \cdot \nabla \Pi ) \Delta
              ( \frac{\partial p}{\partial\pi} )_k
         ] R\sum_{\ell=1}^K H_{k\ell}{T_v}_\ell
     \} \Delta\eta_k } \\
   & & \mbox{} = \sum_{k=1}^K
     \{ R {T_v}_k \frac{\Delta p_k}{\Delta\eta_k}
         \sum_{\ell=1}^K C_{k\ell}
            [ \delta_\ell \Delta p_\ell
               + \pi ( {\mathbf{V}}_\ell \cdot \nabla \Pi )
                 \Delta ( \frac{\partial p}{\partial\pi}
                 )_\ell
            ] \} \Delta\eta_k , \label{3.a.72}

where we have used the relation

.. math::

   \nabla\cdot {\mathbf{V}}(\partial p/\partial\eta )_k =[\delta_k\Delta
   p_k + \\
   \pi( {\mathbf{V}}_k\cdot\nabla\Pi ) \Delta (\partial
       p/\partial\pi)_k ]/\Delta\eta_k

(see [3.a.22]). We can now combine the sums in ([3.a.72]) and simplify
to give

.. math::

   \nonumbereqn{\sum_{k=1}^K \sum_{\ell=1}^K
     \{
         [ \delta_k \Delta p_k
            + \pi ({\mathbf{V}}_k \cdot \nabla \Pi ) \Delta
              ( \frac{\partial p}{\partial\pi} )_k
         ] H_{k\ell}{T_v}_\ell
     \}} \\
   & & \mbox{} = \sum_{k=1}^K \sum_{\ell=1}^K
     \{
         [ \delta_\ell \Delta p_\ell
            + \pi ( {\mathbf{V}}_\ell \cdot \nabla \Pi )
              \Delta ( \frac{\partial p}{\partial\pi} )_\ell
         ] \Delta p_k C_{k\ell}{T_v}_k
     \} . \label{3.a.73}

Interchanging the indexes on the left-hand side of ([3.a.73]) will
obviously result in identical expressions if we require that

.. math:: H_{k\ell} = C_{\ell k} \Delta p_\ell .\label{3.a.74}

Given the definitions of vertical integrals in ([3.a.70]) and ([3.a.71])
and of vertical advection in ([3.a.61]) and ([3.a.62]) the model will
conserve energy as long as we require that :math:`\mathbf{C}` and
:math:`\mathbf{H}` satisfy ([3.a.74]). We are, of course, still
neglecting lack of conservation due to the truncation of the horizontal
spherical harmonic expansions.

Horizontal diffusion
~~~~~~~~~~~~~~~~~~~~

contains a horizontal diffusion term for :math:`T, \zeta`, and
:math:`\delta` to prevent spectral blocking and to provide reasonable
kinetic energy spectra. The horizontal diffusion operator in is also
used to ensure that the CFL condition is not violated in the upper
layers of the model. The horizontal diffusion is a linear
:math:`\nabla^2` form on :math:`\eta` surfaces in the top three levels
of the model and a linear :math:`\nabla^4` form with a partial
correction to pressure surfaces for temperature elsewhere. The
:math:`\nabla^2` diffusion near the model top is used as a simple sponge
to absorb vertically propagating planetary wave energy and also to
control the strength of the stratospheric winter jets. The
:math:`\nabla^2` diffusion coefficient has a vertical variation which
has been tuned to give reasonable Northern and Southern Hemisphere polar
night jets.

In the top three model levels, the :math:`\nabla^2` form of the
horizontal diffusion is given by

.. math::
   :label: 3.a.75-77 

   F_{\zeta_H}   & = K^{(2)} [ \nabla^2 (\zeta + f ) +   2(\zeta + f )/a^2 ] , \\ 
   F_{\delta_H}  & = K^{(2)} [ \nabla^2 \delta + 2(\delta/a^2)] ,  \\
   F_{T_H}       & = K^{(2)} \nabla^2T . 

Since these terms are linear, they are easily calculated in spectral
space. The undifferentiated correction term is added to the vorticity
and divergence diffusion operators to prevent damping of uniform
:math:`(n=1)` rotations . The :math:`\nabla^2` form of the horizontal
diffusion is applied *only* to pressure surfaces in the standard model
configuration.

The horizontal diffusion operator is better applied to pressure surfaces
than to terrain-following surfaces (applying the operator on isentropic
surfaces would be still better). Although the governing system of
equations derived above is designed to reduce to pressure surfaces above
some level, problems can still occur from diffusion along the lower
surfaces. Partial correction to pressure surfaces of harmonic horizontal
diffusion (:math:`\partial\xi/\partial t =
K\nabla^2\xi`) can be included using the relations:

.. math::

   \nabla_p\xi & = \nabla_\eta\xi - p \frac{\partial\xi}{\partial p}
   \nabla_\eta \ln p \nonumber \\ 
   \nabla^2_p\xi & = \nabla^2_\eta\xi
   - p \frac{\partial\xi}{\partial p} \nabla^2_\eta \ln p
   - 2\nabla_\eta( \frac{\partial\xi}{\partial p}
   - )\cdot\nabla_\eta p
   + \frac{\partial^2\xi}{\partial^2 p} \nabla^2_\eta p \,
   . \label{3.a.78}

Retaining only the first two terms above gives a correction to the
:math:`\eta` surface diffusion which involves only a vertical derivative
and the Laplacian of log surface pressure,

.. math::

   \nabla^2_p\xi = \nabla^2_\eta\xi
     - \pi \frac{\partial\xi}{\partial p} \frac{\partial p}{\partial\pi}
         \nabla^2 \Pi + \ldots \label{3.a.79}

Similarly, biharmonic diffusion can be partially corrected to pressure
surfaces as:

.. math::

   \nabla^4_p\xi = \nabla^4_\eta\xi
     - \pi \frac{\partial\xi}{\partial p} \frac{\partial p}{\partial\pi}
         \nabla^4 \Pi + \ldots \label{3.a.80}

The bi-harmonic :math:`\nabla^4` form of the diffusion operator is
applied at all other levels (generally throughout the troposphere) as

.. math::

   F_{\zeta_H} & = -K^{(4)} [\nabla^4 (\zeta + f ) -
   (\zeta + f ) ( 2/a^2)^2 ] , \label{3.a.81} \\ 
   F_{\delta_H} & = -K^{(4)} [\nabla^4 \delta - \delta (2/a^2)^2 ],
    \label{3.a.82} \\
   F_{T_H} & = -K^{(4)} [ \nabla^4 T - \pi \frac{\partial
   T}{\partial p} \frac{\partial p}{\partial \pi} \nabla^4 \Pi ]
   . \label{3.a.83}

The second term in :math:`F_{T_H}` consists of the leading term in the
transformation of the :math:`\nabla^4` operator to pressure surfaces. It
is included to offset partially a spurious diffusion of :math:`T` over
mountains. As with the :math:`\nabla^2` form, the :math:`\nabla^4`
operator can be conveniently calculated in spectral space. The
correction term is then completed after transformation of :math:`T` and
:math:`\nabla^4 \Pi` back to grid–point space. As with the
:math:`\nabla^2` form, an undifferentiated term is added to the
vorticity and divergence diffusion operators to prevent damping of
uniform rotations.

Finite difference equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The governing equations are solved using the spectral method in the
horizontal, so that only the vertical and time differences are presented
here. The dynamics includes horizontal diffusion of :math:`T,
(\zeta + f)`, and :math:`\delta`. Only :math:`T` has the leading term
correction to pressure surfaces. Thus, equations that include the terms
in this time split sub-step are of the form

.. math::

   \frac{\partial \psi}{\partial t} = {\rm Dyn} ( \psi ) -
   (-1)^i K^{(2i)} \nabla^{2i}_\eta \psi \, , \label{3.a.84}

for :math:`(\zeta + f)` and :math:`\delta`, and

.. math::

   \frac{\partial T}{\partial t} = {\rm Dyn} ( T ) - (-1)^i
   K^{(2i)} \{ \nabla^{2i}_\eta T - \pi \frac{\partial T} {\partial
   p} \, \frac{\partial p}{\partial \pi} \nabla^{2i} \Pi \}\, ,
   \label{3.a.85}

where :math:`i = 1` in the top few model levels and :math:`i = 2`
elsewhere (generally within the troposphere). These equations are
further subdivided into time split components:

.. math::

   \psi^{n+1} & = \psi^{n-1} + 2\Delta t \ {\rm Dyn} ( \psi^{n+1},
       \psi^n, \psi^{n-1} ) \ , \label{3.a.86} \\
   \psi^* & = \psi^{n+1} - 2\Delta t \ (-1)^i K^{(2i)} \nabla^{2i}_\eta
       ( {\psi}^{*n+1} ) \ , \label{3.a.87} \\
   \hat \psi^{n+1} & = \psi^* \ , \label{3.a.88}

for :math:`( \zeta + f )` and :math:`\delta`, and

.. math::
   :label:

   T^{n+1}      & = T^{n-1} + 2\Delta t \ {\rm Dyn} ( T^{n+1}, T^n,  T^{n-1} ) \,  \\ 
   T^*          & = T^{n+1} - 2\Delta t \ ( -1 )^i K^{(2i)} \nabla^{2i}\eta ( T^* ) \ ,   \\ 
   \hat T^{n+1} & = T^* + 2\Delta t \ ( -1 )^i K^{(2i)} \pi \, \frac{\partial T^*}{\partial p} \,  \frac{\partial p}{\partial \pi} \, \nabla^{2i} \, \Pi \ , 

for :math:`T`, where in the standard model :math:`i` only takes the
value 2 in ([3.a.91]). The first step from
:math:`( \hskip 5pt )^{n-1}` to
:math:`( \hskip 5pt )^{n+1}` includes the transformation to
spectral coefficients. The second step from :math:`( \hskip 5pt
)^{n+1}` to :math:`( \hat{\hskip 5pt} )^{n+1}` for
:math:`\delta` and :math:`\zeta \,`, or
:math:`( \hskip 5pt )^{n+1}` to :math:`( \hskip
5pt )^*` for :math:`T`, is done on the spectral coefficients, and
the final step from :math:`( \hskip 5pt )^*` to
:math:`( \hat{\hskip
5pt} )^{n+1}` for :math:`T` is done after the inverse transform to
the grid point representation.

The following finite-difference description details only the forecast
given by ([3.a.86]) and ([3.a.89]). The finite-difference form of the
forecast equation for water vapor will be presented later in Section 3c.
The general structure of the complete finite difference equations is
determined by the semi-implicit time differencing and the energy
conservation properties described above. In order to complete the
specification of the finite differencing, we require a definition of the
vertical coordinate. The actual specification of the generalized
vertical coordinate takes advantage of the structure of the equations
([3.a.33])-([3.a.42]). The equations can be finite-differenced in the
vertical and, in time, without having to know the value of :math:`\eta`
anywhere. The quantities that must be known are :math:`p` and
:math:`\partial p/\partial\pi` at the grid points. Therefore the
coordinate is defined implicitly through the relation:

.. math:: p(\eta,\pi) = A(\eta)p_0 + B(\eta)\pi \, , \label{3.a.92}

which gives

.. math:: \frac{\partial p}{\partial\pi} = B(\eta) \, . \label{3.a.93}

A set of levels :math:`\eta_k` may be specified by specifying
:math:`A_k` and :math:`B_k`, such that :math:`\eta_k \equiv A_k + B_k`,
and difference forms of ([3.a.33])-([3.a.42]) may be derived.

The finite difference forms of the Dyn operator ([3.a.33])-([3.a.42]),
including semi-implicit time integration are:

.. math::
   :label:

   \underline{\zeta}^{n+1}
   & =  \underline{\zeta}^{n-1} + 2\Delta t
    {\mathbf{k}\cdot\nabla\times}(\underline{\mathbf{n}}^n/\cos\phi) ,  \\
   \underline{\delta}^{n+1} 
   & =  \underline{\delta}^{n-1} + 2\Delta t [ {\nabla\cdot
          (\underline{\mathbf{n}}^n/\cos\phi)} - \nabla^2 ( \underline{E}^n + \Phi_s \underline{1} +
          R{\mathbf{H}}^n (\underline{T_v}^{'})^n ) ] \nonumber \\
   & \phantom{=} - 2\Delta t R{\mathbf{H}}^r \nabla^2 (
        \frac{(\underline{T}^\prime)^{n-1} +
        (\underline{T}^\prime)^{n+1}}{2}
              - (\underline{T}^\prime)^{n} ) \nonumber \\
   & \phantom{=} - 2\Delta t R( \underline{b}^r + \underline{h}^r
         ) \nabla^2
         ( \frac{\Pi^{n-1} + \Pi^{n+1}}{2} - \Pi^{n}
         ) ,  \\
   (\underline{T}^{'})^{n+1} & =  (\underline{T}^{'})^{n-1}
    - 2 \Delta t [ \frac{1}{a\cos^2\phi}
           \frac{\partial}{\partial\lambda} ( \underline{UT}^\prime
           )^n
         + \frac{1}{a\cos\phi} \frac{\partial}{\partial\phi} (
           \underline{VT}^\prime )^n
         - \underline{\Gamma}^n ] \\ 
	   
    \nonumber & \phantom{=} - 2\Delta t {\mathbf{D}}^r (
          \frac{\underline{\delta}^{n-1} + \underline{\delta}^{n+1}}{2}
               - \underline{\delta}^n ) \,  \\
   \Pi^{n+1} & =  \Pi^{n-1}
    - 2\Delta t \frac{1}{\pi^n} ( (\underline{\delta}^n
         )^T \underline{\Delta p}^n
       + (\underline{\mathbf{V}}^n )^T \cdot \nabla \Pi^n
       \pi^n \underline{\Delta B}
   ) \nonumber \\
   & \phantom{=} - 2\Delta t ( \frac{\underline{\delta}^{n-1} +
          \underline{\delta}^{n+1}}{2}
             - \underline{\delta}^n )^T \frac{1}{\pi^r}
      \underline{\Delta p}^r
   , \\
   ( n_U )_k & =  ( \zeta_k + f ) V_k
    - R{T_v}_k ( \frac{1}{p} \frac{\partial p}{\partial\pi}
      )_k \pi \frac{1}{a} \frac{\partial \Pi}{\partial\lambda}
      \nonumber \\
   & \phantom{=} - \frac{1}{2\Delta p_k}
      [ ( \dot\eta \frac{\partial p}{\partial\eta}
            )_{k+1/2} ( U_{k+1} - U_k )
          + ( \dot\eta \frac{\partial p}{\partial\eta}
            )_{k-1/2} ( U_k - U_{k-1} )
      ] \nonumber \\ & \phantom{=} + ( F_U )_k \ ,
    \\
   ( n_V )_k & =  - ( \zeta_k + f ) U_k
    - R{T_v}_k ( \frac{1}{p} \frac{\partial p}{\partial\pi}
      )_k \pi \frac{\cos \phi}{a} \frac{\partial \Pi}{\partial\phi}
      \nonumber \\
   & \phantom{=} - \frac{1}{2\Delta p_k}
      [ ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k+1/2} ( V_{k+1} - V_k )
         + ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k-1/2} ( V_k - V_{k-1} )
      ]\nonumber \\ 
      & \phantom{=} + ( F_V )_k \ ,  \\
   \Gamma_k  & =  T^{\prime}_k \delta_k + \frac{R{T_v}_k}{(c_p^*)_k} (
    \frac{\omega}{p} )_k - Q \nonumber \\
   & \phantom{=} - \frac{1}{2\Delta p_k}
      [ ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k+1/2} ( T_{k+1} - T_k )
         + ( \dot\eta \frac{\partial p}{\partial\eta}
           )_{k-1/2} ( T_k - T_{k-1} )
      ] , 

.. math::
   :label:

   E_k & = (u_k )^2 + (v_k )^2 ,  \\
   \frac{R {T_v}_{k}}{(c^*_p)_k}
   & = \frac{R}{c_p} ( \frac{T^r_k + {T_v}_k^\prime} {1 + (\frac{c_{p_v}}{c_p}- 1 ) q_k}) ,  \\
       ( \dot\eta \frac{\partial p}{\partial\eta} )_{k+1/2}
   & = B_{k+1/2} \sum^K_{\ell=1}
       [ \delta_\ell \Delta p_\ell + {\mathbf{V}}_\ell \cdot \pi \nabla \Pi \Delta B_\ell ] \nonumber \\
   & \phantom{=} - \sum^k_{\ell=1} [ \delta_\ell \Delta p_\ell + {\mathbf{V}}_\ell \cdot \pi \nabla \Pi \Delta B_\ell] ,  \\
   ( \frac{\omega}{p} )_k
   & = ( \frac{1}{p} \frac{\partial p}{\partial\pi} )_k {\mathbf{V}}_k \cdot \pi \nabla\Pi - \sum^k_{\ell=1} C_{k\ell}
         [ \delta_\ell \Delta p_\ell + {\mathbf{V}}_\ell \cdot \pi \nabla \Pi \Delta B_\ell ]  \\
   C_{k\ell}
   & = \{ \begin{array}{ll} \frac{1}{p_k}, \ell < k \\[6pt]
       \frac{1}{2 p_k}, & \ell = k , \end{array} .   \\
   H_{k\ell}   & = C_{\ell k}\Delta p_\ell, \\
   D^r_{k\ell} & = {\Delta p^r_\ell} \frac{R}{c_p } T^r_k C^r_{\ell k} + \frac{\Delta p^r_\ell}{2\Delta p^r_k}
                   ( T^r_k - T^r_{k-1} ) (\epsilon_{k\ell+1}-B_{k-1/2}) \nonumber \\
               &  \phantom{=} + \frac{\Delta p^r_\ell}{2\Delta p^r_k} ( T^r_{k+1} - T^r_k ) (\epsilon_{k\ell}-B_{k+1/2}) , } \\
   \frac{ \epsilon_{k\ell}}{R}
               & = \{\begin{array}{ll} 1, \ell \leq k \\
                   0, \ell > k, \end{array} .  

where notation such as :math:`( \underline{UT}^\prime )^n`
denotes a column vector with components :math:`( U_k T_k^\prime
)^n`. In order to complete the system, it remains to specify the
reference vector :math:`\underline{h}^r`, together with the term
:math:`(1/p
\, \partial p/\partial\pi)`, which results from the pressure gradient
terms and also appears in the semi-implicit reference vector
:math:`\underline{b}^r`:

.. math::
   :label: 

   ( \frac{1}{p} \frac{\partial p}{\partial\pi} )_k  
                   & = ( \frac{1}{p} )_k ( \frac{\partial p}{\partial\pi} )_k = \frac{B_k}{p_k} ,  \\
   \underline{b}^r & = \underline{T}^r , \\
   \underline{h}^r & = 0 . 

The matrices :math:`{\mathbf{C}}^n` and
:math:`{\mathbf{H}}^n` ( with components :math:`C_{k\ell}` and
:math:`H_{k \ell}`) must be evaluated at each time step and each point
in the horizontal. It is more efficient computationally to substitute
the definitions of these matrices into ([3.a.95]) and ([3.a.104]) at the
cost of some loss of generality in the code. The finite difference
equations have been written in the form ([3.a.94])-([3.a.111]) because
this form is quite general. For example, the equations solved by at
ECMWF can be obtained by changing only the vectors and hydrostatic
matrix defined by ([3.a.108])-([3.a.111]).

Time filter
~~~~~~~~~~~

The time step is completed by applying a recursive time filter
originally designed by and later studied by .

.. math::

   \overline \psi^n = \psi^n + \alpha ( \overline{\psi}^{n-1}
    - 2 \psi^n + \psi^{n+1} )

Spectral transform
~~~~~~~~~~~~~~~~~~

The spectral transform method is used in the horizontal exactly as in
CCM1. As shown earlier, the vertical and temporal aspects of the model
are represented by finite–difference approximations. The horizontal
aspects are treated by the spectral–transform method, which is described
in this section. Thus, at certain points in the integration, the
prognostic variables :math:`(\zeta + f ),
\delta, T,` and :math:`\Pi` are represented in terms of coefficients of
a truncated series of spherical harmonic functions, while at other
points they are given by grid–point values on a corresponding Gaussian
grid. In general, physical parameterizations and nonlinear operations
are carried out in grid–point space. Horizontal derivatives and linear
operations are performed in spectral space. Externally, the model
appears to the user to be a grid–point model, as far as data required
and produced by it. Similarly, since all nonlinear parameterizations are
developed and carried out in grid–point space, the model also appears as
a grid–point model for the incorporation of physical parameterizations,
and the user need not be too concerned with the spectral aspects. For
users interested in diagnosing the balance of terms in the evolution
equations, however, the details are important and care must be taken to
understand which terms have been spectrally truncated and which have
not. The algebra involved in the spectral transformations has been
presented in several publications . In this report, we present only the
details relevant to the model code; for more details and general
philosophy, the reader is referred to these earlier papers.

Spectral algorithm overview
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The horizontal representation of an arbitrary variable :math:`\psi`
  consists of a truncated series of spherical harmonic functions,

  .. math::
     :label:

     \psi(\lambda, \mu) = \sum^M_{m=-M} \; \; \sum^{{\cal N} ( m )}_{n=|m|} \psi^m_n P^m_n (\mu) e^{im \lambda} , 

  where :math:`\mu = \sin \phi, \ M` is the highest Fourier wavenumber
  included in the east–west representation, and :math:`{\cal N}( m
  )` is the highest degree of the associated Legendre polynomials
  for longitudinal wavenumber :math:`m`. The properties of the spherical
  harmonic functions used in the representation can be found in the
  review by . The model is coded for a general pentagonal truncation,
  illustrated in Figure [figure:2], defined by three parameters:
  :math:`M, K`, and :math:`N`, where :math:`M` is defined above,
  :math:`K` is the highest degree of the associated Legendre
  polynomials, and :math:`N` is the highest degree of the Legendre
  polynomials for :math:`m= 0`. The common truncations are subsets of
  this pentagonal case:

  .. math::
     :label:

     \mathrm{Triangular:}  \; & M = N = K , \nonumber \\ 
     \mathrm{Rhomboidal:}  \; & K = N + M ,  \\ 
     \mathrm{Trapezoidal:} \; & N = K > M . \nonumber

The quantity :math:`{\cal N} ( m )` in ([3.b.1])
represents an arbitrary limit on the two-dimensional wavenumber
:math:`n`, and for the pentagonal truncation described above is
simply given by :math:`{\cal N} ( m ) = \min(N + \vert m \vert , K )`.

.. todo:: put in figure3-2 Pentagonal truncation parameters

The associated Legendre polynomials used in the model are normalized
such that

.. math:: \int^1_{-1} [P^m_n(\mu)]^2 d\mu = 1 .  \label{3.b.3}

With this normalization, the Coriolis parameter :math:`f` is

.. math:: f = \frac{\Omega}{\sqrt {0.375}} P^o_1 , \label{3.b.4}

 which is required for the absolute vorticity.

The coefficients of the spectral representation ([3.b.1]) are given by

.. math::

   \psi^m_n = \int^1_{-1} \frac{1}{2 \pi} \int^{2 \pi}_0 \psi (\lambda,
   \mu) e^{-im\lambda} d \lambda P^m_n (\mu) d \mu . \label{3.b.5}

The inner integral represents a Fourier transform,

.. math::

   \psi^m (\mu) = \frac{1}{2 \pi} \int^{2 \pi}_0 \psi (\lambda, \mu)
   e^{-im\lambda} d \lambda , \label{3.b.6}

which is performed by a Fast Fourier Transform (FFT) subroutine. The
outer integral is performed via Gaussian quadrature,

.. math::

   \psi^m_n = \sum^J_{j=1} \psi^m (\mu_j) P^m_n (\mu_j) w_j ,
   \label{3.b.7}

where :math:`\mu_j` denotes the Gaussian grid points in the meridional
direction, :math:`w_j` the Gaussian weight at point :math:`\mu_j`, and
:math:`J` the number of Gaussian grid points from pole to pole. The
Gaussian grid points (:math:`\mu_j`) are given by the roots of the
Legendre polynomial :math:`P_J(\mu)`, and the corresponding weights are
given by

.. math::

   w_j = \frac{2(1 - \mu_j^2)}{[J\; P_{J-1} (\mu_j) ]^2}
   . \label {3.b.8}

The weights themselves satisfy

.. math:: \sum^J_{j=1} w_j = 2.0 \ . \label{3.b.9}

The Gaussian grid used for the north–south transformation is generally
chosen to allow un-aliased computations of quadratic terms only. In this
case, the number of Gaussian latitudes :math:`J` must satisfy

.. math::
   :label: 3.b.11

   J & \geq (2N + K + M + 1)/2  \quad \hbox{for}\> M  \leq 2(K - N)\>, \\ 
   J & \geq (3K + 1)/2  \quad \hbox{for}\> M  \geq 2(K - N)\> . 

For the common truncations, these become

.. math::
   :label: 3.b.13

   J & \geq (3K + 1)/2  \quad  \mathrm{for \; triangular \; and \;  trapezoidal} ,  \\ 
   J & \geq (3N + 2M + 1)/2  \quad \mathrm{for \; rhomboidal} . 

In order to allow exact Fourier transform of quadratic terms, the
number of points :math:`P` in the east–west direction must satisfy

.. math:: P \geq 3M + 1 \ . \label{3.b.14}

The actual values of :math:`J` and :math:`P` are often not set equal to
the lower limit in order to allow use of more efficient transform
programs.

Although in the next section of this model description, we continue to
indicate the Gaussian quadrature as a sum from pole to pole, the code
actually deals with the symmetric and antisymmetric components of
variables and accumulates the sums from equator to pole only. The model
requires an even number of latitudes to easily use the symmetry
conditions. This may be slightly inefficient for some spectral
resolutions. We define a new index, which goes from :math:`-I` at the
point next to the south pole to :math:`+I` at the point next to the
north pole and not including 0 (there are no points at the equator or
pole in the Gaussian grid), *i.e.,* let :math:`I = J/2` and
:math:`i = j - J/2` for :math:`j
\geq J/2+1` and :math:`i = j - J/2 - 1` for :math:`j \leq J/2`; then the
summation in ([3.b.7]) can be rewritten as

.. math::

   \psi^m_n = \sum \limits^{I}_{i = -I, \;i \neq 0} \psi^m (\mu_i) P^m_n
   (\mu_i) w_i .
   \label{3.b.15}

The symmetric (even) and antisymmetric (odd) components of
:math:`\psi^m` are defined by

.. math::

   (\psi_E)^m_i = \frac{1}{2} ( \psi^m_i + \psi^m_{-i}
   ),

.. math::

   (\psi_O)^m_i = \frac{1}{2} (\psi^m_i - \psi^m_{-i}
   ) .
   \label{3.b.16}

Since :math:`w_i` is symmetric about the equator, ([3.b.15]) can be
rewritten to give formulas for the coefficients of even and odd
spherical harmonics:

.. math::

   \psi^{m}_{n} =
   \begin{cases}
    \sum \limits ^I_{i=1} (\psi_E)^m_i(\mu_i) P^m_n(\mu_i)
    2w_i & \text{for $n-m$ even},\\ \sum \limits ^I_{i=1}
    (\psi_O)^m_i (\mu_i) P^m_n (\mu_i) 2w_i & \text{for $n-m$
    odd}.
   \end{cases}
   \label{3.b.17}

The model uses the spectral transform method for all nonlinear terms.
However, the model can be thought of as starting from grid–point values
at time :math:`t` (consistent with the spectral representation) and
producing a forecast of the grid–point values at time
:math:`t + \Delta t` (again, consistent with the spectral resolution).
The forecast procedure involves computation of the nonlinear terms
including physical parameterizations at grid points; transformation via
Gaussian quadrature of the nonlinear terms from grid–point space to
spectral space; computation of the spectral coefficients of the
prognostic variables at time :math:`t + \Delta t` (with the implied
spectral truncation to the model resolution); and transformation back to
grid–point space. The details of the equations involved in the various
transformations are given in the next section.

Combination of terms
~~~~~~~~~~~~~~~~~~~~

In order to describe the transformation to spectral space, for each
equation we first group together all undifferentiated explicit terms,
all explicit terms with longitudinal derivatives, and all explicit terms
with meridional derivatives appearing in the Dyn operator. Thus, the
vorticity equation ([3.a.94]) is rewritten

.. math::

   \underline{(\zeta + f )}^{n+1} = \underline{{\hbox{\sffamily\slshape V}}} +
   \frac{1} {a(1 - \mu^2)} [ \frac{\partial}{\partial \lambda}
   (\underline{{\hbox{\sffamily\slshape V}}}_\lambda) - (1 - \mu^2) \frac{\partial}{\partial
   \mu} (\underline{{\hbox{\sffamily\slshape V}}}_\mu) ] , \label{3.b.18}

where the explicit forms of the vectors
:math:`\underline{{\hbox{\sffamily\slshape V}}},
\underline{{\hbox{\sffamily\slshape V}}}_\lambda,` and
:math:`\underline{{\hbox{\sffamily\slshape V}}}_\mu` are given as

.. math::
   :label: A.1-3

   \underline{{\hbox{\sffamily\slshape V}}} & = \underline{(\zeta+f)}^{n-1} ,  \\
   \underline{{\hbox{\sffamily\slshape V}}}_\lambda & = 2 \Delta t\, \underline{n}^{n}_{V} ,  \\
   \underline{{\hbox{\sffamily\slshape V}}}_{\mu} & = 2 \Delta t\,\underline{n}^{n}_{U}. 

The divergence equation ([3.a.95]) is

.. math::

   \underline{\delta}^{n+1} & = \underline{{\hbox{\sffamily\slshape D}}} + \frac{1}{a(1 -
   \mu^2)} [ \frac{\partial}{\partial \lambda}
   (\underline{{\hbox{\sffamily\slshape D}}}_\lambda) + (1 - \mu^2) \frac{\partial}{\partial
   \mu} (\underline{{\hbox{\sffamily\slshape D}}}_\mu) ] - \nabla^2
   \underline{{\hbox{\sffamily\slshape D}}}_\nabla \nonumber \\ & \phantom{=} - \Delta t
   \nabla^2 (R {\mathbf{H}}^r \underline{T}^{\prime \, n+1} + R
   (\underline{b}^r + \underline{h}^r ) \Pi^{n+1}) .
   \label{3.b.19}

The mean component of the temperature is not included in the
next–to–last term since the Laplacian of it is zero. The thermodynamic
equation ([3.a.96]) is

.. math::

   \underline{T}^{\prime \, n+1} = \underline{{\hbox{\sffamily\slshape T}}} - \frac{1}{a (1 -
   \mu^2)} [ \frac{\partial}{\partial \lambda}
   (\underline{{\hbox{\sffamily\slshape T}}}_\lambda) + (1
   - \mu^2) \frac{\partial}{\partial \mu} (\underline{{\hbox{\sffamily\slshape T}}}_\mu)
   - ] -
   \Delta t {\mathbf{D}}^r \; \underline{\delta}^{n+1} . \label{3.b.20}

The surface–pressure tendency ([3.a.97]) is

.. math::

   \Pi^{n+1} = {{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}} - \frac{\Delta t}{\pi^r} (
   \underline{\Delta p}^r )^T \underline{\delta}^{n+1}
   . \label{3.b.21}

The grouped explicit terms in ([3.b.19])–([3.b.21]) are given as
follows. The terms of ([3.b.19]) are

.. math::
   :label: A.4-6

   \underline{{\hbox{\sffamily\slshape D}}} & = \underline{\delta}^{n-1} ,  \\
   \underline{{\hbox{\sffamily\slshape D}}}_\lambda & = 2 \Delta t \, \underline{n}^{n}_{U} ,  \\
   \underline{{\hbox{\sffamily\slshape D}}}_\mu & = 2 \Delta t \, \underline{n}^{n}_{V} , 

.. math::

   \nonumbereqn{\underline{{\hbox{\sffamily\slshape D}}}_\nabla = 2 \Delta t [ \underline{E}^n +
   \Phi_s \underline{1} + R {\mathbf{H}}^{r} \underline{{\cal T}}^{'n} ]} \\[1ex]
   & & \mbox{} + \Delta t [ R {\mathbf{H}}^{r} ( {( \underline{T}^{'} )}^{n-1} 
   - 2 {(\underline{T}^\prime)}^n )
   + R ( \underline{b}^{r} + \underline{h}^{r}
   ) ( {\Pi}^{n-1} - 2 {\Pi}^{n} ) ] \ . \label{A.7}

The terms of ([3.b.20]) are

.. math::

   \underline{{\hbox{\sffamily\slshape T}}} & = (\underline{T}')^{n-1} + 2 \Delta t \,
   \underline{\Gamma}^{n} \, - \Delta t {\mathbf{D}}^{r}
   [\underline{\delta}^{n-1} - 2\underline{\delta}^{n} ] \ ,
   \label{A.8} \\
   \underline{{\hbox{\sffamily\slshape T}}}_{\lambda} & = 2 \Delta t \underline{(UT')}^n , 
   \label{A.9} \\
   \underline{{\hbox{\sffamily\slshape T}}}_\mu & = 2 \Delta t  \underline{(VT')}^n . 
   \label{A.10}

The nonlinear term in ([3.b.21]) is

.. math::

   & \nonumber {{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}} = \Pi^{n-1} - 2 \Delta t \frac{1}{\pi^n} [ (
   \underline{\delta}^n )^T ( \underline{\Delta p}^n ) +
   ( \underline{\mathbf{V}}^n )^T \nabla \Pi^n \pi^n
   \underline{\Delta B} ]& \\
   & \mbox{} - \Delta t [ (\underline{\Delta p}^r )^T \frac{1}{\pi^r}
   ] [ \underline{\delta}^{n-1} -
   2 \underline{\delta}^n ] \ .& 
   \label{A.11}

Transformation to spectral space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Formally, Equations ([3.b.18])-([3.b.21]) are transformed to spectral
space by performing the operations indicated in ([3.b.22]) to each term.
We see that the equations basically contain three types of terms, for
example, in the vorticity equation the undifferentiated term
:math:`\underline{{\hbox{\sffamily\slshape V}}}`, the longitudinally
differentiated term
:math:`\underline{{\hbox{\sffamily\slshape V}}}_\lambda`, and the
meridionally differentiated term
:math:`\underline{{\hbox{\sffamily\slshape V}}}_\mu`. All terms in the
original equations were grouped into one of these terms on the Gaussian
grid so that they could be transformed at once.

Transformation of the undifferentiated term is obtained by
straightforward application of ([3.b.5])-([3.b.7]),

.. math::

   \{ \underline{{\hbox{\sffamily\slshape V}}} \}^m_n = \sum^J_{j=1}
   \underline{{\hbox{\sffamily\slshape V}}}^m(\mu_j) P^m_n(\mu_j) w_j , \label{3.b.22}

where :math:`\underline{{\hbox{\sffamily\slshape V}}}^m(\mu_j)` is the
Fourier coefficient of :math:`\underline{{\hbox{\sffamily\slshape V}}}`
with wavenumber :math:`m` at the Gaussian grid line :math:`\mu_j`. The
longitudinally differentiated term is handled by integration by parts,
using the cyclic boundary conditions,

.. math::

   \{ \frac{\partial}{\partial \lambda}
   (\underline{{\hbox{\sffamily\slshape V}}}_\lambda) \}^m & = \frac{1}{2 \pi} \int
   ^{2\pi}_o \frac{\partial \underline{{\hbox{\sffamily\slshape V}}}_\lambda} {\partial
   \lambda} e^{-im\lambda} d \lambda ,\\ & = im \frac{1}{2 \pi} \int^{2
   \pi}_o \underline{{\hbox{\sffamily\slshape V}}}_\lambda e^{-im\lambda} d\lambda , \\
   \label{3.b.23}

so that the Fourier transform is performed first, then the
differentiation is carried out in spectral space. The transformation to
spherical harmonic space then follows ([3.b.25]):

.. math::

    \{ \frac{1}{a(1 - \mu^2)} \frac{\partial}{\partial \lambda}
   (\underline{{\hbox{\sffamily\slshape V}}}_\lambda)  \}^m_n = im \sum^J_{j=1}
   \underline{{\hbox{\sffamily\slshape V}}}_\lambda^m (\mu_j) \frac{P^m_n(\mu_j)}{a(1 -
   \mu^2_j)} w_j , \label{3.b.24}

where
:math:`\underline{{\hbox{\sffamily\slshape V}}}_\lambda^m (\mu_j)` is
the Fourier coefficient of
:math:`\underline{{\hbox{\sffamily\slshape V}}}_\lambda` with wavenumber
:math:`m` at the Gaussian grid line :math:`\mu_j`.

The latitudinally differentiated term is handled by integration by parts
using zero boundary conditions at the poles:

.. math::

    \{ \frac{1}{a(1 - \mu^2)} (1 - \mu^2) \frac{\partial}{\partial
   \mu} (\underline{{\hbox{\sffamily\slshape V}}}_\mu)  \}^m_n & = \int^1_{-1}
   \frac{1}{a(1 - \mu^2)} (1 - \mu^2) \frac{\partial} {\partial \mu}
   (\underline{{\hbox{\sffamily\slshape V}}}_\mu)^m P^m_n d\mu ,\\ & = - \int^1_{-1}
   \frac{1}{a(1 - \mu^2)} (\underline{{\hbox{\sffamily\slshape V}}}_\mu)^m (1 - \mu^2)
   \frac{dP^m_n}{d\mu} d\mu .
      \label{3.b.25}

Defining the derivative of the associated Legendre polynomial by

.. math:: H^m_n = (1 - \mu^2) \frac{dP^m_n}{d\mu} , \label{3.b.26}

 ([3.b.28]) can be written

.. math::

    \{ \frac{1}{a(1 - \mu^2)} (1 - \mu^2) \frac{\partial}{\partial
   \mu} (\underline{{\hbox{\sffamily\slshape V}}}_\mu)  \}^m_n = - \sum^J_{j=1}
   (\underline{{\hbox{\sffamily\slshape V}}}_\mu)^m \frac{H^m_n(\mu_j)}{a(1 - \mu^2_j)} w_j
   . \label {3.b.27}

Similarly, the :math:`\nabla^2` operator in the divergence equation can
be converted to spectral space by sequential integration by parts and
then application of the relationship

.. math::

   \nabla^2 P^m_n( \mu) e^{im \lambda} = \frac{-n (n+1)}{a^2} P^m_n(\mu)
   e^{im \lambda}, \label{3.b.28}

 to each spherical harmonic function individually so that

.. math::

    \{ \nabla^2 \underline{{\hbox{\sffamily\slshape D}}}_\nabla  \}^m_n =
   \frac{-n(n+1)} {a^2} \sum^J_{j=1} \underline{{\hbox{\sffamily\slshape D}}}_\nabla^m
   (\mu_j) P^m_n(\mu_j) w_j ,
   \label{3.b.29}

where :math:`\underline{{\hbox{\sffamily\slshape D}}}_\nabla^m (\mu)`
is the Fourier coefficient of the original grid variable
:math:`\underline{{\hbox{\sffamily\slshape D}}}_\nabla`.

Solution of semi-implicit equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The prognostic equations can be converted to spectral form by summation
over the Gaussian grid using ([3.b.22]), ([3.b.24]), and ([3.b.27]). The
resulting equation for absolute vorticity is

.. math::

   \underline{(\zeta + f )}^m_n = \underline
   {{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_n^m ,\label{3.b.30}

where :math:`\underline{(\zeta + f )}_n^m` denotes a
spherical harmonic coefficient of
:math:`\underline{(\zeta + f )}^{n+1}`, and the form of
:math:`\underline{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_n^m`,
as a summation over the Gaussian grid, is given as

.. math::

   \underline{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}^m_n = \sum^J_{j=1}  [ \underline {{\hbox{\sffamily\slshape V}}}^m(\mu_j) P^m_n 
   (\mu_j) + im \underline{{\hbox{\sffamily\slshape V}}}^m_\lambda (\mu_j) \frac{P^m_n (\mu_j)}{a(1 - 
   \mu^2_j)} + \underline{{\hbox{\sffamily\slshape V}}}^m_\mu (\mu_j) \frac{H^m_n (\mu_j)}{a(1 - \mu^2_j)}  
   ] w_j . \label{A.12}

The spectral form of the divergence equation ([3.b.19]) becomes

.. math::

   \underline{\delta}^m_n = \underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}^m_n + \Delta t
   \frac{n(n+1)} {a^2} [ R {\mathbf{H}}^r \underline{T}^{\prime \,
   m}_n + R ( \underline{b}^r + \underline{h}^r ) \Pi^m_n
   ] ,
   \label{3.b.31}

where :math:`\underline{\delta}^m_n, \;\underline{T}'^m_n, \;` and
:math:`\Pi^m_n` are spectral coefficients of
:math:`\underline{\delta}^{n+1}, \;
\underline{T}^{\prime \, n+1}`, and :math:`\Pi^{n+1}`. The Laplacian of
the total temperature in ([3.b.19]) is replaced by the equivalent
Laplacian of the perturbation temperature in ([3.b.31]).
:math:`\underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}^m_n`
is given by

.. math::

   \nonumber\underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}^m_n = \sum^J_{j=1} \Biggl \{ [ \underline {{\hbox{\sffamily\slshape D}}}^m 
   (\mu_j) + \frac{n(n+1)}{a^2} \underline {{\hbox{\sffamily\slshape D}}}^m_\nabla (\mu_j) ] 
   P^m_n (\mu_j) \\
   + im \underline {{\hbox{\sffamily\slshape D}}}^m_\lambda (\mu_j) \frac{P^m_n (\mu_j)} 
   {a(1 - \mu_j^2)} - \underline {{\hbox{\sffamily\slshape D}}}^m_\mu (\mu_j) \frac{H^m_n (\mu_j)} 
   {a(1 - \mu^2_j)} \Biggr \} w_j . \label{A.13}

The spectral thermodynamic equation is

.. math::

   \underline{T}'^m_n = \underline{{\hbox{\sffamily\slshape T}}{\hbox{\sffamily\slshape S}}}^m_n - \Delta t
   {\mathbf{D}}^r \underline{\delta}^m_n , \label{3.b.32}

with :math:`\underline{{\hbox{\sffamily\slshape T}}{\hbox{\sffamily\slshape S}}}^m_n`
defined as

.. math::

   \underline{{\hbox{\sffamily\slshape T}}{\hbox{\sffamily\slshape S}}}^m_n = \sum^J_{j=1} [ \underline{{\hbox{\sffamily\slshape T}}}^m (\mu_j) P^m_n
   (\mu_j) - im \underline{{\hbox{\sffamily\slshape T}}}^m_\lambda (\mu_j) \frac{P^m_n (\mu_j)}{a (1 
   - \mu^2_j)} + \underline{{\hbox{\sffamily\slshape T}}}^m_\mu (\mu_j) \frac{H^m_n (\mu_j)}{a (1 -
   \mu^2_j)}  ] w_j , \label{A.14}

while the surface pressure equation is

.. math::

   \Pi^m_n = {{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}}^m_n - \underline{\delta}^m_n
   (\underline{\Delta p}^r )^T \frac{\Delta t}{\pi^r} ,
   \label{3.b.33}

where :math:`{{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}}^m_n`
is given by

.. math:: {{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}}^m_n = \sum^J_{j=1} {{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}}^m (\mu_j) P^m_n (\mu_j) w_j . \label{A.15}

Equation ([3.b.30]) for vorticity is explicit and complete at this
point. However, the remaining equations ([3.b.31])–([3.b.33]) are
coupled. They are solved by eliminating all variables except
:math:`\underline{\delta}^m_n`:

.. math::

   {\mathbf{A}}_n \underline{\delta}^m_n & =
   \underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}^m_n + \Delta t \frac{n(n+1)}{a^2} [
   R {\mathbf{H}}^r (\underline{{\hbox{\sffamily\slshape T}}{\hbox{\sffamily\slshape S}}})^m_n + R (
   \underline{b}^r + \underline{h}^r ) ({{\hbox{\sffamily\slshape P}}{\hbox{\sffamily\slshape S}}})^m_n
    ] , \label {3.b.34} \\[-1.0em]
   \intertext{where}\nonumber\\[-2.0em] {\mathbf{A}}_n & = {\mathbf{I}}
   + \Delta t^2 \frac{n(n+1)}{a^2}  [ R {\mathbf{H}}^r
   {\mathbf{D}}^r + R ( \underline{b}^r + \underline{h}^r )
   ( ( \underline{\Delta p^r} )^T \, \frac{1}{\pi^r}
   )  ] , \label{3.b.35}

which is simply a set of :math:`K` simultaneous equations for the
coefficients with given wavenumbers (:math:`m,n`) at each level and is
solved by inverting :math:`{\mathbf{A}}_n`. In order to prevent
the accumulation of round–off error in the global mean divergence (which
if exactly zero initially, should remain exactly zero)
:math:`({\mathbf{A}}_o)^{-1}` is set to the null matrix
rather than the identity, and the formal application of ([3.b.34]) then
always guarantees :math:`\underline{\delta}^o_o = 0`. Once
:math:`\delta^m_n` is known, :math:`\underline{T}'^m_n` and
:math:`\Pi^m_n` can be computed from ([3.b.32]) and ([3.b.33]),
respectively, and all prognostic variables are known at time
:math:`n\!+\!1` as spherical harmonic coefficients. Note that the mean
component :math:`\underline{T}'^o_o` is not necessarily zero since the
perturbations are taken with respect to a specified
:math:`\underline{T}^r`.

Horizontal diffusion
~~~~~~~~~~~~~~~~~~~~

As mentioned earlier, the horizontal diffusion in ([3.a.87]) and
([3.a.90]) is computed implicitly via time splitting after the
transformations into spectral space and solution of the semi-implicit
equations. In the following, the :math:`\zeta` and :math:`\delta`
equations have a similar form, so we write only the :math:`\delta`
equation:

.. math::

   (\delta^*)^m_n & = (\delta^{n+1})^m_n - (
   -1
   )^i 2 \Delta t K^{(2i)} [ \nabla^{2i}
       (\delta^*)^m_n - ( - 1 )^i
       (\delta^*)^m_n (2/a^2 )^i ] ,
      \label{3.b.36} \\
   (T^*)^m_n & = (T^{n+1})^m_n - ( -1)^i
   2 \Delta t K^{(2i)} [ \nabla^{2i} (T^*)^m_n ] \
       .
      \label{3.b.37} 

The extra term is present in ([3.b.36]), ([3.b.40]) and ([3.b.42]) to
prevent damping of uniform rotations. The solutions are just

.. math::

   (\delta^*)^m_n & = K^{(2i)}_n (\delta)
   (\delta^{n+1})^m_n , \label{3.b.38} \\ (T^*)^m_n
   & = K^{(2i)}_n (T ) (T^{n+1})^m_n ,
   \label{3.b.39} \\ K^{(2)}_n (\delta ) & = \{ 1 + 2
   \Delta t D_n K^{(2)} [ ( \frac{n(n + 1)}{a^2} ) -
   \frac{2}{a^2} ] \}^{-1} \ , \label{3.b.40} \\ K^{(2)}_n
   ( T ) & = \{ 1 + 2\Delta t D_n K^{(2)} (
   \frac{n(n + 1)}{a^2} ) \}^{-1} \ , \label{3.b.41} \\
   K^{(4)}_n (\delta) & = \{ 1 + 2 \Delta t D_n K^{(4)}
       [ ( \frac{n(n+1)}{a^2} )^2 - \frac{4}{a^4}
       ] \}^{-1} , \label{3.b.42} \\
   K^{(4)}_n (T ) & = \{ 1 + 2 \Delta t D_n K^{(4)} (
       \frac{n(n+1)}{a^2} )^2 \}^{-1} . \label{3.b.43}

:math:`K^{(2)}_n (\delta )` and
:math:`K^{(4)}_n (\delta )` are both set to 1 for :math:`n` =
0. The quantity :math:`D_n` represents the “Courant number limiter”,
normally set to 1. However, :math:`D_n` is modified to ensure that the
CFL criterion is not violated in selected upper levels of the model. If
the maximum wind speed in any of these upper levels is sufficiently
large, then :math:`D_n = 1000` in that level for all :math:`n > n_c`,
where :math:`n_c = a \Delta t \big/ \max \vert
{\mathbf{V}} \vert`. This condition is applied whenever the wind
speed is large enough that :math:`n_c < K`, the truncation parameter in
([3.b.2]), and temporarily reduces the effective resolution of the model
in the affected levels. The number of levels at which this “Courant
number limiter” may be applied is user-selectable, but it is only used
in the top level of the 26 level control runs.

The diffusion of :math:`T` is not complete at this stage. In order to
make the partial correction from :math:`\eta` to :math:`p` in ([3.a.81])
local, it is not included until grid–point values are available. This
requires that :math:`\nabla^4 \Pi` also be transformed from spectral to
grid–point space. The values of the coefficients :math:`K^{(2)}` and
:math:`K^{(4)}` for the standard T42 resolution are
:math:`2.5 \times 10^5`\ m\ :math:`^2`\ sec\ :math:`^{-1}` and
:math:`1.0
\times 10^{16}`\ m\ :math:`^4`\ sec\ :math:`^{-1}`, respectively.

Initial divergence damping
~~~~~~~~~~~~~~~~~~~~~~~~~~

Occasionally, with poorly balanced initial conditions, the model
exhibits numerical instability during the beginning of an integration
because of excessive noise in the solution. Therefore, an optional
divergence damping is included in the model to be applied over the first
few days. The damping has an initial e-folding time of :math:`\Delta t`
and linearly decreases to 0 over a specified number of days,
:math:`t_D`, usually set to be 2. The damping is computed implicitly via
time splitting after the horizontal diffusion.

.. math::

   r & = \max [ \frac{1}{\Delta t} (t_D - t) / t_D , ~ 0 ]
   \label{3.b.44} \\
   (\delta^*)^m_n & = \frac{1}{1 + 2\Delta t r}
   (\delta^*)^m_n
   \label{3.b.45}

Transformation from spectral to physical space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the prognostic variables are completed at time :math:`n+1` in
spectral space :math:`(\underline{(\zeta + f)^*})^m_n`,
:math:`(\underline{\delta^*})^m_n`,
:math:`(\underline{T^*})^m_n`,
:math:`(\Pi^{n+1})^m_n` they are transformed to grid space.
For a variable :math:`\psi`, the transformation is given by

.. math::

   \psi( \lambda, \mu) = \sum^M_{m=-M}  [ \sum^{{\cal N}(m)}_{n=|m|}
   \psi^m_n P^m_n (\mu)  ] e^{im\lambda} . \label{3.b.46}

The inner sum is done essentially as a vector product over :math:`n`,
and the outer is again performed by an FFT subroutine. The term needed
for the remainder of the diffusion terms, :math:`\nabla^4 \Pi`, is
calculated from

.. math::

   \nabla^4 \Pi^{n+1} = \sum^M_{m=-M}  [\sum^{{\cal N}(m)}_{n=|m|}
   ( \frac{n(n+1)}{a^2}  )^2 (\Pi^{n+1} )^m_n P^m_n
   (\mu)  ] e^{im\lambda} . \label{3.b.47}

In addition, the derivatives of :math:`\Pi` are needed on the grid for
the terms involving :math:`\nabla \Pi` and
:math:`{\mathbf{V}} \cdot \nabla \Pi`,

.. math::

   {\mathbf{V}} \cdot \nabla \Pi = \frac{U}{a(1 - \mu^2)} \frac{\partial
   \Pi} {\partial \lambda} + \frac{V}{a(1 - \mu^2)} (1 - \mu^2)
   \frac{\partial \Pi}{\partial \mu}. \label{3.b.48}

These required derivatives are given by

.. math::

   \frac{\partial \Pi}{\partial \lambda} & = \sum^M_{m=-M} im [
   \sum^{{\cal N}(m)}_{n=\vert m \vert} \Pi^m_n P^m_n (\mu) ]
   e^{im\lambda} , \label{3.b.49} \\ \intertext{and using (\ref{3.b.26}),
   }\nonumber\\[-2.0em] (1 - \mu^2) \frac{\partial \Pi}{\partial \mu} & =
   \sum^M_{m=-M}  [ \sum^{{\cal N}(m)}_{n=\vert m \vert} \Pi^m_n
   H^m_n (\mu)  ] e^{im \lambda} , \label{3.b.50}

which involve basically the same operations as ([3.b.47]). The other
variables needed on the grid are :math:`U` and :math:`V`. These can be
computed directly from the absolute vorticity and divergence
coefficients using the relations

.. math::

   (\zeta + f )^m_n & = - \frac{n(n+1)}{a^2} \psi^m_n + f^m_n
   ,
   \label{3.b.51}  \\
   \delta^m_n & = - \frac{n(n+1)}{a^2} \chi^m_n , \label{3.b.52}

in which the only nonzero :math:`f^m_n` is
:math:`f^o_1 = \Omega/ \sqrt{.375},` and

.. math::

   U & = \frac{1}{a} \frac{\partial \chi}{\partial \lambda} - \frac{(1 -
   \mu^2)}{a} \frac{\partial \psi}{\partial \mu} , \label{3.b.53} \\ V &
   = \frac{1}{a} \frac{\partial \psi}{\partial \lambda} + \frac{(1 -
   \mu^2)} {a} \frac{\partial \chi}{\partial \mu} . \label{3.b.54}

 Thus, the direct transformation is

.. math::

   U & = - \sum^M_{m=-M} a \sum^{{\cal N}(m)}_{n=|m|}  [ \frac{im}
   {n(n+1)} \delta^m_n P^m_n (\mu) - \frac{1}{n(n+1)} (\zeta + f)^m_n
   H^m_n (\mu) ] e^{im \lambda} \nonumber \\ & \phantom{=} -
   \>\frac{a}{2} \frac{\Omega}{\sqrt{0.375}} H^o_1 , \label{3.b.55} \\ V
   & = - \sum^M_{m=-M} a \sum^{{\cal N}(m)}_{n=|m|} \bigg
   [\frac{im}{n(n+1)} (\zeta + f)^m_n P^m_n (\mu) + \frac{1}{n(n+1)}
   \delta^m_n H^m_n (\mu) \bigg ] e^{im\lambda}.  \label{3.b.56}

The horizontal diffusion tendencies are also transformed back to grid
space. The spectral coefficients for the horizontal diffusion tendencies
follow from ([3.b.36]) and ([3.b.37]):

.. math::

   F_{T_H}(T^*)^m_n & = (-1)^{i+1} K^{2i} [
       \nabla^{2i} (T^*) ]^m_n , \label{3.b.57} \\
   F_{\zeta_H} ((\zeta + f )^* )^m_n & =
   (-1)^{i+1} K^{2i}  \{\nabla^{2i} (\zeta +
   f)^* - (-1)^i (\zeta + f )^* (2/a^2
   )^i  \} , \label{3.b.58} \\ F_{\delta_H}
   (\delta^*)^m_n & = (-1) K^{2i} \{
   \nabla^{2i} (\delta ^*) - (-1)^i \delta^*
   (2/a^2 )^i \}, \label{3.b.59}

using :math:`i = 1` or 2 as appropriate for the :math:`\nabla^2` or
:math:`\nabla^4` forms. These coefficients are transformed to grid space
following ([3.b.1]) for the :math:`T` term and ([3.b.55]) and ([3.b.56])
for vorticity and divergence. Thus, the vorticity and divergence
diffusion tendencies are converted to equivalent :math:`U` and :math:`V`
diffusion tendencies.

Horizontal diffusion correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After grid–point values are calculated, frictional heating rates are
determined from the momentum diffusion tendencies and are added to the
temperature, and the partial correction of the :math:`\nabla^4`
diffusion from :math:`\eta` to :math:`p` surfaces is applied to
:math:`T`. The frictional heating rate is calculated from the kinetic
energy tendency produced by the momentum diffusion

.. math::

   F_{F_H} = -u^{n-1} F_{u_H} (u^*)/c_p^* -v^{n-1} F_{v_H} (v^*)/c^*_p ,
   \label{3.b.60}

where :math:`F_{u_H}`, and :math:`F_{v_H}` are the momentum equivalent
diffusion tendencies, determined from :math:`F_{\zeta_H}` and
:math:`F_{\delta_H}` just as :math:`U` and :math:`V` are determined from
:math:`\zeta` and :math:`\delta`, and

.. math::

   c^*_p = c_p [ 1 + (\frac{c_{p_v}}{c_p} -1 ) q^{n+1}
    ] . \label{3.b.61}

 These heating rates are then combined with the correction,

.. math::

   \hat {T}^{n+1}_k = T^*_k + ( 2 \Delta t F_{F_H})_k + 2
   \Delta t ( \pi B \frac{\partial T^*}{\partial p} )_k
   K^{(4)} \nabla^4 \Pi^{n+1} .  \label{3.b.62}

The vertical derivatives of :math:`T^*` (where the :math:`^*` notation
is dropped for convenience) are defined by

.. math::

   ( \pi B \frac{\partial T}{\partial p} )_1 & = \frac{\pi}
       {2\Delta p_1} [ B_{1+\frac{1}{2}} ( T_2 - T_1
       ) ] \ , \\
   (\pi B \frac{\partial T}{\partial p} )_k & = \frac{\pi}
       {2\Delta p_k} [ B_{k + \frac{1}{2}} ( T_{k+1} - T_k
       ) + B_{k - \frac{1}{2}} ( T_k - T_{k-1} )
       ] \ , \\
   (\pi B \frac{\partial T}{\partial p} )_K & = \frac{\pi}
       {2\Delta p_K} [ B_{K - \frac{1}{2}} ( T_K - T_{K-1}
       ) ] .  \label{3.b.63}

The corrections are added to the diffusion tendencies calculated earlier
([3.b.57]) to give the total temperature tendency for diagnostic
purposes:

.. math::

   \hat{F}_{T_H}(T^*)_k = F_{T_H}(T^*)_k + ( 2 \Delta t F_{F_H}
   )_k + 2 \Delta t B_k (\pi \frac{\partial T^*}{\partial p}
   )_k K^{(4)} \nabla^4 \Pi^{n+1} . \label{3.b.64}

Semi-Lagrangian Tracer Transport
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The forecast equation for water vapor specific humidity and constituent
mixing ratio in the :math:`\eta` system is from ([3.a.36]) excluding
sources and sinks.

.. math::

   \frac{dq}{dt} & = \frac{\partial q}{\partial t} + {\mathbf{V}} \cdot
       \nabla q + \dot\eta \, \frac{\partial p}{\partial \eta} \,
       \frac{\partial q} {\partial p} = 0 \label{3.c.1} \\[-1.0em]
   \intertext{or}\nonumber\\[-2.0em] \frac{dq}{dt}& = \frac{\partial
   q}{\partial t} + {\mathbf{V}} \cdot \nabla q + \dot\eta \,
   \frac{\partial q}{\partial \eta} = 0 . \label{3.c.2}

Equation ([3.c.2]) is more economical for the semi-Lagrangian vertical
advection, as :math:`\Delta \eta` does not vary in the horizontal, while
:math:`\Delta p` does. Written in this form, the :math:`\eta` advection
equations look exactly like the :math:`\sigma` equations.

The parameterizations are time-split in the moisture equation. The
tendency sources have already been added to the time level
:math:`( n -
1 )`. The semi-Lagrangian advection step is subdivided into
horizontal and vertical advection sub-steps, which, in an Eulerian form,
would be written

.. math::

   q^* & = q^{n-1} + 2\Delta t ( {\mathbf{V}} \cdot \nabla q
   )^n
   \label{3.c.3} \\[-1.0em]
   \intertext{and}\nonumber\\[-2.0em] q^{n+1} & = q^* + 2\Delta t (
   \dot\eta \frac{\partial q}{\partial n} )^n . \label{3.c.4}

In the semi-Lagrangian form used here, the general form is

.. math::

   q^* & = {\rm L}_{\lambda \varphi} ( q^{n-1} ) ,
   \label{3.c.5} \\ q^{n+1} & = {\rm L}_\eta ( q^* ) .
   \label{3.c.6}

Equation ([3.c.5]) represents the horizontal interpolation of
:math:`q^{n-1}` at the departure point calculated assuming
:math:`\dot\eta = 0`. Equation ([3.c.6]) represents the vertical
interpolation of :math:`q^*` at the departure point, assuming
:math:`{\mathbf{V}} = 0`.

The horizontal departure points are found by first iterating for the
mid-point of the trajectory, using winds at time :math:`n`, and a first
guess as the location of the mid-point of the previous time step

.. math::

   \lambda^{k+1}_M & = \lambda_A - \Delta t u^n ( \lambda^k_M,
       \varphi^k_M ) \big/a \cos \varphi^k_M , \label{3.c.7} \\
   \varphi^{k+1}_M & = \varphi_A - \Delta t v^n ( \lambda^k_M,
   \varphi^k_M )/a , \label{3.c.8}

where subscript :math:`A` denotes the arrival (Gaussian grid) point and
subscript :math:`M` the midpoint of the trajectory. The velocity
components at :math:`( \lambda_M^k, \varphi^k_M )` are
determined by Lagrange cubic interpolation. For economic reasons, the
equivalent Hermite cubic interpolant with cubic derivative estimates is
used at some places in this code. The equations will be presented later.

Once the iteration of ([3.c.7]) and ([3.c.8]) is complete, the departure
point is given by

.. math::

   \lambda_D & = \lambda_A - 2 \Delta t u^n ( \lambda_M, \varphi_M
       ) \big/a \cos \varphi_M , \label{3.c.9} \\
   \varphi_D & = \lambda_A - 2 \Delta t v^n ( \lambda_M, \varphi_M
       )/a , \label{3.c.10}

where the subscript :math:`D` denotes the departure point.

The form given by ([3.c.7])-([3.c.10]) is inaccurate near the poles and
thus is only used for arrival points equatorward of 70\ :math:`^\circ`
latitude. Poleward of 70\ :math:`^\circ` we transform to a local
geodesic coordinate for the calculation at each arrival point. The local
geodesic coordinate is essentially a rotated spherical coordinate system
whose equator goes through the arrival point. Details are provided in .
The transformed system is rotated about the axis through
:math:`( \lambda_A - \frac{\pi}{2}, 0
)` and :math:`( \lambda_A + \frac{\pi}{2}, 0 )`, by an
angle :math:`\varphi_A` so the equator goes through
:math:`( \lambda_A,
\varphi_A )`. The longitude of the transformed system is chosen to
be zero at the arrival point. If the local geodesic system is denoted by
:math:`( \lambda^\prime, \varphi^\prime )`, with velocities
:math:`( u^\prime, v^\prime )`, the two systems are related
by

.. math::

   \sin \phi^\prime & =  \sin \phi \cos \phi_A - \cos \phi \sin \phi_A
       \cos ( \lambda_A - \lambda ) , \label{3.c.11} \\
   \sin \phi & =  \sin \phi^\prime \cos \phi_A + \cos \phi^\prime \sin
   \prime_A \cos \lambda^\prime \ , \label{3.c.12} \\ \sin \lambda^\prime
   \cos \phi^\prime & =  -\sin ( \lambda_A
   -\lambda ) \cos \phi \ , \label{3.c.13} \\
   v^\prime \cos \phi^\prime & =  v [ \cos \phi \cos \phi_A + \sin
   \phi \sin \phi_A \cos ( \lambda_A - \lambda ) ]
   \nonumber \\ & \phantom{=} - u \sin \phi_A \sin ( \lambda_A -
   \lambda ), \label{3.c.14} \\ u^\prime \cos \lambda^\prime -
   v^\prime \sin \lambda^\prime \sin \phi^\prime & =  u \cos (
   \lambda_A - \lambda ) + v \sin \phi \sin ( \lambda_A -
   \lambda ) \ . \label{3.c.15}

The calculation of the departure point in the local geodesic system is
identical to ([3.c.7])-([3.c.10]) with all variables carrying a prime.
The equations can be simplified by noting that :math:`(
\lambda^\prime_A, \varphi_A^\prime ) = ( 0, 0 )` by
design and
:math:`u^\prime ( \lambda^\prime_A, \varphi_A^\prime )
= u( \lambda_A, \varphi_A )` and :math:`v^\prime (
\lambda^\prime_A, \varphi_A^\prime ) = v( \lambda_A,
\varphi_A )`. The interpolations are always done in global
spherical coordinates.

The interpolants are most easily defined on the interval 0 :math:`\leq
\theta \leq 1`. Define

.. math::

   \theta = ( x_D - x_i ) \big/ ( x_{i+1} - x_i ),
   \label{3.c.16}

where :math:`x` is either :math:`\lambda` or :math:`\varphi` and the
departure point :math:`x_D` falls within the interval
:math:`( x_i, x_{i+1} )`. Following (23) of with
:math:`r_i=3` the Hermite cubic interpolant is given by

.. math::

   q_D & =  q_{i+1} [ 3-2\theta ]\theta^2
    - d_{i+1} [ h_i \theta^2 ( 1-\theta ) ]
      \nonumber \\
   & \phantom{=} + q_i [ 3-2(1-\theta) ]
                                   (1-\theta )^2
    + d_i [ h_i \theta ( 1-\theta )^2 ]
   \label{3.c.17}

where :math:`q_i` is the value at the grid point :math:`x_i`,
:math:`d_i` is the derivative estimate given below, and
:math:`h_i = x_{i+1} - x_i`.

Following (3.2.12) and (3.2.13) of , the Lagrangian cubic polynomial
interpolant used for the velocity interpolation, is given by

.. math::

   f_{D\phantom{\dot D}} = \sum^2_{j=-1} \ell_j ( x_D )
   f_{i+j}
          \label{3.c.18}

where

.. math::

   \ell_j ( x_D ) =
          \frac{ ( x_D - x_{i-1} ) \ldots ( x_D -
              x_{i+j-1} ) ( x_D - x_{i+j+1} ) \ldots
              ( x_D - x_{i+2} ) }
            { ( x_{i+j} - x_{i-1} ) \ldots ( x_{i+j} -
              x_{i+j-1} ) ( x_{i+j} - x_{i+j+1} ) \ldots
              ( x_{i+j} - x_{i+2} ) } \label{3.c.19}

where :math:`f` can represent either :math:`u` or :math:`v`, or their
counterparts in the geodesic coordinate system.

The derivative approximations used in ([3.c.17]) for :math:`q` are
obtained by differentiating ([3.c.18]) with respect to :math:`x_D`,
replacing :math:`f` by :math:`q` and evaluating the result at
:math:`x_D` equal :math:`x_{i}` and :math:`x_{i+1}`. With these
derivative estimates, the Hermite cubic interpolant ([3.c.17]) is
equivalent to the Lagrangian ([3.c.18]). If we denote the four point
stencil :math:`(
x_{i-1},x_{i},x_{i+1},x_{i+2} )` by :math:`( x_1,x_2,x_3,x_4,
)` the cubic derivative estimates are

.. math::

   d_2 & = [ \frac{(x_2 - x_3)( x_2 - x_4 )} {(x_1 - x_2)(x_1 -
                      x_3)(x_1 - x_4)} ] q_1 \\
       & - [ \frac{1}{(x_1 - x_2)} - \frac{1}{(x_2 - x_3)}
               - \frac{1}{(x_2 - x_4)} ] q_2 \\
       & + [ \frac{(x_2 - x_1)(x_2 - x_4)} {(x_1 - x_3)(x_2 -
                       x_3)(x_3 - x_4)}] q_3 \\
       & - [ \frac{(x_2 - x_1)(x_2 - x_3)} {(x_1 - x_4)(x_2 -
                       x_4)(x_3 - x_4)} ] q_4
   \label{3.c.20}\\[-1.0em]
   \intertext{and}\nonumber\\[-2.0em]
   d_3 & = [ \frac{(x_3 - x_2)(x_3 - x_4)} {(x_1 - x_2)(x_1 -
                       x_3)(x_1 - x_4)} ] q_1 \\
       & - [ \frac{(x_3 - x_1)(x_3 - x_4)} {(x_1 - x_2)(x_2 -
                       x_3)(x_2 - x_4)} ] q_2 \\
       & - [ \frac{1}{(x_1 - x_3)} + \frac{1}{(x_2 - x_3)} -
               \frac{1}{(x_3 - x_4)} ] q_3 \\
       & - [ \frac{(x_3 - x_1)(x_3 - x_2)} {(x_1 - x_4)(x_2 -
                       x_4)(x_3 - x_4)} ] q_4
   \label{3.c.21}

The two dimensional :math:`( \lambda, \varphi )` interpolant
is obtained as a tensor product application of the one-dimensional
interpolants, with :math:`\lambda` interpolations done first. Assume the
departure point falls in the grid box :math:`(\lambda_i,\lambda_{i+1})`
and :math:`(\varphi_i,\varphi_{i+1})`. Four :math:`\lambda`
interpolations are performed to find :math:`q` values at
:math:`(\lambda_D,\varphi_{j-1})`, :math:`(\lambda_D,\varphi_j)`,
:math:`(\lambda_D,\varphi_{j+1})`, and :math:`(\lambda_D,
\varphi_{j+2})`. This is followed by one interpolation in
:math:`\varphi` using these four values to obtain the value at
:math:`(\lambda_D,
\varphi_D)`. Cyclic continuity is used in longitude. In latitude, the
grid is extended to include a pole point (row) and one row across the
pole. The pole row is set equal to the average of the row next to the
pole for :math:`q` and to wavenumber 1 components for :math:`u` and
:math:`v`. The row across the pole is filled with the values from the
first row below the pole shifted :math:`\pi` in longitude for :math:`q`
and minus the value shifted by :math:`\pi` in longitude for :math:`u`
and :math:`v`.

Once the departure point is known, the constituent value of :math:`q^* =
q^{n-1}_D` is obtained as indicated in ([3.c.5]) by Hermite cubic
interpolation ([3.c.17]), with cubic derivative estimates ([3.c.18]) and
([3.c.19]) modified to satisfy the Sufficient Condition for Monotonicity
with C\ :math:`^\circ` continuity (SCMO) described below. Define
:math:`\Delta_i q` by

.. math:: \Delta_i q = \frac{q_{i+1} - q_i}{x_{i+1} - x_i} \ . \label{3.c.28}

First, if :math:`\Delta_i q= 0` then

.. math:: d_i = d_{i+1} = 0 \ . \label{3.c.29}

 Then, if either

.. math:: 0 \leq \frac{d_i}{\Delta_i q} \leq 3 \label{3.c.30}

 or

.. math:: 0 \leq \frac{d_{i+1}}{\Delta_i q} \leq 3 \label{3.c.31}

is violated, :math:`d_i` or :math:`d_{i+1}` is brought to the
appropriate bound of the relationship. These conditions ensure that the
Hermite cubic interpolant is monotonic in the interval
:math:`[ x_i, x_{i+1}
]`.

The horizontal semi-Lagrangian sub-step ([3.c.5]) is followed by the
vertical step ([3.c.6]). The vertical velocity :math:`\dot \eta` is
obtained from that diagnosed in the dynamical calculations ([3.a.93]) by

.. math::

   ( \dot\eta )_{k+\frac{1}{2}} = ( \dot\eta \,
   \frac{\partial p}{\partial \eta} )_{k+\frac{1}{2}} \Bigg/
   (\frac{p_{k+1}
   -p_k}{\eta_{k+1} - \eta_k} ) , \label{3.c.32}

with :math:`\eta_k= A_k + B_k`. Note, this is the only place that the
model actually requires an explicit specification of :math:`\eta`. The
mid-point of the vertical trajectory is found by iteration

.. math::

   \eta^{k+1}_M = \eta_A - \Delta t \dot\eta^n ( \eta^k_M ) .
   \label{3.c.33}

Note, the arrival point :math:`\eta_A` is a mid-level point where
:math:`q` is carried, while the :math:`\dot\eta` used for the
interpolation to mid-points is at interfaces. We restrict :math:`\eta_M`
by

.. math:: \eta_1 \leq \eta_M \leq \eta_K , \label{3.c.34}

which is equivalent to assuming that :math:`q` is constant from the
surface to the first model level and above the top :math:`q` level. Once
the mid-point is determined, the departure point is calculated from

.. math::

   \eta_D = \eta_A - 2 \Delta t \dot\eta^n ( \eta_M ) , \label
   {3.c.35}

 with the restriction

.. math:: \eta_1 \leq \eta_D \leq \eta_K . \label{3.c.36}

The appropriate values of :math:`\dot\eta` and :math:`q` are determined
by interpolation ([3.c.17]), with the derivative estimates given by
([3.c.18]) and ([3.c.19]) for :math:`i = 2` to :math:`K - 1`. At the top
and bottom we assume a zero derivative (which is consistent with
([3.c.34]) and ([3.c.36])), :math:`d_i = 0` for the interval
:math:`k=1`, and :math:`\delta_{i+1} = 0` for the interval
:math:`k = K - 1`. The estimate at the interior end of the first and
last grid intervals is determined from an uncentered cubic
approximation; that is :math:`d_{i+1}` at the :math:`k=1` interval is
equal to :math:`d_i` from the :math:`k=2` interval, and :math:`d_i` at
the :math:`k=K-1` interval is equal to :math:`d_{i+1}` at the
:math:`k = K -2` interval. The monotonic conditions ([3.c.30]) to
([3.c.31]) are applied to the :math:`q` derivative estimates.

Mass fixers
~~~~~~~~~~~

This section describes original and modified fixers used for the
Eulerian and semi-Lagrangian dynamical cores.

Let :math:`\pi^0`, :math:`\Delta p^0` and :math:`q^0` denote the values
of air mass, pressure intervals, and water vapor specific humidity at
the beginning of the time step (which are the same as the values at the
end of the previous time step.)

:math:`\pi^+`, :math:`\Delta p^+` and :math:`q^+` are the values after
fixers are applied at the end of the time step.

:math:`\pi^-`, :math:`\Delta p^-` and :math:`q^-` are the values after
the parameterizations have updated the moisture field and tracers.

Since the physics parameterizations do not change the surface pressure,
:math:`\pi^-` and :math:`\Delta p^-` are also the values at the
beginning of the time step.

The fixers which ensure conservation are applied to the dry atmospheric
mass, water vapor specific humidity and constituent mixing ratios. For
water vapor and atmospheric mass the desired discrete relations,
following are

.. math::

   \int\limits_2 \, \pi^+ - \int\limits_3 \, q^+ \Delta p^+ & =
   {\mathbf{P}} , \label{3.d.1} \\ \int\limits_3 \, q^+ \Delta p^+ & =
   \int\limits_3 \, q^- \Delta p^- ,
    \label{3.d.2} 

where :math:`\mathbf{P}` is the dry mass of the atmosphere. From
the definition of the vertical coordinate,

.. math:: \Delta p = p_0 \Delta A + \pi \Delta B, \label{3.d.3}

and the integral :math:`\int\limits_2` denotes the normal Gaussian
quadrature while :math:`\int\limits_3` includes a vertical sum followed
by Gaussian quadrature. The actual fixers are chosen to have the form

.. math::

   \pi^+ ( \lambda, \varphi ) = {\mathbf{M}} \hat \pi^+
   ( \lambda, \varphi ) , \label{3.d.4}

preserving the horizontal gradient of :math:`\Pi`, which was calculated
earlier during the inverse spectral transform, and

.. math::

   q^+ ( \lambda, \varphi, \eta ) = \hat q^+ + \alpha\eta \hat
   q^+ \vert \hat q^+ - q^- \vert . \label{3.d.5}

In ([3.d.4]) and ([3.d.5]) the :math:`\hat{(\hskip 10pt
)}` denotes the provisional value before adjustment. The form
([3.d.5]) forces the arbitrary corrections to be small when the mixing
ratio is small and when the change made to the mixing ratio by the
advection is small. In addition, the :math:`\eta` factor is included to
make the changes approximately proportional to mass per unit volume .
Satisfying ([3.d.1]) and ([3.d.2]) gives

.. math::

   \alpha = \frac{\int\limits_3 \, q^- \Delta p^- - \int\limits_3 \, \hat
   q^+ p_0\Delta A - M\int\limits_3\hat q^+ \hat\pi^+ \Delta B}
   {\int\limits_3 \eta \hat q^+ \vert \hat q^+ - q^- \vert \, p_0\Delta A
   + M \int\limits_3\eta\hat q^+ \vert \hat q^+ - q^- \vert \hat\pi^+
   \Delta B} \label{3.d.6}

 and

.. math::

   {\mathbf{M}} = ({\mathbf{P}} + \int\limits_3 \, q^- \Delta p^-
   ) \Bigg/ \int\limits_2 \, \hat \pi^+ \ . \label{3.d.7}

Note that water vapor and dry mass are corrected simultaneously.
Additional advected constituents are treated as mixing ratios normalized
by the mass of dry air. This choice was made so that as the water vapor
of a parcel changed, the constituent mixing ratios would not change.
Thus the fixers which ensure conservation involve the dry mass of the
atmosphere rather than the moist mass as in the case of the specific
humidity above. Let :math:`\chi` denote the mixing ratio of
constituents. Historically we have used the following relationship for
conservation:

.. math::

   \int\limits_3\chi^+ (1-q^+) \Delta p^+ = \int\limits_3 \chi^-
   (1-q^-)\Delta p^- \ .
   \label{3.d.8}

The term :math:`(1-q)\Delta p` defines the dry air mass in a layer.
Following the change made by the fixer has the same form as ([3.d.5])

.. math::

   \chi^+ ( \lambda, \varphi, \eta ) = \hat \chi^+ +
   \alpha_\chi \eta \hat \chi^+ \vert \hat \chi^+ - \chi^- \vert
   \label{3.d.9} \ .

Substituting ([3.d.9]) into ([3.d.8]) and using ([3.d.4]) through
([3.d.7]) gives

.. math::

   \alpha_\chi = \frac{ \int\limits_3 \chi^- (1-q^-) \Delta p^- -
           \int\limits_{A,B}\hat\chi^+ (1-\hat q^+)\Delta\hat p^+
           + \alpha\int\limits_{A,B}\hat\chi^+ \eta \hat q^+
               \vert\hat q^+ - q^-\vert\Delta p}
           {\int\limits_{A,B} \eta\hat\chi^+ \vert\hat\chi^+ -
               \chi^-\vert(1-\hat q^+)\Delta p
           - \alpha \int\limits_{A,B}\eta\hat\chi^+
           \vert\hat\chi^+ - \chi^-\vert \eta\hat q^+ \vert \hat
           q^+ - q^-\vert\Delta p}
   \label{3.d.10} \ ,

where the following shorthand notation is adopted:

.. math::

   \int\limits_{A,B} (~~)\Delta p = \int\limits_3 (~~) p_0 \Delta A +
   M\int\limits_3 (~~) p_s \Delta B
   \label{3.d.11} \ .

We note that there is a small error in ([3.d.8]). Consider a situation
in which moisture is transported by a physical parameterization, but
there is no source or sink of moisture. Under this circumstance
:math:`q^- \ne q^0`, but the surface pressure is not allowed to change.
Since :math:`(1- q^-)\Delta p^- \ne
(1-q^0)\Delta p^0`, there is an implied change of dry mass of dry air in
the layer, and even in circumstances where there is no change of dry
mixing ratio :math:`\chi` there would be an implied change in mass of
the tracer. The solution to this inconsistency is to define a dry air
mass *only once* within the model time step, and use it consistently
throughout the model. In this revision, we have chosen to fix the dry
air mass in the model time step where the surface pressure is updated,
e.g. at the end of the model time step. Therefore, we now replace
([3.d.8]) with

.. math::

   \int\limits_3\chi^+ (1-q^+) \Delta p^+ = \int\limits_3 \chi^-
   (1-q^0)\Delta p^0
   \label{3.d.8a} \ .

There is a corresponding change in the first term of the numerator of
([3.d.10]) in which :math:`q^-` is replace by :math:`q^0`. uses
([3.d.10]) for water substances and constituents affecting the
temperature field to prevent changes to the IPCC simulations. In the
future, constituent fields may use a *corrected* version of ([3.d.10]).

Energy Fixer
~~~~~~~~~~~~

Following notation in section [massfixers], the total energy integrals
are

.. math::

   & \int\limits_3 {1 \over g} [ c_p T^+ + \Phi_s + {1 \over 2}
     ( {u^+}^2 + {v^+}^2 ) ] \Delta p^+
   = {\mathbf{E}} \\
    {\mathbf{E}} = \int\limits_3 {1 \over g} [ c_p T^- + \Phi_s +
        {1 \over 2} ( {u^-}^2 + {v^-}^2 ) ] \Delta p^-
   + {\mathbf{S}}

.. math::

   {\mathbf{S}} = \int\limits_2 [ ( FSNT - FLNT )
     -        ( FSNS - FLNS - SHFLX - \rho_{H_2O} L_v PRECT )
     -        ] \Delta t

where :math:`{\mathbf{S}}` is the net source of energy from the
parameterizations. :math:`FSNT` is the net downward solar flux at the
model top, :math:`FLNT` is the net upward longwave flux at the model
top, :math:`FSNS` is the net downward solar flux at the surface,
:math:`FLNS` is the net upward longwave flux at the surface,
:math:`SHFLX` is the surface sensible heat flux, and :math:`PRECT` is
the total precipitation during the time step. From equation ([3.d.4])

.. math::

   \pi^+ ( \lambda, \varphi ) = {\mathbf{M}} \hat \pi^+
   ( \lambda, \varphi )

 and from ([3.d.3])

.. math:: \Delta p = p_0 \Delta A + \pi \Delta B

 The energy fixer is chosen to have the form

.. math::

   T^+ ( \lambda, \varphi, \eta ) & = \hat T^+ + \beta \\ 
   u^+ ( \lambda, \varphi, \eta ) & = \hat u^+ \\ 
   v^+ ( \lambda, \varphi, \eta ) & = \hat v^+

Then

.. math::

   \beta = { g{\mathbf{E}}
    - \int\limits_3 \,
    [ c_p \hat T^+ + \Phi_s + {1 \over 2} ( {\hat u}^{+^2} +
     {\hat v}^{+^2} ) ]
    p_0\Delta A - {\mathbf{M}}\int\limits_3
    [ c_p \hat T^+ + \Phi_s + {1 \over 2} ( {\hat u}^{+^2} +
     {\hat v}^{+^2} ) ]
    \hat\pi^+ \Delta B
   \over \int\limits_3 c_p \, p_0\Delta A + {\mathbf{M}} \int\limits_3 c_p \hat\pi^+ \Delta B }

Statistics Calculations
~~~~~~~~~~~~~~~~~~~~~~~

At each time step, selected global average statistics are computed for
diagnostic purposes when the model is integrated with the Eulerian and
semi-Lagrangian dynamical cores. Let :math:`\int_3` denote a global and
vertical average and :math:`\int_2` a horizontal global average. For an
arbitrary variable :math:`\psi`, these are defined by

.. math::

   \int_3 \psi d V & = \sum^K_{k=1} \sum^J_{j=1} \sum^I_{i=1} 
   \psi_{ijk} w_j (\frac{\Delta p_k}{\pi} ) \bigg/ 2I , \label{8.a.1}\\*[-1.0em]
   \intertext{and}\nonumber\\*[-2.0em]
   \int_2 \psi dA & = \sum^J_{j=1} \sum^I_{i=1} \psi_{ijk} w_j/2I , \label{8.a.2}

where recall that

.. math:: \sum^J_{j=1} w_j = 2 . \label{8.a.3}

The quantities monitored are:

.. math::

   \text{global rms} \; (\zeta+f) (\mathrm{s}^{-1}) & = [ \int_3
   (\zeta^n + f)^2 dV ]^{1/2} , \label{8.a.4} \\
   \text{global rms} \; \delta (\mathrm{s}^{-1}) & =  [\int_3 (\delta^n)^2 dV
   ]^{1/2} , \label{8.a.5} \\
   \text{global rms} \; T \; (\mathrm{K}) \; & = [ \int_3 (T^r + T'^n)^2 dV
   ]^{1/2} , \label{8.a.6} \\
   \text{global average mass times} \; g \; (\mathrm{Pa}) & = \int_2
   \pi^{n} dA , \label{8.a.7} \\
   \text{global average mass of moisture} \; (\mathrm{kg \; m}^{-2}) 
   & = \int_3 \pi^{n} q^{n}/g dV . \label{8.a.8}

Reduced grid
~~~~~~~~~~~~

The Eulerian core and semi-Lagrangian tracer transport can be run on
reduced grids. The term reduced grid generally refers to a grid based on
latitude and longitude circles in which the longitudinal grid increment
increases at latitudes approaching the poles so that the longitudinal
distance between grid points is reasonably constant. Details are
provided in . This option provides a saving of computer time of up to
25%.

Semi-Lagrangian Dynamical Core
------------------------------

Introduction
~~~~~~~~~~~~

The two-time-level semi-implicit semi-Lagrangian spectral transform
dynamical core in evolved from the three-time-level CCM2 semi-Lagrangian
version detailed in hereafter referred to as W&O94. As a first
approximation, to convert from a three-time-level scheme to a
two-time-level scheme, the time level index n-1 becomes n, the time
level index n becomes n+\ :math:`\frac{1}{ 2}`, and :math:`2\Delta t`
becomes :math:`\Delta t`. Terms needed at n+\ :math:`\frac{1}{ 2}` are
extrapolated in time using time n and n-1 terms, except the Coriolis
term which is implicit as the average of time n and n+1. This leads to a
more complex semi-implicit equation to solve. Additional changes have
been made in the scheme to incorporate advances in semi-Lagrangian
methods developed since W&O94. In the following, reference is made to
changes from the scheme developed in W&O94. The reader is referred to
that paper for additional details of the derivation of basic aspects of
the semi-Lagrangian approximations. Only the details of the
two-time-level approximations are provided here.

Vertical coordinate and hydrostatic equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The semi-Lagrangian dynamical core adopts the same hybrid vertical
coordinate (:math:`\eta`) as the Eulerian core defined by

.. math:: p(\eta, p_s) = A (\eta)p_o + B(\eta) p_s \; , \label{sld.1}

where :math:`p` is pressure, :math:`p_s` is surface pressure, and
:math:`p_o` is a specified constant reference pressure. The coefficients
:math:`A` and :math:`B` specify the actual coordinate used. As mentioned
by and implemented by and , the coefficients :math:`A` and :math:`B` are
defined only at the discrete model levels. This has implications in the
continuity equation development which follows.

In the :math:`\eta` system the hydrostatic equation is approximated in a
general way by

.. math:: \Phi_k = \Phi_s + R \sum^K_{l=k} H_{kl}\, (p)\, T_{vl} \label{sld.2}

where k is the vertical grid index running from 1 at the top of the
model to :math:`K` at the first model level above the surface,
:math:`\Phi_k` is the geopotential at level :math:`k`, :math:`\Phi_s` is
the surface geopotential, :math:`T_v` is the virtual temperature, and R
is the gas constant. The matrix :math:`H`, referred to as the
hydrostatic matrix, represents the discrete approximation to the
hydrostatic integral and is left unspecified for now. It depends on
pressure, which varies from horizontal point to point.

Semi-implicit reference state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The semi-implicit equations are linearized about a reference state with
constant :math:`T^r` and :math:`p_s^r`. We choose

.. math:: T^r = 350 {\rm K},~~~~~ p_s^r = 10^5 {\rm Pa}

Perturbation surface pressure prognostic variable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ameliorate the mountain resonance problem, introduce a perturbation
:math:`{\rm ln}\,p_{s}` surface pressure prognostic variable

.. math::

   {\rm ln}\,p'_{s} & = {\rm ln}\,p_{s} - {\rm ln}\,p_{s}^* \\ {\rm
   ln}\,p_{s}^* & = - \frac{\Phi_s}{R T^r}

The perturbation surface pressure, :math:`{\rm ln}\,p'_{s}`, is never
actually used as a grid point variable in the code. It is only used for
the semi-implicit development and solution. The total :math:`{\rm
ln}\,p_{s}` is reclaimed in spectral space from the spectral
coefficients of :math:`\Phi_s` immediately after the semi-implicit
equations are solved, and transformed back to spectral space along with
its derivatives. This is in part because
:math:`\nabla ^4{\rm ln}\,p_{s}` is needed for the horizontal diffusion
correction to pressure surfaces. However the semi-Lagrangian default is
to run with no horizontal diffusion.

Extrapolated variables
~~~~~~~~~~~~~~~~~~~~~~

Variables needed at time (:math:`n+\frac{1}{2}`) are obtained by
extrapolation

.. math::

   (~~~)^{n+\frac{1}{2}} = \frac{3}{2} (~~~)^{n}
                      - \frac{1}{2} (~~~)^{n-1}

Interpolants
~~~~~~~~~~~~

Lagrangian polynomial quasi-cubic interpolation is used in the
prognostic equations for the dynamical core. Monotonic Hermite
quasi-cubic interpolation is used for tracers. Details are provided in
the Eulerian Dynamical Core description. The trajectory calculation uses
tri-linear interpolation of the wind field.

Continuity Equation
~~~~~~~~~~~~~~~~~~~

The discrete semi-Lagrangian, semi-implicit continuity equation is
obtained from (16) of W&O94 modified to be spatially uncentered by a
fraction :math:`\epsilon`, and to predict :math:`{\rm ln}\,p_{s}'`

.. math::

   \Delta B_{_{l}} & \{ ({\rm
             ln}\,p_{s{_{l}}}')^{n+1}_{_{A}}
     - [ ({\rm ln}\,p_{s{_{l}}} )^{n}
            + { \Phi_s \over RT^r } ] _{D_{2}} \} \, \Bigg/
                    \, \Delta t = \nonumber \\
       & - {1\over 2} \bigg\lbrace [ (1 + \epsilon )
                    \Delta ({1 \over p_s}
       \dot\eta {\partial p \over
         \partial\eta})_l]^{n+1}_A + [ (1 -
         \epsilon )\Delta ({1 \over p_s}
       \dot\eta {\partial p \over
       \partial\eta})_l]^{n}_{D_{2}}\bigg\rbrace \\
   &-({1 \over p_s} \delta_{_{l}} \Delta p_{_{l}} )^{n+{1\over
       2}}_{M_{2}} 
    +{\Delta B_{_{l}} \over RT^r} ( {\mathbf{V}}_{_{l}} \cdot
         \nabla \;\Phi_s )^{n+{1\over 2}}_{M_2} \nonumber\\
   &-\bigg\lbrace{1 \over 2}
         [ (1 + \epsilon ) ({1 \over p^r_s}
                \delta_{_{l}} \Delta p_{_{l}}^r
       )^{n+1}_A
       + (1 - \epsilon ) ({1 \over p^r_s} \delta_{_{l}}
                \Delta p_{_{l}}^r
       )^{n}_{D_{2}}]
       - ({1 \over p^r_s} \delta_{_{l}} \Delta p_{_{l}}^r
       )^{n+{1\over 2}}_{M_{2}}\bigg\rbrace \nonumber

where

.. math::

   \Delta (~~~~)_l \,& =\, (~~~~)_{l+ {1 \over 2}} -\, (~~~~)_{l - {1
   \over 2}}\\[-1.0em]
   \intertext{and}\nonumber\\[-2.0em]
    (~~~~ )^{n+{1\over 2}}_{M_{2}} & =
        {1\over 2} [ (1 + \epsilon )(~~~~
                          )^{n+{1\over 2}}_{A}
                      + (1 - \epsilon )(~~~~
                        )^{n+{1\over 2}}_{D_{2}} ]
   \label{sld.3}

:math:`\Delta (~~~~)_l` denotes a vertical difference, :math:`l`
denotes the vertical level, :math:`A` denotes the arrival point,
:math:`D_2` the departure point from horizontal (two-dimensional)
advection, and :math:`M_2` the midpoint of that trajectory.

The surface pressure forecast equation is obtained by summing over all
levels and is related to (18) of W&O94 but is spatially uncentered and
uses :math:`{\rm ln}\,p_{s}'`

.. math::

    ({\rm ln}\,p'_s )^{n+1}_{_{A}} & = \sum^K_{l=1} \Delta
       B_l
       [ ({\rm ln}\,p_{s{_{l}}} )^{n} + { \Phi_s
          \over RT^r } ] _{D_{2}}
        - {1 \over 2} \Delta t \sum^K_{l=1}[ (1 - \epsilon
               )\Delta ({1 \over p_s}
       \dot\eta {\partial p \over
       \partial\eta})_l]^{n}_{D_{2}} \nonumber \\
   & \phantom{=} - \Delta t \sum^K_{l=1}({1 \over p_s} \delta_l
       \Delta p_{_{l}} )^{n+{1\over 2}}_{M_{2}}
   + \Delta t \sum^K_{l=1} {\Delta B_{_{l}} \over RT^r} (
         {\mathbf{V}}_{_{l}} \cdot \nabla \;\Phi_s )^{n+{1\over
         2}}_{M_2} \\
       & \phantom{=}- \Delta t \sum^K_{l=1} {1 \over p^r_s}
       \bigg\lbrace{1 \over 2}
       [ (1 + \epsilon )
                   (\delta_l)^{n+1}_{_{A}}
               + (1 - \epsilon )
                   (\delta_l)^{n}_{_{D_{2}}}
       ]
                  - (\delta_l)^{n+{1\over 2}}_{_{M_{2}}}
             \bigg\rbrace \Delta p^r_{_{l}} \nonumber

The corresponding :math:`({1 \over p_s} \dot\eta {\partial p \over
\partial \eta})` equation for the semi-implicit development
follows and is related to (19) of W&O94, again spatially uncentered and
using :math:`{ln}\,p_{s}'`.

.. math::

   (1 + \epsilon ) ({1 \over p_s} \dot\eta {\partial p
   \over \partial
   \eta})^{n+1}_{k+{1 \over 2}} = - {2 \over \Delta t}
       \bigg\lbrace B_{k+ {1\over 2}} ({\rm ln}\; p'_s
       )^{n+1}_{_A} -
       \sum^k_{l=1} \Delta B_l [ ({\rm ln}\; p_{s_l}
              )^{n} + { \Phi_s \over RT^r } ] _{D_2}
              \bigg\rbrace \nonumber \\
            &- \sum^k_{l=1}[ (1 - \epsilon )\Delta
               ({1 \over p_s}
           \dot\eta {\partial p \over \partial\eta})_l
                 ]^{n}_{D_{2}}\\
   &- 2 \sum^k_{l=1} ({1 \over p_s} \delta_l \Delta p_l
       )^{n+{1\over 2}}_{M_2}
   + 2 \sum^k_{l=1}{\Delta B_{_{l}} \over RT^r} (
         {\mathbf{V}}_{_{l}} \cdot \nabla \;\Phi_s )^{n+{1\over
         2}}_{M_2} \nonumber \\
       &- 2 \sum^k_{l=1} {1 \over p^r_s} \bigg\lbrace{1 \over 2}
           [ (1 + \epsilon )
                   (\delta_l)^{n+1}_{_{A}}
               + (1 - \epsilon )
                   (\delta_l)^{n}_{_{D_{2}}}
           ]
                   - (\delta_l)^{n+{1\over 2}}_{_{M_{2}}}
             \bigg\rbrace \Delta p^r_{_{l}} \nonumber

This is not the actual equation used to determine
:math:`({1 \over p_s}
\dot\eta {\partial p \over \partial \eta})` in the code. The
equation actually used in the code to calculate
:math:`({1 \over p_s}
\dot\eta {\partial p \over \partial \eta})` involves only the
divergence at time (:math:`n+1`) with
:math:`({\rm ln}\; p'_s )^{n+1}` eliminated.

.. math::

   (1 + \epsilon ) &({1 \over p_s} \dot\eta {\partial p
   \over \partial
   \eta})^{n+1}_{k+{1 \over 2}} = \nonumber \\
     {2 \over \Delta t}
       &[ \sum^k_{l=1} ~-~ B_{k+{1\over 2}}\sum^K_{l=1}]
              \Delta B_l [ ({\rm ln}\; p_{s_l} )^{n} + {
              \Phi_s \over RT^r } ] _{D_2} \nonumber \\
    - &[ \sum^k_{l=1} ~-~ B_{k+{1\over 2}}\sum^K_{l=1}]
           [ (1 - \epsilon ) \Delta ({1 \over p_s}
              \dot\eta {\partial p \over \partial
       \eta})_l]^{n}_{D_2} \nonumber \\
    - 2 &[ \sum^k_{l=1} ~-~ B_{k+{1\over 2}}\sum^K_{l=1}]
           ({1
       \over p_s} \delta_l \Delta p_l )^{n+{1\over 2}}_{M_2} \\
    + 2 &[ \sum^k_{l=1} ~-~ B_{k+{1\over 2}}\sum^K_{l=1}]
         {\Delta B_{_{l}} \over RT^r} ( {\mathbf{V}}_{_{l}} \cdot
         \nabla \;\Phi_s )^{n+{1\over 2}}_{M_2} \nonumber \\
    - 2 &[ \sum^k_{l=1} ~-~ B_{k+{1\over 2}}\sum^K_{l=1}]
          {1 \over p_s^r}\bigg\lbrace {1 \over 2}
             [ (1 + \epsilon )
                (\delta_l)^{n+1}_{_A}
              + (1 - \epsilon ) (\delta_l
             )^{n}_{_{D_2}} ]
                                    - (\delta_l )^{n+{1\over
       2}}_{_{M_2}} \bigg\rbrace \Delta p^r_l \nonumber

The combination
:math:`[ ({\rm ln}\,p_{s{_{l}}} )^{n} + { \Phi_s \over RT^r} + {1 \over 2} {\Delta t \over RT^r} ( {\mathbf{V}} \cdot \nabla  \;\Phi_s )^{n+{1\over 2}}  ] _{D_{2}}` is treated as a unit, and follows from ([sld.3]).

Thermodynamic Equation
~~~~~~~~~~~~~~~~~~~~~~

The thermodynamic equation is obtained from (25) of W&O94 modified to be
spatially uncentered and to use :math:`{\rm ln}\,p_{s}'`. In addition
Hortal’s modification is included, in which

.. math::

   {d \over dt} [
   -(p_s B {\partial T \over \partial p})_{ref}
   { \Phi_s \over RT^r } ]

is subtracted from both sides of the temperature equation. This is akin
to horizontal diffusion which includes the first order term converting
horizontal derivatives from eta to pressure coordinates, with
:math:`({\rm ln}\;p_s)` replaced by
:math:`-{ \Phi_s \over RT^r}`, and :math:`(p_s B {\partial T \over \partial p})_{ref}`
taken as a global average so it is invariant with time and can commute
with the differential operators.

.. math::

   {T_A^{n+1} - T_D^{n} \over \Delta t} & =
   \{ \{[
   -(p_s B(\eta) {\partial T \over \partial p})_{ref}
   { \Phi_s \over RT^r } ]_A^{n+1} - [
   -(p_s B(\eta) {\partial T \over \partial p})_{ref}
   { \Phi_s \over RT^r } ]_D^{n} \} \Bigg/ \Delta t .
   \nonumber \\
   & \phantom{=} . +{ 1 \over RT^r }[
   (p_s B(\eta) {\partial T \over \partial p})_{ref}
    {\mathbf{V}} \cdot \nabla \;\Phi_s
   + \Phi_s \dot\eta {\partial \over \partial \eta} (p_s B(\eta)
   {\partial T \over \partial p})_{ref} ]_M^{n+{1\over 2}}
   \}
     \nonumber \\
   & \phantom{=}+({RT_v \over c_p^*} {\omega \over p}
       )^{n+{1\over 2}}_M + Q^{n}_M \nonumber \\
   & \phantom{=}+ {RT^r \over c_p} {p_s^r \over p^r}[ B(\eta)
       {d_2\;{\rm ln}\;p_s' \over dt} + \overline{({1 \over p_s}
       \dot \eta {\partial p \over \partial \eta} )}^t ]
       \\
   & \phantom{=}- {RT^r \over c_p} {p_s^r \over p^r} [({p \over
           p_s}) ( {\omega \over p} ) ]^{n+{1\over
           2}}_M \nonumber \\
   & \phantom{=}- {RT^r \over c_p} {p_s^r \over p^r} B(\eta) [ {1
             \over RT^r} {\mathbf{V}} \cdot \nabla \;\Phi_s
             ]^{n+{1\over 2}}_{M_2} \nonumber

Note that :math:`Q^{n}` represents the heating calculated to advance
from time :math:`n` to time :math:`n+1` and is valid over the interval.

The calculation of :math:`(p_s B {\partial T \over \partial
p})_{ref}` follows that of the ECMWF (Research Manual 3, ECMWF
Forecast Model, Adiabatic Part, ECMWF Research Department, 2nd edition,
1/88, pp 2.25-2.26) Consider a constant lapse rate atmosphere

.. math::

   T & = T_0({p\over p_0})^{R \gamma / g} \\ {\partial T \over
   \partial p} & = {1 \over p}{R \gamma \over g} T_0({p\over
   p_0})^{R \gamma / g} \\ p_s B {\partial T \over \partial p} & =
   B {p_s \over p}{R \gamma \over g} T \\
   ( p_s B {\partial T \over \partial p} )_{ref} & = B_k
    {(p_s)_{ref} \over (p_k)_{ref}}
   {R \gamma \over g} (T_k)_{ref} ~~\hbox{for}~~ (T_k)_{ref} > T_C \\
   ( p_s B {\partial T \over \partial p} )_{ref} & = 0 ~~for~~
    (T_k)_{ref} \leq T_C \\
   (p_k)_{ref} & = A_k p_0 + B_k (p_s)_{ref} \\ (T_k)_{ref} & =
   T_0({(p_k)_{ref} \over (p_s)_{ref} })^{R \gamma / g} \\
   (p_s)_{ref} & = 1013.25 {\rm mb} \\ T_0 & = 288 {\rm K} \\ p_0 & =
   1000 {\rm mb} \\ \gamma & = 6.5 {\rm K / km} \\ T_C & = 216.5 {\rm K}

Momentum equations
~~~~~~~~~~~~~~~~~~

The momentum equations follow from (3) of W&O94 modified to be spatially
uncentered, to use :math:`{\rm ln}\,p_{s}'`, and with the Coriolis term
implicit following and . The semi-implicit, semi-Lagrangian momentum
equation at level :math:`k` (but with the level subscript :math:`k`
suppressed) is

.. math::

   {{\mathbf{V}}^{n+1}_{_{A}} - {\mathbf{V}}^{n}_{_{D}} \over \Delta t} & =  
           -{1 \over 2}\bigg\lbrace \,
              (1 + \epsilon ) [f{\mathbf{{\hat{k}}}}
                  \times {\mathbf{V}}
           ]^{n+1}_{_{A}}
           + (1 - \epsilon ) [f{\mathbf{{\hat{k}}}}
                  \times {\mathbf{V}}
           ]^{n}_{_{D}}\bigg\rbrace + {\mathbf{F}}^n_M \nonumber
           \\
           & \phantom{=} -{1 \over 2}\bigg\lbrace\, (1 + \epsilon
             )
            [ \nabla (\Phi_s + R {\mathbf{H}}_k \cdot
           {\mathbf{T}}_v )
           + RT_v {B \over p} p_s \nabla\,{\rm ln}\,p_s
             ]^{n+{1\over 2}}_{_{A}} \nonumber \\
           & \phantom{=} \phantom{-{1 \over 2}} + (1 - \epsilon
              )
            [ \nabla (\Phi_s + R {\mathbf{H}}_k \cdot
           {\mathbf{T}}_v )
           + RT_v {B \over p} p_s \nabla\,{\rm ln}\,p_s
             ]^{n+{1\over 2}}_{_{D}}\bigg\rbrace
            \nonumber \\ & \phantom{=} - {1 \over 2} \bigg\lbrace\,
          (1 + \epsilon ) \nabla [ R {\mathbf{H}}_k^r
                  \cdot {\mathbf{T}}
              + RT^r\,{\rm ln}\,p_s' ]^{n+1}_{_{A}} \\
          & \phantom{=} \phantom{- {1 \over 2}} - (1 + \epsilon
             )
             \nabla [\Phi_s + R {\mathbf{H}}_k^r \cdot
                  {\mathbf{T}} + RT^r\,{\rm ln}\,p_s]^{n+{1\over
                  2}}_{_{A}} \nonumber \\
           & \phantom{=} \phantom{- {1 \over 2}} + (1 - \epsilon
             )
             \nabla[\Phi_s + R {\mathbf{H}}_k^r \cdot {\mathbf{T}}
                  + RT^r\,{\rm ln}\,p_s]^{n}_{_{D}}\nonumber \\
          & \phantom{=} \phantom{- {1 \over 2}}- (1 - \epsilon
             )
             \nabla[\Phi_s + R {\mathbf{H}}_k^r \cdot {\mathbf{T}}
                  + RT^r\,{\rm ln}\,p_s]^{n+{1\over
                  2}}_{_{D}}\bigg\rbrace \nonumber

The gradient of the geopotential is more complex than in the
:math:`\sigma` system because the hydrostatic matrix
:math:`\mathbf{H}` depends on the local pressure:

.. math::

   \nabla({\mathbf{H}}_k \cdot {\mathbf{T}}_v ) =
   {\mathbf{H}}_k \cdot [(1 + \epsilon_v
        {\mathbf{q}})\nabla {\mathbf{T}} +
        \epsilon_v {\mathbf{T}}\nabla {\mathbf{q}} ] +
            {\mathbf{T}}_v \cdot \nabla {\mathbf{H}}_k

where :math:`\epsilon_v` is :math:`(R_v/R - 1)` and :math:`R_v` is the
gas constant for water vapor. The gradient of :math:`T` is calculated
from the spectral representation and that of :math:`q` from a discrete
cubic approximation that is consistent with the interpolation used in
the semi-Lagrangian water vapor advection. In general, the elements of
:math:`\mathbf{H}` are functions of pressure at adjacent discrete
model levels

.. math:: H_{kl} = f_{kl}(p_{l+1/2},p_l,p_{l-1/2})

The gradient is then a function of pressure and the pressure gradient

.. math::

   \nabla H_{kl} = g_{kl}(p_{_{l+1/2}},\; p_{_{l}},\; p_{_{l-1/2}},
         \nabla p_{_{l+1/2}}, \nabla p_{_{l}}, \nabla p_{_{l-1/2}})

The pressure gradient is available from ([sld.1]) and the surface
pressure gradient calculated from the spectral representation

.. math:: \nabla p_{_{l}} = B_l\nabla p_s = B_l p_s \nabla\,{\rm ln}\,p_s

Development of semi-implicit system equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The momentum equation can be written as

.. math::

   {{\mathbf{V}}^{n+1}_{_{A}} - {\mathbf{V}}^{n}_{_{D}} \over \Delta t}
   =
           -{1 \over 2}\bigg\lbrace
             & (1 + \epsilon ) [f{\mathbf{{\hat{k}}}}
                  \times {\mathbf{V}}
           ]^{n+1}_{_{A}}
           + (1 - \epsilon ) [f{\mathbf{{\hat{k}}}}
                  \times {\mathbf{V}}
           ]^{n}_{_{D}}\bigg\rbrace \nonumber \\
            - {1 \over 2} \bigg\lbrace
          & (1 + \epsilon ) \nabla [ R {\mathbf{H}}_k^r
                  \cdot {\mathbf{T}}
              + RT^r\,{\rm ln}\,p_s' ]^{n+1}_{_{A}}\bigg\rbrace
          + RHS_{\mathbf{V}}
   \label{vmeq}
   ~,

where :math:`RHS_{\mathbf{V}}` contains known terms at times
(:math:`{n+{1\over 2}}`) and (:math:`n`).

By combining terms, [vmeq] can be written in general as

.. math::

   {\cal U}^{^{n+1}}_{_{A}} {\bf {\hat{i}}}_{_{A}} + {\cal
     V}^{^{n+1}}_{_{A}} {\bf {\hat{j}}}_{_{A}} = {\cal U}_{_{A}} {\bf
     {\hat{i}}}_{_{A}} + {\cal V}_{_{A}} {\bf {\hat{j}}}_{_{A}} + {\cal U}_{_{D}}
     {\bf {\hat{i}}}_{_{D}} + {\cal V}_{_{D}} {\bf {\hat{j}}}_{_{D}}
   \label{dw8}
   ~,

where :math:`\bf{{\hat{i}}}` and :math:`\bf{{\hat{j}}}` denote the
spherical unit vectors in the longitudinal and latitudinal directions,
respectively, at the points indicated by the subscripts, and
:math:`\cal U` and :math:`\cal V` denote the appropriate combinations of
terms in [vmeq]. Note that :math:`{\cal U}^{^{n+1}}_{_{A}}` is distinct
from the :math:`{\cal U}_{_{A}}`. Following , equations for the
individual components are obtained by relating the unit vectors at the
departure points
(:math:`{\bf {\hat{i}}}_{_{D}}`,\ :math:`{\bf {\hat{j}}}_{_{D}}`) to
those at the arrival points
(:math:`{\bf {\hat{i}}}_{_{A}}`,\ :math:`{\bf {\hat{j}}}_{_{A}}`):

.. math::

     {\bf {\hat{i}}}_{_{D}} = \alpha^u_{_{A}}{\bf {\hat{i}}}_{_{A}} +
     \beta^u_{_{A}}{\bf {\hat{j}}}_{_{A}} \\ {\bf {\hat{j}}}_{_{D}} =
     \alpha^v_{_{A}}{\bf {\hat{i}}}_{_{A}} + \beta^v_{_{A}}{\bf {\hat{j}}}_{_{A}}
   ~,

in which the vertical components (:math:`{\bf {\hat{k}}}`) are ignored.
The dependence of :math:`\alpha`\ ’s and :math:`\beta`\ ’s on the
latitudes and longitudes of the arrival and departure points is given in
the Appendix of .

W&O94 followed which ignored rotating the vector to remain parallel to
the earth’s surface during translation. We include that factor by
keeping the length of the vector written in terms of
:math:`({\mathbf{{\hat{i}}}}_{_{A}},{\mathbf{{\hat{j}}}}_{_{A}})`
the same as the length of the vector written in terms of
:math:`({\mathbf{{\hat{i}}}}_{_{D}},{\mathbf{{\hat{j}}}}_{_{D}})`.
Thus, (10) of W&O94 becomes

.. math::

     {\cal U}^{^{n+1}}_{_{A}} & = {\cal U}_{_{A}} +
     \gamma\alpha^u_{_{A}}{\cal U}_{_{D}} + \gamma\alpha^v_{_{A}}{\cal
     V}_{_{D}} \nonumber \\ {\cal V}^{^{n+1}}_{_{A}} & = {\cal V}_{_{A}} +
     \gamma\beta^u_{_{A}}{\cal U}_{_{D}} + \gamma\beta^v_{_{A}}{\cal
     V}_{_{D}}

where

.. math::

   \gamma = [
       { {\cal U}_{_{D}}^2 + {\cal V}_{_{D}}^2 \over ({\cal
         U}_{_{D}}\alpha^u_{_{A}} + {\cal
         V}_{_{D}}\alpha^v_{_{A}})^2 + ({\cal
         U}_{_{D}}\beta^u_{_{A}} + {\cal V}_{_{D}}\beta^v_{_{A}})^2
         }
       ]^{1 \over 2}

 After the momentum equation is written in a common set of unit vectors

.. math::

   {\mathbf{V}}^{n+1}_{_{A}} + ({1 + \epsilon \over 2})
     \Delta t
     [f{\mathbf{{\hat{k}}}} \times {\mathbf{V}}
       ]^{n+1}_{_{A}}
     + ({1 + \epsilon \over 2}) \Delta t
     \nabla [ R {\mathbf{H}}_k^r \cdot {\mathbf{T}} + RT^r\,{\rm
       ln}\,p_s' ]^{n+1}_{_{A}}
     = {{\cal R}}_{\mathbf{V}}^*

Drop the :math:`(~~~)^{n+1}_A` from the notation, define

.. math:: \alpha = (1 + \epsilon) \Delta t \Omega

and transform to vorticity and divergence

.. math::

     \zeta + \alpha \sin \varphi \delta + {\alpha \over a} v \cos\varphi
       & = {1 \over a\cos\varphi} [ {\partial{{\cal R}}_v^* \over
       \partial \lambda}
         - {\partial \over \partial\varphi}({{\cal R}}_u^*
       \cos\varphi)] \\ \delta - \alpha \sin \varphi
       \zeta + {\alpha \over a} u \cos\varphi
       & +& ( {1 + \epsilon \over 2} ) \Delta t \nabla^2
           [ R {\mathbf{H}}_k^r \cdot {\mathbf{T}}
         + RT^r\,{\rm ln}\,p_s' ]^{n+1}_{_{A}} \nonumber \\ & =
       {1 \over a\cos\varphi}
       [ {\partial{{\cal R}}_u^* \over \partial \lambda} + {\partial
         \over \partial\varphi}({{\cal R}}_v^* \cos\varphi)]

Note that

.. math::

     u\cos\varphi & = {1 \over a} {{\partial}\over{\partial}\lambda}(\nabla^{-2}
     \delta ) - {\cos\varphi \over a}
     {{\partial}\over{\partial}\varphi}(\nabla^{-2} \zeta ) \\ v\cos\varphi
     & = {1 \over a} {{\partial}\over{\partial}\lambda}(\nabla^{-2} \zeta ) +
     {\cos\varphi \over a} {{\partial}\over{\partial}\varphi}(\nabla^{-2} \delta
     )

Then the vorticity and divergence equations become

.. math::

     \zeta + \alpha \sin \varphi \delta + {\alpha \over
     a^2}{{\partial}\over{\partial}\lambda}(\nabla^{-2} \zeta )
     &+&{\alpha\cos\varphi \over a^2}{{\partial}\over{\partial}\varphi}(\nabla^{-2}
     \delta ) \nonumber \\ & ={1 \over a\cos\varphi} [
     {{\partial}{{\cal R}}_v^* \over {\partial}\lambda}
       - {{\partial}\over {\partial}\varphi}({{\cal R}}_u^* \cos\varphi)] =
     {{\cal L}}\\ \delta - \alpha \sin \varphi \zeta + {\alpha \over
     a^2}{{\partial}\over{\partial}\lambda}(\nabla^{-2} \delta )
     &-&{\alpha\cos\varphi \over a^2}{{\partial}\over{\partial}\varphi}(\nabla^{-2}
     \zeta ) + ( {1 + \epsilon \over 2}) \Delta t
     \nabla^2 [ R {\mathbf{H}}_k^r \cdot {\mathbf{T}} + RT^r\,{\rm
       ln}\,p_s' ]^{n+1}_{_{A}} \nonumber \\
     & = {1 \over a\cos\varphi}
     [ {{\partial}{{\cal R}}_u^* \over {\partial}\lambda} + {{\partial}\over
       {\partial}\varphi}({{\cal R}}_v^* \cos\varphi)]
     = {{\cal M}}

Transform to spectral space as described in the description of the
Eulerian spectral transform dynamical core. Note, from (4.5b) and (4.6)
on page 177 of

.. math::

     \mu P_n^m & = D_{n+1}^m P_{n+1}^m + D_{n}^m P_{n-1}^m \\ D_{n}^m & =
     ({n^2-m^2 \over 4n^2-1})^{1 \over 2}

and from (4.5a) on page 177 of

.. math::

   ( 1 - \mu^2) {{\partial}\over{\partial}\mu}P_n^m = -nD_{n+1}^m P_{n+1}^m
     +  (n+1 ) D_{n}^m P_{n-1}^m

Then the equations for the spectral coefficients at time :math:`n+1` at
each vertical level are

.. math::

     \zeta_n^m ( 1 - {im\alpha \over n(n+1)}) +
     \delta_{n+1}^m \alpha ({n\over n+1}) D_{n+1}^m +
     \delta_{n-1}^m \alpha ({n+1\over n}) D_{n}^m & = {{\cal L}}_n^m
     \\ \delta_n^m ( 1 - {im\alpha \over n(n+1)})
     - \zeta_{n+1}^m \alpha ({n\over n+1}) D_{n+1}^m
     - \zeta_{n-1}^m \alpha ({n+1\over n}) D_{n}^m && \\
     - ({1 + \epsilon \over 2}) \Delta t
     {n(n+1) \over a^2} [ R {\mathbf{H}}_k^r \cdot
       {\mathbf{T}_n^m} + RT^r\,{\rm ln}\,{p_s'}_n^m ]
     & = {{\cal M}}_n^m \nonumber

.. math::

     {{\rm ln} p'_s}^m_n & = {\rm PS}^m_n - ({ 1 + \epsilon \over
     2}) {\Delta t\over p^r_s} (\underline{\Delta p^r}
     )^T \underline{\delta}^m_n \\ \underline{T}^m_n & =
     \underline{\rm TS}^m_n - ({ 1 + \epsilon \over 2}) \Delta
     t {\mathbf{D}}^r \underline{\delta}^m_n

The underbar denotes a vector over vertical levels. Rewrite the
vorticity and divergence equations in terms of vectors over vertical
levels.

.. math::

     \underline{\delta}_n^m ( 1 - {im\alpha \over
     n(n+1)})
     - \underline{\zeta}_{n+1}^m \alpha ({n\over n+1})
     - D_{n+1}^m \underline{\zeta}_{n-1}^m \alpha ({n+1\over
     - n}) D_{n}^m && \\
     - ({1 + \epsilon \over 2}) \Delta t
     {n(n+1) \over a^2} [ R {\mathbf{H}}^r
       \underline{T}_n^m + R\underline{T}^r\,{\rm ln}\,{p_s'}_n^m ]
     & = \underline{DS}_n^m \nonumber \\ \underline{\zeta}_n^m ( 1 -
     {im\alpha \over n(n+1)}) +
     \underline{\delta}_{n+1}^m \alpha ({n\over n+1})
     D_{n+1}^m + \underline{\delta}_{n-1}^m \alpha ({n+1\over
     n}) D_{n}^m & = \underline{VS}_n^m

Define :math:`\underline{h}_n^m` by

.. math::

     g\underline{h}_n^m & = R {\mathbf{H}}^r T_n^m +
              R\underline{T}^r\,{\rm ln}\,{p_s'}_n^m\\[-1.0em]
   \intertext{and}\nonumber\\[-2.0em] {{\cal A}}_n^m & = 1 - {im\alpha \over
   n(n+1)} \\ {{{\cal B}}^+}_n^m & = \alpha({n\over n+1})
   D_{n+1}^m \\ {{{\cal B}}^-}_n^m & = \alpha({n+1\over n}) D_{n}^m

Then the vorticity and divergence equations are

.. math::

   {{\cal A}}_n^m \underline{\zeta}_n^m + {{{\cal B}}^+}_n^m \underline{\delta}_{n+1}^m +
   {{{\cal B}}^-}_n^m \underline{\delta}_{n-1}^m & =
   \underline{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_n^m \\ {{\cal A}}_n^m \underline{\delta}_n^m
   - {{{\cal B}}^+}_n^m \underline{\zeta}_{n+1}^m {{{\cal B}}^-}_n^m
   - \underline{\zeta}_{n-1}^m
   -({1 + \epsilon \over 2}) \Delta t {n(n+1) \over
      a^2} g\underline{h}_n^m
   & = \underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}_n^m

Note that these equations are uncoupled in the vertical, i.e. each
vertical level involves variables at that level only. The equation for
:math:`\underline{h}_n^m` however couples all levels.

.. math::

   g\underline{h}_n^m = -({1 + \epsilon \over 2}) \Delta t
   [
   R {\mathbf{H}}^r {\mathbf{D}}^r + R\underline{T}^r
     {(\underline{\Delta p^r} )^T \over p^r_s} ]
   \underline{\delta}_n^m + R {\mathbf{H}}^r
   \underline{{\hbox{\sffamily\slshape T}}{\hbox{\sffamily\slshape S}}}_n^m + R\underline{T}^r {\rm PS}_n^m

Define :math:`{\mathbf{C}}^r` and
:math:`\underline{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}_n^m`
so that

.. math::

   g\underline{h}_n^m = -({1 + \epsilon \over 2}) \Delta t
   {{\hbox{\sffamily\slshape C}}}^r\underline{\delta}_n^m + \underline{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}_n^m

Let :math:`g{\rm D}_\ell` denote the eigenvalues of
:math:`{\mathbf{C}}^r` with corresponding eigenvectors
:math:`\underline{\Phi}_\ell` and :math:`{\mathbf{\Phi}}` is the
matrix with columns :math:`\underline{\Phi}_\ell`

.. math::

   {\mathbf{\Phi}}= ( \begin{array}{*{3}{c@{\:}}c} \underline{\Phi}_1 &
      \underline{\Phi}_2 & \dots & \underline{\Phi}_L \\
    \end{array})

and :math:`g{\mathbf{D}}` the diagonal matrix of corresponding
eigenvalues

.. math::

   g{\mathbf{D}} & = g (
    \begin{array}{*{3}{c@{\:}}c}
   {\rm D}_1 & 0 & \cdots & 0 \\ 0 & {\rm D}_2 & \cdots & 0 \\ \vdots &
   \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & {\rm D}_L \\
   \end{array}
   ) \\ {\mathbf{C}}^r {{\mathbf{\Phi}}} & = {{\mathbf{\Phi}}} g {\mathbf{D}} \\
   {{\mathbf{\Phi}}}^{-1}{\mathbf{C}}^r {{\mathbf{\Phi}}} & = g {\mathbf{D}}

Then transform

.. math::

   {2}
   \underline{\tilde\zeta}_n^m & = {\mathbf{\Phi}}^{-1} \underline{\zeta}_n^m
   & \quad,\quad \underline{\widetilde{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}}_n^m & =
   {\mathbf{\Phi}}^{-1} \underline{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_n^m\\
   \underline{\tilde\delta}_n^m & = {\mathbf{\Phi}}^{-1} \underline{\delta}_n^m
   & \quad,\quad \underline{\widetilde{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}}_n^m & =
   {\mathbf{\Phi}}^{-1} \underline{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}_n^m\\ \underline{\tilde
   h}_n^m & = {\mathbf{\Phi}}^{-1} \underline{ h}_n^m & \quad,\quad
   \underline{\widetilde{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}}_n^m & = {\mathbf{\Phi}}^{-1}
   \underline{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}_n^m

.. math::

   {{\cal A}}_n^m \underline{\tilde\zeta}_n^m + {{{\cal B}}^+}_n^m
   \underline{\tilde\delta}_{n+1}^m + {{{\cal B}}^-}_n^m
   \underline{\tilde\delta}_{n-1}^m & =
   \underline{\widetilde{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}}_n^m \\ {{\cal A}}_n^m
   \underline{\tilde\delta}_n^m
   - {{{\cal B}}^+}_n^m \underline{\tilde\zeta}_{n+1}^m {{{\cal B}}^-}_n^m
   - \underline{\tilde\zeta}_{n-1}^m
   -({1 + \epsilon \over 2}) \Delta t {n(n+1) \over
      a^2} g\underline{\tilde h}_n^m
   & = \underline{\widetilde{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}}_n^m \\ g\underline{\tilde
   h}_n^m + ({1 + \epsilon \over 2}) \Delta t
   {\mathbf{\Phi}}^{-1}{\mathbf{C}}^r{\mathbf{\Phi}}{\mathbf{\Phi}}^{-1}\underline{\delta}_n^m & =
   \underline{\widetilde{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}}_n^m \\ \underline{\tilde
   h}_n^m + ({1 + \epsilon \over 2}) \Delta t
   {\mathbf{D}}\underline{\tilde \delta}_n^m & = {1 \over g}
   \underline{\widetilde{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}}_n^m

Since :math:`{\mathbf{D}}` is diagonal, all equations are now
uncoupled in the vertical.

For each vertical mode, i.e. element of
:math:`(\underline{\tilde{~~}})_n^m`, and for each Fourier wavenumber
:math:`m` we have a system of equations in :math:`n` to solve. In
following we drop the Fourier index :math:`m` and the modal element
index :math:`(~~)_\ell` from the notation.

.. math::

   {{\cal A}}_n {\tilde\zeta}_n + {{{\cal B}}^+}_n {\tilde\delta}_{n+1} + {{{\cal B}}^-}_n
   {\tilde\delta}_{n-1} & = {\widetilde{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}}_n \\ {{\cal A}}_n
   {\tilde\delta}_n
   - {{{\cal B}}^+}_n {\tilde\zeta}_{n+1} {{{\cal B}}^-}_n {\tilde\zeta}_{n-1}
   -({1 + \epsilon \over 2}) \Delta t {n(n+1) \over
      a^2} g{\tilde h}_n
   & = {\widetilde{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}}_n \\ {\tilde h}_n + ({1 +
   \epsilon \over 2}) \Delta t {\rm D}_\ell{\tilde \delta}_n & = {1
   \over g} {\widetilde{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}}_n

The modal index :math:`(~~)_\ell` was included in the above equation on
:math:`{\rm D}` only as a reminder, but will also be dropped in the
following.

Substitute :math:`{\tilde\zeta}` and :math:`{\tilde h}` into the
:math:`{\tilde\delta}` equation.

.. math::

   [ {{\cal A}}_n + ( {1 + \epsilon \over 2})^2 (\Delta t) ^2 {n(n+1) \over a^2} g{\rm D}
   + {{{\cal B}}^+}_n {{\cal A}}_{n+1}^{-1} {{{\cal B}}^-}_{n+1} + {{{\cal B}}^-}_n {{\cal A}}_{n-1}^{-1}
   {{{\cal B}}^+}_{n-1} ] & {\tilde\delta}_{n} \cr + ( {{{\cal B}}^+}_n
   {{\cal A}}_{n+1}^{-1} {{{\cal B}}^+}_{n+1} ) {\tilde\delta}_{n+2} ~~+~~ (
   {{{\cal B}}^-}_n {{\cal A}}_{n-1}^{-1} {{{\cal B}}^-}_{n-1} ) {\tilde\delta}_{n-2}
   ~~~~~~~~~~~~~~~~~~~ & \\
   = \widetilde{{\hbox{\sffamily\slshape D}}{\hbox{\sffamily\slshape S}}}_n
   + ({1 + \epsilon \over 2}) \Delta t {n(n+1)
                  \over a^2} \widetilde{{\hbox{\sffamily\slshape H}}{\hbox{\sffamily\slshape S}}}_n
   + {{{\cal B}}^+}_n {{\cal A}}_{n+1}^{-1} \widetilde{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_{n+1} + {{{\cal B}}^-}_n
   {{\cal A}}_{n-1}^{-1} \widetilde{{\hbox{\sffamily\slshape V}}{\hbox{\sffamily\slshape S}}}_{n-1} \nonumber

which is just two tri-diagonal systems of equations, one for the even
and one for the odd :math:`n`\ ’s, and :math:`m \le n \le N`

At the end of the system, the boundary conditions are

.. math::

   {2}n & = m,   \quad {{{\cal B}}^-}_n     = {{{\cal B}}^-}_m^m = 0 \\ 
      n & = m+1, \quad {{{\cal B}}^-}_{n-1} = {{{\cal B}}^-}_m^m = {{{\cal B}}^-}_{(m+1)-1}^m = 0 \nonumber

the :math:`{\tilde\delta}_{n-2}` term is not present, and from the underlying truncation

.. math:: \tilde\delta_{N+1}^m = \tilde\delta_{N+2}^m = 0

For each :math:`m` and :math:`\ell` we have the general systems of
equations

.. math::

   -A_n {\tilde\delta}_{n+2} + B_n {\tilde\delta}_{n} -C_n
   -{\tilde\delta}_{n-2} & = D_n \;,
   \{
   \begin{array}{c}
     n=m,m+2,..., \{
   \begin{array}{c}
      N+1 \\ {\rm or} \\ N+2 \\
   \end{array}
   .~~~~~ \\ n=m+1,m+3,..., \{
   \begin{array}{c}
      N+1 \\ {\rm or} \\ N+2 \\
   \end{array}
   . \\
   \end{array}
   . \\ C_m = C_{m+1} & = 0 \\ \tilde\delta_{N+1} =
   \tilde\delta_{N+2} & = 0

Assume solutions of the form

.. math:: \tilde\delta_n = E_n \tilde\delta_{n+2} + F_n

then

.. math::

   E_m & = {A_m \over B_m} \\ F_M & = {D_m \over B_m}

.. math::

   E_n & = {A_n \over B_n - C_nE_{n-2}}~~,~~~ n=m+2,m+4,..., \{
                               \begin{array}{c}
                                 N-2 \\ {\rm or} \\ N-3 \\
                               \end{array} 
                            . \\
   F_n & = {D_n + C_nF_{n-2} \over B_n - C_nE_{n-2}}~~,~~~ n=m+2,m+4,...,
                            \{
                               \begin{array}{c}
                                 N \\ {\rm or} \\ N-1 \\
                               \end{array} 
                            .

.. math:: \tilde\delta_N = F_N~~~~{\rm or}~~~~\tilde\delta_{N-1} = F_{N-1} \;,

.. math::

   \tilde\delta_n = E_n\tilde\delta_{n+2} + F_n \;,~~ \{
        \begin{array}{c}
           n=N-2,N-4,..., \{
                        \begin{array}{c}
                             m \\ {\rm or} \\ m+1 \\
                        \end{array} .
            \\ n=N-3,N-5,..., \{
                         \begin{array}{c}
                             m+1 \\ {\rm or} \\ m \\
                         \end{array} .
            \\
        \end{array}
    .

Divergence in physical space is obtained from the vertical mode
coefficients by

.. math:: \underline{\delta}_n^m = {\mathbf{\Phi}}\underline{\tilde\delta}_n^m

The remaining variables are obtained in physical space by

.. math::

   \zeta_n^m ( 1 - {im\alpha \over n(n+1)}) & =
   {{\cal L}}_n^m
   - \delta_{n+1}^m \alpha ({n\over n+1}) D_{n+1}^m
   - \delta_{n-1}^m \alpha ({n+1\over n}) D_{n}^m \\
   \underline{T}^m_n & = \underline{\rm TS}^m_n
          - ({ 1 + \epsilon \over 2}) \Delta t {\mathbf{D}}^r
              \underline{\delta}^m_n \\
   {{\rm ln} p'_s}^m_n & = {\rm PS}^m_n
                       - ({ 1 + \epsilon \over 2})
                            {\Delta t\over p^r_s} (\underline{\Delta
                             p^r} )^T \underline{\delta}^m_n

Trajectory Calculation
~~~~~~~~~~~~~~~~~~~~~~

The trajectory calculation follows Let :math:`{\mathbf{R}}` denote
the position vector of the parcel,

.. math:: {d{\mathbf{R}} \over dt} = {\mathbf{V}}

which can be approximated in general by

.. math::

   {\mathbf{R}}_D^n = {\mathbf{R}}_A^{n+1} - \Delta t
   {\mathbf{V}}_M^{n+{1 \over 2}}

Hortal’s method is based on a Taylor’s series expansion

.. math::

   {\mathbf{R}}_A^{n+1} = {\mathbf{R}}_D^n + \Delta t ( {d
   {\mathbf{R}} \over d t} )_D^n +{\Delta t^2 \over 2} ( {d^2
   {\mathbf{R}} \over d t^2} )_D^{n} + \dots

or substituting for :math:`{d {\mathbf{R}} / d t}`

.. math::

   {\mathbf{R}}_A^{n+1} = {\mathbf{R}}_D^n + \Delta t {\mathbf{V}}_D^n
   +{\Delta t^2 \over 2} ( {d {\mathbf{V}} \over d t} )_D^{n}
   + \dots

Approximate

.. math::

   ( {d {\mathbf{V}} \over d t} )_D^{n} \approx {{\mathbf{V}}_A^n - {\mathbf{V}}_D^{n-1} \over \Delta t}\\[-1.0em] \intertext{giving}\nonumber\\[-2.0em] {\mathbf{V}}_M^{n+{1 \over 2}} 
   = {1 \over 2}[(2{\mathbf{V}}^n-{\mathbf{V}}^{n-1})_D + {\mathbf{V}}_A^n]

for the trajectory equation.

Mass and energy fixers and statistics calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The semi-Lagrangian dynamical core applies the same mass and energy
fixers and statistical calculations as the Eulerian dynamical core.
These are described in sections [massfixers], [energyfixer], and
[statcalc].

[f:sphere4] Tiling the surface of the sphere with quadrilaterals. An inscribed cube is projected to the surface of the sphere. The faces of the cubed sphere are further subdivided to form a quadrilateral grid of the desired resolution. Coordinate lines from the gnomonic equal-angle projection are shown. | image:: figures/sphere4-thick
[f:GLLnodes] A :math:`4 \times 4` tensor product grid of GLL nodes used within each element, for a degree :math:`d=3` discretization. Nodes on the boundary are shared by neighboring elements. | image:: figures/quad4
