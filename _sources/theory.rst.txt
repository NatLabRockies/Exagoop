Numerical Methods
=================

This section presents the continuum mechanics governing equations that underpin ExaGOOP,
followed by their discretisation within the Material Point Method (MPM) framework.

Governing equations of continuum mechanics
------------------------------------------

The fundamental equations that underpin the Material Point Method (MPM) are the conservation of mass and linear momentum. The energy equation is handled separately through the optional temperature module (see :ref:`temperature_module`). Additionally, material elements with internal torques are not considered here, so the angular momentum conservation equation is omitted. In the context of MPM, the laws of mass and momentum conservation are expressed within a Lagrangian framework.

The Lagrangian description of the motion of a particle is expressed in terms of its material coordinates :math:`X` and time :math:`t`. Material coordinates refer to a coordinate system attached to the particle under consideration in its initial configuration (time :math:`t=0`). In the Lagrangian description, the particle is assumed to move with the local velocity of the medium, and other continuum properties are studied in this coordinate system.

Hence, the motion of a particle is expressed in the Lagrangian description as,

.. math::
	
	\begin{aligned}
	\mathbf{x} = \mathbf{x}(\mathbf{X},t)
	\end{aligned}

It is to be noted that by definition, the above particle has the coordinates defined by :math:`\mathbf{X}` at time :math:`t=0`. The displacement :math:`\mathbf{u}` of the particle with respect to the initial configuration is then expressed as,

.. math::
	
	\begin{aligned}
	\mathbf{u} = \mathbf{x}(\mathbf{X},t)-\mathbf{X}
	\end{aligned}

The definition of the velocity of the particle :math:`\mathbf{v}` then follows as,

.. math::
	
	\begin{aligned}
	\mathbf{v}  &= \dot{\mathbf{x}} =  \frac{d \mathbf{u}(\mathbf{X},t)}{dt}\\
	&= \frac{\partial \mathbf{u}(\mathbf{X},t)}{\partial t}
	\end{aligned}

Similarly, the time rate of change of any property of the particle expressed in the Lagrangian framework is its partial time derivative, and :math:`\mathbf{X}` simply serves as a parameter. Similar to the expression of velocity given above, the acceleration of the particle can also be expressed as,

.. math::
	
	\begin{aligned}
	\mathbf{a}  &= \dot{\mathbf{v}} =  \frac{\partial \mathbf{v}(\mathbf{X},t)}{\partial t}\\
	&= \frac{\partial^2 \mathbf{u}(\mathbf{X},t)}{\partial t^2}
	\end{aligned}

Without going deep into the details, some of the mathematical terms used in the governing equations are defined in the following paragraphs.

Deformation gradient tensor
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The deformation gradient tensor :math:`\mathbf{F}` is defined as,

.. math::
	
	\begin{aligned}
	\mathbf{F}  = \frac{\partial \mathbf{x}(\mathbf{X},t)}{\partial \mathbf{X}}
	\end{aligned}

and is a symmetric, second-order tensor which describes the stretch and rotation of a material element. Mathematically, it is a linear operator that maps the current configuration of a continuum body to its initial configuration.

Velocity gradient tensor
~~~~~~~~~~~~~~~~~~~~~~~~

The velocity gradient tensor :math:`\mathbf{L}` is defined as the spatial gradient of velocity.

.. math::
	
	\begin{aligned}
	\mathbf{L} = \frac{\partial \mathbf{v}}{\partial x}
	\end{aligned}

This second-order tensor can be decomposed into a symmetric part (rate of deformation tensor) and an anti-symmetric part (spin tensor) as shown below.

.. math::
	
	\begin{aligned}
	\mathbf{L} = \mathbf{D}+{\Omega}
	\end{aligned}

where,

.. math::
	
	\begin{aligned}
	\mathbf{D} &= \frac{1}{2} (\mathbf{L}+\mathbf{L}^T)\\
	\Omega &= \frac{1}{2} (\mathbf{L}-\mathbf{L}^T)\\
	\end{aligned}

The rate of deformation tensor :math:`D` indicates the rate of strain suffered by a material element and is used to find the stresses through a constitutive model. The spin tensor :math:`\Omega` refers to the rotation the material element undergoes. The velocity gradient tensor is related to the deformation gradient tensor through the following expression,

.. math::
	
	\begin{aligned}
	\mathbf{L} = \dot{\mathbf{F}}\mathbf{F}^{-1}
	\end{aligned}

Jacobian
~~~~~~~~

The jacobian (:math:`\mathbf{J}`) is defined as the determinant of the deformation gradient tensor (:math:`\mathbf{F}`).

.. math::
	
	\begin{aligned}
	J=\left|\frac{\partial \mathbf{x}}{\partial \mathbf{X}}\right| = \left| \mathbf{F} \right|
	\end{aligned}

A necessary and sufficient condition for the motion to be invertible is to have a non-zero Jacobian at all times. The Jacobian also relates the volume of an infinitesimal body at time :math:`t` to its volume at the initial time through the relation

.. math::
	
	\begin{aligned}
	\mathrm{d} V=\left|\begin{array}{ccc}\mathrm{d} x_1 & \mathrm{~d} x_2 & \mathrm{~d} x_3 \\ \delta x_1 & \delta x_2 & \delta x_3 \\ \Delta x_1 & \Delta x_2 & \Delta x_3\end{array}\right|=J \mathrm{~d} V_0
	\end{aligned}

where :math:`\mathrm{d} V` and :math:`\mathrm{d} V_0` are the volumes at current and initial time.

Stress tensor and Constitutive models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stress tensor :math:`\mathbf{\sigma}` is a symmetric, second-order tensor that defines the state of stress at a point. The traction or force per unit area acting at the point on an imaginary surface with normal :math:`\mathbf{n}` is related to the stress tensor at the point as,

.. math::
	
	\begin{aligned}
	\mathbf{t} = \mathbf{n}.\mathbf{\sigma}
	\end{aligned}

Stress tensor at a point is related to the rate of deformation tensor :math:`\mathbf{D}` through a constitutive relation. Depending on the material considered, a multitude of constitutive relations, such as linear elastic, plastic, and Newtonian fluids, exist.

Equation of conservation of mass
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a continuum body occupying a region :math:`\Omega` in space and bounded by a surface :math:`\Gamma=\Gamma_u \cup \Gamma_t` as shown in :numref:`compdom`, the total mass :math:`m` of the body is given by,

.. figure:: ../landing/_images/MPMEq.png
	:name: compdom
	:align: center
	:figwidth: 50%

	: Computational domain, boundary and various forces acting on it

.. math::
	
	\begin{aligned}
	m=\int_{\Omega} \rho(\mathbf{x}, t) d V
	\end{aligned}

where :math:`\rho` is the density of the material.

Since the mass contained in the region :math:`\Omega` and moving with local material velocity is constant, the total time derivative of the total mass  is zero. Hence,

.. math::
	
	\begin{aligned}
	\frac{\mathrm{D}}{\mathrm{D} t} \int_{\Omega} \rho(\mathbf{x}, t) \mathrm{d} V=\int_{\Omega}(\dot{\rho}+\rho \nabla \cdot \boldsymbol{v}) \mathrm{d} V=0
	\end{aligned}

which leads to the mass conservation equation as,

.. math::
	
	\begin{aligned}
	\rho J-\rho_0=0
	\end{aligned}

Conservation of momentum
~~~~~~~~~~~~~~~~~~~~~~~~

Newton’s second law of motion states that the rate of change of momentum of a body is equal to the sum of the volume and surface forces acting on it. Consider the same body as shown in :numref:`compdom`, with a body force per unit mass :math:`\mathbf{b}` and traction :math:`\mathbf{t}` acting on it’s surface :math:`\Gamma_t`. The law of conservation of momentum is expressed as,

.. math::
	
	\begin{aligned}
	\frac{\mathrm{D}}{\mathrm{D} t} \int_{\Omega} \rho \mathbf{v} \mathrm{d} V=\int_{\Omega} \rho \mathbf{b}(\mathbf{x},t) \mathrm{d} V + \int_{\Gamma_t} \mathbf{t}(\mathbf{x},t).\mathbf{n} \mathrm{d} A	
	\end{aligned}

By invoking the Reynolds transport theorem and upon simplifying, one obtains,

.. math::
	
	\begin{aligned}
	\rho \dot{\mathbf{v}} = \rho \mathbf{b} + \nabla . \sigma
	\end{aligned}

Initial and Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two types of boundary conditions are commonly considered, namely boundaries with specified velocity and specified traction, respectively.

.. math::
	
	\begin{aligned}
	& \left\{\begin{array}{l}
	\left.(\boldsymbol{n} \cdot \boldsymbol{\sigma})\right|_{\Gamma_t}=\overline{\boldsymbol{t}} \\
	\left.\boldsymbol{v}\right|_{\Gamma_u}=\overline{\boldsymbol{v}}
	\end{array}\right. \\
	\end{aligned}

Initial conditions involve specifying the displacement of material points and velocities at time t=0.

.. math::
	
	\begin{aligned}
	& \mathbf{v}(\mathbf{x}, 0)=\mathbf{v}_0(\mathbf{x}), \quad \mathbf{u}(\mathbf{x}, 0)=\mathbf{u}_0(\mathbf{x})
	\end{aligned}

.. note::

   The equations above define the strong form of the governing equations.
   The following section derives their weak form and the MPM discretisation
   used in ExaGOOP.


Governing Equations of MPM
--------------------------

The complete set of governing equations for a non-isothermal problem involving non-polar materials, as discussed in the previous section, is summarized as,

.. math::
   :label: Eq_mce

   \rho J - \rho_0 = 0
 
   
 
.. math::
	
   \begin{aligned}
   \rho \dot{\mathbf{v}} = \rho \mathbf{b} + \nabla . \sigma
   \end{aligned}

The equation of conservation of mass is not solved explicitly. Instead it is used to calculate the updated density field in the computations through Eq. :eq:`Eq_mce`. Like in many finite-element formulations, the method of weighted residuals is used to reduce the residual to zero in an average sense. By taking the virtual
displacement :math:`\delta u_j \in \Re_0, \Re_0=\left\{\delta u_j\left|\delta u_j \in C^0, \delta u_j\right|_{\Gamma_u}=0\right\}` as the test function, one obtains the weak form of the governing equations and the traction boundary condition as,

.. math::

   \begin{aligned}
       \begin{array}{r}
   \int_{\Omega} \delta u_i\left(\sigma_{i j, j}+\rho b_i-\rho \ddot{u}_i\right) \mathrm{d} V=0, \\
   \int_{\Gamma_t} \delta u_i\left(\sigma_{i j} n_j-\bar{t}_i\right) \mathrm{d} V=0,
   \end{array}
   \end{aligned}

which, upon simplification, is shown below.

.. math:: :label: Eq_mom

   \begin{aligned}
   \int_{\Omega} \rho \ddot{u}_i \delta u_i \mathrm{~d} V+\int_{\Omega} \rho \sigma_{i j}^s \delta u_{i, j} \mathrm{~d} V-\int_{\Omega} \rho b_i \delta u_i \mathrm{~d} V-\int_{\Gamma_t} \rho \bar{t}_i^s \delta u_i \mathrm{~d} A=0
   \label{Eq:GovEqMPM}
   \end{aligned}

In the standard formulation of MPM used in ExaGOOP, the density field in the domain is approximated as,

.. math:: :label: Eq_rho

   \begin{aligned}
   \rho(\boldsymbol{x})=\sum_{p=1}^{n_p} m_p \delta\left(\boldsymbol{x}-\boldsymbol{x}_p\right)   
   \end{aligned}

In the above equation, :math:`{n_p}` is the number of material points and :math:`{m_p}` is the mass of each material point.
Substituting Eq. :eq:`Eq_rho` in Eq. :eq:`Eq_mom`  and invoking particle quadrature, one obtains the numerical governing equation of MPM.

.. math:: :label: Eq_NumEqMPM

   \begin{aligned}
   \sum_{p=1}^{n_p} m_p \ddot{u}_{i p} \delta u_{i p}+\sum_{p=1}^{n_p} m_p \sigma_{i j p}^s \delta u_{i p, j}-\sum_{p=1}^{n_p} m_p b_{i p} \delta u_{i p}-\sum_{p=1}^{n_p} m_p \bar{t}_{i p}^s h^{-1} \delta u_{i p}=0   
   \end{aligned}

The subscipts :math:`p` and :math:`i` refer to the particle and spatial dimension, respectively. The solution of the numerical governing equation
above is carried out in four stages, as briefly discussed in the following section.

.. _mpmsect:

MPM discretisation steps
--------------------------
The MPM solution procedure in ExaGOOP begins by generating or reading a collection of material point data from an input file. This material point data includes information on each point's location, velocity, and flags for the constitutive model. Additionally, a background grid is generated based on user-defined inputs, which specify the entire computational domain. The time integration process then advances through each time step, with each step comprising four distinct sub-steps, as outlined in the following subsections.

.. raw:: html

   <div class="figure-grid" style="display: flex; flex-wrap: wrap; gap: 1em; justify-content: center;">
     <div style="flex: 0 0 45%;">
       <img src="./_images/MPM_Step1.png" alt="step 1" style="width: 100%;">
       <p style="text-align: center;"><strong>(a)</strong> particle to grid (P2G) operation</p>
     </div>
     <div style="flex: 0 0 45%;">
       <img src="./_images/MPM_Step2.png" alt="step 2" style="width: 100%;">
       <p style="text-align: center;"><strong>(b)</strong> nodal velocity update </p>
     </div>
     <div style="flex: 0 0 45%;">
       <img src="./_images/MPM_Step4.png" alt="step 3" style="width: 100%;">
       <p style="text-align: center;"><strong>(c)</strong> grid to particle (G2P) operation</p>
     </div>
     <div style="flex: 0 0 45%;">
       <img src="./_images/MPM_Step3.png" alt="step 4" style="width: 100%;">
       <p style="text-align: center;"><strong>(d)</strong> particle position update</p>
     </div>
   </div>
.. figure:: ../landing/_images/none.png
   :name: fig-grid
   :height: 0
   :width: 0   
   :figwidth: 100%
   :align: center   
   :alt: Simulation Snapshots Overview

   The four steps involved in one step of MPM time integration. The material points are shown as circles in red color. The nodes are shown as squares with black outline

Particle to Grid Interpolation (P2G)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
In this step, the material points are assumed to be attached to the background grid as shown in :numref:`fig-grid` (a). The background grid is then considered similar to a finite element grid and based on the shape function defined at the grid node :math:`I` the unknown quantities in Eq. :eq:`Eq_NumEqMPM` are calculated as,

.. math:: :label: Eq_shape

   \begin{aligned}
   \begin{array}{r}
   u_{i p}  =N_{I p} u_{i I} \\
   u_{i p, j}  =N_{I p, j} u_{i I}\\
   \delta u_{i p} =N_{I p} \delta u_{i I}
   \end{array}   
   \end{aligned}

In the equations above, subscript :math:`I` is used to denote the grid node and :math:`N_I` indicate the shape function defined at node :math:`I`. ExaGOOP supports linear, quadratic B-spline and cubic B-spline shape functions. Substituting equations Eq. :eq:`Eq_shape` in Eq. :eq:`Eq_NumEqMPM` and cancelling the common virtual displacement term :math:`\delta u_{i I}`, one obtains,

.. math::

   \begin{aligned}
   m_{I J} \dot{u}_{i J}=f_{i I}^{\mathrm{int}}+f_{i I}^{\mathrm{ext}}, \quad x_I \notin \Gamma_u
   \end{aligned}

where :math:`m_{I J}` is the elements of the mass matrix defined as,

.. math::

   \begin{aligned}
   m_{I J}=\sum_{p=1}^{n_p} m_p N_{I p} N_{J p}\\
   \end{aligned}

and :math:`f_{i I}^{\mathrm{int}}` and :math:`f_{i I}^{\mathrm{ext}}`
are the internal and external forces respectively and given by,

.. math::

   \begin{aligned}
   f_{i I}^{\mathrm{int}}=-\sum_{p=1}^{n_p} N_{I p, j} \sigma_{i j p} \frac{m_p}{\rho_p} \\
   f_{i I}^{\mathrm{ext}}=\sum_{p=1}^{n_p} m_p N_{I p} b_{i p}+\sum_{p=1}^{n_p} N_{I p} \bar{t}_{i p} h^{-1} \frac{m_p}{\rho_p}
   \end{aligned}

Hence, this step of the MPM solution procedure involves ’projecting’
properties from material points to grid nodes and is shown schematically
in :numref:`fig-grid` (a).

Temporal integration at grid nodes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the grid nodal properties are calculated in the P2G operation, the updated velocity at grid nodes are calculated. In ExaGOOP, an explicit, Euler time integration procedure is used. Since this procedure in its original form involves costly inversion of the mass matrix :math:`m_{I J}`, the following mass-lumping approximation is made,

.. math::

   \begin{aligned}
   m_I=\sum_{J=1}^{n_g} m_{I J}=\sum_{p=1}^{n_p} m_p N_{I p}
   \end{aligned}

The velocity components at the nodes are then calculated as,

.. math::

   \begin{aligned}
   \mathbf{v_{I}}^{t+\Delta t}=\mathbf{v_{I}}^{t}+\frac{1}{m_I} \left(\mathbf{f_{i I}}^{\mathrm{int}}+\mathbf{f_{i I}}^{\mathrm{ext}}\right)
   \end{aligned}

where :math:`\Delta t` is the time step used in time integration and is
calculated from the following equation,

.. math::

   \begin{aligned}
   \Delta t= CFL \min \left(\frac{h_x}{c_x}, \frac{h_y}{c_y}, \frac{h_z}{c_z}\right)
   \end{aligned}

Here, :math:`c_{()}` and :math:`h_{()}` refer to characteristic velocity and
grid sizes in different directions respectively.

Grid to Particle (G2P) Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the updated velocities at grid nodes are obtained, the velocities
and their gradients at the material points are obtained in this step as,

.. math::

   \begin{aligned}
   \mathbf{v}_p^{t+\Delta t}=\alpha_{P-F}\left(\mathbf{v}_p^t+\sum_I N_I \left[{\mathbf{v}}_I^{t+\Delta t}-\mathbf{v}_I^t\right]\right)+(1-\alpha_{P-F}) \sum_I N_I {\mathbf{v}}_I^{t+\Delta t}\\
   \nabla \mathbf{v}_p^{t+\Delta t}=\sum_I^{ng} \nabla N_I \mathbf{v}_I^{t+\Delta t}
   \end{aligned}

The term :math:`\alpha_{P-F}` used in the material point velocity update
step above determines the level of blending between Particle-in-Cell
(PIC) and Fluid Implicit Particle Method (FLIP) like updates. The
velocity gradient thus calculated at the material point is used to
compute the stress tensor through the user-provided constitutive
relation.

Material point position update and grid reset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At this step, the updated velocity at the material point is already
obtained and is used to update the material point position as,

.. math::

   \begin{aligned}
   \mathbf{x}_{p}^{t+\Delta t}=\mathbf{x}_{p}^{t} +\Delta t \: \mathbf{v}_p^{t+\Delta t}
   \end{aligned}

The background grid in MPM is used only as a scratch pad to calculate
gradients and for time integration and hence is often reset or
regenerated at the end of each MPM step.

.. container:: float
   :name: Fig:MPM_Steps

   | 
