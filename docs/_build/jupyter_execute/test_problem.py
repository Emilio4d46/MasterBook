#!/usr/bin/env python
# coding: utf-8

# Test Problem
# =======================
# 
# Now that we have gotten a firm understandting of the PML we can start implementing a forward pass: this is to generate training data for our model. We first, however, implement a test case using the same tools as we hope to in the main problem, NGSolve within Python. The test problem will be to model a plan wave in a simple domain.

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *
import matplotlib.pyplot as plt


# We can model the problem as a set of concentric circles with the center most being our region of interest, $\Omega$. This is then sourrounded by a 'normal' region and then a PLM region.

# In[2]:


geo = SplineGeometry()
geo.AddCircle( (0,0), 1.7, leftdomain=3, bc="outerbnd")
geo.AddCircle( (0,0), 1.25, leftdomain=2, rightdomain=3, bc="midbnd")
geo.AddCircle( (0,0), 1, leftdomain=1, rightdomain=2, bc="innerbnd")
geo.SetMaterial(1, "inner")
geo.SetMaterial(2, "mid")
geo.SetMaterial(3, "pmlregion")
mesh = Mesh(geo.GenerateMesh (maxh=0.05))
mesh.Curve(3)

mesh.SetPML(pml.Radial(rad=1.25,alpha=5j,origin=(0,0)), "pmlregion") #Alpha is the strenth of the PML.


# We can set
# 
# $$w(x) = \begin{cases}
#    w_0 &\text{: } x\in\mathbb{R}^3/\Omega \\
#    \tilde{w}(x) &\text{: } x\in\Omega
# \end{cases}$$
# 
# to be the wave number. Here we have set our $\tilde{w}(x)$ to a multivariate Gaussian with equation:
# 
# $$\tilde{w}(x) = 50\exp{\left(-\frac{x^2 + y^2}{\log\left({\frac{5}{2}}\right)}\right)}$$
# 
# and for $w_0(x) = 20$

# In[12]:


import numpy as np

omega_0 = 20
omega_tilde = 50*exp(-(((x*x)*np.log(5/2))+((y*y)*np.log(5/2)))) #Gaussian function for our test Omega.

domain_values = {'inner': omega_tilde, 'mid': omega_0, 'pmlregion': omega_0}
values_list = [domain_values[mat] for mat in mesh.GetMaterials()]
omega = CoefficientFunction(values_list)

Draw(omega, mesh, 'piecewise')


# We can also define the wave as
# 
# $$u(x) = u^\text{inc}(x) + u^{\text{s}}(x)$$
# 
# with $u^\text{inc}(x)$ being the incident wave and $u^{\text{s}}(x)$ being the scattered wave. With
# 
# $$u^\text{inc}(x) = e^{iw_0d\cdot x}$$
# 
# With $d$ being a direction vector and $\|d\|_2 = 1$, or in other words, $d\in\mathbb{S}^1$. We then apply $u(x)$ to a Helmholtz Equation:
# 
# $$\Delta u(x) + w^2(x)u(x) = 0$$
# 
# We wish to solve an inhomogenous Helmholtz euqation, $\Delta u(x) + \omega^2u(x) = f(x)$, so that we can use classiacal finite element methods. We can do this by expanding it via $u(x) = u^\text{inc}(x) + u^{\text{s}}(x)$. This gives us:
# 
# $$\Delta\left(u^\text{inc}(x) + u^{\text{s}}(x)\right) + w^2(x)\left(u^\text{inc}(x) + u^{\text{s}}(x)\right) = 0$$
# 
# Rewriting this with the terms of the scattered wave on the left hand side and with the incident wave on the right hand side we get:
# 
# $$\Delta u^\text{s}(x) + w^2(x)u^{\text{s}}(x) = \underbrace{-\left(\Delta u^\text{inc}(x) + w^2(x)u^{\text{inc}}(x)\right)}_{f(x)}$$
# 
# As we know that $u^\text{inc}(x) = e^{iw_0d\cdot x}$, so, $\Delta u^\text{inc}(x) = -w_0^2(x)u^{\text{inc}}(x)$. This allows us to re-write $f(x)$ and the problem as:
# 
# $$\Delta u^\text{s}(x) + w^2(x)u^{\text{s}}(x) = -u^\text{inc}(x)\left(\omega^2 - \omega_0^2\right)$$
# 
# Which we can solve using classiacal finite element methods.
# 
# _**Note**: We can observe that the right hand side only has support in $\Omega$ as outside of $\Omega$ we have that $w(x) = w_0$, hence $f(x) = 0$._

# In[ ]:


fes = H1(mesh, complex=True, order=5)

u_in = exp(1j*omega_0*(3/5*x + 4/5*y)) #Can use any vector as long as it is on the unit circle. 

#Defining our test and solution functions.
u = fes.TrialFunction()
v = fes.TestFunction()

#Defining our LHS of the problem.
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx - omega**2*u*v*dx
a += -omega*1j*u*v * ds("innerbnd")
a.Assemble()

#Defining the RHS of our problem.
f = LinearForm(fes)
f += u_in * (omega**2 - omega_0**2) * v * dx
f.Assemble()

#Solving our problem.
u_s = GridFunction(fes, name="u")
u_s.vec.data = a.mat.Inverse() * f.vec


# In[14]:


Draw (u_in, mesh, "u_in")


# In[15]:


Draw(u_s, mesh, "u_s")


# In[16]:


Draw(u_in + u_s, mesh, "u_tot")

