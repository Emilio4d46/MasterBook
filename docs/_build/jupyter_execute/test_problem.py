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
import numpy as np


# We can model the problem as a set of concentric circles with the center most being our region of interest, $\Omega$. This is then sourrounded by a 'normal' region and then a PLM region.

# In[2]:


geo = SplineGeometry()
geo.AddCircle( (0,0), 2.25, leftdomain=3, bc="outerbnd")
geo.AddCircle( (0,0), 1.75, leftdomain=2, rightdomain=3, bc="midbnd")
geo.AddCircle( (0,0), 1, leftdomain=1, rightdomain=2, bc="innerbnd")
geo.SetMaterial(1, "inner")
geo.SetMaterial(2, "mid")
geo.SetMaterial(3, "pmlregion")
mesh = Mesh(geo.GenerateMesh (maxh=0.1))
mesh.Curve(3)

mesh.SetPML(pml.Radial(rad=1.75,alpha=1j,origin=(0,0)), "pmlregion") #Alpha is the strenth of the PML.


# We can set
# 
# $$c(x) = \begin{cases}
#    1 &\text{: } x\in\mathbb{R}^2/\Omega \\
#    100 &\text{: } x\in\Omega
# \end{cases}$$
# 
# and for $\omega(x) = 20$

# In[3]:


omega_0 = 20
omega_tilde = 0.2 #-18.5*exp(-(x**2+y**2)) + 20 #70*exp(-(((x*x)*np.log(7/2))+((y*y)*np.log(7/2)))) #Gaussian function for our test Omega.

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

# In[4]:


fes = H1(mesh, complex=True, order=5)

def forward_pass(theta):

    u_in =exp(1j*omega_0*(cos(theta)*x + sin(theta)*y)) #Can use any vector as long as it is on the unit circle. 

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
    f += -u_in * (omega**2 - omega_0**2) * v * dx
    f.Assemble()

    #Solving our problem.
    u_s = GridFunction(fes, name="u")
    u_s.vec.data = a.mat.Inverse() * f.vec

    u_tot = u_in + u_s
    
    return [u_in, u_s, u_tot]

s = 1*np.pi

calc = forward_pass(s)
u_in = calc[0]
u_s = calc[1]
u_tot = calc[2]


# In[5]:


Draw(u_in, mesh, "u_in")
Draw(u_s, mesh, "u_s")
Draw(u_tot, mesh, "u_tot")


# We can now start to compare between the far field of the wave as computed via the PML and the first order approximation that was used in the original paper. This far field approximation is namely
# 
# $$\hat{u}^s(r) = \int_{\Omega}e^{-i\omega_0(\hat{r} - \hat{s})\cdot x}\eta(x)\,dx$$
# 
# Which is over $\Omega$ as $\eta(x)$ is only supported in $\Omega$. While for the full far field wave we can use the following equation
# 
# $$u^s_{\infty}(\hat{r}) = \dfrac{e^{\frac{\pi i}{4}}}{\sqrt{8\pi k}} \int_{\partial\Omega}u(x)\dfrac{\partial e^{-i\omega_0\hat{r}\cdot x}}{\partial n} - \dfrac{\partial u(x)}{\partial n}e^{-i\omega_0\hat{r}\cdot x}\,dS(x)$$
# 
# We can use NGSolve to compute both of these and comapre the result.

# In[8]:


def ff_field(r):
    
    n = specialcf.normal(2)
    us_n = BoundaryFromVolumeCF(Grad(u_s)*n)
    ecomp = exp(-1j*omega_0*(cos(r)*x + sin(r)*y))
    ecomp_n = CoefficientFunction((ecomp.Diff(x),ecomp.Diff(y)))*n
    
    temp = Integrate(u_s*ecomp_n - us_n*ecomp, mesh,definedon=mesh.Boundaries("innerbnd"))
    
    return (exp(1j*r)/np.sqrt(8*np.pi*omega_0))*temp

theta = np.arange(0, 2*np.pi, 0.01)

mag1 = []

for r in theta:
    
    mag1.append(abs(ff_field(r)))
    
    
fig = plt.figure(figsize=(9, 5))

ax = plt.subplot(1, 1, 1, projection='polar')

plt.title('Full Far-field pattern for s = pi')

ax.plot(theta, mag1)
ax.set_rmax(abs(ff_field(s)))
ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

plt.show()


# In[ ]:





# In[ ]:




