Derivation of Foundational Equations
=======================


In order to fully illuminate the methodology and structure behind the model that Yuwei Fan et al. decided to use we must first understand the underlying mathematics behind the scenario they wished to model. I will be following their derivation here and expanding upon it where I find components of their derivation confusing.

  

**Definitions and Assumptions**

Let us define $\Omega$ as a compact domain of interest. We can define the inhomogeneous media scattering problem with a fixed frequency wave, $\omega$, to be modelled by the Helmholtz equation

$$Lu:= \left(-\Delta- \frac{\omega^2}{c^2(x)}\right)$$

with $c^2(x)$ being an unknown velocity field. Assume there exists a known background velocity, $c_0(x)$ such that $c(x)$ is identical to $c_0(x)$ outside the domain $\Omega$.

We define the scatterer, $\eta(x)$ as

$$\eta(x) := \left(\frac{\omega^2}{c^2(x)} - \frac{\omega^2}{c_0^2(x)}\right)$$

which we can see is compactly supported in $\Omega$ as $c^2(x)$ is equal to $c_0^2(x)$ outside $\Omega$. 

We can now define the operator

$$L_0:= -\Delta - \frac{\omega^2}{c_0^2(x)} =  -\Delta - \omega^2 $$

with out loss of generality, we can set $c_0(x)$ to be constant, for convenience we set $c_0(x)\equiv1$. We may then also redefine 

$$L:=L_0-\eta(x)$$

The sources can be seen to be incident anywhere on the boundary of $\Omega$ and denoted as $s$ and be parameterised by $s\in[0,2\pi)$.  Furthermore, we can represent each source wave as an incoming plane wave

$$e^{i\omega\hat{s}x}$$

with $\hat{s} = (\cos(x), \sin(x))$ or $\hat{s} \in \mathbb{S}^1$ being the direction of incidence. We denote the scattered wave by $u^s(x)$ which satisfies the following equation.

$$\left(L_0 - \eta(x)\right) \left(e^{i\omega\hat{s}x} + u^s(x)\right) = 0 \tag{1}$$ 

As well as the Somerfield radiation boundary condition

$$\lim_{|x|\to\infty}|x|^{\frac{1}{2}}\left(\frac{\partial}{\partial|x|} - ik\right)u(x)=0$$

 This is true as TBD.[^1]

Similar to the sources, we can denote the 'recorded' far field pattern $r$, and let them  be parameterised by $r \in [0,2\pi)$. Once again, we will denote $\hat{r} = (\cos(x), \sin(x))$ or $\hat{r} \in \mathbb{S}^1$ being the direction of emission. With this we can define $\hat{u}^s(x)$ to be the far field pattern recorded in direction $\hat{r}$ as

$$\hat{u}^s(x)\equiv\lim_{\rho\to\infty}\sqrt{\rho}e^{-i\omega\rho}u^s(\rho\cdot\hat{r})$$

We can then denote $d(s,r)$ to be the set of all far field patterns from a source $s\in S$ to a direction $r\in R$ by $d(s,r)\equiv\hat{u}^s(x)$

**Relationship between source waves and scattered wave for far field waves**

In order to solve either the forward map or inverse problem we must better understand the relationship between $\eta(x)$ and $d(s,r)$. We can do this by doing a perturbative analysis for a small $\eta(x)$. We can do this by expanding (1)

 $$\left(L_0 - \eta(x)\right) \left(e^{i\omega\hat{s}x} + u^s(x)\right) = 0$$
 $$\Leftrightarrow \left(L_0 - \eta(x)\right) u^s(x) = \left(L_0 - \eta(x)\right)e^{i\omega\hat{s}x}$$
 $$\Leftrightarrow L_0 \cdot u^s(x)-\eta(x)\cdot u^s(x) = L_o \cdot e^{i\omega\hat{s}x}-\eta(x)\cdot e^{i\omega\hat{s}x}$$
 As $u^s(x)$ is dependent on $\eta(x)$ and we are looking only at small perturbations from $\eta(x) = 0$ we can combine them into higher order $\eta(x)$ terms
  $$L_0 \cdot u^s(x)+o(\eta)= L_o \cdot e^{i\omega\hat{s}x}-\eta(x)\cdot e^{i\omega\hat{s}x}$$
  
And we know that $L_0 \cdot e^{i\omega\hat{s}x} = 0$, so we arrive at

$$\left(L_0 \cdot u^s\right)(x)= \eta(x)\cdot e^{i\omega\hat{s}x}+o(\eta)\tag{2}$$

We can now let the inverse of the operator, $L_0^{-1} = G_0$, to be the integral kernel of Green functions of the free space Helmholtz operator, $L_0$. Doing so, as well as allowing $o(\eta(x))\to0$ as we are looking at $\eta(x)$ close to 0, can rewrite (2) as

$$u^s(x) = \int G_0(y-x)\eta(x)e^{i\omega\hat{s}\cdot x}\, dx\tag{3}$$

We can now apply the Green function of the free space Helmholtz operator expansion at infinity

$$G_0(z) = \frac{e^{iwz}+o(1)}{\sqrt{|z|}}$$

Applying this to (3) and combining with equation (1), we can compute $\hat{u}^s(x)$

$$\hat{u}^s(x)=\lim_{\rho\to\infty}\sqrt{\rho}e^{-i\omega\rho}u^s(\rho\cdot\hat{r})\approx\lim_{\rho\to\infty}\sqrt{\rho}e^{i\omega\rho}
\int\left(\frac{e^{iw(\rho-\hat{r}\cdot x)}}{\sqrt{\rho}}\eta(x)e^{-i\omega\hat{s}\cdot x}\right)\, dx$$

Which we can expand it out and simplify to

$$=\int e^{-iw(\hat{r}-\hat{s})\cdot x}\eta(x)\, dx = d_1(s,r)$$

With $d_1(s,r)$ denoting the first order approximation to $d(s,r)$ with respect to $\eta(x)$.

**Problem Setup**

As we are dealing with the far field pattern problem, we can approximate $\Omega$ as the unit disk centred at the origin and can treat $\hat{s}$ and $\hat{r}$ as uniformly sampled from $\mathbb{S}^1$. We can denote this by setting $s=\frac{2\pi j}{N_s}$ and $r=\frac{2\pi k}{N_r}$ with setting $N_s = N_r$. As we are now assuming that $\Omega$ is the unit disc, it will now be useful to use polar coordinates. We may define $x=\left(\rho\cos(\theta),\rho\sin(\theta)\right)$ with $\rho \in [0,1]$ and $\theta \in [0,2\pi)$. As we also have the direction of the source, $s\in[0,2\pi)$ and the direction of our scattering measurement, $r\in[0,2\pi)$; it is suggested to use a change of variables[^2]

$$m=\frac{r+s}{2},h=\frac{r-s}{2}$$

With $m,h,r,s$ are all modulo $2\pi$. We can also trivially redefine

$$d(s,r) \text{ to } d(m,h)$$

and

$$\eta(x)\text{ to }\eta(\rho,\theta)$$

in our new polar coordinates. As $d_1(m,h)$ is linear in terms of $\eta(\rho,\theta)$ there exists a linear operator mapping

$$d_1(m,h) \to \eta(\rho,\theta)$$

As this is a mapping from a function space to a function space there exists a kernel operator $K(m,h,\rho,\theta)$ such that

$$d_1(m,h) = \int_{0}^{1} \int_{0}^{2\pi} K(m,h,\rho,\theta)\cdot \eta(\rho,\theta)\,\, d\theta \, d\rho$$

As the external velocity field is constant, $c_0=1$, and the problem is defined on the unit disc centred at the origin; the system is symmetric with respect to $m$ and $h$.[^3]
This allows us to simplify the system significantly.

[^1]: This is related to the boundary conditions and the fact that the interior and exterior solution should be continuous at the boundary when we expand. However, I am still unsure as to how this relates to the expression (1). I am still making my way through Taylor but if there is a more elegant explanation or resource please do link it to me.

[^2]: Why do we change our variables in this manner? What is the reason for not just using $s$ and $r$ as they are defined above?

[^3]: I understand that the system is fully symmetric with respect to $s$ as there is not preference to the incidence wave. But surely due to the internal structure of $\Omega$ would vary with respect to $r$. Or in a more to phrase it more mathematically, the variation of $\eta(\rho,\theta)$ within $\Omega$ would change $\hat{u}^s(\rho,\theta)$ so that at a given point $r$, we would measure a different wave?

