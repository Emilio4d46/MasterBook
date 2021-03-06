
Perfectly Matched Layers
=======================

Before we can start to build and train our model, we first need to get training data that we will train out model with. We can do this either by experimentation or by modelling the forward pass of the wave interacting with the scattering field. The former is likely out of the scope of the budget for this masters project and, most likely, far more difficult. Modelling the forward pass, however, will not only be free but also provide us a greater appreciation and understanding of computational modelling and the problem at hand. 

One of the largest issues we will face will be modelling the far field behaviours of the waves as, since we cannot model an infinite area for the waves to move through, we will have to have some form of hard boundary which will cause reflection and interreference with the scattering that we would wish to observe. Typical solutions to problems of this ilk are usually along the lines of either remapping an infinite domain to a finite one, typically via $\tanh$ or by creating a large enough modelling space in which the waves will decay naturally by the boundary. However both of these solutions pose issues when it comes to oscillating waves as remapping them will cause them to oscillate infinitely at the far field and as they only decay with respect to $\frac{1}{r}$ they decay too slowly to model in a reasonable domain. Thankfully, there is a solution proposed in "A perfectly matched layer for the absorption of electromagnetic waves" by Jean-Pierre Berenger in 1994 [^1] dubbed the 'Perfectly Matched Layer' (PML). He proposes a method to model waves going off to $\infty$ independent of test set up. This method proposes applying a totally absorbing boundary layer at the boundary of the modelling domain and the waves in this layer exponentially decay; an issue being that there is almost always some reflection when we change media. However, we may analytically extending the wave equation into the complex plane via a form of coordinate stretching.

**Plan for derivation**

1. Analytically continue solution to the wave equation into $\mathbb{C}$.

2. Coordinate Transfer this extended solution back into $\mathbb{R}^d$.

We first start with a prototypical wave.

$$\nabla\left(a\nabla\cdot u \right) = \frac{1}{b}\left(\frac{\partial^2u}{\partial t^2}\right) = \frac{\ddot{u}}{b} \tag{1}$$

With $u(x,t)$ being a wave of scalar amplitude, $c = \sqrt{ab}$ is the phase velocity for some $a(x)$ and $b(x)$ being properties of the medium.

We can then split this into a pair of coupled first order equations via an auxiliary field, $v(x,t)$ with

$$\dot{u} = b\nabla\cdot v\tag{2} \text{ and } v = a\nabla\cdot u$$

If we have some solution of a wave equation in infinite space, $w(x,t)$ and without loss of generality, a region near the origin that we wish to study - in our case this is $\eta(x)$. In order to apply a PML, we first need to analytically extend our solution to $\mathbb{C}$ within a region near the boundary. We will, for simplicity sake, we will focus on the $\pm x$ direction but it is identical for any other direction.

We are able to make one large simplifying assumption from the statement of the problem from the paper - that being that the far field space is homogeneous, linear and time invariant. Under this assumption, radiation sources must take the form of a superposition of plane waves. [We can see this from the form of our scattered waves $u(x)^s$, but this is a general proof/verification of the PML]

$$w(x,t) = \sum_{k,\omega}W_{k,\omega}e^{i(k\cdot x - \omega t)}\tag{3}$$

For some constant amplitudes $W_{k,\omega}$ we can decompose and simplify this further

$$w(x,t) = W(y,z)e^{i(k\cdot x - \omega t)}\tag{4}$$

With $\frac{\omega}{k}$ being the phase velocity (which can be different from the group velocity, $\frac{d\omega}{dk}$), which we take to be positive as we are moving in the positive direction (in a homogeneous medium, these *in general* have the same sign), so we may take $k$ to be positive.

**Analytical Continuation**

Thankfully, we do not need to construct an analytic continuation into $\mathbb{C}$ for $(4)$ - our wave equation already has the form of $e^{ix}$ which is famously analytic! So all we need to do is evaluate this $(6)$ along a complex contour so that we can evaluate along some complex numbers. We can construct a very simple contour, where for some point Re$(x)$, for example Re$(x) \geq 3$, we start to add linearly growing complex component - i.e. Im$(x) = it\cdot$Re$(x)$, $t > 0$ for $x \geq 3$. from this we can clearly see that when we evaluate $(6)$ along this contour, $(6)$ exponentially decays, so along this contour the wave appears to be passing through an absorbing material. However, the critical thing to notice is that for Re$(x) < 3$ it acts like a reflectionless absorbing material - the PML. Moreover, the analytical solution satisfies the exact same differential equation, now all we need to do is to solve the differential equation along this new complex contour.

**Coordinate Transform back to $\mathbb{R}$**

We can denote $\tilde{x}$ as the complex contour and write it as $\tilde{x} = x + if(x)$ with $f(x)$ being some linear function - this is how we have deformed the contour and hence how fast the wave will decay. Since $x\in\mathbb{R}$, $\tilde{x}$ is a real function of $x$ and so all we need to do in order to transform our coordinates back into real values of $x$ is to do a change of variables. We can do this both quickly and easily by applying the following change of variables: 

$$d\tilde{x} \mapsto \left(1+f'(x)\right)dx$$

This is all that is required as we have assumed homogeneity for far field $x$. for future convenience we denote $f'(x) = \frac{\sigma_x(x)}{\omega}$ - with $\sigma_x(x)$ being some, not necessary continuous, function which describes the contour deformation. This is done to make our attenuation rate not depend on the frequency of the wave. With this notation we can define the PML as a single transformation of our original equation:

$$\frac{\partial}{\partial x} \mapsto \frac{1}{1 + i\frac{\sigma_x(x)}{w}}\frac{\partial}{\partial x}\tag{5}$$

Here, we can see that for the regions where $\sigma_x(x) > 0$ that our solutions of the form (4) will decay exponentially by the contour deformation. While for $\sigma_x(x) = 0$ we get that $\frac{\partial}{\partial x}  = \frac{\partial}{\partial x}$. Thanks to this, we do not get any reflections at the boundary of the PML as it is just an analytic continuation of the function into the complex plane, from $x$ to $\tilde{x}$, and so when $x = \tilde{x}$, our original solution remains the same.

Once we have applied $(5)$ we may then truncate the computational region at some sufficiently large $x$ outside of our domain of interest - such as with a hard wall or a Dirichlet boundary. This is as, thanks to the PML, if a wave is reflected by the boundary it will be attenuated by the PML [^2] a second time before being allowed into the region of interest - so it will have a negligible amplitude. 

**Other Directions**

The final direction we wish to check is $-x$ as all other directions are similar to the derivation to the $x$ direction. This is also fairly simple however, apply $(5)$ but with $\sigma_x(x) > 0$ for sufficiently negative $x$ as well. This still works as, for $-x$, we have $k<0$  and so within our PML region our wave solution still attenuates. Using similar conditions we can now generalise to $\pm y, \pm z$.

As we have done the above solution within the frequency domain, there may issues in generalising this to a time domain, as we chose a transformation of the form of $(5)$. This was done to make out attenuation frequency independent, however expressing $\frac{1}{\omega}$ in the time domain may cause some challenges as we don't have $\omega$ separate in the time-domain since the wave function may superimpose multiple frequencies at once. This is not an issue for us as we are dealing with a fixed$\omega$ - however the more general case may be solved via the use of another auxiliary differential equation. This full derivation can be seen from the notes from Steven G. Johnson[^3], the foundation of my understanding of the PML derivation.

Now we wish to apply this idea of a PML to the forward pass that discussed in a previous section, namely the inhomogeneous equation:

$$Lu:= \left(-\Delta - \frac{\omega^2}{c^2(x)}\right)$$

We define the scatterer, $\eta(x)$ as

$$\eta(x) := \left(\frac{\omega^2}{c^2(x)} - \frac{\omega^2}{c_0^2(x)}\right)$$

which we can see is compactly supported in $\Omega$ as $c^2(x)$ is equal to $c_0^2(x)$ outside $\Omega$.

So we can define the remainder of the operator as such  

$$L_0:= -\Delta - \frac{\omega^2}{c_0^2(r)} = -\Delta - \omega^2 $$

Without loss of generality, we can set $c_0(r)$ to be constant, for convenience we set $c_0(r)\equiv1$. so we finally get  

$$Lu:=\left(L_0-\eta(r)\right)u = f$$ 

We can re-write this into an eaiser form to apply the PML too

$$-\nabla\cdot\left(\nabla u\right) - \frac{\omega^2}{c^2(r)}u  - \eta(r)u = f \tag{6}$$

We can then construct and apply a similar change of variables as the ones seen above

$$x' = \begin{cases}  
x + \frac{ic^2(x)}{\omega}\int_{x_0}^{x}\sigma_{0,x}\left(|x|-x_0\right)^n \, dx &\text{if } |x| \geq x_0 \\  
x &\text{if } |x| < x_0
\end{cases}$$

Similarly for $y$. We will also simplify the notation by denoting the following:

$$\sigma_{\xi}(\xi) = \sigma_{0,\xi}\left(|\xi|-\xi_0\right)^n$$

With $\xi$ being either $x$ or $y$ and $\sigma_{0,\xi}$ being a constant, and n $\in \mathbb{Z}$. This simplifies the definition of the derivative

$$d_x(x)= \begin{cases}  
1 + \frac{ic^2(x)}{\omega}\sigma_{x}(x) &\text{if } |x| \geq x_0 \\  
1 &\text{if } |x| < x_0
\end{cases}$$

And once again similarly for $d_y(y)$

We can clearly see that $\frac{\partial x'}{\partial x} = d_x(x)$ and  $\frac{\partial y'}{\partial y} = d_y(y)$. As per our derivation above we now need to apply a substitution to our original equation in order to incorporate the PML, we can do this by the following:

$$\frac{\partial}{\partial x} \to \frac{\partial}{\partial x'} = \frac{1}{d_x(x)}\frac{\partial}{\partial x}$$

This changes $(6)$ into the following:

$$\frac{1}{d_x(x)}\frac{\partial}{\partial x}\left(\frac{1}{d_x(x)}\frac{\partial u}{\partial x}\right) + \frac{1}{d_y(y)}\frac{\partial}{\partial y}\left(\frac{1}{d_y(y)}\frac{\partial u}{\partial y}\right) - \frac{\omega^2}{c^2(r)}u - \eta(r)u = f $$

We will be writing $d_{\xi}(\xi)$ as $d_{\xi}$ to be more concise. If we group the $d_{\xi}$ terms together we arrive at :

$$\frac{1}{d_x d_y}\left(\frac{1}{\partial x}\left(\frac{d_y}{d_x}\frac{\partial u}{\partial x}\right)+\frac{1}{\partial y}\left(\frac{d_x}{d_y}\frac{\partial u}{\partial y}\right)\right) - \frac{\omega^2}{c^2(r)}u - \eta(r)u = f$$

As $d_x$ is independent of y and $d_y$ is independent of x so we can pass them through the derivatives. Furthermore, defining $\zeta = d_x d_y$ and $A = \begin{pmatrix} \frac{d_y}{d_x} & 0\\ 0 & \frac{d_x}{d_y} \end{pmatrix}$ allows us to radically simplify the expression to:

$$\nabla\cdot\left(A\cdot\nabla\cdotp u\right) -\zeta\left(\frac{\omega^2}{c^2(r)}u + \eta(r)u\right) = \zeta f$$

We can clearly see that setting $\zeta = 1$ and $A = I$ we recover the original equation without the PML. 

[^1]: https://www.sciencedirect.com/science/article/pii/S0021999184711594

[^2]: We may also make the PML of arbitrary size as long as we, correspondingly, make sufficient changes to $\sigma_x(x)$. In general, we make $\sigma_x(x)$ $O(x^2)$ or $O(x^3)$ and the PML about $\frac{\lambda}{2}$ thick. 

[^3]: http://math.mit.edu/~stevenj/18.369/pml.pdf

[^4]: https://liveuclac-my.sharepoint.com/personal/ucahtbe_ucl_ac_uk/Documents/Microsoft%20Teams%20Chat%20Files/monk.pdf