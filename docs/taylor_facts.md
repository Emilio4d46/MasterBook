Key Facts and Theorems
=======================

To get a good foundation in the mathematics behind wave scattering and their inverse problems, we will be reading through the book "Partial Differential Equations 2 by Michael E. Taylor" [^1] , as well as notes  from a course on Linear Partial Differential Equations taught by Dr I Kamotski during term 1, 2019 at UCL.

 We will most notably follow chapter 9, Scattering by Objects, and creating a resource of important and useful theorems and facts in order to build up a base of useful knowledge.

We will first lay down the fundamental equations relating to scattering.

Let $K \subset \mathbb{R}^3$ with a smooth boundary and $\Omega := K^c$ being a connected complement. Let $f(x) \in H^s(K)$ be a given function from the Sobolev Space of positive index $s$ and let $k>0$, we wish to solve the following:

$$(\Delta + k^2)v = 0 \text{ on }\Omega\tag{1}$$ 

$$v = f \text{ on }\partial K\tag{2}$$ 

With the typical Somerfield radiation boundary condition applied

$$\left.\lvert{rv(x) \leq C}\right.\rvert \text{ , }\lim_{r\to\infty}r\left(\frac{\partial}{\partial r} - ik\right)v=0 \tag{3}$$

With $r = |x|$, with $v$ satisfying 

$$v(x)=\int_{\partial K} \left[f(y)\frac{\partial g}{\partial v_y}(x,t,k) - g(x,y,k)\frac{\partial v}{\partial v}(y)\right]dS(y)\tag{4}$$

For $x \in \Omega$ and $g(x,y,k)$ being the Green Function.

**Fact 1:** Given $k>0$, if $v$ satisfies $(1)-(3)$ with $f=0$, then $v=0$

**Fact 2:** If $v$ satisfies $(\Delta - k^2)v=0$ for $|x|\geq  R_0$ and

$$v(x)=\int_{S_R} |v|^2dS \to 0 \text{, as } R \to \infty$$

holds, then for $v(x)=0$ for $|x|\geq R_0$. With $S_R$ being the sphere of radius $R$.

**Theorem 3:** For $\epsilon > 0$, if we perturbate equation $(1)$ by $k \to k + i\epsilon$, with $v_{\epsilon}$ being the unique solution to the perturbated equation and $f_{\epsilon}$ being the boundary condition of the perturbated equation.
If $s \geq \frac{3}{2}$, and suppose that as $\epsilon \to 0$

$$f_{\epsilon}\to f\text{,     in } H^s(\partial K)$$

Then we have the unique limit[^2]

$$v_{\epsilon}\to v$$





[^1]: Taylor, M. (2011) Partial Differential Equations II, 2nd ed. New York: Springer

[^2]: Is this limit contained in the Sobolev Space, $H^s(\partial K)$?