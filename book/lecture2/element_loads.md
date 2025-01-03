# Element loads

As the matrix method is a discrete approach, nodal loads were treated with ease. However, what to do with continuous loads or loads which are not applied at the nodes?

::::::{topic} Learning objective
You'll look into how to model element loads using differential equations and conservation of work and how these are combined in the matrix formulation.
::::::

## Force-displacement relations using differential equations
As [before](../lecture1/single_element.md), we can derive the force-displacement relations of a single extension element. However, now let's include the loads, for example a continuous load $q$

```{figure} extensionelementq.svg
:name: extensionelementq
:align: center

Single extension element with distributed load $q$
```

The same approach is used as in [](../lecture1/recap.ipynb). This results in:

- $C_1 = \class{cB}{\cfrac{q\ell}{2EA}}+\cfrac{u_2-u_1}{\ell}$
- $C_2 = u_1$

The continuous distributions for the displacement and section force can be evaluated too:
- $u(x) = \class{cB}{\cfrac{q}{2EA}\left(\ell x - x^2\right)} + u_1\left(1-\cfrac{x}{\ell}\right) + u_2\cfrac{x}{\ell}$
- $N(x) = \class{cB}{\cfrac{q}{2}\left(\ell-2x \right)} -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell}u_2$

These are extended results in comparison to [before](../lecture1/single_element.md)

## Combine elements
As before, we can glue elements together by applying nodal equilibrium:

This leads to:
- $ F_1 = - N = \cfrac{EA}{\ell}u_1-\cfrac{EA}{\ell}u_2 \class{cB}{ -\cfrac{q\ell}{2}}$
- $ F_2 = N = -\cfrac{EA}{\ell}u_1+\cfrac{EA}{\ell} u_2 \class{cB}{-\cfrac{q\ell}{2}}$

or in matrix notation:

$$\cfrac{EA}{\ell}\begin{bmatrix} 1&-1\\-1&1 \end{bmatrix}\begin{bmatrix}u_1\\u_2\end{bmatrix} \class{cB}{- \begin{bmatrix} \cfrac{q\ell}{2}\\ \cfrac{q\ell}{2}\end{bmatrix}}= \begin{bmatrix}F_1\\F_2\end{bmatrix}$$

Effectively, we converted the continuous load to an equivalent nodal load.

This influences the nodal equilibrium (in the global coordinate system) too:

$$\begin{align} -\sum_e\mathbf{f}^e + \mathbf{f}_\text{nodal}& = \mathbf{0}\\
-\sum_e\left(\mathbf{K}^e\mathbf{u}^e \class{cA}{-\mathbf{f}_\mathrm{eq}^e}\right) + \mathbf{f}_\mathrm{nodal}& = \mathbf{0}\\
\sum_e\mathbf{f}^e& = \mathbf{f}_\text{nodal} + \class{cA}{\sum_e\mathbf{f}_\mathrm{eq}^e} \end{align}$$

This means that we can calculate all the equivalent nodal loads separately and add them to the nodal loads. Please note that all these forces should be in the global coordinate system: $\class{cA}{\mathbf{f}_\mathrm{eq}} = \mathbf{T}^\mathrm{T}\class{cB}{\bar{\mathbf{f}}_\mathrm{eq}}$

## Force-displacement relations using conservation of work

### Point load
Let's consider another example in which the application of the differential equations are not trivial: we'll introduce a discontinuous force in the form of a point load halfway the element. We'll compare this with the same element loaded by vertical forces and bending moments at the ends:

```{figure} eqwork3.svg
:name: eqwork3
:align: center

Single extension element with point load $P$ halfway in comparison to element loaded by vertical forces and bending moments at the ends 
```

The conservation of work states that the work done by the force $P$ should be equal to the force done by the forces $\class{cA}{F_1^{eq}}$, $\class{cE}{F_2^{eq}}$, $\class{cB}{T_1^{eq}}$ and $\class{cI}{T_2^{eq}}$: $W_P$ = $W_{eq}$

To ease the calculation, we'll split the displacement in four cubic shape functions (cubic is justified because with $q=0$ a cubic displacement function is found). Each shape function will have a displacement of $1$ in direction of each of the four edge forces ($\class{cA}{w_1}$, $\class{cE}{\varphi_1}$, $\class{cB}{w_2}$ and $\class{cI}{\varphi_2}$):

```{figure} shape_functions.svg
:name: shape_functions
:align: center

Four shape functions
```

The shape functions have the following function:

$$
w(x) = 
	    \class{cA}{\underbrace{\left(\cfrac{2x^3}{\ell^3}-\cfrac{3x^2}{\ell^2}+1\right)}_{s_1}w_1} +
	    \class{cE}{\underbrace{\left(-\cfrac{x^3}{\ell^2}+\cfrac{2x^2}{\ell}-x\right)}_{s_2}\varphi_1} +
	    \class{cB}{\underbrace{\left(-\cfrac{2x^3}{\ell^3}+\cfrac{3x^2}{\ell^2}\right)}_{s_3}w_2} +
	    \class{cI}{\underbrace{\left(-\cfrac{x^3}{\ell^2}+\cfrac{x^2}{\ell}\right)}_{s_4}\varphi_2}
$$

Consequently, the work $W_{eq}$ can be splitted in four independent terms easing the calculation, which can be compared to the same terms in $W_P$.

The work performed by the edge forces equals:

$$ W_{eq} = 
	    \class{cA}{F^\mrm{eq}_1w_1} +
	    \class{cE}{T^\mrm{eq}_1\varphi_1} +
	    \class{cB}{F^\mrm{eq}_2w_2} + 
	    \class{cI}{T^\mrm{eq}_2\varphi_2}$$

The work performed by $P$ (under the same displacement) is:

$$W_P = P\ w\left(\frac{\ell}{2}\right) = 
	  \class{cA}{P\ s_1\left(\frac{\ell}{2}\right) \ w_1} +
	  \class{cE}{P\ s_2\left(\frac{\ell}{2}\right) \ \varphi_1} +
	  \class{cB}{P\ s_3\left(\frac{\ell}{2}\right) \ w_2} +
	  \class{cI}{P\ s_4\left(\frac{\ell}{2}\right) \ \varphi_2}$$

Enforcing $W_F = W_q$ and isolating terms gives:

$$
\mathbf{f}_\mrm{eq}
	    =
	    \begin{bmatrix}
	      \class{cA}{F^\mrm{eq}_1}\\[10pt]
	      \class{cE}{T^\mrm{eq}_1}\\[10pt]
	      \class{cB}{F^\mrm{eq}_2}\\[10pt]
	      \class{cI}{T^\mrm{eq}_2}\\
	    \end{bmatrix}
	    =
	    \begin{bmatrix}
	      \class{cA}{\cfrac{P}{2}}     \\
	      \class{cE}{-\cfrac{P\ell}{8}}\\
	      \class{cB}{\cfrac{P}{2}}     \\
	      \class{cI}{\cfrac{P\ell}{8}} \\
	    \end{bmatrix}
$$

In the local coordinate system.

### Distributed load
The same can be done for a distributed load:

$$W_q = \int_\ell q \ w(x)\mrm{d}x = 
	  \class{cA}{\int_\ell q \ s_1(x)\mrm{d}x \ w_1} +
	  \class{cE}{\int_\ell q \ s_2(x)\mrm{d}x \ \varphi_1} +
	  \class{cB}{\int_\ell q \ s_3(x)\mrm{d}x \ w_2} +
	  \class{cI}{\int_\ell q \ s_4(x)\mrm{d}x \ \varphi_2}$$

Leading to:

$$
\mathbf{f}_\mrm{eq}
	    =
	    \begin{bmatrix}
	      \class{cA}{F^\mrm{eq}_1}\\[10pt]
	      \class{cE}{T^\mrm{eq}_1}\\[10pt]
	      \class{cB}{F^\mrm{eq}_2}\\[10pt]
	      \class{cI}{T^\mrm{eq}_2}\\
	    \end{bmatrix}
	    =
	    \begin{bmatrix}
	      \class{cA}{\cfrac{q\ell}{2}}     \\
	      \class{cE}{-\cfrac{q\ell^2}{12}} \\
	      \class{cB}{\cfrac{q\ell}{2}}   \\
	      \class{cI}{\cfrac{q\ell^2}{12}} \\
	    \end{bmatrix}
$$
