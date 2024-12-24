---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.6
---

# Recap displacement method

In the previous [chapter](./recap.ipynb) you've seen how solving for integration constants because a labour-intensive process for more complicated structures. A way of circumventing that is solving for nodal displacements! You might have seen that before when solving statically indeterminate structures using the displacement method!

The displacement method for statically indeterminate structures works by defining a single nodal displacement of the equivalent statically determinate structure which defines the displacement of the full structure. The nodal displacement can be solved for by equilibrium relations of the external forces and the force corresponding to the displacement. Solving this equation allows you to find the full displacement- and force distribution.

Let's look at an statically indeterminate example

```{figure} extension2fieldsdisp.svg
:name: extension2fieldsq
:align: center

Statically indeterminate extension bar
```

A statically determinate equivalent structure is i.e. the same structure with the middle connection replaced by a displacement $u_2$ and its corresponding reaction force $F_1$ and $F_2$. This leads to two parts. If $F_1 = F_2$, the structure is equivalent to the statically indeterminate structure. 

```{figure} extension2fieldsdispdet.svg
:name: extension2fieldispdet
:align: center

Equivalent statically determinate extension bar if $F_1 = F_2$
```

## Displacements of all parts statically equivalent structure
Due to the nodal load, a constant section force is present in both fields. Using the constitutive equations $\Delta \ell = \cfrac{N \ell}{EA}$, the corresponding displacement $u_2$ can be found. For the load $q$, the displacement can be calculated using the kinematic, constitutive and kinematic relations. This leads to $u_2$ as a function of $F_1$ and $F_2$:

- $u_2 = \cfrac{\ell_1}{EA_1}F_1  + \cfrac{\ell_1^2}{2EA_1}q$
- $u_2 = -\cfrac{\ell_2}{EA_1}F_2  + \cfrac{\ell_2^2}{2EA_2}q$

These relations can be rewritten as $F_1$ and $F_2$ in terms of $u_2$ so that the force equilibrium can be solved for:

- $F_1 = \cfrac{EA_1}{\ell_1} u_2 - \cfrac{\ell_1^2}{2EA_1}q$
- $F_2 = -\cfrac{EA_2}{\ell_2} u_2 - \cfrac{\ell_2^2}{2EA_2}q$

## Solve for displacements

Using $F_1 = F_2$ the displacement $u_2$ can now be solved for:

$$ \cfrac{EA_1}{\ell_1} u_2 - \cfrac{\ell_1^2}{2EA_1}q = -\cfrac{EA_2}{\ell_2} u_2 - \cfrac{\ell_2^2}{2EA_2}q \\
\left(\cfrac{EA_1}{\ell_1} + \cfrac{EA_2}{\ell_2}\right) u_2 = \cfrac{\ell_1^2}{2EA_1}q - \cfrac{\ell_2^2}{2EA_2}q \\
u_2 = \cfrac{\cfrac{\ell_1^2}{2EA_1}q - \cfrac{\ell_2^2}{2EA_2}q }{\cfrac{EA_1}{\ell_1} + \cfrac{EA_2}{\ell_2}}$$

## Equivalence with matrix method

The matrix method applies exactly the same principle as the displacement method; both are solving force equilibrium of nodal forces to find nodal displacements.

However, the displacement method becomes difficult to apply if multiple nodal displacements are taken into account, as nodal forces have effect on multiple nodal displacements. Furthermore, the calculation of the displacements of each part can become tedious because they're problem-dependent and external forces have to be taken into account in the full derivation.

The matrix method addresses these issues by splitting the structure in mostly identical elements for which the force-displacement relations for all potential nodal displacements are evaluated once and can be reused over and over again. The same approach is taken for external forces, of which the resulting relations can be added afterwards. Finally, the calculations are structured in matrices to allow for easy implementation in software.

The similarities and differences are shown in the table below.

:::{table} Equivalence displacement method and matrix method
:widths: auto
:align: center

|Displacement method|Matrix method|
|:-:|:-:|
|Convert structure in statically determinate parts|Convert structure in mostly identical elements|
|Evaluate one or a few nodal displacements for each parts|Evaluate all nodal displacement for default elements a priori|
|Solve nodal equilibrium $\left(\cfrac{EA_1}{\ell_1} + \cfrac{EA_2}{\ell_2}\right) u_2 = \cfrac{\ell_1^2}{2EA_1}q - \cfrac{\ell_2^2}{2EA_2}q $|Solve nodal equilibrium in matrix form $\mathbf{K}\mathbf{u}=\mathbf{f}$|

:::