---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
kernelspec:
  display_name: base
  language: python
  name: python3
---

#Apply

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)

After the workshop, the solution will be added to this preview and to the [GitHub-repository](https://github.com/CIEM5000-2025/practice-assignments)
::::::

In this notebook you will solve a 2-element frame at the end of the notebook.

+++

Our matrix method implementation is now completely stored in a local package, consisting of three classes.

+++

```{custom_download_link} ./Workshop_2_Apply_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./Workshop_2_Apply.md
:text: ".md:myst"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

## Two-element frame

+++

<figure>
  <IMG SRC="https://raw.githubusercontent.com/ibcmrocha/public/main/twoelemframe.png" WIDTH=300 ALIGN="center">
</figure>

With:
- $EI = 1500$
- $EA = 1000$
- $q = 9$
- $L = 5$
- $\bar\varphi = 0.15$

```{exercise-start}
:label: exercise_ws_2
:nonumber: true

The final example for the workshops is the two-element frame above. Here you should make use of all the new code you implemented:
    
- Set up the problem and compute a solution for `u_free`. Remember to consider the prescribed horizontal displacement $\bar{u}$ at the right end of the structure.
- Compute and plot bending moment lines for both elements (in the local and global coordinate systems)
- Compute reactions at both supports

```{code-cell} ipython3
import numpy as np
import matrixmethod as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
#YOUR CODE HERE
```

```{code-cell} ipython3
for elem in elems:
    u_elem = con.full_disp(#YOUR CODE HERE)[#YOUR CODE HERE.global_dofs()]
    elem.plot_displaced #YOUR CODE HERE
```

```{exercise-end}
```

::::::{hint}
:class: dropdown

For the given parameter values, if your implementation is fully correct, you should get the following nodal displacements and support reactions:
$$
\mathbf{u}_\mathrm{free} = \left[-0.09274451, -0.13310939,  0.51159348, -0.01644455\right]
$$

$$
\mathbf{f}_\mathrm{cons} = \left[27.35024439, -63.82451092,  17.64975561, -71.17548908, -36.75489076\right]
$$

You should also get the following moment lines for the two elements:

- in local coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/moments_local.png)

- in global coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/moments_global.png)

And the following displacements:
- in local coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/displacements_local.png)

- in global coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/displacements_global.png)

::::::
