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

# Apply

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)
::::::

::::::{versionadded} v1.2.0 After workshop 2
Solutions in text and downloads 
::::::

In this notebook you will solve a 2-element frame at the end of the notebook.

+++

Our matrix method implementation is now completely stored in a local package, consisting of three classes.

+++

```{custom_download_link} ./Workshop_2_Apply_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./Workshop_2_Apply_stripped_sol.ipynb
:text: ".ipynb solution"
:replace_default: "False"
```

```{custom_download_link} ./Workshop_2_Apply.md
:text: ".md:myst"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments/tree/solution_workshop_2
:text: "All files practice assignments with solutions workshop 2"
:replace_default: "False"
```


## Two-element frame

+++

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/twoelemframe.png
:align: center
:width: 400
```


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
```

```{code-cell} ipython3
:tags: [thebe-remove-input-init]

import matplotlib as plt
import numpy as np
sys.path.insert(1, '/matrixmethod_solution')
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [remove-cell]

import matplotlib as plt
import numpy as np
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [disable-execution-cell]

import numpy as np
import matplotlib as plt
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
![](https://raw.githubusercontent.com/ibcmrocha/public/main/moments_local.svg)

- in global coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/moments_global.svg)

And the following displacements:
- in local coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/displacements_local.svg)

- in global coordinate system:
![](https://raw.githubusercontent.com/ibcmrocha/public/main/displacements_global.svg)
::::::

+++

```{solution-start} exercise_ws_2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-init]

EI = 1500
EA = 1000
q = 9
L = 5
phibar = 0.15

mm.Node.clear()
mm.Element.clear()

nodes = []

nodes.append(mm.Node(0,0))
nodes.append(mm.Node(L,-L))
nodes.append(mm.Node(2*L,0))

elems = []

elems.append(mm.Element(nodes[0], nodes[1]))
elems.append(mm.Element(nodes[1], nodes[2]))

section = {}
section['EI'] = EI
section['EA'] = EA

for elem in elems:
    elem.set_section(section)
    print(elem)

con = mm.Constrainer()

con.fix_dof (nodes[0], 0)
con.fix_dof (nodes[0], 1)
con.fix_dof (nodes[2], 0)
con.fix_dof (nodes[2], 1)
con.fix_dof (nodes[2], 2, phibar)

elems[0].add_distributed_load([0,q])
elems[1].add_distributed_load([0,2*q])

print(con)

global_k = np.zeros ((3*len(nodes), 3*len(nodes)))
global_f = np.zeros (3*len(nodes))

for e in elems:
    elmat = e.stiffness()
    idofs = e.global_dofs()
    
    global_k[np.ix_(idofs,idofs)] += elmat

for n in nodes:
    global_f[n.dofs] += n.p

Kc, Fc = con.constrain ( global_k, global_f )
u_free = np.matmul ( np.linalg.inv(Kc), Fc )
print(u_free)

print(con.support_reactions(global_k,u_free,global_f))
```

```{code-cell} ipython3
:tags: [thebe-init]

for elem in elems:
    u_elem = con.full_disp(u_free)[elem.global_dofs()]
    elem.plot_displaced(u_elem,num_points=51,global_c=False,scale=1)
```

```{code-cell} ipython3
:tags: [thebe-init]

for elem in elems:
    u_elem = con.full_disp(u_free)[elem.global_dofs()]
    elem.plot_displaced(u_elem,num_points=51,global_c=True,scale=1)
```

```{code-cell} ipython3
:tags: [thebe-init]

for elem in elems:
    u_elem = con.full_disp(u_free)[elem.global_dofs()]
    elem.plot_moment_diagram(u_elem,num_points=20,global_c=False)
```

```{code-cell} ipython3
:tags: [thebe-init]

for elem in elems:
    u_elem = con.full_disp(u_free)[elem.global_dofs()]
    elem.plot_moment_diagram(u_elem,num_points=20,global_c=True,scale=0.05)
```

```{solution-end}
```

```{exercise-start} 2
:label: exercise_ws_2_2
:nonumber: true

You've now implemented the matrix method yourself! However, there are many implementation our there. Therefore, let's implement another package which you'll use during the next weeks when using the matrix method for dynamics: PyDSMSM.

The library is provided here: **to be provided**

```


```{exercise-end}
```

```{solution-start} exercise_ws_2_2
:class: dropdown
```

```{code-cell} ipython3
:tags: [thebe-remove-input-init]

import matplotlib as plt
import numpy as np
sys.path.insert(1, '/pydynsm')
import numpy.linalg as LA
import scipy.optimize as opt
import pydynsm as PDM
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
:tags: [disable-execution-cell]

import numpy as np
import matplotlib as plt
import pydynsm as PDM
import numpy.linalg as LA
import scipy.optimize as opt
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell} ipython3
# Initialize the Assembler

Assembler = PDM.Assembler  #Import the Assembler class from the PDM module

s1 = Assembler('1D EB',analysis_type='new')
```

```{code-cell} ipython3
E = 1000
I = 1.5
A = 1
q = 9
L = 5
phibar = 0.15
```

```{code-cell} ipython3
node1 = s1.CreateNode(0,0)
node2 = s1.CreateNode(L,L)
node3 = s1.CreateNode(2*L,0)
elem1 = s1.CreateElement([node1,node2])
elem2 = s1.CreateElement([node2,node3])
```

```{code-cell} ipython3
s1.PlotStructure(plot_elements=True)
```

```{code-cell} ipython3
node1.fix_node('x','z')
node1.apply_dof_change_to_elements('phi_y',phibar)
node3.fix_node('x','z')
```

```{code-cell} ipython3
# Can this be defined without the dynamic properties?
elem1.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'Ib':I})
elem2.SetSection('EulerBernoulli Beam', {'E': E, 'A':A, 'Ib':I})
```

```{code-cell} ipython3
# How to add a distributed load to an element in the local coordinate system?
q_r = lambda ?
q_b = lambda ?
elem.AddDistributedLoad(x=q_r, z=q_b)
```

```{code-cell} ipython3
s1.run_connectivity()
```

```{code-cell} ipython3
# could the package be defined to a static case when omega is not provided?
K_global = s1.GlobalStiffness()
F_global = s1.GlobalForce()
```

```{code-cell} ipython3
Kc_global = s1.GlobalConstrainedStiffness()
Fc_global = s1.GlobalConstrainedForce()
```

```{code-cell} ipython3
u_free = s1.SolveUfree(Kc_global, Fc_global)
```

```{code-cell} ipython3
u_elem = s1.FullDisplacement(u_free)
print(f'u_elem = \n{u_elem}\n')
```

```{code-cell} ipython3
disp = s1.ElementDisplacements(u_elem)
```

```{code-cell} ipython3
s1.PlotElementDisplacements(disp,scale=1.0)
```

```{solution-end}
```
