---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
---

# Hinged beam

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)

After the second workshop, the solution will be added to this preview and to the [GitHub-repository](https://github.com/CIEM5000-2025/practice-assignments)
::::::

```{custom_download_link} ./hinged_beam_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./hinged_beam.md
:text: ".md:myst"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

Given is the following beam {cite:p}`additional_Hans`:

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/hingedbeam.png
:align: center
:width: 400
```

With:
- $EI = 8000$


```{exercise-start}
:label: exercise_hinged_beam
:nonumber: true
```

Solve this problem. How do you deal with the hinges?

```{code-cell}
:tags: [thebe-remove-input-init]

import matplotlib as plt
import numpy as np
sys.path.insert(1, '/matrixmethod_solution')
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell}
:tags: [remove-cell]

import matplotlib as plt
import numpy as np
import matrixmethod_solution as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell}
:tags: [disable-execution-cell]

import numpy as np
import matplotlib as plt
import matrixmethod as mm
%config InlineBackend.figure_formats = ['svg']
```

```{code-cell}
#YOUR_CODE_HERE
```

```{exercise-end}
```
