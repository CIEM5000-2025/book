---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.2
---

# Beam

::::::{attention}
This page shows a preview of the assignment. Please fork and clone the assignment to work on it locally from [GitHub](https://github.com/CIEM5000-2025/practice-assignments)

After the second workshop, the solution will be added to this preview and to the [GitHub-repository](https://github.com/CIEM5000-2025/practice-assignments)
::::::

```{custom_download_link} ./beam_stripped.ipynb
:text: ".ipynb"
:replace_default: "True"
```

```{custom_download_link} ./beam.md
:text: ".md:myst"
:replace_default: "False"
```

```{custom_download_link} https://github.com/CIEM5000-2025/practice-assignments
:text: "All files practice assignments"
:replace_default: "False"
```

Given is the following beam {cite:p}`additional_Hans`:

```{figure} https://raw.githubusercontent.com/ibcmrocha/public/main/beam.png
:align: center
:width: 400
```

With:
- $l_1 = 5.5$
- $l_2 = 5.0$
- $EI_1 = 5000$
- $EI_2 = 8000$
- $q = 6$
- $F = 40$
- $T = 50$

```{exercise-start} Beam
:label: exercise_beam
:nonumber: true
```

Solve this problem.

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
