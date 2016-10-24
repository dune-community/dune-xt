```
This file is part of the dune-pybindx1 project:
  https://github.com/dune-community/dune-pybindx1
The copyright lies with the authors of this file (see below).
License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
Authors:
  Felix Schindler (2016)
```

This module wraps [pybind11](https://github.com/pybind/pybind11) as a
DUNE module, to facilitate the build process with dunecontrol.

To create Python bindings you should

* define the cmake variable `DUNE_PYBINDXI_PYTHON_VERSION` (optional)
* `include(DunePybindxiUtils)` in your CMakeLists.txt
* call `dune_pybindxi_add_module(module)`, where module is one or several C++
  source files

In your sources, you may include the pybind11 header `pybind11/foo.h` as
```
#include <dune/pybindxi/foo.h>
```

For further documentation on pybind1 and Python bindings,
[read the docs](http://pybind11.readthedocs.io/en/latest/).

