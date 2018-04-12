# ~~~
# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
# ~~~

dune_pybindxi_add_module(_grid EXCLUDE_FROM_ALL bindings.cc boundaryinfo.bindings.cc)
add_dune_alberta_flags(GRIDDIM 2 _grid)
