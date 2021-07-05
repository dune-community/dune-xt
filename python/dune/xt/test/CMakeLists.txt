# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2020)
# ~~~

dune_pybindxi_add_module(_test_grid_element grid_element.cc)
dune_pybindxi_add_module(_test_grid_intersection grid_intersection.cc)
dune_pybindxi_add_module(_test_grid_types grid_types.cc)
