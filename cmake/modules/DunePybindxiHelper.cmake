# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2018, 2020)
# ~~~

function(dune_pybindxi_add_helper_lib target_name)
  message(WARNING "dune_pybindxi_add_helper_lib is deprecated. Add source to bindings lib directly")
  add_library(${target_name} EXCLUDE_FROM_ALL ${ARGN})
  target_include_directories(${target_name} PRIVATE ${PYBIND11_INCLUDE_DIR} PRIVATE ${PYTHON_INCLUDE_DIRS})
endfunction()
