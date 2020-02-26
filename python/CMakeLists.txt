# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2018)
# ~~~

add_subdirectory(dune)
add_subdirectory(scripts)
add_subdirectory(wrapper)
if (DUNE_XT_WITH_PYTHON_BINDINGS)
  add_custom_target(bindings DEPENDS common_bindings grid_bindings functions_bindings la_bindings)
endif()
dxt_add_python_tests()
