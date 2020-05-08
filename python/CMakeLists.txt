# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2018, 2020)
# ~~~

# Currently, the python libraries do not have access to the config.h variables, see
# https://github.com/dune-community/dune-xt-common/issues/120
# The following silences all the resulting -Wundef warnings in this directory
add_compile_options(-Wno-undef)
add_subdirectory(dune)
add_subdirectory(scripts)
add_subdirectory(wrapper)
dxt_add_python_tests()
