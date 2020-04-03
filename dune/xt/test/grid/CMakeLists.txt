# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2016 - 2019)
#   Tobias Leibner  (2016, 2018, 2020)
# ~~~

enable_testing()

dxt_exclude_from_headercheck(dd_subdomains_cube.hh)

# we want to discard one of the tests via cmake_guard if UG is parallel, so define a suitable variable here
if("${UG_PARALLEL}" STREQUAL "yes")
  set(UG_IS_NOT_PARALLEL FALSE)
else()
  set(UG_IS_NOT_PARALLEL TRUE)
endif()
