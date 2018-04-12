# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016, 2018)
#   Tobias Leibner  (2016 - 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

set(lib_dune_xt_la_sources container/pattern.cc eigen-solver/internal/numpy.cc)

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  list(APPEND lib_dune_xt_la_sources
              container/common.cc
              container/eigen/dense.cc
              container/eigen/sparse.cc
              container/istl.cc
              solver/common.cc
              solver/eigen.cc
              solver/istl.cc)
endif()

find_package(LAPACKE)
dune_library_add_sources(dunextla SOURCES ${lib_dune_xt_la_sources})

add_subdirectory(test EXCLUDE_FROM_ALL)

