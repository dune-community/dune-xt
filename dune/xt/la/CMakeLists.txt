# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2016, 2018 - 2019)
#   Tobias Leibner  (2016 - 2020)
# ~~~

set(_lib_dune_xt_la_sources container/pattern.cc eigen-solver/internal/numpy.cc)

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  list(
    APPEND
    _lib_dune_xt_la_sources
    container/common.cc
    container/eigen/dense.cc
    container/eigen/sparse.cc
    container/istl.cc
    solver/common.cc
    solver/eigen.cc
    solver/istl.cc)
endif()

dune_library_add_sources(dunext SOURCES ${_lib_dune_xt_la_sources})
