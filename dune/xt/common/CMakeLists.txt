# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017 - 2019)
#   René Fritze     (2017 - 2019)
#   Tobias Leibner  (2017 - 2020)
# ~~~

set(lib_dune_xt_common_sources
    bisect.cc
    cblas.cc
    color.cc
    configuration.cc
    convergence-study.cc
    exceptions.cc
    filesystem.cc
    fix-ambiguous-std-math-overloads.cc
    lapacke.cc
    localization-study.cc
    logging.cc
    logstreams.cc
    math.cc
    memory.cc
    misc.cc
    mkl.cc
    parallel/communicator.cc
    parallel/helper.cc
    parallel/mpi_comm_wrapper.cc
    parallel/threadmanager.cc
    parameter.cc
    python.cc
    signals.cc
    string.cc
    timedlogging.cc
    timings.cc)

dune_library_add_sources(dunext SOURCES ${lib_dune_xt_common_sources})

add_analyze(${lib_dune_xt_common_sources})
