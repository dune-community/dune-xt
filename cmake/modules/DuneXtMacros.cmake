# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013 - 2014, 2016 - 2017, 2020)
#   René Fritze     (2013 - 2016, 2018 - 2020)
#   Sven Kaulmann   (2014)
#   Tobias Leibner  (2016, 2018 - 2021)
# ~~~

# enables "IN_LIST operator
cmake_policy(SET CMP0057 NEW)

# For some reason, the minimum required version is set to 2.8.3 by the find_package(Vc ...) call in
# DuneCommonMacros.cmake in dune-common. This causes some warnings, so we reset it here.
cmake_minimum_required(VERSION 3.8)

include(XtCompilerSupport)
include(XtTooling)
include(DuneXTHints)

set(DXT_DONT_LINK_PYTHON_LIB
    ${DXT_DONT_LINK_PYTHON_LIB}
    CACHE STRING "wheelbuilders want to set this to 1")

# library checks  #########################################################################
find_package(PkgConfig)

set(DS_REQUIRED_BOOST_LIBS
    atomic
    chrono
    date_time
    filesystem
    system
    thread
    timer)
set(_boost_root_hints "$ENV{BOOST_ROOT}" ${root_hints})
# check if any hints are provided by user
if(DEFINED BOOST_ROOT
   OR DEFINED BOOOST_INCLUDEDIR
   OR DEFINED BOOST_LIBRARYDIR)
  find_package(Boost 1.48.0 REQUIRED COMPONENTS ${DS_REQUIRED_BOOST_LIBS})
else(
  DEFINED BOOST_ROOT
  OR DEFINED BOOOST_INCLUDEDIR
  OR DEFINED BOOST_LIBRARYDIR)
  # FindBoost can only take a single hint directory from BOOST_ROOT, so we loop over all hints
  foreach(boost_root_hint ${_boost_root_hints})
    set(BOOST_ROOT ${boost_root_hint})
    find_package(Boost 1.48.0 COMPONENTS ${DS_REQUIRED_BOOST_LIBS})
    if(${Boost_FOUND})
      break()
    endif()
  endforeach(boost_root_hint ${_boost_root_hints})
  # check for Boost again with REQUIRED keyword to make boost mandatory
  find_package(Boost 1.48.0 REQUIRED COMPONENTS ${DS_REQUIRED_BOOST_LIBS})
endif(
  DEFINED BOOST_ROOT
  OR DEFINED BOOOST_INCLUDEDIR
  OR DEFINED BOOST_LIBRARYDIR)
if(${Boost_INCLUDE_DIRS})
  dune_register_package_flags(INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
endif(${Boost_INCLUDE_DIRS})
# if imported targets are available, use them
if(TARGET Boost::headers)
  dune_register_package_flags(LIBRARIES Boost::headers)
endif(TARGET Boost::headers)
foreach(boost_lib ${DS_REQUIRED_BOOST_LIBS})
  set(_boost_lib "")
  string(TOUPPER "${boost_lib}" _boost_lib)
  if(TARGET Boost::${boost_lib})
    dune_register_package_flags(LIBRARIES Boost::${boost_lib})
  else(TARGET Boost::${boost_lib})
    dune_register_package_flags(LIBRARIES ${Boost_${_boost_lib}_LIBRARY})
  endif(TARGET Boost::${boost_lib})
endforeach(boost_lib ${DS_REQUIRED_BOOST_LIBS})

find_package(Eigen3 3.2.0)
if(EIGEN3_FOUND)
  dune_register_package_flags(INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} COMPILE_DEFINITIONS "ENABLE_EIGEN=1")
  set(HAVE_EIGEN 1)
else(EIGEN3_FOUND)
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_EIGEN=0")
  set(HAVE_EIGEN 0)
endif(EIGEN3_FOUND)

find_package(MKL)
if(MKL_FOUND)
  set(HAVE_LAPACKE 0)
else(MKL_FOUND)
  find_package(LAPACKE)
endif(MKL_FOUND)

find_package(Spe10Data)

# intel mic and likwid don't mix
if(NOT CMAKE_SYSTEM_PROCESSOR STREQUAL "k1om")
  include(FindLIKWID)
  find_package(LIKWID)
  if(LIKWID_FOUND)
    dune_register_package_flags(INCLUDE_DIRS ${LIKWID_INCLUDE_DIR} LIBRARIES ${LIKWID_LIBRARY})
  endif(LIKWID_FOUND)
else()
  set(HAVE_LIKWID 0)
endif()

include(DuneTBB)

if(HAVE_MPI)
  include(FindMPI4PY)
  if(MPI4PY_FOUND)

  else()
    execute_process(
      COMMAND ${CMAKE_BINARY_DIR}/run-in-dune-env pip install mpi4py
      ERROR_VARIABLE shell_error
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    myfind_mpi4py()
  endif()
  if(MPI4PY_FOUND)
    # this only works in dependent modules
    dune_register_package_flags(INCLUDE_DIRS "${MPI4PY_INCLUDE_DIR}")
    # this only works in dune-xt itself
    include_directories(SYSTEM "${MPI4PY_INCLUDE_DIR}" ${PYTHON_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR kaput)
  endif()
endif()
# end library checks  #####################################################################

# misc vars  #########################################################################

set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(DS_MAX_MIC_THREADS
    120
    CACHE STRING "")
set(DUNE_XT_COMMON_TEST_DIR ${dune-xt_SOURCE_DIR}/dune/xt/common/test)
set(ENABLE_PERFMON
    0
    CACHE STRING "enable likwid performance monitoring API usage")
if(NOT DS_HEADERCHECK_DISABLE)
  set(ENABLE_HEADERCHECK 1)
endif(NOT DS_HEADERCHECK_DISABLE)
set(DXT_TEST_TIMEOUT
    180
    CACHE STRING "per-test timeout in seconds")
set(DXT_TEST_PROCS
    1
    CACHE STRING "run N tests in parallel")

set(DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS TRUE)

include(DunePybindxiUtils)
