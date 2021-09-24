# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2017 - 2018, 2020)
# ~~~

include(Hints)

find_package(BLAS)
if(BLAS_FOUND)
  list(APPEND LAPACKE_LIBRARIES ${BLAS_LIBRARIES})
  list(APPEND LAPACKE_LIBRARIES gfortran)
endif()
message("-- checking for lapacke library")
find_library(LAPACKE_LIBRARY lapacke HINTS ${lib_hints})
if("${LAPACKE_LIBRARY}" MATCHES "LAPACKE_LIBRARY-NOTFOUND")
  message("--   library 'LAPACKE' not found, make sure you have both LAPACK and LAPACKE installed")
else("${LAPACKE_LIBRARY}" MATCHES "LAPACKE_LIBRARY-NOTFOUND")
  message("--   found LAPACKE library")
  set(LAPACKE_LIBRARIES "${LAPACKE_LIBRARY}")
endif("${LAPACKE_LIBRARY}" MATCHES "LAPACKE_LIBRARY-NOTFOUND")

message("-- checking for lapacke.h header")
find_path(LAPACKE_INCLUDE_DIRS lapacke.h HINTS ${include_hints})
if("${LAPACKE_INCLUDE_DIRS}" MATCHES "LAPACKE_INCLUDE_DIRS-NOTFOUND")
  message("--   lapacke.h header not found")
else("${LAPACKE_INCLUDE_DIRS}" MATCHES "LAPACKE_INCLUDE_DIRS-NOTFOUND")
  message("--   found lapacke.h header")
  include_sys_dir("${LAPACKE_INCLUDE_DIRS}")
endif("${LAPACKE_INCLUDE_DIRS}" MATCHES "LAPACKE_INCLUDE_DIRS-NOTFOUND")

if(NOT DEFINED LAPACKE_FOUND)
  set(LAPACKE_FOUND 0)
endif(NOT DEFINED LAPACKE_FOUND)
if(NOT "${LAPACKE_LIBRARY}" MATCHES "LAPACKE_LIBRARY-NOTFOUND")
  if(NOT "${LAPACKE_INCLUDE_DIRS}" MATCHES "LAPACKE_INCLUDE_DIRS-NOTFOUND")
    set(LAPACKE_FOUND 1)
  endif(NOT "${LAPACKE_INCLUDE_DIRS}" MATCHES "LAPACKE_INCLUDE_DIRS-NOTFOUND")
endif(NOT "${LAPACKE_LIBRARY}" MATCHES "LAPACKE_LIBRARY-NOTFOUND")

set(HAVE_LAPACKE ${LAPACKE_FOUND})

# register all LAPACKE related flags
if(LAPACKE_FOUND)
  dune_register_package_flags(LIBRARIES "${LAPACKE_LIBRARIES}" INCLUDE_DIRS "${LAPACKE_INCLUDE_DIRS}")
endif()
