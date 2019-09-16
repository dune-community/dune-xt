# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018)
#
# dune-xt-common developers and contributors. All rights reserved. License: Dual licensed as BSD 2-Clause License
# (http://opensource.org/licenses/BSD-2-Clause) Authors: Felix Schindler (2017) Rene Milk       (2018) Tobias Leibner
# (2017)
#
# project.org/license.html)
# ~~~

include(Hints)

message("-- checking for cblas library")
find_library(CBLAS_LIBRARY cblas HINTS ${lib_hints})
if("${CBLAS_LIBRARY}" MATCHES "CBLAS_LIBRARY-NOTFOUND")
  message("--   CBLAS library not found, make sure you have CBLAS installed")
else("${CBLAS_LIBRARY}" MATCHES "CBLAS_LIBRARY-NOTFOUND")
  message("--   found CBLAS library")
  set(CBLAS_LIBRARIES "${CBLAS_LIBRARY}")
endif("${CBLAS_LIBRARY}" MATCHES "CBLAS_LIBRARY-NOTFOUND")

message("-- checking for cblas.h header")
find_path(CBLAS_INCLUDE_DIRS cblas.h HINTS ${include_hints})
if("${CBLAS_INCLUDE_DIRS}" MATCHES "CBLAS_INCLUDE_DIRS-NOTFOUND")
  message("--   cblas.h header not found")
else("${CBLAS_INCLUDE_DIRS}" MATCHES "CBLAS_INCLUDE_DIRS-NOTFOUND")
  message("--   found cblas.h header")
  include_sys_dir("${CBLAS_INCLUDE_DIRS}")
endif("${CBLAS_INCLUDE_DIRS}" MATCHES "CBLAS_INCLUDE_DIRS-NOTFOUND")

if(NOT DEFINED CBLAS_FOUND)
  set(CBLAS_FOUND 0)
endif(NOT DEFINED CBLAS_FOUND)
if(NOT "${CBLAS_LIBRARY}" MATCHES "CBLAS_LIBRARY-NOTFOUND")
  if(NOT "${CBLAS_INCLUDE_DIRS}" MATCHES "CBLAS_INCLUDE_DIRS-NOTFOUND")
    set(CBLAS_FOUND 1)
  endif(NOT "${CBLAS_INCLUDE_DIRS}" MATCHES "CBLAS_INCLUDE_DIRS-NOTFOUND")
endif(NOT "${CBLAS_LIBRARY}" MATCHES "CBLAS_LIBRARY-NOTFOUND")

set(HAVE_CBLAS ${CBLAS_FOUND})

# register all CBLAS related flags
if(CBLAS_FOUND)
  dune_register_package_flags(LIBRARIES "${CBLAS_LIBRARIES}" INCLUDE_DIRS "${CBLAS_INCLUDE_DIRS}")
endif()
