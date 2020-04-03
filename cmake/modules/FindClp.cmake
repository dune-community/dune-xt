# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2018, 2020)
# ~~~

include(Hints)

message("-- checking for 'Clp' library")
find_library(CLP_LIBRARY NAMES Clp HINTS ${lib_hints})
find_library(COINUTILS_LIBRARY NAMES libCoinUtils.a HINTS ${lib_hints})
if("${CLP_LIBRARY}" MATCHES "CLP_LIBRARY-NOTFOUND")
  message("--   library 'Clp' not found")
else("${CLP_LIBRARY}" MATCHES "CLP_LIBRARY-NOTFOUND")
  message("--   found Clp library")
  set(Clp_LIBRARIES "${COINUTILS_LIBRARY}")
  list(APPEND Clp_LIBRARIES "${CLP_LIBRARY}")
endif("${CLP_LIBRARY}" MATCHES "CLP_LIBRARY-NOTFOUND")

message("-- checking for coin/ClpSimplex.hpp header")
find_path(Clp_INCLUDE_DIRS coin/ClpSimplex.hpp HINTS ${include_hints})
if("${Clp_INCLUDE_DIRS}" MATCHES "Clp_INCLUDE_DIRS-NOTFOUND")
  message("--   Clp header not found")
else("${Clp_INCLUDE_DIRS}" MATCHES "Clp_INCLUDE_DIRS-NOTFOUND")
  message("--   found Clp header")
  include_sys_dir(${Clp_INCLUDE_DIRS})
endif("${Clp_INCLUDE_DIRS}" MATCHES "Clp_INCLUDE_DIRS-NOTFOUND")

set(Clp_FOUND 0)
if(NOT "${CLP_LIBRARY}" MATCHES "CLP_LIBRARY-NOTFOUND")
  if(NOT "${Clp_INCLUDE_DIRS}" MATCHES "Clp_INCLUDE_DIRS-NOTFOUND")
    set(Clp_FOUND 1)
  endif(NOT "${Clp_INCLUDE_DIRS}" MATCHES "Clp_INCLUDE_DIRS-NOTFOUND")
endif(NOT "${CLP_LIBRARY}" MATCHES "CLP_LIBRARY-NOTFOUND")

set(HAVE_CLP ${Clp_FOUND})

# register all Clp related flags
if(Clp_FOUND)
  dune_register_package_flags(LIBRARIES "${Clp_LIBRARIES}" INCLUDE_DIRS "${Clp_INCLUDE_DIRS}")
endif()
