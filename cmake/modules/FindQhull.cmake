# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2017 - 2018)
# ~~~

include(Hints)

message("-- checking for 'qhullstatic_r' library")
find_library(QHULL_LIBRARY qhullstatic_r HINTS ${lib_hints})
if("${QHULL_LIBRARY}" MATCHES "QHULL_LIBRARY-NOTFOUND")
  message("--   library 'qhullstatic_r' not found")
else("${QHULL_LIBRARY}" MATCHES "QHULL_LIBRARY-NOTFOUND")
  message("--   found qhullstatic_r library")
  set(Qhull_LIBRARIES "${QHULL_LIBRARY}")
endif("${QHULL_LIBRARY}" MATCHES "QHULL_LIBRARY-NOTFOUND")

find_library(QHULLCPP_LIBRARY qhullcpp HINTS ${lib_hints})
if("${QHULLCPP_LIBRARY}" MATCHES "QHULLCPP_LIBRARY-NOTFOUND")
  message("--   library 'qhullcpp' not found")
else("${QHULLCPP_LIBRARY}" MATCHES "QHULLCPP_LIBRARY-NOTFOUND")
  message("--   found qhullcpp library")
  list(APPEND Qhull_LIBRARIES "${QHULLCPP_LIBRARY}")
endif("${QHULLCPP_LIBRARY}" MATCHES "QHULLCPP_LIBRARY-NOTFOUND")

message("-- checking for qhullcpp/Qhull.h header")
find_path(Qhull_INCLUDE_DIRS libqhullcpp/Qhull.h HINTS ${include_hints})
if("${Qhull_INCLUDE_DIRS}" MATCHES "Qhull_INCLUDE_DIRS-NOTFOUND")
  message("--   qhull header not found")
else("${Qhull_INCLUDE_DIRS}" MATCHES "Qhull_INCLUDE_DIRS-NOTFOUND")
  message("--   found qhull header")
  include_sys_dir(${Qhull_INCLUDE_DIRS})
endif("${Qhull_INCLUDE_DIRS}" MATCHES "Qhull_INCLUDE_DIRS-NOTFOUND")

set(Qhull_FOUND 0)
if(NOT "${QHULL_LIBRARY}" MATCHES "QHULL_LIBRARY-NOTFOUND")
  if(NOT "${QHULLCPP_LIBRARY}" MATCHES "QHULLCPP_LIBRARY-NOTFOUND")
    if(NOT "${Qhull_INCLUDE_DIRS}" MATCHES "Qhull_INCLUDE_DIRS-NOTFOUND")
      set(Qhull_FOUND 1)
    endif(NOT "${Qhull_INCLUDE_DIRS}" MATCHES "Qhull_INCLUDE_DIRS-NOTFOUND")
  endif(NOT "${QHULLCPP_LIBRARY}" MATCHES "QHULLCPP_LIBRARY-NOTFOUND")
endif(NOT "${QHULL_LIBRARY}" MATCHES "QHULL_LIBRARY-NOTFOUND")

set(HAVE_QHULL ${Qhull_FOUND})

# register all Qhull related flags
if(Qhull_FOUND)
  dune_register_package_flags(LIBRARIES "${Qhull_LIBRARIES}" INCLUDE_DIRS "${Qhull_INCLUDE_DIRS}")
endif()
