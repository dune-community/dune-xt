# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013 - 2014, 2016 - 2017)
#   René Fritze     (2013, 2015, 2018 - 2019)
#   Tobias Leibner  (2020)
# ~~~

find_path(LIKWID_INCLUDE_DIR NAMES likwid.h PATHS ${LIKWID_ROOT} PATH_SUFFIXES include)
find_library(LIKWID_LIBRARY NAMES likwid PATHS ${LIKWID_ROOT} PATH_SUFFIXES lib)

mark_as_advanced(LIKWID_LIBRARY LIKWID_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set LIKWID_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIKWID LIKWID_LIBRARY LIKWID_INCLUDE_DIR)

set(HAVE_LIKWID 0)
if(LIKWID_FOUND)
  set(LIKWID_INCLUDE_DIRS ${LIKWID_INCLUDE_DIR})
  set(LIKWID_LIBRARIES ${LIKWID_LIBRARY})
  set(HAVE_LIKWID 1)
endif()
