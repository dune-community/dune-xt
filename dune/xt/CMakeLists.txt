# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# The copyright lies with the authors of this file (see below).
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
add_subdirectory(grid/test EXCLUDE_FROM_ALL )

FILE( GLOB_RECURSE _header "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
add_format(${_header})
