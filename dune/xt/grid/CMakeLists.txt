# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# The copyright lies with the authors of this file (see below).
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
add_subdirectory(test EXCLUDE_FROM_ALL )

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(_grid bindings.cc)
  target_link_dune_default_libraries(_grid)
  if(DUNE_XT_WITH_PYTHON_BINDINGS)
    add_custom_target(bindings ALL DEPENDS _grid)
  else()
    add_custom_target(bindings DEPENDS _grid)
  endif()
endif()

