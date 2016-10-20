# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

add_subdirectory(functions)

add_subdirectory(functions/test EXCLUDE_FROM_ALL)

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(functions functions/bindings.cc)
  target_link_dune_default_libraries(functions)
endif()

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
