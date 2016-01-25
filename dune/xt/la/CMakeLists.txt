# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

set(lib_dune_xt_la_sources
    pattern.cc)

dune_library_add_sources(dunextla SOURCES ${lib_dune_la_common_sources})

add_subdirectory(test EXCLUDE_FROM_ALL)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
