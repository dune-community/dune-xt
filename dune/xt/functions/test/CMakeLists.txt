# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016 - 2018)
#   Tobias Leibner  (2016)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)

enable_testing()

begin_testcases()

end_testcases()

# load binning setup from file
if(DEFINED ENV{TRAVIS})
  include("builder_definitions.cmake")
endif(DEFINED ENV{TRAVIS})

foreach (test ${dxt_test_binaries})
  if (${test} MATCHES 2d_simplex_alberta)
    add_dune_alberta_flags(GRIDDIM 2 ${test})
  elseif (${test} MATCHES 3d_simplex_alberta)
    add_dune_alberta_flags(GRIDDIM 3 ${test})
  endif()
endforeach()

# link spe10 data file if present
if(NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")
  set(SPE10MODEL2DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case2a.dat")
  if(NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL2DATA}" "${SPE10MODEL2DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL2DATA_TARGET_FILENAME})
else(NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")
  file(DOWNLOAD https://www.uni-muenster.de/AMM/mirror/spe10/model2/por_perm_case2a.zip
                ${CMAKE_CURRENT_BINARY_DIR}/por_perm_case2a.zip
       SHOW_PROGRESS)
  execute_process(COMMAND unzip -o por_perm_case2a.zip WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set(SPE10MODEL2DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case2a.dat"
      CACHE FILEPATH "spe10 model 2 data file")
endif(NOT ${SPE10MODEL2DATA} STREQUAL "SPE10MODEL2DATA-NOTFOUND")

# link spe10 data file if present
if(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set(SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
else(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  file(DOWNLOAD https://www.uni-muenster.de/AMM/mirror/spe10/model1/perm_case1.zip
                ${CMAKE_CURRENT_BINARY_DIR}/perm_case1.zip
       SHOW_PROGRESS)
  execute_process(COMMAND unzip -o perm_case1.zip WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set(SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat"
      CACHE FILEPATH "spe10 model 1 data file")
endif(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
