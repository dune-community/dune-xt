# ~~~
# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2016 - 2018)
#   TiKeil          (2018)
#   Tobias Leibner  (2016, 2018)
# ~~~

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)

enable_testing()

begin_testcases()

end_testcases()

# load binning setup from file
if(DEFINED ENV{TRAVIS})
  include("builder_definitions.cmake")
endif(DEFINED ENV{TRAVIS})

foreach(test ${dxt_test_binaries})
  if(${test} MATCHES 2d_simplex_alberta)
    add_dune_alberta_flags(GRIDDIM 2 ${test})
  elseif(${test} MATCHES 3d_simplex_alberta)
    add_dune_alberta_flags(GRIDDIM 3 ${test})
  endif()
endforeach()
