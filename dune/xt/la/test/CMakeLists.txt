# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)

enable_testing()

BEGIN_TESTCASES(dunextla)

END_TESTCASES()

# load binning setup from file
if(DEFINED ENV{TRAVIS})
    include("builder_definitions.cmake")
endif(DEFINED ENV{TRAVIS})
