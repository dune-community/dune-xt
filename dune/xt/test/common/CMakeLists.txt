# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2014 - 2017)
#   René Fritze     (2012 - 2016, 2018)
#   Tobias Leibner  (2015 - 2016)
# ~~~

dependencycheck(${xtcommon})

include_directories(SYSTEM ${DUNE_XT_COMMON_TEST_DIR}/gtest)
dune_library_add_sources(gtest_dune_xt_common SOURCES gtest/gtest-all.cxx)

install(DIRECTORY ${DUNE_XT_COMMON_TEST_DIR}/gtest
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/xt/test/
        FILES_MATCHING
        PATTERN "*")
install(FILES ${DUNE_XT_COMMON_TEST_DIR}/main.hxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/xt/test/)
