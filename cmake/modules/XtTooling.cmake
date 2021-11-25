# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017, 2020)
#   René Fritze     (2016 - 2020)
#   Tobias Leibner  (2016, 2018 - 2020)
# ~~~

macro(ADD_FORMAT)
  message(NOTICE "clang-format support was moved to a `pre-commit` hook. See dune-xt/.pre-commit-config.yaml.")
  add_custom_target("format" echo 'clang-format is now a pre-commit hook' && false)
endmacro(ADD_FORMAT)

macro(ADD_TIDY)
  add_custom_target("tidy" echo 'clang-tidy support was removed' && false)
endmacro(ADD_TIDY)

macro(ADD_FORCED_DOXYGEN_TARGET)
  add_doxygen_target()
  if(TARGET doxygen_${ProjectName})
    add_custom_target(doxygen_${ProjectName}_pre_build COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/html)
    add_dependencies(doxygen_${ProjectName} doxygen_${ProjectName}_pre_build)
  endif()
endmacro(ADD_FORCED_DOXYGEN_TARGET)

macro(DEPENDENCYCHECK)
  add_custom_target(dependencycheck SOURCES ${ARGN})
  foreach(header ${ARGN})
    string(REPLACE "/" "_" fn ${header})
    set(TEST_NAME "dependencycheck_${fn}")
    to_list_spaces(CMAKE_CXX_FLAGS TEST_NAME_FLAGS)
    set(XARGS ${TEST_NAME_FLAGS} -DHAVE_CONFIG_H -H -c ${header} -w)
    add_custom_target(${TEST_NAME} + ${dune-xt_SOURCE_DIR}/cmake/dependencyinfo.py ${CMAKE_CXX_COMPILER} ${XARGS}
                                   ${CMAKE_CURRENT_SOURCE_DIR} ${fn}.dep)
    add_dependencies(dependencycheck ${TEST_NAME})
  endforeach(header)
endmacro(DEPENDENCYCHECK)
