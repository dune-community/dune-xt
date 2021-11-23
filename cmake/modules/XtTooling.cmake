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

macro(add_analyze)
  find_program(ANALYZER clang-check)
  if(EXISTS ${ANALYZER})
    message(STATUS "adding analyze target")
    add_custom_target(analyze SOURCES ${ARGN})
    foreach(_file ${ARGN})
      string(REPLACE "/"
                     "_"
                     fn
                     ${_file})
      add_custom_target("analyze_${fn}"
                        ${ANALYZER}
                        -fixit
                        -p=${CMAKE_CURRENT_BINARY_DIR}
                        -analyze
                        ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
      add_dependencies(analyze "analyze_${fn}")
    endforeach(_file)
  else()
    message(WARNING "not adding analyze target because clang-check is missing")
  endif(EXISTS ${ANALYZER})
endmacro(add_analyze)

macro(add_format)
  message(NOTICE clang-format support was moved to a `pre-commit` hook. See dune-xt/.pre-commit-config.yaml)
  add_custom_target("format" echo 'clang-format is now a pre-commit hook' && false)
endmacro(add_format)

find_package(ClangTidy 8)
macro(add_tidy)
  if(ClangTidy_FOUND)
    dune_symlink_to_source_files(FILES .clang-tidy)
    message(STATUS "adding tidy target")
    add_custom_target(tidy)
    add_custom_target(fix_tidy)
    add_tidy_subdir(common)
    add_tidy_subdir(grid)
    add_tidy_subdir(la)
    add_tidy_subdir(functions)
  else()
    message(WARNING "not adding tidy target because clang-tidy is missing or"
                    "wrong version: ${ClangTidy_EXECUTABLE} ${ClangTidy_VERSION}")
  endif(ClangTidy_FOUND)
endmacro()

macro(add_tidy_subdir _dxt_subdir)
  set(BASE ${PROJECT_SOURCE_DIR}/dune/xt/${_dxt_subdir})
  file(GLOB_RECURSE _files
                    "${BASE}/*.hh"
                    "${BASE}/*.h"
                    "${BASE}/*.cc"
                    "${BASE}/*.cxx"
                    "${BASE}/*.cpp"
                    "${BASE}/*.c")
  set(BASE ${PROJECT_SOURCE_DIR}/python/dune/xt/${_dxt_subdir})
  file(GLOB_RECURSE _pyfiles
                    "${BASE}/*.hh"
                    "${BASE}/*.h"
                    "${BASE}/*.cc"
                    "${BASE}/*.cxx"
                    "${BASE}/*.cpp"
                    "${BASE}/*.c")
  list(APPEND _files ${_pyfiles})
  list(REMOVE_DUPLICATES _files)
  set(TIDY_ARGS
      -config=
      -format-style=file
      -extra-arg-before='-includeconfig.h'
      -p=${CMAKE_CURRENT_BINARY_DIR}
      -header-filter=\".*/dune/xt/${_dxt_subdir}.*\")
  add_custom_target(tidy_${_dxt_subdir}
                    COMMAND ${ClangTidy_EXECUTABLE} ${TIDY_ARGS}
                            -export-fixes=${CMAKE_CURRENT_BINARY_DIR}/clang-tidy.fixes ${_files}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  add_custom_target(fix_tidy_${_dxt_subdir}
                    COMMAND ${ClangTidy_EXECUTABLE} ${TIDY_ARGS} -fix ${_files}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  add_dependencies(tidy tidy_${_dxt_subdir})
  add_dependencies(fix_tidy fix_tidy_${_dxt_subdir})
endmacro(add_tidy)

macro(add_forced_doxygen_target)
  add_doxygen_target()
  if(TARGET doxygen_${ProjectName})
    add_custom_target(doxygen_${ProjectName}_pre_build COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/html)
    add_dependencies(doxygen_${ProjectName} doxygen_${ProjectName}_pre_build)
  endif()
endmacro(add_forced_doxygen_target)

macro(DEPENDENCYCHECK)
  add_custom_target(dependencycheck SOURCES ${ARGN})
  foreach(HEADER ${ARGN})
    string(REPLACE "/"
                   "_"
                   fn
                   ${HEADER})
    set(TEST_NAME "dependencycheck_${fn}")
    to_list_spaces(CMAKE_CXX_FLAGS TEST_NAME_FLAGS)
    set(XARGS ${TEST_NAME_FLAGS} -DHAVE_CONFIG_H -H -c ${HEADER} -w)
    add_custom_target(${TEST_NAME}
                      +
                      ${dune-xt_SOURCE_DIR}/cmake/dependencyinfo.py
                      ${CMAKE_CXX_COMPILER}
                      ${XARGS}
                      ${CMAKE_CURRENT_SOURCE_DIR}
                      ${fn}.dep)
    add_dependencies(dependencycheck ${TEST_NAME})
  endforeach(HEADER)
endmacro(DEPENDENCYCHECK)

macro(ADD_CPPCHECK)
  find_program(CPPCHECK_BINARY NAMES cppcheck)
  if(EXISTS ${CPPCHECK_BINARY})
    set(CPPINLINST ${CMAKE_CURRENT_BINARY_DIR}/cppcheck.files)
    if(EXISTS ${CPPINLINST})
      file(REMOVE ${CPPINLINST})
    endif()
    foreach(SOURCEFILE ${ARGN})
      file(APPEND ${CPPINLINST} "${SOURCEFILE}\n")
    endforeach(SOURCEFILE)
    to_list_spaces(CPPCHECK_FLAGS_SPLIT ${CMAKE_CXX_FLAGS})
    add_custom_target(cppcheck
                      cppcheck
                      --enable=all
                      --xml
                      --report-progress
                      --file-list=${CPPINLINST}
                      ${CPPCHECK_FLAGS_SPLIT}
                      2>cppcheck.xml)
  else(EXISTS ${CPPCHECK_BINARY})
    message(STATUS "Not adding cppcheck target because cppcheck executable not found")
  endif(EXISTS ${CPPCHECK_BINARY})
endmacro(ADD_CPPCHECK)
