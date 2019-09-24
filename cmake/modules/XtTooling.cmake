# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2016 - 2019)
#   Tobias Leibner  (2016, 2018 - 2019)
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

find_package(ClangFormat 6 EXACT)
macro(add_format glob_dir)
  if(${ARGC} GREATER 1)
    message(WARNING "'add_format' API has changed. Please provide a single "
                    "search directory instead of multiple filenames")
  endif()
  if(NOT TARGET format)
    add_custom_target(format)
  endif(NOT TARGET format)
  string(REPLACE "/"
                 "_"
                 fn
                 ${glob_dir})
  message(STATUS "adding format target")
  if(ClangFormat_FOUND)
    file(GLOB_RECURSE _files
                      "${glob_dir}/*.hh"
                      "${glob_dir}/*.cc"
                      "${glob_dir}/*.cxx"
                      "${glob_dir}/*.hxx"
                      "${glob_dir}/*.cpp"
                      "${glob_dir}/*.hpp"
                      "${glob_dir}/*.h"
                      "${glob_dir}/*.c"
                      "${glob_dir}/*.pbh")
    # Formatting gtest.h causes extremely high memory usage and takes a long time
    list(FILTER _files EXCLUDE REGEX ".*gtest.h$")
    add_custom_target("format_${fn}" ${ClangFormat_EXECUTABLE} -i -style=file -fallback-style=none ${_files})
    add_dependencies(format "format_${fn}")
  else()
    message(WARNING "not adding format target because clang-format is missing or "
                    "wrong version: ${ClangFormat_EXECUTABLE} ${ClangFormat_VERSION}")
  endif(ClangFormat_FOUND)
  file(GLOB_RECURSE _pyfiles "${glob_dir}/*.py")
  add_custom_target("pyformat_${fn}"
                    ${CMAKE_CURRENT_BINARY_DIR}/run-in-dune-env
                    yapf
                    -i
                    --style=${CMAKE_CURRENT_SOURCE_DIR}/python/.style.yapf
                    ${_pyfiles})
  add_dependencies(format "pyformat_${fn}")
  file(GLOB_RECURSE _files "${glob_dir}/*.cmake" "${glob_dir}/CMakeLists.txt")
  file(GLOB_RECURSE _exclude_files "${glob_dir}/*builder_definitions.cmake")
  list(REMOVE_ITEM _files "${glob_dir}/config.h.cmake")
  list(REMOVE_ITEM _files "${_exclude_files}")
  add_custom_target("format_${fn}_cmake"
                    ${CMAKE_BINARY_DIR}/run-in-dune-env
                    cmake-format
                    -i
                    -c
                    ${glob_dir}/.cmake_format.py
                    ${_files})
  add_dependencies(format "format_${fn}_cmake")
endmacro(add_format)

find_package(ClangTidy 3.7)
macro(add_tidy glob_dir)
  if(ClangTidy_FOUND)
    message(STATUS "adding tidy target")
    if(NOT TARGET tidy)
      add_custom_target(tidy)
    endif(NOT TARGET tidy)
    string(REPLACE "/"
                   "_"
                   fn
                   ${glob_dir})
    file(GLOB_RECURSE _files "${glob_dir}/*.cc" "${glob_dir}/*.c")
    add_custom_target("tidy_${fn}" ${ClangTidy_EXECUTABLE} -p=${CMAKE_CURRENT_BINARY_DIR} ${_files})
    add_dependencies(tidy "tidy_${fn}")
  else()
    message(WARNING "not adding tidy target because clang-tidy is missing or"
                    "wrong version: ${ClangTidy_EXECUTABLE} ${ClangTidy_VERSION}")
  endif(ClangTidy_FOUND)
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
                      ${dune-xt-common_SOURCE_DIR}/cmake/dependencyinfo.py
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
