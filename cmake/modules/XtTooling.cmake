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

# Converts the path to a source file to a unique name that can be used as a target name. If arg points to a header file,
# the corresponding headercheck target is headercheck_${arg} after calling this function. Example:
# /home/user/dune-xt-super/dune-xt/dune/xt/common/algorithm.hh becomes _dune_xt_common_algorithm.hh
function(dxt_path_to_headercheck_name arg)
  string(REPLACE ${PROJECT_SOURCE_DIR} "" rel ${${arg}})
  string(REGEX REPLACE "/" "_" final ${rel})
  set(${arg}
      ${final}
      PARENT_SCOPE)
endfunction()

# Creates targets tidy_* and fix_tidy_* for all source files. These targets apply clang-tidy and clang-tidy -fix,
# respectively, to that file. This macro also creates unqualified tidy and fix_tidy targets which apply clang-tidy to
# all source files. In addition, a target fix_tidy_parallel is created which may be used to run several clang-tidy -fix
# targets in parallel (Parallelity is controlled by the build system, e.g. ninja -j8 fix_tidy_parallel will run 8
# targets in parallel). Note that fix_tidy_parallel may apply the fixes several times to each header (once for each
# compilation unit in which the header is included) and thus often produces broken source files. So if you expect more
# than a few fixes you probably want to avoid clang_tidy_parallel and just do something else while waiting for the
# fix_tidy command to be finished.
macro(ADD_TIDY)
  find_package(ClangTidy 12)
  if(ClangTidy_FOUND)
    dune_symlink_to_source_files(FILES .clang-tidy)
    message(STATUS "adding tidy target")
    set(BASE ${PROJECT_SOURCE_DIR}/dune/xt/)
    file(
      GLOB_RECURSE
      _files
      "${BASE}/*.hh"
      "${BASE}/*.h"
      "${BASE}/*.cc"
      "${BASE}/*.cxx"
      "${BASE}/*.cpp"
      "${BASE}/*.c")
    set(BASE ${PROJECT_SOURCE_DIR}/dune/gdt/)
    file(
      GLOB_RECURSE
      _files
      "${BASE}/*.hh"
      "${BASE}/*.h"
      "${BASE}/*.cc"
      "${BASE}/*.cxx"
      "${BASE}/*.cpp"
      "${BASE}/*.c")
    set(BASE ${PROJECT_SOURCE_DIR}/python/)
    file(
      GLOB_RECURSE
      _pyfiles
      "${BASE}/*.hh"
      "${BASE}/*.h"
      "${BASE}/*.cc"
      "${BASE}/*.cxx"
      "${BASE}/*.cpp"
      "${BASE}/*.c")
    list(APPEND _files ${_pyfiles})
    list(REMOVE_DUPLICATES _files)
    set(TIDY_ARGS --quiet --config-file=${CMAKE_SOURCE_DIR}/.clang-tidy -extra-arg-before='-includeconfig.h'
                  -p=${CMAKE_BINARY_DIR})
    set(fix_tidy_commands)
    add_custom_target(tidy)
    add_custom_target(
      fix_tidy_parallel
      COMMENT "If your fixes have been applied several times to each file, run this command sequentially (-j1)")
    foreach(file ${_files})
      if(${file} MATCHES ".*/functions/expression/mathexpr.*" OR ${file} MATCHES ".*/test/gtest/gtest.*")
        continue()
      endif()
      set(targname ${file})
      dxt_path_to_headercheck_name(targname)
      # Add targets for individual files
      add_custom_target(
        tidy_${targname}
        ${ClangTidy_EXECUTABLE} ${TIDY_ARGS} ${file}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
      add_custom_target(
        fix_tidy_${targname}
        ${ClangTidy_EXECUTABLE} ${TIDY_ARGS} -fix ${file}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
      # Add these targets as dependency to the global targets
      add_dependencies(tidy tidy_${targname})
      add_dependencies(fix_tidy_parallel fix_tidy_${targname})
      # In addition, store the commands in a list to apply them all sequentially in the fix_tidy target)
      list(
        APPEND
        fix_tidy_commands
        COMMAND
        ${ClangTidy_EXECUTABLE}
        ${TIDY_ARGS}
        -fix
        ${file})
    endforeach()
    add_custom_target(
      fix_tidy
      ${fix_tidy_commands}
      COMMENT "Running clang-tidy -fix for all files, this will take a very long time..."
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  else()
    message(WARNING "not adding tidy target because clang-tidy is missing or"
                    "wrong version: ${ClangTidy_EXECUTABLE} ${ClangTidy_VERSION}")
  endif(ClangTidy_FOUND)
endmacro()

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
