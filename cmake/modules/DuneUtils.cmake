# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2019)
#   René Fritze     (2010 - 2019)
#   Sven Kaulmann   (2013)
#   Tobias Leibner  (2015 - 2020)
# ~~~

include(CheckCXXSourceCompiles)
include(DuneXtMacros)
include(CTest)
include(DunePybindxiInstallPythonPackage)

function(TO_LIST_SPACES _LIST_NAME OUTPUT_VAR)
  set(NEW_LIST_SPACE)
  foreach(ITEM ${${_LIST_NAME}})
    set(NEW_LIST_SPACE ${NEW_LIST_SPACE} ${ITEM})
  endforeach()
  set(${OUTPUT_VAR}
      "${NEW_LIST_SPACE}"
      PARENT_SCOPE)
endfunction()

if(DEFINED dune-xt_DIR)
  set(dune-xt-path ${dune-xt_DIR})
else(DEFINED dune-xt_DIR)
  set(dune-xt-path ${dune-xt_SOURCE_DIR})
endif(DEFINED dune-xt_DIR)
if(DEFINED dune-xt_MODULE_PATH) # dependent modules
  set(dune-xt-module-path ${dune-xt_MODULE_PATH})
else(DEFINED dune-xt_MODULE_PATH) # dune-xt itself
  set(dune-xt-module-path ${PROJECT_SOURCE_DIR}/cmake/modules)
endif(DEFINED dune-xt_MODULE_PATH)

enable_testing()

include(DuneXTTesting)

macro(add_header_listing) # header
  file(
    GLOB_RECURSE
    current_module_header
    "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.pbh"
    "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/python/*.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/python/*.pbh"
    "${CMAKE_CURRENT_SOURCE_DIR}/python/*.hxx")
  set(COMMON_HEADER ${current_module_header} ${DUNE_HEADERS})

  # add header of dependent modules for header listing
  foreach(_mod ${ALL_DEPENDENCIES})
    file(
      GLOB_RECURSE
      HEADER_LIST
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/dune/*.hh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/dune/*.pbh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/dune/*.hxx"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/python/*.hh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/python/*.pbh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/python/*.hxx")
    list(APPEND COMMON_HEADER ${HEADER_LIST})
  endforeach(_mod DEPENDENCIES)
  set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)
  add_custom_target(all_header SOURCES ${COMMON_HEADER})
endmacro(add_header_listing)

macro(make_dependent_modules_sys_included) # disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES})
    if(${_mod}_INCLUDE_DIRS)
      foreach(_idir ${${_mod}_INCLUDE_DIRS})
        include_sys_dir(${_idir})
      endforeach(_idir)
    endif(${_mod}_INCLUDE_DIRS)
  endforeach(_mod DEPENDENCIES)
endmacro(make_dependent_modules_sys_included)

macro(add_pylicense)
  file(GLOB configs ${PROJECT_SOURCE_DIR}/.pylicense*.py)
  foreach(cfg ${configs})
    string(REGEX REPLACE ".*/([^/]*)" "\\1" cfg_target ${cfg})
    string(REPLACE ${PROJECT_SOURCE_DIR} "" cfg_target ${cfg_target})
    string(REGEX REPLACE "(.*)/[^/]*" "\\1" cfg_target ${cfg_target})
    string(REGEX REPLACE "/" "_" cfg_target ${cfg_target})
    list(APPEND cfg_targets ${cfg_target})
    add_custom_target(
      ${cfg_target}
      ${CMAKE_BINARY_DIR}/run-in-dune-env pylicense "--cfg=${cfg}" "${PROJECT_SOURCE_DIR}"
      VERBATIM USES_TERMINAL)
  endforeach(cfg ${configs})
  add_custom_target(license DEPENDS ${cfg_targets})
endmacro(add_pylicense)

function(dump_cmake_variables)
  get_cmake_property(_variableNames VARIABLES)
  list(SORT _variableNames)
  foreach(_variableName ${_variableNames})
    if(ARGV0)
      unset(MATCHED)
      string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
      if(NOT MATCHED)
        continue()
      endif()
    endif()
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endfunction()
