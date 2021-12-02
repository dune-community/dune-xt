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

# cmake-lint: disable=C0103
include(CheckCXXSourceCompiles)
include(DuneXtMacros)
include(CTest)
include(DunePybindxiInstallPythonPackage)

function(to_list_spaces list_name output_var)
  set(new_list_space)
  foreach(item ${${list_name}})
    set(new_list_space ${new_list_space} ${item})
  endforeach()
  set(${output_var}
      "${new_list_space}"
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

macro(ADD_HEADER_LISTING)
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
  foreach(mod ${ALL_DEPENDENCIES})
    file(
      GLOB_RECURSE
      HEADER_LIST
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/dune/*.hh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/dune/*.pbh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/dune/*.hxx"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/python/*.hh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/python/*.pbh"
      "${CMAKE_CURRENT_SOURCE_DIR}/../${mod}/python/*.hxx")
    list(APPEND COMMON_HEADER ${HEADER_LIST})
  endforeach(mod DEPENDENCIES)
  set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)
  add_custom_target(all_header SOURCES ${COMMON_HEADER})
endmacro(ADD_HEADER_LISTING)

macro(MAKE_DEPENDENT_MODULES_SYS_INCLUDED) # disable most warnings from dependent modules
  foreach(mod ${ALL_DEPENDENCIES})
    if(${mod}_INCLUDE_DIRS)
      foreach(idir ${${mod}_INCLUDE_DIRS})
        include_sys_dir(${idir})
      endforeach(idir)
    endif(${mod}_INCLUDE_DIRS)
  endforeach(mod DEPENDENCIES)
endmacro(MAKE_DEPENDENT_MODULES_SYS_INCLUDED)

macro(ADD_PYLICENSE)
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
endmacro(ADD_PYLICENSE)

function(dump_cmake_variables)
  get_cmake_property(_variableNames VARIABLES)
  list(SORT _variableNames)
  foreach(variable_name ${_variableNames})
    if(ARGV0)
      unset(MATCHED)
      string(REGEX MATCH ${ARGV0} MATCHED ${variable_name})
      if(NOT MATCHED)
        continue()
      endif()
    endif()
    message(STATUS "${variable_name}=${${variable_name}}")
  endforeach()
endfunction()
