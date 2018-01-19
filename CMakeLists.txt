# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

cmake_minimum_required(VERSION 3.1)
project(dune-xt-functions CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# dune-common
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${PROJECT_SOURCE_DIR}/cmake/modules"
     "${dune-common_MODULE_PATH}")
include(DuneMacros)

# dune-xt-common
find_package(dune-xt-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${dune-xt-common_MODULE_PATH}")
include(DuneUtils)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages(MODULE_LIBRARIES dunextfunctions
                         VERBOSE)

add_header_listing()
make_dependent_modules_sys_included()

add_subdirectory(dune)
add_subdirectory(doc)

# indicate that the python subdir is a package to install
dune_install_python_package(PATH python)

if(dune-pybindxi_FOUND)
  foreach(_file
          dune/__init__.py
          dune/xt/__init__.py
          dune/xt/functions/__init__.py
         )
  execute_process(COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/python/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})
add_pylicense()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
