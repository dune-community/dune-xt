# ~~~
# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2016 - 2018)
#   Tobias Leibner  (2016)
# ~~~

cmake_minimum_required(VERSION 3.1)
project(dune-xt-la CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE ${CMAKE_PROJECT_NAME}
                 dune-common
                 dune-common_DIR
                 ${PROJECT_BINARY_DIR})
endif()

# dune-common
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" "${dune-common_MODULE_PATH}")
include(DuneMacros)

# dune-xt-common
find_package(dune-xt-common REQUIRED)
find_package(dune-pybindxi REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${dune-xt-common_MODULE_PATH}" "${dune-pybindxi_MODULE_PATH}")
include(DuneUtils)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages(INCLUDE_DIRS ${dune-xt-la_SOURCE_DIR}/dune MODULE_LIBRARIES dunextla)

add_header_listing()
make_dependent_modules_sys_included()

add_subdirectory(dune)
add_subdirectory(doc)

include(DuneXTInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
dune_python_install_package(PATH python)
include_dependent_binary_python_dirs()
# do not change order here
add_subdirectory(python)
dxt_add_make_dependent_bindings(dune-xt-common)

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})
add_pylicense()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
