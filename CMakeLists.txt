# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2018)
#   René Fritze     (2010 - 2018)
#   Sven Kaulmann   (2013)
#   Tobias Leibner  (2016, 2018)
# ~~~

project("dune-xt" CXX)

# local environment
set(ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}")
set(ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

# cmake specific
cmake_minimum_required(VERSION 3.1)

set(CMAKE_COLOR_MAKEFILE ON)

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE ${CMAKE_PROJECT_NAME}
                 dune-common
                 dune-common_DIR
                 ${PROJECT_BINARY_DIR})
endif()

# find dune-common and set the module path
find_package(dune-common REQUIRED)
find_package(dune-pybindxi REQUIRED)
list(APPEND CMAKE_MODULE_PATH
            ${dune-common_MODULE_PATH}
            ${dune-pybindxi_MODULE_PATH}
            ${PROJECT_SOURCE_DIR}/cmake/modules)

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# dune-xt-common cmake includes
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----------------------------------------------------------------------------------------------------
# General Settings
# ----------------------------------------------------------------------------------------------------
execute_process(COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
                OUTPUT_VARIABLE COMMIT
                ERROR_VARIABLE shell_error
                OUTPUT_STRIP_TRAILING_WHITESPACE)

dune_enable_all_packages(INCLUDE_DIRS
                         ${dune-xt-common_SOURCE_DIR}/dune
                         MODULE_LIBRARIES
                         dunextcommon
                         dunextgrid
                         dunextfunctions
                         dunextla
                         gtest_dune_xt_common)

include(DunePybindxiInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
dune_pybindxi_install_python_package(PATH python)
include_dependent_binary_python_dirs()
# do not change order here
include(DuneUtils)
add_subdirectory(python)

include(GridUtils)

add_header_listing()
make_dependent_modules_sys_included()

add_subdirectory(dune)
add_subdirectory("cmake/modules")
add_subdirectory(doc/doxygen)

add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})
add_pylicense()

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

# patch config.h to remove warnings about unused HAVE_GRIDTYPE
execute_process(COMMAND sed -i -e
                        "s,^\\([ ]*\\)\#if HAVE_GRIDTYPE$,\\1\#if defined(HAVE_GRIDTYPE) \\&\\& HAVE_GRIDTYPE,g"
                        ${CMAKE_CURRENT_BINARY_DIR}/config.h)
