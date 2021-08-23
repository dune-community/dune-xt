# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2018)
#   René Fritze     (2010 - 2019)
#   Sven Kaulmann   (2013)
#   Tim Keil        (2018)
#   Tobias Leibner  (2016, 2018 - 2020)
# ~~~

cmake_minimum_required(VERSION 3.8)

project("dune-xt" CXX)

# local environment
set(ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}")
set(ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}")
set(ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

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
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# dune-xt cmake includes
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
add_subdirectory(pybind11)
dune_enable_all_packages(INCLUDE_DIRS
                         ${dune-xt_SOURCE_DIR}/dune
                         ${Python3_INCLUDE_DIRS}
                         MODULE_LIBRARIES
                         dunepybindxi
                         gtest_dune_xt
                         dunext)

target_link_libraries(dunepybindxi PRIVATE pybind11::module)
target_link_libraries(dunepybindxi PRIVATE ${Python3_LIBRARIES})
target_link_libraries(dunext INTERFACE ${Python3_LIBRARIES})
target_include_directories(dunepybindxi PRIVATE pybind11::module)

# ----------------------------------------------------------------------------------------------------
# General Settings
# ----------------------------------------------------------------------------------------------------
execute_process(COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
                OUTPUT_VARIABLE COMMIT
                ERROR_VARIABLE shell_error
                OUTPUT_STRIP_TRAILING_WHITESPACE)

# silence most warnings we are not responsible for
include(DuneUtils)
make_dependent_modules_sys_included()

include(DunePybindxiInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
dune_pybindxi_install_python_package(PATH python)
include_dependent_binary_python_dirs()

add_subdirectory(python)

include(GridUtils)

add_header_listing()

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
