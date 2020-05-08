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
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# dune-xt cmake includes
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

dune_enable_all_packages(INCLUDE_DIRS
                         ${dune-xt-common_SOURCE_DIR}/dune
                         MODULE_LIBRARIES
                         dunepybindxiinterface
                         dunepybindximodule
                         dunepybindxiembed
                         dunepybindxi
                         gtest_dune_xt
                         dunext)
# The -fsized-deallocation flag is not needed anymore in newer versions of pybind11. TODO: Remove once we update
# pybind11.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(dunepybindxi INTERFACE -fsized-deallocation)
  target_compile_options(dunepybindxiinterface INTERFACE -fsized-deallocation)
  target_compile_options(dunepybindximodule INTERFACE -fsized-deallocation)
  target_compile_options(dunepybindxiembed INTERFACE -fsized-deallocation)
endif()
# Build an interface library target:
add_library(pybind11::pybind11 ALIAS dunepybindxiinterface) # to match exported target
target_include_directories(dunepybindxiinterface
                           INTERFACE $<BUILD_INTERFACE:${PYBIND11_INCLUDE_DIR}>
                                     $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>
                                     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_options(dunepybindxiinterface INTERFACE $<BUILD_INTERFACE:${PYBIND11_CPP_STANDARD}>)

add_library(pybind11::module ALIAS dunepybindximodule)
if(NOT MSVC)
  target_compile_options(dunepybindximodule PRIVATE -fvisibility=hidden)
endif()
target_link_libraries(dunepybindximodule INTERFACE pybind11::pybind11)
if(WIN32 OR CYGWIN)
  target_link_libraries(dunepybindximodule INTERFACE $<BUILD_INTERFACE:${PYTHON_LIBRARIES}>)
elseif(APPLE)
  target_link_libraries(dunepybindximodule INTERFACE "-undefined dynamic_lookup")
endif()

target_link_libraries(dunepybindxiembed INTERFACE pybind11::pybind11 $<BUILD_INTERFACE:${PYTHON_LIBRARIES}>)
add_library(pybind11::embed ALIAS dunepybindxiembed)
target_link_libraries(dunepybindxi dunepybindxiembed)

# ----------------------------------------------------------------------------------------------------
# General Settings
# ----------------------------------------------------------------------------------------------------
execute_process(COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
                OUTPUT_VARIABLE COMMIT
                ERROR_VARIABLE shell_error
                OUTPUT_STRIP_TRAILING_WHITESPACE)

include(DunePybindxiInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
make_dependent_modules_sys_included()
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
