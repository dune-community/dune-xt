# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

cmake_minimum_required(VERSION 3.1)
project(dune-pybindxi CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})
#message(FATAL_ERROR "${CMAKE_MODULE_PATH}")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

if(NOT (CMAKE_VERSION VERSION_LESS 3.0))  # CMake >= 3.0
  dune_enable_all_packages(INCLUDE_DIRS ${dune-pybindxi_SOURCE_DIR}/dune
                           MODULE_LIBRARIES dunepybindxiinterface dunepybindximodule dunepybindxiembed dunepybindxi)
  # The -fsized-deallocation flag is not needed anymore in newer versions of pybind11. 
  # TODO: Remove once we update pybind11.
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(dunepybindxiinterface INTERFACE -fsized-deallocation)
  endif()
  # Build an interface library target:
  add_library(pybind11::pybind11 ALIAS dunepybindxiinterface)  # to match exported target
  target_include_directories(dunepybindxiinterface INTERFACE $<BUILD_INTERFACE:${PYBIND11_INCLUDE_DIR}>
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
else()
  dune_enable_all_packages(INCLUDE_DIRS ${dune-pybindxi_SOURCE_DIR}/dune)
endif()

add_subdirectory(dune)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

include(DunePybindxiInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
dune_pybindxi_install_python_package(PATH python)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
