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

# add header of this module for header listing
file(GLOB_RECURSE xtgrid "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh")
set(COMMON_HEADER ${xtgrid} ${DUNE_HEADERS})

# add header of dependent modules for header listing
foreach(_mod ${ALL_DEPENDENCIES})
        file(GLOB_RECURSE HEADER_LIST "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.hh")
        list(APPEND COMMON_HEADER ${HEADER_LIST})
endforeach(_mod DEPENDENCIES)
set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)

add_subdirectory(dune)
add_subdirectory(doc)

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
