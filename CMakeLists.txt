# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

cmake_minimum_required(VERSION 3.1)
project(dune-xt-la CXX)

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
dune_enable_all_packages(INCLUDE_DIRS ${dune-xt-la_SOURCE_DIR}/dune
                         MODULE_LIBRARIES dunextla
                         VERBOSE)


# add header of this module for header listing
file(GLOB_RECURSE xtlaheader "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh")
set(COMMON_HEADER ${xtlaheader} ${DUNE_HEADERS})

# add header of dependent modules for header listing
foreach(_mod ${ALL_DEPENDENCIES})
        file(GLOB_RECURSE HEADER_LIST "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.hh")
        list(APPEND COMMON_HEADER ${HEADER_LIST})
endforeach(_mod DEPENDENCIES)
set_source_files_properties(${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1)

#disable most warnings from dependent modules
foreach(_mod ${ALL_DEPENDENCIES})
    dune_module_to_uppercase(_upper_case "${_mod}")
    if(${_mod}_INCLUDE_DIRS)
        foreach( _idir ${${_mod}_INCLUDE_DIRS} )
            add_definitions("-isystem ${_idir}")
        endforeach( _idir )
    endif(${_mod}_INCLUDE_DIRS)
endforeach(_mod DEPENDENCIES)

add_subdirectory(dune)
add_subdirectory(doc)
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
