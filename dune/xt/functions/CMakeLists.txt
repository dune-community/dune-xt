# This file is part of the dune-xt-functions project:
#   https://github.com/dune-community/dune-xt-functions
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

set(lib_dune_xt_functions_sources
    expression/mathexpr.cc)

FILE( GLOB_RECURSE _header "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
add_format(${lib_dune_xt_functions_sources} ${_header})

dune_library_add_sources(dunextfunctions SOURCES ${lib_dune_xt_functions_sources})

