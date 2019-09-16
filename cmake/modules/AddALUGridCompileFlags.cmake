# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2015 - 2017)
#   René Fritze     (2015, 2018 - 2019)
#
# Module providing a convenience method to add ALUGrid compile flags.
#
# Provides add_dune_alugrid_compile_flags(target1 target2 ...), which
# adds ALUGrid flags to the targets for compilation only (see
# add_dune_alugrid_flags in dune-grid for compilation and linking).
# This function is essentially copied from AddALUGridFlags.cmake
#
# This is needed if you define an object like
#
#     add_library(obj_foo
#                 OBJECT
#                 EXCLUDE_FROM_ALL
#                 foo.cc)
#
# to be used in a library
#
#     add_library(lib_foo
#                 STATIC
#                 EXCLUDE_FROM_ALL
#                 $<TARGET_OBJECTS:obj_foo>
#
# that is used by a target
#
#     target_link_libraries(bar lib_foo)
#
# becaus while we ca use add_dune_alugrid_flags on lib_foo we can not use it
# on obj_foo. Instead we can use
#
#     add_dune_alugrid_compile_flags(obj_foo)
#
# to make it compile.
# ~~~

function(add_dune_alugrid_compile_flags)
  if(ALUGRID_FOUND)
    cmake_parse_arguments(ADD_ALU "SOURCE_ONLY;OBJECT" "" "" ${ARGN})
    if(ADD_ALU_SOURCE_ONLY)
      set(_prefix SOURCE)
      set(_source_only SOURCE_ONLY)
    else()
      if(ADD_ALU_OBJECT)
        set(_object OBJECT)
      endif(ADD_ALU_OBJECT)
      set(_prefix TARGET)
    endif()
    include_directories(${ALUGRID_INCLUDES})
    set_property(${_prefix} ${ADD_ALU_UNPARSED_ARGUMENTS} APPEND PROPERTY COMPILE_DEFINITIONS ENABLE_ALUGRID=1)
  endif(ALUGRID_FOUND)
endfunction(add_dune_alugrid_compile_flags)
