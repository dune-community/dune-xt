# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   René Fritze     (2018 - 2020)
#   Tobias Leibner  (2020 - 2021)
#
# The dune_pybind11_add_module function is a renamed copy of pybind11_add_module from ../../pybind11/tools/pybind11Tools.cmake, see ../../pybind11/LICENSE for license
# information.
# ~~~

# Build a Python extension module: dune_pybindxi_add_module(<name> [MODULE | SHARED] [EXCLUDE_FROM_ALL] [NO_EXTRAS]
# [THIN_LTO] source1 [source2 ...])
# Renamed copy of pybind11_add_module, added code blocks are marked with dune-pybindxi START/END
function(dune_pybindxi_add_module target_name)
  set(options "MODULE;SHARED;EXCLUDE_FROM_ALL;NO_EXTRAS;SYSTEM;THIN_LTO;OPT_SIZE")
  # the next two lines were added/modified compared to the original function
  set(oneValueArgs LIBNAME)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "" ${ARGN})

  if(ARG_MODULE AND ARG_SHARED)
    message(FATAL_ERROR "Can't be both MODULE and SHARED")
  elseif(ARG_SHARED)
    set(lib_type SHARED)
  else()
    set(lib_type MODULE)
  endif()

  # dune-pybindxi START
  if(NOT ARG_LIBNAME)
    set(lib_name bindings)
  else()
    set(lib_name ${ARG_LIBNAME})
  endif()
  # dune-pybindxi END

  if(ARG_EXCLUDE_FROM_ALL)
    set(exclude_from_all EXCLUDE_FROM_ALL)
  else()
    set(exclude_from_all "")
  endif()

  add_library(${target_name} ${lib_type} ${exclude_from_all} ${ARG_UNPARSED_ARGUMENTS})
  # dune-pybindxi START
  dune_target_link_libraries(${target_name} "${DUNE_LIB_ADD_LIBS}")
  add_dune_all_flags(${target_name})

  target_link_libraries(${target_name} dunepybindxi)
  target_include_directories(${target_name} PRIVATE ${PYBIND11_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
  # dune-pybindxi END

  if(ARG_SYSTEM)
    message(STATUS "Warning: this does not have an effect - use NO_SYSTEM_FROM_IMPORTED if using imported targets")
  endif()

  pybind11_extension(${target_name})

  # -fvisibility=hidden is required to allow multiple modules compiled against different pybind versions to work
  # properly, and for some features (e.g. py::module_local).  We force it on everything inside the `pybind11`
  # namespace; also turning it on for a pybind module compilation here avoids potential warnings or issues from having
  # mixed hidden/non-hidden types.
  if(NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
    set_target_properties(${target_name} PROPERTIES CXX_VISIBILITY_PRESET "hidden")
  endif()

  if(NOT DEFINED CMAKE_CUDA_VISIBILITY_PRESET)
    set_target_properties(${target_name} PROPERTIES CUDA_VISIBILITY_PRESET "hidden")
  endif()

  # dune-pybindxi START
  if(TARGET ${lib_name})
    add_dependencies(${lib_name} ${target_name})
    add_dependencies(${lib_name}_no_ext ${target_name})
  else()
    if(DUNE_XT_WITH_PYTHON_BINDINGS)
      add_custom_target(${lib_name} ALL DEPENDS ${target_name})
      add_custom_target(${lib_name}_no_ext ALL DEPENDS ${target_name})
    else()
      add_custom_target(${lib_name} DEPENDS ${target_name})
      add_custom_target(${lib_name}_no_ext DEPENDS ${target_name})
    endif()
  endif()
  # dune-pybindxi END

  if(ARG_NO_EXTRAS)
    return()
  endif()

  if(NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
    if(ARG_THIN_LTO)
      target_link_libraries(${target_name} pybind11::thin_lto)
    else()
      target_link_libraries(${target_name} pybind11::lto)
    endif()
  endif()

  if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    pybind11_strip(${target_name})
  endif()

  if(MSVC)
    target_link_libraries(${target_name} pybind11::windows_extras)
  endif()

  if(ARG_OPT_SIZE)
    target_link_libraries(${target_name} pybind11::opt_size)
  endif()
endfunction()

macro(dxt_add_make_dependent_bindings)
  add_custom_target(dependent_bindings)
  if(TARGET bindings AND NOT DXT_NO_AUTO_BINDINGS_DEPENDS)
    add_dependencies(bindings dependent_bindings)
  endif()
  foreach(_mod ${ARGN})
    dune_module_path(MODULE ${_mod} RESULT ${_mod}_binary_dir BUILD_DIR)
    set(tdir ${${_mod}_binary_dir})
    if(IS_DIRECTORY ${tdir})
      add_custom_target(${_mod}_bindings
                        COMMAND ${CMAKE_COMMAND}
                                --build ${tdir}
                                --target bindings_no_ext
                                -- -j1)
      add_dependencies(dependent_bindings ${_mod}_bindings)
    endif()
  endforeach()
endmacro()
