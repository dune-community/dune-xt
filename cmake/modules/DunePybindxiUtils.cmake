# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#
# The code below is a renamed copy of parts of ../../pybind11/CMakeLists.txt,
# see ../../pybind11/LICENSE for license information.

# Build a Python extension module:
# dune_pybindxi_add_module(<name> source1 [source2 ...])
#
function(dune_pybindxi_add_module target_name)
  add_library(${target_name} MODULE ${ARGN})
  target_include_directories(${target_name}
    PRIVATE ${PYBIND11_INCLUDE_DIR}
    PRIVATE ${PYTHON_INCLUDE_DIRS})

  # The prefix and extension are provided by FindPythonLibsNew.cmake
  set_target_properties(${target_name} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
  set_target_properties(${target_name} PROPERTIES SUFFIX "${PYTHON_MODULE_EXTENSION}")

  if(WIN32 OR CYGWIN)
    # Link against the Python shared library on Windows
    target_link_libraries(${target_name} PRIVATE ${PYTHON_LIBRARIES})
  elseif(APPLE)
    # It's quite common to have multiple copies of the same Python version
    # installed on one's system. E.g.: one copy from the OS and another copy
    # that's statically linked into an application like Blender or Maya.
    # If we link our plugin library against the OS Python here and import it
    # into Blender or Maya later on, this will cause segfaults when multiple
    # conflicting Python instances are active at the same time (even when they
    # are of the same version).

    # Windows is not affected by this issue since it handles DLL imports
    # differently. The solution for Linux and Mac OS is simple: we just don't
    # link against the Python library. The resulting shared library will have
    # missing symbols, but that's perfectly fine -- they will be resolved at
    # import time.

    target_link_libraries(${target_name} PRIVATE "-undefined dynamic_lookup")
  endif()

  if(NOT MSVC)
    # Make sure C++11/14 are enabled
    target_compile_options(${target_name} PUBLIC ${PYBIND11_CPP_STANDARD})

    # Enable link time optimization and set the default symbol
    # visibility to hidden (very important to obtain small binaries)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      # Check for Link Time Optimization support (GCC/Clang)
      check_cxx_compiler_flag("-flto" HAS_LTO_FLAG)
      if(HAS_LTO_FLAG AND NOT CYGWIN)
        target_compile_options(${target_name} PRIVATE -flto)
      endif()

      # Intel equivalent to LTO is called IPO
      if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        check_cxx_compiler_flag("-ipo" HAS_IPO_FLAG)
        if(HAS_IPO_FLAG)
          target_compile_options(${target_name} PRIVATE -ipo)
        endif()
      endif()

      # Default symbol visibility
      target_compile_options(${target_name} PRIVATE "-fvisibility=hidden")

      # Strip unnecessary sections of the binary on Linux/Mac OS
      if(CMAKE_STRIP)
        if(APPLE)
          add_custom_command(TARGET ${target_name} POST_BUILD
                             COMMAND ${CMAKE_STRIP} -u -r $<TARGET_FILE:${target_name}>)
        else()
          add_custom_command(TARGET ${target_name} POST_BUILD
                             COMMAND ${CMAKE_STRIP} $<TARGET_FILE:${target_name}>)
        endif()
      endif()
    endif()
  elseif(MSVC)
    # /MP enables multithreaded builds (relevant when there are many files), /bigobj is
    # needed for bigger binding projects due to the limit to 64k addressable sections
    target_compile_options(${target_name} PRIVATE /MP /bigobj)

    # Enforce link time code generation on MSVC, except in debug mode
    target_compile_options(${target_name} PRIVATE $<$<NOT:$<CONFIG:Debug>>:/GL>)

    # Fancy generator expressions don't work with linker flags, for reasons unknown
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_RELEASE /LTCG)
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_MINSIZEREL /LTCG)
    set_property(TARGET ${target_name} APPEND_STRING PROPERTY LINK_FLAGS_RELWITHDEBINFO /LTCG)
  endif()
endfunction()


#set(PYBIND11_HEADERS
#  include/pybind11/attr.h
#  include/pybind11/cast.h
#  include/pybind11/chrono.h
#  include/pybind11/common.h
#  include/pybind11/complex.h
#  include/pybind11/descr.h
#  include/pybind11/eigen.h
#  include/pybind11/eval.h
#  include/pybind11/functional.h
#  include/pybind11/numpy.h
#  include/pybind11/operators.h
#  include/pybind11/pybind11.h
#  include/pybind11/pytypes.h
#  include/pybind11/stl.h
#  include/pybind11/stl_bind.h
#  include/pybind11/typeid.h
#)
#string(REPLACE "include/" "${CMAKE_CURRENT_SOURCE_DIR}/include/"
#       PYBIND11_HEADERS "${PYBIND11_HEADERS}")

