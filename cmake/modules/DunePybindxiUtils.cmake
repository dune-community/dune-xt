# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#
# The code below is a renamed copy of parts of ../../pybind11/CMakeLists.txt,
# see ../../pybind11/LICENSE for license information.

# Checks whether the given CXX/linker flags can compile and link a cxx file.  cxxflags and
# linkerflags are lists of flags to use.  The result variable is a unique variable name for each set
# of flags: the compilation result will be cached base on the result variable.  If the flags work,
# sets them in cxxflags_out/linkerflags_out internal cache variables (in addition to ${result}).
function(_dune_pybind11_return_if_cxx_and_linker_flags_work result cxxflags linkerflags cxxflags_out linkerflags_out)
  set(CMAKE_REQUIRED_LIBRARIES ${linkerflags})
  check_cxx_compiler_flag("${cxxflags}" ${result})
  if (${result})
    set(${cxxflags_out} "${cxxflags}" CACHE INTERNAL "" FORCE)
    set(${linkerflags_out} "${linkerflags}" CACHE INTERNAL "" FORCE)
  endif()
endfunction()

# Internal: find the appropriate link time optimization flags for this compiler
function(_dune_pybind11_add_lto_flags target_name prefer_thin_lto)
  if (NOT DEFINED PYBIND11_LTO_CXX_FLAGS)
    set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
    set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      set(cxx_append "")
      set(linker_append "")
      if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND NOT APPLE)
        # Clang Gold plugin does not support -Os; append -O3 to MinSizeRel builds to override it
        set(linker_append ";$<$<CONFIG:MinSizeRel>:-O3>")
      elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(cxx_append ";-fno-fat-lto-objects")
      endif()

      if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND prefer_thin_lto)
        _dune_pybind11_return_if_cxx_and_linker_flags_work(HAS_FLTO_THIN
          "-flto=thin${cxx_append}" "-flto=thin${linker_append}"
          PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
      endif()

      if (NOT HAS_FLTO_THIN)
        _dune_pybind11_return_if_cxx_and_linker_flags_work(HAS_FLTO
          "-flto${cxx_append}" "-flto${linker_append}"
          PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
      endif()
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      # Intel equivalent to LTO is called IPO
      _dune_pybind11_return_if_cxx_and_linker_flags_work(HAS_INTEL_IPO
      "-ipo" "-ipo" PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    elseif(MSVC)
      # cmake only interprets libraries as linker flags when they start with a - (otherwise it
      # converts /LTCG to \LTCG as if it was a Windows path).  Luckily MSVC supports passing flags
      # with - instead of /, even if it is a bit non-standard:
      _dune_pybind11_return_if_cxx_and_linker_flags_work(HAS_MSVC_GL_LTCG
        "/GL" "-LTCG" PYBIND11_LTO_CXX_FLAGS PYBIND11_LTO_LINKER_FLAGS)
    endif()

    if (PYBIND11_LTO_CXX_FLAGS)
      message(STATUS "LTO enabled")
    else()
      message(STATUS "LTO disabled (not supported by the compiler and/or linker)")
    endif()
  endif()

  # Enable LTO flags if found, except for Debug builds
  if (PYBIND11_LTO_CXX_FLAGS)
    target_compile_options(${target_name} PRIVATE "$<$<NOT:$<CONFIG:Debug>>:${PYBIND11_LTO_CXX_FLAGS}>")
  endif()
  if (PYBIND11_LTO_LINKER_FLAGS)
    target_link_libraries(${target_name} "$<$<NOT:$<CONFIG:Debug>>:${PYBIND11_LTO_LINKER_FLAGS}>")
  endif()
endfunction()

# Build a Python extension module:
# dune_pybind11_add_module(<name> [MODULE | SHARED] [EXCLUDE_FROM_ALL]
     #                     [NO_EXTRAS] [THIN_LTO] source1 [source2 ...])
#
function(dune_pybindxi_add_module target_name)
  set(options MODULE SHARED EXCLUDE_FROM_ALL NO_EXTRAS THIN_LTO)
  set(oneValueArgs LIBNAME)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "" ${ARGN})

  if(ARG_MODULE AND ARG_SHARED)
    message(FATAL_ERROR "Can't be both MODULE and SHARED")
  elseif(ARG_SHARED)
    set(lib_type SHARED)
  else()
    set(lib_type MODULE)
  endif()

  if (NOT ARG_LIBNAME)
    set(lib_name bindings)
  else()
    set(lib_name ${ARG_LIBNAME})
  endif()

  if(ARG_EXCLUDE_FROM_ALL)
    set(exclude_from_all EXCLUDE_FROM_ALL)
  endif()

  add_library(${target_name} ${lib_type} ${exclude_from_all} ${ARG_UNPARSED_ARGUMENTS})
  dune_target_link_libraries(${target_name} "${DUNE_LIB_ADD_LIBS}")

  if(ARG_SYSTEM)
    set(inc_isystem SYSTEM)
  endif()

  target_include_directories(${target_name} ${inc_isystem}
    PRIVATE ${PYBIND11_INCLUDE_DIR}  # from project CMakeLists.txt
    PRIVATE ${pybind11_INCLUDE_DIR}  # from pybind11Config
    PRIVATE ${PYTHON_INCLUDE_DIRS})

  # Python debug libraries expose slightly different objects
  # https://docs.python.org/3.6/c-api/intro.html#debugging-builds
  # https://stackoverflow.com/questions/39161202/how-to-work-around-missing-pymodule-create2-in-amd64-win-python35-d-lib
  if(PYTHON_IS_DEBUG)
    target_compile_definitions(${target_name} PRIVATE Py_DEBUG)
  endif()

  # The prefix and extension are provided by FindPythonLibsNew.cmake
  set_target_properties(${target_name} PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}")
  set_target_properties(${target_name} PROPERTIES SUFFIX "${PYTHON_MODULE_EXTENSION}")

  # -fvisibility=hidden is required to allow multiple modules compiled against
  # different pybind versions to work properly, and for some features (e.g.
  # py::module_local).  We force it on everything inside the `pybind11`
  # namespace; also turning it on for a pybind module compilation here avoids
  # potential warnings or issues from having mixed hidden/non-hidden types.
  set_target_properties(${target_name} PROPERTIES CXX_VISIBILITY_PRESET "hidden")
  set_target_properties(${target_name} PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)


  if(WIN32 OR CYGWIN)
    # Link against the Python shared library on Windows
    target_link_libraries(${target_name} ${PYTHON_LIBRARIES})
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

    target_link_libraries(${target_name} "-undefined dynamic_lookup")

    if(ARG_SHARED)
      # Suppress CMake >= 3.0 warning for shared libraries
      set_target_properties(${target_name} PROPERTIES MACOSX_RPATH ON)
    endif()
  endif()

  # Make sure C++11/14 are enabled
  target_compile_options(${target_name} PUBLIC ${PYBIND11_CPP_STANDARD} -g)
  add_dune_all_flags(${target_name})
  target_include_directories(${target_name} PUBLIC ${MPI4PY_INCLUDE_DIR})

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

  if(ARG_NO_EXTRAS)
    return()
  endif()

  if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8)
    message(WARNING "not enabling LTO although it was requested due to buggy gcc version")
  else()
    _dune_pybind11_add_lto_flags(${target_name} ${ARG_THIN_LTO})
  endif()

  if(MSVC)
    # /MP enables multithreaded builds (relevant when there are many files), /bigobj is
    # needed for bigger binding projects due to the limit to 64k addressable sections
    target_compile_options(${target_name} PRIVATE /MP /bigobj)
  endif()


endfunction()


macro(dxt_add_make_dependent_bindings)
    add_custom_target(dependent_bindings)
    if(TARGET bindings AND NOT DXT_NO_AUTO_BINDINGS_DEPENDS)
      add_dependencies(bindings dependent_bindings)
    endif()
    foreach(_mod ${ARGN} )
      dune_module_path(MODULE ${_mod}
                     RESULT ${_mod}_binary_dir
                     BUILD_DIR)
      set(tdir ${${_mod}_binary_dir})
      if(IS_DIRECTORY ${tdir})
        add_custom_target( ${_mod}_bindings
                            COMMAND ${CMAKE_COMMAND} --build ${tdir} --target bindings_no_ext -- -j1)
        add_dependencies(dependent_bindings ${_mod}_bindings)
      endif()
    endforeach()
endmacro()
