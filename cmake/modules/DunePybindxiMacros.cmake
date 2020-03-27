# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#
# The code below is a renamed copy of parts of ../../pybind11/CMakeLists.txt,
# see ../../pybind11/LICENSE for license information.

# dune-python's way of forcing a version
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add a CMake parameter for choosing a desired Python version
set(DUNE_PYBINDXI_PYTHON_VERSION "" CACHE STRING "Python version to use for dune-pybindxi")

set(PYTHON_EXECUTABLE ${DUNE_PYTHON_VIRTUALENV_INTERPRETER})
unset(PYTHON_LIBRARY)
unset(PYTHON_LIBRARY CACHE)
find_package(PythonLibsNew ${PYBIND11_PYTHON_VERSION} REQUIRED)

set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} CACHE INTERNAL "")
set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES} CACHE INTERNAL "")
set(PYTHON_MODULE_PREFIX ${PYTHON_MODULE_PREFIX} CACHE INTERNAL "")
set(PYTHON_MODULE_EXTENSION ${PYTHON_MODULE_EXTENSION} CACHE INTERNAL "")

include(CheckCXXCompilerFlag)
if(NOT MSVC AND NOT DUNE_PYBINDXI_CPP_STANDARD)
    check_cxx_compiler_flag("-std=c++14" HAS_CPP14_FLAG)
    check_cxx_compiler_flag("-std=c++11" HAS_CPP11_FLAG)

    if (HAS_CPP14_FLAG)
        set(DUNE_PYBINDXI_CPP_STANDARD -std=c++14)
    elseif (HAS_CPP11_FLAG)
        set(DUNE_PYBINDXI_CPP_STANDARD -std=c++11)
    else()
        message(FATAL_ERROR "Unsupported compiler -- pybind11 requires C++11 support!")
    endif()

    set(DUNE_PYBINDXI_CPP_STANDARD ${DUNE_PYBINDXI_CPP_STANDARD} CACHE STRING
          "C++ standard flag, e.g. -std=c++11 or -std=c++14. Defaults to C++11." FORCE)
endif()

include(DunePybindxiUtils)
include(DunePybindxiHelper)

