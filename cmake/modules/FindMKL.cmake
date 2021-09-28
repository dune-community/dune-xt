# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018 - 2019)
#   Tobias Leibner (2017 - 2021)
# ~~~

include(Hints)

set(mkl_hints $ENV{MKLROOT} ${MKLROOT} ${root_hints})
set(mkl_lib_hints "")
set(mkl_include_hints "")
list(APPEND mkl_hints "/opt/intel/mkl/" "$ENV{HOME}/intel/mkl/" "/opt/intel/oneapi/mkl/latest/")
append_to_each("${mkl_hints}" "/lib/" mkl_lib_hints)
append_to_each("${mkl_hints}" "/lib/intel64/" mkl_lib_hints_intel)
list(APPEND mkl_lib_hints ${mkl_lib_hints_intel})
append_to_each("${mkl_hints}" "/include/" mkl_include_hints)

find_library(MKL_LP64_LIBRARY mkl_intel_lp64 HINTS ${mkl_lib_hints})
find_library(MKL_SEQUENTIAL_LIBRARY mkl_sequential HINTS ${mkl_lib_hints})
find_library(MKL_CORE_LIBRARY mkl_core HINTS ${mkl_lib_hints})

set(MKL_FOUND 0)
if("${MKL_LP64_LIBRARY}" MATCHES "MKL_LP64_LIBRARY-NOTFOUND")
  message("--   library 'mkl_intel_lp64' not found, make sure you have Intel MKL installed")
else("${MKL_LP64_LIBRARY}" MATCHES "MKL_LP64_LIBRARY-NOTFOUND")
  message("--   found mkl_intel_lp64 library")
  if("${MKL_SEQUENTIAL_LIBRARY}" MATCHES "MKL_SEQUENTIAL_LIBRARY-NOTFOUND")
    message("--   library 'mkl_sequential' not found, make sure you have Intel MKL installed")
  else("${MKL_SEQUENTIAL_LIBRARY}" MATCHES "MKL_SEQUENTIAL_LIBRARY-NOTFOUND")
    message("--   found mkl_sequential library")
    if("${MKL_CORE_LIBRARY}" MATCHES "MKL_CORE_LIBRARY-NOTFOUND")
      message("--   library 'mkl_core' not found, make sure you have Intel MKL installed")
    else("${MKL_CORE_LIBRARY}" MATCHES "MKL_CORE_LIBRARY-NOTFOUND")
      message("--   found mkl_core library")
      set(MKL_LIBRARIES "${MKL_LP64_LIBRARY}")
      list(APPEND MKL_LIBRARIES ${MKL_SEQUENTIAL_LIBRARY})
      list(APPEND MKL_LIBRARIES ${MKL_CORE_LIBRARY})
    endif("${MKL_CORE_LIBRARY}" MATCHES "MKL_CORE_LIBRARY-NOTFOUND")
  endif("${MKL_SEQUENTIAL_LIBRARY}" MATCHES "MKL_SEQUENTIAL_LIBRARY-NOTFOUND")
endif("${MKL_LIBRARY}" MATCHES "MKL_LP64_LIBRARY-NOTFOUND")

message("-- checking for mkl.h header")
# message(FATAL_ERROR ${mkl_include_hints})
find_path(MKL_INCLUDE_DIRS mkl.h HINTS ${mkl_include_hints})
if("${MKL_INCLUDE_DIRS}" MATCHES "MKL_INCLUDE_DIRS-NOTFOUND")
  message("--   mkl.h header not found")
else("${MKL_INCLUDE_DIRS}" MATCHES "MKL_INCLUDE_DIRS-NOTFOUND")
  message("--   found mkl.h header")
  include_sys_dir("${MKL_INCLUDE_DIRS}")
endif("${MKL_INCLUDE_DIRS}" MATCHES "MKL_INCLUDE_DIRS-NOTFOUND")

if(NOT "${MKL_LP64_LIBRARY}" MATCHES "MKL_LP64_LIBRARY-NOTFOUND")
  if(NOT "${MKL_SEQUENTIAL_LIBRARY}" MATCHES "MKL_SEQUENTIAL_LIBRARY-NOTFOUND")
    if(NOT "${MKL_CORE_LIBRARY}" MATCHES "MKL_CORE_LIBRARY-NOTFOUND")
      if(NOT "${MKL_INCLUDE_DIRS}" MATCHES "MKL_INCLUDE_DIRS-NOTFOUND")
        set(MKL_FOUND 1)
      endif(NOT "${MKL_INCLUDE_DIRS}" MATCHES "MKL_INCLUDE_DIRS-NOTFOUND")
    endif(NOT "${MKL_CORE_LIBRARY}" MATCHES "MKL_CORE_LIBRARY-NOTFOUND")
  endif(NOT "${MKL_SEQUENTIAL_LIBRARY}" MATCHES "MKL_SEQUENTIAL_LIBRARY-NOTFOUND")
endif(NOT "${MKL_LP64_LIBRARY}" MATCHES "MKL_LP64_LIBRARY-NOTFOUND")

set(HAVE_MKL ${MKL_FOUND})

# register all MKL related flags
if(MKL_FOUND)
  dune_register_package_flags(LIBRARIES "${MKL_LIBRARIES}" INCLUDE_DIRS "${MKL_INCLUDE_DIRS}")
endif()
