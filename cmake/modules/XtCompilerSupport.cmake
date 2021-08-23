# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   René Fritze     (2016 - 2019)
#   Tobias Leibner  (2018, 2020)
# ~~~

macro(ADD_IF_SUPPORTED dest)
  foreach(flag ${ARGN})
    check_cxx_accepts_flag("${flag}" has_${flag})
    if(has_${flag})
      set(${dest} "${${dest}} ${flag}")
    else(has_${flag})
      message("compiler doesn't support: ${flag}")
    endif(has_${flag})
  endforeach(flag ${ARGN})
endmacro(ADD_IF_SUPPORTED)

macro(INCLUDE_SYS_DIR)
  foreach(ARG ${ARGN})
    if(IS_DIRECTORY ${ARG})
      include_directories(SYSTEM ${ARG}) # due to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=70129  we have to filter
                                         # what to sys-include includes
    else(IS_DIRECTORY ${ARG})
      message(STATUS "Include directory ${ARG} does not exist")
    endif(IS_DIRECTORY ${ARG})
  endforeach(ARG)
endmacro(INCLUDE_SYS_DIR)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #include <map>
    int main(void)
    {
      std::map<int, int> a;
      a.emplace(2, 2);
      return 0;
    };
" HAVE_MAP_EMPLACE)

check_cxx_source_compiles("
    void foo([[maybe_unused]] bool arg) {}
    int main(void){};
" HAS_WORKING_UNUSED_ATTRIBUTE)
