# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
#
# ~~~

macro(include_dependent_binary_python_dirs) # disable most warnings from dependent modules
  foreach(_mod ${ALL_DEPENDENCIES} ${PROJECT_NAME})
    dune_module_path(MODULE ${_mod} RESULT ${_mod}_binary_dir BUILD_DIR)
    set(tdir ${${_mod}_binary_dir})
    if(IS_DIRECTORY ${tdir})
      dune_register_package_flags(INCLUDE_DIRS ${tdir})
    endif()
  endforeach(_mod DEPENDENCIES)
endmacro(include_dependent_binary_python_dirs)

