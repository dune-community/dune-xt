# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2014, 2016 - 2017)
#   René Fritze     (2012 - 2013, 2015 - 2016, 2018)
# ~~~

install(FILES AddALUGridCompileFlags.cmake
              DuneTBB.cmake
              DuneUtils.cmake
              DuneXtMacros.cmake
              FindClangFormat.cmake
              FindClangTidy.cmake
              FindEigen3.cmake
              FindFASP.cmake
              FindLIKWID.cmake
              FindTBB.cmake
              GridUtils.cmake
              XtCompilerSupport.cmake
              XtTooling.cmake
        DESTINATION ${DUNE_INSTALL_MODULEDIR})
