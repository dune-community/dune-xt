# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2014, 2016 - 2017)
#   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2020)
#   Tobias Leibner  (2020)
# ~~~

install(
  FILES DunePybindxiInstallPythonPackage.cmake
        DunePybindxiUtils.cmake
        DuneTBB.cmake
        DuneUtils.cmake
        DuneXtMacros.cmake
        DuneXTTesting.cmake
        DuneXtVersionHelper.cmake
        FindCBLAS.cmake
        FindClangTidy.cmake
        FindEigen3.cmake
        FindFASP.cmake
        FindLAPACKE.cmake
        FindLIKWID.cmake
        FindMKL.cmake
        FindMPI4PY.cmake
        FindSpe10Data.cmake
        GridUtils.cmake
        Hints.cmake
        XtCompilerSupport.cmake
        XtTooling.cmake
        versioneer.py
  DESTINATION ${DUNE_INSTALL_MODULEDIR})
set(_modules DuneGdtMacros.cmake)
install(FILES "${_modules}" DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/modules)
