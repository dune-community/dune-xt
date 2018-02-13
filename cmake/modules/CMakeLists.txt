# This file is part of the dune-pybindx1 project:
#   https://github.com/dune-community/dune-pybindx1
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

set(modules
    DunePybindxiMacros.cmake
    DunePybindxiUtils.cmake
    FindPythonLibsNew.cmake)

install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
