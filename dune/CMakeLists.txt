# This file is part of the dune-xt-la project:
#   https://github.com/dune-community/dune-xt-la
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

add_subdirectory(xt)

install(DIRECTORY xt
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/
        FILES_MATCHING PATTERN "*.hh")
