# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
#   Tobias Leibner  (2016)

enable_testing()

# link grid files
dune_symlink_to_source_files(FILES gmsh_2d_simplices.msh dgf_1d_interval.dgf dgf_2d_interval.dgf dgf_3d_interval.dgf dgf_4d_interval.dgf)

BEGIN_TESTCASES(dunextgrid)

END_TESTCASES()

# load binning setup from file
if(DEFINED ENV{TRAVIS})
    include("builder_definitions.cmake")
endif(DEFINED ENV{TRAVIS})

add_dune_alberta_flags(GRIDDIM 2 test_dd_glued_2d)
add_dune_alberta_flags(GRIDDIM 3 test_dd_glued_3d)

foreach (test ${dxt_test_binaries})
  if (${test} MATCHES alberta_1d)
    add_dune_alberta_flags(GRIDDIM 1 ${test})
  elseif (${test} MATCHES alberta_2d)
    add_dune_alberta_flags(GRIDDIM 2 ${test})
  elseif (${test} MATCHES alberta_3d)
    add_dune_alberta_flags(GRIDDIM 3 ${test})
  endif()
endforeach()

exclude_from_headercheck(dd_subdomains_cube.hh)
