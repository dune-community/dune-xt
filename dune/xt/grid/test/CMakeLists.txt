# This file is part of the dune-xt-grid project:
#   https://github.com/dune-community/dune-xt-grid
# The copyright lies with the authors of this file (see below).
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

enable_testing()

# link grid files
dune_symlink_to_source_files(FILES gmsh_2d_simplices.msh dgf_1d_interval.dgf dgf_2d_interval.dgf dgf_3d_interval.dgf dgf_4d_interval.dgf)

BEGIN_TESTCASES()

add_format(${test_sources})
 
END_TESTCASES()
