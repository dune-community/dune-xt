# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017, 2019)
#   René Fritze     (2018 - 2019)
#   Tim Keil        (2018)
#   Tobias Leibner  (2019 - 2020)
# ~~~

from dune.xt import guarded_import

for mod_name in (
    '_functions_function_interface_1d',
    '_functions_function_interface_2d',
    '_functions_checkerboard',
    '_functions_constant',
    '_functions_expression',
    '_functions_function_as_grid_function',
    '_functions_function_interface_3d',
    '_functions_gridfunction',
    '_functions_gridfunction_interface_1d',
    '_functions_gridfunction_interface_2d',
    '_functions_gridfunction_interface_3d',
    '_functions_indicator',
    '_functions_spe10',
    ):
    guarded_import(globals(), 'dune.xt.functions', mod_name)

