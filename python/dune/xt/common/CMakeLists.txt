# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2018 - 2019)
# ~~~

dune_pybindxi_add_module(logging EXCLUDE_FROM_ALL logging.cc)
dune_pybindxi_add_module(timedlogging EXCLUDE_FROM_ALL timedlogging.cc)

dune_pybindxi_add_module(_common_empty EXCLUDE_FROM_ALL empty.cc)
dune_pybindxi_add_module(_common_exceptions EXCLUDE_FROM_ALL exceptions.cc)
dune_pybindxi_add_module(_common_mpi EXCLUDE_FROM_ALL mpi.cc)
dune_pybindxi_add_module(_common_timings EXCLUDE_FROM_ALL timings.cc)
