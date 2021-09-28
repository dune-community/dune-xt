# ~~~
# This file is part of the dune-xt project:
#   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
# Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2020)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2019 - 2020)
# ~~~

import pytest
from dune.xt.test.base import load_all_submodule


def test_mpi4py():
    try:
        from mpi4py import MPI
    except ImportError:
        pytest.skip('optional mpi4py is missing')
        return

    mpi_comm = MPI.COMM_WORLD
    from dune.xt.common import CollectiveCommunication
    comm_def = CollectiveCommunication()
    comm = CollectiveCommunication(mpi_comm)
    assert type(comm) is type(comm_def)
    assert comm.size > 0
    assert comm.rank < comm.size
    assert comm.sum(1) == comm.size


def test_wrapper():
    try:
        from mpi4py import MPI
    except ImportError:
        pytest.skip('optional mpi4py is missing')
        return

    mpi_comm = MPI.COMM_WORLD


if __name__ == '__main__':
    from dune.xt.common.test import runmodule
    runmodule(__file__)
