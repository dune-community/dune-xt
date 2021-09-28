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
#   Tobias Leibner  (2018 - 2019)
# ~~~

import pytest


def test_tags():
    from dune.xt.grid import (
        Simplex,
        Cube,
        Pyramid,
        Prism,
        Dim,
    )
    Simplex()
    Cube()
    Pyramid()
    Prism()
    for dd in range(10):
        Dim(dd)
