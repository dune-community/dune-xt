# ~~~
# This file is part of the dune-xt project:
#   https://github.com/dune-community/dune-xt
# Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2019)
#   René Fritze     (2018 - 2019)
#   Tobias Leibner  (2019 - 2020)
# ~~~

import pytest
from dune.xt.test.base import load_all_submodule
from dune.xt.test.grid_types import types


def test_load_all():
    import dune.xt.grid as xtc
    load_all_submodule(xtc)


def test_empty():
    from dune.xt.common._empty import Dog, Pet, Terrier

    dog = Dog('Susi')
    pet = Pet('Bello')
    ter = Terrier()

    assert ter.getName() == 'Berti'
    assert pet.getName() == 'Bello'
    assert ter.bark() == 'woof!'


def test_logging():
    import dune.xt.common.logging as lg
    lg.create(lg.log_max)
    lg.info('log info test')
    lg.error('log error test')
    lg.debug('log debug test')


def test_timings():
    from dune.xt.common.timings import instance
    timings = instance()
    timings.reset()
    timings.start("foo.bar")
    timings.stop()
    timings.output_simple()


if __name__ == '__main__':
    from dune.xt.test.base import runmodule
    runmodule(__file__)


def test_types():
    cache = {}
    rt = types.all_types(cache=cache, dims=(2, 3))
    assert len(rt) == 0
