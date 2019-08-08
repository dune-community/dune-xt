# ~~~
# This file is part of the dune-xt-common project:
#   https://github.com/dune-community/dune-xt-common
# Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze (2018 - 2019)
# ~~~

import pytest
from dune.xt.common.test import load_all_submodule


def test_load_all():
    import dune.xt.common as xtc
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
    from dune.xt.common.test import runmodule
    runmodule(__file__)
