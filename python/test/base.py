
import pytest
from dune.xt.common.test import load_all_submodule


def test_load_all():
    import dune.xt.grid as xtc
    load_all_submodule(xtc)


