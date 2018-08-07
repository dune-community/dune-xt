from dune.xt._provider import *


def default_options_cube_grid(type_str):
    factory = 'GridProviderFactory__{}'.format(type_str)
    try:
        fac = globals()[factory]
    except KeyError:
        raise TypeError('no GridProviderFactory available for Type {}'.format(type_str))
    return fac.default_config('xt.grid.gridprovider.cube')
