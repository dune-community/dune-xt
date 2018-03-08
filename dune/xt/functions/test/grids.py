import itertools
import sys
from dune.xt.codegen import typeid_to_typedef_name

dimDomain = [1, 2, 3]

def _make_alu():
    tmpl = 'Dune::ALUGrid<{d},{d},Dune::{g},Dune::{r}>'
    alu_dims = (2, 3)
    grids = []
    for d in alu_dims:
        grids.append((tmpl.format(d=d, g='cube', r='nonconforming'), d))
    for d, r in itertools.product(alu_dims, ('nonconforming', 'conforming')):
        grids.append((tmpl.format(d=d, g='simplex', r=r), d))
    return grids

def _if_active(val, guard, cache):
    return val if cache[guard] else []

def type_and_dim(cache):
    grid_alu = _if_active(_make_alu(), 'dune-alugrid', cache)
    grid_yasp = [('Dune::YaspGrid<{d},Dune::EquidistantOffsetCoordinates<double,{d}>>'.format(d=d), d) for d in dimDomain]
    grid_ug = _if_active([('Dune::UGGrid<{}>'.format(d), d) for d in dimDomain if d != 1], 'dune-uggrid', cache)
    grid_alberta = _if_active([('Dune::AlbertaGrid<{d},{d}>'.format(d=d), d) for d in dimDomain if d != 1], 'ALBERTA_FOUND', cache)
    grid_oneD = [('Dune::OneDGrid',1)]
    return grid_alberta + grid_alu + grid_yasp + grid_ug + grid_oneD

def pretty_print(grid_name, dim):
    if 'Alberta' in grid_name:
        return '{}d_simplex_albertagrid'.format(dim)
    elif 'ALU' in grid_name and 'simplex' in grid_name and 'nonconforming' in grid_name:
        return '{}d_simplex_alunonconformgrid'.format(dim)
    elif 'ALU' in grid_name and 'simplex' in grid_name and 'conforming' in grid_name:
        return '{}d_simplex_aluconformgrid'.format(dim)
    elif 'ALU' in grid_name and 'cube' in grid_name and 'nonconforming' in grid_name:
        return '{}d_cube_alunonconformgrid'.format(dim)
    elif 'Yasp' in grid_name:
        return '{}d_cube_yaspgrid'.format(dim)
    elif 'UG' in grid_name:
        return '{}d_simplex_uggrid'.format(dim)
    elif 'OneD' in grid_name:
        return '{}d_cube_onedgrid'.format(dim)
    else:
        raise RuntimeError('unknown type: {}'.format(grid_name))
