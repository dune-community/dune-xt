from itertools import product

from dune.xt import codegen

import matrices


def n_cimp(c, f):
  v = 'Dune::XT::LA::{}<{}>'.format(c, f)
  return codegen.typeid_to_typedef_name(v), v


conts = matrices.matrices(cache) + matrices.vectors(cache)
container = [n_cimp(c,f)
             for c,f in product(conts, matrices.fieldtype)
             if matrices.vector_filter(c,f)]

