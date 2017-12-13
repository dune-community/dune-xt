from itertools import product

from dune.xt import codegen

import matrices


conts = matrices.matrices(cache) + matrices.vectors(cache)
container = [matrices.name_type_tuple(c,f)
             for c,f in product(conts, matrices.fieldtypes(cache))
             if matrices.vector_filter(c,f)]
