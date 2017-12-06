from itertools import product
from matrices import matrices, latype, vectors, fieldtypes, vector_filter
from dune.xt.codegen import typeid_to_typedef_name as safe_name

testtypes = [(safe_name('{}_{}_{}'.format(*mv,f)), latype(mv[0],f), latype(mv[1],f))
             for mv,f in product(zip(matrices(cache), vectors(cache)), fieldtypes(cache))
             if vector_filter(mv[1], f)]