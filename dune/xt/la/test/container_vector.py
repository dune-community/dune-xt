from itertools import product
from matrices import matrices, latype, vectors, fieldtypes, vector_filter
from dune.xt.codegen import typeid_to_typedef_name as safe_name

testtypes = [(safe_name('{}_{}'.format(mv,f)), latype(mv,f))
             for mv,f in product(vectors(cache), fieldtypes(cache))
             if vector_filter(mv, f)]