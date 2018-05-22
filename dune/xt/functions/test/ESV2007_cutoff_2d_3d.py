import grids
import itertools
from dune.xt.codegen import typeid_to_typedef_name

dimDomain = [2, 3]

multi_out = {grids.pretty_print(g[0], g[1]) : g[0] for g in grids.type_and_dim(cache, dimDomain)}

multi_out = {filename + '.cc': {'types': [(filename, grid)]}
             for filename, grid in multi_out.items()}
