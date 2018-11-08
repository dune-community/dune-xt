import grids
import itertools
from dune.xt.codegen import typeid_to_typedef_name

dim_range = [1]
dim_range_cols = [1]
dimDomain = [1, 2, 3]

multi_out = {grids.pretty_print(g[0], g[1]) : g[0] for g in grids.type_and_dim(cache, dimDomain)}

multi_out = {filename + '.cc': {'types': [(filename, grid, r, rC)
                             for r, rC in itertools.product(dim_range, dim_range_cols)]}
             for filename, grid in multi_out.items()}
