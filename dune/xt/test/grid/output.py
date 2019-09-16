import dune.xt.grid.types as grid_types
from dune.xt.codegen import typeid_to_typedef_name as safe_name

# alberta needs manual flag adding in cmake, so we skip it here
all_grids = ((safe_name(g), g) for g in grid_types.all_types(cache, list(range(1, 4))) if 'Alberta' not in g)
