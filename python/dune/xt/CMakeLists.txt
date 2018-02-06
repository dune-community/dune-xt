
dune_pybindxi_add_module(_grid EXCLUDE_FROM_ALL bindings.cc boundaryinfo.bindings.cc)
add_dune_alberta_flags(GRIDDIM 2 _grid)
