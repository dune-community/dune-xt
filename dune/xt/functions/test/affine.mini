__name = affine
__exec_suffix = {gridname}_{dimDomain}d_r{dimRange}_rc{dimRangeCols}

dimRange = 1, 3 | expand
dimRangeCols = 1, 3 | expand

include grids.mini

[__static]
GRIDTYPE = {grid}
TESTFUNCTIONTYPE = Dune::XT::Functions::AffineFunction<{entity_type}, double, {dimDomain}, double, {dimRange}, {dimRangeCols}>
