__name = la_container_matrix
__exec_suffix = {matrix}_{vector}_{fieldtype_short}

include vectors.mini

include matrices.mini

[__static]
TESTMATRIXTYPE = Dune::XT::LA::{matrix}<{fieldtype}>
TESTVECTORTYPE = Dune::XT::LA::{vector}<{fieldtype}>



