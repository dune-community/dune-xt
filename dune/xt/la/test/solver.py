__name = la_solver
__exec_suffix = {matrix}_{vector}_{vector2}_{fieldtype_short}

include vectors.mini
include matrices.mini

__local.vector_eigen2 = EigenDenseVector, EigenMappedDenseVector | expand
vector2 = {__local.vector_common}, {__local.vector_eigen2}, {__local.vector_istl} | expand types

('{vector}' == 'EigenMappedDenseVector' or '{vector2}' == 'EigenMappedDenseVector') and '{matrix}' == 'EigenRowMajorSparseMatrix'  | exclude
'{vector2}' == 'EigenMappedDenseVector' and '{fieldtype_short}' == 'complex'  | exclude
'{vector}' == 'CommonSparseVector' or '{vector2}' == 'CommonSparseVector' | exclude
'{matrix}' == 'CommonSparseMatrixCsr'  | exclude
'{matrix}' == 'CommonSparseMatrixCsc'  | exclude
'{matrix}' == 'CommonSparseOrDenseMatrixCsr'  | exclude
'{matrix}' == 'CommonSparseOrDenseMatrixCsc'  | exclude
'{matrix}' == 'IstlRowMajorSparseMatrix' and '{fieldtype_short}' == 'complex'  | exclude

[__static]
TESTMATRIXTYPE = Dune::XT::LA::{matrix}<{fieldtype}>
TESTRHSVECTORTYPE = Dune::XT::LA::{vector}<{fieldtype}>
TESTSOLUTIONVECTORTYPE = Dune::XT::LA::{vector2}<{fieldtype}>
