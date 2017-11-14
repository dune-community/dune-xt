__name = la_eigensolver
__exec_suffix = {matrix_short}

matrix = EigenDenseMatrix<double>, FieldMatrix<double\, 5\, 5> | expand types
matrix_short = eigendense_double, fieldmatrix_double | expand types
complex_matrix = EigenDenseMatrix<std::complex<double>>, FieldMatrix<std::complex<double>\, 5\, 5> | expand types
real_matrix = EigenDenseMatrix<double>, FieldMatrix<double\, 5\, 5> | expand types
field = double, double | expand types

EIGEN3_FOUND, 1 | expand types | cmake_guard

[__static]
TESTMATRIXSIZE = 5
TESTMATRIXTYPE = {matrix}
TESTFIELDTYPE = {field}
TESTCOMPLEXMATRIXTYPE = {complex_matrix}
TESTREALMATRIXTYPE = {real_matrix}
