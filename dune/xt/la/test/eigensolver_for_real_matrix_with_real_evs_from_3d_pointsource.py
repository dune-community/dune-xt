
from dune.xt.codegen import typeid_to_typedef_name as safe_name, have_eigen

matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 32, 32>']
field = ['double', 'double']
complex_matrix = ['EigenDenseMatrix<std::complex<double>>', 'FieldMatrix<std::complex<double>, 32, 32>']
real_matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 32, 32>']

def _ok(ft):
    if 'Eigen' in ft[0]:
        return have_eigen(cache)
    return True

testtypes = [(safe_name('_'.join(ft)), *ft) for ft in zip(matrix, field, complex_matrix, real_matrix) if _ok(ft)]
