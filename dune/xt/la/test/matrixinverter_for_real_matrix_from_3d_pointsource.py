
from dune.xt.codegen import typeid_to_typedef_name as safe_name, have_eigen

matrix = ['EigenDenseMatrix<double>', 'FieldMatrix<double, 6, 6>']

def _ok(ft):
    if 'Eigen' in ft[0]:
        return have_eigen(cache)
    return True

testtypes = [(safe_name('_'.join(ft)), ft) for ft in matrix if _ok(ft)]
