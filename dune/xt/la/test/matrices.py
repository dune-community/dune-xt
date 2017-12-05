from dune.xt import codegen

def matrices(cache):
  mat = ['CommonDenseMatrix', 'CommonSparseMatrixCsr', 'CommonSparseMatrixCsc',
    'CommonSparseOrDenseMatrixCsr', 'CommonSparseOrDenseMatrixCsc']

  if codegen.have_eigen(cache):
      mat.extend(['EigenRowMajorSparseMatrix', 'EigenDenseMatrix'])
  if codegen.have_eigen(cache):
      mat.append('IstlRowMajorSparseMatrix')
  return mat

def vectors(cache):
  vecs = ['CommonDenseVector', 'CommonSparseVector']
  if codegen.have_eigen(cache):
      vecs.extend(['EigenDenseVector', 'EigenMappedDenseVector'])
  if codegen.have_istl(cache):
      vecs.append('IstlDenseVector')
  return vecs

fieldtype = ['double', 'std::complex<double>']
#fieldtype_short = double, complex | expand field


#'{vector}' ==  and '{fieldtype_short}' == 'complex'  | exclude
def vector_filter(vector, field):
    return not (vector == 'EigenMappedDenseVector' and 'complex' in fieldtype)
