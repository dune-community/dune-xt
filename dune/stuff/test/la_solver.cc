#include "test_common.hh"

#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/stuff/common/tuple.hh>

using namespace Dune;
using namespace Dune::Stuff;

static const size_t dimension = 5;


template <class MatrixType>
MatrixType* createIdentityMatrix(const MatrixType& /*emptyDummy*/);

#if HAVE_EIGEN
template <class ElementType>
LA::EigenDenseMatrix<ElementType>* createIdentityMatrix(const LA::EigenDenseMatrix<ElementType>& /*emptyDummy*/)
{
  typedef LA::EigenDenseMatrix<ElementType> MatrixType;
  MatrixType* matrix = new MatrixType(dimension, dimension);
  for (size_t ii = 0; ii < dimension; ++ii)
    matrix->set(ii, ii, ElementType(1));
  return matrix;
}

template <class ElementType>
LA::EigenRowMajorSparseMatrix<ElementType>*
createIdentityMatrix(const LA::EigenRowMajorSparseMatrix<ElementType>& /*emptyDummy*/)
{
  typedef LA::EigenRowMajorSparseMatrix<ElementType> MatrixType;
  typename LA::SparsityPatternDefault pattern(dimension);
  for (size_t ii = 0; ii < dimension; ++ii)
    pattern.inner(ii).insert(ii);
  MatrixType* matrix = new MatrixType(dimension, dimension, pattern);
  for (size_t ii = 0; ii < dimension; ++ii)
    matrix->set(ii, ii, ElementType(1));
  return matrix;
}
#endif // HAVE_EIGEN

template <class VectorType>
VectorType* createVector(const VectorType& /*emptyDummy*/);

#if HAVE_EIGEN
template <class ElementType>
LA::EigenDenseVector<ElementType>* createVector(const LA::EigenDenseVector<ElementType>& /*emptyDummy*/)
{
  typedef LA::EigenDenseVector<ElementType> VectorType;
  VectorType* vector = new VectorType(dimension);
  for (size_t ii = 0; ii < dimension; ++ii)
    vector->set(ii, ElementType(1));
  return vector;
}
#endif // HAVE_EIGEN


struct SolverBaseTest
{
  template <class MatrixType, class VectorType>
  static void run()
  {
    const MatrixType* A = createIdentityMatrix(MatrixType());
    const VectorType* f = createVector(VectorType());
    VectorType x;
    for (std::string solverType : LA::solverTypes()) {
      const LA::SolverInterface<MatrixType, VectorType>* solver = LA::createSolver<MatrixType, VectorType>(solverType);
      const auto solverSettings = LA::solverDefaultSettings<MatrixType, VectorType>(solverType);
      solver->apply(*A, *f, x, solverSettings);
      delete solver;
    }
    delete A;
    delete f;
  }
};


#if HAVE_EIGEN
template <class TestFunctor>
struct SolverTest : public ::testing::Test
{
  typedef boost::mpl::vector</*LA::EigenDenseMatrix< double >,*/
                             LA::EigenRowMajorSparseMatrix<double>> MatrixTypes;
  typedef boost::mpl::vector<LA::EigenDenseVector<double>> VectorTypes;
  typedef typename DSC::TupleProduct::Combine<MatrixTypes, VectorTypes, TestFunctor>::template Generate<>
      base_generator_type;
  void run()
  {
    base_generator_type::Run();
  }
}; // struct FunctionTest
#endif // HAVE_EIGEN


typedef ::testing::Types<SolverBaseTest> SolverTestTypes;

#if HAVE_EIGEN
TYPED_TEST_CASE(SolverTest, SolverTestTypes);
TYPED_TEST(SolverTest, All)
{
  this->run();
}
#endif // HAVE_EIGEN


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
