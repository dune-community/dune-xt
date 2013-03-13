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
LA::Container::EigenDenseMatrix<ElementType>*
createIdentityMatrix(const LA::Container::EigenDenseMatrix<ElementType>& /*emptyDummy*/)
{
  typedef LA::Container::EigenDenseMatrix<ElementType> MatrixType;
  MatrixType* matrix = new MatrixType(dimension, dimension);
  for (size_t ii = 0; ii < dimension; ++ii)
    matrix->set(ii, ii, ElementType(1));
  return matrix;
}

template <class ElementType>
LA::Container::EigenRowMajorSparseMatrix<ElementType>*
createIdentityMatrix(const LA::Container::EigenRowMajorSparseMatrix<ElementType>& /*emptyDummy*/)
{
  typedef LA::Container::EigenRowMajorSparseMatrix<ElementType> MatrixType;
  typename LA::Container::SparsityPatternDefault pattern(dimension);
  for (size_t ii = 0; ii < dimension; ++ii)
    pattern.set(ii).insert(ii);
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
LA::Container::EigenDenseVector<ElementType>*
createVector(const LA::Container::EigenDenseVector<ElementType>& /*emptyDummy*/)
{
  typedef LA::Container::EigenDenseVector<ElementType> VectorType;
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
    for (std::string solverType : LA::Solver::types()) {
      const LA::Solver::Interface<MatrixType, VectorType>* solver =
          LA::Solver::create<MatrixType, VectorType>(solverType);
      solver->apply(*A, *f, x);
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
  typedef boost::mpl::vector<LA::Container::EigenDenseMatrix<double>, LA::Container::EigenRowMajorSparseMatrix<double>>
      MatrixTypes;
  typedef boost::mpl::vector<LA::Container::EigenDenseVector<double>> VectorTypes;
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
