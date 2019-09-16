// This file is part of the dune-xt-common project:
//   https://github.com/dune-community/dune-xt-common
// Copyright 2009-2018 dune-xt-common developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   René Fritze     (2012 - 2018)
//   Tobias Leibner  (2014, 2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <array>
#include <ostream>

#include <boost/assign/list_of.hpp>
#include <boost/array.hpp>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/random.hh>
#include <dune/xt/common/tuple.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/validation.hh>

#include <dune/xt/common/test/float_cmp.hh>

// uncomment this for output
// std::ostream& test_out = std::cout;
std::ostream& test_out = DXTC_LOG.devnull();

using namespace Dune::XT::Common;
using Dune::XT::Common::Exceptions::results_are_not_as_expected;
using namespace Dune::XT::Common::FloatCmp;

struct CreateByOperator
{
  static Configuration create()
  {
    Configuration config;
    config["string"] = "string";
    config["sub1.int"] = "1";
    config["sub2.size_t"] = "1";
    config["sub2.subsub1.vector"] = "[0 1]";
    config["sub2.subsub1.matrix"] = "[0 1; 1 2]";
    return config;
  }
};

struct CreateByOperatorAndAssign
{
  static Configuration create()
  {
    Configuration config;
    config["string"] = "string";
    config["sub1.int"] = "1";
    config["sub2.size_t"] = "1";
    config["sub2.subsub1.vector"] = "[0 1]";
    config["sub2.subsub1.matrix"] = "[0 1; 1 2]";
    Configuration config2;
    config2 = config;
    return config2;
  }
};

struct CreateByKeyAndValueAndAddConfiguration
{
  static Configuration create()
  {
    Configuration config({"string"}, {"string"});
    config.set("sub1.int", "1");
    config.set("sub2.size_t", 1);
    config.add(Configuration({"vector"}, {"[0 1]"}), "sub2.subsub1");
    config.add(Configuration({"matrix"}, {"[0 1; 1 2]"}), "sub2.subsub1");
    Configuration config_1({"bool"}, {true});
    Configuration config_2({"int"}, {int(1)});
    return config;
  }
};

struct CreateByKeyAndValueAndAddParameterTree
{
  static Configuration create()
  {
    const auto str = std::string("string");
    Configuration config({str}, {str});
    config.set("sub1.int", "1");
    config.set("sub2.size_t", 1);
    Dune::ParameterTree paramtree;
    paramtree["vector"] = "[0 1]";
    paramtree["matrix"] = "[0 1; 1 2]";
    config.add(paramtree, "sub2.subsub1");
    return config;
  }
};

struct CreateByKeyAndValueVectorsAndAddParameterTree
{
  static Configuration create()
  {
    using namespace std;
    const auto str = string("string");
    Configuration config(vector<string>({str, "sub1.int"}), {str, string("1")});
    config.set("sub2.size_t", 1);
    Dune::ParameterTree paramtree;
    paramtree["vector"] = "[0 1]";
    paramtree["matrix"] = "[0 1; 1 2]";
    config.add(paramtree, "sub2.subsub1");
    return config;
  }
};

struct CreateByKeysAndValues
{
  static Configuration create()
  {
    return Configuration({"string", "sub1.int", "sub2.size_t", "sub2.subsub1.vector", "sub2.subsub1.matrix"},
                         {"string", "1", "1", "[0 1]", "[0 1; 1 2]"});
  }
};

struct CreateByParameterTree
{
  static Configuration create()
  {
    Dune::ParameterTree paramtree;
    paramtree["string"] = "string";
    paramtree["sub1.int"] = "1";
    paramtree["sub2.size_t"] = "1";
    paramtree["sub2.subsub1.vector"] = "[0 1]";
    paramtree["sub2.subsub1.matrix"] = "[0 1; 1 2]";
    return Configuration(paramtree);
  }
};

typedef testing::
    Types<double, float, std::string, std::complex<double>, int, unsigned int, unsigned long, long long, char>
        TestTypes;

typedef testing::Types<CreateByOperator,
                       CreateByKeyAndValueAndAddConfiguration,
                       CreateByKeyAndValueAndAddParameterTree,
                       CreateByKeyAndValueVectorsAndAddParameterTree,
                       CreateByKeysAndValues,
                       CreateByParameterTree,
                       CreateByOperatorAndAssign>
    ConfigurationCreators;

constexpr const auto SEED = std::random_device::result_type(0);

template <class T>
static DefaultRNG<T> rng_setup()
{
  return DefaultRNG<T>(std::numeric_limits<T>::min(), std::numeric_limits<T>::max(), T(SEED));
}

template <>
DefaultRNG<std::string> rng_setup()
{
  return DefaultRNG<std::string>(12, SEED);
}

template <>
DefaultRNG<std::complex<double>> rng_setup()
{
  return DefaultRNG<std::complex<double>>(-2, 2, SEED);
}

template <>
DefaultRNG<double> rng_setup()
{
  return DefaultRNG<double>(-2, 2, SEED);
}

template <class T>
static void val_compare_eq(const T& aa, const T& bb)
{
  DXTC_EXPECT_FLOAT_EQ(aa, bb);
}

static void val_compare_eq(const std::string& aa, const std::string& bb)
{
  EXPECT_EQ(aa, bb);
}

template <class T>
struct ConfigTest : public testing::Test
{
  static const int count = 2;
  DefaultRNG<T> rng;
  RandomStrings key_gen;
  // std::array is not assignable from list_of it seems. Don't make the following two arrays const (triggers boost/intel
  // compiler bug, last tested with icpc version 14.0.3).
  boost::array<T, count> values;
  boost::array<std::string, count> keys;
  ConfigTest()
    : rng(rng_setup<T>())
    , key_gen(8)
    , values(boost::assign::list_of<T>().repeat_fun(values.size() - 1, rng))
    , keys(boost::assign::list_of<std::string>().repeat_fun(values.size() - 1, key_gen))
  {}

  virtual ~ConfigTest() {}

  void get()
  {
    std::set<std::string> uniq_keys;
    for (T val : values) {
      const auto key = key_gen();
      const auto got_val = DXTC_CONFIG_GET(key, val);
      // since the value invariably goes through string conversion, we need to adjust the expected value as well
      const T adjusted_val = from_string<T>(to_string(val));
      val_compare_eq(adjusted_val, got_val);
      //! TODO add a float compare check that makes sure introduced error is only due to limited precision in str conv
      uniq_keys.insert(key);
    }
    EXPECT_EQ(values.size(), uniq_keys.size());
  }

  void set()
  {
    for (T val : values) {
      auto key = key_gen();
      // since the value invariably goes through string conversion, we need to adjust the expected value as well
      const T adjusted_val = from_string<T>(to_string(val));
      DXTC_CONFIG.set(key, val);
      // get with default diff from expected
      auto re = DXTC_CONFIG.get(key, T(val + Epsilon<T>::value));
      val_compare_eq(re, adjusted_val);
    }
  }

  void other()
  {
    auto key = this->key_gen();
    DXTC_CONFIG.set(key, T());
    EXPECT_THROW(DXTC_CONFIG.get(key, T(), ValidateNone<T>()), Exceptions::configuration_error);
  }

  void issue_42()
  {
    using namespace Dune::XT::Common;
    using namespace std;
    Configuration empty;
    Configuration to_add(vector<string>{"subsection.key"}, {0l});
    empty.add(to_add, "", true);
    EXPECT_TRUE(empty.has_sub("subsection"));
  }

  void add_overwrite(bool do_overwrite)
  {
    using namespace Dune::XT::Common;
    using namespace std;
    Configuration base(vector<string>{"subsection.key"}, {1l});
    Configuration to_add(vector<string>{"subsection.otherkey"}, {0l});
    base.add(to_add, "", do_overwrite);
    EXPECT_TRUE(base.has_sub("subsection"));
    EXPECT_EQ(base.get<long>("subsection.key"), 1l);
    EXPECT_EQ(base.get<long>("subsection.otherkey"), 0l);
  }

  void subtree_ctor()
  {
    using namespace Dune::XT::Common;
    using namespace std;
    Configuration to_add(vector<string>{"key"}, {1l});
    Configuration base(to_add, string("subsection"));
    EXPECT_TRUE(base.has_sub("subsection"));
    EXPECT_EQ(base.get<long>("subsection.key"), 1l);
  }
}; // struct ConfigTest

struct StaticCheck
{
  typedef boost::mpl::vector<Int<1>, Int<2>> Ints;

  template <class MatrixType>
  static void check_matrix_static_size(const Configuration& config)
  {
    typedef MatrixAbstraction<MatrixType> MT;
    const auto r = MT::rows(MatrixType());
    const auto c = MT::cols(MatrixType());

    const auto check = [&r, &c](const MatrixType& mat) {
      for (size_t cc = 0; cc < c; ++cc) {
        for (size_t rr = 0; rr < r; ++rr) {
          val_compare_eq(MT::get_entry(mat, rr, cc), double(rr + cc));
        }
      }
    };

    check(config.get("matrix", MatrixType(), r, c));
    check(config.get("matrix", MatrixType()));
    check(config.get<MatrixType>("matrix", r, c));
    check(config.get<MatrixType>("matrix"));
  }

  template <class U, class V>
  static void run(const Configuration& config)
  {
    const auto rows = U::value;
    const auto cols = V::value;
    check_matrix_static_size<FieldMatrix<double, rows, cols>>(config);
    check_matrix_static_size<Dune::FieldMatrix<double, rows, cols>>(config);
  }
};

template <class ConfigurationCreator>
struct ConfigurationTest : public ::testing::Test
{
  template <class VectorType>
  static void check_vector(const Configuration& config)
  {
    VectorType vec = config.get("vector", VectorType(), 1);
    EXPECT_EQ(1, vec.size());
    EXPECT_FLOAT_EQ(0.0, vec[0]);
    vec = config.get("vector", VectorType(), 2);
    EXPECT_EQ(2, vec.size());
    for (auto ii : {0.0, 1.0})
      EXPECT_FLOAT_EQ(ii, vec[ii]);
    vec = config.get<VectorType>("vector", 1);
    EXPECT_EQ(1, vec.size());
    EXPECT_FLOAT_EQ(0.0, vec[0]);
    vec = config.get<VectorType>("vector", 2);
    EXPECT_EQ(2, vec.size());
    for (auto ii : {0.0, 1.0})
      EXPECT_FLOAT_EQ(ii, vec[ii]);
    vec = config.get<VectorType>("vector");
    EXPECT_EQ(2, vec.size());
    for (auto ii : {0.0, 1.0})
      EXPECT_FLOAT_EQ(ii, vec[ii]);
  } // ... check_vector< ... >(...)

  template <class K, int d>
  static void check_field_vector(const Configuration& config)
  {
    typedef FieldVector<K, d> VectorType;
    VectorType vec = config.get("vector", VectorType(), d);
    EXPECT_EQ(d, vec.size());
    for (size_t ii = 0; ii < d; ++ii)
      EXPECT_TRUE(FloatCmp::eq(vec[ii], double(ii)));

    vec = config.get<VectorType>("vector", d);
    EXPECT_EQ(vec.size(), d);

    for (size_t ii = 0; ii < d; ++ii)
      EXPECT_TRUE(FloatCmp::eq(vec[ii], double(ii)));

    vec = config.get<VectorType>("vector");
    EXPECT_EQ(vec.size(), d);

    for (size_t ii = 0; ii < d; ++ii)
      EXPECT_TRUE(FloatCmp::eq(vec[ii], double(ii)));
  } // ... check_field_vector< ... >(...)

  template <class MatrixType>
  static void check_matrix(const Configuration& config)
  {
    MatrixType mat = config.get("matrix", MatrixType(), 1, 1);
    typedef MatrixAbstraction<MatrixType> MT;
    EXPECT_FALSE(MT::rows(mat) != 1 || MT::cols(mat) != 1);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    mat = config.get("matrix", MatrixType(), 1, 2);
    EXPECT_FALSE(MT::rows(mat) != 1 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    mat = config.get("matrix", MatrixType(), 2, 1);
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 1);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    mat = config.get("matrix", MatrixType(), 2, 2);
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 1), 2.0));
    mat = config.get("matrix", MatrixType());
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 1), 2.0));

    mat = config.get<MatrixType>("matrix", 1, 1);
    EXPECT_FALSE(MT::rows(mat) != 1 || MT::cols(mat) != 1);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    mat = config.get<MatrixType>("matrix", 1, 2);
    EXPECT_FALSE(MT::rows(mat) != 1 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    mat = config.get<MatrixType>("matrix", 2, 1);
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 1);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    mat = config.get<MatrixType>("matrix", 2, 2);
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 1), 2.0));
    mat = config.get<MatrixType>("matrix");
    EXPECT_FALSE(MT::rows(mat) != 2 || MT::cols(mat) != 2);
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 0), 0.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 0, 1), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 0), 1.0));
    EXPECT_TRUE(FloatCmp::eq(MT::get_entry(mat, 1, 1), 2.0));
  } // ... check_matrix< ... >(...)

  static void behaves_correctly()
  {
    const Configuration config = ConfigurationCreator::create();
    config.report(); // <- this works as well but will produce output
    config.report(test_out);
    config.report(test_out, "'prefix '");
    test_out << config << std::endl;
    std::string report_str DUNE_UNUSED = config.report_string();
    std::string str = config.get("string", std::string("foo"));
    EXPECT_EQ(str, "string");
    str = config.get("foo", std::string("string"));
    EXPECT_EQ(str, "string");
    str = config.get<std::string>("string");
    EXPECT_EQ(str, "string");

    EXPECT_TRUE(config.has_sub("sub1"));
    Configuration sub1_config = config.sub("sub1");
    int nt = sub1_config.get("int", int(0));
    EXPECT_EQ(nt, 1);
    nt = sub1_config.get("intt", int(1));
    EXPECT_EQ(nt, 1);
    nt = sub1_config.get<int>("int");
    EXPECT_EQ(nt, 1);
    size_t st = config.get("sub2.size_t", size_t(0));
    EXPECT_EQ(st, 1);
    st = config.get("sub2.size_tt", size_t(1));
    EXPECT_EQ(st, 1);
    st = config.get<size_t>("sub2.size_t");
    EXPECT_EQ(st, 1);

    const auto subsub1 = config.sub("sub2.subsub1");
    check_vector<std::vector<double>>(subsub1);
    check_field_vector<double, 1>(subsub1);
    check_field_vector<double, 2>(subsub1);
    check_vector<Dune::DynamicVector<double>>(subsub1);
    check_matrix<Dune::DynamicMatrix<double>>(subsub1);

    TupleProduct::Combine<StaticCheck::Ints, StaticCheck::Ints, StaticCheck>::Generate<>::Run(subsub1);

  } // ... behaves_correctly(...)
}; // struct ConfigurationTest

TYPED_TEST_CASE(ConfigTest, TestTypes);
TYPED_TEST(ConfigTest, Get)
{
  this->get();
}
TYPED_TEST(ConfigTest, Set)
{
  this->set();
}
TYPED_TEST(ConfigTest, Other)
{
  this->other();
  this->issue_42();
  this->add_overwrite(true);
  this->add_overwrite(false);
  this->subtree_ctor();
}

TYPED_TEST_CASE(ConfigurationTest, ConfigurationCreators);
TYPED_TEST(ConfigurationTest, behaves_correctly)
{
  this->behaves_correctly();
}
