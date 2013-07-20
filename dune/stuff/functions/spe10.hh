#ifndef DUNE_STUFF_FUNCTION_SPE10_HH
#define DUNE_STUFF_FUNCTION_SPE10_HH

#include <sstream>
#include <iostream>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/string.hh>

#include "interfaces.hh"


namespace Dune {
namespace Stuff {


// default, to allow for specialization
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class FunctionSpe10Model1
{
public:
  FunctionSpe10Model1() = delete;
}; // class FunctionSpe10Model1


template <class DomainFieldImp, class RangeFieldImp>
class FunctionSpe10Model1<DomainFieldImp, 2, RangeFieldImp, 1, 1>
    : public FunctionInterface<DomainFieldImp, 2, RangeFieldImp, 1, 1>
{
  typedef FunctionInterface<DomainFieldImp, 2, RangeFieldImp, 1, 1> BaseType;

public:
  typedef FunctionSpe10Model1<DomainFieldImp, 2, RangeFieldImp, 1, 1> ThisType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".spe10.model1";
  }

private:
  static const size_t numXelements = 100;
  static const size_t numYelements = 1;
  static const size_t numZelements = 20;
  static constexpr double minValue = 0.001;
  static constexpr double maxValue = 998.915;

public:
  FunctionSpe10Model1(const std::string filename, const DomainType& _lowerLeft, const DomainType& _upperRight,
                      //                      const std::vector< size_t >& _numElements,
                      const std::string _name = static_id(), const int _order = 0, const double _minValue = minValue,
                      const double _maxValue = maxValue)
    : lowerLeft_(_lowerLeft)
    , upperRight_(_upperRight)
    //    , numElements_(_numElements)
    , name_(_name)
    , order_(_order)
    , minValue_(_minValue)
  {
    // sanity checks
    std::stringstream msg;
    size_t throw_up = 0;
    for (int dd = 0; dd < dimDomain; ++dd) {
      if (!(lowerLeft_[dd] < upperRight_[dd])) {
        ++throw_up;
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " lowerLeft[" << dd << "] (" << lowerLeft_[dd]
            << ") has to be smaller than upperRight[" << dd << "] (" << upperRight_[dd] << ")!";
      }
      //      if (!(numElements_[dd] > 0)) {
      //        ++throw_up;
      //        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
      //            << " numElements[" << dd << "] has to be positive (is " << numElements_[dd] << ")!";
      //      }
    } // for (int dd = 0; dd < dimDomain; ++dd)
    //    if (!(numElements_[0] <= numXelements)) {
    //      ++throw_up;
    //      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
    //          << " numElements[0] has to be smaller than " << numXelements << " (is " << numElements_[0] << ")!";
    //    }
    //    if (!(numElements_[1] <= numZelements)) {
    //      ++throw_up;
    //      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
    //          << " numElements[1] has to be smaller than " << numZelements << " (is " << numElements_[1] << ")!";
    //    }
    if (!(_maxValue > minValue_)) {
      ++throw_up;
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " maxValue_ (is " << _maxValue
          << ") has to be larger than minValue_ (is " << _minValue << ")!";
    }
    scale_ = (_maxValue - minValue_) / (maxValue - minValue);
    shift_ = minValue_ - scale_ * minValue;
    // read all the data from the file
    std::ifstream datafile(filename);
    if (datafile.is_open()) {
      static const size_t entriesPerDim = numXelements * numYelements * numZelements;
      // create storage (there should be exactly 6000 values in the file)
      data_           = new double[3 * entriesPerDim];
      double tmpValue = 0.0;
      size_t counter = 0;
      while (datafile >> tmpValue && counter < 3 * entriesPerDim) {
        data_[counter] = tmpValue;
        ++counter;
      }
      datafile.close();
      if (counter != 3 * entriesPerDim) {
        ++throw_up;
        msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " wrong number of entries in '" << filename
            << "' (are " << counter << ", should be " << 3 * entriesPerDim << ")!";
      }
    } else { // if (datafile)
      ++throw_up;
      msg << "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " could not open '" << filename << "'!";
    } // if (datafile)
    // throw up, if needed
    if (throw_up)
      DUNE_THROW(Dune::RangeError, msg.str());
    //    // keep the data we need
    //    data_ = new double[numElements_[0] * numElements_[1]];


    //    // clean up
    //    delete data;
  } // FunctionSpe10Model1()

  ~FunctionSpe10Model1()
  {
    delete data_;
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["filename"]   = "perm_case1.dat";
    description["lowerLeft"]  = "[0.0; 0.0]";
    description["upperRight"] = "[762.0; 15.24]";
    description["minValue"]   = "0.001"; // Dune::Stuff::Common::toString(minValue);
    description["maxValue"]   = "998.915"; // Dune::Stuff::Common::toString(maxValue);
    //    description["numElements"] = "[100; 20]";
    description["name"]  = static_id();
    description["order"] = "0";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const Dune::Stuff::Common::ExtendedParameterTree settings)
  {
    // get data
    const std::string filenameIn                    = settings.get<std::string>("filename");
    const std::vector<DomainFieldType> lowerLeftIn  = settings.getVector<DomainFieldType>("lowerLeft", dimDomain);
    const std::vector<DomainFieldType> upperRightIn = settings.getVector<DomainFieldType>("upperRight", dimDomain);
    //    const std::vector< size_t >          numElements  = settings.getVector< size_t >(         "numElements",
    //                                                                                                 dimDomain);
    double minValueIn = minValue;
    if (settings.hasKey("minValue"))
      minValueIn      = settings.get<double>("minValue");
    double maxValueIn = maxValue;
    if (settings.hasKey("maxValue"))
      maxValueIn             = settings.get<double>("maxValue");
    const std::string nameIn = settings.get<std::string>("name", static_id());
    const int orderIn        = settings.get<int>("order", 0);
    // convert and leave the checks to the constructor
    DomainType lowerLeft;
    DomainType upperRight;
    for (int dd = 0; dd < dimDomain; ++dd) {
      lowerLeft[dd]  = lowerLeftIn[dd];
      upperRight[dd] = upperRightIn[dd];
    }
    // create and return
    return new ThisType(filenameIn, lowerLeft, upperRight, /*numElements,*/ nameIn, orderIn, minValueIn, maxValueIn);
  } // ... create(...)

  const DomainType& lowerLeft() const
  {
    return lowerLeft_;
  }

  const DomainType& upperRight() const
  {
    return upperRight_;
  }

  //  const std::vector< size_t >& numElements() const
  //  {
  //    return numElements_;
  //  }

  virtual std::string name() const
  {
    return name_;
  }

  virtual int order() const
  {
    return order_;
  }

  virtual void evaluate(const DomainType& x, RangeType& ret) const
  {
    // decide on the interval x belongs to
    Dune::FieldVector<size_t, dimDomain> interval;
    interval[0] = std::floor(numXelements * ((x[0] - lowerLeft_[0]) / (upperRight_[0] - lowerLeft_[0])));
    interval[1] = std::floor(numZelements * ((x[1] - lowerLeft_[1]) / (upperRight_[1] - lowerLeft_[1])));
    //    for (size_t dd = 0; dd < dimDomain; ++dd) {
    //      interval[dd] = std::floor(numElements_[dd]*(           (x[dd] - lowerLeft_[dd])
    //                                                  /(upperRight_[dd] - lowerLeft_[dd])));
    //    }
    //    if (interval[0] >= numElements_[0] || interval[1] >= numElements_[1]) {
    if (interval[0] >= numXelements || interval[1] >= numZelements) {
      ret[0] = minValue_;
    } else {
      const size_t index = interval[0] + numXelements * 0 + numXelements * numYelements * interval[1];
      ret[0]             = (scale_ * data_[index]) + shift_;
    }
  } // virtual void evaluate(const DomainType& x, RangeType& ret) const

private:
  const std::string filename_;
  const DomainType lowerLeft_;
  const DomainType upperRight_;
  //  const std::vector< size_t > numElements_;
  const std::string name_;
  const int order_;
  const double minValue_;
  double* data_;
  double scale_;
  double shift_;
}; // class FunctionSpe10Model1< DomainFieldImp, 2, RangeFieldImp, 1, 1 >


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_SPE10_HH
