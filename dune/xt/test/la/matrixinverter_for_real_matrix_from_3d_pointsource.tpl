// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#include <dune/xt/test/main.hxx> // <- has to come first (includes the config.h)!

#include <dune/xt/test/la/matrixinverter.hh>

{% for T_NAME, TESTMATRIXTYPE in config.testtypes %}
struct MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}
: public MatrixInverterTest<{{TESTMATRIXTYPE}}>
{
  using BaseType = MatrixInverterTest;
  using typename BaseType::MatrixType;

  /**
     The matrix below is a real non-symmetric invertible matrix from a moment model for the kinetic Pointsource test case.
       "[0.1756212892262638 0.18004482126181995 -0.49348712464381461 0.49938830949606494 -0.7073160963417815 1.0595994834402057e-06;"
       " 0.17562806606385517 0.18005184462676252 -0.49354113600539418 0.50059575375120657 0.70689735319270453 -3.769499436967368e-07;"
       " 0.17562307226079987 0.1800466692525447 -0.49350050991711036 -0.5000065175076156 0.00018887507812282846 -0.70710715811504954;"
       " 0.17562308446070105 0.18004668189625178 -0.49350060714612815 -0.50000869003275417 0.00019031361405394119 0.70710640425695015;"
       " -0.0072214111281474463 0.93288324029450198 -0.11009998093332186 -1.7482015044681947e-06 -2.35420746900079e-06 -4.2380607559371285e-09;"
       " 0.93625470097440933 -0.0077746247590777659 -0.11696151733678119 -1.8717676241478393e-06 -2.5225363177584535e-06 -4.5410877139483271e-09]"
   Matlab (this will be our reference below) gives the following inverse
       "[-0.069956619842954  -0.069956322880040 -0.069956501823745 -0.069956501289142 0.063349638850509 1.121064161778902;"
       " -0.066113473123754 -0.066113223084417 -0.066113362249636 -0.066113361799508 1.123470950632021 0.058271943290769;"
       " -0.555587502096003 -0.555615651279932 -0.555585807267011 -0.555585857939820 0.432422844944552 0.420211281044740;"
       "  0.499710573383257 0.500274796075355 -0.500006831431901 -0.500007846623773 0.000003909674199 0.000003817686226;"
       " -0.707554041861306 0.706659150542343 0.000405628342406 0.000407065756770 0.000010628642550 0.000010383891450;"
       "  0.000001450379141 0.000000012708409 -0.707107586716496 0.707105975654669 0.000000019133995 0.000000018693387]"
   In contrast, the invert() method of the Dune::FieldMatrix gives the following inverse
       "[-0.069965348200504485 -0.069779548837673677 -0.30508339894419151 -0.17022243546397972 0.063349643874073922 1.1210641697581696;"
       " -0.066589121801199844 -0.065891332461474883 0.29630032422042812 -0.17642735288715061 1.1234709478464535 0.058271937555510793;"
       " -0.55576414032886623 -0.55547179000514746 -0.50703985852133382 -0.63151350394304429 0.43242284571601636 0.42021128179190376;"
       " 0.49971057178651246 0.50027479737875802 -0.50000639672394975 -0.50000853465992079 3.9096742059732438e-06 3.8176862328221191e-06;"
       " -0.7075540462020502 0.70665915408649083 0.00040680890049179673 0.00040519484866008535 1.0628642569406787e-05 1.0383891468237634e-05;"
       " 1.4503713262541104e-06 1.271478900101307e-08 -0.70710758459122069 0.70710597228660532 1.9133995394523902e-08 1.8693386844274778e-08]"
   which is not an inverse as inverse * matrix gives
       "[0.94112680036257423 -0.060356101812430632 0.16543360359954729 0.16778295289847314 6.7643227058435479e-05 0.095361227759531528;"
       " 0.044229952510739405 1.0453440202592406 -0.12428610216989937 -0.12617770697826253 0.00054074451885857743 -0.33426904206670099;"
       " -0.0048146133851407757 -0.0049358839470007923 1.0135290560056918 0.013674997318695366 0.00022135327484569504 -0.088016112845982833;"
       " -4.4541786474508619e-08 -4.5663705791181705e-08 1.25162349731269e-07 1.0000001265223264 2.001938470488048e-09 -7.938998719977453e-07;"
       " -1.2138127276098885e-07 -1.2443862642877074e-07 3.410811850005737e-07 3.4479004458057095e-07 1.0000000054425457 -2.1577121855966705e-06;"
       " -2.1851414040509861e-10 -2.2401807192658647e-10 6.1402437727322846e-10 6.2070116263273302e-10 9.7978374586421261e-12 0.99999999611562296]"
   */
  MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}()
  {
    matrix_ = XT::Common::from_string<MatrixType>(
       "[0.1756212892262638 0.18004482126181995 -0.49348712464381461 0.49938830949606494 -0.7073160963417815 1.0595994834402057e-06;"
       " 0.17562806606385517 0.18005184462676252 -0.49354113600539418 0.50059575375120657 0.70689735319270453 -3.769499436967368e-07;"
       " 0.17562307226079987 0.1800466692525447 -0.49350050991711036 -0.5000065175076156 0.00018887507812282846 -0.70710715811504954;"
       " 0.17562308446070105 0.18004668189625178 -0.49350060714612815 -0.50000869003275417 0.00019031361405394119 0.70710640425695015;"
       " -0.0072214111281474463 0.93288324029450198 -0.11009998093332186 -1.7482015044681947e-06 -2.35420746900079e-06 -4.2380607559371285e-09;"
       " 0.93625470097440933 -0.0077746247590777659 -0.11696151733678119 -1.8717676241478393e-06 -2.5225363177584535e-06 -4.5410877139483271e-09]");
    expected_inverse_ = XT::Common::from_string<MatrixType>(
       "[-0.069956619842954  -0.069956322880040 -0.069956501823745 -0.069956501289142 0.063349638850509 1.121064161778902;"
       " -0.066113473123754 -0.066113223084417 -0.066113362249636 -0.066113361799508 1.123470950632021 0.058271943290769;"
       " -0.555587502096003 -0.555615651279932 -0.555585807267011 -0.555585857939820 0.432422844944552 0.420211281044740;"
       "  0.499710573383257 0.500274796075355 -0.500006831431901 -0.500007846623773 0.000003909674199 0.000003817686226;"
       " -0.707554041861306 0.706659150542343 0.000405628342406 0.000407065756770 0.000010628642550 0.000010383891450;"
       "  0.000001450379141 0.000000012708409 -0.707107586716496 0.707105975654669 0.000000019133995 0.000000018693387]");
    all_matrices_and_inverse_matrices_are_computed_ = true;
  }

  using BaseType::all_matrices_and_inverse_matrices_are_computed_;
  using BaseType::matrix_;
  using BaseType::expected_inverse_;
}; // struct MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}

TEST_F(MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}, exports_correct_types)
{
  exports_correct_types();
}

TEST_F(MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}, has_types_and_options)
{
  has_types_and_options();
}

TEST_F(MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}, throws_on_broken_matrix_construction)
{
  throws_on_broken_matrix_construction();
}

TEST_F(MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}, is_constructible)
{
  is_constructible();
}

TEST_F(MatrixInverterForMatrixFrom3dPointsource_{{T_NAME}}, gives_correct_inverse)
{
  gives_correct_inverse({ {"direct", "1e-13"} });
}

{% endfor %}