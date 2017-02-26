__name = ESV2007
__exec_suffix = {testcase_short}_{gridname}_{dimDomain}d_r{dimRange}_rc{dimRangeCols}

include functions.mini

testcase1 = Dune::XT::Functions::ESV2007::Testcase1Force<{entity_type}, double, 2, double, 1, 1>
testcase2 = Dune::XT::Functions::ESV2007::Testcase1ExactSolution<{entity_type}, double, 2, double, 1, 1>
testcase_short = force, exact | expand testcase

dimDomain = 2

[__static]
TESTFUNCTIONTYPE = {testcase1}, {testcase2} | expand testcase
