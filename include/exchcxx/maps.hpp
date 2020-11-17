#include <exchcxx/xc_functional.hpp>
#include <map>

namespace ExchCXX {

std::map<std::string, XCFunctional::Functional> str2func = {
    {"SVWN3", XCFunctional::Functional::SVWN3},
    {"BLYP", XCFunctional::Functional::BLYP},
    {"PBE0", XCFunctional::Functional::PBE0},
    {"PBE_X", XCFunctional::Functional::PBE_X},
    {"SlaterExchange", XCFunctional::Functional::SlaterExchange}};

} // namespace ExchCXX
