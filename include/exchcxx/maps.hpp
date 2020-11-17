#include <exchcxx/xc_functional.hpp>
#include <map>

namespace ExchCXX {

template <typename T> class MapStr : public std::map<std::string, T> {
  std::map<std::string, T> str2val_;
  std::map<T, std::string> val2str_;

public:
  MapStr(std::initializer_list<std::pair<const std::string, T>> strval)
      : str2val_(strval) {

    for (auto &&v : str2val_) {
      val2str_.insert(std::make_pair(v.second, v.first));
    }
  }

  T to_val(std::string str) { return str2val_.at(str); }

  std::string to_str(T val) { return val2str_.at(val); }
};

MapStr<XCFunctional::Functional>
    FuncMap({{"SVWN3", XCFunctional::Functional::SVWN3},
             {"BLYP", XCFunctional::Functional::BLYP},
             {"PBE0", XCFunctional::Functional::PBE0},
             {"PBE_X", XCFunctional::Functional::PBE_X},
             {"SlaterExchange", XCFunctional::Functional::SlaterExchange}});

MapStr<Kernel> KerMap({{"SlaterExchange", Kernel::SlaterExchange},
                       {"PBE_X", Kernel::PBE_X},
                       {"PBE_C", Kernel::PBE_C},
                       {"LYP", Kernel::LYP},
                       {"B3LYP", Kernel::B3LYP},
                       {"PBE0", Kernel::PBE0},
                       {"VWN3", Kernel::VWN3},
                       {"VWN5", Kernel::VWN5},
                       {"PZ81", Kernel::PZ81},
                       {"PZ81_MOD", Kernel::PZ81_MOD},
                       {"PW91_LDA", Kernel::PW91_LDA},
                       {"PW91_LDA_MOD", Kernel::PW91_LDA_MOD},
                       {"PW91_LDA_RPA", Kernel::PW91_LDA_RPA},
                       {"B88", Kernel::B88}});

inline std::ostream &operator<<(std::ostream &out,
                                XCFunctional::Functional functional) {
  out << FuncMap.to_str(functional);
  return out;
}

} // namespace ExchCXX
