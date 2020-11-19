#pragma once

#include <map>

namespace ExchCXX {

enum class Functional {
  SVWN3,
  SVWN5,
  BLYP,
  B3LYP,
  PBE0,
};

class FuncMap {
  std::map<std::string, Functional> str2val_;
  std::map<Functional, std::string> val2str_;

public:
  FuncMap()
      : str2val_({{"SVWN3", Functional::SVWN3},
                  {"SVWN5", Functional::SVWN5},
                  {"BLYP", Functional::BLYP},
                  {"B3LYP", Functional::B3LYP},
                  {"PBE0", Functional::PBE0}}) {

    for (auto &&v : str2val_) {
      val2str_.insert(std::make_pair(v.second, v.first));
    }
  }

  Functional to_val(std::string str) { return str2val_.at(str); }

  std::string to_str(Functional val) { return val2str_.at(val); }
};

extern FuncMap str2func;

std::ostream &operator<<(std::ostream &out, Functional functional);

}
