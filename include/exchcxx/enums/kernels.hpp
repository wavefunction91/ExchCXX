#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <iostream>
#include <map>

namespace ExchCXX {

enum class Kernel {
  // LDA Functionals
  SlaterExchange,
  VWN3,
  VWN5,
  PZ81,
  PZ81_MOD,
  PW91_LDA,
  PW91_LDA_MOD,
  PW91_LDA_RPA,
/*
  Wigner,
  RPA,
  HedinLundqvist,
  GunnarsonLundqvist,
  XAlpha,
*/
/*
  VWN_RPA,
  PerdewZunger,
  PerdewZungerMod,
*/
  // GGA functionals
  PBE_X,
  PBE_C,
  B88,
  LYP,

  // Hybrid GGA functionals
  B3LYP,
  PBE0,
};

class KerMap {
  std::map<std::string, Kernel> str2val_;
  std::map<Kernel, std::string> val2str_;

public:
  KerMap()
      : str2val_({{"SlaterExchange", Kernel::SlaterExchange},
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
                  {"B88", Kernel::B88}}) {

    for (auto &&v : str2val_) {
      val2str_.insert(std::make_pair(v.second, v.first));
    }
  }

  Kernel to_val(std::string str) { return str2val_.at(str); }

  std::string to_str(Kernel val) { return val2str_.at(val); }
};

std::ostream& operator<<( std::ostream& out, Kernel kern );

}
