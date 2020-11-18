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

class KerMap : public std::map<std::string, Kernel> {
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

inline std::ostream& operator<<( std::ostream& out, Kernel kern ) {

  switch(kern) {
    case Kernel::SlaterExchange:
      out << "SlaterExchange";
      break;
    case Kernel::VWN3:
      out << "VWN3";
      break;
    case Kernel::VWN5:
      out << "VWN5";
      break;
    case Kernel::PBE_X:
      out << "PBE_X";
      break;
    case Kernel::PBE_C:
      out << "PBE_C";
      break;
    case Kernel::B88:
      out << "B88";
      break;
    case Kernel::LYP:
      out << "LYP";
      break;
    case Kernel::B3LYP:
      out << "B3LYP";
      break;
    case Kernel::PBE0:
      out << "PBE0";
      break;
    default:
      out << "Unknown Kernel";
  }
  return out;

}


}
