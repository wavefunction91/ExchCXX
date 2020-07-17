#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <iostream>

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
