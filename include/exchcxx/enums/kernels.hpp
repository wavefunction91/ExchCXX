#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

namespace ExchCXX {

enum class Kernel {
  // LDA Functionals
  SlaterExchange,
/*
  Wigner,
  RPA,
  HedinLundqvist,
  GunnarsonLundqvist,
  XAlpha,
*/
  VWN3,
  VWN5,
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

}
