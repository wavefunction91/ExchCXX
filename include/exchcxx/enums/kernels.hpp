#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <iostream>
#include <exchcxx/bidirectional_map.hpp>

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

extern BidirectionalMap<std::string, Kernel> kernel_map;

std::ostream& operator<<( std::ostream& out, Kernel kern );

}
