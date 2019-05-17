#ifndef __INCLUDED_KERNELS_KERNELS_HPP__
#define __INCLUDED_KERNELS_KERNELS_HPP__

// This file is meant to be included only from xc_kernel.hpp

enum Kernel {
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
  B88,
  LYP,

  // Hybrid GGA functionals
  B3LYP,
  PBE0,
};

#endif
