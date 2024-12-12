/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <exchcxx/util/bidirectional_map.hpp>
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
  revPBE_X,
  B88,
  LYP,

  // MGGA functionals
  SCAN_C,
  SCAN_X,
  SCANL_C,
  SCANL_X,
  R2SCAN_C,
  R2SCAN_X,
  R2SCANL_C,
  R2SCANL_X,
  FT98_X,
  M062X_X,
  M062X_C,
  PKZB_X,
  PKZB_C,
  

  // KEDFs
  PC07_K,
  PC07OPT_K,

  // Hybrid GGA functionals
  B3LYP,
  PBE0,

  // NEO LDA Functionals
  EPC17_1,
  EPC17_2,
  EPC18_1,
  EPC18_2,
};

inline static bool supports_unpolarized(ExchCXX::Kernel kern) {
  switch (kern) {
  case ExchCXX::Kernel::EPC17_1:
  case ExchCXX::Kernel::EPC17_2:
  case ExchCXX::Kernel::EPC18_1:
  case ExchCXX::Kernel::EPC18_2:
    return false;
  default:
    return true;
  }
}


extern BidirectionalMap<std::string, Kernel> kernel_map;

std::ostream& operator<<( std::ostream& out, Kernel kern );

}
