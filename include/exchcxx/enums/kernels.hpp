/**
 * ExchCXX 
 *
 * Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). 
 *
 * Portions Copyright (c) Microsoft Corporation.
 *
 * All rights reserved.
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
  VWN,
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
  B97_D,
  ITYH_X,
  ITYH_X_033,
  ITYH_X_015,
  P86_C,
  P86VWN_FT_C,
  PW91_C,
  PBE_SOL_C,
  BMK_C,
  N12_C,
  N12_SX_C,
  SOGGA11_X_C,
  PW91_X,
  MPW91_X,
  OPTX_X,
  RPBE_X,
  SOGGA11_X_X,
  PW86_X,
  HJS_PBE_X,
  LCwPBE_wPBEh_X,
  LRCwPBE_HJS_PBE_X,
  LRCwPBEh_HJS_PBE_X,
  wPBEh_X_default0,
  HSE03_wPBEh_X,
  HSE06_wPBEh_X,


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
  TPSS_X,
  revTPSS_X,
  M06_L_X,
  M06_X,
  revM06_L_X,
  M06_HF_X,
  M06_SX_X,
  M06_L_C,
  M06_C,
  revM06_L_C,
  M06_HF_C,
  M06_SX_C,
  M05_2X_C,
  M05_C,
  M08_HX_C,
  M08_SO_C,
  CF22D_C,
  M11_C,
  MN12_L_C,
  MN12_SX_C,
  MN15_C,
  MN15_L_C,
  TPSS_C,
  revTPSS_C,
  RSCAN_C,
  BC95_C,
  mBEEF_X,
  RSCAN_X,
  BMK_X,
  M08_HX_X,
  M08_SO_X,
  MN12_L_X,
  MN15_L_X,
  MN15_X,
  CF22D_X,
  MN12_SX_X,
  M11_X,
  M05_X,
  M05_2X_X,
  wB97_XC,
  wB97X_XC,
  wB97X_V_XC,
  wB97X_D_XC,
  wB97X_D3_XC,


  // KEDFs
  PC07_K,
  PC07OPT_K,


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
