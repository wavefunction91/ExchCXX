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


namespace ExchCXX {

class BuiltinKernel;

struct BuiltinSlaterExchange;
struct BuiltinVWN3;
struct BuiltinVWN_RPA;
struct BuiltinVWN;
struct BuiltinPW91_LDA;
struct BuiltinPW91_LDA_MOD;
struct BuiltinPW91_LDA_RPA;
struct BuiltinPZ81;
struct BuiltinPZ81_MOD;

struct BuiltinB88;
struct BuiltinLYP;
struct BuiltinPBE_X;
struct BuiltinRevPBE_X;
struct BuiltinPBE_C;
struct BuiltinB97_D;
struct BuiltinITYH_X;
struct BuiltinITYH_X_033;
struct BuiltinITYH_X_015;
struct BuiltinP86_C;
struct BuiltinP86VWN_FT_C;
struct BuiltinPW91_C;
struct BuiltinPBE_SOL_C;
struct BuiltinBMK_C;
struct BuiltinN12_C;
struct BuiltinN12_SX_C;
struct BuiltinSOGGA11_X_C;
struct BuiltinPW91_X;
struct BuiltinMPW91_X;
struct BuiltinOPTX_X;
struct BuiltinRPBE_X;
struct BuiltinSOGGA11_X_X;
struct BuiltinPW86_X;
struct BuiltinWB97_XC;
struct BuiltinWB97X_XC;
struct BuiltinWB97X_V_XC;
struct BuiltinWB97X_D_XC;
struct BuiltinWB97X_D3_XC;
struct BuiltinHJS_PBE_X;
struct BuiltinLCwPBE_wPBEh_X;
struct BuiltinLRCwPBE_HJS_PBE_X;
struct BuiltinLRCwPBEh_HJS_PBE_X;
struct BuiltinWPBEh_X_default0;
struct BuiltinHSE03_wPBEh_X;
struct BuiltinHSE06_wPBEh_X;


struct BuiltinSCAN_X;
struct BuiltinSCAN_C;
struct BuiltinR2SCAN_X;
struct BuiltinR2SCAN_C;
struct BuiltinM062X_X;
struct BuiltinM062X_C;
struct BuiltinPKZB_X;
struct BuiltinPKZB_C;
struct BuiltinFT98_X;
struct BuiltinTPSS_X;
struct BuiltinRevTPSS_X;
struct BuiltinM06_L_X;
struct BuiltinM06_X;
struct BuiltinM06_HF_X;
struct BuiltinRevM06_L_X;
struct BuiltinM06_SX_X;
struct BuiltinM06_L_C;
struct BuiltinM06_C;
struct BuiltinM06_HF_C;
struct BuiltinRevM06_L_C;
struct BuiltinM06_SX_C;
struct BuiltinM05_2X_C;
struct BuiltinM05_C;
struct BuiltinM08_HX_C;
struct BuiltinM08_SO_C;
struct BuiltinCF22D_C;
struct BuiltinM11_C;
struct BuiltinMN12_L_C;
struct BuiltinMN12_SX_C;
struct BuiltinMN15_C;
struct BuiltinMN15_L_C;
struct BuiltinTPSS_C;
struct BuiltinRevTPSS_C;
struct BuiltinRSCAN_C;
struct BuiltinBC95_C;
struct BuiltinMBEEF_X;
struct BuiltinRSCAN_X;
struct BuiltinBMK_X;
struct BuiltinM08_HX_X;
struct BuiltinM08_SO_X;
struct BuiltinMN12_L_X;
struct BuiltinMN15_L_X;
struct BuiltinMN15_X;
struct BuiltinCF22D_X;
struct BuiltinMN12_SX_X;
struct BuiltinM11_X;
struct BuiltinM05_X;
struct BuiltinM05_2X_X;


struct BuiltinPC07_K;
struct BuiltinPC07OPT_K;

template <typename XCEF, typename KEDF>
struct Deorbitalized;
struct BuiltinSCANL_X;
struct BuiltinSCANL_C;
struct BuiltinR2SCANL_X;
struct BuiltinR2SCANL_C;

struct BuiltinEPC17_1;
struct BuiltinEPC17_2;
struct BuiltinEPC18_1;
struct BuiltinEPC18_2;

template <typename KernelType>
struct kernel_traits;




namespace detail {
  template <typename KernelType, typename = void>
  struct BuiltinKernelImpl;

  class BuiltinKernelInterface;
}

}
