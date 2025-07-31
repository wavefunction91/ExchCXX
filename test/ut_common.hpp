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

#include "catch2/catch_all.hpp"
#include <exchcxx/exchcxx.hpp>
#include <xc.h>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <random>

#include "reference_values.hpp"

using Catch::Approx;

enum class TestInterface {
  EXC,
  EXC_VXC,

  EXC_INC,
  EXC_VXC_INC,
  
  FXC,
  VXC_FXC,

  FXC_INC,
  VXC_FXC_INC,
};

enum class EvalType {
  Regular,
  Small,
  Zero
};


static std::vector<ExchCXX::Kernel> lda_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5,
  ExchCXX::Kernel::VWN,
  ExchCXX::Kernel::PZ81,
  ExchCXX::Kernel::PZ81_MOD,
  ExchCXX::Kernel::PW91_LDA,
  ExchCXX::Kernel::PW91_LDA_MOD,
  ExchCXX::Kernel::PW91_LDA_RPA
};

static std::vector<ExchCXX::Kernel> gga_kernels = {
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,
  ExchCXX::Kernel::revPBE_X,
  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::B97_D,
  ExchCXX::Kernel::ITYH_X,
  ExchCXX::Kernel::P86_C,
  ExchCXX::Kernel::P86VWN_FT_C,
  ExchCXX::Kernel::PW91_C,
  ExchCXX::Kernel::PBE_SOL_C,
  ExchCXX::Kernel::BMK_C,
  ExchCXX::Kernel::N12_C,
  ExchCXX::Kernel::N12_SX_C,
  ExchCXX::Kernel::SOGGA11_X_C,
  ExchCXX::Kernel::PW91_X,
  ExchCXX::Kernel::MPW91_X,
  ExchCXX::Kernel::OPTX_X,
  ExchCXX::Kernel::RPBE_X,
  ExchCXX::Kernel::SOGGA11_X_X,
  ExchCXX::Kernel::PW86_X,
  ExchCXX::Kernel::wB97_XC,
  ExchCXX::Kernel::wB97X_XC,
  ExchCXX::Kernel::wB97X_V_XC,
  ExchCXX::Kernel::wB97X_D_XC,
  ExchCXX::Kernel::wB97X_D3_XC,
  ExchCXX::Kernel::HJS_PBE_X,
  ExchCXX::Kernel::wPBEh_X_default0,

};

static std::vector<ExchCXX::Kernel> mgga_kernels = {
  ExchCXX::Kernel::SCAN_X,
  ExchCXX::Kernel::SCAN_C,
  ExchCXX::Kernel::SCANL_C,
  ExchCXX::Kernel::R2SCAN_X,
  ExchCXX::Kernel::R2SCAN_C,
  ExchCXX::Kernel::R2SCANL_X,
  ExchCXX::Kernel::R2SCANL_C,
  ExchCXX::Kernel::M062X_X,
  ExchCXX::Kernel::M062X_C,
  ExchCXX::Kernel::FT98_X,
  ExchCXX::Kernel::PC07_K,
  ExchCXX::Kernel::PC07OPT_K,
  ExchCXX::Kernel::TPSS_X,
  ExchCXX::Kernel::revTPSS_X,
  ExchCXX::Kernel::M06_L_X,
  ExchCXX::Kernel::M06_X,
  ExchCXX::Kernel::revM06_L_X,
  ExchCXX::Kernel::M06_HF_X,
  ExchCXX::Kernel::M06_SX_X,
  ExchCXX::Kernel::M06_L_C,
  ExchCXX::Kernel::M06_C,
  ExchCXX::Kernel::revM06_L_C,
  ExchCXX::Kernel::M06_HF_C,
  ExchCXX::Kernel::M06_SX_C,
  ExchCXX::Kernel::M05_2X_C,
  ExchCXX::Kernel::M05_C,
  ExchCXX::Kernel::M08_HX_C,
  ExchCXX::Kernel::M08_SO_C,
  ExchCXX::Kernel::CF22D_C,
  ExchCXX::Kernel::M11_C,
  ExchCXX::Kernel::MN12_L_C,
  ExchCXX::Kernel::MN12_SX_C,
  ExchCXX::Kernel::MN15_C,
  ExchCXX::Kernel::MN15_L_C,
  ExchCXX::Kernel::TPSS_C,
  ExchCXX::Kernel::revTPSS_C,
  ExchCXX::Kernel::BC95_C,
  ExchCXX::Kernel::BMK_X,
  ExchCXX::Kernel::M08_HX_X,
  ExchCXX::Kernel::M08_SO_X,
  ExchCXX::Kernel::MN12_L_X,
  ExchCXX::Kernel::MN15_L_X,
  ExchCXX::Kernel::MN15_X,
  ExchCXX::Kernel::CF22D_X,
  ExchCXX::Kernel::MN12_SX_X,
  ExchCXX::Kernel::M11_X,
  ExchCXX::Kernel::M05_X,
  ExchCXX::Kernel::M05_2X_X,
};

static std::vector<ExchCXX::Kernel> epc_lda_kernels = {
  ExchCXX::Kernel::EPC17_1,
  ExchCXX::Kernel::EPC17_2,
  ExchCXX::Kernel::EPC18_1,
  ExchCXX::Kernel::EPC18_2
};

static std::vector<ExchCXX::Kernel> builtin_supported_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5,
  ExchCXX::Kernel::VWN,
  ExchCXX::Kernel::PZ81,
  ExchCXX::Kernel::PZ81_MOD,
  ExchCXX::Kernel::PW91_LDA,
  ExchCXX::Kernel::PW91_LDA_MOD,
  ExchCXX::Kernel::PW91_LDA_RPA,

  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,
  ExchCXX::Kernel::revPBE_X,
  ExchCXX::Kernel::B97_D,
  ExchCXX::Kernel::ITYH_X,
  ExchCXX::Kernel::P86_C,
  ExchCXX::Kernel::P86VWN_FT_C,
  ExchCXX::Kernel::PW91_C,
  ExchCXX::Kernel::PBE_SOL_C,
  ExchCXX::Kernel::BMK_C,
  ExchCXX::Kernel::N12_C,
  ExchCXX::Kernel::N12_SX_C,
  ExchCXX::Kernel::SOGGA11_X_C,
  ExchCXX::Kernel::PW91_X,
  ExchCXX::Kernel::MPW91_X,
  ExchCXX::Kernel::OPTX_X,
  ExchCXX::Kernel::RPBE_X,
  ExchCXX::Kernel::SOGGA11_X_X,
  ExchCXX::Kernel::PW86_X,
  ExchCXX::Kernel::wB97_XC,
  ExchCXX::Kernel::wB97X_XC,
  ExchCXX::Kernel::wB97X_V_XC,
  ExchCXX::Kernel::wB97X_D_XC,
  ExchCXX::Kernel::wB97X_D3_XC,
  ExchCXX::Kernel::HJS_PBE_X,
  ExchCXX::Kernel::wPBEh_X_default0,


  ExchCXX::Kernel::SCAN_X, 
  ExchCXX::Kernel::SCAN_C,
  ExchCXX::Kernel::SCANL_C,
  ExchCXX::Kernel::SCANL_X,
  ExchCXX::Kernel::R2SCAN_X, 
  ExchCXX::Kernel::R2SCAN_C,
  ExchCXX::Kernel::R2SCANL_X, 
  ExchCXX::Kernel::R2SCANL_C,
  ExchCXX::Kernel::FT98_X,
  ExchCXX::Kernel::M062X_X,
  ExchCXX::Kernel::M062X_C,
  ExchCXX::Kernel::PKZB_X,
  ExchCXX::Kernel::PKZB_C,
  ExchCXX::Kernel::TPSS_X,
  ExchCXX::Kernel::revTPSS_X,
  ExchCXX::Kernel::M06_L_X,
  ExchCXX::Kernel::M06_X,
  ExchCXX::Kernel::revM06_L_X,
  ExchCXX::Kernel::M06_HF_X,
  ExchCXX::Kernel::M06_SX_X,
  ExchCXX::Kernel::M06_L_C,
  ExchCXX::Kernel::M06_C,
  ExchCXX::Kernel::revM06_L_C,
  ExchCXX::Kernel::M06_HF_C,
  ExchCXX::Kernel::M06_SX_C,
  ExchCXX::Kernel::M05_2X_C,
  ExchCXX::Kernel::M05_C,
  ExchCXX::Kernel::M08_HX_C,
  ExchCXX::Kernel::M08_SO_C,
  ExchCXX::Kernel::CF22D_C,
  ExchCXX::Kernel::M11_C,
  ExchCXX::Kernel::MN12_L_C,
  ExchCXX::Kernel::MN12_SX_C,
  ExchCXX::Kernel::MN15_C,
  ExchCXX::Kernel::MN15_L_C,
  ExchCXX::Kernel::TPSS_C,
  ExchCXX::Kernel::revTPSS_C,
  ExchCXX::Kernel::BC95_C,
  ExchCXX::Kernel::BMK_X,
  ExchCXX::Kernel::M08_HX_X,
  ExchCXX::Kernel::M08_SO_X,
  ExchCXX::Kernel::MN12_L_X,
  ExchCXX::Kernel::MN15_L_X,
  ExchCXX::Kernel::MN15_X,
  ExchCXX::Kernel::CF22D_X,
  ExchCXX::Kernel::MN12_SX_X,
  ExchCXX::Kernel::M11_X,
  ExchCXX::Kernel::M05_X,
  ExchCXX::Kernel::M05_2X_X,

  ExchCXX::Kernel::PC07_K,
  ExchCXX::Kernel::PC07OPT_K,

  ExchCXX::Kernel::EPC17_1,
  ExchCXX::Kernel::EPC17_2,
  ExchCXX::Kernel::EPC18_1,
  ExchCXX::Kernel::EPC18_2
};

static std::vector<ExchCXX::Kernel> deorbitalized_kernels = {

  ExchCXX::Kernel::SCANL_C,
  ExchCXX::Kernel::SCANL_X,
  ExchCXX::Kernel::R2SCANL_X,
  ExchCXX::Kernel::R2SCANL_C
};

static std::vector<ExchCXX::Kernel> unstable_small_kernels = {

  // Unstable on CPU and GPU
  ExchCXX::Kernel::SCANL_C, // vrho, vsigma (GPU only)
  ExchCXX::Kernel::ITYH_X, //vsigma
  ExchCXX::Kernel::PKZB_C, // exc, vrho, vsigma
  ExchCXX::Kernel::TPSS_C, // exc, vrho, vsigma
  ExchCXX::Kernel::revTPSS_C, // vrho, vsigma
  ExchCXX::Kernel::RSCAN_C,
  ExchCXX::Kernel::BMK_X,
  ExchCXX::Kernel::M05_2X_C,
  ExchCXX::Kernel::M05_C,
  ExchCXX::Kernel::PC07_K,
  ExchCXX::Kernel::PC07OPT_K,
  ExchCXX::Kernel::SOGGA11_X_X,
  
  // The following kernels is unstable on GPU only
  // vsigma unstable only and report -0.0 == Approx( -0.0 ): B97 series, BMK_C, N12_C, OPTX_X
  ExchCXX::Kernel::B97_D,
  ExchCXX::Kernel::wB97_XC,
  ExchCXX::Kernel::wB97X_XC,
  ExchCXX::Kernel::wB97X_V_XC,
  ExchCXX::Kernel::wB97X_D3_XC,
  ExchCXX::Kernel::wB97X_D_XC,
  ExchCXX::Kernel::BMK_C,
  ExchCXX::Kernel::N12_C,
  ExchCXX::Kernel::OPTX_X,
  // vrho and vsigma unstable and report -0.0 == Approx( -0.0 ): PW91_X, MPW91_X, OPTX_X
  ExchCXX::Kernel::PW91_X,
  ExchCXX::Kernel::MPW91_X,

  ExchCXX::Kernel::SCANL_X, // vrho, vsigma, vrho sometimes big error
  ExchCXX::Kernel::mBEEF_X, // vrho, vsigma, vtau, vrho huge error

  // vsigma unstable only with small error
  ExchCXX::Kernel::MN12_L_C,
  ExchCXX::Kernel::MN15_C,
  ExchCXX::Kernel::M08_HX_C,
  ExchCXX::Kernel::CF22D_C,
  ExchCXX::Kernel::N12_SX_C,
  ExchCXX::Kernel::MN12_SX_C,
  ExchCXX::Kernel::M11_C,
  ExchCXX::Kernel::M08_SO_C,
  ExchCXX::Kernel::MN15_L_C
};

static std::vector<ExchCXX::Kernel> unstable_small_polarized_fxc_exchange_due_to_libxc_bug = {

  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::M05_2X_X,
  ExchCXX::Kernel::M05_X,
  ExchCXX::Kernel::M062X_X,
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::RPBE_X,
  ExchCXX::Kernel::revPBE_X,
  ExchCXX::Kernel::PW86_X,
  ExchCXX::Kernel::wPBEh_X_default0,
  ExchCXX::Kernel::SCAN_X,
  ExchCXX::Kernel::PKZB_X,
  ExchCXX::Kernel::TPSS_X,
  ExchCXX::Kernel::revTPSS_X,
  ExchCXX::Kernel::M06_L_X,
  ExchCXX::Kernel::M06_X,
  ExchCXX::Kernel::revM06_L_X,
  ExchCXX::Kernel::M06_HF_X,
  ExchCXX::Kernel::M06_SX_X,
  ExchCXX::Kernel::BMK_X,
  ExchCXX::Kernel::M08_HX_X,
  ExchCXX::Kernel::M08_SO_X,
  ExchCXX::Kernel::MN12_L_X,
  ExchCXX::Kernel::MN15_L_X,
  ExchCXX::Kernel::MN15_X,
  ExchCXX::Kernel::CF22D_X,
  ExchCXX::Kernel::MN12_SX_X,
};

static std::vector<ExchCXX::Kernel> unstable_small_2nd_deriv_device ={

  ExchCXX::Kernel::B97_D,
  ExchCXX::Kernel::BMK_C,
  ExchCXX::Kernel::BMK_X,
  ExchCXX::Kernel::N12_C,
  ExchCXX::Kernel::N12_SX_C,
  ExchCXX::Kernel::OPTX_X,
  ExchCXX::Kernel::PW91_X,
  ExchCXX::Kernel::MPW91_X,
  ExchCXX::Kernel::wB97_XC,
  ExchCXX::Kernel::wB97X_XC,
  ExchCXX::Kernel::wB97X_V_XC,
  ExchCXX::Kernel::wB97X_D3_XC,
  ExchCXX::Kernel::wB97X_D_XC,
  ExchCXX::Kernel::M08_HX_C,
  ExchCXX::Kernel::M08_SO_C,
  ExchCXX::Kernel::CF22D_C,
  ExchCXX::Kernel::M11_C,
  ExchCXX::Kernel::MN12_L_C,
  ExchCXX::Kernel::MN12_SX_C,
  ExchCXX::Kernel::MN15_C,
  ExchCXX::Kernel::MN15_L_C,
  ExchCXX::Kernel::TPSS_C, 
  ExchCXX::Kernel::PC07_K, 
  ExchCXX::Kernel::PC07OPT_K, 
  ExchCXX::Kernel::mBEEF_X, // vrho, vsigma, vtau, vrho huge error

  //Polarized unstable only
  ExchCXX::Kernel::PKZB_C, 
  ExchCXX::Kernel::revTPSS_C, 
  ExchCXX::Kernel::M05_2X_C,
  ExchCXX::Kernel::M05_C,

};

inline bool is_unstable_small(ExchCXX::Kernel kern) {
  return std::find(unstable_small_kernels.begin(), unstable_small_kernels.end(),
                   kern) != unstable_small_kernels.end();
}
inline bool is_unstable_small_polarized_fxc_exchange_due_to_libxc_bug(ExchCXX::Kernel kern) {
  return std::find(unstable_small_polarized_fxc_exchange_due_to_libxc_bug.begin(), unstable_small_polarized_fxc_exchange_due_to_libxc_bug.end(),
                   kern) != unstable_small_polarized_fxc_exchange_due_to_libxc_bug.end();
}
inline bool is_unstable_small_2nd_deriv_device(ExchCXX::Kernel kern) {
  return std::find(unstable_small_2nd_deriv_device.begin(), unstable_small_2nd_deriv_device.end(),
                   kern) != unstable_small_2nd_deriv_device.end();
}
inline bool is_deorbitalized(ExchCXX::Kernel kern) {
  return std::find(deorbitalized_kernels.begin(), deorbitalized_kernels.end(),
                   kern) != deorbitalized_kernels.end();
}
inline bool is_epc(ExchCXX::Kernel kern) {
  return std::find(epc_lda_kernels.begin(), epc_lda_kernels.end(),
                   kern) != epc_lda_kernels.end();
}

static constexpr std::array string_kernal_pairs = {
    std::pair("SlaterExchange", ExchCXX::Kernel::SlaterExchange),
    std::pair("PBE_X",ExchCXX::Kernel::PBE_X),
    std::pair("revPBE_X", ExchCXX::Kernel::revPBE_X),
    std::pair("PBE_C", ExchCXX::Kernel::PBE_C),
    std::pair("B97_D", ExchCXX::Kernel::B97_D),
    std::pair("ITYH_X", ExchCXX::Kernel::ITYH_X),
    std::pair("P86_C", ExchCXX::Kernel::P86_C),
    std::pair("P86VWN_FT_C", ExchCXX::Kernel::P86VWN_FT_C),
    std::pair("PW91_C", ExchCXX::Kernel::PW91_C),
    std::pair("PBE_SOL_C", ExchCXX::Kernel::PBE_SOL_C),
    std::pair("BMK_C", ExchCXX::Kernel::BMK_C),
    std::pair("N12_C", ExchCXX::Kernel::N12_C),
    std::pair("N12_SX_C", ExchCXX::Kernel::N12_SX_C),
    std::pair("SOGGA11_X_C", ExchCXX::Kernel::SOGGA11_X_C),
    std::pair("PW91_X", ExchCXX::Kernel::PW91_X),
    std::pair("MPW91_X", ExchCXX::Kernel::MPW91_X),
    std::pair("OPTX_X", ExchCXX::Kernel::OPTX_X),
    std::pair("RPBE_X", ExchCXX::Kernel::RPBE_X),
    std::pair("SOGGA11_X_X", ExchCXX::Kernel::SOGGA11_X_X),
    std::pair("PW86_X", ExchCXX::Kernel::PW86_X),
    std::pair("wB97_XC", ExchCXX::Kernel::wB97_XC),
    std::pair("wB97X_XC", ExchCXX::Kernel::wB97X_XC),
    std::pair("wB97X_V_XC", ExchCXX::Kernel::wB97X_V_XC),
    std::pair("wB97X_D_XC", ExchCXX::Kernel::wB97X_D_XC),
    std::pair("wB97X_D3_XC", ExchCXX::Kernel::wB97X_D3_XC),
    std::pair("HJS_PBE_X", ExchCXX::Kernel::HJS_PBE_X),
    std::pair("wPBEh_X_default0", ExchCXX::Kernel::wPBEh_X_default0),

    std::pair("SCAN_X",ExchCXX::Kernel::SCAN_X),
    std::pair("SCAN_C", ExchCXX::Kernel::SCAN_C),
    std::pair("SCANL_C", ExchCXX::Kernel::SCANL_C),
    std::pair("FT98_X",ExchCXX::Kernel::FT98_X),
    std::pair("PC07_K",ExchCXX::Kernel::PC07_K),
    std::pair("PC07OPT_K",ExchCXX::Kernel::PC07OPT_K),
    std::pair("R2SCANL_X",ExchCXX::Kernel::R2SCANL_X),
    std::pair("R2SCANL_C", ExchCXX::Kernel::R2SCANL_C),
    std::pair("R2SCAN_X",ExchCXX::Kernel::R2SCAN_X),
    std::pair("R2SCAN_C", ExchCXX::Kernel::R2SCAN_C),
    std::pair("M062X_X",ExchCXX::Kernel::M062X_X),
    std::pair("M062X_C", ExchCXX::Kernel::M062X_C),
    std::pair("PKZB_X",ExchCXX::Kernel::PKZB_X),
    std::pair("PKZB_C", ExchCXX::Kernel::PKZB_C),
    std::pair("TPSS_X", ExchCXX::Kernel::TPSS_X),
    std::pair("revTPSS_X", ExchCXX::Kernel::revTPSS_X),
    std::pair("M06_L_X", ExchCXX::Kernel::M06_L_X),
    std::pair("M06_X", ExchCXX::Kernel::M06_X),
    std::pair("revM06_L_X", ExchCXX::Kernel::revM06_L_X),
    std::pair("M06_HF_X", ExchCXX::Kernel::M06_HF_X),
    std::pair("M06_SX_X", ExchCXX::Kernel::M06_SX_X),
    std::pair("M06_L_C", ExchCXX::Kernel::M06_L_C),
    std::pair("M06_C", ExchCXX::Kernel::M06_C),
    std::pair("revM06_L_C", ExchCXX::Kernel::revM06_L_C),
    std::pair("M06_HF_C", ExchCXX::Kernel::M06_HF_C),
    std::pair("M06_SX_C", ExchCXX::Kernel::M06_SX_C),
    std::pair("M05_2X_C", ExchCXX::Kernel::M05_2X_C),
    std::pair("M05_C", ExchCXX::Kernel::M05_C),
    std::pair("M08_HX_C", ExchCXX::Kernel::M08_HX_C),
    std::pair("M08_SO_C", ExchCXX::Kernel::M08_SO_C),
    std::pair("CF22D_C", ExchCXX::Kernel::CF22D_C),
    std::pair("M11_C", ExchCXX::Kernel::M11_C),
    std::pair("MN12_L_C", ExchCXX::Kernel::MN12_L_C),
    std::pair("MN12_SX_C", ExchCXX::Kernel::MN12_SX_C),
    std::pair("MN15_C", ExchCXX::Kernel::MN15_C),
    std::pair("MN15_L_C", ExchCXX::Kernel::MN15_L_C),
    std::pair("TPSS_C", ExchCXX::Kernel::TPSS_C),
    std::pair("revTPSS_C", ExchCXX::Kernel::revTPSS_C),
    std::pair("RSCAN_C", ExchCXX::Kernel::RSCAN_C),
    std::pair("BC95_C", ExchCXX::Kernel::BC95_C),
    std::pair("mBEEF_X", ExchCXX::Kernel::mBEEF_X),
    std::pair("RSCAN_X", ExchCXX::Kernel::RSCAN_X),
    std::pair("BMK_X", ExchCXX::Kernel::BMK_X),
    std::pair("M08_HX_X", ExchCXX::Kernel::M08_HX_X),
    std::pair("M08_SO_X", ExchCXX::Kernel::M08_SO_X),
    std::pair("MN12_L_X", ExchCXX::Kernel::MN12_L_X),
    std::pair("MN15_L_X", ExchCXX::Kernel::MN15_L_X),
    std::pair("MN15_X", ExchCXX::Kernel::MN15_X),
    std::pair("CF22D_X", ExchCXX::Kernel::CF22D_X),
    std::pair("MN12_SX_X", ExchCXX::Kernel::MN12_SX_X),
    std::pair("M11_X", ExchCXX::Kernel::M11_X),
    std::pair("M05_X", ExchCXX::Kernel::M05_X),
    std::pair("M05_2X_X", ExchCXX::Kernel::M05_2X_X),

    std::pair("LYP", ExchCXX::Kernel::LYP),
    std::pair("VWN3", ExchCXX::Kernel::VWN3),
    std::pair("VWN5", ExchCXX::Kernel::VWN5),
    std::pair("VWN", ExchCXX::Kernel::VWN),
    std::pair("PZ81", ExchCXX::Kernel::PZ81),
    std::pair("PZ81_MOD", ExchCXX::Kernel::PZ81_MOD),
    std::pair("PW91_LDA", ExchCXX::Kernel::PW91_LDA),
    std::pair("PW91_LDA_MOD", ExchCXX::Kernel::PW91_LDA_MOD),
    std::pair("PW91_LDA_RPA", ExchCXX::Kernel::PW91_LDA_RPA),
    std::pair("B88", ExchCXX::Kernel::B88),
    std::pair("EPC17_1", ExchCXX::Kernel::EPC17_1),
    std::pair("EPC17_2", ExchCXX::Kernel::EPC17_2),
    std::pair("EPC18_1", ExchCXX::Kernel::EPC18_1),
    std::pair("EPC18_2", ExchCXX::Kernel::EPC18_2)
};

static constexpr std::array string_functional_pairs = {
    std::pair("SVWN3", ExchCXX::Functional::SVWN3),
    std::pair("SVWN5", ExchCXX::Functional::SVWN5),
    std::pair("BLYP", ExchCXX::Functional::BLYP),
    std::pair("B3LYP", ExchCXX::Functional::B3LYP),
    std::pair("PBE", ExchCXX::Functional::PBE),
    std::pair("SCAN", ExchCXX::Functional::SCAN),
    std::pair("R2SCANL", ExchCXX::Functional::R2SCANL),
    std::pair("M062X", ExchCXX::Functional::M062X),
    std::pair("PKZB", ExchCXX::Functional::PKZB),
    std::pair("PBE0", ExchCXX::Functional::PBE0),
    std::pair("EPC17_1", ExchCXX::Functional::EPC17_1),
    std::pair("EPC17_2", ExchCXX::Functional::EPC17_2),
    std::pair("EPC18_1", ExchCXX::Functional::EPC18_1),
    std::pair("EPC18_2", ExchCXX::Functional::EPC18_2),
    
    std::pair("B97D", ExchCXX::Functional::B97D),
    std::pair("B97D3ZERO", ExchCXX::Functional::B97D3ZERO),
    std::pair("CAMB3LYP", ExchCXX::Functional::CAMB3LYP),
    std::pair("LDA", ExchCXX::Functional::LDA),
    std::pair("M06L", ExchCXX::Functional::M06L),
    std::pair("SCAN0", ExchCXX::Functional::SCAN0),
    std::pair("SPW92", ExchCXX::Functional::SPW92),
    std::pair("TPSS", ExchCXX::Functional::TPSS),
    std::pair("TPSSH", ExchCXX::Functional::TPSSh),
    std::pair("TPSS0", ExchCXX::Functional::TPSS0),
    std::pair("VWN3", ExchCXX::Functional::VWN3),
    std::pair("VWN5", ExchCXX::Functional::VWN5),
    std::pair("LRCWPBE", ExchCXX::Functional::LRCwPBE),
    std::pair("LRCWPBEH", ExchCXX::Functional::LRCwPBEh),
    std::pair("BP86", ExchCXX::Functional::BP86),
    std::pair("HSE03", ExchCXX::Functional::HSE03),
    std::pair("HSE06", ExchCXX::Functional::HSE06),
    std::pair("REVB3LYP", ExchCXX::Functional::revB3LYP),
    std::pair("REVPBE0", ExchCXX::Functional::revPBE0),
    std::pair("REVTPSS", ExchCXX::Functional::revTPSS),
    std::pair("REVTPSSH", ExchCXX::Functional::revTPSSh),
    std::pair("PW91", ExchCXX::Functional::PW91),
    std::pair("MBEEF", ExchCXX::Functional::mBEEF),
    std::pair("B3PW91", ExchCXX::Functional::B3PW91),
    std::pair("O3LYP", ExchCXX::Functional::O3LYP),
    std::pair("OLYP", ExchCXX::Functional::OLYP),
    std::pair("OPBE", ExchCXX::Functional::OPBE),
    std::pair("MPW1K", ExchCXX::Functional::MPW1K),
    std::pair("RPBE", ExchCXX::Functional::RPBE),
    std::pair("B88", ExchCXX::Functional::B88),
    std::pair("MPW91", ExchCXX::Functional::MPW91),
    std::pair("RSCAN", ExchCXX::Functional::RSCAN),
    std::pair("TUNEDCAMB3LYP", ExchCXX::Functional::TUNEDCAMB3LYP),
    std::pair("WB97", ExchCXX::Functional::wB97),
    std::pair("WB97X", ExchCXX::Functional::wB97X),
    std::pair("WB97XD", ExchCXX::Functional::wB97XD),
    std::pair("WB97XD3", ExchCXX::Functional::wB97XD3),
    std::pair("LCWPBE", ExchCXX::Functional::LCwPBE),
    std::pair("X3LYP", ExchCXX::Functional::X3LYP),
    std::pair("XLYP", ExchCXX::Functional::XLYP),
    std::pair("BHANDH", ExchCXX::Functional::BHANDH),
    std::pair("BMK", ExchCXX::Functional::BMK),
    std::pair("BP86VWN", ExchCXX::Functional::BP86VWN),
    std::pair("PW86B95", ExchCXX::Functional::PW86B95),
    std::pair("PW86PBE", ExchCXX::Functional::PW86PBE),
    std::pair("R2SCAN0", ExchCXX::Functional::R2SCAN0),
    std::pair("R2SCANH", ExchCXX::Functional::R2SCANh),
    std::pair("R2SCAN50", ExchCXX::Functional::R2SCAN50),
    std::pair("M05", ExchCXX::Functional::M05),
    std::pair("M06", ExchCXX::Functional::M06),
    std::pair("M08HX", ExchCXX::Functional::M08HX),
    std::pair("M08SO", ExchCXX::Functional::M08SO),
    std::pair("M052X", ExchCXX::Functional::M052X),
    std::pair("M06SX", ExchCXX::Functional::M06SX),
    std::pair("CF22D", ExchCXX::Functional::CF22D),
    std::pair("SOGGA11X", ExchCXX::Functional::SOGGA11X),
    std::pair("M06HF", ExchCXX::Functional::M06HF),
    std::pair("M11", ExchCXX::Functional::M11),
    std::pair("MN12L", ExchCXX::Functional::MN12L),
    std::pair("MN12SX", ExchCXX::Functional::MN12SX),
    std::pair("MN15", ExchCXX::Functional::MN15),
    std::pair("MN15L", ExchCXX::Functional::MN15L),
    std::pair("REVM06L", ExchCXX::Functional::revM06L),
};
