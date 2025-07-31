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

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX
{

  BidirectionalMap<std::string, Kernel> kernel_map{
      {{"SlaterExchange", Kernel::SlaterExchange},
       {"PBE_X", Kernel::PBE_X},
       {"PBE_C", Kernel::PBE_C},
       {"B97_D", Kernel::B97_D},
       {"ITYH_X", Kernel::ITYH_X},
       {"ITYH_X", Kernel::ITYH_X_033},
       {"ITYH_X", Kernel::ITYH_X_015},
       {"wB97_XC", Kernel::wB97_XC},
       {"wB97X_XC", Kernel::wB97X_XC},
       {"wB97X_V_XC", Kernel::wB97X_V_XC},
       {"wB97X_D_XC", Kernel::wB97X_D_XC},
       {"wB97X_D3_XC", Kernel::wB97X_D3_XC},

       {"SCAN_X", Kernel::SCAN_X},
       {"SCAN_C", Kernel::SCAN_C},
       {"SCANL_C", Kernel::SCANL_C},
       {"SCANL_X", Kernel::SCANL_X},
       {"FT98_X", Kernel::FT98_X},
       {"PC07_K", Kernel::PC07_K},
       {"PC07OPT_K", Kernel::PC07OPT_K},
       {"R2SCAN_X", Kernel::R2SCAN_X},
       {"R2SCAN_C", Kernel::R2SCAN_C},
       {"R2SCANL_X", Kernel::R2SCANL_X},
       {"R2SCANL_C", Kernel::R2SCANL_C},
       {"M062X_X", Kernel::M062X_X},
       {"M062X_C", Kernel::M062X_C},
       {"PKZB_X", Kernel::PKZB_X},
       {"PKZB_C", Kernel::PKZB_C},
       {"TPSS_X", Kernel::TPSS_X},
       {"revTPSS_X", Kernel::revTPSS_X},
       {"M06_L_X", Kernel::M06_L_X},
       {"M06_X", Kernel::M06_X},
       {"revM06_L_X", Kernel::revM06_L_X},
       {"M06_HF_X", Kernel::M06_HF_X},
       {"M06_SX_X", Kernel::M06_SX_X},
       {"M06_L_C", Kernel::M06_L_C},
       {"M06_C", Kernel::M06_C},
       {"revM06_L_C", Kernel::revM06_L_C},
       {"M06_HF_C", Kernel::M06_HF_C},
       {"M06_SX_C", Kernel::M06_SX_C},
       {"M05_2X_C", Kernel::M05_2X_C},
       {"M05_C", Kernel::M05_C},
       {"M08_HX_C", Kernel::M08_HX_C},
       {"M08_SO_C", Kernel::M08_SO_C},
       {"CF22D_C", Kernel::CF22D_C},
       {"M11_C", Kernel::M11_C},
       {"MN12_L_C", Kernel::MN12_L_C},
       {"MN12_SX_C", Kernel::MN12_SX_C},
       {"MN15_C", Kernel::MN15_C},
       {"MN15_L_C", Kernel::MN15_L_C},
       {"M05_2X_C", Kernel::M05_2X_C},
       {"M05_C", Kernel::M05_C},
       {"M08_HX_C", Kernel::M08_HX_C},
       {"M08_SO_C", Kernel::M08_SO_C},
       {"CF22D_C", Kernel::CF22D_C},
       {"M11_C", Kernel::M11_C},
       {"MN12_L_C", Kernel::MN12_L_C},
       {"MN12_SX_C", Kernel::MN12_SX_C},
       {"MN15_C", Kernel::MN15_C},
       {"MN15_L_C", Kernel::MN15_L_C},
       {"P86_C", Kernel::P86_C},
       {"P86VWN_FT_C", Kernel::P86VWN_FT_C},
       {"PW91_C", Kernel::PW91_C},
       {"PBE_SOL_C", Kernel::PBE_SOL_C},
       {"BMK_C", Kernel::BMK_C},
       {"N12_C", Kernel::N12_C},
       {"N12_SX_C", Kernel::N12_SX_C},
       {"SOGGA11_X_C", Kernel::SOGGA11_X_C},
       {"TPSS_C", Kernel::TPSS_C},
       {"revTPSS_C", Kernel::revTPSS_C},
       {"RSCAN_C", Kernel::RSCAN_C},
       {"BC95_C", Kernel::BC95_C},
       {"PW91_X", Kernel::PW91_X},
       {"MPW91_X", Kernel::MPW91_X},
       {"OPTX_X", Kernel::OPTX_X},
       {"RPBE_X", Kernel::RPBE_X},
       {"SOGGA11_X_X", Kernel::SOGGA11_X_X},
       {"PW86_X", Kernel::PW86_X},
       {"HJS_PBE_X", Kernel::HJS_PBE_X},
       {"LCwPBE_wPBEh_X", Kernel::LCwPBE_wPBEh_X},
       {"LRCwPBE_HJS_PBE_X", Kernel::LRCwPBE_HJS_PBE_X},
       {"LRCwPBEh_HJS_PBE_X", Kernel::LRCwPBEh_HJS_PBE_X},
       {"wPBEh_X_default0", Kernel::wPBEh_X_default0},
       {"HSE03_wPBEh_X", Kernel::HSE03_wPBEh_X},
       {"HSE06_wPBEh_X", Kernel::HSE06_wPBEh_X},

       {"mBEEF_X", Kernel::mBEEF_X},
       {"RSCAN_X", Kernel::RSCAN_X},
       {"BMK_X", Kernel::BMK_X},
       {"M08_HX_X", Kernel::M08_HX_X},
       {"M08_SO_X", Kernel::M08_SO_X},
       {"MN12_L_X", Kernel::MN12_L_X},
       {"MN15_L_X", Kernel::MN15_L_X},
       {"MN15_X", Kernel::MN15_X},
       {"CF22D_X", Kernel::CF22D_X},
       {"MN12_SX_X", Kernel::MN12_SX_X},
       {"M11_X", Kernel::M11_X},
       {"M05_X", Kernel::M05_X},
       {"M05_2X_X", Kernel::M05_2X_X},

       {"revPBE_X", Kernel::revPBE_X},
       {"LYP", Kernel::LYP},
       {"VWN3", Kernel::VWN3},
       {"VWN5", Kernel::VWN5},
       {"VWN", Kernel::VWN},
       {"PZ81", Kernel::PZ81},
       {"PZ81_MOD", Kernel::PZ81_MOD},
       {"PW91_LDA", Kernel::PW91_LDA},
       {"PW91_LDA_MOD", Kernel::PW91_LDA_MOD},
       {"PW91_LDA_RPA", Kernel::PW91_LDA_RPA},
       {"B88", Kernel::B88},
       {"EPC17_1", Kernel::EPC17_1},
       {"EPC17_2", Kernel::EPC17_2},
       {"EPC18_1", Kernel::EPC18_1},
       {"EPC18_2", Kernel::EPC18_2}}};

  std::ostream &operator<<(std::ostream &out, Kernel kern)
  {
    out << kernel_map.key(kern);
    return out;
  }

  XCKernel::XCKernel(
      const Backend backend,
      const Kernel kern,
      const Spin polar) : XCKernel(kernel_factory(backend, kern, polar)) {}

#ifdef EXCHCXX_ENABLE_LIBXC
  XCKernel::XCKernel(
      const libxc_name_string &xc_name,
      const Spin polar) : XCKernel(libxc_kernel_factory(xc_name.get(), polar)) {}
#endif

  XCKernel::XCKernel(impl_ptr &&ptr) : pimpl_(std::move(ptr)) {}

  XCKernel::XCKernel(const XCKernel &other) : pimpl_(other.pimpl_->clone()) {}

  XCKernel::XCKernel(XCKernel &&other) noexcept = default;

  XCKernel &XCKernel::operator=(XCKernel &&other) noexcept = default;

  XCKernel &XCKernel::operator=(const XCKernel &other)
  {
    return *this = XCKernel(other);
  }

  XCKernel::~XCKernel() noexcept = default;

}
