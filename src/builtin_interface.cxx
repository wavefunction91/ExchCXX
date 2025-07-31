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

#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/interface.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>
#include <exchcxx/exceptions/exchcxx_exception.hpp>


namespace ExchCXX {
namespace detail  {





std::unique_ptr<BuiltinKernel> 
  gen_from_kern( Kernel kern, Spin polar ) {

  // Bail if polarized eval is requested and not supported
  EXCHCXX_BOOL_CHECK(kernel_map.key(kern) + " Needs to be Spin-Polarized!",
                     supports_unpolarized(kern) or polar == Spin::Polarized);

  if( kern == Kernel::SlaterExchange )
    return std::make_unique<BuiltinSlaterExchange>( polar );
  else if( kern == Kernel::VWN3 )
    return std::make_unique<BuiltinVWN3>( polar );
  else if( kern == Kernel::VWN5 )
    return std::make_unique<BuiltinVWN_RPA>( polar );
  else if( kern == Kernel::VWN )
    return std::make_unique<BuiltinVWN>( polar );
  else if( kern == Kernel::PW91_LDA )
    return std::make_unique<BuiltinPW91_LDA>( polar );
  else if( kern == Kernel::PW91_LDA_MOD )
    return std::make_unique<BuiltinPW91_LDA_MOD>( polar );
  else if( kern == Kernel::PW91_LDA_RPA )
    return std::make_unique<BuiltinPW91_LDA_RPA>( polar );
  else if( kern == Kernel::PZ81 )
    return std::make_unique<BuiltinPZ81>( polar );
  else if( kern == Kernel::PZ81_MOD )
    return std::make_unique<BuiltinPZ81_MOD>( polar );


  else if( kern == Kernel::B88 )
    return std::make_unique<BuiltinB88>( polar );
  else if( kern == Kernel::LYP )
    return std::make_unique<BuiltinLYP>( polar );
  else if( kern == Kernel::PBE_X )
    return std::make_unique<BuiltinPBE_X>( polar );
  else if( kern == Kernel::revPBE_X )
    return std::make_unique<BuiltinRevPBE_X>( polar );
  else if( kern == Kernel::PBE_C )
    return std::make_unique<BuiltinPBE_C>( polar );
  else if( kern == Kernel::B97_D )
    return std::make_unique<BuiltinB97_D>( polar );
  else if( kern == Kernel::ITYH_X )
    return std::make_unique<BuiltinITYH_X>( polar );
  else if( kern == Kernel::ITYH_X_033 )
    return std::make_unique<BuiltinITYH_X_033>( polar );
  else if( kern == Kernel::ITYH_X_015 )
    return std::make_unique<BuiltinITYH_X_015>( polar );
  else if( kern == Kernel::P86_C )
  return std::make_unique<BuiltinP86_C>( polar );
  else if( kern == Kernel::P86VWN_FT_C )
    return std::make_unique<BuiltinP86VWN_FT_C>( polar );
  else if( kern == Kernel::PW91_C )
    return std::make_unique<BuiltinPW91_C>( polar );
  else if( kern == Kernel::PBE_SOL_C )
    return std::make_unique<BuiltinPBE_SOL_C>( polar );
  else if( kern == Kernel::BMK_C )
    return std::make_unique<BuiltinBMK_C>( polar );
  else if( kern == Kernel::N12_C )
    return std::make_unique<BuiltinN12_C>( polar );
  else if( kern == Kernel::N12_SX_C )
    return std::make_unique<BuiltinN12_SX_C>( polar );
  else if( kern == Kernel::SOGGA11_X_C )
    return std::make_unique<BuiltinSOGGA11_X_C>( polar );
  else if( kern == Kernel::PW91_X )
    return std::make_unique<BuiltinPW91_X>( polar );
  else if( kern == Kernel::MPW91_X )
    return std::make_unique<BuiltinMPW91_X>( polar );
  else if( kern == Kernel::OPTX_X )
    return std::make_unique<BuiltinOPTX_X>( polar );
  else if( kern == Kernel::RPBE_X )
    return std::make_unique<BuiltinRPBE_X>( polar );
  else if( kern == Kernel::SOGGA11_X_X )
    return std::make_unique<BuiltinSOGGA11_X_X>( polar );
  else if( kern == Kernel::PW86_X )
    return std::make_unique<BuiltinPW86_X>( polar );
  else if( kern == Kernel::wB97_XC )
    return std::make_unique<BuiltinWB97_XC>( polar );
  else if( kern == Kernel::wB97X_XC )
    return std::make_unique<BuiltinWB97X_XC>( polar );
  else if( kern == Kernel::wB97X_V_XC )
    return std::make_unique<BuiltinWB97X_V_XC>( polar );
  else if( kern == Kernel::wB97X_D_XC )
    return std::make_unique<BuiltinWB97X_D_XC>( polar );
  else if( kern == Kernel::wB97X_D3_XC )
    return std::make_unique<BuiltinWB97X_D3_XC>( polar );
  else if( kern == Kernel::HJS_PBE_X )
    return std::make_unique<BuiltinHJS_PBE_X>( polar );
  else if( kern == Kernel::LCwPBE_wPBEh_X )
    return std::make_unique<BuiltinLCwPBE_wPBEh_X>( polar );
  else if( kern == Kernel::LRCwPBE_HJS_PBE_X )
    return std::make_unique<BuiltinLRCwPBE_HJS_PBE_X>( polar );
  else if( kern == Kernel::LRCwPBEh_HJS_PBE_X )
    return std::make_unique<BuiltinLRCwPBEh_HJS_PBE_X>( polar );
  else if( kern == Kernel::wPBEh_X_default0 )
    return std::make_unique<BuiltinWPBEh_X_default0>( polar );
  else if( kern == Kernel::HSE03_wPBEh_X )
    return std::make_unique<BuiltinHSE03_wPBEh_X>( polar );
  else if( kern == Kernel::HSE06_wPBEh_X )
    return std::make_unique<BuiltinHSE06_wPBEh_X>( polar );

  else if( kern == Kernel::SCAN_X )
    return std::make_unique<BuiltinSCAN_X>( polar );
  else if( kern == Kernel::SCAN_C )
    return std::make_unique<BuiltinSCAN_C>( polar );
  else if( kern == Kernel::SCANL_C )
    return std::make_unique<BuiltinSCANL_C>( polar );
  else if( kern == Kernel::SCANL_X )
    return std::make_unique<BuiltinSCANL_X>( polar );
  else if( kern == Kernel::R2SCAN_X )
    return std::make_unique<BuiltinR2SCAN_X>( polar );
  else if( kern == Kernel::R2SCAN_C )
    return std::make_unique<BuiltinR2SCAN_C>( polar );
  else if( kern == Kernel::R2SCANL_X )
    return std::make_unique<BuiltinR2SCANL_X>( polar );
  else if( kern == Kernel::R2SCANL_C )
    return std::make_unique<BuiltinR2SCANL_C>( polar );
  else if( kern == Kernel::FT98_X )
    return std::make_unique<BuiltinFT98_X>( polar );
  else if( kern == Kernel::M062X_X )
    return std::make_unique<BuiltinM062X_X>( polar );
  else if( kern == Kernel::M062X_C )
    return std::make_unique<BuiltinM062X_C>( polar );
  else if( kern == Kernel::PKZB_X )
    return std::make_unique<BuiltinPKZB_X>( polar );
  else if( kern == Kernel::PKZB_C )
    return std::make_unique<BuiltinPKZB_C>( polar );
  else if( kern == Kernel::TPSS_X )
    return std::make_unique<BuiltinTPSS_X>( polar );
  else if( kern == Kernel::revTPSS_X )
    return std::make_unique<BuiltinRevTPSS_X>( polar );
  else if( kern == Kernel::M06_L_X )
    return std::make_unique<BuiltinM06_L_X>( polar );
  else if( kern == Kernel::M06_X )
    return std::make_unique<BuiltinM06_X>( polar );
  else if( kern == Kernel::M06_HF_X )
    return std::make_unique<BuiltinM06_HF_X>( polar );
  else if( kern == Kernel::revM06_L_X )
    return std::make_unique<BuiltinRevM06_L_X>( polar );
  else if( kern == Kernel::M06_SX_X )
    return std::make_unique<BuiltinM06_SX_X>( polar );
    else if( kern == Kernel::M06_L_C )
    return std::make_unique<BuiltinM06_L_C>( polar );
  else if( kern == Kernel::M06_C )
    return std::make_unique<BuiltinM06_C>( polar );
  else if( kern == Kernel::M06_HF_C )
    return std::make_unique<BuiltinM06_HF_C>( polar );
  else if( kern == Kernel::revM06_L_C )
    return std::make_unique<BuiltinRevM06_L_C>( polar );
  else if( kern == Kernel::M06_SX_C )
    return std::make_unique<BuiltinM06_SX_C>( polar );
  else if( kern == Kernel::M05_2X_C )
    return std::make_unique<BuiltinM05_2X_C>( polar );
  else if( kern == Kernel::M05_C )
    return std::make_unique<BuiltinM05_C>( polar );
  else if( kern == Kernel::M08_HX_C )
    return std::make_unique<BuiltinM08_HX_C>( polar );
  else if( kern == Kernel::M08_SO_C )
    return std::make_unique<BuiltinM08_SO_C>( polar );
  else if( kern == Kernel::CF22D_C )
    return std::make_unique<BuiltinCF22D_C>( polar );
  else if( kern == Kernel::M11_C )
    return std::make_unique<BuiltinM11_C>( polar );
  else if( kern == Kernel::MN12_L_C )
    return std::make_unique<BuiltinMN12_L_C>( polar );
  else if( kern == Kernel::MN12_SX_C )
    return std::make_unique<BuiltinMN12_SX_C>( polar );
  else if( kern == Kernel::MN15_C )
    return std::make_unique<BuiltinMN15_C>( polar );
  else if( kern == Kernel::MN15_L_C )
    return std::make_unique<BuiltinMN15_L_C>( polar );
  else if( kern == Kernel::TPSS_C )
    return std::make_unique<BuiltinTPSS_C>( polar );
  else if( kern == Kernel::revTPSS_C )
    return std::make_unique<BuiltinRevTPSS_C>( polar );
  else if( kern == Kernel::RSCAN_C )
    return std::make_unique<BuiltinRSCAN_C>( polar );
  else if( kern == Kernel::BC95_C )
    return std::make_unique<BuiltinBC95_C>( polar );
  else if( kern == Kernel::MN12_L_C )
    return std::make_unique<BuiltinMN12_L_C>( polar );
  else if( kern == Kernel::MN12_SX_C )
    return std::make_unique<BuiltinMN12_SX_C>( polar );
  else if( kern == Kernel::MN15_C )
    return std::make_unique<BuiltinMN15_C>( polar );
  else if( kern == Kernel::MN15_L_C )
    return std::make_unique<BuiltinMN15_L_C>( polar );
  else if( kern == Kernel::M06_L_C )
    return std::make_unique<BuiltinM06_L_C>( polar );
  else if( kern == Kernel::M06_C )
    return std::make_unique<BuiltinM06_C>( polar );
  else if( kern == Kernel::mBEEF_X )
    return std::make_unique<BuiltinMBEEF_X>( polar );
  else if( kern == Kernel::RSCAN_X )
    return std::make_unique<BuiltinRSCAN_X>( polar );
  else if( kern == Kernel::BMK_X )
    return std::make_unique<BuiltinBMK_X>( polar );
  else if( kern == Kernel::M08_HX_X )
    return std::make_unique<BuiltinM08_HX_X>( polar );
  else if( kern == Kernel::M08_SO_X )
    return std::make_unique<BuiltinM08_SO_X>( polar );
  else if( kern == Kernel::MN12_L_X )
    return std::make_unique<BuiltinMN12_L_X>( polar );
  else if( kern == Kernel::MN15_L_X )
    return std::make_unique<BuiltinMN15_L_X>( polar );
  else if( kern == Kernel::MN15_X )
    return std::make_unique<BuiltinMN15_X>( polar );
  else if( kern == Kernel::CF22D_X )
    return std::make_unique<BuiltinCF22D_X>( polar );
  else if( kern == Kernel::MN12_SX_X )
    return std::make_unique<BuiltinMN12_SX_X>( polar );
  else if( kern == Kernel::M11_X )
    return std::make_unique<BuiltinM11_X>( polar );
  else if( kern == Kernel::M05_X )
    return std::make_unique<BuiltinM05_X>( polar );
  else if( kern == Kernel::M05_2X_X )
    return std::make_unique<BuiltinM05_2X_X>( polar );
  

  else if( kern == Kernel::PC07_K )
    return std::make_unique<BuiltinPC07_K>( polar );
  else if( kern == Kernel::PC07OPT_K )
    return std::make_unique<BuiltinPC07OPT_K>( polar );
  
  else if( kern == Kernel::EPC17_1) {
    return std::make_unique<BuiltinEPC17_1>( polar );
  } else if( kern == Kernel::EPC17_2) {
    return std::make_unique<BuiltinEPC17_2>( polar );
  } else if( kern == Kernel::EPC18_1) {
    return std::make_unique<BuiltinEPC18_1>( polar );
  } else if( kern == Kernel::EPC18_2) {
    return std::make_unique<BuiltinEPC18_2>( polar );

  } else
    throw std::runtime_error("Specified kernel does not have a builtin implementation");

}

BuiltinKernelInterface::~BuiltinKernelInterface() noexcept = default;

BuiltinKernelInterface::BuiltinKernelInterface( Kernel kern, 
  Spin polar ) : 
    whatami_(kern), impl_(gen_from_kern( kern, polar )) { }

BuiltinKernelInterface::BuiltinKernelInterface( 
  const BuiltinKernelInterface& other 
) : BuiltinKernelInterface( other.whatami_, other.impl_->polar() ) { }

std::unique_ptr<XCKernelImpl> BuiltinKernelInterface::clone_() const {
  return std::make_unique<BuiltinKernelInterface>( *this );
}







bool BuiltinKernelInterface::is_lda_()       const noexcept {
  return impl_->is_lda();
}
bool BuiltinKernelInterface::is_gga_()       const noexcept {
  return impl_->is_gga();
}
bool BuiltinKernelInterface::is_mgga_()      const noexcept {
  return impl_->is_mgga();
}
bool BuiltinKernelInterface::is_polarized_() const noexcept {
  return impl_->is_polarized();
}
bool BuiltinKernelInterface::needs_laplacian_()       const noexcept {
  return impl_->needs_laplacian();
}
bool BuiltinKernelInterface::needs_tau_()       const noexcept {
  return impl_->needs_tau();
}
bool BuiltinKernelInterface::is_epc_()       const noexcept {
  return impl_->is_epc();
}

bool BuiltinKernelInterface::supports_inc_interface_() const noexcept {
  return true;
}


#define FORWARD_FOR_BUILTIN( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_DEVICE( APPROX, TYPE, func ) \
  FORWARD_XC_ARGS_DEVICE( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_INC( APPROX, TYPE, func ) \
  FORWARD_XC_INC_ARGS( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )
#define FORWARD_FOR_BUILTIN_INC_DEVICE( APPROX, TYPE, func ) \
  FORWARD_XC_INC_ARGS_DEVICE( APPROX, TYPE, BuiltinKernelInterface:: func ## _, \
                   impl_->func, const )



FORWARD_FOR_BUILTIN( LDA,  EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( LDA,  EXC_VXC, eval_exc_vxc )
FORWARD_FOR_BUILTIN( LDA,  FXC,     eval_fxc     )
FORWARD_FOR_BUILTIN( LDA,  VXC_FXC, eval_vxc_fxc )
FORWARD_FOR_BUILTIN( GGA,  EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( GGA,  EXC_VXC, eval_exc_vxc )
FORWARD_FOR_BUILTIN( GGA,  FXC,     eval_fxc     )
FORWARD_FOR_BUILTIN( GGA,  VXC_FXC, eval_vxc_fxc )
FORWARD_FOR_BUILTIN( MGGA, EXC,     eval_exc     )
FORWARD_FOR_BUILTIN( MGGA, EXC_VXC, eval_exc_vxc )
FORWARD_FOR_BUILTIN( MGGA, FXC,     eval_fxc     )
FORWARD_FOR_BUILTIN( MGGA, VXC_FXC, eval_vxc_fxc )

FORWARD_FOR_BUILTIN_INC( LDA,  EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( LDA,  EXC_VXC, eval_exc_vxc_inc )
FORWARD_FOR_BUILTIN_INC( LDA,  FXC,     eval_fxc_inc     )
FORWARD_FOR_BUILTIN_INC( LDA,  VXC_FXC, eval_vxc_fxc_inc )
FORWARD_FOR_BUILTIN_INC( GGA,  EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( GGA,  EXC_VXC, eval_exc_vxc_inc )
FORWARD_FOR_BUILTIN_INC( GGA,  FXC,     eval_fxc_inc     )
FORWARD_FOR_BUILTIN_INC( GGA,  VXC_FXC, eval_vxc_fxc_inc )
FORWARD_FOR_BUILTIN_INC( MGGA, EXC,     eval_exc_inc     )
FORWARD_FOR_BUILTIN_INC( MGGA, EXC_VXC, eval_exc_vxc_inc )
FORWARD_FOR_BUILTIN_INC( MGGA, FXC,     eval_fxc_inc     )
FORWARD_FOR_BUILTIN_INC( MGGA, VXC_FXC, eval_vxc_fxc_inc )


#ifdef EXCHCXX_ENABLE_DEVICE

FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( LDA,  EXC_VXC, eval_exc_vxc_device )
FORWARD_FOR_BUILTIN_DEVICE( LDA,  FXC,     eval_fxc_device     )
FORWARD_FOR_BUILTIN_DEVICE( LDA,  VXC_FXC, eval_vxc_fxc_device )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  EXC_VXC, eval_exc_vxc_device )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  FXC,     eval_fxc_device     )
FORWARD_FOR_BUILTIN_DEVICE( GGA,  VXC_FXC, eval_vxc_fxc_device )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC,     eval_exc_device     )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_device )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, FXC,     eval_fxc_device     )
FORWARD_FOR_BUILTIN_DEVICE( MGGA, VXC_FXC, eval_vxc_fxc_device )


FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  EXC_VXC, eval_exc_vxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  FXC,     eval_fxc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( LDA,  VXC_FXC, eval_vxc_fxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  EXC_VXC, eval_exc_vxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  FXC,     eval_fxc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( GGA,  VXC_FXC, eval_vxc_fxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, EXC,     eval_exc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_inc_device )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, FXC,     eval_fxc_inc_device     )
FORWARD_FOR_BUILTIN_INC_DEVICE( MGGA, VXC_FXC, eval_vxc_fxc_inc_device )

#endif















}
}

