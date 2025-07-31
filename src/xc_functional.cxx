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

#include <exchcxx/xc_functional.hpp>

namespace ExchCXX {


BidirectionalMap<std::string, Functional> functional_map{{
    {"SVWN3", Functional::SVWN3},
    {"SVWN5", Functional::SVWN5},
    {"BLYP", Functional::BLYP},
    {"B3LYP", Functional::B3LYP},
    {"PBE", Functional::PBE},
    {"SCAN", Functional::SCAN},
    {"R2SCAN", Functional::R2SCAN},
    {"R2SCANL", Functional::R2SCANL},
    {"M062X", Functional::M062X},
    {"PKZB", Functional::PKZB},
    {"REVPBE", Functional::revPBE},
    {"PBE0", Functional::PBE0},
    {"EPC17_1", Functional::EPC17_1},
    {"EPC17_2", Functional::EPC17_2},
    {"EPC18_1", Functional::EPC18_1},
    {"EPC18_2", Functional::EPC18_2},
    
    {"B97D", Functional::B97D},
    {"B97D3ZERO", Functional::B97D3ZERO},
    {"CAMB3LYP", Functional::CAMB3LYP},
    {"LDA", Functional::LDA},
    {"M06L", Functional::M06L},
    {"SCAN0", Functional::SCAN0},
    {"SPW92", Functional::SPW92},
    {"TPSS", Functional::TPSS},
    {"TPSSH", Functional::TPSSh},
    {"TPSS0", Functional::TPSS0},
    {"VWN3", Functional::VWN3},
    {"VWN5", Functional::VWN5},
    {"LRCWPBE", Functional::LRCwPBE},
    {"LRCWPBEH", Functional::LRCwPBEh},
    {"BP86", Functional::BP86},
    {"HSE03", Functional::HSE03},
    {"HSE06", Functional::HSE06},
    {"REVB3LYP", Functional::revB3LYP},
    {"REVPBE0", Functional::revPBE0},
    {"REVTPSS", Functional::revTPSS},
    {"REVTPSSH", Functional::revTPSSh},
    {"PW91", Functional::PW91},
    {"MBEEF", Functional::mBEEF},
    {"B3PW91", Functional::B3PW91},
    {"O3LYP", Functional::O3LYP},
    {"OLYP", Functional::OLYP},
    {"OPBE", Functional::OPBE},
    {"MPW1K", Functional::MPW1K},
    {"RPBE", Functional::RPBE},
    {"B88", Functional::B88},
    {"MPW91", Functional::MPW91},
    {"RSCAN", Functional::RSCAN},
    {"TUNEDCAMB3LYP", Functional::TUNEDCAMB3LYP},
    {"WB97", Functional::wB97},
    {"WB97X", Functional::wB97X},
    {"WB97XD", Functional::wB97XD},
    {"WB97XD3", Functional::wB97XD3},
    {"LCWPBE", Functional::LCwPBE},
    {"X3LYP", Functional::X3LYP},
    {"XLYP", Functional::XLYP},
    {"BHANDH", Functional::BHANDH},
    {"BMK", Functional::BMK},
    {"BP86VWN", Functional::BP86VWN},
    {"PW86B95", Functional::PW86B95},
    {"PW86PBE", Functional::PW86PBE},
    {"R2SCAN0", Functional::R2SCAN0},
    {"R2SCANH", Functional::R2SCANh},
    {"R2SCAN50", Functional::R2SCAN50},
    {"M05", Functional::M05},
    {"M06", Functional::M06},
    {"M08HX", Functional::M08HX},
    {"M08SO", Functional::M08SO},
    {"M052X", Functional::M052X},
    {"M06SX", Functional::M06SX},
    {"CF22D", Functional::CF22D},
    {"SOGGA11X", Functional::SOGGA11X},
    {"M06HF", Functional::M06HF},
    {"M11", Functional::M11},
    {"MN12L", Functional::MN12L},
    {"MN12SX", Functional::MN12SX},
    {"MN15", Functional::MN15},
    {"MN15L", Functional::MN15L},
    {"REVM06L", Functional::revM06L}}
  };

std::ostream &operator<<(std::ostream &out, Functional functional) {
  out << functional_map.key(functional);
  return out;
}

std::pair<std::vector<std::pair<double, XCKernel>>, HybCoeffs> functional_factory( 
  const Backend        backend,
  const Functional func,
  const Spin          polar
) {

  std::vector< std::pair<double, XCKernel> > kerns;
  HybCoeffs hyb_coefs = {0.0, 0.0, 0.0}; // Default initialization

  if( func == Functional::SVWN3 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {1.0, XCKernel( backend, Kernel::VWN3,           polar )}
    };
  else if( func == Functional::SVWN5 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {1.0, XCKernel( backend, Kernel::VWN5,           polar )}
    };
  else if( func == Functional::BLYP )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B88, polar )},
      {1.0, XCKernel( backend, Kernel::LYP, polar )}
    };
  else if( func == Functional::B3LYP ){
    hyb_coefs = {0.2, 0.0, 0.0};
    kerns = { 
      {0.08, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.72, XCKernel( backend, Kernel::B88,           polar )},
      {0.19, XCKernel( backend, Kernel::VWN5,          polar )},
      {0.81, XCKernel( backend, Kernel::LYP,           polar )}
    };
  }
  else if( func == Functional::PBE )
    kerns = {
        {1.0, XCKernel( backend, Kernel::PBE_X, polar )},
        {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  else if( func == Functional::SCAN )
    kerns = {
        {1.0, XCKernel( backend, Kernel::SCAN_X, polar )},
        {1.0, XCKernel( backend, Kernel::SCAN_C, polar )}
    };
  else if( func == Functional::R2SCAN )
    kerns = {
        {1.0, XCKernel( backend, Kernel::R2SCAN_X, polar )},
        {1.0, XCKernel( backend, Kernel::R2SCAN_C, polar )}
    };
  else if( func == Functional::PKZB )
    kerns = {
        {1.0, XCKernel( backend, Kernel::PKZB_X, polar )},
        {1.0, XCKernel( backend, Kernel::PKZB_C, polar )}
    };
  else if( func == Functional::M062X ){
    hyb_coefs = {0.54, 0.0, 0.0};
    kerns = {
        {1.0, XCKernel( backend, Kernel::M062X_X, polar )},
        {1.0, XCKernel( backend, Kernel::M062X_C, polar )}
    };
  }
  else if( func == Functional::R2SCANL )
    kerns = {
        {1.0, XCKernel( backend, Kernel::R2SCANL_X, polar )},
        {1.0, XCKernel( backend, Kernel::R2SCANL_C, polar )}
    };
  else if( func == Functional::revPBE )
    kerns = {
        {1.0, XCKernel( backend, Kernel::revPBE_X, polar )},
        {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  else if( func == Functional::PBE0 ){
    hyb_coefs = {0.25, 0.0, 0.0};
    kerns = { {0.75, XCKernel( backend, Kernel::PBE_X, polar )},
              {1.0, XCKernel( backend, Kernel::PBE_C, polar )} };
  }
  else if( func == Functional::EPC17_1 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::EPC17_1, polar )}
    };
  else if( func == Functional::EPC17_2 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::EPC17_2, polar )}
    };
  else if( func == Functional::EPC18_1 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::EPC18_1, polar )}
    };
  else if( func == Functional::EPC18_2 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::EPC18_2, polar )}
    };
  else if( func == Functional::B97D )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B97_D, polar )}
    };
  else if( func == Functional::B97D3ZERO )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B97_D, polar )}
    };
  else if( func == Functional::CAMB3LYP ){
    hyb_coefs = {0.65, -0.46, 0.33};
    kerns = {
      {0.35, XCKernel( backend, Kernel::B88, polar)},
      {0.46, XCKernel( backend, Kernel::ITYH_X_033, polar)},
      {0.19, XCKernel( backend, Kernel::VWN, polar)},
      {0.81, XCKernel( backend, Kernel::LYP, polar)}
    };
  }
  else if( func == Functional::LDA )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::SlaterExchange, polar )}
    };
  else if( func == Functional::M06L )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M06_L_X, polar )},
      {1.0, XCKernel( backend, Kernel::M06_L_C, polar )}
    };
  else if( func == Functional::SCAN0 ){
    hyb_coefs = {0.25, 0.0, 0.0};
    kerns = { 
      {0.75, XCKernel( backend, Kernel::SCAN_X, polar )},
      {1.0, XCKernel( backend, Kernel::SCAN_C, polar )}
    };
  }
  else if( func == Functional::SPW92 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {1.0, XCKernel( backend, Kernel::PW91_LDA, polar )}
    };
  else if( func == Functional::TPSS )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::TPSS_X, polar )},
      {1.0, XCKernel( backend, Kernel::TPSS_C, polar )}
    };
  else if( func == Functional::TPSSh ){
    hyb_coefs = {0.1, 0.0, 0.0};
    kerns = { 
      {0.9, XCKernel( backend, Kernel::TPSS_X, polar )},
      {1.0, XCKernel( backend, Kernel::TPSS_C, polar )}
    };
  }
  else if( func == Functional::TPSS0 ){
    hyb_coefs = {0.25, 0.0, 0.0};
    kerns = { 
      {0.75, XCKernel( backend, Kernel::TPSS_X, polar )},
      {1.0, XCKernel( backend, Kernel::TPSS_C, polar )}
    };
  }
  else if( func == Functional::VWN3 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::VWN3,           polar )}
    };
  else if( func == Functional::VWN5 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::VWN5,           polar )}
    };
  else if( func == Functional::LRCwPBE ){
    hyb_coefs = {1.0, -1.0, 0.3};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::LRCwPBE_HJS_PBE_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if( func == Functional::LRCwPBEh ){
    hyb_coefs = {1.0, -0.8, 0.2};
    kerns = { 
      {0.8, XCKernel( backend, Kernel::LRCwPBEh_HJS_PBE_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if( func == Functional::BP86 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B88, polar )},
      {1.0, XCKernel( backend, Kernel::P86_C, polar )}
    };
  else if( func == Functional::HSE03 ){
    hyb_coefs = {0.0, 0.25, 0.10606601717798213};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wPBEh_X_default0, polar )},
      {-0.25, XCKernel( backend, Kernel::HSE03_wPBEh_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if( func == Functional::HSE06 ){
    hyb_coefs = {0.0, 0.25, 0.11};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wPBEh_X_default0, polar )},
      {-0.25, XCKernel( backend, Kernel::HSE06_wPBEh_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if( func == Functional::revB3LYP ){
    hyb_coefs = {0.2, 0.0, 0.0};
    kerns = { 
      {0.13, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.67, XCKernel( backend, Kernel::B88,           polar )},
      {0.16, XCKernel( backend, Kernel::VWN5,          polar )},
      {0.84, XCKernel( backend, Kernel::LYP,           polar )}
    };
  }
  else if( func == Functional::revPBE0 ){
    hyb_coefs = {0.25, 0.0, 0.0};
    kerns = { 
      {0.75, XCKernel( backend, Kernel::revPBE_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if( func == Functional::revTPSS )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::revTPSS_X, polar )},
      {1.0, XCKernel( backend, Kernel::revTPSS_C, polar )}
    };
  else if( func == Functional::revTPSSh ){
    hyb_coefs = {0.1, 0.0, 0.0};
    kerns = { 
      {0.9, XCKernel( backend, Kernel::revTPSS_X, polar )},
      {1.0, XCKernel( backend, Kernel::revTPSS_C, polar )}
    };
  }
  else if ( func == Functional::PW91 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::PW91_X, polar )},
      {1.0, XCKernel( backend, Kernel::PW91_C, polar )}
    };
  else if ( func == Functional::mBEEF )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::mBEEF_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_SOL_C, polar )}
    };
  else if ( func == Functional::B3PW91 ){
    hyb_coefs = {0.2, 0.0, 0.0};
    // 0.08, 0.72, 0.19, 0.81
    kerns = { 
      {0.08, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.72, XCKernel( backend, Kernel::B88, polar )},
      {0.19, XCKernel( backend, Kernel::PW91_LDA, polar )},
      {0.81, XCKernel( backend, Kernel::PW91_C, polar )}
    };
  }
  else if ( func == Functional::O3LYP ){
    hyb_coefs = {0.1161, 0.0, 0.0};
    kerns = { 
      {0.071006917, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.8133, XCKernel( backend, Kernel::OPTX_X,           polar )},
      {0.19, XCKernel( backend, Kernel::VWN,          polar )},
      {0.81, XCKernel( backend, Kernel::LYP,          polar )}
    };
  }
  else if ( func == Functional::OLYP )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::OPTX_X, polar )},
      {1.0, XCKernel( backend, Kernel::LYP,           polar )}
    };
  else if ( func == Functional::OPBE )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::OPTX_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C,           polar )}
    };
  else if ( func == Functional::MPW1K ){
    hyb_coefs = {0.428, 0.0, 0.0};
    kerns = { 
      {0.572, XCKernel( backend, Kernel::MPW91_X, polar )},
      {1.0, XCKernel( backend, Kernel::PW91_C, polar )}
    };
  }
  else if ( func == Functional::RPBE )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::RPBE_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  else if ( func == Functional::B88 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B88, polar )}
    };
  else if ( func == Functional::MPW91 )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::MPW91_X, polar )},
      {1.0, XCKernel( backend, Kernel::PW91_C, polar )}
    };
  else if ( func == Functional::RSCAN )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::RSCAN_X, polar )},
      {1.0, XCKernel( backend, Kernel::RSCAN_C, polar )}
    };
  else if ( func == Functional::TUNEDCAMB3LYP ){
    hyb_coefs = {1.0, -0.9201, 0.15};
    kerns = {
      {0.9201, XCKernel( backend, Kernel::ITYH_X_015, polar)},
      {0.19, XCKernel( backend, Kernel::VWN, polar)},
      {0.81, XCKernel( backend, Kernel::LYP, polar)}
    };
  }
  else if ( func == Functional::wB97 ){
    hyb_coefs = {1.0, -1.0, 0.4};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wB97_XC, polar )},
    };
  }
  else if ( func == Functional::wB97X ){
    hyb_coefs = {1.0, -0.842294, 0.3};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wB97X_XC, polar )},
    };
  }
  else if ( func == Functional::wB97XD ){
    hyb_coefs = {1.0, -0.777964, 0.2};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wB97X_D_XC, polar )},
    };
  }
  else if ( func == Functional::wB97XD3 ){
    hyb_coefs = {1.0, -0.804272, 0.25};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::wB97X_D3_XC, polar )},
    };
  }
  else if ( func == Functional::LCwPBE ){
    hyb_coefs = {1.0, -1.0, 0.4};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::LCwPBE_wPBEh_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if ( func == Functional::X3LYP ){
    hyb_coefs = {0.218, 0.0, 0.0};
    kerns = { 
      {0.073, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.166615, XCKernel( backend, Kernel::PW91_X, polar )},
      {0.542385, XCKernel( backend, Kernel::B88, polar )},
      {0.129, XCKernel( backend, Kernel::VWN5, polar )},
      {0.871, XCKernel( backend, Kernel::LYP, polar )}
    };
  }
  else if ( func == Functional::XLYP ){
    kerns = { 
      {-0.069, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {0.722, XCKernel( backend, Kernel::B88, polar )},
      {0.347, XCKernel( backend, Kernel::PW91_X, polar )},
      {1.0, XCKernel( backend, Kernel::LYP, polar )}
    };
  }
  else if ( func == Functional::BHANDH ){
    hyb_coefs = {0.5, 0.0, 0.0};
    kerns = { 
      {0.5, XCKernel( backend, Kernel::SlaterExchange, polar )},
      {1.0, XCKernel( backend, Kernel::LYP, polar )}
    };  
  }
  else if ( func == Functional::BMK ){
    hyb_coefs = {0.42, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::BMK_X, polar )},
      {1.0, XCKernel( backend, Kernel::BMK_C, polar )}
    };  
  }
  else if ( func == Functional::BP86VWN ){
    kerns = { 
      {1.0, XCKernel( backend, Kernel::B88, polar )},
      {1.0, XCKernel( backend, Kernel::P86VWN_FT_C, polar )}
    };
  }
  else if ( func == Functional::PW86B95 ){
    hyb_coefs = {0.29, 0.0, 0.0};
    kerns = { 
      {0.71, XCKernel( backend, Kernel::PW86_X, polar )},
      {1.0, XCKernel( backend, Kernel::BC95_C, polar )}
    };
  }
  else if ( func == Functional::PW86PBE ){
    kerns = { 
      {1.0, XCKernel( backend, Kernel::PW86_X, polar )},
      {1.0, XCKernel( backend, Kernel::PBE_C, polar )}
    };
  }
  else if ( func == Functional::R2SCAN0 ){
    hyb_coefs = {0.25, 0.0, 0.0};
    kerns = { 
      {0.75, XCKernel( backend, Kernel::R2SCAN_X, polar )},
      {1.0, XCKernel( backend, Kernel::R2SCAN_C, polar )}
    };
  }
  else if ( func == Functional::R2SCANh ){
    hyb_coefs = {0.1, 0.0, 0.0};
    kerns = { 
      {0.9, XCKernel( backend, Kernel::R2SCAN_X, polar )},
      {1.0, XCKernel( backend, Kernel::R2SCAN_C, polar )}
    };
  }
  else if ( func == Functional::R2SCAN50 ){
    hyb_coefs = {0.5, 0.0, 0.0};
    kerns = { 
      {0.5, XCKernel( backend, Kernel::R2SCAN_X, polar )},
      {1.0, XCKernel( backend, Kernel::R2SCAN_C, polar )}
    };
  }
  else if ( func == Functional::M05 ){
    hyb_coefs = {0.28, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M05_X, polar )},
      {1.0, XCKernel( backend, Kernel::M05_C, polar )}
    };
  }
  else if ( func == Functional::M06 ){
    hyb_coefs = {0.27, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M06_X, polar )},
      {1.0, XCKernel( backend, Kernel::M06_C, polar )}
    };
  }
  else if ( func == Functional::M08HX ){
    hyb_coefs = {0.5223, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M08_HX_X, polar )},
      {1.0, XCKernel( backend, Kernel::M08_HX_C, polar )}
    };
  }
  else if ( func == Functional::M08SO ){
    hyb_coefs = {0.5679, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M08_SO_X, polar )},
      {1.0, XCKernel( backend, Kernel::M08_SO_C, polar )}
    };
  }
  else if ( func == Functional::M052X ){
    hyb_coefs = {0.56, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M05_2X_X, polar )},
      {1.0, XCKernel( backend, Kernel::M05_2X_C, polar )}
    };
  }
  else if ( func == Functional::M06SX ){
    hyb_coefs = {0.0, 0.335, 0.1};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M06_SX_X, polar )},
      {1.0, XCKernel( backend, Kernel::M06_SX_C, polar )}
    };
  }
  else if ( func == Functional::CF22D ){
    hyb_coefs = {0.462806, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::CF22D_X, polar )},
      {1.0, XCKernel( backend, Kernel::CF22D_C, polar )}
    };
  }
  else if ( func == Functional::SOGGA11X ){
    hyb_coefs = {0.4015, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::SOGGA11_X_X, polar )},
      {1.0, XCKernel( backend, Kernel::SOGGA11_X_C, polar )}
    };
  }
  else if ( func == Functional::M06HF ){
    hyb_coefs = {1.0, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M06_HF_X, polar )},
      {1.0, XCKernel( backend, Kernel::M06_HF_C, polar )}
    };
  }
  else if ( func == Functional::M11 ){
    hyb_coefs = {1.0, -0.572, 0.25};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::M11_X, polar )},
      {1.0, XCKernel( backend, Kernel::M11_C, polar )}
    };
  }
  else if ( func == Functional::MN12L )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::MN12_L_X, polar )},
      {1.0, XCKernel( backend, Kernel::MN12_L_C, polar )}
    };
  else if ( func == Functional::MN12SX ){
    hyb_coefs = {0.0, 0.25, 0.11};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::MN12_SX_X, polar )},
      {1.0, XCKernel( backend, Kernel::MN12_SX_C, polar )}
    };
  }
  else if ( func == Functional::MN15 ){
    hyb_coefs = {0.44, 0.0, 0.0};
    kerns = { 
      {1.0, XCKernel( backend, Kernel::MN15_X, polar )},
      {1.0, XCKernel( backend, Kernel::MN15_C, polar )}
    };
  }
  else if ( func == Functional::MN15L )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::MN15_L_X, polar )},
      {1.0, XCKernel( backend, Kernel::MN15_L_C, polar )}
    };
  else if ( func == Functional::revM06L )
    kerns = { 
      {1.0, XCKernel( backend, Kernel::revM06_L_X, polar )},
      {1.0, XCKernel( backend, Kernel::revM06_L_C, polar )}
    };  
  else {
    EXCHCXX_BOOL_CHECK( "Functional NYI Through Builtin Backend", false );
  }


  return {kerns, hyb_coefs};
}

XCFunctional::XCFunctional() = default;

XCFunctional::XCFunctional( const std::vector< XCKernel > &ks ) {
  hyb_coefs_ = {0.0, 0.0, 0.0};
  for(const auto& k : ks )
    kernels_.push_back( { 1., k } );
}
XCFunctional::XCFunctional( const std::vector< XCKernel > &ks, const HybCoeffs& hyb ) :
  hyb_coefs_(hyb) {
  for(const auto& k : ks )
    kernels_.push_back( { 1., k } );
}

XCFunctional::XCFunctional( const std::initializer_list< value_type >& list ) : kernels_{ list } { }
XCFunctional::XCFunctional( const std::initializer_list< value_type >& list , const HybCoeffs& hyb ) :
  kernels_{ list }, hyb_coefs_(hyb) { }

XCFunctional::XCFunctional( const std::vector<value_type>& ks ) :
  kernels_(ks) { }
  XCFunctional::XCFunctional( const std::vector<value_type>& ks, const HybCoeffs& hyb ) :
  kernels_(ks), hyb_coefs_(hyb) { }
XCFunctional::XCFunctional( std::vector<value_type>&& ks ) :
  kernels_(std::move(ks)) { }
  XCFunctional::XCFunctional( std::vector<value_type>&& ks, HybCoeffs&& hyb ) :
  kernels_(std::move(ks)), hyb_coefs_(std::move(hyb)) { }


  XCFunctional::XCFunctional( 
    const Backend        backend, 
    const Functional func,
    const Spin           polar
  ) {
    auto [kernels, hyb_coefs] = functional_factory(backend, func, polar);
    hyb_coefs_ = std::move(hyb_coefs);
    for(const auto& k : kernels )
      kernels_.push_back( k );
  }



XCFunctional& XCFunctional::operator=( const XCFunctional& ) = default;
XCFunctional& XCFunctional::operator=( XCFunctional&&      ) noexcept = default;

XCFunctional::XCFunctional( const XCFunctional& )       = default;
XCFunctional::XCFunctional( XCFunctional&& )  noexcept  = default;





void _scal( size_t len, double alpha, double* v ) {
  for( auto k = 0ul; k < len; ++k ) v[k] *= alpha;
}
void _addscal( size_t len, double alpha, double* v, const double* w ) {
  for( auto k = 0ul; k < len; ++k ) v[k] += alpha * w[k];
}



LDA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_inc( 
        kernels_[i].first, N, rho, eps
      );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;
      kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );
  
    }
  }

}


LDA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_exc_buffer = exc_buffer_len(N);
  const size_t len_vxc_buffer = vrho_buffer_len(N);

  std::vector<double> eps_scr, vxc_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vxc_scr.resize( len_vxc_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vxc, len_vxc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_vxc_inc(
        kernels_[i].first, N, rho, eps, vxc
      );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;
      double* vxc_eval = i ? vxc_scr.data() : vxc;
      kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vxc_eval);

      if( i ) {

        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
        _addscal( len_vxc_buffer, kernels_[i].first, vxc, vxc_eval );

      } else {

        _scal( len_exc_buffer, kernels_[i].first, eps );
        _scal( len_vxc_buffer, kernels_[i].first, vxc );

      }

    }

  }

}
LDA_FXC_GENERATOR( XCFunctional::eval_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_fxc_buffer = v2rho2_buffer_len(N);

  std::vector<double> fxc_scr;
  if( kernels_.size() > 1 && !supports_inc_interface() )
    fxc_scr.resize( len_fxc_buffer );

  std::fill_n( fxc, len_fxc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {
      kernels_[i].second.eval_fxc_inc(
        kernels_[i].first, N, rho, fxc
      );
    } else {
      double* fxc_eval = i ? fxc_scr.data() : fxc;
      kernels_[i].second.eval_fxc(N, rho, fxc_eval);

      if( i )
        _addscal( len_fxc_buffer, kernels_[i].first, fxc, fxc_eval );
      else
        _scal( len_fxc_buffer, kernels_[i].first, fxc );
    }
  }
}

LDA_VXC_FXC_GENERATOR( XCFunctional::eval_vxc_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_vxc_buffer = vrho_buffer_len(N);
  const size_t len_fxc_buffer = v2rho2_buffer_len(N);

  std::vector<double> vxc_scr, fxc_scr;
  if( kernels_.size() > 1 && !supports_inc_interface() ) {
    vxc_scr.resize( len_vxc_buffer );
    fxc_scr.resize( len_fxc_buffer );
  }

  std::fill_n( vxc, len_vxc_buffer, 0. );
  std::fill_n( fxc, len_fxc_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {
      kernels_[i].second.eval_vxc_fxc_inc(
        kernels_[i].first, N, rho, vxc, fxc
      );
    } else {
      double* vxc_eval = i ? vxc_scr.data() : vxc;
      double* fxc_eval = i ? fxc_scr.data() : fxc;
      kernels_[i].second.eval_vxc_fxc(N, rho, vxc_eval, fxc_eval);

      if( i ) {
        _addscal( len_vxc_buffer, kernels_[i].first, vxc, vxc_eval );
        _addscal( len_fxc_buffer, kernels_[i].first, fxc, fxc_eval );
      } else {
        _scal( len_vxc_buffer, kernels_[i].first, vxc );
        _scal( len_fxc_buffer, kernels_[i].first, fxc );
      }
    }
  }
}


// GGA Interfaces

GGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, eps
        );
      else
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, eps
        );

    } else {

      double* eps_eval = i ? eps_scr.data() : eps;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
      else
        kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );

    }
  
  }

}


GGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, eps, vrho
        );

    } else {

      double* eps_eval    = i ? eps_scr.data()    : eps;
      double* vrho_eval   = i ? vrho_scr.data()   : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, 
          vsigma_eval );
      else
        kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

      if( i ) {

        _addscal( len_exc_buffer,    kernels_[i].first, eps,    eps_eval  );
        _addscal( len_vrho_buffer,   kernels_[i].first, vrho,   vrho_eval );
   
        if( kernels_[i].second.is_gga() )
          _addscal( len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval );

      } else {

        _scal( len_exc_buffer,    kernels_[i].first, eps  );
        _scal( len_vrho_buffer,   kernels_[i].first, vrho );

        if( kernels_[i].second.is_gga() )
          _scal( len_vsigma_buffer, kernels_[i].first, vsigma );

      }

    }
  
  }

}
GGA_FXC_GENERATOR( XCFunctional::eval_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_v2rho2_buffer = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2sigma2_buffer = v2sigma2_buffer_len(N);

  std::vector<double> v2sigma2_scr, v2rhosigma_scr, v2rho2_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    v2rho2_scr.resize( len_v2rho2_buffer );
    v2rhosigma_scr.resize( len_v2rhosigma_buffer );
    v2sigma2_scr.resize( len_v2sigma2_buffer );
  }

  std::fill_n( v2rho2, len_v2rho2_buffer, 0. );
  std::fill_n( v2rhosigma, len_v2rhosigma_buffer, 0. );
  std::fill_n( v2sigma2, len_v2sigma2_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_inc(
          kernels_[i].first, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2
        );
      else
        kernels_[i].second.eval_fxc_inc(
          kernels_[i].first, N, rho, v2rho2
        );

    } else {

      double* v2rho2_eval    = i ? v2rho2_scr.data()    : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr.data() : v2rhosigma;
      double* v2sigma2_eval  = i ? v2sigma2_scr.data()  : v2sigma2;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc(N, rho, sigma, v2rho2_eval, 
          v2rhosigma_eval, v2sigma2_eval );
      else
        kernels_[i].second.eval_fxc(N, rho, v2rho2_eval);

      if( i ) {
        _addscal( len_v2rho2_buffer,    kernels_[i].first, v2rho2,    v2rho2_eval  );
        if( kernels_[i].second.is_gga() ){
          _addscal( len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma, v2rhosigma_eval );
          _addscal( len_v2sigma2_buffer,   kernels_[i].first, v2sigma2,   v2sigma2_eval );
        }

      } else {
        _scal( len_v2rho2_buffer,    kernels_[i].first, v2rho2 );
        if( kernels_[i].second.is_gga() ){
          _scal( len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma );
          _scal( len_v2sigma2_buffer,   kernels_[i].first, v2sigma2 );
        }
      }
    }
  }
}

GGA_VXC_FXC_GENERATOR( XCFunctional::eval_vxc_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_v2rho2_buffer = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);

  std::vector<double> vrho_scr, vsigma_scr, v2rho2_scr, v2rhosigma_scr, v2sigma2_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
    v2rho2_scr.resize( len_v2rho2_buffer );
    v2rhosigma_scr.resize( len_v2rhosigma_buffer );
    v2sigma2_scr.resize( len_v2sigma2_buffer );
  }

  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );
  std::fill_n( v2rho2, len_v2rho2_buffer, 0. );
  std::fill_n( v2rhosigma, len_v2rhosigma_buffer, 0. );
  std::fill_n( v2sigma2, len_v2sigma2_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_vxc_fxc_inc(
          kernels_[i].first, N, rho, sigma, vrho, vsigma,
          v2rho2, v2rhosigma, v2sigma2
        );
      else
        kernels_[i].second.eval_vxc_fxc_inc(
          kernels_[i].first, N, rho, vrho, v2rho2
        );

    } else {

      double* vrho_eval = i ? vrho_scr.data() : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;
      double* v2rho2_eval = i ? v2rho2_scr.data() : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr.data() : v2rhosigma;
      double* v2sigma2_eval = i ? v2sigma2_scr.data() : v2sigma2;

      if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc(N, rho, sigma, vrho_eval, vsigma_eval,
                                        v2rho2_eval, v2rhosigma_eval, v2sigma2_eval);
      } else {
        kernels_[i].second.eval_vxc_fxc(N, rho, vrho_eval, v2rho2_eval);
      }

      if (i) {
        _addscal(len_vrho_buffer, kernels_[i].first, vrho, vrho_eval);
        _addscal(len_v2rho2_buffer, kernels_[i].first, v2rho2, v2rho2_eval);

        if (kernels_[i].second.is_gga()) {
          _addscal(len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval);
          _addscal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma, v2rhosigma_eval);
          _addscal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2, v2sigma2_eval);
        }
      } else {
        _scal(len_vrho_buffer, kernels_[i].first, vrho);
        _scal(len_v2rho2_buffer, kernels_[i].first, v2rho2);

        if (kernels_[i].second.is_gga()) {
          _scal(len_vsigma_buffer, kernels_[i].first, vsigma);
          _scal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma);
          _scal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2);
        }
      }
    }
  }
}


// mGGA Interfaces

MGGA_EXC_GENERATOR( XCFunctional::eval_exc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_exc_buffer = exc_buffer_len(N);

  std::vector<double> eps_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr.resize( len_exc_buffer );

  std::fill_n( eps, len_exc_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, sigma, eps
        );
      else
        kernels_[i].second.eval_exc_inc(
          kernels_[i].first, N, rho, eps
        );

    } else { 

      double* eps_eval = i ? eps_scr.data() : eps;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc(N, rho, sigma, lapl, tau, eps_eval);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc(N, rho, sigma, eps_eval);
      else
        kernels_[i].second.eval_exc(N, rho, eps_eval);

      if( i ) 
        _addscal( len_exc_buffer, kernels_[i].first, eps, eps_eval );
      else
        _scal( len_exc_buffer, kernels_[i].first, eps );

    }
  
  }

}


MGGA_EXC_VXC_GENERATOR( XCFunctional::eval_exc_vxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_exc_buffer    = exc_buffer_len(N);
  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_vlapl_buffer  = vlapl_buffer_len(N);
  const size_t len_vtau_buffer   = vtau_buffer_len(N);

  std::vector<double> eps_scr, vrho_scr, vsigma_scr, vlapl_scr, vtau_scr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr.resize( len_exc_buffer );
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
    vlapl_scr.resize( len_vlapl_buffer );
    vtau_scr.resize( len_vtau_buffer );
  }

  std::fill_n( eps, len_exc_buffer, 0. );
  std::fill_n( vrho, len_vrho_buffer, 0. );
  std::fill_n( vsigma, len_vsigma_buffer, 0. );
  std::fill_n( vlapl, len_vlapl_buffer, 0. );
  std::fill_n( vtau, len_vtau_buffer, 0. );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps, 
          vrho, vsigma, vlapl, vtau 
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc(
          kernels_[i].first, N, rho, eps, vrho
        );
    
    } else {

      double* eps_eval    = i ? eps_scr.data()    : eps;
      double* vrho_eval   = i ? vrho_scr.data()   : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;
      double* vlapl_eval  = i ? vlapl_scr.data()  : vlapl;
      double* vtau_eval   = i ? vtau_scr.data()   : vtau;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, lapl, tau, eps_eval, 
          vrho_eval, vsigma_eval, vlapl_eval, vtau_eval );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc(N, rho, sigma, eps_eval, vrho_eval, 
          vsigma_eval );
      else
        kernels_[i].second.eval_exc_vxc(N, rho, eps_eval, vrho_eval);

      if( i ) {

        _addscal( len_exc_buffer,    kernels_[i].first, eps,    eps_eval  );
        _addscal( len_vrho_buffer,   kernels_[i].first, vrho,   vrho_eval );
   
        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          _addscal( len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval );

        if( kernels_[i].second.needs_laplacian() ) 
          _addscal( len_vlapl_buffer, kernels_[i].first, vlapl, vlapl_eval );

        if( kernels_[i].second.is_mgga() ) 
          _addscal( len_vtau_buffer,  kernels_[i].first, vtau,  vtau_eval  );

      } else {

        _scal( len_exc_buffer,    kernels_[i].first, eps  );
        _scal( len_vrho_buffer,   kernels_[i].first, vrho );

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          _scal( len_vsigma_buffer, kernels_[i].first, vsigma );

        if( kernels_[i].second.needs_laplacian() ) 
          _scal( len_vlapl_buffer, kernels_[i].first, vlapl );

        if( kernels_[i].second.is_mgga() ) 
          _scal( len_vtau_buffer,  kernels_[i].first, vtau  );

      }

    }
  
  }

}

MGGA_FXC_GENERATOR( XCFunctional::eval_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_v2rho2_buffer    = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2rholapl_buffer  = v2rholapl_buffer_len(N);
  const size_t len_v2rhotau_buffer   = v2rhotau_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);
  const size_t len_v2sigmalapl_buffer = v2sigmalapl_buffer_len(N);
  const size_t len_v2sigmatau_buffer  = v2sigmatau_buffer_len(N);
  const size_t len_v2lapl2_buffer    = v2lapl2_buffer_len(N);
  const size_t len_v2lapltau_buffer  = v2lapltau_buffer_len(N);
  const size_t len_v2tau2_buffer     = v2tau2_buffer_len(N);

  std::vector<double> v2rho2_scr, v2rhosigma_scr, v2rholapl_scr, v2rhotau_scr,
    v2sigma2_scr, v2sigmalapl_scr, v2sigmatau_scr, v2lapl2_scr, v2lapltau_scr, v2tau2_scr;
  if( kernels_.size() > 1 && !supports_inc_interface() ) {
    v2rho2_scr.resize( len_v2rho2_buffer );
    v2rhosigma_scr.resize( len_v2rhosigma_buffer );
    v2rholapl_scr.resize( len_v2rholapl_buffer );
    v2rhotau_scr.resize( len_v2rhotau_buffer );
    v2sigma2_scr.resize( len_v2sigma2_buffer );
    v2sigmalapl_scr.resize( len_v2sigmalapl_buffer );
    v2sigmatau_scr.resize( len_v2sigmatau_buffer );
    v2lapl2_scr.resize( len_v2lapl2_buffer );
    v2lapltau_scr.resize( len_v2lapltau_buffer );
    v2tau2_scr.resize( len_v2tau2_buffer );
  }

  std::fill_n( v2rho2, len_v2rho2_buffer, 0. );
  std::fill_n( v2rhosigma, len_v2rhosigma_buffer, 0. );
  std::fill_n( v2rholapl, len_v2rholapl_buffer, 0. );
  std::fill_n( v2rhotau, len_v2rhotau_buffer, 0. );
  std::fill_n( v2sigma2, len_v2sigma2_buffer, 0. );
  std::fill_n( v2sigmalapl, len_v2sigmalapl_buffer, 0. );
  std::fill_n( v2sigmatau, len_v2sigmatau_buffer, 0. );
  std::fill_n( v2lapl2, len_v2lapl2_buffer, 0. );
  std::fill_n( v2lapltau, len_v2lapltau_buffer, 0. );
  std::fill_n( v2tau2, len_v2tau2_buffer, 0. );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_fxc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau,
          v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_inc(
          kernels_[i].first, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2
        );
      else
        kernels_[i].second.eval_fxc_inc(
          kernels_[i].first, N, rho, v2rho2
        );

    } else {

      double* v2rho2_eval    = i ? v2rho2_scr.data()    : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr.data() : v2rhosigma;
      double* v2rholapl_eval  = i ? v2rholapl_scr.data()  : v2rholapl;
      double* v2rhotau_eval   = i ? v2rhotau_scr.data()   : v2rhotau;
      double* v2sigma2_eval   = i ? v2sigma2_scr.data()   : v2sigma2;
      double* v2sigmalapl_eval = i ? v2sigmalapl_scr.data() : v2sigmalapl;
      double* v2sigmatau_eval  = i ? v2sigmatau_scr.data()  : v2sigmatau;
      double* v2lapl2_eval    = i ? v2lapl2_scr.data()    : v2lapl2;
      double* v2lapltau_eval  = i ? v2lapltau_scr.data()  : v2lapltau;
      double* v2tau2_eval     = i ? v2tau2_scr.data()     : v2tau2;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_fxc(N, rho, sigma, lapl, tau, v2rho2_eval, 
          v2rhosigma_eval, v2rholapl_eval, v2rhotau_eval, v2sigma2_eval, v2sigmalapl_eval, 
          v2sigmatau_eval, v2lapl2_eval, v2lapltau_eval, v2tau2_eval);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc(N, rho, sigma, v2rho2_eval, v2rhosigma_eval, v2sigma2_eval );
      else
        kernels_[i].second.eval_fxc(N, rho, v2rho2_eval);

      if (i) {
        _addscal(len_v2rho2_buffer, kernels_[i].first, v2rho2, v2rho2_eval);

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() ){
          _addscal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma, v2rhosigma_eval);
          _addscal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2, v2sigma2_eval);
        }

        if( kernels_[i].second.needs_laplacian() ) {
          _addscal(len_v2rholapl_buffer, kernels_[i].first, v2rholapl, v2rholapl_eval);
          _addscal(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl, v2sigmalapl_eval);
          _addscal(len_v2lapl2_buffer, kernels_[i].first, v2lapl2, v2lapl2_eval);
        }
        if( kernels_[i].second.is_mgga() ) {
          _addscal(len_v2rhotau_buffer, kernels_[i].first, v2rhotau, v2rhotau_eval);
          _addscal(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau, v2sigmatau_eval);
          _addscal(len_v2tau2_buffer, kernels_[i].first, v2tau2, v2tau2_eval);
        }
        if ( kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga() ) {
          _addscal(len_v2lapltau_buffer, kernels_[i].first, v2lapltau, v2lapltau_eval);
        }

      } else{

        _scal(len_v2rho2_buffer, kernels_[i].first, v2rho2);

        if (kernels_[i].second.is_gga() or kernels_[i].second.is_mgga()) {
          _scal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma);
          _scal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2);
        }

        if (kernels_[i].second.needs_laplacian()) {
          _scal(len_v2rholapl_buffer, kernels_[i].first, v2rholapl);
          _scal(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl);
          _scal(len_v2lapl2_buffer, kernels_[i].first, v2lapl2);
        }

        if (kernels_[i].second.is_mgga()) {
          _scal(len_v2rhotau_buffer, kernels_[i].first, v2rhotau);
          _scal(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau);
          _scal(len_v2tau2_buffer, kernels_[i].first, v2tau2);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          _scal(len_v2lapltau_buffer, kernels_[i].first, v2lapltau);
        }
      }
    }
  }
}


MGGA_VXC_FXC_GENERATOR( XCFunctional::eval_vxc_fxc ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_vlapl_buffer  = vlapl_buffer_len(N);
  const size_t len_vtau_buffer   = vtau_buffer_len(N);
  const size_t len_v2rho2_buffer    = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2rholapl_buffer  = v2rholapl_buffer_len(N);
  const size_t len_v2rhotau_buffer   = v2rhotau_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);
  const size_t len_v2sigmalapl_buffer = v2sigmalapl_buffer_len(N);
  const size_t len_v2sigmatau_buffer  = v2sigmatau_buffer_len(N);
  const size_t len_v2lapl2_buffer    = v2lapl2_buffer_len(N);
  const size_t len_v2lapltau_buffer  = v2lapltau_buffer_len(N);
  const size_t len_v2tau2_buffer     = v2tau2_buffer_len(N);

  std::vector<double> vrho_scr, vsigma_scr, vlapl_scr, vtau_scr;
  std::vector<double> v2rho2_scr, v2rhosigma_scr, v2rholapl_scr, v2rhotau_scr,
    v2sigma2_scr, v2sigmalapl_scr, v2sigmatau_scr, v2lapl2_scr, v2lapltau_scr, v2tau2_scr;
  
  if( kernels_.size() > 1 && !supports_inc_interface() ) {
    vrho_scr.resize( len_vrho_buffer );
    vsigma_scr.resize( len_vsigma_buffer );
    vlapl_scr.resize( len_vlapl_buffer );
    vtau_scr.resize( len_vtau_buffer );
    v2rho2_scr.resize( len_v2rho2_buffer );
    v2rhosigma_scr.resize( len_v2rhosigma_buffer );
    v2rholapl_scr.resize( len_v2rholapl_buffer );
    v2rhotau_scr.resize(len_v2rhotau_buffer);
    v2sigma2_scr.resize(len_v2sigma2_buffer);
    v2sigmalapl_scr.resize(len_v2sigmalapl_buffer);
    v2sigmatau_scr.resize(len_v2sigmatau_buffer);
    v2lapl2_scr.resize(len_v2lapl2_buffer);
    v2lapltau_scr.resize(len_v2lapltau_buffer);
    v2tau2_scr.resize(len_v2tau2_buffer);
  }

  std::fill_n(vrho, len_vrho_buffer, 0.);
  std::fill_n(vsigma, len_vsigma_buffer, 0.);
  std::fill_n(vlapl, len_vlapl_buffer, 0.);
  std::fill_n(vtau, len_vtau_buffer, 0.);
  std::fill_n(v2rho2, len_v2rho2_buffer, 0.);
  std::fill_n(v2rhosigma, len_v2rhosigma_buffer, 0.);
  std::fill_n(v2rholapl, len_v2rholapl_buffer, 0.);
  std::fill_n(v2rhotau, len_v2rhotau_buffer, 0.);
  std::fill_n(v2sigma2, len_v2sigma2_buffer, 0.);
  std::fill_n(v2sigmalapl, len_v2sigmalapl_buffer, 0.);
  std::fill_n(v2sigmatau, len_v2sigmatau_buffer, 0.);
  std::fill_n(v2lapl2, len_v2lapl2_buffer, 0.);
  std::fill_n(v2lapltau, len_v2lapltau_buffer, 0.);
  std::fill_n(v2tau2, len_v2tau2_buffer, 0.);

  for (auto i = 0ul; i < kernels_.size(); ++i) {

    if( supports_inc_interface() ) {

      if (kernels_[i].second.is_mgga()) {
        kernels_[i].second.eval_vxc_fxc_inc(
          kernels_[i].first, N, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau,
          v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2
        );
      } else if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc_inc(
          kernels_[i].first, N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2
        );
      } else {
        kernels_[i].second.eval_vxc_fxc_inc(
          kernels_[i].first, N, rho, vrho, v2rho2
        );
      }

    } else {

      double* vrho_eval = i ? vrho_scr.data() : vrho;
      double* vsigma_eval = i ? vsigma_scr.data() : vsigma;
      double* vlapl_eval = i ? vlapl_scr.data() : vlapl;
      double* vtau_eval = i ? vtau_scr.data() : vtau;
      double* v2rho2_eval = i ? v2rho2_scr.data() : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr.data() : v2rhosigma;
      double* v2rholapl_eval = i ? v2rholapl_scr.data() : v2rholapl;
      double* v2rhotau_eval = i ? v2rhotau_scr.data() : v2rhotau;
      double* v2sigma2_eval = i ? v2sigma2_scr.data() : v2sigma2;
      double* v2sigmalapl_eval = i ? v2sigmalapl_scr.data() : v2sigmalapl;
      double* v2sigmatau_eval = i ? v2sigmatau_scr.data() : v2sigmatau;
      double* v2lapl2_eval = i ? v2lapl2_scr.data() : v2lapl2;
      double* v2lapltau_eval = i ? v2lapltau_scr.data() : v2lapltau;
      double* v2tau2_eval = i ? v2tau2_scr.data() : v2tau2;

      if (kernels_[i].second.is_mgga()) {
        kernels_[i].second.eval_vxc_fxc(
          N, rho, sigma, lapl, tau, vrho_eval, vsigma_eval, vlapl_eval, vtau_eval,
          v2rho2_eval, v2rhosigma_eval, v2rholapl_eval, v2rhotau_eval,
          v2sigma2_eval, v2sigmalapl_eval, v2sigmatau_eval, v2lapl2_eval,
          v2lapltau_eval, v2tau2_eval);
      } else if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc(
          N, rho, sigma, vrho_eval, vsigma_eval, v2rho2_eval, v2rhosigma_eval,
          v2sigma2_eval);
      } else {
        kernels_[i].second.eval_vxc_fxc(N, rho, vrho_eval, v2rho2_eval);
      }

      if (i) {
        _addscal(len_vrho_buffer, kernels_[i].first, vrho, vrho_eval);
        _addscal(len_v2rho2_buffer, kernels_[i].first, v2rho2, v2rho2_eval);

        if (kernels_[i].second.is_gga() || kernels_[i].second.is_mgga()) {
          _addscal(len_vsigma_buffer, kernels_[i].first, vsigma, vsigma_eval);
          _addscal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma, v2rhosigma_eval);
          _addscal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2, v2sigma2_eval);
        }

        if (kernels_[i].second.needs_laplacian()) {
          _addscal(len_vlapl_buffer, kernels_[i].first, vlapl, vlapl_eval);
          _addscal(len_v2rholapl_buffer, kernels_[i].first, v2rholapl, v2rholapl_eval);
          _addscal(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl, v2sigmalapl_eval);
          _addscal(len_v2lapl2_buffer, kernels_[i].first, v2lapl2, v2lapl2_eval);
        }

        if (kernels_[i].second.is_mgga()) {
          _addscal(len_vtau_buffer, kernels_[i].first, vtau, vtau_eval);
          _addscal(len_v2rhotau_buffer, kernels_[i].first, v2rhotau, v2rhotau_eval);
          _addscal(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau, v2sigmatau_eval);
          _addscal(len_v2tau2_buffer, kernels_[i].first, v2tau2, v2tau2_eval);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          _addscal(len_v2lapltau_buffer, kernels_[i].first, v2lapltau, v2lapltau_eval);
        }
      } else {
        _scal(len_vrho_buffer, kernels_[i].first, vrho);
        _scal(len_v2rho2_buffer, kernels_[i].first, v2rho2);

        if (kernels_[i].second.is_gga() || kernels_[i].second.is_mgga()) {
          _scal(len_vsigma_buffer, kernels_[i].first, vsigma);
          _scal(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma);
          _scal(len_v2sigma2_buffer, kernels_[i].first, v2sigma2);
        }

        if (kernels_[i].second.needs_laplacian()) {
          _scal(len_vlapl_buffer, kernels_[i].first, vlapl);
          _scal(len_v2rholapl_buffer, kernels_[i].first, v2rholapl);
          _scal(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl);
          _scal(len_v2lapl2_buffer, kernels_[i].first, v2lapl2);
        }

        if (kernels_[i].second.is_mgga()) {
          _scal(len_vtau_buffer, kernels_[i].first, vtau);
          _scal(len_v2rhotau_buffer, kernels_[i].first, v2rhotau);
          _scal(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau);
          _scal(len_v2tau2_buffer, kernels_[i].first, v2tau2);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          _scal(len_v2lapltau_buffer, kernels_[i].first, v2lapltau);
        }
      }
    }
  }
}


} // ExchCXX


