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

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX {

BidirectionalMap<std::string, Kernel> kernel_map{
    {{"SlaterExchange", Kernel::SlaterExchange},
     {"PBE_X", Kernel::PBE_X},
     {"PBE_C", Kernel::PBE_C},
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
     {"revPBE_X", Kernel::revPBE_X},
     {"LYP", Kernel::LYP},
     {"B3LYP", Kernel::B3LYP},
     {"PBE0", Kernel::PBE0},
     {"VWN3", Kernel::VWN3},
     {"VWN5", Kernel::VWN5},
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

std::ostream& operator<<( std::ostream& out, Kernel kern ) {
  out << kernel_map.key(kern);
  return out;
}

XCKernel::XCKernel( 
  const Backend backend, 
  const Kernel kern, 
  const Spin polar) : 
XCKernel( kernel_factory( backend, kern, polar ) ) { }

XCKernel::XCKernel(
  const libxc_name_string& xc_name,
  const Spin polar) :
XCKernel( libxc_kernel_factory( xc_name.get(), polar ) ) { }


XCKernel::XCKernel( impl_ptr&& ptr ) :
  pimpl_(std::move(ptr)) { }

XCKernel::XCKernel( const XCKernel& other ) :
  pimpl_(other.pimpl_->clone()) { }

XCKernel::XCKernel( XCKernel&& other ) noexcept = default;

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept =default;


XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = XCKernel(other);
}

XCKernel::~XCKernel() noexcept = default;

}
