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


#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/impl/builtin/fwd.hpp>


namespace ExchCXX {
namespace detail {

class BuiltinKernelInterface : public XCKernelImpl {

  using unique_me = XCKernelImpl::unique_me;

  unique_me clone_() const override;

  Kernel whatami_;

  std::unique_ptr<BuiltinKernel> impl_;

  bool is_lda_()       const noexcept override;
  bool is_gga_()       const noexcept override;
  bool is_mgga_()      const noexcept override;
  bool is_hyb_()       const noexcept override;
  bool needs_laplacian_()       const noexcept override;
  bool is_polarized_() const noexcept override;
  double hyb_exx_()    const noexcept override;

  bool supports_inc_interface_() const noexcept override;

  // LDA interface
  LDA_EXC_GENERATOR( eval_exc_ )                 const override;
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const override;
  LDA_EXC_INC_GENERATOR( eval_exc_inc_ )         const override;
  LDA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const override;

  // GGA interface
  GGA_EXC_GENERATOR( eval_exc_ )                 const override;
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const override;
  GGA_EXC_INC_GENERATOR( eval_exc_inc_ )         const override;
  GGA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const override;

  // MGGA interface
  MGGA_EXC_GENERATOR( eval_exc_ )                 const override;
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const override;
  MGGA_EXC_INC_GENERATOR( eval_exc_inc_ )         const override;
  MGGA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const override;

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  LDA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const override;
  LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const override;
  LDA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const override;
  LDA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const override;

  // GGA interface
  GGA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const override;
  GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const override;
  GGA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const override;
  GGA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const override;

  // MGGA interface
  MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const override;
  MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const override;
  MGGA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const override;
  MGGA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const override;

#endif

public:

  BuiltinKernelInterface() = delete;
  
  BuiltinKernelInterface( Kernel kern, Spin p);
  BuiltinKernelInterface( const BuiltinKernelInterface& );
  BuiltinKernelInterface( BuiltinKernelInterface&& ) noexcept = delete;

  // Destroy interal Builtin data
  ~BuiltinKernelInterface() noexcept;
    
};

}
}
