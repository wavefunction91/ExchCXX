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

#include <exchcxx/enums/enums.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/util/exchcxx_macros.hpp>

// Standard Libs
#include <vector>
#include <algorithm>
#include <memory>


namespace ExchCXX {

/**
 *  \brief A class which manages the lifetime and evaluation
 *  of exchange-correlation (XC) kernels for density functional
 *  theory.
 */
class XCKernel {

  using impl_ptr = std::unique_ptr< detail::XCKernelImpl >;

  impl_ptr pimpl_; ///< Pointer to implementation

public:

  // Avoid stateless kernel
  XCKernel() = delete;
  
  XCKernel( const Backend backend, const Kernel kern, 
    const Spin polar );
  XCKernel( const libxc_name_string& xc_name, 
    const Spin polar );
  XCKernel( const Kernel kern, const Spin polar ) : 
    XCKernel( Backend::libxc, kern, polar ){ };

  XCKernel( impl_ptr&& ptr )                     ;
  XCKernel( const XCKernel& )                    ;
  XCKernel( XCKernel&&      )            noexcept;
  XCKernel& operator=( const XCKernel& )         ;
  XCKernel& operator=( XCKernel&&      ) noexcept;

  // Destroy interal Libxc data
  ~XCKernel() noexcept;



  bool is_lda()       const noexcept { return pimpl_->is_lda();       }
  bool is_gga()       const noexcept { return pimpl_->is_gga();       }
  bool is_mgga()      const noexcept { return pimpl_->is_mgga();      }
  bool is_hyb()       const noexcept { return pimpl_->is_hyb();       }
  bool is_polarized() const noexcept { return pimpl_->is_polarized(); }

  bool needs_laplacian() const noexcept { return pimpl_->needs_laplacian(); }
  bool needs_tau()    const noexcept { return pimpl_->needs_tau(); }
  
  double hyb_exx() const noexcept { return pimpl_->hyb_exx(); }

  bool supports_inc_interface() const noexcept {
    return pimpl_->supports_inc_interface();
  }

  inline size_t rho_buffer_len( size_t npts ) const noexcept {
    return pimpl_->rho_buffer_len( npts );
  }
  inline size_t sigma_buffer_len( size_t npts ) const noexcept {
    return pimpl_->sigma_buffer_len( npts );
  }
  inline size_t lapl_buffer_len( size_t npts ) const noexcept {
    return pimpl_->lapl_buffer_len( npts );
  }
  inline size_t tau_buffer_len( size_t npts ) const noexcept {
    return pimpl_->tau_buffer_len( npts );
  }

  inline size_t exc_buffer_len( size_t npts ) const noexcept {
    return pimpl_->exc_buffer_len( npts );
  }
  inline size_t vrho_buffer_len( size_t npts ) const noexcept {
    return pimpl_->vrho_buffer_len( npts );
  }
  inline size_t vsigma_buffer_len( size_t npts ) const noexcept {
    return pimpl_->vsigma_buffer_len( npts );
  }
  inline size_t vlapl_buffer_len( size_t npts ) const noexcept {
    return pimpl_->vlapl_buffer_len( npts );
  }
  inline size_t vtau_buffer_len( size_t npts ) const noexcept {
    return pimpl_->vtau_buffer_len( npts );
  }


  template <typename... Args>
  void eval_exc( Args&&... args ) const {
    pimpl_->eval_exc( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc( Args&&... args ) const {
    pimpl_->eval_exc_vxc( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_inc( Args&&... args ) const {
    pimpl_->eval_exc_inc( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_inc( Args&&... args ) const {
    pimpl_->eval_exc_vxc_inc( std::forward<Args>(args)... );
  }

  // Device code
#ifdef EXCHCXX_ENABLE_DEVICE
  
  template <typename... Args>
  void eval_exc_device( Args&&... args ) const {
    pimpl_->eval_exc_device( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_device( Args&&... args ) const {
    pimpl_->eval_exc_vxc_device( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_inc_device( Args&&... args ) const {
    pimpl_->eval_exc_inc_device( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_inc_device( Args&&... args ) const {
    pimpl_->eval_exc_vxc_inc_device( std::forward<Args>(args)... );
  }

#endif



};









}

