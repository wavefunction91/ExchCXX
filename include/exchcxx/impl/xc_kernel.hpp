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


#include <exchcxx/util/exchcxx_macros.hpp>
#include <string>
#include <memory>

namespace ExchCXX {
namespace detail {

struct XCKernelImpl {

  using unique_me = std::unique_ptr< XCKernelImpl >;

  XCKernelImpl() = default;

  XCKernelImpl( XCKernelImpl&& )                 = delete;
  XCKernelImpl& operator=( const XCKernelImpl& ) = delete;
  XCKernelImpl& operator=( XCKernelImpl&& )      = delete;
  
  virtual ~XCKernelImpl() = default;

  unique_me clone() const { return clone_(); };

  bool is_lda()          const noexcept { return is_lda_();       }
  bool is_gga()          const noexcept { return is_gga_();       }
  bool is_mgga()         const noexcept { return is_mgga_();      }
  bool is_hyb()          const noexcept { return is_hyb_();       }
  bool is_polarized()    const noexcept { return is_polarized_(); }
  bool needs_laplacian() const noexcept { return needs_laplacian_();       }
  bool needs_tau() const noexcept { return needs_tau_();       }

  double hyb_exx() const noexcept { return hyb_exx_(); };

  bool supports_inc_interface() const noexcept {
    return supports_inc_interface_();
  }




  inline size_t rho_buffer_len( size_t npts ) const noexcept {
    return is_polarized() ? 2*npts : npts;
  }
  inline size_t sigma_buffer_len( size_t npts ) const noexcept {
    return is_lda() ? 0 : is_polarized() ? 3*npts : npts;
  }
  inline size_t lapl_buffer_len( size_t npts ) const noexcept {
    return needs_laplacian() ? rho_buffer_len(npts) : 0;
  }
  inline size_t tau_buffer_len( size_t npts ) const noexcept {
    return is_mgga() ? rho_buffer_len(npts) : 0;
  }

  inline size_t exc_buffer_len( size_t npts ) const noexcept {
    return npts;
  }
  inline size_t vrho_buffer_len( size_t npts ) const noexcept {
    return rho_buffer_len( npts );
  }
  inline size_t vsigma_buffer_len( size_t npts ) const noexcept {
    return sigma_buffer_len( npts );
  }
  inline size_t vlapl_buffer_len( size_t npts ) const noexcept {
    return lapl_buffer_len( npts );
  }
  inline size_t vtau_buffer_len( size_t npts ) const noexcept {
    return tau_buffer_len( npts );
  }

  template <typename... Args>
  void eval_exc( Args&&... args ) {
    eval_exc_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_inc( Args&&... args ) {
    eval_exc_inc_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc( Args&&... args ) {
    eval_exc_vxc_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_inc( Args&&... args ) {
    eval_exc_vxc_inc_( std::forward<Args>(args)... );
  }
 
       
    
  // Device Code
#ifdef EXCHCXX_ENABLE_DEVICE

  template <typename... Args>
  void eval_exc_device( Args&&... args ) {
    eval_exc_device_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_inc_device( Args&&... args ) {
    eval_exc_inc_device_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_device( Args&&... args ) {
    eval_exc_vxc_device_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_inc_device( Args&&... args ) {
    eval_exc_vxc_inc_device_( std::forward<Args>(args)... );
  }

#endif

protected:

  XCKernelImpl( const XCKernelImpl& ) = default;

private:

  virtual unique_me clone_() const = 0;

  virtual bool is_lda_()        const noexcept = 0;
  virtual bool is_gga_()        const noexcept = 0;
  virtual bool is_mgga_()       const noexcept = 0;
  virtual bool is_hyb_()        const noexcept = 0;
  virtual bool is_polarized_()  const noexcept = 0;
  virtual bool needs_laplacian_() const noexcept = 0;
  virtual bool needs_tau_() const noexcept = 0;

  virtual double hyb_exx_() const noexcept = 0;

  virtual bool supports_inc_interface_() const noexcept = 0;

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc_ )                 const = 0;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const = 0;
  virtual LDA_EXC_INC_GENERATOR( eval_exc_inc_ )         const = 0;
  virtual LDA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc_ )                 const = 0;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const = 0;
  virtual GGA_EXC_INC_GENERATOR( eval_exc_inc_ )         const = 0;
  virtual GGA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc_ )                 const = 0;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_ )         const = 0;
  virtual MGGA_EXC_INC_GENERATOR( eval_exc_inc_ )         const = 0;
  virtual MGGA_EXC_VXC_INC_GENERATOR( eval_exc_vxc_inc_ ) const = 0;

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  virtual LDA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const = 0;
  virtual LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const = 0;
  virtual LDA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const = 0;
  virtual LDA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const = 0;
  virtual GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const = 0;
  virtual GGA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const = 0;
  virtual GGA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_ )                 const = 0;
  virtual MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_ )         const = 0;
  virtual MGGA_EXC_INC_GENERATOR_DEVICE( eval_exc_inc_device_ )         const = 0;
  virtual MGGA_EXC_VXC_INC_GENERATOR_DEVICE( eval_exc_vxc_inc_device_ ) const = 0;

#endif
    
    


};

} // detail
} // ExchCXX

