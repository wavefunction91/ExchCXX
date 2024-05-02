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
#include <exchcxx/util/exchcxx_macros.hpp>
#include <exchcxx/exceptions/exchcxx_exception.hpp>

#include <numeric>

namespace ExchCXX {

class XCFunctional {

  using value_type = std::pair<double, XCKernel>;
private:

  std::vector< value_type > kernels_;

  inline bool sanity_check() const {

    // Must have one kernel
    if( not kernels_.size() ) return false;

    // Polarization is all or nothing
    int polar_one = kernels_.at(0).second.is_polarized();
    bool polar_all = std::all_of(
      kernels_.begin(), kernels_.end(),
      [&](const auto& a){ 
        return (int)a.second.is_polarized() == polar_one; 
      }
    ); 

    if( not polar_all ) return false;

    // If we made it, kernel is sane
    return true;

  }


  void throw_if_not_sane() const { EXCHCXX_BOOL_CHECK("Functional Not Sane", sanity_check() ); }


  inline bool supports_inc_interface() const noexcept {
    return std::all_of( 
      kernels_.begin(), kernels_.end(),
      [&](const auto& a){ 
        return a.second.supports_inc_interface();
      }
    ); 
  }

public:

  XCFunctional();
  explicit XCFunctional( const std::vector< XCKernel >& );
  explicit XCFunctional( const std::initializer_list< value_type >& list );
  explicit XCFunctional( const std::vector<value_type>& ks );
  explicit XCFunctional( std::vector<value_type>&& ks );

  XCFunctional( const Backend, const Functional, const Spin );
  XCFunctional( const Functional func, const Spin polar) :
    XCFunctional( Backend::libxc, func, polar) { };

  XCFunctional( const XCFunctional& )                    ;
  XCFunctional( XCFunctional&& )                 noexcept;
  XCFunctional& operator=( const XCFunctional& )         ;
  XCFunctional& operator=( XCFunctional&&      ) noexcept;



  inline bool is_lda() const {
    throw_if_not_sane();
    return std::all_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_lda(); }
    );
  }

  inline bool is_gga() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_gga(); }
    ) and not is_mgga();
  }

  inline bool is_mgga() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_mgga(); }
    );
  }

  inline bool is_polarized() const {
    throw_if_not_sane();
    // Polarization is all or nothing
    return kernels_.at(0).second.is_polarized(); 
  }

  inline bool is_hyb() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_hyb(); }
    );
  }

  inline bool is_epc() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.is_epc(); }
    );
  }

  inline bool needs_laplacian() const {
    throw_if_not_sane();
    return std::any_of( 
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.needs_laplacian(); }
    );
  }

  inline bool needs_tau() const {
    throw_if_not_sane();
    return std::any_of(
      kernels_.begin(), kernels_.end(),
      [](const auto& x) { return x.second.needs_tau(); }
    );
  }


  inline double hyb_exx() const {
    throw_if_not_sane();
    return std::accumulate( 
      kernels_.begin(), kernels_.end(), 0.,
      [](const auto& x, const auto &y) { 
        return x + y.second.hyb_exx(); 
      }
    );
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




  // LDA Interfaces
  LDA_EXC_GENERATOR(     eval_exc     ) const;
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // GGA Interfaces
  GGA_EXC_GENERATOR(     eval_exc     ) const;
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // mGGA interface
  MGGA_EXC_GENERATOR(     eval_exc     ) const;
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;


  // Device code
#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA Interfaces
  LDA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

  // GGA Interfaces
  GGA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

  // mGGA interface
  MGGA_EXC_GENERATOR_DEVICE(     eval_exc_device     ) const;
  MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device ) const;

#endif

};

} // namespace ExchCXX

