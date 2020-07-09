#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <exchcxx/enums/enums.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/util/exchcxx_macros.hpp>

// Standard Libs
#include <vector>
#include <cassert>
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
  XCKernel( const Kernel kern, const Spin polar ) : 
    XCKernel( Backend::libxc, kern, polar ){ };

  XCKernel( impl_ptr&& ptr )                     ;
  XCKernel( const XCKernel& )                    ;
  XCKernel( XCKernel&&      )            noexcept;
  XCKernel& operator=( const XCKernel& )         ;
  XCKernel& operator=( XCKernel&&      ) noexcept;

  // Destroy interal Libxc data
  ~XCKernel() noexcept;



  bool is_lda()       const noexcept { return pimpl_->is_lda();       };
  bool is_gga()       const noexcept { return pimpl_->is_gga();       };
  bool is_mgga()      const noexcept { return pimpl_->is_mgga();      };
  bool is_hyb()       const noexcept { return pimpl_->is_hyb();       };
  bool is_polarized() const noexcept { return pimpl_->is_polarized(); };
  
  double hyb_exx() const noexcept { return pimpl_->hyb_exx(); }

  bool supports_inc_interface() const noexcept {
    return pimpl_->supports_inc_interface();
  }

  inline size_t rho_buffer_len( size_t npts ) const noexcept {
    return is_polarized() ? 2*npts : npts;
  }
  inline size_t sigma_buffer_len( size_t npts ) const noexcept {
    return is_lda() ? 0 : is_polarized() ? 3*npts : npts;
  }
  inline size_t lapl_buffer_len( size_t npts ) const noexcept {
    return is_mgga() ? rho_buffer_len(npts) : 0;
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

