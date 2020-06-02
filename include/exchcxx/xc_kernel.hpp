#ifndef __INCLUDED_XC_KERNEL_HPP__
#define __INCLUDED_XC_KERNEL_HPP__

//#include "impl/xc_kernel_fwd.hpp" // XCKernelImpl
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/device/cuda_type_fwd.hpp> // cuda_stream_t*

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

  enum Backend {
    libxc,
    builtin
  };

  enum Spin {
    Polarized,
    Unpolarized
  };

  #include "kernels/kernels.hpp"

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



  template <typename... Args>
  void eval_exc( Args&&... args ) const {
    pimpl_->eval_exc( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc( Args&&... args ) const {
    pimpl_->eval_exc_vxc( std::forward<Args>(args)... );
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
  void eval_exc_device_async( Args&&... args ) const {
    pimpl_->eval_exc_device_async( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_device_async( Args&&... args ) const {
    pimpl_->eval_exc_vxc_device_async( std::forward<Args>(args)... );
  }

#endif



};









}

#endif
