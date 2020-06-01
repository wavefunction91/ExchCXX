#pragma once

#include <exchcxx/device/cuda_type_fwd.hpp> // cuda_stream_t*
#include <exchcxx/util/exchcxx_macros.hpp>
//#include <exchcxx/xc_kernel.hpp>
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

  bool is_lda()        const noexcept { return is_lda_();       };
  bool is_gga()        const noexcept { return is_gga_();       };
  bool is_mgga()       const noexcept { return is_mgga_();      };
  bool is_hyb()        const noexcept { return is_hyb_();       };
  bool is_polarized()  const noexcept { return is_polarized_(); };

  double hyb_exx() const noexcept { return hyb_exx_(); };


  template <typename... Args>
  void eval_exc( Args&&... args ) {
    eval_exc_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc( Args&&... args ) {
    eval_exc_vxc_( std::forward<Args>(args)... );
  }
 
       
    
  // Device Code
#ifdef EXCHCXX_ENABLE_DEVICE

  template <typename... Args>
  void eval_exc_device( Args&&... args ) {
    eval_exc_device_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_device( Args&&... args ) {
    eval_exc_vxc_device_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_device_async( Args&&... args ) {
    eval_exc_device_async_( std::forward<Args>(args)... );
  }

  template <typename... Args>
  void eval_exc_vxc_device_async( Args&&... args ) {
    eval_exc_vxc_device_async_( std::forward<Args>(args)... );
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

  virtual double hyb_exx_() const noexcept = 0;

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc_ )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc_ ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc_ )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc_ ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc_ )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_ ) const = 0;

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc_device_ )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device_ ) const = 0;
  virtual LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async_ )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async_ ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc_device_ )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc_device_ ) const = 0;
  virtual GGA_EXC_GENERATOR_DEVICE( eval_exc_device_async_ )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async_ ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc_device_ )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_device_ ) const = 0;
  virtual MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_async_ )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async_ ) const = 0;

#endif
    
    


};

}; // detail
}; // ExchCXX

