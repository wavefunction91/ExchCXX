#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif


#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/exceptions/exchcxx_exception.hpp>
#include <xc.h> // Libxc

namespace ExchCXX {
namespace detail {

class LibxcKernelImpl : public XCKernelImpl {

  using unique_me = XCKernelImpl::unique_me;

  unique_me clone_() const override;

  bool is_lda_()       const noexcept override;
  bool is_gga_()       const noexcept override;
  bool is_mgga_()      const noexcept override;
  bool is_hyb_()       const noexcept override;
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

protected:

  int          polar_;  ///< Spin polarization
  xc_func_type kernel_; ///< Libxc kernel definition

  bool initialized_ = false;
  void throw_if_uninitialized() const { EXCHCXX_BOOL_CHECK("Kernel Uninitialized", initialized_ ); }


  auto xc_info() const { 
    throw_if_uninitialized();
    return kernel_.info; 
  };

  LibxcKernelImpl( xc_func_type kern, const int spin_polar, 
    const bool init );

public:

  LibxcKernelImpl() = delete;
  
  LibxcKernelImpl( const Kernel kern, const Spin polar);
  LibxcKernelImpl( const int kern, const int spin_polar );
  LibxcKernelImpl( const LibxcKernelImpl& );

  // Destroy interal Libxc data
  ~LibxcKernelImpl() noexcept;




    
};

} // namespace detail

} // namespace ExchCXX

