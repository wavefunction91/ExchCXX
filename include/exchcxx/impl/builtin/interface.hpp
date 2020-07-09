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
  bool is_polarized_() const noexcept override;
  double hyb_exx_()    const noexcept override;

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
  BuiltinKernelInterface( BuiltinKernelInterface&& ) noexcept = default;

  // Destroy interal Builtin data
  ~BuiltinKernelInterface() noexcept;
    
};

}
}
