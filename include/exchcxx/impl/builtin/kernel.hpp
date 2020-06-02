#pragma once

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/builtin/fwd.hpp>
#include <type_traits>

namespace ExchCXX {



class BuiltinKernel {

  XCKernel::Spin polar_;

public:

  explicit BuiltinKernel( XCKernel::Spin p ) : polar_(p) { }
  virtual ~BuiltinKernel() = default;

  virtual bool is_lda()       const noexcept = 0;
  virtual bool is_gga()       const noexcept = 0;
  virtual bool is_mgga()      const noexcept = 0;
  virtual bool is_hyb()       const noexcept = 0;
  virtual double hyb_exx()    const noexcept = 0;

  inline bool is_polarized() const noexcept { 
    return polar_ == XCKernel::Spin::Polarized; 
  };

  inline auto polar() const noexcept { return polar_; }

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc )         const;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc )         const;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc )         const;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const;

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc_device )         const;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const;
  virtual LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const;
  virtual LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc_device )         const;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const;
  virtual GGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const;
  virtual GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc_device )         const;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const;
  virtual MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const;
  virtual MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const;

#endif


};



namespace detail {

template <typename KernelType>
struct BuiltinKernelImpl<
  KernelType,
  std::enable_if_t<kernel_traits<KernelType>::is_lda>
> : public BuiltinKernel {

  using traits = kernel_traits<KernelType>;

  BuiltinKernelImpl( XCKernel::Spin p ) : BuiltinKernel(p) { }
  virtual ~BuiltinKernelImpl() noexcept = default;

  inline bool is_hyb()  const noexcept override { return traits::is_hyb;  }
  inline bool is_lda()  const noexcept override { return traits::is_lda;  }
  inline bool is_gga()  const noexcept override { return traits::is_gga;  }
  inline bool is_mgga() const noexcept override { return traits::is_mgga; }

  double hyb_exx() const noexcept override {
    return traits::exx_coeff;
  }

  LDA_EXC_GENERATOR( eval_exc ) const override {

    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_unpolar( rho[i], eps[i] );

  }

  LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override {

    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_vxc_unpolar( rho[i], eps[i], vxc[i] );

  }

};



template <typename KernelType>
struct BuiltinKernelImpl<
  KernelType,
  std::enable_if_t<kernel_traits<KernelType>::is_gga>
> : public BuiltinKernel {

  using traits = kernel_traits<KernelType>;

  BuiltinKernelImpl( XCKernel::Spin p ) : BuiltinKernel(p) { }
  virtual ~BuiltinKernelImpl() noexcept = default;

  inline bool is_hyb()  const noexcept override { return traits::is_hyb;  }
  inline bool is_lda()  const noexcept override { return traits::is_lda;  }
  inline bool is_gga()  const noexcept override { return traits::is_gga;  }
  inline bool is_mgga() const noexcept override { return traits::is_mgga; }

  double hyb_exx() const noexcept override {
    return traits::exx_coeff;
  }

  GGA_EXC_GENERATOR( eval_exc ) const override {

    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_unpolar( rho[i], sigma[i], eps[i] );

  }

  GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override {

    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_vxc_unpolar( rho[i], sigma[i], eps[i], vrho[i], vsigma[i] );

  }

};

}





}
