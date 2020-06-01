#pragma once

#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <cmath>

#include <exchcxx/impl/builtin/constants.hpp>

namespace ExchCXX {


inline static void disabled_lda_interface() {
  throw std::runtime_error("LDA Interface is disabled for the specified kernel");
}

inline static void disabled_gga_interface() {
  throw std::runtime_error("GGA Interface is disabled for the specified kernel");
}

inline static void disabled_mgga_interface() {
  throw std::runtime_error("MGGA Interface is disabled for the specified kernel");
}











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
  virtual LDA_EXC_GENERATOR( eval_exc )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const = 0;

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  virtual LDA_EXC_GENERATOR( eval_exc_device )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const = 0;
  virtual LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const = 0;
  virtual LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const = 0;

  // GGA interface
  virtual GGA_EXC_GENERATOR( eval_exc_device )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const = 0;
  virtual GGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const = 0;
  virtual GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const = 0;

  // MGGA interface
  virtual MGGA_EXC_GENERATOR( eval_exc_device )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const = 0;
  virtual MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const = 0;
  virtual MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const = 0;

#endif

};


template <typename T>
struct PureKernel : virtual T {

  bool is_hyb()     const noexcept override { return false; };
  double hyb_exx()  const noexcept override { return 0.; };

  virtual ~PureKernel() = default;
};


template <typename T>
class HybKernel : virtual T {

  double exx_;

public:

  bool is_hyb()     const noexcept override { return true; };
  double hyb_exx()  const noexcept override { return exx_; };

  HybKernel( double exx ) : exx_(exx) { }
  virtual ~HybKernel() = default;

};

template <typename T>
struct LDAKernel : virtual T {

  bool is_lda()  const noexcept override { return true;  };
  bool is_gga()  const noexcept override { return false; };
  bool is_mgga() const noexcept override { return false; };

  virtual ~LDAKernel() = default;

  // DISABLE GGA / MGGA interfaces

  // GGA interface
  GGA_EXC_GENERATOR( eval_exc )         const override { disabled_gga_interface(); };
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_gga_interface(); };

  // MGGA interface
  MGGA_EXC_GENERATOR( eval_exc )         const override { disabled_mgga_interface(); };
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_mgga_interface(); };

#ifdef EXCHCXX_ENABLE_DEVICE

  // GGA interface
  GGA_EXC_GENERATOR( eval_exc_device )         const override { disabled_gga_interface(); }
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_gga_interface(); }
  GGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_gga_interface(); }
  GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_gga_interface(); }

  // MGGA interface
  MGGA_EXC_GENERATOR( eval_exc_device )         const override { disabled_mgga_interface(); }
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_mgga_interface(); }
  MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_mgga_interface(); }
  MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_mgga_interface(); }

#endif

};

template <typename T>
struct GGAKernel : virtual T {

  bool is_lda()  const noexcept override { return false; };
  bool is_gga()  const noexcept override { return true;  };
  bool is_mgga() const noexcept override { return false; };

  virtual ~GGAKernel() = default;

  // DISABLE LDA / MGGA interfaces

  // LDA interface
  LDA_EXC_GENERATOR( eval_exc )         const override { disabled_lda_interface(); };
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_lda_interface(); };

  // MGGA interface
  MGGA_EXC_GENERATOR( eval_exc )         const override { disabled_mgga_interface(); };
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_mgga_interface(); };

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  LDA_EXC_GENERATOR( eval_exc_device )         const override { disabled_lda_interface(); }
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_lda_interface(); }
  LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_lda_interface(); }
  LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_lda_interface(); }

  // MGGA interface
  MGGA_EXC_GENERATOR( eval_exc_device )         const override { disabled_mgga_interface(); }
  MGGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_mgga_interface(); }
  MGGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_mgga_interface(); }
  MGGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_mgga_interface(); }

#endif

};

template <typename T>
struct MGGAKernel : virtual T {

  bool is_lda()  const noexcept override { return false; };
  bool is_gga()  const noexcept override { return false; };
  bool is_mgga() const noexcept override { return true;  };

  virtual ~MGGAKernel() = default;

  // DISABLE LDA / GGA interfaces

  // LDA interface
  LDA_EXC_GENERATOR( eval_exc )         const override { disabled_lda_interface(); };
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_lda_interface(); };

  // GGA interface
  GGA_EXC_GENERATOR( eval_exc )         const override { disabled_gga_interface(); };
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override { disabled_gga_interface(); };

#ifdef EXCHCXX_ENABLE_DEVICE

  // LDA interface
  LDA_EXC_GENERATOR( eval_exc_device )         const override { disabled_lda_interface(); }
  LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_lda_interface(); }
  LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_lda_interface(); }
  LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_lda_interface(); }

  // GGA interface
  GGA_EXC_GENERATOR( eval_exc_device )         const override { disabled_gga_interface(); }
  GGA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override { disabled_gga_interface(); }
  GGA_EXC_GENERATOR_DEVICE( eval_exc_device_async )         const override { disabled_gga_interface(); }
  GGA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override { disabled_gga_interface(); }

#endif

};



template <typename KernelType>
struct kernel_traits;


struct BuiltinSlaterExchange;

template <>
struct kernel_traits<BuiltinSlaterExchange> {

  static void eval_exc_unpolar( double rho, double& eps ) {

    constexpr double alpha = 1.;

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;

    double t11 = std::pow( rho, constants::m_third );
    double t13 = t5 * t10 * t11;

    eps = - 3. / 16. * t13;

  }



  static void eval_exc_vxc_unpolar( double rho, double& eps, double& vxc ) {

    constexpr double alpha = 1.;

    constexpr double t1 = constants::m_cbrt_3;
    constexpr double t6 = constants::m_cbrt_4;
    constexpr double t8 = constants::m_cbrt_2;
    constexpr double t4 = constants::m_cbrt_one_ov_pi;
    constexpr double t5 = alpha * t1 * t4;
    constexpr double t7 = t6 * t6;
    constexpr double t9 = t8 * t8;
    constexpr double t10 = t7 * t9;

    double t11 = std::pow( rho, constants::m_third );
    double t13 = t5 * t10 * t11;

    eps = - 3. / 16. * t13;
    vxc = -t13 / 4.;

  }

};





struct BuiltinSlaterExchange : 
  public virtual BuiltinKernel, 
  public virtual PureKernel<BuiltinKernel>,
  public virtual LDAKernel<BuiltinKernel> {

  using traits = kernel_traits<BuiltinSlaterExchange>;

  BuiltinSlaterExchange( XCKernel::Spin p ) :
    BuiltinKernel(p) { }
  
  virtual ~BuiltinSlaterExchange() = default;

  inline LDA_EXC_GENERATOR( eval_exc ) const override {
    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_unpolar( rho[i], eps[i] );
  }

  inline LDA_EXC_VXC_GENERATOR( eval_exc_vxc ) const override {
    for( size_t i = 0; i < N; ++i )
      traits::eval_exc_vxc_unpolar( rho[i], eps[i], vxc[i] );
  }





  inline LDA_EXC_GENERATOR( eval_exc_device ) const override {
    throw std::runtime_error("NYI");
  }
  inline LDA_EXC_VXC_GENERATOR( eval_exc_vxc_device ) const override {
    throw std::runtime_error("NYI");
  }

  inline LDA_EXC_GENERATOR_DEVICE( eval_exc_device_async ) const override {
    throw std::runtime_error("NYI");
  }
  inline LDA_EXC_VXC_GENERATOR_DEVICE( eval_exc_vxc_device_async ) const override {
    throw std::runtime_error("NYI");
  }
};



namespace detail {

class BuiltinKernelInterface : public XCKernelImpl {

  using unique_me = XCKernelImpl::unique_me;

  unique_me clone_() const override;

  XCKernel::Kernel whatami_;

  std::unique_ptr<BuiltinKernel> impl_;

  inline bool is_lda_()       const noexcept override {
    return impl_->is_lda();
  };
  inline bool is_gga_()       const noexcept override {
    return impl_->is_gga();
  }
  inline bool is_mgga_()      const noexcept override {
    return impl_->is_mgga();
  }
  inline bool is_hyb_()       const noexcept override {
    return impl_->is_hyb();
  }
  inline bool is_polarized_() const noexcept override {
    return impl_->is_polarized();
  }
  inline double hyb_exx_()    const noexcept override {
    return impl_->hyb_exx();
  }

  inline FORWARD_XC_ARGS( LDA, EXC, eval_exc_, impl_->eval_exc, const override );
  inline FORWARD_XC_ARGS( LDA, EXC_VXC, eval_exc_vxc_, impl_->eval_exc_vxc, const override );

  inline FORWARD_XC_ARGS( GGA, EXC, eval_exc_, impl_->eval_exc, const override );
  inline FORWARD_XC_ARGS( GGA, EXC_VXC, eval_exc_vxc_, impl_->eval_exc_vxc, const override );

  inline FORWARD_XC_ARGS( MGGA, EXC, eval_exc_, impl_->eval_exc, const override );
  inline FORWARD_XC_ARGS( MGGA, EXC_VXC, eval_exc_vxc_, impl_->eval_exc_vxc, const override );

#ifdef EXCHCXX_ENABLE_DEVICE

  inline FORWARD_XC_ARGS( LDA, EXC, eval_exc_device_, impl_->eval_exc_device, const override );
  inline FORWARD_XC_ARGS( LDA, EXC_VXC, eval_exc_vxc_device_, impl_->eval_exc_vxc_device, const override );

  inline FORWARD_XC_ARGS( GGA, EXC, eval_exc_device_, impl_->eval_exc_device, const override );
  inline FORWARD_XC_ARGS( GGA, EXC_VXC, eval_exc_vxc_device_, impl_->eval_exc_vxc_device, const override );

  inline FORWARD_XC_ARGS( MGGA, EXC, eval_exc_device_, impl_->eval_exc_device, const override );
  inline FORWARD_XC_ARGS( MGGA, EXC_VXC, eval_exc_vxc_device_, impl_->eval_exc_vxc_device, const override );

  inline FORWARD_XC_ARGS_DEVICE( LDA, EXC, eval_exc_device_async_, impl_->eval_exc_device_async, const override );
  inline FORWARD_XC_ARGS_DEVICE( LDA, EXC_VXC, eval_exc_vxc_device_async_, impl_->eval_exc_vxc_device_async, const override );

  inline FORWARD_XC_ARGS_DEVICE( GGA, EXC, eval_exc_device_async_, impl_->eval_exc_device_async, const override );
  inline FORWARD_XC_ARGS_DEVICE( GGA, EXC_VXC, eval_exc_vxc_device_async_, impl_->eval_exc_vxc_device_async, const override );

  inline FORWARD_XC_ARGS_DEVICE( MGGA, EXC, eval_exc_device_async_, impl_->eval_exc_device_async, const override );
  inline FORWARD_XC_ARGS_DEVICE( MGGA, EXC_VXC, eval_exc_vxc_device_async_, impl_->eval_exc_vxc_device_async, const override );

#endif

  //BuiltinKernelInterface( std::unique_ptr<BuiltinKernel>&& ptr ) noexcept ;

public:

  BuiltinKernelInterface() = delete;
  
  BuiltinKernelInterface( XCKernel::Kernel kern, XCKernel::Spin p);
  BuiltinKernelInterface( const BuiltinKernelInterface& );
  BuiltinKernelInterface( BuiltinKernelInterface&& ) noexcept = default;

  // Destroy interal Builtin data
  ~BuiltinKernelInterface() noexcept;
    
};

}
}
