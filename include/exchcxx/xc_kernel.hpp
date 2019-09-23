#ifndef __INCLUDED_XC_KERNEL_HPP__
#define __INCLUDED_XC_KERNEL_HPP__

#include "impl/xc_kernel_fwd.hpp" // XCKernelImpl
#include "device/cuda_type_fwd.hpp" // cuda_stream_t*

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
    libxc
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



  /**
   *  \brief Determine if kernel uses the local density 
   *  approximation (LDA).
   *
   *  No throw guarantee.
   *
   *  \return true if kernel is LDA, false otherwise
   */ 
  bool is_lda() const noexcept;

  /**
   *  \brief Determine if kernel uses the generalized
   *  gradient approximation (GGA)
   *
   *  No throw guarantee.
   *
   *  \return true if kernel is GGA, false otherwise
   */ 
  bool is_gga() const noexcept;

  /**
   *  \brief Determine if kernel uses the meta generalized
   *  gradient approximation (mGGA)
   *
   *  No throw guarantee.
   *
   *  \return true if kernel is mGGA, false otherwise
   */ 
  bool is_mgga() const noexcept;

  /**
   *  \brief Determine if kernel is a hybrid kernel. 
   *
   *  No throw guarantee.
   *
   *  \return true if kernel is hybrid, false otherwise
   */ 
  bool is_hyb() const noexcept;

  /**
   *  \brief Determine if kernel was initialized as spin
   *  polarized
   *
   *  No throw guarantee.
   *
   *  \return true if kernel is spin polarized, false otherwise
   */ 
  bool is_polarized() const noexcept;

  /**
   *  \brief Determine the exact (HF) exchange coefficient
   *  for the kernel
   *
   *  No throw guarantee.
   *
   *  \return the HF exchange coefficient for hybrid functionals,
   *  0. otherwise.
   */ 
  double hyb_exx() const noexcept;



  // LDA interfaces
    
  /**
   *  \brief Evaluate the XC energy density for an LDA kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] eps  Evalaution of the XC energy density for the
   *                   specified density points. An array of length
   *                   N.
   */
  void eval_exc( 
    const int     N, 
    const double* rho, 
    double*       eps 
  ) const; 

  /**
   *  \brief Evaluate the XC potential for an LDA kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] vxc  Evalaution of the XC potential for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 2N if initialized as spin
   *                   polarized.
   */
  void eval_vxc( 
    const int     N, 
    const double* rho, 
    double*       vxc 
  ) const; 

  /**
   *  \brief Evaluate the XC energy density and the XC potential 
   *  for an LDA kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] eps  Evalaution of the XC energy density for the
   *                   specified density points. An array of length
   *                   N.
   *  \param[out] vxc  Evalaution of the XC potential for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 2N if initialized as spin
   *                   polarized.
   */
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const; 

  /**
   *  \brief Evaluate the functional second derivative for
   *  an LDA XC kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] fxc  Evalaution of the functional second 
   *                   derivative of the XC kernel for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 3N if initialized as spin
   *                   polarized.
   */
  void eval_fxc( 
    const int     N, 
    const double* rho, 
    double*       fxc 
  ) const; 

  /**
   *  \brief Evaluate the functional third derivative for
   *  an LDA XC kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] kxc  Evalaution of the functional third
   *                   derivative of the XC kernel for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 4N if initialized as spin
   *                   polarized.
   */
  void eval_kxc( 
    const int     N, 
    const double* rho, 
    double*       kxc 
  ) const; 


  /**
   *  \brief Evaluate the XC energy density, XC potential,
   *  and the functional second and third derivatives
   *  for an LDA kernel.
   *
   *  Only valid for LDA kernels, throws for non-LDA kernels.
   *  Strong throw guarantee.
   *
   *  \param[in]  N    Number of density points to evaluate kernel
   *  \parap[in]  rho  Density points to evaluate the kernel. An 
   *                   array of N density evaluations if 
   *                   initialized as spin unpolarized, 2N if
   *                   initialized as spin polarized.
   *  \param[out] eps  Evalaution of the XC energy density for the
   *                   specified density points. An array of length
   *                   N.
   *  \param[out] vxc  Evalaution of the XC potential for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 2N if initialized as spin
   *                   polarized.
   *  \param[out] fxc  Evalaution of the functional second 
   *                   derivative of the XC kernel for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 3N if initialized as spin
   *                   polarized.
   *  \param[out] kxc  Evalaution of the functional third
   *                   derivative of the XC kernel for the
   *                   specified density points. An array of length
   *                   N if kernel was initialized as spin 
   *                   unpolarized, 4N if initialized as spin
   *                   polarized.
   */
  void eval( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc 
  ) const; 
    
    
  // GGA interface
  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const; 

  void eval_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma
  ) const;
  
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const;

  void eval_fxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const;

  void eval( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const;


  // mGGA interface
  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const; 

  
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       vlapl,
    double*       vtau
  ) const;

  // TODO: mGGA fxc/kxc interfaces  
    

  // Device code
#ifdef EXCHCXX_ENABLE_DEVICE
  
  // LDA interfaces
  void eval_exc_device( 
    const int     N, 
    const double* rho, 
    double*       eps 
  ) const; 
    
  void eval_vxc_device( 
    const int     N, 
    const double* rho, 
    double*       vxc 
  ) const; 

  void eval_exc_vxc_device( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const; 
    
  void eval_fxc_device( 
    const int     N, 
    const double* rho, 
    double*       fxc 
  ) const; 

  void eval_kxc_device( 
    const int     N, 
    const double* rho, 
    double*       kxc 
  ) const; 

  void eval_device( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc 
  ) const; 

  // GGA interface
  void eval_exc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const; 

  void eval_vxc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma
  ) const;
  
  void eval_exc_vxc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const;

  void eval_fxc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const;

  void eval_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const;


  // mGGA interface
  void eval_exc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const; 

  
  void eval_exc_vxc_device( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       vlapl,
    double*       vtau
  ) const;

  // TODO: mGGA fxc/kxc interfaces  
    
    
  // LDA interfaces
  void eval_exc_device_async( 
    const int     N, 
    const double* rho, 
    double*       eps,
    device::cuda_stream_t* stream
  ) const; 
    
  void eval_vxc_device_async( 
    const int     N, 
    const double* rho, 
    double*       vxc, 
    device::cuda_stream_t* stream
  ) const; 

  void eval_exc_vxc_device_async( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc, 
    device::cuda_stream_t* stream
  ) const; 
    
  void eval_fxc_device_async( 
    const int     N, 
    const double* rho, 
    double*       fxc, 
    device::cuda_stream_t* stream
  ) const; 

  void eval_kxc_device_async( 
    const int     N, 
    const double* rho, 
    double*       kxc,
    device::cuda_stream_t* stream
  ) const; 

  void eval_device_async( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc, 
    device::cuda_stream_t* stream
  ) const; 

  // GGA interface
  void eval_exc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    device::cuda_stream_t* stream
  ) const; 

  void eval_vxc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma,
    device::cuda_stream_t* stream
  ) const;
  
  void eval_exc_vxc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    device::cuda_stream_t* stream
  ) const;

  void eval_fxc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2,
    device::cuda_stream_t* stream
  ) const;

  void eval_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2,
    device::cuda_stream_t* stream
  ) const;


  // mGGA interface
  void eval_exc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps,
    device::cuda_stream_t* stream
  ) const; 

  
  void eval_exc_vxc_device_async( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       vlapl,
    double*       vtau,
    device::cuda_stream_t* stream
  ) const;

  // TODO: mGGA fxc/kxc interfaces  
#endif
};









};

#endif
