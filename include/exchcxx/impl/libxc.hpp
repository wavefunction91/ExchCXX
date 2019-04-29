#ifndef __INCLUDED_IMPL_LIBXC_HPP__
#define __INCLUDED_IMPL_LIBXC_HPP__

#include <exchcxx/impl/xc_kernel.hpp>
#include <cassert>
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


  // LDA interfaces
  void eval_exc_( 
    const int     N, 
    const double* rho, 
    double*       eps 
  ) const override; 

/*
  void eval_vxc_( 
    const int     N, 
    const double* rho, 
    double*       vxc 
  ) const override; 
*/

  void eval_exc_vxc_( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const override; 

/*
  void eval_fxc_( 
    const int     N, 
    const double* rho, 
    double*       fxc 
  ) const override; 

  void eval_kxc_( 
    const int     N, 
    const double* rho, 
    double*       kxc 
  ) const override; 

  void eval_( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc 
  ) const override; 
*/
    
    
  // GGA interface
  void eval_exc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const override; 

/*
  void eval_vxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma
  ) const override;
*/
  
  void eval_exc_vxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const override;

/*
  void eval_fxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const override;

  void eval_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const override;
*/


  // mGGA interface
  void eval_exc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const override; 

  
  void eval_exc_vxc_( 
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
  ) const override;

  // TODO: mGGA fxc/kxc interfaces  
protected:

  int          polar_;  ///< Spin polarization
  xc_func_type kernel_; ///< Libxc kernel definition

  bool initialized_ = false;
  void throw_if_uninitialized() const { assert( initialized_ ); }


  auto xc_info() const { 
    throw_if_uninitialized();
    return kernel_.info; 
  };

  LibxcKernelImpl( xc_func_type kern, const int spin_polar, 
    const bool init );

public:

  LibxcKernelImpl() = delete;
  
  LibxcKernelImpl( const std::string& kname, const bool spin_polar );
  LibxcKernelImpl( const int kern, const int spin_polar );
  LibxcKernelImpl( const LibxcKernelImpl& );

  // Destroy interal Libxc data
  ~LibxcKernelImpl() noexcept;




    
};

}; // namespace detail
}; // namespace ExchCXX

#endif
