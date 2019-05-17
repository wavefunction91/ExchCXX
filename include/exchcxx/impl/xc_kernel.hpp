#ifndef __INCLUDED_IMPL_XC_KERNEL_HPP__
#define __INCLUDED_IMPL_XC_KERNEL_HPP__

#include <exchcxx/xc_kernel.hpp>
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


  // LDA interfaces
    
  void eval_exc( 
    const int N, 
    const double* rho, 
    double* eps 
  ) const { eval_exc_(N, rho, eps); }

/*
  void eval_vxc( 
    const int     N, 
    const double* rho, 
    double*       vxc 
  ) const { eval_vxc_(N, rho, vxc); }; 
*/

  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const { eval_exc_vxc_(N, rho, eps, vxc); }; 

/*
  void eval_fxc( 
    const int     N, 
    const double* rho, 
    double*       fxc 
  ) const { eval_fxc_(N, rho, fxc); }; 

  void eval_kxc( 
    const int     N, 
    const double* rho, 
    double*       kxc 
  ) const { eval_fxc_(N, rho, kxc); }; 

  void eval( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc 
  ) const { eval_(N,rho,eps,vxc,fxc,kxc); }; 
*/

 
  // GGA Interfaces

  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const { eval_exc_(N,rho,sigma,eps); }; 

/*
  void eval_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma
  ) const { eval_vxc_(N,rho,sigma,vrho,vsigma); }; 
*/
  
  void eval_exc_vxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const { eval_exc_vxc_( N, rho, sigma, eps, vrho, vsigma ); }; 

/*
  void eval_fxc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const { eval_fxc_(N,rho,sigma,v2rho2,v2rhosigma,v2sigma2); }; 

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
  ) const { eval_(N,rho,sigma,eps,vrho,vsigma,v2rho2,v2rhosigma,
                  v2sigma2); }; 
*/

  // TODO: GGA kxc interface
       
  // mGGA interface
    
  void eval_exc( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const { eval_exc_(N,rho,sigma,lapl,tau,eps); } 

  
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
  ) const { eval_exc_vxc_(N,rho,sigma,lapl,tau,eps,vrho,vsigma,
                          vlapl,vtau); } 

  // TODO: mGGA fxc/kxc interfaces  
    
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

  // LDA interfaces
  virtual void eval_exc_( 
    const int     N, 
    const double* rho, 
    double*       eps 
  ) const = 0; 

/*
  virtual void eval_vxc_( 
    const int     N, 
    const double* rho, 
    double*       vxc 
  ) const = 0; 
*/

  virtual void eval_exc_vxc_( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc 
  ) const = 0; 

/*
  virtual void eval_fxc_( 
    const int     N, 
    const double* rho, 
    double*       fxc 
  ) const = 0; 

  virtual void eval_kxc_( 
    const int     N, 
    const double* rho, 
    double*       kxc 
  ) const = 0; 

  virtual void eval_( 
    const int     N, 
    const double* rho, 
    double*       eps, 
    double*       vxc,
    double*       fxc,
    double*       kxc 
  ) const = 0; 
*/
    
    
  // GGA interface
  virtual void eval_exc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps
  ) const = 0; 

/*
  virtual void eval_vxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       vrho,
    double*       vsigma
  ) const = 0;
*/
  
  virtual void eval_exc_vxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma
  ) const = 0;

/*
  virtual void eval_fxc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const = 0;

  virtual void eval_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    double*       eps,
    double*       vrho,
    double*       vsigma,
    double*       v2rho2, 
    double*       v2rhosigma, 
    double*       v2sigma2
  ) const = 0;
*/


  // mGGA interface
  virtual void eval_exc_( 
    const int     N, 
    const double* rho, 
    const double* sigma, 
    const double* lapl, 
    const double* tau, 
    double*       eps
  ) const = 0; 

  
  virtual void eval_exc_vxc_( 
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
  ) const = 0;

  // TODO: mGGA fxc/kxc interfaces  
};

}; // detail
}; // ExchCXX

#endif
