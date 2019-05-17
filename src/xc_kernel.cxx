#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>
#include <exchcxx/factory/xc_kernel.hpp>

namespace ExchCXX {

XCKernel::XCKernel( 
  const Backend backend, 
  const std::string& kname, 
  const Spin polar) : 
XCKernel( std::move(libxc_kernel_factory( kname, polar )) ) { }
  

XCKernel::XCKernel( impl_ptr&& ptr ) :
  pimpl_(std::move(ptr)) { };

XCKernel::XCKernel( const XCKernel& other ) :
  pimpl_(other.pimpl_->clone()) { };

XCKernel::XCKernel( XCKernel&& other ) noexcept = default;

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept =default;


XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = std::move( XCKernel(other) );
}


XCKernel::~XCKernel() = default;

bool XCKernel::is_lda()       const noexcept { return pimpl_->is_lda();       };
bool XCKernel::is_gga()       const noexcept { return pimpl_->is_gga();       };
bool XCKernel::is_mgga()      const noexcept { return pimpl_->is_mgga();      };
bool XCKernel::is_hyb()       const noexcept { return pimpl_->is_hyb();       };
bool XCKernel::is_polarized() const noexcept { return pimpl_->is_polarized(); };

double XCKernel::hyb_exx() const noexcept { return pimpl_->hyb_exx(); }


// LDA interfaces
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  pimpl_->eval_exc( N, rho, eps );  

}


void XCKernel::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc 
) const {

  pimpl_->eval_exc_vxc( N, rho, eps, vxc );  

}

// TODO: LDA kxc interfaces

// GGA interface
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps
) const {

  pimpl_->eval_exc( N, rho, sigma, eps );  

}


void XCKernel::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  double*       vrho,
  double*       vsigma
) const {

  pimpl_->eval_exc_vxc( N, rho, sigma, eps, vrho, vsigma );  

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps
) const {

  pimpl_->eval_exc( N, rho, sigma, lapl, tau, eps );  

}


void XCKernel::eval_exc_vxc( 
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
) const {

  pimpl_->eval_exc_vxc( N, rho, sigma, lapl, tau, 
    eps, vrho, vsigma, vlapl, vtau );  

}
};
