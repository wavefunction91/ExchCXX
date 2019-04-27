#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {

XCKernel::XCKernel( 
  const int kern, 
  const int spin_polar
) : polar_(spin_polar) {

  // Initialize XC Kernel using Libxc
  int info = xc_func_init( &kernel_, kern, spin_polar );

  assert( info == 0 );

  initialized_ = true;
} 

XCKernel::XCKernel( const XCKernel& other ) :
  XCKernel( other.xc_info()->number, other.polar_ ){ };

XCKernel::XCKernel( XCKernel&& other ) noexcept :
  XCKernel( other.kernel_, other.polar_, other.initialized_) { 
  other.initialized_ = false; // Avoid double destruction
};

XCKernel& XCKernel::operator=( XCKernel&& other ) noexcept {
  kernel_     = other.kernel_;
  polar_      = other.polar_;
  initialized_= other.initialized_;

  other.initialized_ = false; // Avoid double destruction

  return *this;
}

XCKernel& XCKernel::operator=( const XCKernel& other ) {
  return *this = std::move( XCKernel(other) );
}


XCKernel::~XCKernel() noexcept {
  if( initialized_ ) xc_func_end( &kernel_ );
}



// LDA interfaces
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc( &kernel_, N, rho, eps );

}


void XCKernel::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc 
) const {

  throw_if_uninitialized();
  assert( is_lda() );
  xc_lda_exc_vxc( &kernel_, N, rho, eps, vxc );

}

// TODO: LDA kxc interfaces

// GGA interface
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps
) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc( &kernel_, N, rho, sigma, eps );

}


void XCKernel::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  double*       vrho,
  double*       vsigma
) const {

  throw_if_uninitialized();
  assert( is_gga() );
  xc_gga_exc_vxc( &kernel_, N, rho, sigma, eps, vrho, vsigma );

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

  throw_if_uninitialized();
  assert( is_mgga() );
  xc_mgga_exc( &kernel_, N, rho, sigma, lapl, tau, eps );

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

  throw_if_uninitialized();
  assert( is_gga() );
  xc_mgga_exc_vxc( &kernel_, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau );

}
};
