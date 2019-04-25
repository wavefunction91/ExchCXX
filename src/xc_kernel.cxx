#include <exchcxx/xc_kernel.hpp>

namespace ExchCXX {


XCKernel::XCKernel( 
  const int kern, 
  const int spin_polar
) : polar_(spin_polar) {

  // Initialize XC Kernel using Libxc
  int info = xc_func_init( &kernel_, kern, spin_polar );

  assert( info == 0 );

} 

XCKernel::XCKernel( const XCKernel& other) :
  XCKernel( other.xc_info()->number, other.polar_ ){ };





// LDA interfaces
void XCKernel::eval_exc( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  assert( is_lda() );
  xc_lda_exc( &kernel_, N, rho, eps );

}


void XCKernel::eval_exc_vxc( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc 
) const {

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

  assert( is_gga() );
  xc_gga_exc_vxc( &kernel_, N, rho, sigma, eps, vrho, vsigma );

}

// TODO: GGA kxc interfaces  
};
