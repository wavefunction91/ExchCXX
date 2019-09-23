#include <exchcxx/xc_kernel.hpp>
#include <exchcxx/impl/xc_kernel.hpp>

namespace ExchCXX {

// LDA interfaces
void XCKernel::eval_exc_device( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  pimpl_->eval_exc_device( N, rho, eps );  

}


void XCKernel::eval_exc_vxc_device( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc 
) const {

  pimpl_->eval_exc_vxc_device( N, rho, eps, vxc );  

}

// TODO: LDA kxc interfaces

// GGA interface
void XCKernel::eval_exc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps
) const {

  pimpl_->eval_exc_device( N, rho, sigma, eps );  

}


void XCKernel::eval_exc_vxc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  double*       vrho,
  double*       vsigma
) const {

  pimpl_->eval_exc_vxc_device( N, rho, sigma, eps, vrho, vsigma );  

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void XCKernel::eval_exc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps
) const {

  pimpl_->eval_exc_device( N, rho, sigma, lapl, tau, eps );  

}


void XCKernel::eval_exc_vxc_device( 
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

  pimpl_->eval_exc_vxc_device( N, rho, sigma, lapl, tau, 
    eps, vrho, vsigma, vlapl, vtau );  

}


// LDA interfaces
void XCKernel::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  device::cuda_stream_t* stream
) const {

  pimpl_->eval_exc_device_async( N, rho, eps, stream );  

}


void XCKernel::eval_exc_vxc_device_async( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc, 
  device::cuda_stream_t* stream
) const {

  pimpl_->eval_exc_vxc_device_async( N, rho, eps, vxc, stream );  

}

// TODO: LDA kxc interfaces

// GGA interface
void XCKernel::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  device::cuda_stream_t* stream
) const {

  pimpl_->eval_exc_device_async( N, rho, sigma, eps, stream );  

}


void XCKernel::eval_exc_vxc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps,
  double*       vrho,
  double*       vsigma,
  device::cuda_stream_t* stream
) const {

  pimpl_->eval_exc_vxc_device_async( N, rho, sigma, eps, vrho, vsigma, stream );  

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void XCKernel::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps,
  device::cuda_stream_t* stream
) const {

  pimpl_->eval_exc_device_async( N, rho, sigma, lapl, tau, eps, stream );  

}


void XCKernel::eval_exc_vxc_device_async( 
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
) const {

  pimpl_->eval_exc_vxc_device_async( N, rho, sigma, lapl, tau, 
    eps, vrho, vsigma, vlapl, vtau, stream );  

}


};
