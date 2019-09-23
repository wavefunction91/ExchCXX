#include <exchcxx/impl/libxc.hpp>
#include <exchcxx/factory/xc_kernel.hpp>
#include <exchcxx/device/cuda_type_wrappers.hpp>

void throw_if_fail( cudaError_t stat, std::string msg ) {
  if( stat != cudaSuccess ) throw std::runtime_error( msg );
}

void recv_from_device( void* dest, const void* src, const size_t len ) {

  auto stat = cudaMemcpy( dest, src, len, cudaMemcpyDeviceToHost );
  throw_if_fail( stat, "recv failed" );

}

void recv_from_device_async( void* dest, const void* src, const size_t len, cudaStream_t& stream ) {

  auto stat = cudaMemcpyAsync( dest, src, len, cudaMemcpyDeviceToHost, stream );
  throw_if_fail( stat, "recv failed" );

}

void send_to_device( void* dest, const void* src, const size_t len ) {

  auto stat = cudaMemcpy( dest, src, len, cudaMemcpyHostToDevice);
  throw_if_fail( stat, "send failed" );

}

void send_to_device_async( void* dest, const void* src, const size_t len, cudaStream_t& stream ) {

  auto stat = cudaMemcpyAsync( dest, src, len, cudaMemcpyHostToDevice, stream);
  throw_if_fail( stat, "send failed" );

}

void stream_sync( cudaStream_t& stream ) {

  auto stat = cudaStreamSynchronize( stream );
  throw_if_fail( stat, "sync failed" );

}

namespace ExchCXX {
namespace detail {

// LDA interfaces
void LibxcKernelImpl::eval_exc_device_( 
  const int     N, 
  const double* rho_device, 
  double*       eps_device 
) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N );

  recv_from_device( rho_host.data(), rho_device, len_rho );

  xc_lda_exc( &kernel_, N, rho_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}


void LibxcKernelImpl::eval_exc_vxc_device_( 
  const int     N, 
  const double* rho_device, 
  double*       eps_device, 
  double*       vxc_device 
) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);
  size_t len_vxc = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), vxc_host( N );

  recv_from_device( rho_host.data(), rho_device, len_rho );

  xc_lda_exc_vxc( &kernel_, N, rho_host.data(), eps_host.data(), vxc_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );
  send_to_device( vxc_device, vxc_host.data(), len_vxc );

}

// TODO: LDA kxc interfaces

// GGA interface
void LibxcKernelImpl::eval_exc_device_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  double*       eps_device
) const {

  throw_if_uninitialized();
  assert( is_gga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );

  xc_gga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}

/*
__global__ void xc_gga_exc_vxc_device(
  const xc_func_type* func_device,
  const int           np,
  const double*       rho_device,
  const double*       sigma_device,
        double*       eps_device,
        double*       vrho_device,
        double*       vsigma_device
) {

  dvc_xc_gga_exc_vxc( func_device, np, rho_device, sigma_device, eps_device, vrho_device, vsigma_device );

}
*/

void LibxcKernelImpl::eval_exc_vxc_device_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  double*       eps_device,
  double*       vrho_device,
  double*       vsigma_device
) const {

  throw_if_uninitialized();
  assert( is_gga() );

#if 1

  size_t len_rho    = N*sizeof(double);
  size_t len_sigma  = N*sizeof(double);
  size_t len_vrho   = N*sizeof(double);
  size_t len_vsigma = N*sizeof(double);
  size_t len_eps    = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), vrho_host( N ), vsigma_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );

  xc_gga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data(), vrho_host.data(), vsigma_host.data() );

  send_to_device( eps_device,    eps_host.data(),    len_eps    );
  send_to_device( vrho_device,   vrho_host.data(),   len_vrho   );
  send_to_device( vsigma_device, vsigma_host.data(), len_vsigma );

#else

  //dim3 threads = 1024;
  //dim3 blocks  = std::ceil( N / 1024. );
  //xc_gga_exc_vxc_device<<< blocks, threads >>>( &kernel_, N, rho_device, 
  //  sigma_device, eps_device, vrho_device, vsigma_device );

  cudaError_t stat;
  stat = cudaMemset( eps_device,    0, N*sizeof(double) ); throw_if_fail( stat, "EPS    ZERO" );
  stat = cudaMemset( vrho_device,   0, N*sizeof(double) ); throw_if_fail( stat, "VRHO   ZERO" );
  stat = cudaMemset( vsigma_device, 0, N*sizeof(double) ); throw_if_fail( stat, "VSIGMA ZERO" );

#endif

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void LibxcKernelImpl::eval_exc_device_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  const double* lapl_device, 
  const double* tau_device, 
  double*       eps_device
) const {

  throw_if_uninitialized();
  assert( is_mgga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);
  size_t len_lapl  = N*sizeof(double);
  size_t len_tau   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), lapl_host( N ), tau_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );
  recv_from_device( lapl_host.data(),  lapl_device,  len_lapl  );
  recv_from_device( tau_host.data(),   tau_device,   len_tau   );

  xc_mgga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), lapl_host.data(), tau_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}


void LibxcKernelImpl::eval_exc_vxc_device_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  const double* lapl_device, 
  const double* tau_device, 
  double*       eps_device,
  double*       vrho_device,
  double*       vsigma_device,
  double*       vlapl_device, 
  double*       vtau_device
) const {

  throw_if_uninitialized();
  assert( is_mgga() );

  size_t len_rho    = N*sizeof(double);
  size_t len_sigma  = N*sizeof(double);
  size_t len_lapl   = N*sizeof(double);
  size_t len_tau    = N*sizeof(double);
  size_t len_eps    = N*sizeof(double);
  size_t len_vrho   = N*sizeof(double);
  size_t len_vsigma = N*sizeof(double);
  size_t len_vlapl  = N*sizeof(double);
  size_t len_vtau   = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), lapl_host( N ), tau_host( N );
  std::vector<double> vrho_host( N ), vsigma_host( N ),  vlapl_host( N ), vtau_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );
  recv_from_device( lapl_host.data(),  lapl_device,  len_lapl  );
  recv_from_device( tau_host.data(),   tau_device,   len_tau   );

  xc_mgga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), lapl_host.data(), tau_host.data(), eps_host.data(), vrho_host.data(), vsigma_host.data(), vlapl_host.data(), vtau_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );
  send_to_device( vrho_device,   vrho_host.data(),   len_vrho   );
  send_to_device( vsigma_device, vsigma_host.data(), len_vsigma );
  send_to_device( vlapl_device,  vlapl_host.data(),  len_vlapl  );
  send_to_device( vtau_device,   vtau_host.data(),   len_vtau   );
}






// LDA interfaces
void LibxcKernelImpl::eval_exc_device_async_( 
  const int     N, 
  const double* rho_device, 
  double*       eps_device,
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N );

  recv_from_device( rho_host.data(), rho_device, len_rho );

  xc_lda_exc( &kernel_, N, rho_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}


void LibxcKernelImpl::eval_exc_vxc_device_async_( 
  const int     N, 
  const double* rho_device, 
  double*       eps_device, 
  double*       vxc_device, 
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);
  size_t len_vxc = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), vxc_host( N );

  recv_from_device( rho_host.data(), rho_device, len_rho );

  xc_lda_exc_vxc( &kernel_, N, rho_host.data(), eps_host.data(), vxc_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );
  send_to_device( vxc_device, vxc_host.data(), len_vxc );

}

// TODO: LDA kxc interfaces

// GGA interface
void LibxcKernelImpl::eval_exc_device_async_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  double*       eps_device,
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_gga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );

  xc_gga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}


void LibxcKernelImpl::eval_exc_vxc_device_async_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  double*       eps_device,
  double*       vrho_device,
  double*       vsigma_device,
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_gga() );

#if 1

  size_t len_rho    = N*sizeof(double);
  size_t len_sigma  = N*sizeof(double);
  size_t len_vrho   = N*sizeof(double);
  size_t len_vsigma = N*sizeof(double);
  size_t len_eps    = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), vrho_host( N ), vsigma_host( N );

  cudaStream_t& st = *stream->stream;
  recv_from_device_async( rho_host.data(),   rho_device,   len_rho  , st );
  recv_from_device_async( sigma_host.data(), sigma_device, len_sigma, st );
 
  stream_sync( st );
  xc_gga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data(), vrho_host.data(), vsigma_host.data() );

  send_to_device( eps_device,    eps_host.data(),    len_eps    );
  send_to_device( vrho_device,   vrho_host.data(),   len_vrho   );
  send_to_device( vsigma_device, vsigma_host.data(), len_vsigma );

#else

  //dim3 threads = 1024;
  //dim3 blocks  = std::ceil( N / 1024. );
  //xc_gga_exc_vxc_device<<< blocks, threads >>>( &kernel_, N, rho_device, 
  //  sigma_device, eps_device, vrho_device, vsigma_device );

  cudaError_t stat;
  stat = cudaMemsetAsync( eps_device,    0, N*sizeof(double), *stream->stream ); throw_if_fail( stat, "EPS    ZERO" );
  stat = cudaMemsetAsync( vrho_device,   0, N*sizeof(double), *stream->stream ); throw_if_fail( stat, "VRHO   ZERO" );
  stat = cudaMemsetAsync( vsigma_device, 0, N*sizeof(double), *stream->stream ); throw_if_fail( stat, "VSIGMA ZERO" );

#endif

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
void LibxcKernelImpl::eval_exc_device_async_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  const double* lapl_device, 
  const double* tau_device, 
  double*       eps_device,
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_mgga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);
  size_t len_lapl  = N*sizeof(double);
  size_t len_tau   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), lapl_host( N ), tau_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );
  recv_from_device( lapl_host.data(),  lapl_device,  len_lapl  );
  recv_from_device( tau_host.data(),   tau_device,   len_tau   );

  xc_mgga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), lapl_host.data(), tau_host.data(), eps_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );

}


void LibxcKernelImpl::eval_exc_vxc_device_async_( 
  const int     N, 
  const double* rho_device, 
  const double* sigma_device, 
  const double* lapl_device, 
  const double* tau_device, 
  double*       eps_device,
  double*       vrho_device,
  double*       vsigma_device,
  double*       vlapl_device, 
  double*       vtau_device,
  device::cuda_stream_t* stream
) const {

  throw_if_uninitialized();
  assert( is_mgga() );

  size_t len_rho    = N*sizeof(double);
  size_t len_sigma  = N*sizeof(double);
  size_t len_lapl   = N*sizeof(double);
  size_t len_tau    = N*sizeof(double);
  size_t len_eps    = N*sizeof(double);
  size_t len_vrho   = N*sizeof(double);
  size_t len_vsigma = N*sizeof(double);
  size_t len_vlapl  = N*sizeof(double);
  size_t len_vtau   = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), lapl_host( N ), tau_host( N );
  std::vector<double> vrho_host( N ), vsigma_host( N ),  vlapl_host( N ), vtau_host( N );

  recv_from_device( rho_host.data(),   rho_device,   len_rho   );
  recv_from_device( sigma_host.data(), sigma_device, len_sigma );
  recv_from_device( lapl_host.data(),  lapl_device,  len_lapl  );
  recv_from_device( tau_host.data(),   tau_device,   len_tau   );

  xc_mgga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), lapl_host.data(), tau_host.data(), eps_host.data(), vrho_host.data(), vsigma_host.data(), vlapl_host.data(), vtau_host.data() );

  send_to_device( eps_device, eps_host.data(), len_eps );
  send_to_device( vrho_device,   vrho_host.data(),   len_vrho   );
  send_to_device( vsigma_device, vsigma_host.data(), len_vsigma );
  send_to_device( vlapl_device,  vlapl_host.data(),  len_vlapl  );
  send_to_device( vtau_device,   vtau_host.data(),   len_vtau   );
}

}
}
