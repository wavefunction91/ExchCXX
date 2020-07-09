#include "libxc_common.hpp"
#include <exchcxx/impl/builtin/util.hpp>
//#include <functionals.cuh>

void throw_if_fail( cudaError_t stat, std::string msg ) {
  if( stat != cudaSuccess ) throw std::runtime_error( msg );
}

void recv_from_device( void* dest, const void* src, const size_t len ) {

  auto stat = cudaMemcpy( dest, src, len, cudaMemcpyDeviceToHost );
  throw_if_fail( stat, "recv failed" );

}

void recv_from_device( void* dest, const void* src, const size_t len, 
  cudaStream_t& stream ) {

  auto stat = cudaMemcpyAsync( dest, src, len, cudaMemcpyDeviceToHost, stream );
  throw_if_fail( stat, "recv failed" );

}

void send_to_device( void* dest, const void* src, const size_t len ) {

  auto stat = cudaMemcpy( dest, src, len, cudaMemcpyHostToDevice);
  throw_if_fail( stat, "send failed" );

}

void send_to_device( void* dest, const void* src, const size_t len, 
  cudaStream_t& stream ) {

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
LDA_EXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_device_ ) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N );

  recv_from_device( rho_host.data(), rho, len_rho, stream );

  stream_sync( stream );
  xc_lda_exc( &kernel_, N, rho_host.data(), eps_host.data() );

  send_to_device( eps, eps_host.data(), len_eps, stream );
  stream_sync( stream ); // Lifetime of host vectors

}


LDA_EXC_VXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_device_ ) const {

  throw_if_uninitialized();
  assert( is_lda() );

  size_t len_rho = N*sizeof(double);
  size_t len_eps = N*sizeof(double);
  size_t len_vxc = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), vxc_host( N );

  recv_from_device( rho_host.data(), rho, len_rho, stream );

  stream_sync( stream );
  xc_lda_exc_vxc( &kernel_, N, rho_host.data(), eps_host.data(), vxc_host.data() );

  send_to_device( eps, eps_host.data(), len_eps, stream );
  send_to_device( vxc, vxc_host.data(), len_vxc, stream );
  stream_sync( stream ); // Lifetime of host vectors

}

// TODO: LDA kxc interfaces

// GGA interface
GGA_EXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_device_ ) const {

  throw_if_uninitialized();
  assert( is_gga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N );

  recv_from_device( rho_host.data(),   rho,   len_rho  , stream );
  recv_from_device( sigma_host.data(), sigma, len_sigma, stream );

  stream_sync( stream );
  xc_gga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data() );

  send_to_device( eps, eps_host.data(), len_eps, stream );
  stream_sync( stream ); // Lifetime of host vectors

}


GGA_EXC_VXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_device_ ) const {

  throw_if_uninitialized();
  assert( is_gga() );


  size_t len_rho    = N*sizeof(double);
  size_t len_sigma  = N*sizeof(double);
  size_t len_vrho   = N*sizeof(double);
  size_t len_vsigma = N*sizeof(double);
  size_t len_eps    = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), 
    vrho_host( N ), vsigma_host( N );

  recv_from_device( rho_host.data(),   rho,   len_rho  , stream );
  recv_from_device( sigma_host.data(), sigma, len_sigma, stream );
 
  stream_sync( stream );
  xc_gga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), eps_host.data(), 
    vrho_host.data(), vsigma_host.data() );

  send_to_device( eps,    eps_host.data(),    len_eps   , stream);
  send_to_device( vrho,   vrho_host.data(),   len_vrho  , stream);
  send_to_device( vsigma, vsigma_host.data(), len_vsigma, stream);
  stream_sync( stream ); // Lifetime of host vectors

}

// TODO: GGA kxc interfaces  
  
  
// mGGA interface
MGGA_EXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_device_ ) const {

  throw_if_uninitialized();
  assert( is_mgga() );

  size_t len_rho   = N*sizeof(double);
  size_t len_sigma = N*sizeof(double);
  size_t len_lapl  = N*sizeof(double);
  size_t len_tau   = N*sizeof(double);
  size_t len_eps   = N*sizeof(double);

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), 
    lapl_host( N ), tau_host( N );

  recv_from_device( rho_host.data(),   rho,   len_rho  , stream );
  recv_from_device( sigma_host.data(), sigma, len_sigma, stream );
  recv_from_device( lapl_host.data(),  lapl,  len_lapl , stream );
  recv_from_device( tau_host.data(),   tau,   len_tau  , stream );

  stream_sync( stream );
  xc_mgga_exc( &kernel_, N, rho_host.data(), sigma_host.data(), lapl_host.data(), 
    tau_host.data(), eps_host.data() );

  send_to_device( eps, eps_host.data(), len_eps, stream );
  stream_sync( stream ); // Lifetime of host vectors

}


MGGA_EXC_VXC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_device_ ) const {

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

  std::vector<double> rho_host( N ), eps_host( N ), sigma_host( N ), 
    lapl_host( N ), tau_host( N );
  std::vector<double> vrho_host( N ), vsigma_host( N ),  vlapl_host( N ), 
    vtau_host( N );

  recv_from_device( rho_host.data(),   rho,   len_rho  , stream );
  recv_from_device( sigma_host.data(), sigma, len_sigma, stream );
  recv_from_device( lapl_host.data(),  lapl,  len_lapl , stream );
  recv_from_device( tau_host.data(),   tau,   len_tau  , stream );

  stream_sync( stream );
  xc_mgga_exc_vxc( &kernel_, N, rho_host.data(), sigma_host.data(), 
    lapl_host.data(), tau_host.data(), eps_host.data(), vrho_host.data(), 
    vsigma_host.data(), vlapl_host.data(), vtau_host.data() );

  send_to_device( eps,    eps_host.data(), len_eps      , stream );
  send_to_device( vrho,   vrho_host.data(),   len_vrho  , stream );
  send_to_device( vsigma, vsigma_host.data(), len_vsigma, stream );
  send_to_device( vlapl,  vlapl_host.data(),  len_vlapl , stream );
  send_to_device( vtau,   vtau_host.data(),   len_vtau  , stream );
  stream_sync( stream ); // Lifetime of host vectors

}







LDA_EXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_inc_device_ ) const {
  disabled_inc_interface();
}

LDA_EXC_VXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_inc_device_ ) const {
  disabled_inc_interface();
}

GGA_EXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_inc_device_ ) const {
  disabled_inc_interface();
}

GGA_EXC_VXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_inc_device_ ) const {
  disabled_inc_interface();
}

MGGA_EXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_inc_device_ ) const {
  disabled_inc_interface();
}

MGGA_EXC_VXC_INC_GENERATOR_DEVICE( LibxcKernelImpl::eval_exc_vxc_inc_device_ ) const {
  disabled_inc_interface();
}

}
}
