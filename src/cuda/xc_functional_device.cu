#include <exchcxx/xc_functional.hpp>


__global__ void scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if( tid < N ) Y_device[tid] = X_device[tid] * fact;

}

__global__ void add_scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if( tid < N ) Y_device[tid] += X_device[tid] * fact;

}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  scal_kernel<<< blocks, threads >>>( N, fact, X_device, Y_device );
}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  scal_kernel<<< blocks, threads, 0, stream >>>( N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  add_scal_kernel<<< blocks, threads >>>( N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  add_scal_kernel<<< blocks, threads, 0, stream >>>( N, fact, X_device, Y_device );
}


template <typename T = double>
T* safe_cuda_malloc( size_t N ) {

  T* ptr = nullptr;
  auto stat = cudaMalloc( &ptr, N*sizeof(T) );
  if( stat != cudaSuccess ) throw std::runtime_error("Alloc Failed");

  return ptr;

}

namespace ExchCXX {






LDA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  assert( is_lda() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;
    kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps, stream );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps, stream );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


LDA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  assert( is_lda() );

  int len_vxc = is_polarized() ? 2*N : N;

  double* eps_scr(nullptr), *vxc_scr(nullptr);
  if( kernels_.size() > 1 ) {
    eps_scr = safe_cuda_malloc( N );
    vxc_scr = safe_cuda_malloc( len_vxc );
  }

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;
    double* vxc_eval = i ? vxc_scr : vxc;
    kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vxc_eval, stream);

    if( i ) {

      add_scal_device( N,       kernels_[i].first, eps_eval, eps, stream );
      add_scal_device( len_vxc, kernels_[i].first, vxc_eval, vxc, stream );

    } else {

      scal_device( N,       kernels_[i].first, eps_eval, eps, stream );
      scal_device( len_vxc, kernels_[i].first, vxc_eval, vxc, stream );

    }
  
  }

  if( eps_scr ) cudaFree( eps_scr );
  if( vxc_scr ) cudaFree( vxc_scr );

}



// GGA Interfaces

GGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  assert( is_gga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, stream);
    else
      kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps, stream );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps, stream );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


GGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  assert( is_gga() );

  int len_vrho   = is_polarized() ? 2*N : N;
  int len_vsigma = is_polarized() ? 3*N : N;

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr);
  if( kernels_.size() > 1 ) {
    eps_scr    = safe_cuda_malloc( N );
    vrho_scr   = safe_cuda_malloc( len_vrho );
    vsigma_scr = safe_cuda_malloc( len_vsigma );
  }

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr    : eps;
    double* vrho_eval   = i ? vrho_scr   : vrho;
    double* vsigma_eval = i ? vsigma_scr : vsigma;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval, 
        vsigma_eval, stream );
    else
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval, stream);

    if( i ) {

      add_scal_device( N, kernels_[i].first, eps_eval, eps, stream );
      add_scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho, stream);
      if( kernels_[i].second.is_gga() )
        add_scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, stream );

    } else {

      scal_device( N, kernels_[i].first, eps_eval, eps, stream );
      scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho, stream );
      if( kernels_[i].second.is_gga() )
        scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, stream );

    }
  
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );

}




// mGGA Interfaces

MGGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  assert( is_mgga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, lapl, tau, eps_eval, stream);
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, stream);
    else
      kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps, stream );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps, stream );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


MGGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  assert( is_gga() );

  int len_vrho   = is_polarized() ? 2*N : N;
  int len_vsigma = is_polarized() ? 3*N : N;
  int len_vlapl  = is_polarized() ? 2*N : N;
  int len_vtau   = is_polarized() ? 2*N : N;

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr), 
    *vlapl_scr(nullptr), *vtau_scr(nullptr);
  if( kernels_.size() > 1 ) {
    eps_scr    = safe_cuda_malloc( N );
    vrho_scr   = safe_cuda_malloc( len_vrho );
    vsigma_scr = safe_cuda_malloc( len_vsigma );
    vlapl_scr  = safe_cuda_malloc( len_vlapl );
    vtau_scr   = safe_cuda_malloc( len_vtau );
  }

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr    : eps;
    double* vrho_eval   = i ? vrho_scr   : vrho;
    double* vsigma_eval = i ? vsigma_scr : vsigma;
    double* vlapl_eval  = i ? vlapl_scr  : vlapl;
    double* vtau_eval   = i ? vtau_scr   : vtau;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, lapl, tau, eps_eval, 
        vrho_eval, vsigma_eval, vlapl_eval, vtau_eval, stream );
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval, 
        vsigma_eval, stream );
    else
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval, stream);

    if( i ) {

      add_scal_device( N, kernels_[i].first, eps_eval, eps, stream );
      add_scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho, stream );

      if( kernels_[i].second.is_gga() )
        add_scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, stream );

      if( kernels_[i].second.is_mgga() ) {
        add_scal_device( len_vlapl, kernels_[i].first, vlapl_eval, vlapl, stream );
        add_scal_device( len_vtau,  kernels_[i].first, vtau_eval,  vtau, stream  );
      }

    } else {

      scal_device( N, kernels_[i].first, eps_eval, eps, stream );
      scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho, stream );

      if( kernels_[i].second.is_gga() )
        scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, stream );

      if( kernels_[i].second.is_mgga() ) {
        scal_device( len_vlapl, kernels_[i].first, vlapl_eval, vlapl, stream );
        scal_device( len_vtau,  kernels_[i].first, vtau_eval,  vtau, stream  );
      }

    }
  
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );
  if( vlapl_scr )  cudaFree( vlapl_scr );
  if( vtau_scr )   cudaFree( vtau_scr );
}

}
