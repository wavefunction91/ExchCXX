#include <exchcxx/xc_functional.hpp>
#include <exchcxx/device/cuda_type_wrappers.hpp>


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

void scal_device_async( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  scal_kernel<<< blocks, threads, 0, stream >>>( N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 1024;
  int blocks  = std::ceil( N / 1024. );
  add_scal_kernel<<< blocks, threads >>>( N, fact, X_device, Y_device );
}

void add_scal_device_async( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
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

void XCFunctional::eval_exc_device( 
  const int     N, 
  const double* rho, 
  double*       eps 
) const {

  throw_if_not_sane();
  assert( is_lda() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;
    kernels_[i].second.eval_exc_device(N, rho, eps_eval);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc
) const {

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
    kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vxc_eval);

    if( i ) {

      add_scal_device( N,       kernels_[i].first, eps_eval, eps );
      add_scal_device( len_vxc, kernels_[i].first, vxc_eval, vxc );

    } else {

      scal_device( N,       kernels_[i].first, eps_eval, eps );
      scal_device( len_vxc, kernels_[i].first, vxc_eval, vxc );

    }
  
  }

  if( eps_scr ) cudaFree( eps_scr );
  if( vxc_scr ) cudaFree( vxc_scr );

}



// GGA Interfaces

void XCFunctional::eval_exc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps 
) const {

  throw_if_not_sane();
  assert( is_gga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc_device(N, rho, eps_eval);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps, 
  double*       vrho,
  double*       vsigma
) const {

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
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval );
    else
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval);

    if( i ) {

      add_scal_device( N, kernels_[i].first, eps_eval, eps );
      add_scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho );
      if( kernels_[i].second.is_gga() )
        add_scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma );

    } else {

      scal_device( N, kernels_[i].first, eps_eval, eps );
      scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho );
      if( kernels_[i].second.is_gga() )
        scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma );

    }
  
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );

}




// mGGA Interfaces

void XCFunctional::eval_exc_device( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps
) const {

  throw_if_not_sane();
  assert( is_mgga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, lapl, tau, eps_eval);
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval);
    else
      kernels_[i].second.eval_exc_device(N, rho, eps_eval);

    if( i ) 
      add_scal_device( N, kernels_[i].first, eps_eval, eps );
    else
      scal_device( N, kernels_[i].first, eps_eval, eps );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device( 
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

  throw_if_not_sane();
  assert( is_gga() );

  int len_vrho   = is_polarized() ? 2*N : N;
  int len_vsigma = is_polarized() ? 3*N : N;
  int len_vlapl  = is_polarized() ? 2*N : N;
  int len_vtau   = is_polarized() ? 2*N : N;

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr), *vlapl_scr(nullptr), *vtau_scr(nullptr);
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
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, lapl, tau, eps_eval, vrho_eval, vsigma_eval, vlapl_eval, vtau_eval );
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval );
    else
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval);

    if( i ) {

      add_scal_device( N, kernels_[i].first, eps_eval, eps );
      add_scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho );

      if( kernels_[i].second.is_gga() )
        add_scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma );

      if( kernels_[i].second.is_mgga() ) {
        add_scal_device( len_vlapl, kernels_[i].first, vlapl_eval, vlapl );
        add_scal_device( len_vtau,  kernels_[i].first, vtau_eval,  vtau  );
      }

    } else {

      scal_device( N, kernels_[i].first, eps_eval, eps );
      scal_device( len_vrho, kernels_[i].first, vrho_eval, vrho );

      if( kernels_[i].second.is_gga() )
        scal_device( len_vsigma, kernels_[i].first, vsigma_eval, vsigma );

      if( kernels_[i].second.is_mgga() ) {
        scal_device( len_vlapl, kernels_[i].first, vlapl_eval, vlapl );
        scal_device( len_vtau,  kernels_[i].first, vtau_eval,  vtau  );
      }

    }
  
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );
  if( vlapl_scr )  cudaFree( vlapl_scr );
  if( vtau_scr )   cudaFree( vtau_scr );
}


















void XCFunctional::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  double*       eps,
  device::cuda_stream_t* stream 
) const {

  throw_if_not_sane();
  assert( is_lda() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );

  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;
    kernels_[i].second.eval_exc_device_async(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
    else
      scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device_async( 
  const int     N, 
  const double* rho, 
  double*       eps, 
  double*       vxc,
  device::cuda_stream_t* stream
) const {

  throw_if_not_sane();
  assert( is_lda() );

  int len_vxc = is_polarized() ? 2*N : N;

  double* eps_scr(nullptr), *vxc_scr(nullptr);
  if( kernels_.size() > 1 ) {
    eps_scr = safe_cuda_malloc( N );
    vxc_scr = safe_cuda_malloc( len_vxc );
  }

  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;
    double* vxc_eval = i ? vxc_scr : vxc;
    kernels_[i].second.eval_exc_vxc_device_async(N, rho, eps_eval, vxc_eval, stream);

    if( i ) {

      add_scal_device_async( N,       kernels_[i].first, eps_eval, eps, st );
      add_scal_device_async( len_vxc, kernels_[i].first, vxc_eval, vxc, st );

    } else {

      scal_device_async( N,       kernels_[i].first, eps_eval, eps, st );
      scal_device_async( len_vxc, kernels_[i].first, vxc_eval, vxc, st );

    }
  
  }

  if( eps_scr ) cudaFree( eps_scr );
  if( vxc_scr ) cudaFree( vxc_scr );

}



// GGA Interfaces

void XCFunctional::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps, 
  device::cuda_stream_t* stream
) const {

  throw_if_not_sane();
  assert( is_gga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device_async(N, rho, sigma, eps_eval, stream);
    else
      kernels_[i].second.eval_exc_device_async(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
    else
      scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  double*       eps, 
  double*       vrho,
  double*       vsigma,
  device::cuda_stream_t* stream
) const {

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

  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr    : eps;
    double* vrho_eval   = i ? vrho_scr   : vrho;
    double* vsigma_eval = i ? vsigma_scr : vsigma;

    if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc_device_async(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval, stream );
    else
      kernels_[i].second.eval_exc_vxc_device_async(N, rho, eps_eval, vrho_eval, stream);

    if( i ) {

      add_scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
      add_scal_device_async( len_vrho, kernels_[i].first, vrho_eval, vrho, st);
      if( kernels_[i].second.is_gga() )
        add_scal_device_async( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, st );

    } else {

      scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
      scal_device_async( len_vrho, kernels_[i].first, vrho_eval, vrho, st );
      if( kernels_[i].second.is_gga() )
        scal_device_async( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, st );

    }
  
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );

}




// mGGA Interfaces

void XCFunctional::eval_exc_device_async( 
  const int     N, 
  const double* rho, 
  const double* sigma, 
  const double* lapl, 
  const double* tau, 
  double*       eps,
  device::cuda_stream_t* stream
) const {

  throw_if_not_sane();
  assert( is_mgga() );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 ) 
    eps_scr = safe_cuda_malloc( N );


  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval = i ? eps_scr : eps;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_device_async(N, rho, sigma, lapl, tau, eps_eval, stream);
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_device_async(N, rho, sigma, eps_eval, stream);
    else
      kernels_[i].second.eval_exc_device_async(N, rho, eps_eval, stream);

    if( i ) 
      add_scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
    else
      scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
  
  }

  if( eps_scr ) cudaFree( eps_scr );

}


void XCFunctional::eval_exc_vxc_device_async( 
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

  throw_if_not_sane();
  assert( is_gga() );

  int len_vrho   = is_polarized() ? 2*N : N;
  int len_vsigma = is_polarized() ? 3*N : N;
  int len_vlapl  = is_polarized() ? 2*N : N;
  int len_vtau   = is_polarized() ? 2*N : N;

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr), *vlapl_scr(nullptr), *vtau_scr(nullptr);
  if( kernels_.size() > 1 ) {
    eps_scr    = safe_cuda_malloc( N );
    vrho_scr   = safe_cuda_malloc( len_vrho );
    vsigma_scr = safe_cuda_malloc( len_vsigma );
    vlapl_scr  = safe_cuda_malloc( len_vlapl );
    vtau_scr   = safe_cuda_malloc( len_vtau );
  }

  cudaStream_t& st = *stream->stream;
  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    double* eps_eval    = i ? eps_scr    : eps;
    double* vrho_eval   = i ? vrho_scr   : vrho;
    double* vsigma_eval = i ? vsigma_scr : vsigma;
    double* vlapl_eval  = i ? vlapl_scr  : vlapl;
    double* vtau_eval   = i ? vtau_scr   : vtau;

    if( kernels_[i].second.is_mgga() )
      kernels_[i].second.eval_exc_vxc_device_async(N, rho, sigma, lapl, tau, eps_eval, vrho_eval, vsigma_eval, vlapl_eval, vtau_eval, stream );
    else if( kernels_[i].second.is_gga() )
      kernels_[i].second.eval_exc_vxc_device_async(N, rho, sigma, eps_eval, vrho_eval, vsigma_eval, stream );
    else
      kernels_[i].second.eval_exc_vxc_device_async(N, rho, eps_eval, vrho_eval, stream);

    if( i ) {

      add_scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
      add_scal_device_async( len_vrho, kernels_[i].first, vrho_eval, vrho, st );

      if( kernels_[i].second.is_gga() )
        add_scal_device_async( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, st );

      if( kernels_[i].second.is_mgga() ) {
        add_scal_device_async( len_vlapl, kernels_[i].first, vlapl_eval, vlapl, st );
        add_scal_device_async( len_vtau,  kernels_[i].first, vtau_eval,  vtau, st  );
      }

    } else {

      scal_device_async( N, kernels_[i].first, eps_eval, eps, st );
      scal_device_async( len_vrho, kernels_[i].first, vrho_eval, vrho, st );

      if( kernels_[i].second.is_gga() )
        scal_device_async( len_vsigma, kernels_[i].first, vsigma_eval, vsigma, st );

      if( kernels_[i].second.is_mgga() ) {
        scal_device_async( len_vlapl, kernels_[i].first, vlapl_eval, vlapl, st );
        scal_device_async( len_vtau,  kernels_[i].first, vtau_eval,  vtau, st  );
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
