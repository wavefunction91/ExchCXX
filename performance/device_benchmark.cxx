#include <chrono>
#include <iostream>
#include <cuda_runtime.h>

#include <exchcxx/xc_kernel.hpp>


template <typename T>
T* safe_cuda_malloc( size_t n ) {

  T* ptr;
  auto stat = cudaMalloc( (void**)&ptr, n*sizeof(T) );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
  return ptr;

}

template <typename T>
void safe_cuda_cpy( T* dest, const T* src, size_t len ) {

  auto stat = cudaMemcpy( dest, src, len*sizeof(T), cudaMemcpyDefault );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

}

void cuda_free_all(){ }
template <typename T, typename... Args>
void cuda_free_all( T* ptr, Args&&... args ) {

  auto stat = cudaFree( ptr );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  cuda_free_all( std::forward<Args>(args)... );


}

void device_synchronize() {
  auto stat = cudaDeviceSynchronize();
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
}

void stream_synchronize(cudaStream_t stream) {
  auto stat = cudaStreamSynchronize( stream );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
}


template <typename Op>
double time_cpu_op( const Op& op ) {

  auto st = std::chrono::high_resolution_clock::now();

  op();

  auto en = std::chrono::high_resolution_clock::now();

  return std::chrono::duration<double,std::milli>(en - st).count();

}

template <typename Op>
double time_gpu_op( const Op& op, cudaStream_t stream ) {

  cudaEvent_t st, en;

  auto stat = cudaEventCreate( &st );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
  stat = cudaEventCreate( &en );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  stream_synchronize( stream );

  stat = cudaEventRecord( st, stream );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  op();

  stream_synchronize( stream );
  stat = cudaEventRecord( en, stream );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  stat = cudaEventSynchronize( en );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  float dur;
  stat = cudaEventElapsedTime( &dur, st, en );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  stat = cudaEventDestroy( st );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));
  stat = cudaEventDestroy( en );
  if( stat != cudaSuccess )
    throw std::runtime_error(cudaGetErrorString( stat ));

  return dur;
  
}

template <typename T>
auto max_diff_vec( const std::vector<T>& a, const std::vector<T>& b ) {

  assert( a.size() == b.size() );

  std::vector<T> diff = a;
  for( auto i = 0; i < a.size(); ++i )
    diff[i] = std::abs( a[i] - b[i] );

  return *std::max_element(diff.begin(), diff.end() );

}

int main() {


  using namespace ExchCXX;
  auto kern  = XCKernel::Kernel::PBE0;
  auto polar = XCKernel::Spin::Unpolarized;

  auto builtin_backend = XCKernel::Backend::builtin;

  XCKernel pbe_libxc( kern, polar );
  XCKernel pbe_built( builtin_backend, kern, polar );

  int npts = 100000;
  std::vector<double> rho( npts, 0.1 );
  std::vector<double> sigma( npts, 0.2 );

  std::vector<double> eps(npts), vrho(npts), vsigma(npts);


  double* rho_device    = safe_cuda_malloc<double>( npts );
  double* sigma_device  = safe_cuda_malloc<double>( npts );
  double* eps_device    = safe_cuda_malloc<double>( npts );
  double* vrho_device   = safe_cuda_malloc<double>( npts );
  double* vsigma_device = safe_cuda_malloc<double>( npts );

  safe_cuda_cpy( rho_device, rho.data(), npts );
  safe_cuda_cpy( sigma_device, sigma.data(), npts );
  device_synchronize();

  cudaStream_t stream = 0;

  int nrep = 5;

  // Warmup CPU
  for( int i = 0; i < nrep; ++i ) {
    pbe_libxc.eval_exc_vxc( npts, rho.data(), sigma.data(), eps.data(), 
      vrho.data(), vsigma.data() );
  }

  auto eps_ref = eps;
  auto vrho_ref = vrho;
  auto vsigma_ref = vsigma;

  // Warmup GPU
  for( int i = 0; i < nrep; ++i ) {
    pbe_built.eval_exc_vxc_device( npts, rho_device, sigma_device, eps_device, 
      vrho_device, vsigma_device, stream );
  }


  // Get timings for CPU
  auto libxc_host_dur = time_cpu_op([&]() {
    for( int i = 0; i < nrep; ++i ) {
      pbe_libxc.eval_exc_vxc( npts, rho.data(), sigma.data(), eps.data(), 
        vrho.data(), vsigma.data() );
    }
  });

  auto builtin_host_dur = time_cpu_op([&]() {
    for( int i = 0; i < nrep; ++i ) {
      pbe_built.eval_exc_vxc( npts, rho.data(), sigma.data(), eps.data(), 
        vrho.data(), vsigma.data() );
    }
  });

  // Get timings for GPU
  auto libxc_device_dur = time_gpu_op([&]() {
    for( int i = 0; i < nrep; ++i ) {
      pbe_libxc.eval_exc_vxc_device( npts, rho_device, sigma_device, eps_device, 
        vrho_device, vsigma_device, stream );
    }
  }, stream);

  auto builtin_device_dur = time_gpu_op([&]() {
    for( int i = 0; i < nrep; ++i ) {
      pbe_built.eval_exc_vxc_device( npts, rho_device, sigma_device, eps_device, 
        vrho_device, vsigma_device, stream );
    }
  }, stream);


  safe_cuda_cpy( eps.data(), eps_device, npts );
  safe_cuda_cpy( vrho.data(), vrho_device, npts );
  safe_cuda_cpy( vsigma.data(), vsigma_device, npts );

  device_synchronize();

  std::cout << max_diff_vec( eps, eps_ref ) << std::endl;


  
  std::cout << "Libxc Host   = " << libxc_host_dur << std::endl;
  std::cout << "Built Host   = " << builtin_host_dur << std::endl;
  std::cout << "Libxc Device = " << libxc_device_dur << std::endl;
  std::cout << "Built Device = " << builtin_device_dur << std::endl;



  cuda_free_all( rho_device, sigma_device, eps_device, vrho_device, vsigma_device );
}
