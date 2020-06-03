#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>

namespace ExchCXX {
namespace detail {

template <typename Integral1, typename Integral2>
int64_t div_ceil(Integral1 x, Integral2 y) {
  int64_t x_ll = x;
  int64_t y_ll = y;

  auto d = std::div(x_ll, y_ll);
  return d.quot + !!d.rem;
}


template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N )
    traits::eval_exc_unpolar( rho[tid], eps[tid] );

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N )
    traits::eval_exc_vxc_unpolar( rho[tid], eps[tid], vxc[tid] );

}

template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N )
    traits::eval_exc_unpolar( rho[tid], sigma[tid], eps[tid] );

}

template <typename KernelType>
__global__ GGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N )
    traits::eval_exc_vxc_unpolar( rho[tid], sigma[tid], eps[tid], 
      vrho[tid], vsigma[tid] );

}




template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper ) {

  dim3 threads(std::min( 1024, N ));
  dim3 blocks( div_ceil( N, threads.x) );
  device_eval_exc_helper_kernel<KernelType><<<blocks,threads>>>(
    N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper ) {

  dim3 threads(std::min( 1024, N ));
  dim3 blocks( div_ceil( N, threads.x) );
  device_eval_exc_vxc_helper_kernel<KernelType><<<blocks,threads>>>(
    N, rho, eps, vxc
  );

}




template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper ) {

  dim3 threads(std::min( 1024, N ));
  dim3 blocks( div_ceil( N, threads.x) );
  device_eval_exc_helper_kernel<KernelType><<<blocks,threads>>>(
    N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper ) {

  dim3 threads(std::min( 1024, N ));
  dim3 blocks( div_ceil( N, threads.x) );
  device_eval_exc_vxc_helper_kernel<KernelType><<<blocks,threads>>>(
    N, rho, sigma, eps, vrho, vsigma
  );

}

#define LDA_GENERATE_DEVICE_HELPERS(KERN) \
  template LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper<KERN> ); \
  template LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper<KERN> ); \

#define GGA_GENERATE_DEVICE_HELPERS(KERN) \
  template GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper<KERN> ); \
  template GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper<KERN> ); \

LDA_GENERATE_DEVICE_HELPERS( BuiltinSlaterExchange );

GGA_GENERATE_DEVICE_HELPERS( BuiltinLYP   );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE_X );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE_C );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE0  );



}
}

