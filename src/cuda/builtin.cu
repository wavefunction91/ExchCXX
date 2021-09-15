#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>
#include <exchcxx/util/div_ceil.hpp>

namespace ExchCXX {
namespace detail {


template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_unpolar( rho_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );

    traits::eval_exc_polar( rho_a_use, rho_b_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, eps[tid], vxc[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    auto vxc_i = vxc + 2*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );

    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, eps[tid], 
      vxc_i[0], vxc_i[1] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e;
  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_unpolar( rho_use, e );
    eps[tid] += scal_fact * e;

  }

}

template <typename KernelType>
__global__ LDA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );

    double e;
    traits::eval_exc_polar( rho_a_use, rho_b_use, e );
    
    eps[tid] += scal_fact * e;

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e,v;
  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, e, v );
    eps[tid] += scal_fact * e;
    vxc[tid] += scal_fact * v;

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    auto vxc_i = vxc + 2*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );

    double v_a, v_b, e;
    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, e, v_a, v_b);
    eps[tid] += scal_fact * e;
    vxc_i[0] += scal_fact * v_a;
    vxc_i[1] += scal_fact * v_b;

  }

}

template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double rho_use   = fmax( rho[tid],   0.    );
    const double sigma_use = fmax( sigma[tid], 1e-40 );
    traits::eval_exc_unpolar( rho_use, sigma_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );
    const double sigma_aa_use = fmax( sigma_i[0], 1e-40 );
    const double sigma_bb_use = fmax( sigma_i[2], 1e-40 );
    const double sigma_ab_use = fmax( 
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );

    traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use, 
      sigma_ab_use, sigma_bb_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double rho_use   = fmax( rho[tid],   0.    );
    const double sigma_use = fmax( sigma[tid], 1e-40 );
    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, eps[tid], 
      vrho[tid], vsigma[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i    = rho   + 2*tid;
    auto* sigma_i  = sigma + 3*tid;
    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );
    const double sigma_aa_use = fmax( sigma_i[0], 1e-40 );
    const double sigma_bb_use = fmax( sigma_i[2], 1e-40 );
    const double sigma_ab_use = fmax( 
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );
                                                         
                                                         
    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, sigma_aa_use, 
      sigma_ab_use, sigma_bb_use, eps[tid], vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2] );

  }

}


template <typename KernelType>
__global__ GGA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e;
  if( tid < N ) {

    const double rho_use   = fmax( rho[tid],   0.    );
    const double sigma_use = fmax( sigma[tid], 1e-40 );
                                      
    traits::eval_exc_unpolar( rho_use, sigma_use, e );
    eps[tid] += scal_fact * e;
     

  }

}

template <typename KernelType>
__global__ GGA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );
    const double sigma_aa_use = fmax( sigma_i[0], 1e-40 );
    const double sigma_bb_use = fmax( sigma_i[2], 1e-40 );
    const double sigma_ab_use = fmax( 
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );

    double e;
    traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use, 
      sigma_ab_use, sigma_bb_use, e );
    eps[tid] += scal_fact * e;
     

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e, vr, vs;
  if( tid < N ) {

    const double rho_use   = fmax( rho[tid],   0.    );
    const double sigma_use = fmax( sigma[tid], 1e-40 );

    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, e, vr, vs );
    eps[tid]    += scal_fact * e;
    vrho[tid]   += scal_fact * vr;
    vsigma[tid] += scal_fact * vs;

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i    = rho   + 2*tid;
    auto* sigma_i  = sigma + 3*tid;
    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;

    const double rho_a_use = fmax( rho_i[0], 0. );
    const double rho_b_use = fmax( rho_i[1], 0. );
    const double sigma_aa_use = fmax( sigma_i[0], 1e-40 );
    const double sigma_bb_use = fmax( sigma_i[2], 1e-40 );
    const double sigma_ab_use = fmax( 
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );
                                                         
                                                         
    double e, vra, vrb, vsaa,vsab,vsbb;
    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, sigma_aa_use, 
      sigma_ab_use, sigma_bb_use, e, vra, vrb, vsaa, vsab, vsbb );

    eps[tid]    += scal_fact * e;
    vrho_i[0]   += scal_fact * vra;
    vrho_i[1]   += scal_fact * vrb;
    vsigma_i[0] += scal_fact * vsaa;
    vsigma_i[1] += scal_fact * vsab;
    vsigma_i[2] += scal_fact * vsbb;

  }

}


template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_vxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_vxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_vxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, eps, vxc
  );

}




template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, eps, vrho, vsigma
  );

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, eps, vrho, vsigma
  );

}


template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, eps, vrho, vsigma
  );

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, eps, vrho, vsigma
  );

}

#define LDA_GENERATE_DEVICE_HELPERS(KERN) \
  template LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar<KERN> ); \
  template LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar<KERN> ); \
  template LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); 

#define GGA_GENERATE_DEVICE_HELPERS(KERN) \
  template GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar<KERN> ); \
  template GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar<KERN> ); \
  template GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); 

LDA_GENERATE_DEVICE_HELPERS( BuiltinSlaterExchange );
LDA_GENERATE_DEVICE_HELPERS( BuiltinVWN3 );
LDA_GENERATE_DEVICE_HELPERS( BuiltinVWN_RPA );
LDA_GENERATE_DEVICE_HELPERS( BuiltinPW91_LDA );
LDA_GENERATE_DEVICE_HELPERS( BuiltinPW91_LDA_MOD );
LDA_GENERATE_DEVICE_HELPERS( BuiltinPW91_LDA_RPA );
LDA_GENERATE_DEVICE_HELPERS( BuiltinPZ81 );
LDA_GENERATE_DEVICE_HELPERS( BuiltinPZ81_MOD );

GGA_GENERATE_DEVICE_HELPERS( BuiltinB88   );
GGA_GENERATE_DEVICE_HELPERS( BuiltinLYP   );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE_X );
GGA_GENERATE_DEVICE_HELPERS( BuiltinRevPBE_X );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE_C );

GGA_GENERATE_DEVICE_HELPERS( BuiltinB3LYP );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE0  );



}
}

