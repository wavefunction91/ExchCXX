/**
 * ExchCXX 
 *
 * Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). 
 *
 * Portions Copyright (c) Microsoft Corporation.
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>
#include <exchcxx/util/div_ceil.hpp>

namespace ExchCXX {
namespace detail {

template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_unpolar( rho_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

  if( tid < N ) {

    const double rho_use = fmax( rho[tid], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, eps[tid], vxc[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

  if( tid < N ) {

    const double rho_use   = fmax( rho[tid],   0.    );
    const double sigma_use = fmax( sigma[tid], 1e-40 );
    traits::eval_exc_unpolar( rho_use, sigma_use, eps[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
  int tid = hipThreadIdx_x + hipBlockIdx_x * hipBlockDim_x; 

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
__global__ MGGA_EXC_GENERATOR( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;
    traits::eval_exc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid], eps[tid] );

  }

}


template <typename KernelType>
__global__ MGGA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* lapl_i  = traits::needs_laplacian ? (lapl + 2*tid) : nullptr;
    auto* tau_i   = tau   + 2*tid;

    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0],
      tau_i[1], eps[tid] );

  }

}

template <typename KernelType>
__global__ MGGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;

    double dummy;
    auto& vlapl_return = traits::needs_laplacian ? vlapl[tid] : dummy; 
    traits::eval_exc_vxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      eps[tid], vrho[tid], vsigma[tid], vlapl_return, vtau[tid] );

  }

}

template <typename KernelType>
__global__ MGGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double dummy_vlapl[2];

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* lapl_i  = traits::needs_laplacian ? (lapl + 2*tid) : lapl;
    auto* tau_i   = tau   + 2*tid;

    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;
    auto* vlapl_i  = traits::needs_laplacian ? vlapl + 2*tid : dummy_vlapl;
    auto* vtau_i   = vtau   + 2*tid;
    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
                                                         
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0],
      tau_i[1], eps[tid], vrho_i[0], vrho_i[1], vsigma_i[0], vsigma_i[1], 
      vsigma_i[2], vlapl_i[0], vlapl_i[1], vtau_i[0], vtau_i[1] );

  }

}


template <typename KernelType>
__global__ MGGA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e;
  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;
    traits::eval_exc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid], e );
    eps[tid] += scal_fact * e;
     

  }

}

template <typename KernelType>
__global__ MGGA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* lapl_i  = traits::needs_laplacian ? (lapl + 2*tid) : lapl;
    auto* tau_i   = tau   + 2*tid;

    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0], 
      tau_i[1], e );
    eps[tid] += scal_fact * e;
     

  }

}

template <typename KernelType>
__global__ MGGA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e, vr, vs, vl, vt;
  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;

    traits::eval_exc_vxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      e, vr, vs, vl, vt );
    eps[tid]    += scal_fact * e;
    vrho[tid]   += scal_fact * vr;
    vsigma[tid] += scal_fact * vs;
    vtau[tid]   += scal_fact * vt;
    if(traits::needs_laplacian) vlapl[tid] += scal_fact * vl;

  }

}

template <typename KernelType>
__global__ MGGA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double dummy_vlapl[2];
  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* lapl_i  = traits::needs_laplacian ? (lapl + 2*tid) : lapl;
    auto* tau_i   = tau   + 2*tid;

    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;
    auto* vlapl_i  = traits::needs_laplacian ? vlapl + 2*tid : dummy_vlapl;
    auto* vtau_i   = vtau   + 2*tid;

    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
                                                         
                                                         
    double e, vra, vrb, vsaa,vsab,vsbb, vla, vlb, vta, vtb;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0],
      tau_i[1], e, vra, vrb, vsaa, vsab, vsbb, vla, vlb, vta, vtb );

    eps[tid]    += scal_fact * e;
    vrho_i[0]   += scal_fact * vra;
    vrho_i[1]   += scal_fact * vrb;
    vsigma_i[0] += scal_fact * vsaa;
    vsigma_i[1] += scal_fact * vsab;
    vsigma_i[2] += scal_fact * vsbb;
    vtau_i[0]   += scal_fact * vta;
    vtau_i[1]   += scal_fact * vtb;
    if(traits::needs_laplacian) {
      vlapl_i[0]   += scal_fact * vla;
      vlapl_i[1]   += scal_fact * vlb;
    }

  }

}

template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(
    HIP_KERNEL_NAME(device_eval_exc_helper_unpolar_kernel<KernelType>),
    blocks, threads, 0, stream,
    N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_inc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_inc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, eps
  );

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, eps, vxc
  );

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, eps, vxc
  );

}




template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, sigma, eps, vrho, vsigma
  );

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    N, rho, sigma, eps, vrho, vsigma
  );

}


template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_inc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_inc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, sigma, eps
  );

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, sigma, eps, vrho, vsigma
  );

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  hipLaunchKernelGGL(HIP_KERNEL_NAME(device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>), blocks, threads, 0, stream, 
    scal_fact, N, rho, sigma, eps, vrho, vsigma
  );

}

template <typename KernelType>
MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, eps
  );

}

template <typename KernelType>
MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, eps
  );

}

template <typename KernelType>
MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau
  );

}

template <typename KernelType>
MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau
  );

}


template <typename KernelType>
MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau, eps
  );

}

template <typename KernelType>
MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_exc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau, eps
  );

}

template <typename KernelType>
MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau
  );

}

template <typename KernelType>
MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );

  device_eval_exc_vxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau, eps, vrho, vsigma, vlapl, vtau
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

#define MGGA_GENERATE_DEVICE_HELPERS(KERN) \
  template MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); 

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


MGGA_GENERATE_DEVICE_HELPERS( BuiltinSCAN_X );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinSCAN_C );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinR2SCAN_X );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinR2SCAN_C );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinFT98_X );

MGGA_GENERATE_DEVICE_HELPERS( BuiltinPC07_K );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinPC07OPT_K );

MGGA_GENERATE_DEVICE_HELPERS( BuiltinSCANL_C );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinSCANL_X );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinR2SCANL_C );
MGGA_GENERATE_DEVICE_HELPERS( BuiltinR2SCANL_X );

LDA_GENERATE_DEVICE_HELPERS( BuiltinEPC17_1 )
LDA_GENERATE_DEVICE_HELPERS( BuiltinEPC17_2 )
LDA_GENERATE_DEVICE_HELPERS( BuiltinEPC18_1 )
LDA_GENERATE_DEVICE_HELPERS( BuiltinEPC18_2 )

}
}

