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
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_exc_unpolar( rho[tid], eps[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    traits::eval_exc_polar( rho_i[0], rho_i[1], eps[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_exc_vxc_unpolar( rho[tid], eps[tid], vxc[tid] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    auto vxc_i = vxc + 2*tid;

    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], eps[tid], 
      vxc_i[0], vxc_i[1] );

  }

}


template <typename KernelType>
__global__ LDA_FXC_GENERATOR( device_eval_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_fxc_unpolar( rho[tid],  fxc[tid] );

  }

}

template <typename KernelType>
__global__ LDA_FXC_GENERATOR( device_eval_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    auto v2rho2_i = fxc + 3*tid;

    traits::eval_fxc_polar( rho_i[0], rho_i[1], v2rho2_i[0], 
      v2rho2_i[1], v2rho2_i[2] );
  }
}

template <typename KernelType>
__global__ LDA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_vxc_fxc_unpolar( rho[tid], vxc[tid], fxc[tid] );

  }

}

template <typename KernelType>
__global__ LDA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;
    auto vxc_i = vxc + 2*tid;
    auto v2rho2_i = fxc + 3*tid;

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], vxc_i[0], vxc_i[1], 
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2] );

  }

}

template <typename KernelType>
__global__ LDA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e;
  if( tid < N ) {

    traits::eval_exc_unpolar( rho[tid], e );
    eps[tid] += scal_fact * e;

  }

}

template <typename KernelType>
__global__ LDA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto rho_i = rho + 2*tid;

    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], e );
    
    eps[tid] += scal_fact * e;

  }

}

template <typename KernelType>
__global__ LDA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e,v;
  if( tid < N ) {

    traits::eval_exc_vxc_unpolar( rho[tid], e, v );
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

    double v_a, v_b, e;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], e, v_a, v_b);
    eps[tid] += scal_fact * e;
    vxc_i[0] += scal_fact * v_a;
    vxc_i[1] += scal_fact * v_b;

  }

}

template <typename KernelType>
__global__ LDA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    double f;
    traits::eval_fxc_unpolar( rho[tid], f );
    fxc[tid] += scal_fact * f;
  }
}

template <typename KernelType>
__global__ LDA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_polar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    auto rho_i = rho + 2*tid;
    auto fxc_i = fxc + 3*tid;
    double f0, f1, f2;
    traits::eval_fxc_polar( rho_i[0], rho_i[1], f0, f1, f2 );
    fxc_i[0] += scal_fact * f0;
    fxc_i[1] += scal_fact * f1;
    fxc_i[2] += scal_fact * f2;
  }
}

template <typename KernelType>
__global__ LDA_VXC_FXC_INC_GENERATOR( device_eval_vxc_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    double v, f;
    traits::eval_vxc_fxc_unpolar( rho[tid], v, f );
    vxc[tid] += scal_fact * v;
    fxc[tid] += scal_fact * f;
  }
}

template <typename KernelType>
__global__ LDA_VXC_FXC_INC_GENERATOR( device_eval_vxc_fxc_inc_helper_polar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    auto rho_i = rho + 2*tid;
    auto vxc_i = vxc + 2*tid;
    auto fxc_i = fxc + 3*tid;
    double v0, v1, f0, f1, f2;
    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], v0, v1, f0, f1, f2 );
    vxc_i[0] += scal_fact * v0;
    vxc_i[1] += scal_fact * v1;
    fxc_i[0] += scal_fact * f0;
    fxc_i[1] += scal_fact * f1;
    fxc_i[2] += scal_fact * f2;
  }
}






template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_exc_unpolar( rho[tid], sigma[tid], eps[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_GENERATOR( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;

    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], eps[tid] );

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_GENERATOR( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_exc_vxc_unpolar( rho[tid], sigma[tid], eps[tid], 
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

    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], eps[tid], vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2] );

  }

}

template <typename KernelType>
__global__ GGA_FXC_GENERATOR( device_eval_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_fxc_unpolar( rho[tid], sigma[tid], v2rho2[tid], v2rhosigma[tid], v2sigma2[tid] );

  }

}

template <typename KernelType>
__global__ GGA_FXC_GENERATOR( device_eval_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* v2rho2_i = v2rho2 + 3*tid;
    auto* v2rhosigma_i = v2rhosigma + 6*tid;
    auto* v2sigma2_i = v2sigma2 + 6*tid;


    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2], 
      v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
      v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
      v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
      v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5] );

  }
}

template <typename KernelType>
__global__ GGA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    traits::eval_vxc_fxc_unpolar( rho[tid], sigma[tid], vrho[tid], vsigma[tid],
      v2rho2[tid], v2rhosigma[tid], v2sigma2[tid] );

  }

}

template <typename KernelType>
__global__ GGA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;
    auto* v2rho2_i = v2rho2 + 3*tid;
    auto* v2rhosigma_i = v2rhosigma + 6*tid;
    auto* v2sigma2_i = v2sigma2 + 6*tid;

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      vrho_i[0], vrho_i[1], vsigma_i[0], vsigma_i[1], vsigma_i[2],
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2], 
      v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
      v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
      v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
      v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5] );

  }

}


template <typename KernelType>
__global__ GGA_EXC_INC_GENERATOR( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e;
  if( tid < N ) {
                                      
    traits::eval_exc_unpolar( rho[tid], sigma[tid], e );
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
    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], e );
    eps[tid] += scal_fact * e;
     

  }

}

template <typename KernelType>
__global__ GGA_EXC_VXC_INC_GENERATOR( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  double e, vr, vs;
  if( tid < N ) {

    traits::eval_exc_vxc_unpolar( rho[tid], sigma[tid], e, vr, vs );
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
                                                         
    double e, vra, vrb, vsaa,vsab,vsbb;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], e, vra, vrb, vsaa, vsab, vsbb );

    eps[tid]    += scal_fact * e;
    vrho_i[0]   += scal_fact * vra;
    vrho_i[1]   += scal_fact * vrb;
    vsigma_i[0] += scal_fact * vsaa;
    vsigma_i[1] += scal_fact * vsab;
    vsigma_i[2] += scal_fact * vsbb;

  }

}


template <typename KernelType>
__global__ GGA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    double f2, f3, f4;
    traits::eval_fxc_unpolar( rho[tid], sigma[tid], f2, f3, f4 );
    v2rho2[tid]    += scal_fact * f2;
    v2rhosigma[tid] += scal_fact * f3;
    v2sigma2[tid]  += scal_fact * f4;
  }
}

template <typename KernelType>
__global__ GGA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_polar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* v2rho2_i = v2rho2 + 3*tid;
    auto* v2rhosigma_i = v2rhosigma + 6*tid;
    auto* v2sigma2_i = v2sigma2 + 6*tid;
    double f2[3], f3[6], f4[6];
    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      f2[0], f2[1], f2[2],
      f3[0], f3[1], f3[2], f3[3], f3[4], f3[5],
      f4[0], f4[1], f4[2], f4[3], f4[4], f4[5] );
    for(int i=0;i<3;++i) v2rho2_i[i] += scal_fact * f2[i];
    for(int i=0;i<6;++i) v2rhosigma_i[i] += scal_fact * f3[i];
    for(int i=0;i<6;++i) v2sigma2_i[i] += scal_fact * f4[i];
  }
}

template <typename KernelType>
__global__ GGA_VXC_FXC_INC_GENERATOR( device_eval_vxc_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    double v, s, f2, f3, f4;
    traits::eval_vxc_fxc_unpolar( rho[tid], sigma[tid], v, s, f2, f3, f4 );
    vrho[tid]   += scal_fact * v;
    vsigma[tid] += scal_fact * s;
    v2rho2[tid]    += scal_fact * f2;
    v2rhosigma[tid] += scal_fact * f3;
    v2sigma2[tid]  += scal_fact * f4;
  }
}

template <typename KernelType>
__global__ GGA_VXC_FXC_INC_GENERATOR( device_eval_vxc_fxc_inc_helper_polar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if( tid < N ) {
    auto* rho_i   = rho   + 2*tid;
    auto* sigma_i = sigma + 3*tid;
    auto* vrho_i   = vrho   + 2*tid;
    auto* vsigma_i = vsigma + 3*tid;
    auto* v2rho2_i = v2rho2 + 3*tid;
    auto* v2rhosigma_i = v2rhosigma + 6*tid;
    auto* v2sigma2_i = v2sigma2 + 6*tid;
    double v[2], s[3], f2[3], f3[6], f4[6];
    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      v[0], v[1], s[0], s[1], s[2],
      f2[0], f2[1], f2[2],
      f3[0], f3[1], f3[2], f3[3], f3[4], f3[5],
      f4[0], f4[1], f4[2], f4[3], f4[4], f4[5] );
    for(int i=0;i<2;++i) vrho_i[i] += scal_fact * v[i];
    for(int i=0;i<3;++i) vsigma_i[i] += scal_fact * s[i];
    for(int i=0;i<3;++i) v2rho2_i[i] += scal_fact * f2[i];
    for(int i=0;i<6;++i) v2rhosigma_i[i] += scal_fact * f3[i];
    for(int i=0;i<6;++i) v2sigma2_i[i] += scal_fact * f4[i];
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
__global__ MGGA_FXC_GENERATOR( device_eval_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;
    double local_v2rholapl, local_v2sigmalapl, local_v2lapl2, local_v2lapltau;

    auto& v2rholapl_return = traits::needs_laplacian ?  v2rholapl[tid] : local_v2rholapl;
    auto& v2sigmalapl_return = traits::needs_laplacian ?  v2sigmalapl[tid] : local_v2sigmalapl;
    auto& v2lapl2_return = traits::needs_laplacian ?  v2lapl2[tid] : local_v2lapl2;
    auto& v2lapltau_return = traits::needs_laplacian ?  v2lapltau[tid] : local_v2lapltau;

    traits::eval_fxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      v2rho2[tid], v2rhosigma[tid], v2rholapl_return, v2rhotau[tid],
      v2sigma2[tid], v2sigmalapl_return, v2sigmatau[tid],
      v2lapl2_return, v2lapltau_return, v2tau2[tid] );
  }

}

template <typename KernelType>
__global__ MGGA_FXC_GENERATOR( device_eval_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 
  double dummy_v2rholapl[4];
  double dummy_v2sigmalapl[6];
  double dummy_v2lapl2[3];
  double dummy_v2lapltau[4];

  if( tid < N ) {

    auto* rho_i           = rho           + 2 * tid;
    auto* sigma_i         = sigma         + 3 * tid;
    auto* tau_i           = tau           + 2 * tid;
    auto* v2rho2_i        = v2rho2        + 3 * tid;
    auto* v2rhosigma_i    = v2rhosigma    + 6 * tid;
    auto* v2rhotau_i      = v2rhotau      + 4 * tid;
    auto* v2sigma2_i      = v2sigma2      + 6 * tid;
    auto* v2sigmatau_i    = v2sigmatau    + 6 * tid;
    auto* v2tau2_i        = v2tau2        + 3 * tid;

    auto* lapl_i          = traits::needs_laplacian ? (lapl + 2 * tid) : lapl;
    auto* v2rholapl_i     = traits::needs_laplacian ? (v2rholapl + 4 * tid) : dummy_v2rholapl;
    auto* v2sigmalapl_i   = traits::needs_laplacian ? (v2sigmalapl + 6 * tid) : dummy_v2sigmalapl;
    auto* v2lapl2_i       = traits::needs_laplacian ? (v2lapl2 + 3 * tid) : dummy_v2lapl2;
    auto* v2lapltau_i     = traits::needs_laplacian ? (v2lapltau + 4 * tid) : dummy_v2lapltau;

    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2],
      v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
      v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
      v2rholapl_i[0], v2rholapl_i[1], v2rholapl_i[2], v2rholapl_i[3],
      v2rhotau_i[0], v2rhotau_i[1], v2rhotau_i[2], v2rhotau_i[3],
      v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
      v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5],
      v2sigmalapl_i[0], v2sigmalapl_i[1], v2sigmalapl_i[2],
      v2sigmalapl_i[3], v2sigmalapl_i[4], v2sigmalapl_i[5],
      v2sigmatau_i[0], v2sigmatau_i[1], v2sigmatau_i[2],
      v2sigmatau_i[3], v2sigmatau_i[4], v2sigmatau_i[5],
      v2lapl2_i[0], v2lapl2_i[1], v2lapl2_i[2],
      v2lapltau_i[0], v2lapltau_i[1], v2lapltau_i[2], v2lapltau_i[3],
      v2tau2_i[0], v2tau2_i[1], v2tau2_i[2] );

    
  }
}

template <typename KernelType>
__global__ MGGA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 

  if( tid < N ) {

    const double lapl_use  = traits::needs_laplacian ? lapl[tid] : 0.0;
    double dummy_v2rholapl, dummy_v2sigmalapl, dummy_v2lapl2, dummy_v2lapltau, dummy_vlapl;
    auto& vlapl_return = traits::needs_laplacian ? vlapl[tid] : dummy_vlapl;
    auto& v2rholapl_return = traits::needs_laplacian ?  v2rholapl[tid] : dummy_v2rholapl;
    auto& v2sigmalapl_return = traits::needs_laplacian ?  v2sigmalapl[tid] : dummy_v2sigmalapl;
    auto& v2lapl2_return = traits::needs_laplacian ?  v2lapl2[tid] : dummy_v2lapl2;
    auto& v2lapltau_return = traits::needs_laplacian ?  v2lapltau[tid] : dummy_v2lapltau;

    traits::eval_vxc_fxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      vrho[tid], vsigma[tid], vlapl_return, vtau[tid],
      v2rho2[tid], v2rhosigma[tid], v2rholapl_return, 
      v2rhotau[tid], v2sigma2[tid], v2sigmalapl_return,
      v2sigmatau[tid], v2lapl2_return, v2lapltau_return,
      v2tau2[tid] );

  }

}

template <typename KernelType>
__global__ MGGA_VXC_FXC_GENERATOR( device_eval_vxc_fxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 
  double dummy_vlapl[2];
  double dummy_v2rholapl[4];
  double dummy_v2sigmalapl[6];
  double dummy_v2lapl2[3];
  double dummy_v2lapltau[4];

  if( tid < N ) {

    auto* rho_i           = rho           + 2 * tid;
    auto* sigma_i         = sigma         + 3 * tid;
    auto* tau_i           = tau           + 2 * tid;
    auto* vrho_i          = vrho          + 2 * tid;
    auto* vsigma_i        = vsigma        + 3 * tid;
    auto* vtau_i          = vtau          + 2 * tid;

    auto* v2rho2_i        = v2rho2        + 3 * tid;
    auto* v2rhosigma_i    = v2rhosigma    + 6 * tid;
    auto* v2rhotau_i      = v2rhotau      + 4 * tid;
    auto* v2sigma2_i      = v2sigma2      + 6 * tid;
    auto* v2sigmatau_i    = v2sigmatau    + 6 * tid;
    auto* v2tau2_i        = v2tau2        + 3 * tid;

    auto* lapl_i          = traits::needs_laplacian ? (lapl + 2 * tid) : lapl;
    auto* vlapl_i         = traits::needs_laplacian ? (vlapl + 2 * tid) : dummy_vlapl;
    auto* v2rholapl_i     = traits::needs_laplacian ? (v2rholapl + 4 * tid) : dummy_v2rholapl;
    auto* v2sigmalapl_i   = traits::needs_laplacian ? (v2sigmalapl + 6 * tid) : dummy_v2sigmalapl;
    auto* v2lapl2_i       = traits::needs_laplacian ? (v2lapl2 + 3 * tid) : dummy_v2lapl2;
    auto* v2lapltau_i     = traits::needs_laplacian ? (v2lapltau + 4 * tid) : dummy_v2lapltau;
    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
      vrho_i[0], vrho_i[1], vsigma_i[0], vsigma_i[1], vsigma_i[2],
      vlapl_i[0], vlapl_i[1], vtau_i[0], vtau_i[1],
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2],
      v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
      v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
      v2rholapl_i[0], v2rholapl_i[1], v2rholapl_i[2], v2rholapl_i[3],
      v2rhotau_i[0], v2rhotau_i[1], v2rhotau_i[2], v2rhotau_i[3],
      v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
      v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5],
      v2sigmalapl_i[0], v2sigmalapl_i[1], v2sigmalapl_i[2],
      v2sigmalapl_i[3], v2sigmalapl_i[4], v2sigmalapl_i[5],
      v2sigmatau_i[0], v2sigmatau_i[1], v2sigmatau_i[2],
      v2sigmatau_i[3], v2sigmatau_i[4], v2sigmatau_i[5],
      v2lapl2_i[0], v2lapl2_i[1], v2lapl2_i[2],
      v2lapltau_i[0], v2lapltau_i[1], v2lapltau_i[2], v2lapltau_i[3],
      v2tau2_i[0], v2tau2_i[1], v2tau2_i[2] );

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
__global__ MGGA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 
  if( tid < N ) {
    const double lapl_use = traits::needs_laplacian ? lapl[tid] : 0.0;
    double f_rho2, f_rhosigma, f_rholapl, f_rhotau, f_sigma2, f_sigmalapl, f_sigmatau, f_lapl2, f_lapltau, f_tau2;
    traits::eval_fxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      f_rho2, f_rhosigma, f_rholapl, f_rhotau,
      f_sigma2, f_sigmalapl, f_sigmatau,
      f_lapl2, f_lapltau, f_tau2 );
    v2rho2[tid]     += scal_fact * f_rho2;
    v2rhosigma[tid] += scal_fact * f_rhosigma;
    v2rhotau[tid]   += scal_fact * f_rhotau;
    v2sigma2[tid]   += scal_fact * f_sigma2;
    v2sigmatau[tid] += scal_fact * f_sigmatau;
    v2tau2[tid]     += scal_fact * f_tau2;
    if(traits::needs_laplacian) {
      v2rholapl[tid]   += scal_fact * f_rholapl;
      v2sigmalapl[tid] += scal_fact * f_sigmalapl;
      v2lapl2[tid]     += scal_fact * f_lapl2;
      v2lapltau[tid]   += scal_fact * f_lapltau;
    }
  }
}

template <typename KernelType>
__global__ MGGA_FXC_INC_GENERATOR( device_eval_fxc_inc_helper_polar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if( tid < N ) {
    auto* rho_i           = rho           + 2 * tid;
    auto* sigma_i         = sigma         + 3 * tid;
    auto* tau_i           = tau           + 2 * tid;
    auto* v2rho2_i        = v2rho2        + 3 * tid;
    auto* v2rhosigma_i    = v2rhosigma    + 6 * tid;
    auto* v2rhotau_i      = v2rhotau      + 4 * tid;
    auto* v2sigma2_i      = v2sigma2      + 6 * tid;
    auto* v2sigmatau_i    = v2sigmatau    + 6 * tid;
    auto* v2tau2_i        = v2tau2        + 3 * tid;

    auto* lapl_i          = traits::needs_laplacian ? (lapl + 2 * tid) : lapl;
    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    double f_rho2[3], f_rhosigma[6], f_rholapl[4], f_rhotau[4], f_sigma2[6], f_sigmalapl[6], f_sigmatau[6], f_lapl2[3], f_lapltau[4], f_tau2[3];
    
    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
      f_rho2[0], f_rho2[1], f_rho2[2],
      f_rhosigma[0], f_rhosigma[1], f_rhosigma[2], f_rhosigma[3], f_rhosigma[4], f_rhosigma[5],
      f_rholapl[0], f_rholapl[1], f_rholapl[2], f_rholapl[3],
      f_rhotau[0], f_rhotau[1], f_rhotau[2], f_rhotau[3],
      f_sigma2[0], f_sigma2[1], f_sigma2[2], f_sigma2[3], f_sigma2[4], f_sigma2[5],
      f_sigmalapl[0], f_sigmalapl[1], f_sigmalapl[2], f_sigmalapl[3], f_sigmalapl[4], f_sigmalapl[5],
      f_sigmatau[0], f_sigmatau[1], f_sigmatau[2], f_sigmatau[3], f_sigmatau[4], f_sigmatau[5],
      f_lapl2[0], f_lapl2[1], f_lapl2[2],
      f_lapltau[0], f_lapltau[1], f_lapltau[2], f_lapltau[3],
      f_tau2[0], f_tau2[1], f_tau2[2] );
    
    for(int i=0;i<3;++i) v2rho2_i[i] += scal_fact * f_rho2[i];
    for(int i=0;i<6;++i) v2rhosigma_i[i] += scal_fact * f_rhosigma[i];
    for(int i=0;i<4;++i) v2rhotau_i[i] += scal_fact * f_rhotau[i];
    for(int i=0;i<6;++i) v2sigma2_i[i] += scal_fact * f_sigma2[i];
    for(int i=0;i<6;++i) v2sigmatau_i[i] += scal_fact * f_sigmatau[i];
    for(int i=0;i<3;++i) v2tau2_i[i] += scal_fact * f_tau2[i];
    
    if(traits::needs_laplacian) {
      auto* v2rholapl_i     = v2rholapl + 4 * tid;
      auto* v2sigmalapl_i   = v2sigmalapl + 6 * tid;
      auto* v2lapl2_i       = v2lapl2 + 3 * tid;
      auto* v2lapltau_i     = v2lapltau + 4 * tid;
      for(int i=0;i<4;++i) v2rholapl_i[i] += scal_fact * f_rholapl[i];
      for(int i=0;i<6;++i) v2sigmalapl_i[i] += scal_fact * f_sigmalapl[i];
      for(int i=0;i<3;++i) v2lapl2_i[i] += scal_fact * f_lapl2[i];
      for(int i=0;i<4;++i) v2lapltau_i[i] += scal_fact * f_lapltau[i];
    }
  }
}


template <typename KernelType>
__global__ MGGA_VXC_FXC_INC_GENERATOR( device_eval_vxc_fxc_inc_helper_unpolar_kernel ) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x; 
  if( tid < N ) {
    const double lapl_use = traits::needs_laplacian ? lapl[tid] : 0.0;
    double f_rho2, f_rhosigma, f_rholapl, f_rhotau, f_sigma2, f_sigmalapl, f_sigmatau, f_lapl2, f_lapltau, f_tau2;
    double vr, vs, vl, vt;
    traits::eval_vxc_fxc_unpolar( rho[tid], sigma[tid], lapl_use, tau[tid],
      vr, vs, vl, vt,
      f_rho2, f_rhosigma, f_rholapl, f_rhotau,
      f_sigma2, f_sigmalapl, f_sigmatau,
      f_lapl2, f_lapltau, f_tau2);
      
    vrho[tid]   += scal_fact * vr;
    vsigma[tid] += scal_fact * vs;
    vtau[tid]   += scal_fact * vt;
    v2rho2[tid]     += scal_fact * f_rho2;
    v2rhosigma[tid] += scal_fact * f_rhosigma;
    v2rhotau[tid]   += scal_fact * f_rhotau;
    v2sigma2[tid]   += scal_fact * f_sigma2;
    v2sigmatau[tid] += scal_fact * f_sigmatau;
    v2tau2[tid]     += scal_fact * f_tau2;

    if(traits::needs_laplacian) {
      vlapl[tid] += scal_fact * vl;
      v2rholapl[tid]   += scal_fact * f_rholapl;
      v2sigmalapl[tid] += scal_fact * f_sigmalapl;
      v2lapl2[tid]     += scal_fact * f_lapl2;
      v2lapltau[tid]   += scal_fact * f_lapltau;
    }
  }
}

template <typename KernelType>
__global__ MGGA_VXC_FXC_INC_GENERATOR(device_eval_vxc_fxc_inc_helper_polar_kernel) {
  using traits = kernel_traits<KernelType>;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if( tid < N ) {
    auto* rho_i           = rho           + 2 * tid;
    auto* sigma_i         = sigma         + 3 * tid;
    auto* tau_i           = tau           + 2 * tid;
    auto* vrho_i          = vrho          + 2 * tid;
    auto* vsigma_i        = vsigma        + 3 * tid;
    auto* vtau_i          = vtau          + 2 * tid;

    auto* v2rho2_i        = v2rho2        + 3 * tid;
    auto* v2rhosigma_i    = v2rhosigma    + 6 * tid;
    auto* v2rhotau_i      = v2rhotau      + 4 * tid;
    auto* v2sigma2_i      = v2sigma2      + 6 * tid;
    auto* v2sigmatau_i    = v2sigmatau    + 6 * tid;
    auto* v2tau2_i        = v2tau2        + 3 * tid;

    auto* lapl_i          = traits::needs_laplacian ? (lapl + 2 * tid) : lapl;
    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    double frho[2], fsigma[3], flapl[2], ftau[2];
    double f_rho2[3], f_rhosigma[6], f_rholapl[4], f_rhotau[4], f_sigma2[6], f_sigmalapl[6], f_sigmatau[6], f_lapl2[3], f_lapltau[4], f_tau2[3];

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
      frho[0], frho[1], fsigma[0], fsigma[1], fsigma[2],
      flapl[0], flapl[1], ftau[0], ftau[1],
      f_rho2[0], f_rho2[1], f_rho2[2],
      f_rhosigma[0], f_rhosigma[1], f_rhosigma[2],
      f_rhosigma[3], f_rhosigma[4], f_rhosigma[5],
      f_rholapl[0], f_rholapl[1], f_rholapl[2], f_rholapl[3],
      f_rhotau[0], f_rhotau[1], f_rhotau[2], f_rhotau[3],
      f_sigma2[0], f_sigma2[1], f_sigma2[2],
      f_sigma2[3], f_sigma2[4], f_sigma2[5],
      f_sigmalapl[0], f_sigmalapl[1], f_sigmalapl[2],
      f_sigmalapl[3], f_sigmalapl[4], f_sigmalapl[5],
      f_sigmatau[0], f_sigmatau[1], f_sigmatau[2],
      f_sigmatau[3], f_sigmatau[4], f_sigmatau[5],
      f_lapl2[0], f_lapl2[1], f_lapl2[2],
      f_lapltau[0], f_lapltau[1], f_lapltau[2], f_lapltau[3],
      f_tau2[0], f_tau2[1], f_tau2[2] );

    for(int i=0;i<2;++i) vrho_i[i] += scal_fact * frho[i];
    for(int i=0;i<3;++i) vsigma_i[i] += scal_fact * fsigma[i];
    for(int i=0;i<2;++i) vtau_i[i] += scal_fact * ftau[i];

    for(int i=0;i<3;++i) v2rho2_i[i] += scal_fact * f_rho2[i];
    for(int i=0;i<6;++i) v2rhosigma_i[i] += scal_fact * f_rhosigma[i];
    for(int i=0;i<4;++i) v2rhotau_i[i] += scal_fact * f_rhotau[i];
    for(int i=0;i<6;++i) v2sigma2_i[i] += scal_fact * f_sigma2[i];
    for(int i=0;i<6;++i) v2sigmatau_i[i] += scal_fact * f_sigmatau[i];
    for(int i=0;i<3;++i) v2tau2_i[i] += scal_fact * f_tau2[i];

    if(traits::needs_laplacian) {
      auto* vlapl_i         = vlapl + 2 * tid;
      auto* v2rholapl_i     = v2rholapl + 4 * tid;
      auto* v2sigmalapl_i   = v2sigmalapl + 6 * tid;
      auto* v2lapl2_i       = v2lapl2 + 3 * tid;
      auto* v2lapltau_i     = v2lapltau + 4 * tid;
      for(int i=0;i<2;++i) vlapl_i[i]   += scal_fact * flapl[i];
      for(int i=0;i<4;++i) v2rholapl_i[i] += scal_fact * f_rholapl[i];
      for(int i=0;i<6;++i) v2sigmalapl_i[i] += scal_fact * f_sigmalapl[i];
      for(int i=0;i<3;++i) v2lapl2_i[i] += scal_fact * f_lapl2[i];
      for(int i=0;i<4;++i) v2lapltau_i[i] += scal_fact * f_lapltau[i];
    }
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
LDA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, fxc
  );

}

template <typename KernelType>
LDA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, fxc
  );
}

template <typename KernelType>
LDA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, vxc, fxc
  );

}

template <typename KernelType>
LDA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, vxc, fxc
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
LDA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, fxc
  );

}

template <typename KernelType>
LDA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, fxc
  );

}

template <typename KernelType>
LDA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, vxc, fxc
  );

}

template <typename KernelType>
LDA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, vxc, fxc
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
GGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, v2rho2, v2rhosigma, v2sigma2
  );

}

template <typename KernelType>
GGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, v2rho2, v2rhosigma, v2sigma2
  );

}

template <typename KernelType>
GGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2
  );

}

template <typename KernelType>
GGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2
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


template <typename KernelType>
GGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2 );
}

template <typename KernelType>
GGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2 );
}

template <typename KernelType>
GGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);
}

template <typename KernelType>
GGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2);

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
MGGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    N, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
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

template <typename KernelType>
MGGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_unpolar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau,
    vrho, vsigma, vlapl, vtau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
  );

}

template <typename KernelType>
MGGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar ) {

  dim3 threads(32);
  dim3 blocks( util::div_ceil( N, threads.x) );
  device_eval_vxc_fxc_inc_helper_polar_kernel<KernelType><<<blocks,threads,0,stream>>>(
    scal_fact, N, rho, sigma, lapl, tau,
    vrho, vsigma, vlapl, vtau,
    v2rho2, v2rhosigma, v2rholapl, v2rhotau,
    v2sigma2, v2sigmalapl, v2sigmatau,
    v2lapl2, v2lapltau, v2tau2
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
  template LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template LDA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar<KERN> ); \
  template LDA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar<KERN> ); \
  template LDA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template LDA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar<KERN> ); \
  template LDA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar<KERN> ); \
  template LDA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar<KERN> ); \
  template LDA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template LDA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar<KERN> );

#define GGA_GENERATE_DEVICE_HELPERS(KERN) \
  template GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar<KERN> ); \
  template GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar<KERN> ); \
  template GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template GGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar<KERN> ); \
  template GGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar<KERN> ); \
  template GGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template GGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar<KERN> ); \
  template GGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar<KERN> ); \
  template GGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar<KERN> ); \
  template GGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template GGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar<KERN> );

#define MGGA_GENERATE_DEVICE_HELPERS(KERN) \
  template MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template MGGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template MGGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_unpolar<KERN> ); \
  template MGGA_FXC_GENERATOR_DEVICE( device_eval_fxc_helper_polar<KERN> ); \
  template MGGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template MGGA_VXC_FXC_GENERATOR_DEVICE( device_eval_vxc_fxc_helper_polar<KERN> ); \
  template MGGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_unpolar<KERN> ); \
  template MGGA_FXC_INC_GENERATOR_DEVICE( device_eval_fxc_inc_helper_polar<KERN> ); \
  template MGGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template MGGA_VXC_FXC_INC_GENERATOR_DEVICE( device_eval_vxc_fxc_inc_helper_polar<KERN> );

}
}

