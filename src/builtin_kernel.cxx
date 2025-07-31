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
#include <exchcxx/util/unused.hpp>


namespace ExchCXX {


// LDA interface
UNUSED_INTERFACE_GENERATOR( LDA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( LDA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( LDA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( LDA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )

// GGA interface
UNUSED_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( GGA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( GGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )

// MGGA interface
UNUSED_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )


// INC interfaces
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )

UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )

UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )


#ifdef EXCHCXX_ENABLE_DEVICE

// LDA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )

// GGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )

// MGGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )


// INC interfaces
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

#endif



namespace detail {

template <typename KernelType>
LDA_EXC_GENERATOR( host_eval_exc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_exc_unpolar( rho[i], eps[i] );

  }

}


template <typename KernelType>
LDA_EXC_GENERATOR( host_eval_exc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto rho_i = rho + 2*i;
    traits::eval_exc_polar( rho_i[0], rho_i[1], eps[i] );

  }

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_exc_vxc_unpolar( rho[i], eps[i], vxc[i] );

  }

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto rho_i = rho + 2*i;
    auto vxc_i = vxc + 2*i;

    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], eps[i], 
      vxc_i[0], vxc_i[1] );

  }

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    double e;
    traits::eval_exc_unpolar( rho[i], e );
    eps[i] += scal_fact * e;

  }

}


template <typename KernelType>
LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto rho_i = rho + 2*i;

    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], e );
    
    eps[i] += scal_fact * e;

  }

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    double v,e;
    traits::eval_exc_vxc_unpolar( rho[i], e, v );
    eps[i] += scal_fact * e;
    vxc[i] += scal_fact * v;

  }

}


template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto rho_i = rho + 2*i;
    auto vxc_i = vxc + 2*i;

    double v_a, v_b, e;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], e, v_a, v_b);
    eps[i]   += scal_fact * e;
    vxc_i[0] += scal_fact * v_a;
    vxc_i[1] += scal_fact * v_b;

  }

}

template <typename KernelType>
LDA_FXC_GENERATOR( host_eval_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_fxc_unpolar( rho[i], fxc[i] );

  }

}

template <typename KernelType>
LDA_FXC_GENERATOR( host_eval_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho    + 2*i;
    auto* v2rho2_i = fxc + 3*i;

    traits::eval_fxc_polar( rho_i[0], rho_i[1], v2rho2_i[0], 
                           v2rho2_i[1], v2rho2_i[2] );

  }

}

template <typename KernelType>
LDA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_vxc_fxc_unpolar( rho[i], vxc[i], fxc[i] );

  }

}

template <typename KernelType>
LDA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho    + 2*i;
    auto* vrho_i   = vxc   + 2*i;
    auto* v2rho2_i = fxc + 3*i;

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], 
                               vrho_i[0], vrho_i[1],
                               v2rho2_i[0], v2rho2_i[1], v2rho2_i[2] );

  }

}

template <typename KernelType>
LDA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    double f;
    traits::eval_fxc_unpolar( rho[i], f );
    fxc[i] += scal_fact * f;
  }
}

template <typename KernelType>
LDA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i = rho + 2*i;
    auto* fxc_i = fxc + 3*i;
    double f_aa, f_ab, f_bb;
    traits::eval_fxc_polar( rho_i[0], rho_i[1], f_aa, f_ab, f_bb );
    fxc_i[0] += scal_fact * f_aa;
    fxc_i[1] += scal_fact * f_ab;
    fxc_i[2] += scal_fact * f_bb;
  }
}

template <typename KernelType>
LDA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    double v, f;
    traits::eval_vxc_fxc_unpolar( rho[i], v, f );
    vxc[i] += scal_fact * v;
    fxc[i] += scal_fact * f;
  }
}

template <typename KernelType>
LDA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i = rho + 2*i;
    auto* vxc_i = vxc + 2*i;
    auto* fxc_i = fxc + 3*i;
    double v_a, v_b, f_aa, f_ab, f_bb;
    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], v_a, v_b, f_aa, f_ab, f_bb );
    vxc_i[0] += scal_fact * v_a;
    vxc_i[1] += scal_fact * v_b;
    fxc_i[0] += scal_fact * f_aa;
    fxc_i[1] += scal_fact * f_ab;
    fxc_i[2] += scal_fact * f_bb;
  }
}

template <typename KernelType>
GGA_EXC_GENERATOR( host_eval_exc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_exc_unpolar( rho[i], sigma[i], eps[i] );

  }

}

template <typename KernelType>
GGA_EXC_GENERATOR( host_eval_exc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i   = rho   + 2*i;
    auto* sigma_i = sigma + 3*i;


    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], eps[i] );

  }

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_exc_vxc_unpolar( rho[i], sigma[i], 
      eps[i], vrho[i], vsigma[i] );

  }

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho   + 2*i;
    auto* sigma_i  = sigma + 3*i;
    auto* vrho_i   = vrho   + 2*i;
    auto* vsigma_i = vsigma + 3*i;

    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], eps[i], vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2] );

  }

}



template <typename KernelType>
GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    double e;
    traits::eval_exc_unpolar( rho[i], sigma[i], e );
    eps[i] += scal_fact * e;

  }

}


template <typename KernelType>
GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i   = rho   + 2*i;
    auto* sigma_i = sigma + 3*i;


    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], e );
    eps[i] += scal_fact * e;

  }

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    double e, vr, vs;
    traits::eval_exc_vxc_unpolar( rho[i], sigma[i], e, vr, vs );
    eps[i]    += scal_fact * e;
    vrho[i]   += scal_fact * vr;
    vsigma[i] += scal_fact * vs;

  }


}
template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho   + 2*i;
    auto* sigma_i  = sigma + 3*i;
    auto* vrho_i   = vrho   + 2*i;
    auto* vsigma_i = vsigma + 3*i;

    double e, vra, vrb, vsaa,vsab,vsbb;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], e, vra, vrb, vsaa, vsab, vsbb );

    eps[i]      += scal_fact * e;
    vrho_i[0]   += scal_fact * vra;
    vrho_i[1]   += scal_fact * vrb;
    vsigma_i[0] += scal_fact * vsaa;
    vsigma_i[1] += scal_fact * vsab;
    vsigma_i[2] += scal_fact * vsbb;

  }

}

template <typename KernelType>
GGA_FXC_GENERATOR( host_eval_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_fxc_unpolar( rho[i], sigma[i], v2rho2[i], v2rhosigma[i], v2sigma2[i] );

  }

}

template <typename KernelType>
GGA_FXC_GENERATOR( host_eval_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i        = rho        + 2*i;
    auto* sigma_i      = sigma      + 3*i;
    auto* v2rho2_i     = v2rho2     + 3*i;
    auto* v2rhosigma_i = v2rhosigma + 6*i;
    auto* v2sigma2_i   = v2sigma2   + 6*i;

    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
                           v2rho2_i[0], v2rho2_i[1], v2rho2_i[2], 
                           v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
                           v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
                           v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
                           v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5] );

  }

}

template <typename KernelType>
GGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    traits::eval_vxc_fxc_unpolar( rho[i], sigma[i], vrho[i], vsigma[i],
                                 v2rho2[i], v2rhosigma[i], v2sigma2[i] );

  }

}

template <typename KernelType>
GGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i        = rho        + 2*i;
    auto* sigma_i      = sigma      + 3*i;
    auto* vrho_i       = vrho       + 2*i;
    auto* vsigma_i     = vsigma     + 3*i;
    auto* v2rho2_i     = v2rho2     + 3*i;
    auto* v2rhosigma_i = v2rhosigma + 6*i;
    auto* v2sigma2_i   = v2sigma2   + 6*i;

    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
                               vrho_i[0], vrho_i[1], 
                               vsigma_i[0], vsigma_i[1], vsigma_i[2],
                               v2rho2_i[0], v2rho2_i[1], v2rho2_i[2], 
                               v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
                               v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
                               v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
                               v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5] );

  }

}

template <typename KernelType>
GGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    double f_rho2, f_rhosigma, f_sigma2;
    traits::eval_fxc_unpolar( rho[i], sigma[i], f_rho2, f_rhosigma, f_sigma2 );
    v2rho2[i]     += scal_fact * f_rho2;
    v2rhosigma[i] += scal_fact * f_rhosigma;
    v2sigma2[i]   += scal_fact * f_sigma2;
  }
}

template <typename KernelType>
GGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i        = rho        + 2*i;
    auto* sigma_i      = sigma      + 3*i;
    auto* v2rho2_i     = v2rho2     + 3*i;
    auto* v2rhosigma_i = v2rhosigma + 6*i;
    auto* v2sigma2_i   = v2sigma2   + 6*i;
    double f_rho2[3], f_rhosigma[6], f_sigma2[6];
    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      f_rho2[0], f_rho2[1], f_rho2[2],
      f_rhosigma[0], f_rhosigma[1], f_rhosigma[2], f_rhosigma[3], f_rhosigma[4], f_rhosigma[5],
      f_sigma2[0], f_sigma2[1], f_sigma2[2], f_sigma2[3], f_sigma2[4], f_sigma2[5] );
    for(int j=0;j<3;++j) v2rho2_i[j]     += scal_fact * f_rho2[j];
    for(int j=0;j<6;++j) v2rhosigma_i[j] += scal_fact * f_rhosigma[j];
    for(int j=0;j<6;++j) v2sigma2_i[j]   += scal_fact * f_sigma2[j];
  }
}

template <typename KernelType>
GGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    double vr, vs, f_rho2, f_rhosigma, f_sigma2;
    traits::eval_vxc_fxc_unpolar( rho[i], sigma[i], vr, vs, f_rho2, f_rhosigma, f_sigma2 );
    vrho[i]    += scal_fact * vr;
    vsigma[i]  += scal_fact * vs;
    v2rho2[i]     += scal_fact * f_rho2;
    v2rhosigma[i] += scal_fact * f_rhosigma;
    v2sigma2[i]   += scal_fact * f_sigma2;
  }
}

template <typename KernelType>
GGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i        = rho        + 2*i;
    auto* sigma_i      = sigma      + 3*i;
    auto* vrho_i       = vrho       + 2*i;
    auto* vsigma_i     = vsigma     + 3*i;
    auto* v2rho2_i     = v2rho2     + 3*i;
    auto* v2rhosigma_i = v2rhosigma + 6*i;
    auto* v2sigma2_i   = v2sigma2   + 6*i;
    double vr[2], vs[3], f_rho2[3], f_rhosigma[6], f_sigma2[6];
    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2],
      vr[0], vr[1], vs[0], vs[1], vs[2],
      f_rho2[0], f_rho2[1], f_rho2[2],
      f_rhosigma[0], f_rhosigma[1], f_rhosigma[2], f_rhosigma[3], f_rhosigma[4], f_rhosigma[5],
      f_sigma2[0], f_sigma2[1], f_sigma2[2], f_sigma2[3], f_sigma2[4], f_sigma2[5] );
    for(int j=0;j<2;++j) vrho_i[j]   += scal_fact * vr[j];
    for(int j=0;j<3;++j) vsigma_i[j] += scal_fact * vs[j];
    for(int j=0;j<3;++j) v2rho2_i[j]     += scal_fact * f_rho2[j];
    for(int j=0;j<6;++j) v2rhosigma_i[j] += scal_fact * f_rhosigma[j];
    for(int j=0;j<6;++j) v2sigma2_i[j]   += scal_fact * f_sigma2[j];
  }
}

template <typename KernelType>
MGGA_EXC_GENERATOR( host_eval_exc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    const double lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    traits::eval_exc_unpolar( rho[i], sigma[i], lapl_i, tau[i], eps[i] );

  }

}

template <typename KernelType>
MGGA_EXC_GENERATOR( host_eval_exc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i   = rho   + 2*i;
    auto* sigma_i = sigma + 3*i;
    auto* tau_i   = tau   + 2*i;
    auto* lapl_i  = lapl  + 2*i;

    const auto lapl_i_a = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_i_b = traits::needs_laplacian ? lapl_i[1] : 0.0;
   
    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_i_a, lapl_i_b, tau_i[0], tau_i[1], 
      eps[i] );

  }

}

template <typename KernelType>
MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    double vl;
    traits::eval_exc_vxc_unpolar( rho[i], sigma[i], lapl_i, tau[i],
      eps[i], vrho[i], vsigma[i], vl, vtau[i] );
    if(traits::needs_laplacian) vlapl[i] = vl;

  }

}

template <typename KernelType>
MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho   + 2*i;
    auto* sigma_i  = sigma + 3*i;
    auto* tau_i    = tau   + 2*i;
    auto* lapl_i   = lapl  + 2*i;
    auto* vrho_i   = vrho   + 2*i;
    auto* vsigma_i = vsigma + 3*i;
    auto* vtau_i   = vtau   + 2*i;
    auto* vlapl_i  = vlapl  + 2*i;


    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
    double vla, vlb;
    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0], tau_i[1], 
      eps[i], vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2], vla, vlb,
      vtau_i[0], vtau_i[1] );
    if(traits::needs_laplacian) {
      vlapl_i[0] = vla;
      vlapl_i[1] = vlb;
    }

  }

}



template <typename KernelType>
MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    double e;
    const double lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    traits::eval_exc_unpolar( rho[i], sigma[i], lapl_i, tau[i], e );
    eps[i] += scal_fact * e;

  }

}


template <typename KernelType>
MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i   = rho   + 2*i;
    auto* sigma_i = sigma + 3*i;
    auto* tau_i   = tau   + 2*i;
    auto* lapl_i  = lapl  + 2*i;

    const auto lapl_i_a = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_i_b = traits::needs_laplacian ? lapl_i[1] : 0.0;
   
    double e;
    traits::eval_exc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_i_a, lapl_i_b, tau_i[0], tau_i[1], 
      e );
    eps[i] += scal_fact * e;

  }

}

template <typename KernelType>
MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    double e, vr, vs, vl, vt;
    
    traits::eval_exc_vxc_unpolar( rho[i], sigma[i], lapl_i, tau[i], 
      e, vr, vs, vl, vt );
    eps[i]    += scal_fact * e;
    vrho[i]   += scal_fact * vr;
    vsigma[i] += scal_fact * vs;
    vtau[i]   += scal_fact * vt;
    if(traits::needs_laplacian) vlapl[i] += scal_fact * vl;

  }


}
template <typename KernelType>
MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i    = rho   + 2*i;
    auto* sigma_i  = sigma + 3*i;
    auto* tau_i    = tau   + 2*i;
    auto* lapl_i   = lapl  + 2*i;
    auto* vrho_i   = vrho   + 2*i;
    auto* vsigma_i = vsigma + 3*i;
    auto* vtau_i   = vtau   + 2*i;
    auto* vlapl_i  = vlapl  + 2*i;

    const double lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const double lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
    double e, vra, vrb, vsaa,vsab,vsbb, vla, vlb, vta, vtb;

    traits::eval_exc_vxc_polar( rho_i[0], rho_i[1], sigma_i[0], 
      sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0], tau_i[1], 
      e, vra, vrb, vsaa, vsab, vsbb, vla, vlb, vta, vtb );

    eps[i]      += scal_fact * e;
    vrho_i[0]   += scal_fact * vra;
    vrho_i[1]   += scal_fact * vrb;
    vsigma_i[0] += scal_fact * vsaa;
    vsigma_i[1] += scal_fact * vsab;
    vsigma_i[2] += scal_fact * vsbb;
    vtau_i[0]   += scal_fact * vta;
    vtau_i[1]   += scal_fact * vtb;
    if(traits::needs_laplacian) {
      vlapl_i[0] += scal_fact * vla;
      vlapl_i[1] += scal_fact * vlb;
    }

  }

}

template <typename KernelType>
MGGA_FXC_GENERATOR( host_eval_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    
    // Use local variables for laplacian-related second derivatives
    double local_v2rholapl, local_v2sigmalapl, local_v2lapl2, local_v2lapltau;
    
    traits::eval_fxc_unpolar( rho[i], sigma[i], lapl_i, tau[i], 
                             v2rho2[i], v2rhosigma[i], local_v2rholapl, v2rhotau[i], 
                             v2sigma2[i], local_v2sigmalapl, v2sigmatau[i], 
                             local_v2lapl2, local_v2lapltau, v2tau2[i] );
    
    // Only update laplacian-related second derivatives if needs_laplacian
    if (traits::needs_laplacian) {
      v2rholapl[i] = local_v2rholapl;
      v2sigmalapl[i] = local_v2sigmalapl;
      v2lapl2[i] = local_v2lapl2;
      v2lapltau[i] = local_v2lapltau;
    }
  }
}

template <typename KernelType>
MGGA_FXC_GENERATOR( host_eval_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;
  
  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i          = rho          + 2*i;
    auto* sigma_i        = sigma        + 3*i;
    auto* lapl_i         = lapl         + 2*i;
    auto* tau_i          = tau          + 2*i;
    auto* v2rho2_i       = v2rho2       + 3*i;
    auto* v2rhosigma_i   = v2rhosigma   + 6*i;
    auto* v2rhotau_i     = v2rhotau     + 4*i;
    auto* v2sigma2_i     = v2sigma2     + 6*i;
    auto* v2sigmatau_i   = v2sigmatau   + 6*i;
    auto* v2tau2_i       = v2tau2       + 3*i;

    const auto lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
    
    // Use local variables for laplacian-related second derivatives
    double local_v2rholapl_aa, local_v2rholapl_ab, local_v2rholapl_ba, local_v2rholapl_bb;
    double local_v2sigmalapl_aa_a, local_v2sigmalapl_aa_b, local_v2sigmalapl_ab_a, local_v2sigmalapl_ab_b, local_v2sigmalapl_bb_a, local_v2sigmalapl_bb_b;
    double local_v2lapl2_aa, local_v2lapl2_ab, local_v2lapl2_bb;
    double local_v2lapltau_aa, local_v2lapltau_ab, local_v2lapltau_ba, local_v2lapltau_bb;

    traits::eval_fxc_polar( rho_i[0], rho_i[1], 
                           sigma_i[0], sigma_i[1], sigma_i[2], 
                           lapl_a_use, lapl_b_use, 
                           tau_i[0], tau_i[1],
                           v2rho2_i[0], v2rho2_i[1], v2rho2_i[2],
                           v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
                           v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
                           local_v2rholapl_aa, local_v2rholapl_ab, local_v2rholapl_ba, local_v2rholapl_bb,
                           v2rhotau_i[0], v2rhotau_i[1], v2rhotau_i[2], v2rhotau_i[3],
                           v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
                           v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5],
                           local_v2sigmalapl_aa_a, local_v2sigmalapl_aa_b, local_v2sigmalapl_ab_a, 
                           local_v2sigmalapl_ab_b, local_v2sigmalapl_bb_a, local_v2sigmalapl_bb_b,
                           v2sigmatau_i[0], v2sigmatau_i[1], v2sigmatau_i[2],
                           v2sigmatau_i[3], v2sigmatau_i[4], v2sigmatau_i[5],
                           local_v2lapl2_aa, local_v2lapl2_ab, local_v2lapl2_bb,
                           local_v2lapltau_aa, local_v2lapltau_ab, local_v2lapltau_ba, local_v2lapltau_bb,
                           v2tau2_i[0], v2tau2_i[1], v2tau2_i[2] );
    
    // Only update laplacian-related second derivatives if needs_laplacian
    if (traits::needs_laplacian) {
      auto* v2rholapl_i    = v2rholapl    + 4*i;
      auto* v2sigmalapl_i  = v2sigmalapl  + 6*i;
      auto* v2lapl2_i      = v2lapl2      + 3*i;
      auto* v2lapltau_i    = v2lapltau    + 4*i;
      
      v2rholapl_i[0] = local_v2rholapl_aa;
      v2rholapl_i[1] = local_v2rholapl_ab;
      v2rholapl_i[2] = local_v2rholapl_ba;
      v2rholapl_i[3] = local_v2rholapl_bb;
      
      v2sigmalapl_i[0] = local_v2sigmalapl_aa_a;
      v2sigmalapl_i[1] = local_v2sigmalapl_aa_b;
      v2sigmalapl_i[2] = local_v2sigmalapl_ab_a;
      v2sigmalapl_i[3] = local_v2sigmalapl_ab_b;
      v2sigmalapl_i[4] = local_v2sigmalapl_bb_a;
      v2sigmalapl_i[5] = local_v2sigmalapl_bb_b;
      
      v2lapl2_i[0] = local_v2lapl2_aa;
      v2lapl2_i[1] = local_v2lapl2_ab;
      v2lapl2_i[2] = local_v2lapl2_bb;
      
      v2lapltau_i[0] = local_v2lapltau_aa;
      v2lapltau_i[1] = local_v2lapltau_ab;
      v2lapltau_i[2] = local_v2lapltau_ba;
      v2lapltau_i[3] = local_v2lapltau_bb;
    }
  }
}

template <typename KernelType>
MGGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar ) {

  using traits = kernel_traits<KernelType>;

  for( int32_t i = 0; i < N; ++i ) {

    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;

    /* First–derivative w.r.t. laplacian */
    double local_vlapl;

    /* Second–derivatives involving the laplacian */
    double local_v2rholapl, local_v2sigmalapl, local_v2lapl2, local_v2lapltau;

    traits::eval_vxc_fxc_unpolar( rho[i], sigma[i], lapl_i, tau[i],
                                  vrho[i], vsigma[i], local_vlapl, vtau[i],
                                  v2rho2[i], v2rhosigma[i], local_v2rholapl,
                                  v2rhotau[i], v2sigma2[i], local_v2sigmalapl,
                                  v2sigmatau[i], local_v2lapl2,
                                  local_v2lapltau, v2tau2[i] );

    if( traits::needs_laplacian ) {
      vlapl[i]        = local_vlapl;
      v2rholapl[i]    = local_v2rholapl;
      v2sigmalapl[i]  = local_v2sigmalapl;
      v2lapl2[i]      = local_v2lapl2;
      v2lapltau[i]    = local_v2lapltau;
    }

  }

}



template <typename KernelType>
MGGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar ) {

  using traits = kernel_traits<KernelType>;

  for( int32_t i = 0; i < N; ++i ) {

    auto* rho_i          = rho          + 2*i;
    auto* sigma_i        = sigma        + 3*i;
    auto* lapl_i         = lapl         + 2*i;
    auto* tau_i          = tau          + 2*i;
    auto* vrho_i         = vrho         + 2*i;
    auto* vsigma_i       = vsigma       + 3*i;
    auto* vlapl_i        = vlapl        + 2*i;
    auto* vtau_i         = vtau         + 2*i;
    auto* v2rho2_i       = v2rho2       + 3*i;
    auto* v2rhosigma_i   = v2rhosigma   + 6*i;
    auto* v2rhotau_i     = v2rhotau     + 4*i;
    auto* v2sigma2_i     = v2sigma2     + 6*i;
    auto* v2sigmatau_i   = v2sigmatau   + 6*i;
    auto* v2tau2_i       = v2tau2       + 3*i;

    const auto lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;

    /* First–derivative w.r.t. laplacian */
    double local_vlapl_a, local_vlapl_b;

    /* Second–derivatives involving the laplacian */
    double local_v2rholapl_aa, local_v2rholapl_ab, local_v2rholapl_ba, local_v2rholapl_bb;

    double local_v2sigmalapl_aa_a, local_v2sigmalapl_aa_b, local_v2sigmalapl_ab_a, local_v2sigmalapl_ab_b, local_v2sigmalapl_bb_a, local_v2sigmalapl_bb_b;

    double local_v2lapl2_aa,  local_v2lapl2_ab,  local_v2lapl2_bb;
    double local_v2lapltau_aa, local_v2lapltau_ab, local_v2lapltau_ba, local_v2lapltau_bb;

    traits::eval_vxc_fxc_polar(
      rho_i[0], rho_i[1],
      sigma_i[0], sigma_i[1], sigma_i[2],
      lapl_a_use, lapl_b_use,
      tau_i[0], tau_i[1],
      vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2],
      local_vlapl_a, local_vlapl_b,
      vtau_i[0], vtau_i[1],
      v2rho2_i[0], v2rho2_i[1], v2rho2_i[2],
      v2rhosigma_i[0], v2rhosigma_i[1], v2rhosigma_i[2],
      v2rhosigma_i[3], v2rhosigma_i[4], v2rhosigma_i[5],
      local_v2rholapl_aa, local_v2rholapl_ab, local_v2rholapl_ba, local_v2rholapl_bb,
      v2rhotau_i[0], v2rhotau_i[1], v2rhotau_i[2], v2rhotau_i[3],
      v2sigma2_i[0], v2sigma2_i[1], v2sigma2_i[2],
      v2sigma2_i[3], v2sigma2_i[4], v2sigma2_i[5],
      local_v2sigmalapl_aa_a, local_v2sigmalapl_aa_b, local_v2sigmalapl_ab_a,
      local_v2sigmalapl_ab_b, local_v2sigmalapl_bb_a, local_v2sigmalapl_bb_b,
      v2sigmatau_i[0], v2sigmatau_i[1], v2sigmatau_i[2],
      v2sigmatau_i[3], v2sigmatau_i[4], v2sigmatau_i[5],
      local_v2lapl2_aa, local_v2lapl2_ab, local_v2lapl2_bb,
      local_v2lapltau_aa, local_v2lapltau_ab, local_v2lapltau_ba, local_v2lapltau_bb,
      v2tau2_i[0], v2tau2_i[1], v2tau2_i[2] );

    if( traits::needs_laplacian ) {
      /* first derivative */
      vlapl_i[0] = local_vlapl_a;
      vlapl_i[1] = local_vlapl_b;

      /* second derivatives */
      auto* v2rholapl_i   = v2rholapl   + 4*i;
      auto* v2sigmalapl_i = v2sigmalapl + 6*i;
      auto* v2lapl2_i     = v2lapl2     + 3*i;
      auto* v2lapltau_i   = v2lapltau   + 4*i;

      v2rholapl_i[0] = local_v2rholapl_aa;
      v2rholapl_i[1] = local_v2rholapl_ab;
      v2rholapl_i[2] = local_v2rholapl_ba;
      v2rholapl_i[3] = local_v2rholapl_bb;

      v2sigmalapl_i[0] = local_v2sigmalapl_aa_a;
      v2sigmalapl_i[1] = local_v2sigmalapl_aa_b;
      v2sigmalapl_i[2] = local_v2sigmalapl_ab_a;
      v2sigmalapl_i[3] = local_v2sigmalapl_ab_b;
      v2sigmalapl_i[4] = local_v2sigmalapl_bb_a;
      v2sigmalapl_i[5] = local_v2sigmalapl_bb_b;

      v2lapl2_i[0] = local_v2lapl2_aa;
      v2lapl2_i[1] = local_v2lapl2_ab;
      v2lapl2_i[2] = local_v2lapl2_bb;

      v2lapltau_i[0] = local_v2lapltau_aa;
      v2lapltau_i[1] = local_v2lapltau_ab;
      v2lapltau_i[2] = local_v2lapltau_ba;
      v2lapltau_i[3] = local_v2lapltau_bb;
    }

  }

}

template <typename KernelType>
MGGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    double f_rho2, f_rhosigma, f_rholapl, f_rhotau, f_sigma2, f_sigmalapl, f_sigmatau, f_lapl2, f_lapltau, f_tau2;
    traits::eval_fxc_unpolar( rho[i], sigma[i], lapl_i, tau[i],
      f_rho2, f_rhosigma, f_rholapl, f_rhotau,
      f_sigma2, f_sigmalapl, f_sigmatau,
      f_lapl2, f_lapltau, f_tau2 );
    v2rho2[i]     += scal_fact * f_rho2;
    v2rhosigma[i] += scal_fact * f_rhosigma;
    v2rhotau[i]   += scal_fact * f_rhotau;
    v2sigma2[i]   += scal_fact * f_sigma2;
    v2sigmatau[i] += scal_fact * f_sigmatau;
    v2tau2[i]     += scal_fact * f_tau2;
    if(traits::needs_laplacian) {
      v2rholapl[i] += scal_fact * f_rholapl;
      v2sigmalapl[i] += scal_fact * f_sigmalapl;
      v2lapl2[i] += scal_fact * f_lapl2;
      v2lapltau[i] += scal_fact * f_lapltau;
    }
  }
}

template <typename KernelType>
MGGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i          = rho          + 2*i;
    auto* sigma_i        = sigma        + 3*i;
    auto* lapl_i         = lapl         + 2*i;
    auto* tau_i          = tau          + 2*i;
    auto* v2rho2_i       = v2rho2       + 3*i;
    auto* v2rhosigma_i   = v2rhosigma   + 6*i;
    auto* v2rhotau_i     = v2rhotau     + 4*i;
    auto* v2sigma2_i     = v2sigma2     + 6*i;
    auto* v2sigmatau_i   = v2sigmatau   + 6*i;
    auto* v2tau2_i       = v2tau2       + 3*i;
    const auto lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
    double f_rho2[3], f_rhosigma[6], f_rholapl[4], f_rhotau[4], f_sigma2[6], f_sigmalapl[6], f_sigmatau[6], f_lapl2[3], f_lapltau[4], f_tau2[3];
    traits::eval_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
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

    for(int j=0;j<3;++j) v2rho2_i[j]     += scal_fact * f_rho2[j];
    for(int j=0;j<6;++j) v2rhosigma_i[j] += scal_fact * f_rhosigma[j];
    for(int j=0;j<4;++j) v2rhotau_i[j]   += scal_fact * f_rhotau[j];
    for(int j=0;j<6;++j) v2sigma2_i[j]   += scal_fact * f_sigma2[j];
    for(int j=0;j<6;++j) v2sigmatau_i[j] += scal_fact * f_sigmatau[j];
    for(int j=0;j<3;++j) v2tau2_i[j]     += scal_fact * f_tau2[j];

    if(traits::needs_laplacian) {
      for(int j=0;j<4;++j) v2rholapl[4*i+j] += scal_fact * f_rholapl[j];
      for(int j=0;j<6;++j) v2sigmalapl[6*i+j] += scal_fact * f_sigmalapl[j];
      for(int j=0;j<3;++j) v2lapl2[3*i+j] += scal_fact * f_lapl2[j];
      for(int j=0;j<4;++j) v2lapltau[4*i+j] += scal_fact * f_lapltau[j];
    }
  }
}

template <typename KernelType>
MGGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    const auto lapl_i = traits::needs_laplacian ? lapl[i] : 0.0;
    double vr, vs, vl, vt, f_rho2, f_rhosigma, f_rholapl, f_rhotau, f_sigma2, f_sigmalapl, f_sigmatau, f_lapl2, f_lapltau, f_tau2;
    traits::eval_vxc_fxc_unpolar( rho[i], sigma[i], lapl_i, tau[i],
      vr, vs, vl, vt,
      f_rho2, f_rhosigma, f_rholapl, f_rhotau,
      f_sigma2, f_sigmalapl, f_sigmatau,
      f_lapl2, f_lapltau, f_tau2 );

    vrho[i]    += scal_fact * vr;
    vsigma[i]  += scal_fact * vs;
    vtau[i]    += scal_fact * vt;
    v2rho2[i]     += scal_fact * f_rho2;
    v2rhosigma[i] += scal_fact * f_rhosigma;
    v2rhotau[i]   += scal_fact * f_rhotau;
    v2sigma2[i]   += scal_fact * f_sigma2;
    v2sigmatau[i] += scal_fact * f_sigmatau;
    v2tau2[i]     += scal_fact * f_tau2;

    if(traits::needs_laplacian) {
      vlapl[i] += scal_fact * vl;
      v2rholapl[i] += scal_fact * f_rholapl;
      v2sigmalapl[i] += scal_fact * f_sigmalapl;
      v2lapl2[i] += scal_fact * f_lapl2;
      v2lapltau[i] += scal_fact * f_lapltau;
    }
  }
}

template <typename KernelType>
MGGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar ) {
  using traits = kernel_traits<KernelType>;
  for( int32_t i = 0; i < N; ++i ) {
    auto* rho_i          = rho          + 2*i;
    auto* sigma_i        = sigma        + 3*i;
    auto* lapl_i         = lapl         + 2*i;
    auto* tau_i          = tau          + 2*i;
    auto* vrho_i         = vrho         + 2*i;
    auto* vsigma_i       = vsigma       + 3*i;
    auto* vlapl_i        = vlapl        + 2*i;
    auto* vtau_i         = vtau         + 2*i;
    auto* v2rho2_i       = v2rho2       + 3*i;
    auto* v2rhosigma_i   = v2rhosigma   + 6*i;
    auto* v2rhotau_i     = v2rhotau     + 4*i;
    auto* v2sigma2_i     = v2sigma2     + 6*i;
    auto* v2sigmatau_i   = v2sigmatau   + 6*i;
    auto* v2tau2_i       = v2tau2       + 3*i;
    const auto lapl_a_use = traits::needs_laplacian ? lapl_i[0] : 0.0;
    const auto lapl_b_use = traits::needs_laplacian ? lapl_i[1] : 0.0;
    double vr[2], vs[3], vl[2], vt[2], f_rho2[3], f_rhosigma[6], f_rholapl[4], f_rhotau[4], f_sigma2[6], f_sigmalapl[6], f_sigmatau[6], f_lapl2[3], f_lapltau[4], f_tau2[3];
    traits::eval_vxc_fxc_polar( rho_i[0], rho_i[1], sigma_i[0], sigma_i[1], sigma_i[2], lapl_a_use, lapl_b_use, tau_i[0], tau_i[1],
      vr[0], vr[1], vs[0], vs[1], vs[2], vl[0], vl[1], vt[0], vt[1],
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
      
    for(int j=0;j<2;++j) vrho_i[j]   += scal_fact * vr[j];
    for(int j=0;j<3;++j) vsigma_i[j] += scal_fact * vs[j];
    for(int j=0;j<2;++j) vtau_i[j]   += scal_fact * vt[j];
    
    for(int j=0;j<3;++j) v2rho2_i[j]     += scal_fact * f_rho2[j];
    for(int j=0;j<6;++j) v2rhosigma_i[j] += scal_fact * f_rhosigma[j];
    for(int j=0;j<4;++j) v2rhotau_i[j]   += scal_fact * f_rhotau[j];
    for(int j=0;j<6;++j) v2sigma2_i[j]   += scal_fact * f_sigma2[j];
    for(int j=0;j<6;++j) v2sigmatau_i[j] += scal_fact * f_sigmatau[j];
    for(int j=0;j<3;++j) v2tau2_i[j]     += scal_fact * f_tau2[j];

    if(traits::needs_laplacian) {
      for(int j=0;j<2;++j) vlapl_i[j]  += scal_fact * vl[j];
      for(int j=0;j<4;++j) v2rholapl[4*i+j] += scal_fact * f_rholapl[j];
      for(int j=0;j<6;++j) v2sigmalapl[6*i+j] += scal_fact * f_sigmalapl[j];
      for(int j=0;j<3;++j) v2lapl2[3*i+j] += scal_fact * f_lapl2[j];
      for(int j=0;j<4;++j) v2lapltau[4*i+j] += scal_fact * f_lapltau[j];
    }
  }
}

#define LDA_GENERATE_HOST_HELPERS(KERN) \
  template LDA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template LDA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template LDA_FXC_GENERATOR( host_eval_fxc_helper_unpolar<KERN> ); \
  template LDA_FXC_GENERATOR( host_eval_fxc_helper_polar<KERN> ); \
  template LDA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template LDA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar<KERN> ); \
  template LDA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar<KERN> ); \
  template LDA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar<KERN> ); \
  template LDA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template LDA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar<KERN> );

#define GGA_GENERATE_HOST_HELPERS(KERN) \
  template GGA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template GGA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template GGA_FXC_GENERATOR( host_eval_fxc_helper_unpolar<KERN> ); \
  template GGA_FXC_GENERATOR( host_eval_fxc_helper_polar<KERN> ); \
  template GGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template GGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar<KERN> ); \
  template GGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar<KERN> ); \
  template GGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar<KERN> ); \
  template GGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template GGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar<KERN> );

#define MGGA_GENERATE_HOST_HELPERS(KERN) \
  template MGGA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template MGGA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); \
  template MGGA_FXC_GENERATOR( host_eval_fxc_helper_unpolar<KERN> ); \
  template MGGA_FXC_GENERATOR( host_eval_fxc_helper_polar<KERN> ); \
  template MGGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_unpolar<KERN> ); \
  template MGGA_VXC_FXC_GENERATOR( host_eval_vxc_fxc_helper_polar<KERN> ); \
  template MGGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_unpolar<KERN> ); \
  template MGGA_FXC_INC_GENERATOR( host_eval_fxc_inc_helper_polar<KERN> ); \
  template MGGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_unpolar<KERN> ); \
  template MGGA_VXC_FXC_INC_GENERATOR( host_eval_vxc_fxc_inc_helper_polar<KERN> );

LDA_GENERATE_HOST_HELPERS( BuiltinSlaterExchange )
LDA_GENERATE_HOST_HELPERS( BuiltinVWN3 )
LDA_GENERATE_HOST_HELPERS( BuiltinVWN_RPA )
LDA_GENERATE_HOST_HELPERS( BuiltinVWN)
LDA_GENERATE_HOST_HELPERS( BuiltinPW91_LDA )
LDA_GENERATE_HOST_HELPERS( BuiltinPW91_LDA_MOD )
LDA_GENERATE_HOST_HELPERS( BuiltinPW91_LDA_RPA )
LDA_GENERATE_HOST_HELPERS( BuiltinPZ81 )
LDA_GENERATE_HOST_HELPERS( BuiltinPZ81_MOD )

GGA_GENERATE_HOST_HELPERS( BuiltinB88   )
GGA_GENERATE_HOST_HELPERS( BuiltinLYP   )
GGA_GENERATE_HOST_HELPERS( BuiltinPBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinRevPBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinPBE_C )
GGA_GENERATE_HOST_HELPERS( BuiltinB97_D )
GGA_GENERATE_HOST_HELPERS( BuiltinITYH_X )
GGA_GENERATE_HOST_HELPERS( BuiltinITYH_X_033 )
GGA_GENERATE_HOST_HELPERS( BuiltinITYH_X_015 )
GGA_GENERATE_HOST_HELPERS( BuiltinP86_C )
GGA_GENERATE_HOST_HELPERS( BuiltinP86VWN_FT_C )
GGA_GENERATE_HOST_HELPERS( BuiltinPW91_C )
GGA_GENERATE_HOST_HELPERS( BuiltinPBE_SOL_C )
GGA_GENERATE_HOST_HELPERS( BuiltinBMK_C )
GGA_GENERATE_HOST_HELPERS( BuiltinN12_C )
GGA_GENERATE_HOST_HELPERS( BuiltinN12_SX_C )
GGA_GENERATE_HOST_HELPERS( BuiltinSOGGA11_X_C )
GGA_GENERATE_HOST_HELPERS( BuiltinPW91_X )
GGA_GENERATE_HOST_HELPERS( BuiltinMPW91_X )
GGA_GENERATE_HOST_HELPERS( BuiltinOPTX_X )
GGA_GENERATE_HOST_HELPERS( BuiltinRPBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinSOGGA11_X_X )
GGA_GENERATE_HOST_HELPERS( BuiltinPW86_X )
GGA_GENERATE_HOST_HELPERS( BuiltinWB97_XC )
GGA_GENERATE_HOST_HELPERS( BuiltinWB97X_XC )
GGA_GENERATE_HOST_HELPERS( BuiltinWB97X_V_XC )
GGA_GENERATE_HOST_HELPERS( BuiltinWB97X_D_XC )
GGA_GENERATE_HOST_HELPERS( BuiltinWB97X_D3_XC )
GGA_GENERATE_HOST_HELPERS( BuiltinHJS_PBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinLCwPBE_wPBEh_X )
GGA_GENERATE_HOST_HELPERS( BuiltinLRCwPBE_HJS_PBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinLRCwPBEh_HJS_PBE_X )
GGA_GENERATE_HOST_HELPERS( BuiltinWPBEh_X_default0 )
GGA_GENERATE_HOST_HELPERS( BuiltinHSE03_wPBEh_X )
GGA_GENERATE_HOST_HELPERS( BuiltinHSE06_wPBEh_X )


MGGA_GENERATE_HOST_HELPERS( BuiltinSCAN_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinSCAN_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCAN_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCAN_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinFT98_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM062X_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM062X_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinPKZB_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinPKZB_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinTPSS_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinRevTPSS_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_L_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_HF_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinRevM06_L_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_SX_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_L_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_HF_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinRevM06_L_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM06_SX_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM05_2X_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM05_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM08_HX_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM08_SO_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinCF22D_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinM11_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN12_L_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN12_SX_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN15_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN15_L_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinTPSS_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinRevTPSS_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinRSCAN_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinBC95_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinMBEEF_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinRSCAN_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinBMK_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM08_HX_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM08_SO_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN12_L_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN15_L_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN15_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinCF22D_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinMN12_SX_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM11_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM05_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinM05_2X_X )


MGGA_GENERATE_HOST_HELPERS( BuiltinPC07_K )
MGGA_GENERATE_HOST_HELPERS( BuiltinPC07OPT_K )

MGGA_GENERATE_HOST_HELPERS( BuiltinSCANL_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinSCANL_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCANL_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCANL_X )

LDA_GENERATE_HOST_HELPERS( BuiltinEPC17_1 )
LDA_GENERATE_HOST_HELPERS( BuiltinEPC17_2 )
LDA_GENERATE_HOST_HELPERS( BuiltinEPC18_1 )
LDA_GENERATE_HOST_HELPERS( BuiltinEPC18_2 )

}
}

