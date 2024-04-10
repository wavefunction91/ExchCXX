/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
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

// GGA interface
UNUSED_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )

// MGGA interface
UNUSED_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )


// INC interfaces
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )

#ifdef EXCHCXX_ENABLE_DEVICE

// LDA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )

// GGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )

// MGGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )


// INC interfaces
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )

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

#define LDA_GENERATE_HOST_HELPERS(KERN) \
  template LDA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template LDA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); 

#define GGA_GENERATE_HOST_HELPERS(KERN) \
  template GGA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template GGA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); 

#define MGGA_GENERATE_HOST_HELPERS(KERN) \
  template MGGA_EXC_GENERATOR( host_eval_exc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_unpolar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_unpolar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_unpolar<KERN> );\
  template MGGA_EXC_GENERATOR( host_eval_exc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper_polar<KERN> ); \
  template MGGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper_polar<KERN> ); \
  template MGGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper_polar<KERN> ); 

LDA_GENERATE_HOST_HELPERS( BuiltinSlaterExchange )
LDA_GENERATE_HOST_HELPERS( BuiltinVWN3 )
LDA_GENERATE_HOST_HELPERS( BuiltinVWN_RPA )
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

GGA_GENERATE_HOST_HELPERS( BuiltinB3LYP  )
GGA_GENERATE_HOST_HELPERS( BuiltinPBE0  )

MGGA_GENERATE_HOST_HELPERS( BuiltinSCAN_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinSCAN_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCAN_X )
MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCAN_C )
MGGA_GENERATE_HOST_HELPERS( BuiltinFT98_X )

MGGA_GENERATE_HOST_HELPERS( BuiltinPC07_K )
MGGA_GENERATE_HOST_HELPERS( BuiltinPC07OPT_K )

MGGA_GENERATE_HOST_HELPERS( BuiltinSCANL_C )
//MGGA_GENERATE_HOST_HELPERS( BuiltinSCANL_X )
//MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCANL_C )
//MGGA_GENERATE_HOST_HELPERS( BuiltinR2SCANL_X )

}
}

