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

#include <exchcxx/xc_functional.hpp>
#include <exchcxx/util/div_ceil.hpp>
#include <string>


__global__ void scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if( tid < N ) Y_device[tid] = X_device[tid] * fact;

}

__global__ void add_scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
  if( tid < N ) Y_device[tid] += X_device[tid] * fact;

}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 1024;
  int blocks  = ExchCXX::util::div_ceil(N,1024);
  hipLaunchKernelGGL(scal_kernel, dim3(blocks), dim3(threads ), 0, 0,  N, fact, X_device, Y_device );
}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device, hipStream_t& stream ) {
  int threads = 1024;
  int blocks  = ExchCXX::util::div_ceil(N,1024);
  hipLaunchKernelGGL(scal_kernel, dim3(blocks), dim3(threads), 0, stream ,  N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 1024;
  int blocks  = ExchCXX::util::div_ceil(N,1024);
  hipLaunchKernelGGL(add_scal_kernel, dim3(blocks), dim3(threads ), 0, 0,  N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device, hipStream_t& stream ) {
  int threads = 1024;
  int blocks  = ExchCXX::util::div_ceil(N,1024);
  hipLaunchKernelGGL(add_scal_kernel, dim3(blocks), dim3(threads), 0, stream ,  N, fact, X_device, Y_device );
}


template <typename T = double>
T* safe_hip_malloc( size_t N ) {

  T* ptr = nullptr;
  auto stat = hipMalloc( &ptr, N*sizeof(T) );
  if( stat != hipSuccess ) throw std::runtime_error("Alloc Failed");

  return ptr;

}

template <typename T>
void safe_zero( size_t len, T* ptr, hipStream_t stream ) {
  auto stat = hipMemsetAsync( ptr, 0, len*sizeof(T), stream );
  if( stat != hipSuccess ) 
    throw std::runtime_error("Memset Failed : " + std::string(hipGetErrorString( stat )));
}

namespace ExchCXX {






LDA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_hip_malloc( len_exc_buffer );

  safe_zero( len_exc_buffer, eps, stream );


  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_inc_device(
        kernels_[i].first, N, rho, eps, stream
      );

    } else {

      double* eps_eval = i ? eps_scr : eps;
      kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

      if( i ) 
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );

    }
  
  }

  if( eps_scr ) hipFree( eps_scr );

}


LDA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );
  size_t len_vxc_buffer = vrho_buffer_len( N );

  double* eps_scr(nullptr), *vxc_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr = safe_hip_malloc( len_exc_buffer );
    vxc_scr = safe_hip_malloc( len_vxc_buffer );
  }

  safe_zero( len_exc_buffer, eps, stream );
  safe_zero( len_vxc_buffer, vxc, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      kernels_[i].second.eval_exc_vxc_inc_device(
        kernels_[i].first, N, rho, eps, vxc, stream
      );

    } else {

      double* eps_eval = i ? eps_scr : eps;
      double* vxc_eval = i ? vxc_scr : vxc;
      kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vxc_eval, stream);

      if( i ) {

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        add_scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, stream );

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, stream );

      }

    }
  
  }

  if( eps_scr ) hipFree( eps_scr );
  if( vxc_scr ) hipFree( vxc_scr );

}



// GGA Interfaces

GGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_hip_malloc( len_exc_buffer );

  safe_zero( len_exc_buffer, eps, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, stream
        );
      else
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, eps, stream
        );

    } else {

      double* eps_eval = i ? eps_scr : eps;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, stream);
      else
        kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

      if( i ) 
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
  
    }
  }

  if( eps_scr ) hipFree( eps_scr );

}


GGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer    = exc_buffer_len(N);
  size_t len_vrho_buffer   = vrho_buffer_len(N);
  size_t len_vsigma_buffer = vsigma_buffer_len(N);

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr    = safe_hip_malloc( len_exc_buffer );
    vrho_scr   = safe_hip_malloc( len_vrho_buffer );
    vsigma_scr = safe_hip_malloc( len_vsigma_buffer );
  }

  safe_zero( len_exc_buffer,    eps,    stream );
  safe_zero( len_vrho_buffer,   vrho,   stream );
  safe_zero( len_vsigma_buffer, vsigma, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma, stream 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, eps, vrho, stream
        );

    } else {

      double* eps_eval    = i ? eps_scr    : eps;
      double* vrho_eval   = i ? vrho_scr   : vrho;
      double* vsigma_eval = i ? vsigma_scr : vsigma;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_device(N, rho, sigma, eps_eval, vrho_eval, 
          vsigma_eval, stream );
      else
        kernels_[i].second.eval_exc_vxc_device(N, rho, eps_eval, vrho_eval, stream);

      if( i ) {

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        add_scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream);
        if( kernels_[i].second.is_gga() )
          add_scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream );
        if( kernels_[i].second.is_gga() )
          scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );

      }

    }
  }

  if( eps_scr )    hipFree( eps_scr );
  if( vrho_scr )   hipFree( vrho_scr );
  if( vsigma_scr ) hipFree( vsigma_scr );

}




// mGGA Interfaces

MGGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_hip_malloc( len_exc_buffer );

  safe_zero( len_exc_buffer,    eps,    stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps, stream
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, stream
        );
      else
        kernels_[i].second.eval_exc_inc_device(
          kernels_[i].first, N, rho, eps, stream
        );

    } else {

      double* eps_eval = i ? eps_scr : eps;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, lapl, tau, eps_eval, stream);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_device(N, rho, sigma, eps_eval, stream);
      else
        kernels_[i].second.eval_exc_device(N, rho, eps_eval, stream);

      if( i ) 
        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
      else
        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
  
    }
  }

  if( eps_scr ) hipFree( eps_scr );

}


MGGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  size_t len_exc_buffer    = exc_buffer_len(N);
  size_t len_vrho_buffer   = vrho_buffer_len(N);
  size_t len_vsigma_buffer = vsigma_buffer_len(N);
  size_t len_vlapl_buffer   = vlapl_buffer_len(N);
  size_t len_vtau_buffer   = vtau_buffer_len(N);

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr), 
    *vlapl_scr(nullptr), *vtau_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr    = safe_hip_malloc( len_exc_buffer );
    vrho_scr   = safe_hip_malloc( len_vrho_buffer );
    vsigma_scr = safe_hip_malloc( len_vsigma_buffer );
    vlapl_scr  = safe_hip_malloc( len_vlapl_buffer );
    vtau_scr   = safe_hip_malloc( len_vtau_buffer );
  }

  safe_zero( len_exc_buffer, eps, stream );
  safe_zero( len_vrho_buffer, vrho, stream );
  safe_zero( len_vsigma_buffer, vsigma, stream );
  safe_zero( len_vlapl_buffer, vlapl, stream );
  safe_zero( len_vtau_buffer, vtau, stream );
  

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( supports_inc_interface() ) {

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, eps, 
          vrho, vsigma, vlapl, vtau, stream 
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, sigma, eps, vrho, 
          vsigma, stream 
        );
      else
        kernels_[i].second.eval_exc_vxc_inc_device(
          kernels_[i].first, N, rho, eps, vrho, stream
        );

    } else {

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

        add_scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        add_scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream );

        if( kernels_[i].second.is_gga() )
          add_scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );

        if( kernels_[i].second.is_mgga() ) {
          add_scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream );
          add_scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, stream  );
        }

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream );

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() )
          scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );

        if( kernels_[i].second.is_mgga() ) {
          scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream );
          scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, stream  );
        }

      }
    }
  }

  if( eps_scr )    hipFree( eps_scr );
  if( vrho_scr )   hipFree( vrho_scr );
  if( vsigma_scr ) hipFree( vsigma_scr );
  if( vlapl_scr )  hipFree( vlapl_scr );
  if( vtau_scr )   hipFree( vtau_scr );
}

}
