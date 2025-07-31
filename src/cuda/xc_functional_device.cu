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

#include <exchcxx/xc_functional.hpp>
#include <exchcxx/util/div_ceil.hpp>
#include <string>


__global__ void scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if( tid < N ) Y_device[tid] = X_device[tid] * fact;

}

__global__ void add_scal_kernel( const int N, const double fact, const double* X_device, double* Y_device ) {

  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if( tid < N ) Y_device[tid] += X_device[tid] * fact;

}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 512;
  int blocks  = ExchCXX::util::div_ceil(N,512);
  scal_kernel<<< blocks, threads >>>( N, fact, X_device, Y_device );
}

void scal_device( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
  int threads = 512;
  int blocks  = ExchCXX::util::div_ceil(N,512);
  scal_kernel<<< blocks, threads, 0, stream >>>( N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device ) {
  int threads = 512;
  int blocks  = ExchCXX::util::div_ceil(N,512);
  add_scal_kernel<<< blocks, threads >>>( N, fact, X_device, Y_device );
}

void add_scal_device( const int N, const double fact, const double* X_device, double* Y_device, cudaStream_t& stream ) {
  int threads = 512;
  int blocks  = ExchCXX::util::div_ceil(N,512);
  add_scal_kernel<<< blocks, threads, 0, stream >>>( N, fact, X_device, Y_device );
}


template <typename T = double>
T* safe_cuda_malloc( size_t N ) {

  T* ptr = nullptr;
  auto stat = cudaMalloc( &ptr, N*sizeof(T) );
  if( stat != cudaSuccess ) throw std::runtime_error("Alloc Failed");

  return ptr;

}

template <typename T>
void safe_zero( size_t len, T* ptr, cudaStream_t stream ) {
  auto stat = cudaMemsetAsync( ptr, 0, len*sizeof(T), stream );
  if( stat != cudaSuccess ) 
    throw std::runtime_error("Memset Failed : " + std::string(cudaGetErrorString( stat )));
}

namespace ExchCXX {






LDA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_cuda_malloc( len_exc_buffer );

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

  if( eps_scr ) cudaFree( eps_scr );

}


LDA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  size_t len_exc_buffer = exc_buffer_len( N );
  size_t len_vxc_buffer = vrho_buffer_len( N );

  double* eps_scr(nullptr), *vxc_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr = safe_cuda_malloc( len_exc_buffer );
    vxc_scr = safe_cuda_malloc( len_vxc_buffer );
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

  if( eps_scr ) cudaFree( eps_scr );
  if( vxc_scr ) cudaFree( vxc_scr );

}


LDA_FXC_GENERATOR_DEVICE( XCFunctional::eval_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_fxc_buffer = v2rho2_buffer_len(N);

  double* fxc_scr = nullptr;
  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) 
    fxc_scr = safe_cuda_malloc( len_fxc_buffer );

  safe_zero( len_fxc_buffer, fxc, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {
    if (use_inc) {
      kernels_[i].second.eval_fxc_inc_device(
        kernels_[i].first, N, rho, fxc, stream
      );
    } else {
      double* fxc_eval = i ? fxc_scr : fxc;
      kernels_[i].second.eval_fxc_device(N, rho, fxc_eval, stream);

      if( i ) 
        add_scal_device( len_fxc_buffer, kernels_[i].first, fxc_eval, fxc, stream );
      else
        scal_device( len_fxc_buffer, kernels_[i].first, fxc_eval, fxc, stream );
    }
  }

  if( fxc_scr ) cudaFree( fxc_scr );
}

LDA_VXC_FXC_GENERATOR_DEVICE( XCFunctional::eval_vxc_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT LDA",  is_lda() );

  const size_t len_vxc_buffer = vrho_buffer_len(N);
  const size_t len_fxc_buffer = v2rho2_buffer_len(N);

  double* vxc_scr(nullptr), *fxc_scr(nullptr);
  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) {
    vxc_scr = safe_cuda_malloc( len_vxc_buffer );
    fxc_scr = safe_cuda_malloc( len_fxc_buffer );
  }

  safe_zero( len_vxc_buffer, vxc, stream );
  safe_zero( len_fxc_buffer, fxc, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {
    if (use_inc) {
      kernels_[i].second.eval_vxc_fxc_inc_device(
        kernels_[i].first, N, rho, vxc, fxc, stream
      );
    } else {
      double* vxc_eval = i ? vxc_scr : vxc;
      double* fxc_eval = i ? fxc_scr : fxc;
      kernels_[i].second.eval_vxc_fxc_device(N, rho, vxc_eval, fxc_eval, stream);

      if( i ) {
        add_scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, stream );
        add_scal_device( len_fxc_buffer, kernels_[i].first, fxc_eval, fxc, stream );
      } else {
        scal_device( len_vxc_buffer, kernels_[i].first, vxc_eval, vxc, stream );
        scal_device( len_fxc_buffer, kernels_[i].first, fxc_eval, fxc, stream );
      }
    }
  }

  if( vxc_scr ) cudaFree( vxc_scr );
  if( fxc_scr ) cudaFree( fxc_scr );
}


// GGA Interfaces

GGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_cuda_malloc( len_exc_buffer );

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

  if( eps_scr ) cudaFree( eps_scr );

}


GGA_EXC_VXC_GENERATOR_DEVICE( XCFunctional::eval_exc_vxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  size_t len_exc_buffer    = exc_buffer_len(N);
  size_t len_vrho_buffer   = vrho_buffer_len(N);
  size_t len_vsigma_buffer = vsigma_buffer_len(N);

  double* eps_scr(nullptr), *vrho_scr(nullptr), *vsigma_scr(nullptr);
  if( kernels_.size() > 1 and not supports_inc_interface() ) {
    eps_scr    = safe_cuda_malloc( len_exc_buffer );
    vrho_scr   = safe_cuda_malloc( len_vrho_buffer );
    vsigma_scr = safe_cuda_malloc( len_vsigma_buffer );
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

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );

}

GGA_FXC_GENERATOR_DEVICE( XCFunctional::eval_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_v2rho2_buffer = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2sigma2_buffer = v2sigma2_buffer_len(N);

  double* v2rho2_scr(nullptr), *v2rhosigma_scr(nullptr), *v2sigma2_scr(nullptr);
  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) {
    v2rho2_scr = safe_cuda_malloc( len_v2rho2_buffer );
    v2rhosigma_scr = safe_cuda_malloc( len_v2rhosigma_buffer );
    v2sigma2_scr = safe_cuda_malloc( len_v2sigma2_buffer );
  }

  safe_zero( len_v2rho2_buffer, v2rho2, stream );
  safe_zero( len_v2rhosigma_buffer, v2rhosigma, stream );
  safe_zero( len_v2sigma2_buffer, v2sigma2, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {
    if (use_inc) {
      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2, stream
        );
      else
        kernels_[i].second.eval_fxc_inc_device(
          kernels_[i].first, N, rho, v2rho2, stream
        );
    } else {
      double* v2rho2_eval    = i ? v2rho2_scr    : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr : v2rhosigma;
      double* v2sigma2_eval  = i ? v2sigma2_scr  : v2sigma2;

      if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_device(N, rho, sigma, v2rho2_eval, 
          v2rhosigma_eval, v2sigma2_eval, stream );
      else
        kernels_[i].second.eval_fxc_device(N, rho, v2rho2_eval, stream);

      if( i ) {
        add_scal_device( len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream );
        if( kernels_[i].second.is_gga() ){
          add_scal_device( len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream );
          add_scal_device( len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream );
        }

      } else {
        scal_device( len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream );
        if( kernels_[i].second.is_gga() ){
          scal_device( len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream );
          scal_device( len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream );
        }
      }
    }
  }

  if( v2rho2_scr ) cudaFree( v2rho2_scr );
  if( v2rhosigma_scr ) cudaFree( v2rhosigma_scr );
  if( v2sigma2_scr ) cudaFree( v2sigma2_scr );
}

GGA_VXC_FXC_GENERATOR_DEVICE( XCFunctional::eval_vxc_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT GGA",  is_gga() );

  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_v2rho2_buffer = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);

  double* vrho_scr(nullptr), *vsigma_scr(nullptr); 
  double* v2rho2_scr(nullptr), *v2rhosigma_scr(nullptr), *v2sigma2_scr(nullptr);
  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) {
    vrho_scr = safe_cuda_malloc( len_vrho_buffer );
    vsigma_scr = safe_cuda_malloc( len_vsigma_buffer );
    v2rho2_scr = safe_cuda_malloc( len_v2rho2_buffer );
    v2rhosigma_scr = safe_cuda_malloc( len_v2rhosigma_buffer );
    v2sigma2_scr = safe_cuda_malloc( len_v2sigma2_buffer );
  }

  safe_zero( len_vrho_buffer, vrho, stream );
  safe_zero( len_vsigma_buffer, vsigma, stream );
  safe_zero( len_v2rho2_buffer, v2rho2, stream );
  safe_zero( len_v2rhosigma_buffer, v2rhosigma, stream );
  safe_zero( len_v2sigma2_buffer, v2sigma2, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {
    if (use_inc) {
      if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2, stream
        );
      } else {
        kernels_[i].second.eval_vxc_fxc_inc_device(
          kernels_[i].first, N, rho, vrho, v2rho2, stream
        );
      }
    } else {
      double* vrho_eval = i ? vrho_scr : vrho;
      double* vsigma_eval = i ? vsigma_scr : vsigma;
      double* v2rho2_eval = i ? v2rho2_scr : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr : v2rhosigma;
      double* v2sigma2_eval = i ? v2sigma2_scr : v2sigma2;

      if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc_device(
          N, rho, sigma, vrho_eval, vsigma_eval, v2rho2_eval, v2rhosigma_eval, v2sigma2_eval, stream);
      } else {
        kernels_[i].second.eval_vxc_fxc_device(N, rho, vrho_eval, v2rho2_eval, stream);
      }

      if (i) {
        add_scal_device(len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream);
        add_scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if (kernels_[i].second.is_gga()) {
          add_scal_device(len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream);
          add_scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          add_scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }
      } else {
        scal_device(len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream);
        scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if (kernels_[i].second.is_gga()) {
          scal_device(len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream);
          scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }
      }
    }
  }

  if( vrho_scr ) cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );
  if( v2rho2_scr ) cudaFree( v2rho2_scr );
  if( v2rhosigma_scr ) cudaFree( v2rhosigma_scr );
  if( v2sigma2_scr ) cudaFree( v2sigma2_scr );
}

// mGGA Interfaces

MGGA_EXC_GENERATOR_DEVICE( XCFunctional::eval_exc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  size_t len_exc_buffer = exc_buffer_len( N );

  double* eps_scr = nullptr;
  if( kernels_.size() > 1 and not supports_inc_interface() ) 
    eps_scr = safe_cuda_malloc( len_exc_buffer );

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

  if( eps_scr ) cudaFree( eps_scr );

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
    eps_scr    = safe_cuda_malloc( len_exc_buffer );
    vrho_scr   = safe_cuda_malloc( len_vrho_buffer );
    vsigma_scr = safe_cuda_malloc( len_vsigma_buffer );
    vtau_scr   = safe_cuda_malloc( len_vtau_buffer );
    if(needs_laplacian())
      vlapl_scr  = safe_cuda_malloc( len_vlapl_buffer );
  }

  safe_zero( len_exc_buffer, eps, stream );
  safe_zero( len_vrho_buffer, vrho, stream );
  safe_zero( len_vsigma_buffer, vsigma, stream );
  safe_zero( len_vtau_buffer, vtau, stream );
  
  if(needs_laplacian()) 
    safe_zero( len_vlapl_buffer, vlapl, stream );

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

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() ) {
          add_scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );
        }

        if( kernels_[i].second.is_mgga() ) {
          add_scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, stream  );
        }

        if( kernels_[i].second.needs_laplacian() ) {
          add_scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream );
        }

      } else {

        scal_device( len_exc_buffer, kernels_[i].first, eps_eval, eps, stream );
        scal_device( len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream );

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() ) {
          scal_device( len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream );
        }

        if( kernels_[i].second.is_mgga() ) {
          scal_device( len_vtau_buffer,  kernels_[i].first, vtau_eval,  vtau, stream  );
        }

        if( kernels_[i].second.needs_laplacian() ) {
          scal_device( len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream );
        }

      }
    }
  }

  if( eps_scr )    cudaFree( eps_scr );
  if( vrho_scr )   cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );
  if( vlapl_scr )  cudaFree( vlapl_scr );
  if( vtau_scr )   cudaFree( vtau_scr );
}

MGGA_FXC_GENERATOR_DEVICE( XCFunctional::eval_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_v2rho2_buffer    = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2rholapl_buffer  = v2rholapl_buffer_len(N);
  const size_t len_v2rhotau_buffer   = v2rhotau_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);
  const size_t len_v2sigmalapl_buffer = v2sigmalapl_buffer_len(N);
  const size_t len_v2sigmatau_buffer  = v2sigmatau_buffer_len(N);
  const size_t len_v2lapl2_buffer    = v2lapl2_buffer_len(N);
  const size_t len_v2lapltau_buffer  = v2lapltau_buffer_len(N);
  const size_t len_v2tau2_buffer     = v2tau2_buffer_len(N);

  double* v2rho2_scr(nullptr), *v2rhosigma_scr(nullptr), *v2rholapl_scr(nullptr), *v2rhotau_scr(nullptr),
    *v2sigma2_scr(nullptr), *v2sigmalapl_scr(nullptr), *v2sigmatau_scr(nullptr), *v2lapl2_scr(nullptr), 
    *v2lapltau_scr(nullptr), *v2tau2_scr(nullptr);

  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) {
    v2rho2_scr = safe_cuda_malloc( len_v2rho2_buffer );
    v2rhosigma_scr = safe_cuda_malloc( len_v2rhosigma_buffer );
    v2rholapl_scr = safe_cuda_malloc( len_v2rholapl_buffer );
    v2rhotau_scr = safe_cuda_malloc( len_v2rhotau_buffer );
    v2sigma2_scr = safe_cuda_malloc( len_v2sigma2_buffer );
    v2sigmalapl_scr = safe_cuda_malloc( len_v2sigmalapl_buffer );
    v2sigmatau_scr = safe_cuda_malloc( len_v2sigmatau_buffer );
    v2lapl2_scr = safe_cuda_malloc( len_v2lapl2_buffer );
    v2lapltau_scr = safe_cuda_malloc( len_v2lapltau_buffer );
    v2tau2_scr = safe_cuda_malloc( len_v2tau2_buffer );
  }

  safe_zero( len_v2rho2_buffer, v2rho2, stream );
  safe_zero( len_v2rhosigma_buffer, v2rhosigma, stream );
  safe_zero( len_v2rholapl_buffer, v2rholapl, stream );
  safe_zero( len_v2rhotau_buffer, v2rhotau, stream );
  safe_zero( len_v2sigma2_buffer, v2sigma2, stream );
  safe_zero( len_v2sigmalapl_buffer, v2sigmalapl, stream );
  safe_zero( len_v2sigmatau_buffer, v2sigmatau, stream );
  safe_zero( len_v2lapl2_buffer, v2lapl2, stream );
  safe_zero( len_v2lapltau_buffer, v2lapltau, stream );
  safe_zero( len_v2tau2_buffer, v2tau2, stream );

  for( auto i = 0ul; i < kernels_.size(); ++i ) {

    if( use_inc ) {
      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, v2rho2, v2rhosigma, v2rholapl, v2rhotau,
          v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2, v2lapltau, v2tau2, stream
        );
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, v2rho2, v2rhosigma, v2sigma2, stream
        );
      else
        kernels_[i].second.eval_fxc_inc_device(
          kernels_[i].first, N, rho, v2rho2, stream
        );
    } else {
      double* v2rho2_eval    = i ? v2rho2_scr    : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr : v2rhosigma;
      double* v2rholapl_eval  = i ? v2rholapl_scr  : v2rholapl;
      double* v2rhotau_eval   = i ? v2rhotau_scr   : v2rhotau;
      double* v2sigma2_eval   = i ? v2sigma2_scr   : v2sigma2;
      double* v2sigmalapl_eval = i ? v2sigmalapl_scr : v2sigmalapl;
      double* v2sigmatau_eval  = i ? v2sigmatau_scr  : v2sigmatau;
      double* v2lapl2_eval    = i ? v2lapl2_scr    : v2lapl2;
      double* v2lapltau_eval  = i ? v2lapltau_scr  : v2lapltau;
      double* v2tau2_eval     = i ? v2tau2_scr     : v2tau2;

      if( kernels_[i].second.is_mgga() )
        kernels_[i].second.eval_fxc_device(N, rho, sigma, lapl, tau, v2rho2_eval, 
          v2rhosigma_eval, v2rholapl_eval, v2rhotau_eval, v2sigma2_eval, v2sigmalapl_eval, 
          v2sigmatau_eval, v2lapl2_eval, v2lapltau_eval, v2tau2_eval, stream);
      else if( kernels_[i].second.is_gga() )
        kernels_[i].second.eval_fxc_device(N, rho, sigma, v2rho2_eval, v2rhosigma_eval, v2sigma2_eval, stream);
      else
        kernels_[i].second.eval_fxc_device(N, rho, v2rho2_eval, stream);

      if (i) {
        add_scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if( kernels_[i].second.is_gga() or kernels_[i].second.is_mgga() ){
          add_scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          add_scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }

        if( kernels_[i].second.needs_laplacian() ) {
          add_scal_device(len_v2rholapl_buffer, kernels_[i].first, v2rholapl_eval, v2rholapl, stream);
          add_scal_device(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl_eval, v2sigmalapl, stream);
          add_scal_device(len_v2lapl2_buffer, kernels_[i].first, v2lapl2_eval, v2lapl2, stream);
        }
        
        if( kernels_[i].second.is_mgga() ) {
          add_scal_device(len_v2rhotau_buffer, kernels_[i].first, v2rhotau_eval, v2rhotau, stream);
          add_scal_device(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau_eval, v2sigmatau, stream);
          add_scal_device(len_v2tau2_buffer, kernels_[i].first, v2tau2_eval, v2tau2, stream);
        }
        
        if ( kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga() ) {
          add_scal_device(len_v2lapltau_buffer, kernels_[i].first, v2lapltau_eval, v2lapltau, stream);
        }

      } else{

        scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if (kernels_[i].second.is_gga() or kernels_[i].second.is_mgga()) {
          scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }

        if (kernels_[i].second.needs_laplacian()) {
          scal_device(len_v2rholapl_buffer, kernels_[i].first, v2rholapl_eval, v2rholapl, stream);
          scal_device(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl_eval, v2sigmalapl, stream);
          scal_device(len_v2lapl2_buffer, kernels_[i].first, v2lapl2_eval, v2lapl2, stream);
        }

        if (kernels_[i].second.is_mgga()) {
          scal_device(len_v2rhotau_buffer, kernels_[i].first, v2rhotau_eval, v2rhotau, stream);
          scal_device(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau_eval, v2sigmatau, stream);
          scal_device(len_v2tau2_buffer, kernels_[i].first, v2tau2_eval, v2tau2, stream);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          scal_device(len_v2lapltau_buffer, kernels_[i].first, v2lapltau_eval, v2lapltau, stream);
        }
      }
    }
  }

  if( v2rho2_scr ) cudaFree( v2rho2_scr );
  if( v2rhosigma_scr ) cudaFree( v2rhosigma_scr );
  if( v2rholapl_scr ) cudaFree( v2rholapl_scr );
  if( v2rhotau_scr ) cudaFree( v2rhotau_scr );
  if( v2sigma2_scr ) cudaFree( v2sigma2_scr );
  if( v2sigmalapl_scr ) cudaFree( v2sigmalapl_scr );
  if( v2sigmatau_scr ) cudaFree( v2sigmatau_scr );
  if( v2lapl2_scr ) cudaFree( v2lapl2_scr );
  if( v2lapltau_scr ) cudaFree( v2lapltau_scr );
  if( v2tau2_scr ) cudaFree( v2tau2_scr );
}

MGGA_VXC_FXC_GENERATOR_DEVICE( XCFunctional::eval_vxc_fxc_device ) const {

  throw_if_not_sane();
  EXCHCXX_BOOL_CHECK("KERNEL IS NOT MGGA",  is_mgga() );

  const size_t len_vrho_buffer   = vrho_buffer_len(N);
  const size_t len_vsigma_buffer = vsigma_buffer_len(N);
  const size_t len_vlapl_buffer  = vlapl_buffer_len(N);
  const size_t len_vtau_buffer   = vtau_buffer_len(N);
  const size_t len_v2rho2_buffer    = v2rho2_buffer_len(N);
  const size_t len_v2rhosigma_buffer = v2rhosigma_buffer_len(N);
  const size_t len_v2rholapl_buffer  = v2rholapl_buffer_len(N);
  const size_t len_v2rhotau_buffer   = v2rhotau_buffer_len(N);
  const size_t len_v2sigma2_buffer   = v2sigma2_buffer_len(N);
  const size_t len_v2sigmalapl_buffer = v2sigmalapl_buffer_len(N);
  const size_t len_v2sigmatau_buffer  = v2sigmatau_buffer_len(N);
  const size_t len_v2lapl2_buffer    = v2lapl2_buffer_len(N);
  const size_t len_v2lapltau_buffer  = v2lapltau_buffer_len(N);
  const size_t len_v2tau2_buffer     = v2tau2_buffer_len(N);

  double* vrho_scr(nullptr), *vsigma_scr(nullptr), *vlapl_scr(nullptr), *vtau_scr(nullptr);
  double* v2rho2_scr(nullptr), *v2rhosigma_scr(nullptr), *v2rholapl_scr(nullptr), *v2rhotau_scr(nullptr),
    *v2sigma2_scr(nullptr), *v2sigmalapl_scr(nullptr), *v2sigmatau_scr(nullptr), *v2lapl2_scr(nullptr), 
    *v2lapltau_scr(nullptr), *v2tau2_scr(nullptr);

  bool use_inc = supports_inc_interface();
  if( kernels_.size() > 1 && !use_inc ) {
    vrho_scr = safe_cuda_malloc( len_vrho_buffer );
    vsigma_scr = safe_cuda_malloc( len_vsigma_buffer );
    vlapl_scr = safe_cuda_malloc( len_vlapl_buffer );
    vtau_scr = safe_cuda_malloc( len_vtau_buffer );
    v2rho2_scr = safe_cuda_malloc( len_v2rho2_buffer );
    v2rhosigma_scr = safe_cuda_malloc( len_v2rhosigma_buffer );
    v2rholapl_scr = safe_cuda_malloc( len_v2rholapl_buffer );
    v2rhotau_scr = safe_cuda_malloc(len_v2rhotau_buffer);
    v2sigma2_scr = safe_cuda_malloc(len_v2sigma2_buffer);
    v2sigmalapl_scr = safe_cuda_malloc(len_v2sigmalapl_buffer);
    v2sigmatau_scr = safe_cuda_malloc(len_v2sigmatau_buffer);
    v2lapl2_scr = safe_cuda_malloc(len_v2lapl2_buffer);
    v2lapltau_scr = safe_cuda_malloc(len_v2lapltau_buffer);
    v2tau2_scr = safe_cuda_malloc(len_v2tau2_buffer);
  }

  safe_zero(len_vrho_buffer, vrho, stream);
  safe_zero(len_vsigma_buffer, vsigma, stream);
  safe_zero(len_vlapl_buffer, vlapl, stream);
  safe_zero(len_vtau_buffer, vtau, stream);
  safe_zero(len_v2rho2_buffer, v2rho2, stream);
  safe_zero(len_v2rhosigma_buffer, v2rhosigma, stream);
  safe_zero(len_v2rholapl_buffer, v2rholapl, stream);
  safe_zero(len_v2rhotau_buffer, v2rhotau, stream);
  safe_zero(len_v2sigma2_buffer, v2sigma2, stream);
  safe_zero(len_v2sigmalapl_buffer, v2sigmalapl, stream);
  safe_zero(len_v2sigmatau_buffer, v2sigmatau, stream);
  safe_zero(len_v2lapl2_buffer, v2lapl2, stream);
  safe_zero(len_v2lapltau_buffer, v2lapltau, stream);
  safe_zero(len_v2tau2_buffer, v2tau2, stream);

  for (auto i = 0ul; i < kernels_.size(); ++i) {
    if( use_inc ) {
      if (kernels_[i].second.is_mgga()) {
        kernels_[i].second.eval_vxc_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, lapl, tau, vrho, vsigma, vlapl, vtau,
          v2rho2, v2rhosigma, v2rholapl, v2rhotau,
          v2sigma2, v2sigmalapl, v2sigmatau, v2lapl2,
          v2lapltau, v2tau2, stream);
      } else if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc_inc_device(
          kernels_[i].first, N, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma,
          v2sigma2, stream);
      } else {
        kernels_[i].second.eval_vxc_fxc_inc_device(
          kernels_[i].first, N, rho, vrho, v2rho2, stream);
      }
    } else {
      double* vrho_eval = i ? vrho_scr : vrho;
      double* vsigma_eval = i ? vsigma_scr : vsigma;
      double* vlapl_eval = i ? vlapl_scr : vlapl;
      double* vtau_eval = i ? vtau_scr : vtau;
      double* v2rho2_eval = i ? v2rho2_scr : v2rho2;
      double* v2rhosigma_eval = i ? v2rhosigma_scr : v2rhosigma;
      double* v2rholapl_eval = i ? v2rholapl_scr : v2rholapl;
      double* v2rhotau_eval = i ? v2rhotau_scr : v2rhotau;
      double* v2sigma2_eval = i ? v2sigma2_scr : v2sigma2;
      double* v2sigmalapl_eval = i ? v2sigmalapl_scr : v2sigmalapl;
      double* v2sigmatau_eval = i ? v2sigmatau_scr : v2sigmatau;
      double* v2lapl2_eval = i ? v2lapl2_scr : v2lapl2;
      double* v2lapltau_eval = i ? v2lapltau_scr : v2lapltau;
      double* v2tau2_eval = i ? v2tau2_scr : v2tau2;

      if (kernels_[i].second.is_mgga()) {
        kernels_[i].second.eval_vxc_fxc_device(
          N, rho, sigma, lapl, tau, vrho_eval, vsigma_eval, vlapl_eval, vtau_eval,
          v2rho2_eval, v2rhosigma_eval, v2rholapl_eval, v2rhotau_eval,
          v2sigma2_eval, v2sigmalapl_eval, v2sigmatau_eval, v2lapl2_eval,
          v2lapltau_eval, v2tau2_eval, stream);
      } else if (kernels_[i].second.is_gga()) {
        kernels_[i].second.eval_vxc_fxc_device(
          N, rho, sigma, vrho_eval, vsigma_eval, v2rho2_eval, v2rhosigma_eval,
          v2sigma2_eval, stream);
      } else {
        kernels_[i].second.eval_vxc_fxc_device(N, rho, vrho_eval, v2rho2_eval, stream);
      }

      if (i) {
        add_scal_device(len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream);
        add_scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if (kernels_[i].second.is_gga() || kernels_[i].second.is_mgga()) {
          add_scal_device(len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream);
          add_scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          add_scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }

        if (kernels_[i].second.needs_laplacian()) {
          add_scal_device(len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream);
          add_scal_device(len_v2rholapl_buffer, kernels_[i].first, v2rholapl_eval, v2rholapl, stream);
          add_scal_device(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl_eval, v2sigmalapl, stream);
          add_scal_device(len_v2lapl2_buffer, kernels_[i].first, v2lapl2_eval, v2lapl2, stream);
        }

        if (kernels_[i].second.is_mgga()) {
          add_scal_device(len_vtau_buffer, kernels_[i].first, vtau_eval, vtau, stream);
          add_scal_device(len_v2rhotau_buffer, kernels_[i].first, v2rhotau_eval, v2rhotau, stream);
          add_scal_device(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau_eval, v2sigmatau, stream);
          add_scal_device(len_v2tau2_buffer, kernels_[i].first, v2tau2_eval, v2tau2, stream);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          add_scal_device(len_v2lapltau_buffer, kernels_[i].first, v2lapltau_eval, v2lapltau, stream);
        }
      } else {
        scal_device(len_vrho_buffer, kernels_[i].first, vrho_eval, vrho, stream);
        scal_device(len_v2rho2_buffer, kernels_[i].first, v2rho2_eval, v2rho2, stream);

        if (kernels_[i].second.is_gga() || kernels_[i].second.is_mgga()) {
          scal_device(len_vsigma_buffer, kernels_[i].first, vsigma_eval, vsigma, stream);
          scal_device(len_v2rhosigma_buffer, kernels_[i].first, v2rhosigma_eval, v2rhosigma, stream);
          scal_device(len_v2sigma2_buffer, kernels_[i].first, v2sigma2_eval, v2sigma2, stream);
        }

        if (kernels_[i].second.needs_laplacian()) {
          scal_device(len_vlapl_buffer, kernels_[i].first, vlapl_eval, vlapl, stream);
          scal_device(len_v2rholapl_buffer, kernels_[i].first, v2rholapl_eval, v2rholapl, stream);
          scal_device(len_v2sigmalapl_buffer, kernels_[i].first, v2sigmalapl_eval, v2sigmalapl, stream);
          scal_device(len_v2lapl2_buffer, kernels_[i].first, v2lapl2_eval, v2lapl2, stream);
        }

        if (kernels_[i].second.is_mgga()) {
          scal_device(len_vtau_buffer, kernels_[i].first, vtau_eval, vtau, stream);
          scal_device(len_v2rhotau_buffer, kernels_[i].first, v2rhotau_eval, v2rhotau, stream);
          scal_device(len_v2sigmatau_buffer, kernels_[i].first, v2sigmatau_eval, v2sigmatau, stream);
          scal_device(len_v2tau2_buffer, kernels_[i].first, v2tau2_eval, v2tau2, stream);
        }

        if (kernels_[i].second.needs_laplacian() && kernels_[i].second.is_mgga()) {
          scal_device(len_v2lapltau_buffer, kernels_[i].first, v2lapltau_eval, v2lapltau, stream);
        }
      }
    }
  }

  if( vrho_scr ) cudaFree( vrho_scr );
  if( vsigma_scr ) cudaFree( vsigma_scr );
  if( vlapl_scr ) cudaFree( vlapl_scr );
  if( vtau_scr ) cudaFree( vtau_scr );
  if( v2rho2_scr ) cudaFree( v2rho2_scr );
  if( v2rhosigma_scr ) cudaFree( v2rhosigma_scr );
  if( v2rholapl_scr ) cudaFree( v2rholapl_scr );
  if( v2rhotau_scr ) cudaFree( v2rhotau_scr );
  if( v2sigma2_scr ) cudaFree( v2sigma2_scr );
  if( v2sigmalapl_scr ) cudaFree( v2sigmalapl_scr );
  if( v2sigmatau_scr ) cudaFree( v2sigmatau_scr );
  if( v2lapl2_scr ) cudaFree( v2lapl2_scr );
  if( v2lapltau_scr ) cudaFree( v2lapltau_scr );
  if( v2tau2_scr ) cudaFree( v2tau2_scr );
}


}
