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

#pragma once

#include <exchcxx/exchcxx_config.hpp>

#include <cmath>

namespace ExchCXX {


#if defined(__CUDACC__) || defined(__HIPCC__)

#define SAFE_CONSTEXPR_INLINE(TYPE) static inline constexpr TYPE __host__ __device__
#define SAFE_INLINE(TYPE)           static inline           TYPE __host__ __device__

template <typename F>
SAFE_INLINE(F) safe_max(F x, F y) { return fmax(x,y); }
template <typename F>
SAFE_INLINE(F) safe_min(F x, F y) { return fmin(x,y); }

#elif defined(__SYCL_DEVICE_ONLY__) && defined(EXCHCXX_ENABLE_SYCL)

#define SAFE_CONSTEXPR_INLINE(TYPE) static __attribute__((always_inline)) constexpr TYPE
#define SAFE_INLINE(TYPE)           static __attribute__((always_inline)) TYPE

template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_max(F x, F y) { return sycl::fmax(x,y); }
template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_min(F x, F y) { return sycl::fmin(x,y); }

#else

#define SAFE_CONSTEXPR_INLINE(TYPE) static inline constexpr TYPE
#define SAFE_INLINE(TYPE)           static inline           TYPE

template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_max(F x, F y) { return std::max(x,y); }
template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_min(F x, F y) { return std::min(x,y); }

#endif

#define BUILTIN_KERNEL_EVAL_RETURN SAFE_INLINE(void)

namespace safe_math {

#ifdef EXCHCXX_ENABLE_SYCL
namespace sm = sycl;
#else
namespace sm = std;
#endif

template <typename T>
SAFE_INLINE(auto) sqrt( T x ) { return sm::sqrt(x); }
template <typename T>
SAFE_INLINE(auto) cbrt( T x ) { return sm::cbrt(x); }
template <typename T>
SAFE_INLINE(auto) log( T x ) { return sm::log(x); }
template <typename T>
SAFE_INLINE(auto) exp( T x ) { return sm::exp(x); }
template <typename T>
SAFE_INLINE(auto) atan( T x ) { return sm::atan(x); }
template <typename T, typename U>
SAFE_INLINE(auto) pow( T x, U e ) { return sm::pow(x,e); }

}

template <typename KernelType, typename... Args>
struct max_dens_tol {
  static constexpr double dens_tol = std::max(
    KernelType::dens_tol, max_dens_tol<Args...>::dens_tol
  );
};

template <typename KernelType>
struct max_dens_tol<KernelType> {
  static constexpr double dens_tol = KernelType::dens_tol;
};


template <typename F>
SAFE_CONSTEXPR_INLINE( F ) square( F x ) { return x*x; }
template <typename F>
SAFE_INLINE( F ) pow_3_2( F x ) { F y = safe_math::sqrt(x); return y*y*y; }

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_3( bool b, F x, F y ) {
  return b ? x : y;
}

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_5( bool b, F x, bool c, F y, F z ) {
  return b ? x : (c ? y : z);
}

}
