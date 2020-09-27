#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif

#include <cmath>

namespace ExchCXX {


#ifdef __CUDACC__

#define SAFE_CONSTEXPR_INLINE(TYPE) static inline constexpr TYPE __host__ __device__
#define SAFE_INLINE(TYPE)           static inline           TYPE __host__ __device__

template <typename F>
SAFE_INLINE(F) safe_max(F x, F y) { return fmax(x,y); }
template <typename F>
SAFE_INLINE(F) safe_min(F x, F y) { return fmin(x,y); }

#else

#define SAFE_CONSTEXPR_INLINE(TYPE) static inline constexpr TYPE
#define SAFE_INLINE(TYPE)           static inline           TYPE

template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_max(F x, F y) { return std::max(x,y); }
template <typename F>
SAFE_CONSTEXPR_INLINE(F) safe_min(F x, F y) { return std::min(x,y); }

#endif

#define BUILTIN_KERNEL_EVAL_RETURN SAFE_INLINE(void)

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
SAFE_INLINE( F ) pow_3_2( F x ) { F y = std::sqrt(x); return y*y*y; }

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_3( bool b, F x, F y ) {
  return b ? x : y;
}

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_5( bool b, F x, bool c, F y, F z ) {
  return b ? x : (c ? y : z);
}

}
