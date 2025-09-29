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

#pragma once

#include <exchcxx/exchcxx_config.hpp>

#include <cmath>
#include <cfloat>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define EXCHCXX_READONLY_TABLE  static __device__
#elif defined(__SYCL_DEVICE_ONLY__)
#define EXCHCXX_READONLY_TABLE  inline constexpr
#else
#define EXCHCXX_READONLY_TABLE  static
#endif

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

// Apparently C / C++ cbrt differ ever-so-slightly
template <typename T>
SAFE_INLINE(double) cbrt( T x ) { return ::cbrt(x); }
template <typename T>
SAFE_INLINE(auto) log( T x ) { return sm::log(x); }
template <typename T>
SAFE_INLINE(auto) exp( T x ) { return sm::exp(x); }
template <typename T>
SAFE_INLINE(auto) atan( T x ) { return sm::atan(x); }
template <typename T>
SAFE_INLINE(auto) erf( T x ) { return sm::erf(x); }
template <typename T, typename U>
SAFE_INLINE(auto) pow( T x, U e ) { return sm::pow(x,e); }
template <typename T>
SAFE_INLINE(auto) xc_erfcx( T x ) { return sm::exp(x*x)*sm::erfc(x); } 



// This function is taken from libxc
SAFE_INLINE(auto) xc_cheb_eval(const double x, const double *cs, const int N)
{
  int i;
  double twox, b0, b1, b2;

  b2 = b1 = b0 = 0.0;

  twox = 2.0*x;
  for(i=N-1; i>=0; i--){
    b2 = b1;
    b1 = b0;
    b0 = twox*b1 - b2 + cs[i];
  }

  return 0.5*(b0 - b2);
}
// The following data is taken from libxc
EXCHCXX_READONLY_TABLE double AE11_data[39] = {
   0.121503239716065790, -0.065088778513550150,  0.004897651357459670, -0.000649237843027216,  0.000093840434587471,
   0.000000420236380882, -0.000008113374735904,  0.000002804247688663,  0.000000056487164441, -0.000000344809174450,
   0.000000058209273578,  0.000000038711426349, -0.000000012453235014, -0.000000005118504888,  0.000000002148771527,
   0.000000000868459898, -0.000000000343650105, -0.000000000179796603,  0.000000000047442060,  0.000000000040423282,
  -0.000000000003543928, -0.000000000008853444, -0.000000000000960151,  0.000000000001692921,  0.000000000000607990,
  -0.000000000000224338, -0.000000000000200327, -0.000000000000006246,  0.000000000000045571,  0.000000000000016383,
  -0.000000000000005561, -0.000000000000006074, -0.000000000000000862,  0.000000000000001223,  0.000000000000000716,
  -0.000000000000000024, -0.000000000000000201, -0.000000000000000082,  0.000000000000000017
};

EXCHCXX_READONLY_TABLE double AE12_data[25] = {
   0.582417495134726740, -0.158348850905782750, -0.006764275590323141,  0.005125843950185725,  0.000435232492169391,
  -0.000143613366305483, -0.000041801320556301, -0.000002713395758640,  0.000001151381913647,  0.000000420650022012,
   0.000000066581901391,  0.000000000662143777, -0.000000002844104870, -0.000000000940724197, -0.000000000177476602,
  -0.000000000015830222,  0.000000000002905732,  0.000000000001769356,  0.000000000000492735,  0.000000000000093709,
   0.000000000000010707, -0.000000000000000537, -0.000000000000000716, -0.000000000000000244, -0.000000000000000058
};

EXCHCXX_READONLY_TABLE double E11_data[19] = {
  -16.11346165557149402600,   7.79407277874268027690,  -1.95540581886314195070,   0.37337293866277945612,  -0.05692503191092901938,
    0.00721107776966009185,  -0.00078104901449841593,   0.00007388093356262168,  -0.00000620286187580820,   0.00000046816002303176,
   -0.00000003209288853329,   0.00000000201519974874,  -0.00000000011673686816,   0.00000000000627627066,  -0.00000000000031481541,
    0.00000000000001479904,  -0.00000000000000065457,   0.00000000000000002733,  -0.00000000000000000108
};

EXCHCXX_READONLY_TABLE double E12_data[16] = {
  -0.03739021479220279500,  0.04272398606220957700, -0.13031820798497005440,  0.01441912402469889073, -0.00134617078051068022,
   0.00010731029253063780, -0.00000742999951611943,  0.00000045377325690753, -0.00000002476417211390,  0.00000000122076581374,
  -0.00000000005485141480,  0.00000000000226362142, -0.00000000000008635897,  0.00000000000000306291, -0.00000000000000010148,
   0.00000000000000000315
};

EXCHCXX_READONLY_TABLE double AE13_data[25] = {
  -0.605773246640603460, -0.112535243483660900,  0.013432266247902779, -0.001926845187381145,  0.000309118337720603,
  -0.000053564132129618,  0.000009827812880247, -0.000001885368984916,  0.000000374943193568, -0.000000076823455870,
   0.000000016143270567, -0.000000003466802211,  0.000000000758754209, -0.000000000168864333,  0.000000000038145706,
  -0.000000000008733026,  0.000000000002023672, -0.000000000000474132,  0.000000000000112211, -0.000000000000026804,
   0.000000000000006457, -0.000000000000001568,  0.000000000000000383, -0.000000000000000094,  0.000000000000000023
};

EXCHCXX_READONLY_TABLE double AE14_data[26] = {
  -0.18929180007530170, -0.08648117855259871,  0.00722410154374659, -0.00080975594575573,  0.00010999134432661,
  -0.00001717332998937,  0.00000298562751447, -0.00000056596491457,  0.00000011526808397, -0.00000002495030440,
   0.00000000569232420, -0.00000000135995766,  0.00000000033846628, -0.00000000008737853,  0.00000000002331588,
  -0.00000000000641148,  0.00000000000181224, -0.00000000000052538,  0.00000000000015592, -0.00000000000004729,
   0.00000000000001463, -0.00000000000000461,  0.00000000000000148, -0.00000000000000048,  0.00000000000000016,
  -0.00000000000000005
};


// This function is taken from libxc
/* implementation for E1, allowing for scaling by exp(x) */
template <typename T>
SAFE_INLINE(auto) xc_E1_scaled(T x){
  int scale = 1;
  const double xmaxt = -log(DBL_MIN);        /* XMAXT = -log (R1MACH(1)) */
  const double xmax  = xmaxt - log(xmaxt);    /* XMAX = XMAXT - log(XMAXT) */

  double e1 = 0.0;

  /* this is a workaround not to have argument errors */
  if(! scale) x = safe_min(x, xmax);

  if(x <= -10.0){
    const double s = 1.0/x;
    e1 = s * (1.0 + xc_cheb_eval(20.0/x + 1.0, AE11_data, 39));
  }else if(x <= -4.0){
    const double s = 1.0/x;
    e1 = s * (1.0 + xc_cheb_eval((40.0/x + 7.0)/3.0, AE12_data, 25));
  }else if(x <= -1.0){
    const double scale_factor = exp(x);
    e1 = scale_factor * (-log(fabs(x)) + xc_cheb_eval((2.0*x + 5.0)/3.0, E11_data, 19));
  }else if(x == 0.0) {
#if defined(__CUDACC__) || defined(__HIPCC__)
    printf("Argument cannot be 0.0 in expint_e1\n");
#endif
  }else if(x <= 1.0){
    const double scale_factor = exp(x);
    e1 = scale_factor*(-log(fabs(x)) - 0.6875 + x + xc_cheb_eval(x, E12_data, 16));
  }else if(x <= 4.0){
    const double s = 1.0/x;
    e1 = s * (1.0 + xc_cheb_eval((8.0/x - 5.0)/3.0, AE13_data, 25));
  }else if(x <= xmax || scale){
    const double s = 1.0/x;
    e1 = s * (1.0 + xc_cheb_eval(8.0/x - 1.0, AE14_data, 26));
  }else{
#if defined(__CUDACC__) || defined(__HIPCC__)
  printf("Argument %14.10le is larger than xmax=%14.10le in expint_e1\n", x, xmax);
#endif
  }

  return e1;
}


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
SAFE_INLINE( F ) pow_1_4( F x ) { F y = safe_math::sqrt(x); return safe_math::sqrt(y); }

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_3( bool b, F x, F y ) {
  return b ? x : y;
}

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) piecewise_functor_5( bool b, F x, bool c, F y, F z ) {
  return b ? x : (c ? y : z);
}



template <typename F>
SAFE_CONSTEXPR_INLINE( F ) enforce_fermi_hole_curvature(F sigma, F rho, F tau) {
  return safe_min(sigma, F(8) * rho * tau);
}

template <typename F>
SAFE_CONSTEXPR_INLINE( F ) enforce_polar_sigma_constraints(F sigma_aa, F sigma_ab, F sigma_bb) {
  const auto s_ave = 0.5 * (sigma_aa + sigma_bb);
  sigma_ab = (sigma_ab >= -s_ave ? sigma_ab : -s_ave);
  sigma_ab = (sigma_ab <=  s_ave ? sigma_ab :  s_ave);
  return sigma_ab;
}


}
