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

namespace ExchCXX {
namespace detail {

template <typename KernelType>
inline LDA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  const double rho_use = sycl::max( rho[idx], 0. );
  traits::eval_exc_unpolar( rho_use, eps[idx] );

}

template <typename KernelType>
inline LDA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto rho_i = rho + 2*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );

  traits::eval_exc_polar( rho_a_use, rho_b_use, eps[idx] );

}

template <typename KernelType>
inline LDA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  const double rho_use = sycl::max( rho[idx], 0. );
  traits::eval_exc_vxc_unpolar( rho_use, eps[idx], vxc[idx] );

}

template <typename KernelType>
inline LDA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto rho_i = rho + 2*idx;
  auto vxc_i = vxc + 2*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );

  traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, eps[idx],
    vxc_i[0], vxc_i[1] );

}

template <typename KernelType>
inline LDA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e;

  const double rho_use = sycl::max( rho[idx], 0. );
  traits::eval_exc_unpolar( rho_use, e );
  eps[idx] += scal_fact * e;

}

template <typename KernelType>
inline LDA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto rho_i = rho + 2*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );

  double e;
  traits::eval_exc_polar( rho_a_use, rho_b_use, e );

  eps[idx] += scal_fact * e;

}

template <typename KernelType>
inline LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e,v;

  const double rho_use = sycl::max( rho[idx], 0. );
  traits::eval_exc_vxc_unpolar( rho_use, e, v );
  eps[idx] += scal_fact * e;
  vxc[idx] += scal_fact * v;

}

template <typename KernelType>
inline LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto rho_i = rho + 2*idx;
  auto vxc_i = vxc + 2*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );

  double v_a, v_b, e;
  traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, e, v_a, v_b);
  eps[idx] += scal_fact * e;
  vxc_i[0] += scal_fact * v_a;
  vxc_i[1] += scal_fact * v_b;

}

template <typename KernelType>
inline GGA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  const double rho_use   = sycl::max( rho[idx],   0.    );
  const double sigma_use = sycl::max( sigma[idx], 1e-40 );
  traits::eval_exc_unpolar( rho_use, sigma_use, eps[idx] );

}

template <typename KernelType>
inline GGA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto* rho_i   = rho   + 2*idx;
  auto* sigma_i = sigma + 3*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );
  const double sigma_aa_use = sycl::max( sigma_i[0], 1e-40 );
  const double sigma_bb_use = sycl::max( sigma_i[2], 1e-40 );
  const double sigma_ab_use = sycl::max(
    sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
  );

  traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use,
    sigma_ab_use, sigma_bb_use, eps[idx] );

}

template <typename KernelType>
inline GGA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  const double rho_use   = sycl::max( rho[idx],   0.    );
  const double sigma_use = sycl::max( sigma[idx], 1e-40 );
  traits::eval_exc_vxc_unpolar( rho_use, sigma_use, eps[idx],
    vrho[idx], vsigma[idx] );

}

template <typename KernelType>
inline GGA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto* rho_i    = rho   + 2*idx;
  auto* sigma_i  = sigma + 3*idx;
  auto* vrho_i   = vrho   + 2*idx;
  auto* vsigma_i = vsigma + 3*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );
  const double sigma_aa_use = sycl::max( sigma_i[0], 1e-40 );
  const double sigma_bb_use = sycl::max( sigma_i[2], 1e-40 );
  const double sigma_ab_use = sycl::max(
    sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
  );


  traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, sigma_aa_use,
    sigma_ab_use, sigma_bb_use, eps[idx], vrho_i[0], vrho_i[1],
    vsigma_i[0], vsigma_i[1], vsigma_i[2] );

}


template <typename KernelType>
inline GGA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e;
  const double rho_use   = sycl::max( rho[idx],   0.    );
  const double sigma_use = sycl::max( sigma[idx], 1e-40 );

  traits::eval_exc_unpolar( rho_use, sigma_use, e );
  eps[idx] += scal_fact * e;

}

template <typename KernelType>
inline GGA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto* rho_i   = rho   + 2*idx;
  auto* sigma_i = sigma + 3*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );
  const double sigma_aa_use = sycl::max( sigma_i[0], 1e-40 );
  const double sigma_bb_use = sycl::max( sigma_i[2], 1e-40 );
  const double sigma_ab_use = sycl::max(
    sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
  );

  double e;
  traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use,
    sigma_ab_use, sigma_bb_use, e );
  eps[idx] += scal_fact * e;

}

template <typename KernelType>
inline GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e, vr, vs;
  const double rho_use   = sycl::max( rho[idx],   0.    );
  const double sigma_use = sycl::max( sigma[idx], 1e-40 );

  traits::eval_exc_vxc_unpolar( rho_use, sigma_use, e, vr, vs );
  eps[idx]    += scal_fact * e;
  vrho[idx]   += scal_fact * vr;
  vsigma[idx] += scal_fact * vs;

}

template <typename KernelType>
inline GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  auto* rho_i    = rho   + 2*idx;
  auto* sigma_i  = sigma + 3*idx;
  auto* vrho_i   = vrho   + 2*idx;
  auto* vsigma_i = vsigma + 3*idx;

  const double rho_a_use = sycl::max( rho_i[0], 0. );
  const double rho_b_use = sycl::max( rho_i[1], 0. );
  const double sigma_aa_use = sycl::max( sigma_i[0], 1e-40 );
  const double sigma_bb_use = sycl::max( sigma_i[2], 1e-40 );
  const double sigma_ab_use = sycl::max(
    sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
  );


  double e, vra, vrb, vsaa,vsab,vsbb;
  traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, sigma_aa_use,
    sigma_ab_use, sigma_bb_use, e, vra, vrb, vsaa, vsab, vsbb );

  eps[idx]    += scal_fact * e;
  vrho_i[0]   += scal_fact * vra;
  vrho_i[1]   += scal_fact * vrb;
  vsigma_i[0] += scal_fact * vsaa;
  vsigma_i[1] += scal_fact * vsab;
  vsigma_i[2] += scal_fact * vsbb;

}






template <typename KernelType> class lda_eval_exc_unpolar;
template <typename KernelType> class lda_eval_exc_polar;
template <typename KernelType> class lda_eval_exc_vxc_unpolar;
template <typename KernelType> class lda_eval_exc_vxc_polar;
template <typename KernelType> class lda_eval_exc_inc_unpolar;
template <typename KernelType> class lda_eval_exc_inc_polar;
template <typename KernelType> class lda_eval_exc_vxc_inc_unpolar;
template <typename KernelType> class lda_eval_exc_vxc_inc_polar;

template <typename KernelType> class gga_eval_exc_unpolar;
template <typename KernelType> class gga_eval_exc_polar;
template <typename KernelType> class gga_eval_exc_vxc_unpolar;
template <typename KernelType> class gga_eval_exc_vxc_polar;
template <typename KernelType> class gga_eval_exc_inc_unpolar;
template <typename KernelType> class gga_eval_exc_inc_polar;
template <typename KernelType> class gga_eval_exc_vxc_inc_unpolar;
template <typename KernelType> class gga_eval_exc_vxc_inc_polar;







template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  queue->parallel_for<lda_eval_exc_unpolar<KernelType>>( sycl::range<1>(N),
    [=](sycl::id<1> idx) {
      device_eval_exc_helper_unpolar_kernel<KernelType>(
          N, rho, eps, idx
          );
    });

}

template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  queue->parallel_for<lda_eval_exc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_helper_polar_kernel<KernelType>(
        N, rho, eps, idx
      );
    });

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  queue->parallel_for<lda_eval_exc_vxc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_helper_unpolar_kernel<KernelType>(
        N, rho, eps, vxc, idx
      );
    });

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  queue->parallel_for<lda_eval_exc_vxc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_helper_polar_kernel<KernelType>(
        N, rho, eps, vxc, idx
      );
    });

}









template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  queue->parallel_for<lda_eval_exc_inc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_inc_helper_unpolar_kernel<KernelType>(
        scal_fact, N, rho, eps, idx
      );
    });

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  queue->parallel_for<lda_eval_exc_inc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_inc_helper_polar_kernel<KernelType>(
        scal_fact, N, rho, eps, idx
      );
    });

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  queue->parallel_for<lda_eval_exc_vxc_inc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>(
        scal_fact, N, rho, eps, vxc, idx
      );
    });

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  queue->parallel_for<lda_eval_exc_vxc_inc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>(
        scal_fact, N, rho, eps, vxc, idx
      );
    });

}










template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  queue->parallel_for<gga_eval_exc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_helper_unpolar_kernel<KernelType>(
        N, rho, sigma, eps, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  queue->parallel_for<gga_eval_exc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_helper_polar_kernel<KernelType>(
        N, rho, sigma, eps, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  queue->parallel_for<gga_eval_exc_vxc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_helper_unpolar_kernel<KernelType>(
        N, rho, sigma, eps, vrho, vsigma, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  queue->parallel_for<gga_eval_exc_vxc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_helper_polar_kernel<KernelType>(
        N, rho, sigma, eps, vrho, vsigma, idx
      );
    });

}









template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  queue->parallel_for<gga_eval_exc_inc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_inc_helper_unpolar_kernel<KernelType>(
        scal_fact, N, rho, sigma, eps, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  queue->parallel_for<gga_eval_exc_inc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_inc_helper_polar_kernel<KernelType>(
        scal_fact, N, rho, sigma, eps, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  queue->parallel_for<gga_eval_exc_vxc_inc_unpolar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>(
        scal_fact, N, rho, sigma, eps, vrho, vsigma, idx
      );
    });

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  queue->parallel_for<gga_eval_exc_vxc_inc_polar<KernelType>>( sycl::range<1>( N ),
    [=](sycl::id<1> idx) {
      device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>(
        scal_fact, N, rho, sigma, eps, vrho, vsigma, idx
      );
    });

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
