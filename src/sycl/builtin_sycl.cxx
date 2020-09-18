#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>

namespace ExchCXX {
namespace detail {

template <typename Integral1, typename Integral2>
int64_t div_ceil(Integral1 x, Integral2 y) {
  int64_t x_ll = x;
  int64_t y_ll = y;

  auto d = std::div(x_ll, y_ll);
  return d.quot + !!d.rem;
}


template <typename KernelType>
inline LDA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    const double rho_use = cl::sycl::max( rho[idx], 0. );
    traits::eval_exc_unpolar( rho_use, eps[idx] );

  }

}

template <typename KernelType>
inline LDA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto rho_i = rho + 2*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );

    traits::eval_exc_polar( rho_a_use, rho_b_use, eps[idx] );

  }

}

template <typename KernelType>
inline LDA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    const double rho_use = cl::sycl::max( rho[idx], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, eps[idx], vxc[idx] );

  }

}

template <typename KernelType>
inline LDA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto rho_i = rho + 2*idx;
    auto vxc_i = vxc + 2*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );

    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, eps[idx],
      vxc_i[0], vxc_i[1] );

  }

}

template <typename KernelType>
inline LDA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e;
  if( idx < N ) {

    const double rho_use = cl::sycl::max( rho[idx], 0. );
    traits::eval_exc_unpolar( rho_use, e );
    eps[idx] += scal_fact * e;

  }

}

template <typename KernelType>
inline LDA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto rho_i = rho + 2*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );

    double e;
    traits::eval_exc_polar( rho_a_use, rho_b_use, e );

    eps[idx] += scal_fact * e;

  }

}

template <typename KernelType>
inline LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e,v;
  if( idx < N ) {

    const double rho_use = cl::sycl::max( rho[idx], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, e, v );
    eps[idx] += scal_fact * e;
    vxc[idx] += scal_fact * v;

  }

}

template <typename KernelType>
inline LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto rho_i = rho + 2*idx;
    auto vxc_i = vxc + 2*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );

    double v_a, v_b, e;
    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, e, v_a, v_b);
    eps[idx] += scal_fact * e;
    vxc_i[0] += scal_fact * v_a;
    vxc_i[1] += scal_fact * v_b;

  }

}

template <typename KernelType>
inline GGA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    const double rho_use   = cl::sycl::max( rho[idx],   0.    );
    const double sigma_use = cl::sycl::max( sigma[idx], 1e-40 );
    traits::eval_exc_unpolar( rho_use, sigma_use, eps[idx] );

  }

}

template <typename KernelType>
inline GGA_EXC_GENERATOR_SYCL_KERNEL( device_eval_exc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto* rho_i   = rho   + 2*idx;
    auto* sigma_i = sigma + 3*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );
    const double sigma_aa_use = cl::sycl::max( sigma_i[0], 1e-40 );
    const double sigma_bb_use = cl::sycl::max( sigma_i[2], 1e-40 );
    const double sigma_ab_use = cl::sycl::max(
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );

    traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use,
      sigma_ab_use, sigma_bb_use, eps[idx] );

  }

}

template <typename KernelType>
inline GGA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    const double rho_use   = cl::sycl::max( rho[idx],   0.    );
    const double sigma_use = cl::sycl::max( sigma[idx], 1e-40 );
    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, eps[idx],
      vrho[idx], vsigma[idx] );

  }

}

template <typename KernelType>
inline GGA_EXC_VXC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto* rho_i    = rho   + 2*idx;
    auto* sigma_i  = sigma + 3*idx;
    auto* vrho_i   = vrho   + 2*idx;
    auto* vsigma_i = vsigma + 3*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );
    const double sigma_aa_use = cl::sycl::max( sigma_i[0], 1e-40 );
    const double sigma_bb_use = cl::sycl::max( sigma_i[2], 1e-40 );
    const double sigma_ab_use = cl::sycl::max(
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );


    traits::eval_exc_vxc_polar( rho_a_use, rho_b_use, sigma_aa_use,
      sigma_ab_use, sigma_bb_use, eps[idx], vrho_i[0], vrho_i[1],
      vsigma_i[0], vsigma_i[1], vsigma_i[2] );

  }

}


template <typename KernelType>
inline GGA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e;
  if( idx < N ) {

    const double rho_use   = cl::sycl::max( rho[idx],   0.    );
    const double sigma_use = cl::sycl::max( sigma[idx], 1e-40 );

    traits::eval_exc_unpolar( rho_use, sigma_use, e );
    eps[idx] += scal_fact * e;


  }

}

template <typename KernelType>
inline GGA_EXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto* rho_i   = rho   + 2*idx;
    auto* sigma_i = sigma + 3*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );
    const double sigma_aa_use = cl::sycl::max( sigma_i[0], 1e-40 );
    const double sigma_bb_use = cl::sycl::max( sigma_i[2], 1e-40 );
    const double sigma_ab_use = cl::sycl::max(
      sigma_i[1], -(sigma_i[0] + sigma_i[1]) / 2.
    );

    double e;
    traits::eval_exc_polar( rho_a_use, rho_b_use, sigma_aa_use,
      sigma_ab_use, sigma_bb_use, e );
    eps[idx] += scal_fact * e;


  }

}

template <typename KernelType>
inline GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_unpolar_kernel ) {

  using traits = kernel_traits<KernelType>;

  double e, vr, vs;
  if( idx < N ) {

    const double rho_use   = cl::sycl::max( rho[idx],   0.    );
    const double sigma_use = cl::sycl::max( sigma[idx], 1e-40 );

    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, e, vr, vs );
    eps[idx]    += scal_fact * e;
    vrho[idx]   += scal_fact * vr;
    vsigma[idx] += scal_fact * vs;

  }

}

template <typename KernelType>
inline GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL( device_eval_exc_vxc_inc_helper_polar_kernel ) {

  using traits = kernel_traits<KernelType>;

  if( idx < N ) {

    auto* rho_i    = rho   + 2*idx;
    auto* sigma_i  = sigma + 3*idx;
    auto* vrho_i   = vrho   + 2*idx;
    auto* vsigma_i = vsigma + 3*idx;

    const double rho_a_use = cl::sycl::max( rho_i[0], 0. );
    const double rho_b_use = cl::sycl::max( rho_i[1], 0. );
    const double sigma_aa_use = cl::sycl::max( sigma_i[0], 1e-40 );
    const double sigma_bb_use = cl::sycl::max( sigma_i[2], 1e-40 );
    const double sigma_ab_use = cl::sycl::max(
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

}














template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {
  std::cout << N << ", " << rho << ", " << eps << std::endl;

  try {
      queue->submit([&](cl::sycl::handler& cgh) {
        cgh.parallel_for( cl::sycl::range<1>(N),
          [=](cl::sycl::id<1> idx) {
            device_eval_exc_helper_unpolar_kernel<KernelType>(
                N, rho, eps, idx
                );
          });
      });
  }
  catch( cl::sycl::exception const& ex ) {
      std::cout << "SYCL: " << ex.what() << std::endl;
    throw;
  }
}

template <typename KernelType>
LDA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_helper_polar_kernel<KernelType>(
          N, rho, eps, idx
        );
      });

  });

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_helper_unpolar_kernel<KernelType>(
          N, rho, eps, vxc, idx
        );
      });

  });

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_helper_polar_kernel<KernelType>(
          N, rho, eps, vxc, idx
        );
      });

  });

}









template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_inc_helper_unpolar_kernel<KernelType>(
          scal_fact, N, rho, eps, idx
        );
      });

  });

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_inc_helper_polar_kernel<KernelType>(
          scal_fact, N, rho, eps, idx
        );
      });

  });

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>(
          scal_fact, N, rho, eps, vxc, idx
        );
      });

  });

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>(
          scal_fact, N, rho, eps, vxc, idx
        );
      });

  });

}










template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_helper_unpolar_kernel<KernelType>(
          N, rho, sigma, eps, idx
        );
      });

  });


}

template <typename KernelType>
GGA_EXC_GENERATOR_DEVICE( device_eval_exc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_helper_polar_kernel<KernelType>(
          N, rho, sigma, eps, idx
        );
      });

  });

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_helper_unpolar_kernel<KernelType>(
          N, rho, sigma, eps, vrho, vsigma, idx
        );
      });

  });

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR_DEVICE( device_eval_exc_vxc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_helper_polar_kernel<KernelType>(
          N, rho, sigma, eps, vrho, vsigma, idx
        );
      });

  });

}









template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_inc_helper_unpolar_kernel<KernelType>(
          scal_fact, N, rho, sigma, eps, idx
        );
      });

  });

}

template <typename KernelType>
GGA_EXC_INC_GENERATOR_DEVICE( device_eval_exc_inc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_inc_helper_polar_kernel<KernelType>(
          scal_fact, N, rho, sigma, eps, idx
        );
      });

  });

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_unpolar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_inc_helper_unpolar_kernel<KernelType>(
          scal_fact, N, rho, sigma, eps, vrho, vsigma, idx
        );
      });

  });

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR_DEVICE( device_eval_exc_vxc_inc_helper_polar ) {

  queue->submit( [&](cl::sycl::handler &cgh) {

    cgh.parallel_for( cl::sycl::range<1>( N ),
      [=](cl::sycl::id<1> idx) {
        device_eval_exc_vxc_inc_helper_polar_kernel<KernelType>(
          scal_fact, N, rho, sigma, eps, vrho, vsigma, idx
        );
      });

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
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE_C );

GGA_GENERATE_DEVICE_HELPERS( BuiltinB3LYP );
GGA_GENERATE_DEVICE_HELPERS( BuiltinPBE0  );



}
}
