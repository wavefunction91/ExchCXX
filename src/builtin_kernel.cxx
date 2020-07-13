#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>


namespace ExchCXX {


// LDA interface
LDA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_lda_interface();
};
LDA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_lda_interface();
};

// GGA interface
GGA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_gga_interface();
};
GGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_gga_interface();
};

// MGGA interface
MGGA_EXC_GENERATOR( BuiltinKernel::eval_exc ) const {
  disabled_mgga_interface();
};
MGGA_EXC_VXC_GENERATOR( BuiltinKernel::eval_exc_vxc ) const {
  disabled_mgga_interface();
};


// INC interfaces
LDA_EXC_INC_GENERATOR( BuiltinKernel::eval_exc_inc ) const {
  disabled_inc_interface();
};
LDA_EXC_VXC_INC_GENERATOR( BuiltinKernel::eval_exc_vxc_inc ) const {
  disabled_inc_interface();
};
GGA_EXC_INC_GENERATOR( BuiltinKernel::eval_exc_inc ) const {
  disabled_inc_interface();
};
GGA_EXC_VXC_INC_GENERATOR( BuiltinKernel::eval_exc_vxc_inc ) const {
  disabled_inc_interface();
};
MGGA_EXC_INC_GENERATOR( BuiltinKernel::eval_exc_inc ) const {
  disabled_inc_interface();
};
MGGA_EXC_VXC_INC_GENERATOR( BuiltinKernel::eval_exc_vxc_inc ) const {
  disabled_inc_interface();
};

#ifdef EXCHCXX_ENABLE_DEVICE

// LDA interface
LDA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device ) const {
  disabled_lda_device_interface();
};
LDA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_lda_device_interface();
};

// GGA interface
GGA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device ) const {
  disabled_gga_device_interface();
};
GGA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_gga_device_interface();
};

// MGGA interface
MGGA_EXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_device ) const {
  disabled_mgga_device_interface();
};
MGGA_EXC_VXC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_device ) const {
  disabled_mgga_device_interface();
};

// INC interface
LDA_EXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_inc_device ) const {
  disabled_inc_device_interface();
};
LDA_EXC_VXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_inc_device ) const {
  disabled_inc_device_interface();
};
GGA_EXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_inc_device ) const {
  disabled_inc_device_interface();
};
GGA_EXC_VXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_inc_device ) const {
  disabled_inc_device_interface();
};
MGGA_EXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_inc_device ) const {
  disabled_inc_device_interface();
};
MGGA_EXC_VXC_INC_GENERATOR_DEVICE( BuiltinKernel::eval_exc_vxc_inc_device ) const {
  disabled_inc_device_interface();
};

#endif



namespace detail {

template <typename KernelType>
LDA_EXC_GENERATOR( host_eval_exc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use = std::max( rho[i], 0. );
    traits::eval_exc_unpolar( rho_use, eps[i] );

  }

}

template <typename KernelType>
LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use = std::max( rho[i], 0. );
    traits::eval_exc_vxc_unpolar( rho_use, eps[i], vxc[i] );

  }

}

template <typename KernelType>
LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use = std::max( rho[i], 0. );
    double e;
    traits::eval_exc_unpolar( rho_use, e );
    eps[i] += scal_fact * e;

  }

}

template <typename KernelType>
LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use = std::max( rho[i], 0. );
    double v,e;
    traits::eval_exc_vxc_unpolar( rho_use, e, v );
    eps[i] += scal_fact * e;
    vxc[i] += scal_fact * v;

  }

}

template <typename KernelType>
GGA_EXC_GENERATOR( host_eval_exc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use   = std::max( rho[i],   0.    );
    const double sigma_use = std::max( sigma[i], 1e-40 );
    traits::eval_exc_unpolar( rho_use, sigma_use, eps[i] );

  }

}

template <typename KernelType>
GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use   = std::max( rho[i],   0.    );
    const double sigma_use = std::max( sigma[i], 1e-40 );
    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, 
      eps[i], vrho[i], vsigma[i] );

  }

}


template <typename KernelType>
GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use   = std::max( rho[i],   0.    );
    const double sigma_use = std::max( sigma[i], 1e-40 );

    double e;
    traits::eval_exc_unpolar( rho_use, sigma_use, e );
    eps[i] += scal_fact * e;

  }

}

template <typename KernelType>
GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper ) {

  using traits = kernel_traits<KernelType>;
  
  for( size_t i = 0; i < N; ++i ) {

    const double rho_use   = std::max( rho[i],   0.    );
    const double sigma_use = std::max( sigma[i], 1e-40 );

    double e, vr, vs;
    traits::eval_exc_vxc_unpolar( rho_use, sigma_use, e, vr, vs );
    eps[i]    += scal_fact * e;
    vrho[i]   += scal_fact * vr;
    vsigma[i] += scal_fact * vs;

  }

}

#define LDA_GENERATE_HOST_HELPERS(KERN) \
  template LDA_EXC_GENERATOR( host_eval_exc_helper<KERN> ); \
  template LDA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper<KERN> ); \
  template LDA_EXC_INC_GENERATOR( host_eval_exc_inc_helper<KERN> ); \
  template LDA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper<KERN> ); 

#define GGA_GENERATE_HOST_HELPERS(KERN) \
  template GGA_EXC_GENERATOR( host_eval_exc_helper<KERN> ); \
  template GGA_EXC_VXC_GENERATOR( host_eval_exc_vxc_helper<KERN> ); \
  template GGA_EXC_INC_GENERATOR( host_eval_exc_inc_helper<KERN> ); \
  template GGA_EXC_VXC_INC_GENERATOR( host_eval_exc_vxc_inc_helper<KERN> ); 

LDA_GENERATE_HOST_HELPERS( BuiltinSlaterExchange );
LDA_GENERATE_HOST_HELPERS( BuiltinVWN3 );
LDA_GENERATE_HOST_HELPERS( BuiltinVWN_RPA );

GGA_GENERATE_HOST_HELPERS( BuiltinB88   );
GGA_GENERATE_HOST_HELPERS( BuiltinLYP   );
GGA_GENERATE_HOST_HELPERS( BuiltinPBE_X );
GGA_GENERATE_HOST_HELPERS( BuiltinPBE_C );

GGA_GENERATE_HOST_HELPERS( BuiltinB3LYP  );
GGA_GENERATE_HOST_HELPERS( BuiltinPBE0  );

}



}

