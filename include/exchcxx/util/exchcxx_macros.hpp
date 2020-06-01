#pragma once

// Device params
#define DEVICE_PARAMS device::cuda_stream_t* stream
#define DEVICE_PARAMS_NOTYPE stream


// LDA Generators
#define LDA_IPARAMS          const int N, const double* rho
#define LDA_OPARAMS_EXC      double* eps
#define LDA_OPARAMS_VXC      double* vxc
#define LDA_OPARAMS_FXC      double* fxc
#define LDA_OPARAMS_KXC      double* kxc
#define LDA_OPARAMS_EXC_VXC  double* eps, double* vxc

#define LDA_IPARAMS_NOTYPE          N, rho
#define LDA_OPARAMS_EXC_NOTYPE      eps
#define LDA_OPARAMS_VXC_NOTYPE      vxc
#define LDA_OPARAMS_FXC_NOTYPE      fxc
#define LDA_OPARAMS_KXC_NOTYPE      kxc
#define LDA_OPARAMS_EXC_VXC_NOTYPE  eps, vxc


#define LDA_EXC_GENERATOR(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_EXC )

#define LDA_VXC_GENERATOR(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_VXC )

#define LDA_FXC_GENERATOR(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_FXC )

#define LDA_KXC_GENERATOR(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_KXC )

#define LDA_EXC_VXC_GENERATOR(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_EXC_VXC )

#define LDA_EXC_GENERATOR_DEVICE(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_EXC, DEVICE_PARAMS )

#define LDA_VXC_GENERATOR_DEVICE(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_VXC, DEVICE_PARAMS )

#define LDA_FXC_GENERATOR_DEVICE(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_FXC, DEVICE_PARAMS )

#define LDA_KXC_GENERATOR_DEVICE(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_KXC, DEVICE_PARAMS )

#define LDA_EXC_VXC_GENERATOR_DEVICE(func) \
  void func( LDA_IPARAMS , LDA_OPARAMS_EXC_VXC, DEVICE_PARAMS )



// GGA Generators 
#define GGA_IPARAMS          const int N, const double* rho, const double* sigma
#define GGA_OPARAMS_EXC      double* eps
#define GGA_OPARAMS_VXC      double* vrho, double* vsigma
#define GGA_OPARAMS_FXC      double *v2rho2, double *v2rhosigma, double *v2sigma2
#define GGA_OPARAMS_KXC      double *v3rho3, double *v3rho2sigma, \
                             double *v3rhosigma2, double *v3sigma3 

#define GGA_OPARAMS_EXC_VXC  GGA_OPARAMS_EXC , GGA_OPARAMS_VXC

#define GGA_IPARAMS_NOTYPE          N, rho, sigma
#define GGA_OPARAMS_EXC_NOTYPE      eps
#define GGA_OPARAMS_VXC_NOTYPE      vrho, vsigma
#define GGA_OPARAMS_FXC_NOTYPE      v2rho2, v2rhosigma, v2sigma2
#define GGA_OPARAMS_KXC_NOTYPE      v3rho3, v3rho2sigma, \
                                    v3rhosigma2, v3sigma3 

#define GGA_OPARAMS_EXC_VXC_NOTYPE  GGA_OPARAMS_EXC_NOTYPE , GGA_OPARAMS_VXC_NOTYPE

#define GGA_EXC_GENERATOR(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_EXC )

#define GGA_VXC_GENERATOR(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_VXC )

#define GGA_FXC_GENERATOR(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_FXC )

#define GGA_KXC_GENERATOR(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_KXC )

#define GGA_EXC_VXC_GENERATOR(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_EXC_VXC )

#define GGA_EXC_GENERATOR_DEVICE(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_EXC, DEVICE_PARAMS )

#define GGA_VXC_GENERATOR_DEVICE(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_VXC, DEVICE_PARAMS )

#define GGA_FXC_GENERATOR_DEVICE(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_FXC, DEVICE_PARAMS )

#define GGA_KXC_GENERATOR_DEVICE(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_KXC, DEVICE_PARAMS )

#define GGA_EXC_VXC_GENERATOR_DEVICE(func) \
  void func( GGA_IPARAMS , GGA_OPARAMS_EXC_VXC, DEVICE_PARAMS )





// MGGA Generators
#define MGGA_IPARAMS          const int N, const double* rho, const double* sigma, \
                              const double* lapl, const double* tau

#define MGGA_OPARAMS_EXC      double* eps

#define MGGA_OPARAMS_VXC      double* vrho, double* vsigma, double* vlapl, \
                              double* vtau

#define MGGA_OPARAMS_EXC_VXC  MGGA_OPARAMS_EXC , MGGA_OPARAMS_VXC



#define MGGA_IPARAMS_NOTYPE           N,  rho,  sigma, lapl,  tau
#define MGGA_OPARAMS_EXC_NOTYPE       eps
#define MGGA_OPARAMS_VXC_NOTYPE       vrho,  vsigma,  vlapl, vtau

#define MGGA_OPARAMS_EXC_VXC_NOTYPE  MGGA_OPARAMS_EXC_NOTYPE , MGGA_OPARAMS_VXC_NOTYPE

#define MGGA_EXC_GENERATOR(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_EXC )

#define MGGA_VXC_GENERATOR(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_VXC )

#define MGGA_EXC_VXC_GENERATOR(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_EXC_VXC )

#define MGGA_EXC_GENERATOR_DEVICE(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_EXC, DEVICE_PARAMS )

#define MGGA_VXC_GENERATOR_DEVICE(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_VXC, DEVICE_PARAMS )

#define MGGA_EXC_VXC_GENERATOR_DEVICE(func) \
  void func( MGGA_IPARAMS , MGGA_OPARAMS_EXC_VXC, DEVICE_PARAMS )







#define FORWARD_XC_ARGS(APPROX, TYPE, func, finternal, qualifiers) \
  APPROX ## _ ## TYPE ## _GENERATOR(func) qualifiers {             \
    finternal( APPROX ## _IPARAMS_NOTYPE , APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE ); \
  }

#define FORWARD_XC_ARGS_DEVICE(APPROX, TYPE, func, finternal, qualifiers) \
  APPROX ## _ ## TYPE ## _GENERATOR_DEVICE(func) qualifiers {             \
    finternal( APPROX ## _IPARAMS_NOTYPE , APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE, DEVICE_PARAMS_NOTYPE ); \
  }




