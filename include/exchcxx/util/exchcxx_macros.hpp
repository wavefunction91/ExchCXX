#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif



#include <exchcxx/util/param_macros.hpp>

// Device params
#ifdef EXCHCXX_ENABLE_CUDA

  #define DEVICE_PARAMS cudaStream_t stream
  #define DEVICE_PARAMS_NOTYPE stream

#endif

#ifdef EXCHCXX_ENABLE_SYCL

  #define DEVICE_PARAMS cl::sycl::queue& queue
  #define DEVICE_PARAMS_NOTYPE queue

  #define SYCL_KERNEL_PARAMS cl::sycl::nd_item<1> item_ct

#endif

#define RET_GENERATOR( APPROX, TYPE, func ) \
  func( APPROX ## _IPARAMS , APPROX ## _OPARAMS_ ## TYPE )

#define RET_INC_GENERATOR( APPROX, TYPE, func ) \
  func( double scal_fact, APPROX ## _IPARAMS , APPROX ## _OPARAMS_ ## TYPE )

#ifdef EXCHCXX_ENABLE_DEVICE

  #define RET_GENERATOR_DEVICE( APPROX, TYPE, func ) \
    func( DEV_ ## APPROX ## _IPARAMS , DEV_ ## APPROX ## _OPARAMS_ ## TYPE, \
          DEVICE_PARAMS )

  #define RET_INC_GENERATOR_DEVICE( APPROX, TYPE, func ) \
    func( double scal_fact, DEV_ ## APPROX ## _IPARAMS , DEV_ ## APPROX ## _OPARAMS_ ## TYPE, \
          DEVICE_PARAMS )

#endif

#ifdef EXCHCXX_ENABLE_SYCL

  #define RET_GENERATOR_SYCL_KERNEL( APPROX, TYPE, func ) \
    func( DEV_ ## APPROX ## _IPARAMS , DEV_ ## APPROX ## _OPARAMS_ ## TYPE, \
          SYCL_KERNEL_PARAMS )

  #define RET_INC_GENERATOR_SYCL_KERNEL( APPROX, TYPE, func ) \
    func( double scal_fact, DEV_ ## APPROX ## _IPARAMS , DEV_ ## APPROX ## _OPARAMS_ ## TYPE, \
          SYCL_KERNEL_PARAMS )

#endif


// LDA Generators

#define RET_LDA_EXC_GENERATOR(func)     RET_GENERATOR( LDA, EXC,     func )
#define RET_LDA_VXC_GENERATOR(func)     RET_GENERATOR( LDA, VXC,     func )
#define RET_LDA_FXC_GENERATOR(func)     RET_GENERATOR( LDA, FXC,     func )
#define RET_LDA_KXC_GENERATOR(func)     RET_GENERATOR( LDA, KXC,     func )
#define RET_LDA_EXC_VXC_GENERATOR(func) RET_GENERATOR( LDA, EXC_VXC, func )

#define LDA_EXC_GENERATOR(func)            void RET_LDA_EXC_GENERATOR(func)
#define LDA_VXC_GENERATOR(func)            void RET_LDA_VXC_GENERATOR(func)
#define LDA_FXC_GENERATOR(func)            void RET_LDA_FXC_GENERATOR(func)
#define LDA_KXC_GENERATOR(func)            void RET_LDA_KXC_GENERATOR(func)
#define LDA_EXC_VXC_GENERATOR(func)        void RET_LDA_EXC_VXC_GENERATOR(func)

#define RET_LDA_EXC_INC_GENERATOR(func)     RET_INC_GENERATOR( LDA, EXC,     func )
#define RET_LDA_VXC_INC_GENERATOR(func)     RET_INC_GENERATOR( LDA, VXC,     func )
#define RET_LDA_FXC_INC_GENERATOR(func)     RET_INC_GENERATOR( LDA, FXC,     func )
#define RET_LDA_KXC_INC_GENERATOR(func)     RET_INC_GENERATOR( LDA, KXC,     func )
#define RET_LDA_EXC_VXC_INC_GENERATOR(func) RET_INC_GENERATOR( LDA, EXC_VXC, func )

#define LDA_EXC_INC_GENERATOR(func)            void RET_LDA_EXC_INC_GENERATOR(func)
#define LDA_VXC_INC_GENERATOR(func)            void RET_LDA_VXC_INC_GENERATOR(func)
#define LDA_FXC_INC_GENERATOR(func)            void RET_LDA_FXC_INC_GENERATOR(func)
#define LDA_KXC_INC_GENERATOR(func)            void RET_LDA_KXC_INC_GENERATOR(func)
#define LDA_EXC_VXC_INC_GENERATOR(func)        void RET_LDA_EXC_VXC_INC_GENERATOR(func)


#ifdef EXCHCXX_ENABLE_DEVICE

  #define RET_LDA_EXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( LDA, EXC, func )
  #define RET_LDA_VXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( LDA, VXC, func )
  #define RET_LDA_FXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( LDA, FXC, func )
  #define RET_LDA_KXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( LDA, KXC, func )
  #define RET_LDA_EXC_VXC_GENERATOR_DEVICE(func) \
    RET_GENERATOR_DEVICE( LDA, EXC_VXC, func )
  

  #define LDA_EXC_GENERATOR_DEVICE(func)     \
    void RET_LDA_EXC_GENERATOR_DEVICE(func)
  #define LDA_VXC_GENERATOR_DEVICE(func)     \
    void RET_LDA_VXC_GENERATOR_DEVICE(func)
  #define LDA_FXC_GENERATOR_DEVICE(func)     \
    void RET_LDA_FXC_GENERATOR_DEVICE(func)
  #define LDA_KXC_GENERATOR_DEVICE(func)     \
    void RET_LDA_KXC_GENERATOR_DEVICE(func)
  #define LDA_EXC_VXC_GENERATOR_DEVICE(func) \
    void RET_LDA_EXC_VXC_GENERATOR_DEVICE(func)

  #define RET_LDA_EXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( LDA, EXC, func )
  #define RET_LDA_VXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( LDA, VXC, func )
  #define RET_LDA_FXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( LDA, FXC, func )
  #define RET_LDA_KXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( LDA, KXC, func )
  #define RET_LDA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    RET_INC_GENERATOR_DEVICE( LDA, EXC_VXC, func )
  

  #define LDA_EXC_INC_GENERATOR_DEVICE(func)     \
    void RET_LDA_EXC_INC_GENERATOR_DEVICE(func)
  #define LDA_VXC_INC_GENERATOR_DEVICE(func)     \
    void RET_LDA_VXC_INC_GENERATOR_DEVICE(func)
  #define LDA_FXC_INC_GENERATOR_DEVICE(func)     \
    void RET_LDA_FXC_INC_GENERATOR_DEVICE(func)
  #define LDA_KXC_INC_GENERATOR_DEVICE(func)     \
    void RET_LDA_KXC_INC_GENERATOR_DEVICE(func)
  #define LDA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    void RET_LDA_EXC_VXC_INC_GENERATOR_DEVICE(func)

#endif

#ifdef EXCHCXX_ENABLE_SYCL

  #define RET_LDA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( LDA, EXC, func )
  #define RET_LDA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( LDA, VXC, func )
  #define RET_LDA_FXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( LDA, FXC, func )
  #define RET_LDA_KXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( LDA, KXC, func )
  #define RET_LDA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    RET_GENERATOR_SYCL_KERNEL( LDA, EXC_VXC, func )
  

  #define LDA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_EXC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_VXC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_FXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_FXC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_KXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_KXC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    void RET_LDA_EXC_VXC_GENERATOR_SYCL_KERNEL(func)

  #define RET_LDA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( LDA, EXC, func )
  #define RET_LDA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( LDA, VXC, func )
  #define RET_LDA_FXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( LDA, FXC, func )
  #define RET_LDA_KXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( LDA, KXC, func )
  #define RET_LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    RET_INC_GENERATOR_SYCL_KERNEL( LDA, EXC_VXC, func )
  

  #define LDA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_EXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_VXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_FXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_FXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_KXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_LDA_KXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    void RET_LDA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func)

#endif


// GGA Generators 

#define RET_GGA_EXC_GENERATOR(func)     RET_GENERATOR( GGA, EXC,     func )
#define RET_GGA_VXC_GENERATOR(func)     RET_GENERATOR( GGA, VXC,     func )
#define RET_GGA_FXC_GENERATOR(func)     RET_GENERATOR( GGA, FXC,     func )
#define RET_GGA_KXC_GENERATOR(func)     RET_GENERATOR( GGA, KXC,     func )
#define RET_GGA_EXC_VXC_GENERATOR(func) RET_GENERATOR( GGA, EXC_VXC, func )

#define GGA_EXC_GENERATOR(func)            void RET_GGA_EXC_GENERATOR(func)
#define GGA_VXC_GENERATOR(func)            void RET_GGA_VXC_GENERATOR(func)
#define GGA_FXC_GENERATOR(func)            void RET_GGA_FXC_GENERATOR(func)
#define GGA_KXC_GENERATOR(func)            void RET_GGA_KXC_GENERATOR(func)
#define GGA_EXC_VXC_GENERATOR(func)        void RET_GGA_EXC_VXC_GENERATOR(func)

#define RET_GGA_EXC_INC_GENERATOR(func)     RET_INC_GENERATOR( GGA, EXC,     func )
#define RET_GGA_VXC_INC_GENERATOR(func)     RET_INC_GENERATOR( GGA, VXC,     func )
#define RET_GGA_FXC_INC_GENERATOR(func)     RET_INC_GENERATOR( GGA, FXC,     func )
#define RET_GGA_KXC_INC_GENERATOR(func)     RET_INC_GENERATOR( GGA, KXC,     func )
#define RET_GGA_EXC_VXC_INC_GENERATOR(func) RET_INC_GENERATOR( GGA, EXC_VXC, func )

#define GGA_EXC_INC_GENERATOR(func)            void RET_GGA_EXC_INC_GENERATOR(func)
#define GGA_VXC_INC_GENERATOR(func)            void RET_GGA_VXC_INC_GENERATOR(func)
#define GGA_FXC_INC_GENERATOR(func)            void RET_GGA_FXC_INC_GENERATOR(func)
#define GGA_KXC_INC_GENERATOR(func)            void RET_GGA_KXC_INC_GENERATOR(func)
#define GGA_EXC_VXC_INC_GENERATOR(func)        void RET_GGA_EXC_VXC_INC_GENERATOR(func)


#ifdef EXCHCXX_ENABLE_DEVICE

  #define RET_GGA_EXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( GGA, EXC, func )
  #define RET_GGA_VXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( GGA, VXC, func )
  #define RET_GGA_FXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( GGA, FXC, func )
  #define RET_GGA_KXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( GGA, KXC, func )
  #define RET_GGA_EXC_VXC_GENERATOR_DEVICE(func) \
    RET_GENERATOR_DEVICE( GGA, EXC_VXC, func )


  #define GGA_EXC_GENERATOR_DEVICE(func)     \
    void RET_GGA_EXC_GENERATOR_DEVICE(func)
  #define GGA_VXC_GENERATOR_DEVICE(func)     \
    void RET_GGA_VXC_GENERATOR_DEVICE(func)
  #define GGA_FXC_GENERATOR_DEVICE(func)     \
    void RET_GGA_FXC_GENERATOR_DEVICE(func)
  #define GGA_KXC_GENERATOR_DEVICE(func)     \
    void RET_GGA_KXC_GENERATOR_DEVICE(func)
  #define GGA_EXC_VXC_GENERATOR_DEVICE(func) \
    void RET_GGA_EXC_VXC_GENERATOR_DEVICE(func)

  #define RET_GGA_EXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( GGA, EXC, func )
  #define RET_GGA_VXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( GGA, VXC, func )
  #define RET_GGA_FXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( GGA, FXC, func )
  #define RET_GGA_KXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( GGA, KXC, func )
  #define RET_GGA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    RET_INC_GENERATOR_DEVICE( GGA, EXC_VXC, func )


  #define GGA_EXC_INC_GENERATOR_DEVICE(func)     \
    void RET_GGA_EXC_INC_GENERATOR_DEVICE(func)
  #define GGA_VXC_INC_GENERATOR_DEVICE(func)     \
    void RET_GGA_VXC_INC_GENERATOR_DEVICE(func)
  #define GGA_FXC_INC_GENERATOR_DEVICE(func)     \
    void RET_GGA_FXC_INC_GENERATOR_DEVICE(func)
  #define GGA_KXC_INC_GENERATOR_DEVICE(func)     \
    void RET_GGA_KXC_INC_GENERATOR_DEVICE(func)
  #define GGA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    void RET_GGA_EXC_VXC_INC_GENERATOR_DEVICE(func)

#endif


#ifdef EXCHCXX_ENABLE_SYCL

  #define RET_GGA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( GGA, EXC, func )
  #define RET_GGA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( GGA, VXC, func )
  #define RET_GGA_FXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( GGA, FXC, func )
  #define RET_GGA_KXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( GGA, KXC, func )
  #define RET_GGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    RET_GENERATOR_SYCL_KERNEL( GGA, EXC_VXC, func )


  #define GGA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_EXC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_VXC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_FXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_FXC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_KXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_KXC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    void RET_GGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func)

  #define RET_GGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( GGA, EXC, func )
  #define RET_GGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( GGA, VXC, func )
  #define RET_GGA_FXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( GGA, FXC, func )
  #define RET_GGA_KXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( GGA, KXC, func )
  #define RET_GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    RET_INC_GENERATOR_SYCL_KERNEL( GGA, EXC_VXC, func )


  #define GGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_FXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_FXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_KXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_GGA_KXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    void RET_GGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func)

#endif


// MGGA Generators

#define RET_MGGA_EXC_GENERATOR(func)     RET_GENERATOR( MGGA, EXC,     func )
#define RET_MGGA_VXC_GENERATOR(func)     RET_GENERATOR( MGGA, VXC,     func )
#define RET_MGGA_EXC_VXC_GENERATOR(func) RET_GENERATOR( MGGA, EXC_VXC, func )

#define MGGA_EXC_GENERATOR(func)            void RET_MGGA_EXC_GENERATOR(func)
#define MGGA_VXC_GENERATOR(func)            void RET_MGGA_VXC_GENERATOR(func)
#define MGGA_EXC_VXC_GENERATOR(func)        void RET_MGGA_EXC_VXC_GENERATOR(func)

#define RET_MGGA_EXC_INC_GENERATOR(func)     RET_INC_GENERATOR( MGGA, EXC,     func )
#define RET_MGGA_VXC_INC_GENERATOR(func)     RET_INC_GENERATOR( MGGA, VXC,     func )
#define RET_MGGA_EXC_VXC_INC_GENERATOR(func) RET_INC_GENERATOR( MGGA, EXC_VXC, func )

#define MGGA_EXC_INC_GENERATOR(func)            void RET_MGGA_EXC_INC_GENERATOR(func)
#define MGGA_VXC_INC_GENERATOR(func)            void RET_MGGA_VXC_INC_GENERATOR(func)
#define MGGA_EXC_VXC_INC_GENERATOR(func)        void RET_MGGA_EXC_VXC_INC_GENERATOR(func)


#ifdef EXCHCXX_ENABLE_DEVICE

  #define RET_MGGA_EXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( MGGA, EXC, func )
  #define RET_MGGA_VXC_GENERATOR_DEVICE(func)     \
    RET_GENERATOR_DEVICE( MGGA, VXC, func )
  #define RET_MGGA_EXC_VXC_GENERATOR_DEVICE(func) \
    RET_GENERATOR_DEVICE( MGGA, EXC_VXC, func )


  #define MGGA_EXC_GENERATOR_DEVICE(func)     \
    void RET_MGGA_EXC_GENERATOR_DEVICE(func)
  #define MGGA_VXC_GENERATOR_DEVICE(func)     \
    void RET_MGGA_VXC_GENERATOR_DEVICE(func)
  #define MGGA_EXC_VXC_GENERATOR_DEVICE(func) \
    void RET_MGGA_EXC_VXC_GENERATOR_DEVICE(func)


  #define RET_MGGA_EXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( MGGA, EXC, func )
  #define RET_MGGA_VXC_INC_GENERATOR_DEVICE(func)     \
    RET_INC_GENERATOR_DEVICE( MGGA, VXC, func )
  #define RET_MGGA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    RET_INC_GENERATOR_DEVICE( MGGA, EXC_VXC, func )


  #define MGGA_EXC_INC_GENERATOR_DEVICE(func)     \
    void RET_MGGA_EXC_INC_GENERATOR_DEVICE(func)
  #define MGGA_VXC_INC_GENERATOR_DEVICE(func)     \
    void RET_MGGA_VXC_INC_GENERATOR_DEVICE(func)
  #define MGGA_EXC_VXC_INC_GENERATOR_DEVICE(func) \
    void RET_MGGA_EXC_VXC_INC_GENERATOR_DEVICE(func)

#endif



#ifdef EXCHCXX_ENABLE_SYCL

  #define RET_MGGA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( MGGA, EXC, func )
  #define RET_MGGA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    RET_GENERATOR_SYCL_KERNEL( MGGA, VXC, func )
  #define RET_MGGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    RET_GENERATOR_SYCL_KERNEL( MGGA, EXC_VXC, func )


  #define MGGA_EXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_MGGA_EXC_GENERATOR_SYCL_KERNEL(func)
  #define MGGA_VXC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_MGGA_VXC_GENERATOR_SYCL_KERNEL(func)
  #define MGGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func) \
    void RET_MGGA_EXC_VXC_GENERATOR_SYCL_KERNEL(func)


  #define RET_MGGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( MGGA, EXC, func )
  #define RET_MGGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    RET_INC_GENERATOR_SYCL_KERNEL( MGGA, VXC, func )
  #define RET_MGGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    RET_INC_GENERATOR_SYCL_KERNEL( MGGA, EXC_VXC, func )


  #define MGGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_MGGA_EXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define MGGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)     \
    void RET_MGGA_VXC_INC_GENERATOR_SYCL_KERNEL(func)
  #define MGGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func) \
    void RET_MGGA_EXC_VXC_INC_GENERATOR_SYCL_KERNEL(func)

#endif






#define FORWARD_XC_ARGS(APPROX, TYPE, func, finternal, qualifiers) \
  APPROX ## _ ## TYPE ## _GENERATOR(func) qualifiers {             \
    finternal( APPROX ## _IPARAMS_NOTYPE ,                         \
               APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE              \
    );                                                             \
  }

#define FORWARD_XC_INC_ARGS(APPROX, TYPE, func, finternal, qualifiers) \
  APPROX ## _ ## TYPE ## _INC_GENERATOR(func) qualifiers {             \
    finternal( scal_fact, APPROX ## _IPARAMS_NOTYPE ,                  \
               APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE                  \
    );                                                                 \
  }


#ifdef EXCHCXX_ENABLE_DEVICE

  #define FORWARD_XC_ARGS_DEVICE(APPROX, TYPE, func, finternal, qualifiers) \
    APPROX ## _ ## TYPE ## _GENERATOR_DEVICE(func) qualifiers {             \
      finternal( APPROX ## _IPARAMS_NOTYPE ,                                \
                 APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE,                    \
                 DEVICE_PARAMS_NOTYPE );                                    \
    }

  #define FORWARD_XC_INC_ARGS_DEVICE(APPROX, TYPE, func, finternal, qualifiers) \
    APPROX ## _ ## TYPE ## _INC_GENERATOR_DEVICE(func) qualifiers {             \
      finternal( scal_fact, APPROX ## _IPARAMS_NOTYPE ,                         \
                 APPROX ## _OPARAMS_ ## TYPE ## _NOTYPE,                        \
                 DEVICE_PARAMS_NOTYPE );                                        \
    }

#endif



