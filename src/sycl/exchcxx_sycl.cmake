set( EXCHCXX_SYCL_SOURCES 
  sycl/xc_functional_device.cxx 
  sycl/libxc_device.cxx 
  sycl/builtin_sycl.cxx 
)


target_sources( exchcxx PRIVATE ${EXCHCXX_SYCL_SOURCES} )

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fno-sycl-early-optimizations" EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
if( EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-early-optimizations>
  )
endif()

#target_link_libraries( exchcxx PUBLIC SYCL::SYCL )
