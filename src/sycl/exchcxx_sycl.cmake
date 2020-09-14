set( EXCHCXX_SYCL_SOURCES 
  sycl/xc_functional_device.cxx 
  sycl/libxc_device.cxx 
  sycl/builtin_sycl.cxx 
)


target_sources( exchcxx PRIVATE ${EXCHCXX_SYCL_SOURCES} )
#target_link_libraries( exchcxx PUBLIC SYCL::SYCL )
