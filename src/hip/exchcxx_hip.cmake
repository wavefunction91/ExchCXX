set( EXCHCXX_HIP_SOURCES 
  hip/xc_functional_device.cxx 
  hip/libxc_device.cxx 
  hip/builtin.cxx
)

target_sources( exchcxx PRIVATE ${EXCHCXX_HIP_SOURCES} )
