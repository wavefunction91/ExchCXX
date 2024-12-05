set( EXCHCXX_HIP_SOURCES 
  hip/xc_functional_device.hip 
  hip/builtin.hip
)
if( EXCHCXX_ENABLE_LIBXC )
  list(APPEND EXCHCXX_HIP_SOURCES hip/libxc_device.hip)
endif()

find_package( hip REQUIRED )

target_sources( exchcxx PRIVATE ${EXCHCXX_HIP_SOURCES} )
target_link_libraries( exchcxx PUBLIC hip::host )
