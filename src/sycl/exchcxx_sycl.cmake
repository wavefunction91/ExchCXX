set( EXCHCXX_SYCL_SOURCES
  sycl/xc_functional_device.cxx
  sycl/libxc_device.cxx
  sycl/builtin_sycl.cxx
)


target_sources( exchcxx PRIVATE ${EXCHCXX_SYCL_SOURCES} )

list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )
find_package( SYCL REQUIRED )
target_link_libraries( exchcxx PUBLIC SYCL::SYCL )

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fno-sycl-id-queries-fit-in-int"     EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT )
check_cxx_compiler_flag("-fsycl-device-code-split=per_kernel" EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
check_cxx_compiler_flag("-fno-sycl-early-optimizations"       EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )


if( EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-id-queries-fit-in-int>
  )
endif()

if( EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fsycl-device-code-split=per_kernel>
  )
endif()

if( EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-early-optimizations>
  )
endif()
