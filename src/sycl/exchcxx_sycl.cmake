set( EXCHCXX_SYCL_SOURCES
  sycl/xc_functional_device.cxx
  sycl/libxc_device.cxx
  sycl/builtin_sycl.cxx
)


target_sources( exchcxx PRIVATE ${EXCHCXX_SYCL_SOURCES} )

#set(SYCL_FLAGS -fsycl -fsycl-unnamed-lambda -fsycl-device-code-split=per_kernel -fsycl-targets=spir64_gen -Xsycl-target-backend "-device 12.1.0,12.4.0")
#set(SYCL_FLAGS -fsycl -fsycl-unnamed-lambda -fsycl-device-code-split=per_kernel -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend --cuda-gpu-arch=sm_80)
#set(SYCL_FLAGS -fsycl -fsycl-unnamed-lambda -fsycl-device-code-split=per_kernel -fsycl-targets=amdgcn-amd-amdhsa -Xsycl-target-backend --offload-arch=gfx90a)

include(CheckCXXCompilerFlag)
target_compile_options( exchcxx
  PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: ${SYCL_FLAGS} > )

check_cxx_compiler_flag("-fno-sycl-early-optimizations" EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
if( EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-early-optimizations>
  )
endif()

target_link_libraries( exchcxx PUBLIC ${SYCL_FLAGS} )
