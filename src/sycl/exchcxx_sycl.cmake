set( EXCHCXX_SYCL_SOURCES
  sycl/xc_functional_device.cxx
  sycl/builtin_sycl.cxx
)
if( EXCHCXX_ENABLE_LIBXC )
  list(APPEND EXCHCXX_SYCL_SOURCES sycl/libxc_device.cxx)
endif()


target_sources( exchcxx PRIVATE ${EXCHCXX_SYCL_SOURCES} )

list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )
find_package( SYCL REQUIRED )
target_link_libraries( exchcxx PUBLIC SYCL::SYCL )

# --- AoT-builds SYCL target alias pass-through ---
set(_EXCHCXX_SYCL_ALLOWED
    intel_gpu_pvc
    nvidia_gpu_sm_80
    nvidia_gpu_sm_90
    amd_gpu_gfx90a
    amd_gpu_gfx942)

if(DEFINED EXCHCXX_SYCL_TARGET AND NOT EXCHCXX_SYCL_TARGET STREQUAL "")
  list(FIND _EXCHCXX_SYCL_ALLOWED "${EXCHCXX_SYCL_TARGET}" _exchcxx_sycl_idx)
  if(_exchcxx_sycl_idx EQUAL -1)
    message(FATAL_ERROR
      "Invalid EXCHCXX_SYCL_TARGET='${EXCHCXX_SYCL_TARGET}'. "
      "Allowed values: ${_EXCHCXX_SYCL_ALLOWED}")
  endif()

  # Apply ONLY to this target (both compile & link)
  target_compile_options( exchcxx PRIVATE -fsycl-device-only -fsycl-targets=${EXCHCXX_SYCL_TARGET} )
  target_link_options( exchcxx PRIVATE -fsycl-device-only -fsycl-targets=${EXCHCXX_SYCL_TARGET} )
  message(STATUS "ExchCXX SYCL AoT enabled for target: ${EXCHCXX_SYCL_TARGET}")

  # target_compile_options( exchcxx PRIVATE -Wno-unused-parameter -Wno-unused-variable -fsycl-device-only -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc" )
  # target_link_options( exchcxx PRIVATE -fsycl-device-only -fsycl-targets=spir64_gen -Xsycl-target-backend "-device pvc" )

endif()

target_compile_options(exchcxx PRIVATE  $<$<COMPILE_LANGUAGE:CXX>:-ffp-model=precise>)
target_link_options(exchcxx PRIVATE -fsycl-max-parallel-link-jobs=20)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fno-sycl-id-queries-fit-in-int"     EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT )
check_cxx_compiler_flag("-fsycl-device-code-split=per_source" EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_SOURCE )
check_cxx_compiler_flag("-fno-sycl-early-optimizations"       EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )


if( EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-id-queries-fit-in-int>
  )
endif()

if( EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_SOURCE )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fsycl-device-code-split=per_source>
  )
endif()

if( EXCHCXX_SYCL_HAS_NO_EARLY_OPTIMIZATIONS )
  target_compile_options( exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -fno-sycl-early-optimizations>
  )
endif()
