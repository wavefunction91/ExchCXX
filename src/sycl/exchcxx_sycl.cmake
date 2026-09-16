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
# User-facing aliases
set(_EXCHCXX_SYCL_ALLOWED
  intel_gpu_pvc
  spir64_x86_64
  nvidia_gpu_sm_80
  nvidia_gpu_sm_90
  amd_gpu_gfx90a
  amd_gpu_gfx942
)

if(DEFINED EXCHCXX_SYCL_TARGET AND NOT EXCHCXX_SYCL_TARGET STREQUAL "")
  list(FIND _EXCHCXX_SYCL_ALLOWED "${EXCHCXX_SYCL_TARGET}" _exchcxx_sycl_idx)
  if(_exchcxx_sycl_idx EQUAL -1)
    message(FATAL_ERROR
      "Invalid EXCHCXX_SYCL_TARGET='${EXCHCXX_SYCL_TARGET}'. "
      "Allowed values: ${_EXCHCXX_SYCL_ALLOWED}")
  endif()

  unset(_exchcxx_sycl_compile_opts)
  unset(_exchcxx_sycl_link_opts)

  if(EXCHCXX_SYCL_TARGET STREQUAL "intel_gpu_pvc")
    list(APPEND _exchcxx_sycl_compile_opts
      -fsycl-default-sub-group-size=32
      -fsycl-targets=spir64_gen
      "SHELL:-Xsycl-target-backend \"-device pvc\""
    )
    list(APPEND _exchcxx_sycl_link_opts
      "SHELL:-ftarget-register-alloc-mode=pvc:large"
      "SHELL:-fsycl-targets=spir64_gen"
      "SHELL:-Xsycl-target-backend \"-device pvc\""
    )
  elseif(EXCHCXX_SYCL_TARGET STREQUAL "spir64_x86_64")
    # CPU AoT via opencl-aot. No -march is passed: opencl-aot then leaves
    # CL_CONFIG_CPU_TARGET_ARCH unset and the resulting binary stays portable
    # across x86_64 hosts, which matters for CI runners whose CPU model is not
    # fixed (and may not even be Intel).
    list(APPEND _exchcxx_sycl_compile_opts
      -fsycl-targets=spir64_x86_64
    )
    list(APPEND _exchcxx_sycl_link_opts
      "SHELL:-fsycl-targets=spir64_x86_64"
    )
  endif()

  # Every allowlisted target must map to real flags. Without this an entry that
  # is accepted above but has no branch here would silently produce a JIT build
  # while still reporting "AoT enabled" below.
  if(NOT _exchcxx_sycl_compile_opts)
    message(FATAL_ERROR
      "EXCHCXX_SYCL_TARGET='${EXCHCXX_SYCL_TARGET}' is allowed but has no "
      "AoT flag mapping in ${CMAKE_CURRENT_LIST_FILE}; refusing to fall back "
      "to a JIT build that would be reported as AoT.")
  endif()

  target_compile_options(exchcxx PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:${_exchcxx_sycl_compile_opts}>
  )
  target_link_options(exchcxx PUBLIC
    ${_exchcxx_sycl_link_opts}
  )

  message(STATUS "ExchCXX SYCL AoT enabled for target: ${EXCHCXX_SYCL_TARGET}")
endif()


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-fno-sycl-id-queries-fit-in-int"     EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT )
check_cxx_compiler_flag("-fsycl-device-code-split=per_kernel" EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL )
check_cxx_compiler_flag("-Xsycl-target-frontend \"-fp-model=precise\"" EXCHCXX_HAVE_SYCL_TARGET_FRONTEND_FP_MODEL_PRECISE )


include(CheckLinkerFlag)
check_linker_flag(CXX "-flink-huge-device-code"          EXCHCXX_SYCL_LINK_HUGE_DEVICE_CODE)
check_linker_flag(CXX "-fsycl-max-parallel-link-jobs=4"  EXCHCXX_SYCL_MAX_PARALLEL_LINK_JOBS)


if(EXCHCXX_SYCL_ID_QUERIES_FIT_IN_INT)
  target_compile_options(exchcxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-fno-sycl-id-queries-fit-in-int>
  )
endif()

if(EXCHCXX_SYCL_DEVICE_CODE_SPLIT_PER_KERNEL)
  target_compile_options(exchcxx PUBLIC
    $<$<COMPILE_LANGUAGE:CXX>:-fsycl-device-code-split=per_kernel>
  )
  target_link_options(exchcxx PUBLIC
    $<$<LINK_LANGUAGE:CXX>:-fsycl-device-code-split=per_kernel>
  )
endif()

# Device-side FP model. -Xsycl-target-frontend reaches only the SYCL device
# frontend; the matching host-side -fp-model=precise is applied project-wide in
# the top-level CMakeLists.txt so that libxc and the tests inherit it too.
if(EXCHCXX_HAVE_SYCL_TARGET_FRONTEND_FP_MODEL_PRECISE)
  target_compile_options(exchcxx PRIVATE
    "$<$<COMPILE_LANGUAGE:CXX>:SHELL:-Xsycl-target-frontend -fp-model=precise>"
  )
endif()

if(EXCHCXX_SYCL_LINK_HUGE_DEVICE_CODE)
  target_link_options(exchcxx PRIVATE
    $<$<LINK_LANGUAGE:CXX>:-flink-huge-device-code>
  )
endif()

if(EXCHCXX_SYCL_MAX_PARALLEL_LINK_JOBS)
  target_link_options(exchcxx PRIVATE
    $<$<LINK_LANGUAGE:CXX>:-fsycl-max-parallel-link-jobs=4>
  )
endif()
