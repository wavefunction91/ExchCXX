cmake_minimum_required( VERSION 3.21 FATAL_ERROR ) # Require CMake 3.21+
include(FetchContent)

# Set up project definition + version information
project( ExchCXX VERSION 1.0.0 LANGUAGES C CXX )

# ExchCXX Options
option( EXCHCXX_ENABLE_BENCHMARK "Enable Performance Benchmark" OFF )
option( EXCHCXX_ENABLE_CUDA      "Enable Device Code (CUDA)"    OFF )
option( EXCHCXX_ENABLE_HIP       "Enable Device Code (HIP)"     OFF )
option( EXCHCXX_ENABLE_SYCL      "Enable Device Code (SYCL)"    OFF )
option( EXCHCXX_ENABLE_LIBXC     "Enable Libxc Backend"         ON  )
option( BUILD_SHARED_LIBS        "Build Shared Libs"            OFF )

# --- SYCL Architecture Options (mutually exclusive) ---
option(EXCHCXX_SYCL_INTEL_PVC_ARCH   "Target Intel PVC (Ponte Vecchio)" OFF)
option(EXCHCXX_SYCL_NV_SM_80_ARCH    "Target NVIDIA SM_80 (A100)"       OFF)
option(EXCHCXX_SYCL_NV_SM_90_ARCH    "Target NVIDIA SM_90 (H100/H200)"  OFF)
option(EXCHCXX_SYCL_AMD_GFX90A_ARCH  "Target AMD gfx90a (MI200)"        OFF)
option(EXCHCXX_SYCL_AMD_GFX942_ARCH  "Target AMD gfx942 (MI300A/X)"     OFF)


# Decided if we're compiling device bindings
if( EXCHCXX_ENABLE_CUDA OR EXCHCXX_ENABLE_SYCL OR EXCHCXX_ENABLE_HIP )
  set( EXCHCXX_ENABLE_DEVICE TRUE CACHE BOOL "Enable Device Code" )
endif()

# Die if both CUDA and SYCL are enabled
if( EXCHCXX_ENABLE_CUDA AND EXCHCXX_ENABLE_SYCL )
  message( FATAL_ERROR "CUDA / SYCL bindings are mutually exclusive" )
endif()
if( EXCHCXX_ENABLE_CUDA AND EXCHCXX_ENABLE_HIP )
  message( FATAL_ERROR "CUDA / HIP bindings are mutually exclusive" )
endif()
if( EXCHCXX_ENABLE_SYCL AND EXCHCXX_ENABLE_HIP )
  message( FATAL_ERROR "SYCL / HIP bindings are mutually exclusive" )
endif()


if(EXCHCXX_ENABLE_SYCL)
  # Allowed alias strings (as consumed by -fsycl-targets)
  set(_EXCHCXX_SYCL_ALLOWED
    intel_gpu_pvc
    nvidia_gpu_sm_80
    nvidia_gpu_sm_90
    amd_gpu_gfx90a
    amd_gpu_gfx942
  )

  # Collect selections from booleans
  set(_EXCHCXX_SYCL_FROM_BOOLEANS "")
  if(EXCHCXX_SYCL_INTEL_PVC_ARCH)
    list(APPEND _EXCHCXX_SYCL_FROM_BOOLEANS intel_gpu_pvc)
  endif()
  if(EXCHCXX_SYCL_NV_SM_80_ARCH)
    list(APPEND _EXCHCXX_SYCL_FROM_BOOLEANS nvidia_gpu_sm_80)
  endif()
  if(EXCHCXX_SYCL_NV_SM_90_ARCH)
    list(APPEND _EXCHCXX_SYCL_FROM_BOOLEANS nvidia_gpu_sm_90)
  endif()
  if(EXCHCXX_SYCL_AMD_GFX90A_ARCH)
    list(APPEND _EXCHCXX_SYCL_FROM_BOOLEANS amd_gpu_gfx90a)
  endif()
  if(EXCHCXX_SYCL_AMD_GFX942_ARCH)
    list(APPEND _EXCHCXX_SYCL_FROM_BOOLEANS amd_gpu_gfx942)
  endif()
  
  # Resolve the selected target (prefer booleans; fall back to user-provided EXCHCXX_SYCL_TARGET)
  set(_EXCHCXX_SELECTED "")
  list(LENGTH _EXCHCXX_SYCL_FROM_BOOLEANS _exchcxx_bool_count)
  if(_exchcxx_bool_count GREATER 1)
    message(FATAL_ERROR
      "Multiple SYCL arch options enabled (${_EXCHCXX_SYCL_FROM_BOOLEANS}). Enable exactly one.")
  elseif(_exchcxx_bool_count EQUAL 1)
    list(GET _EXCHCXX_SYCL_FROM_BOOLEANS 0 _EXCHCXX_SELECTED)
  elseif(DEFINED EXCHCXX_SYCL_TARGET AND NOT EXCHCXX_SYCL_TARGET STREQUAL "")
    set(_EXCHCXX_SELECTED "${EXCHCXX_SYCL_TARGET}")  # backward compat
  endif()

  # Validate & cache the resolved target (or leave empty for JIT)
  if(NOT _EXCHCXX_SELECTED STREQUAL "")
    list(FIND _EXCHCXX_SYCL_ALLOWED "${_EXCHCXX_SELECTED}" _exchcxx_sycl_idx)
    if(_exchcxx_sycl_idx EQUAL -1)
      message(FATAL_ERROR
        "Invalid SYCL AoT target '${_EXCHCXX_SELECTED}'. "
        "Allowed values: ${_EXCHCXX_SYCL_ALLOWED}")
    endif()
    set(EXCHCXX_SYCL_TARGET "${_EXCHCXX_SELECTED}" CACHE STRING "Resolved SYCL target alias" FORCE)
  else()
    set(EXCHCXX_SYCL_TARGET "" CACHE STRING "Resolved SYCL target alias (empty => JIT)" FORCE)
    message(STATUS "ExchCXX SYCL AoT disabled (no target selected) — will JIT at runtime.")
  endif()
  
endif(EXCHCXX_ENABLE_SYCL)


# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()

# Populate BUILD_TESTING prior to dependencies to avoid clash
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()



# Enable CUDA if desired
if( EXCHCXX_ENABLE_CUDA )
  enable_language( CUDA )
endif()

if( EXCHCXX_ENABLE_HIP )
  enable_language( HIP )
endif()


if(EXCHCXX_ENABLE_LIBXC)

## Find LibXC
find_package( Libxc 7.0.0 CONFIG QUIET )

if( ${Libxc_FOUND} )

  message( STATUS "Found: Libxc Version ${Libxc_VERSION}" )
  if( Libxc_VERSION VERSION_GREATER_EQUAL "8.0.0" )
    message( FATAL_ERROR "Libxc version 8+ breaks the API currently used in ExchCXX" )
  endif()

else()

  option( FETCHCONTENT_LIBXC_GIT_SHALLOW      "Whether to use GIT_SHALLOW for FetchContent'ing libxc"     ON )
  
  FetchContent_Declare(
    libxc
    GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
    # if pinning to specific SHA change the FETCHCONTENT_LIBXC_GIT_SHALLOW default to OFF, https://cmake.org/cmake/help/latest/module/ExternalProject.html#git
    GIT_TAG 7.0.0
    GIT_SHALLOW ${FETCHCONTENT_LIBXC_GIT_SHALLOW}
    PATCH_COMMAND sed -i -e "s/p->info->family != XC_KINETIC/p->info->kind != XC_KINETIC/g" src/work_mgga_inc.c
  )
  set( Libxc_VERSION 7.0.0 )

  set( OLD_BUILD_TESTING ${BUILD_TESTING} )
  set( BUILD_TESTING OFF CACHE BOOL "" FORCE )

  FetchContent_MakeAvailable( libxc )
  add_library( Libxc::xc ALIAS xc )
  target_include_directories( xc 
    PUBLIC 
      $<BUILD_INTERFACE:${libxc_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}>
  )
  message( STATUS "Libxc Source: ${libxc_SOURCE_DIR}" )

  # disable unity builds for libxc
  if (CMAKE_UNITY_BUILD)
    set_target_properties(xc PROPERTIES UNITY_BUILD OFF)
    message(STATUS "Will disable unity-build for Libxc::xc")
  endif()
  
  set( BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" FORCE )

endif()

else( EXCHCXX_ENABLE_LIBXC )
  set( Libxc_FOUND FALSE )
endif( EXCHCXX_ENABLE_LIBXC )
        
add_subdirectory( src )

# Testing
if( NOT DEFINED EXCHCXX_ENABLE_TESTS )
  set( EXCHCXX_ENABLE_TESTS ON )
endif()

if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND EXCHCXX_ENABLE_TESTS AND BUILD_TESTING )
  if(${Libxc_FOUND} OR NOT EXCHCXX_ENABLE_LIBXC)
    message(WARNING "ExchCXX Unit Tests Require Local Patch of Libxc - Disabling")
  else()
    add_subdirectory( test )
  endif()
endif()

if( EXCHCXX_ENABLE_BENCHMARK AND EXCHCXX_ENABLE_LIBXC )
  add_subdirectory( performance )
endif()
