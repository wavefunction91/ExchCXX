cmake_minimum_required( VERSION 3.21 FATAL_ERROR ) # Require CMake 3.21+
include(FetchContent)

# Set up project definition + version information
project( ExchCXX VERSION 0.1.0 LANGUAGES C CXX )

# ExchCXX Options
option( EXCHCXX_ENABLE_TESTS     "Enable Unit Tests"            ON  )
option( EXCHCXX_ENABLE_BENCHMARK "Enable Performance Benchmark" OFF )
option( EXCHCXX_ENABLE_CUDA      "Enable Device Code (CUDA)"    OFF )
option( EXCHCXX_ENABLE_HIP       "Enable Device Code (HIP)"     OFF )
option( EXCHCXX_ENABLE_SYCL      "Enable Device Code (SYCL)"    OFF )
option( BUILD_SHARED_LIBS        "Build Shared Libs"            OFF )


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



## Find LibXC
find_package( Libxc 6.1.0 CONFIG QUIET )

if( ${Libxc_FOUND} )

  message( STATUS "Found: Libxc Version ${Libxc_VERSION}" )
  if( Libxc_VERSION VERSION_GREATER_EQUAL "7.0.0" )
    message( FATAL_ERROR "Libxc version 7+ breaks the API currently used in ExchCXX" )
  endif()

else()

  FetchContent_Declare(
    libxc
    GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
    GIT_TAG 6.1.0  # v6.1.0
  )
  set( Libxc_VERSION 6.1.0 )

  set( OLD_BUILD_TESTING ${BUILD_TESTING} )
  set( BUILD_TESTING OFF CACHE BOOL "" FORCE )

  FetchContent_MakeAvailable( libxc )
  add_library( Libxc::xc ALIAS xc )
  target_include_directories( xc 
    PUBLIC 
      $<BUILD_INTERFACE:${libxc_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/gen_funcidx>
  )

  # disable unity builds for libxc
  if (CMAKE_UNITY_BUILD)
    set_target_properties(xc PROPERTIES UNITY_BUILD OFF)
    message(STATUS "Will disable unity-build for Libxc::xc")
  endif()
  
  set( BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" FORCE )

endif()
        
add_subdirectory( src )

# Testing
if( NOT DEFINED EXCHCXX_ENABLE_TESTS )
  set( EXCHCXX_ENABLE_TESTS ON )
endif()

if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND EXCHCXX_ENABLE_TESTS AND BUILD_TESTING )
  add_subdirectory( test )
endif()

if( EXCHCXX_ENABLE_BENCHMARK )
  add_subdirectory( performance )
endif()
