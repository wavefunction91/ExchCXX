cmake_minimum_required( VERSION 3.17 FATAL_ERROR ) # Require CMake 3.17+
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



# Enable CUDA if desired
if( EXCHCXX_ENABLE_CUDA )
  enable_language( CUDA )
endif()

# Find SYCL if requested
if( EXCHCXX_ENABLE_SYCL )
#  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/oneAPI" )
#  find_package( SYCL REQUIRED )
endif()


## Find LibXC
find_package( Libxc 5.0.0 CONFIG QUIET )

if( ${Libxc_FOUND} )

  message( STATUS "Found: Libxc Version ${Libxc_VERSION}" )
  if( Libxc_VERSION VERSION_GREATER_EQUAL "6.0.0" )
    message( FATAL_ERROR "Libxc version 6+ breaks the API currently used in ExchCXX" )
  endif()

else()

  FetchContent_Declare(
    libxc
    GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
    GIT_TAG        5.0.0
  )

  set( OLD_BUILD_TESTING ${BUILD_TESTING} )
  set( BUILD_TESTING OFF CACHE BOOL "" )

  FetchContent_MakeAvailable( libxc )
  add_library( Libxc::xc ALIAS xc )
  target_include_directories( xc 
    PUBLIC 
      $<BUILD_INTERFACE:${libxc_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/gen_funcidx>
  )

  set( BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" )

endif()
        
add_subdirectory( src )

if( EXCHCXX_ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif()

if( EXCHCXX_ENABLE_BENCHMARK )
  add_subdirectory( performance )
endif()
