cmake_minimum_required( VERSION 3.17 FATAL_ERROR ) # Require CMake 3.17+
include(FetchContent)

# Set up project definition + version information
project( ExchCXX VERSION 0.1.0 LANGUAGES C CXX )

# ExchCXX Options
option( EXCHCXX_ENABLE_TESTS     "Enable Unit Tests"            ON  )
option( EXCHCXX_ENABLE_BENCHMARK "Enable Performance Benchmark" OFF )
option( EXCHCXX_ENABLE_CUDA      "Enable Device Code (CUDA)"    OFF )


# Decided if we're compiling device bindings
if( EXCHCXX_ENABLE_CUDA )
  set( EXCHCXX_ENABLE_DEVICE TRUE CACHE BOOL "Enable Device Code" )
endif()



# Enable CUDA if desired
if( EXCHCXX_ENABLE_CUDA )
  enable_language( CUDA )
endif()






# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
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
    GIT_REPOSITORY https://gitlab.com/eduard1/libxc.git
    GIT_TAG        master
  )

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

endif()
        
add_subdirectory( src  )

if( EXCHCXX_ENABLE_TESTS )
  enable_testing()
  add_subdirectory( test )
endif()

if( EXCHCXX_ENABLE_BENCHMARK )
  add_subdirectory( performance )
endif()
