set( EXCHCXX_SOURCES
     xc_kernel.cxx
     xc_functional.cxx
     boilerplate.cxx
     builtin.cxx
     builtin_interface.cxx
     builtin_kernel.cxx
)
if(EXCHCXX_ENABLE_LIBXC)
  list(APPEND EXCHCXX_SOURCES libxc.cxx)
endif()
message(STATUS ${EXCHCXX_SOURCES})

add_library( exchcxx ${EXCHCXX_SOURCES} )

# TARGET properties

target_compile_features( exchcxx PUBLIC cxx_std_17 )
if(MSVC)
  target_compile_definitions( exchcxx PUBLIC _USE_MATH_DEFINES )
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")  # clang-cl
    target_compile_options( exchcxx PUBLIC
      -Wno-c++98-compat
      -Wno-c++98-compat-local-type-template-args
      -Wno-c++98-compat-pedantic
      -Wno-deprecated-declarations
      -Wno-exit-time-destructors
      -Wno-extra-semi
      -Wno-extra-semi-stmt
      -Wno-float-conversion
      -Wno-float-equal
      -Wno-global-constructors
      -Wno-inconsistent-missing-destructor-override
      -Wno-missing-prototypes
      -Wno-missing-variable-declarations
      -Wno-old-style-cast
      -Wno-pre-c++14-compat
      -Wno-pre-c++17-compat
      -Wno-reserved-macro-identifier
      -Wno-sign-conversion
      -Wno-suggest-destructor-override
      -Wno-switch-enum
      -Wno-undefined-func-template
      -Wno-unsafe-buffer-usage
      -Wno-unsafe-buffer-usage-in-libc-call
      -Wno-unused-macros
      -Wno-unused-parameter
      -Wno-unused-variable
    )
  else()  # cl
    target_compile_options( exchcxx PUBLIC 
      /wd4003  # not enough actual parameters for macro 'identifier'
      /wd4061  # enumerator 'enumerator' in switch of enum 'enumeration' is not explicitly handled by a case label
      /wd4244  # conversion from 'type1' to 'type2', possible loss of data
      /wd4266  # no override available for virtual member function from base 'class'; function is hidden
      /wd4267  # conversion from 'size_t' to 'type', possible loss of data
      /wd4365  # conversion from 'type1' to 'type2', signed/unsigned mismatch
      /wd4514  # unreferenced inline function has been removed
      /wd4267  # conversion from 'size_t' to 'type', possible loss of data
      /wd4710  # function 'function' not inlined
      /wd4711  # function 'function' selected for automatic inline expansion
      /wd4800  # forcing value to bool 'true' or 'false' (performance warning)
      /wd4820  # 'bytes' bytes padding added after construct 'member_name'
      /wd5246  # C-style struct initialization may not be portable
    )
  endif()
endif()
if( EXCHCXX_ENABLE_LIBXC )
  target_link_libraries( exchcxx PUBLIC Libxc::xc )
endif()

# Generate exchcxx_config.hpp
configure_file(
  ${PROJECT_SOURCE_DIR}/include/exchcxx/exchcxx_config.hpp.in
  ${PROJECT_BINARY_DIR}/include/exchcxx/exchcxx_config.hpp
)
target_compile_definitions( exchcxx PUBLIC "EXCHCXX_HAS_CONFIG_H=1" )



# Device specific
if( EXCHCXX_ENABLE_CUDA )
  include( cuda/exchcxx_cuda.cmake )
endif()

if( EXCHCXX_ENABLE_HIP )
  include( hip/exchcxx_hip.cmake )
endif()

if( EXCHCXX_ENABLE_SYCL )
  include( sycl/exchcxx_sycl.cmake )
endif()



target_include_directories( exchcxx
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
)

include( CheckCXXCompilerFlag )
check_cxx_compiler_flag( -Wall              EXCHCXX_CXX_HAS_WALL              )
check_cxx_compiler_flag( -Wextra            EXCHCXX_CXX_HAS_WEXTRA            )
check_cxx_compiler_flag( -Wpedantic         EXCHCXX_CXX_HAS_WPEDANTIC         )
check_cxx_compiler_flag( -Wnon-virtual-dtor EXCHCXX_CXX_HAS_WNON_VIRTUAL_DTOR )
check_cxx_compiler_flag( -Wshadow           EXCHCXX_CXX_HAS_WSHADOW           )

if( EXCHCXX_CXX_HAS_WALL )
  target_compile_options( exchcxx PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wall> )
endif()

if( EXCHCXX_CXX_HAS_WEXTRA )
  target_compile_options( exchcxx PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wextra> )
endif()

if( EXCHCXX_CXX_HAS_WPEDANTIC )
  target_compile_options( exchcxx PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wpedantic> )
endif()

if( EXCHCXX_CXX_HAS_WNON_VIRTUAL_DTOR )
  target_compile_options( exchcxx PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wnon-virtual-dtor -Werror=non-virtual-dtor> )
endif()

if( EXCHCXX_CXX_HAS_WSHADOW )
  target_compile_options( exchcxx PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -Wshadow -Werror=shadow> )
endif()

# INSTALL rules
add_library( ExchCXX::ExchCXX ALIAS exchcxx )


include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/ExchCXX )

# Targets
install(TARGETS exchcxx
    EXPORT exchcxx-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)


set_target_properties( exchcxx PROPERTIES EXPORT_NAME ExchCXX )

# Install Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hpp" )
install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hpp" )

# Export target to script
install(EXPORT exchcxx-targets
  FILE
    ExchCXXTargets.cmake
  NAMESPACE
    ExchCXX::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

# Export build-tree targets also (to be usable by e.g. FetchContent)
export(EXPORT exchcxx-targets
    NAMESPACE ExchCXX::
    FILE "${PROJECT_BINARY_DIR}/ExchCXXTargets.cmake")

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Setup ExchCXXConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/ExchCXXConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
#install(DIRECTORY
#  ${PROJECT_SOURCE_DIR}/cmake/
#  DESTINATION ${INSTALL_CONFIGDIR}
#  FILES_MATCHING PATTERN "*.cmake"
#)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
