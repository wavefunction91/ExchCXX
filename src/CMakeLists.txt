set( EXCHCXX_SOURCES 
     xc_kernel.cxx 
     xc_functional.cxx 
     libxc.cxx 
     boilerplate.cxx 
     builtin.cxx 
     builtin_interface.cxx 
     builtin_kernel.cxx 
)

if( EXCHCXX_ENABLE_DEVICE )
  list( APPEND EXCHCXX_SOURCES xc_functional_device.cu libxc_device.cxx builtin.cu )
endif()

add_library( exchcxx STATIC ${EXCHCXX_SOURCES} )


# TARGET properties

target_compile_features( exchcxx PUBLIC cxx_std_14 cuda_std_14 )
target_link_libraries( exchcxx PUBLIC Libxc::xc )

if( EXCHCXX_ENABLE_DEVICE )
  find_package( CUDAToolkit REQUIRED )
  target_compile_definitions( exchcxx PUBLIC EXCHCXX_ENABLE_DEVICE=1 )
  target_link_libraries( exchcxx PUBLIC CUDA::cudart )
endif()


target_include_directories( exchcxx
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

target_compile_options( exchcxx INTERFACE "-Wall" "-Wextra" "-Wpedantic" "-Wnon-virtual-dtor")


# INSTALL rules
add_library( ExchCXX::exchcxx ALIAS exchcxx )


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
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

if( ExchCXX_ENABLE_FORTRAN )
  install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
endif()

# Export target to script
install(EXPORT exchcxx-targets
  FILE
    ExchCXXTargets.cmake
  NAMESPACE
    ExchCXX::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)


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
install(DIRECTORY 
  ${PROJECT_SOURCE_DIR}/cmake/ 
  DESTINATION ${INSTALL_CONFIGDIR} 
  FILES_MATCHING PATTERN "*.cmake"
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ExchCXXConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
