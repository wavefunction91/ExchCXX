#   FindLibXC.cmake
#
#   Finds the LibXC library
#
#   This module will define the following variables:
#
#     LIBXC_FOUND         - System has found LibXC installation
#     LIBXC_INCLUDE_DIR   - Location of LibXC headers
#     LIBXC_LIBRARIES     - LibXC libraries
#
#   This module will export the following targets if LIBXC_FOUND
#
#     LibXC::libxc
#
#    Proper usage:
#
#      project( TEST_FIND_LIBXC C )
#      find_package( LibXC )
#
#      if( LIBXC_FOUND )
#        add_executable( test test.cxx )
#        target_link_libraries( test LibXC::libxc )
#      endif()
#
#   This module will use the following variables to change
#   default behaviour if set
#
#     libxc_PREFIX
#     libxc_INCLUDE_DIR
#     libxc_LIBRARY_DIR
#     libxc_LIBRARIES
#
#==================================================================
#   Copyright (c) 2018 The Regents of the University of California,
#   through Lawrence Berkeley National Laboratory.  
#
#   Author: David Williams-Young
#   
#   This file is part of cmake-modules. All rights reserved.
#   
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#   
#   (1) Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#   (2) Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#   (3) Neither the name of the University of California, Lawrence Berkeley
#   National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#   
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#   
#   You are under no obligation whatsoever to provide any bug fixes, patches, or
#   upgrades to the features, functionality or performance of the source code
#   ("Enhancements") to anyone; however, if you choose to make your Enhancements
#   available either publicly, or directly to Lawrence Berkeley National
#   Laboratory, without imposing a separate written license agreement for such
#   Enhancements, then you hereby grant the following license: a non-exclusive,
#   royalty-free perpetual license to install, use, modify, prepare derivative
#   works, incorporate into other computer software, distribute, and sublicense
#   such enhancements or derivative works thereof, in binary and source code form.
#
#==================================================================

cmake_minimum_required( VERSION 3.11 ) # Require CMake 3.11+
include(FindPackageHandleStandardArgs)

# Set up some auxillary vars if hints have been set

if( libxc_PREFIX AND NOT libxc_INCLUDE_DIR )
  set( libxc_INCLUDE_DIR ${libxc_PREFIX}/include )
endif()


if( libxc_PREFIX AND NOT libxc_LIBRARY_DIR )
  set( libxc_LIBRARY_DIR 
    ${libxc_PREFIX}/lib 
    ${libxc_PREFIX}/lib32 
    ${libxc_PREFIX}/lib64 
  )
endif()







# Try to find the header
find_path( LIBXC_INCLUDE_DIR 
  NAMES xc.h
  HINTS ${libxc_PREFIX}
  PATHS ${libxc_INCLUDE_DIR}
  PATH_SUFFIXES include
  DOC "Location of LibXC header"
)



# Try to serial find libraries if not already set
if( NOT libxc_LIBRARIES )

  find_library( LIBXC_LIBRARIES
    NAMES xc 
    HINTS ${libxc_PREFIX}
    PATHS ${libxc_LIBRARY_DIR}
    PATH_SUFFIXES lib lib64 lib32
    DOC "LibXC Library"
  )

else()

  # fiXME
  set( LIBXC_LIBRARIES ${libxc_LIBRARIES} )

endif()



# Check version
if( EXISTS ${LIBXC_INCLUDE_DIR}/xc_version.h )
  set( version_pattern 
  "^#define[\t ]+XC_(MAJOR|MINOR|MICRO)_VERSION[\t ]+([0-9\\.]+)$"
  )
  file( STRINGS ${LIBXC_INCLUDE_DIR}/xc_version.h libxc_version
        REGEX ${version_pattern} )
  
  foreach( match ${libxc_version} )
  
    if(LIBXC_VERSION_STRING)
      set(LIBXC_VERSION_STRING "${LIBXC_VERSION_STRING}.")
    endif()
  
    string(REGEX REPLACE ${version_pattern} 
      "${LIBXC_VERSION_STRING}\\2" 
      LIBXC_VERSION_STRING ${match}
    )
  
    set(LIBXC_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
  
  endforeach()
  
  unset( libxc_version )
  unset( version_pattern )
endif()



# Determine if we've found LibXC
mark_as_advanced( LIBXC_FOUND LIBXC_INCLUDE_DIR LIBXC_LIBRARIES )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( LIBXC
  REQUIRED_VARS LIBXC_LIBRARIES LIBXC_INCLUDE_DIR
  VERSION_VAR LIBXC_VERSION_STRING
)

# Export target
if( LIBXC_FOUND AND NOT TARGET LibXC::libxc )

  add_library( LibXC::libxc INTERFACE IMPORTED )
  set_target_properties( LibXC::libxc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES      "${LIBXC_LIBRARIES}" 
  )

endif()
