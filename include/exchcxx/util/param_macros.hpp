/**
 * ExchCXX 
 *
 * Copyright (c) 2020-2024, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). 
 *
 * Portions Copyright (c) Microsoft Corporation.
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * (1) Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 
 * (2) Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * (3) Neither the name of the University of California, Lawrence Berkeley
 * National Laboratory, U.S. Dept. of Energy nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * 
 * You are under no obligation whatsoever to provide any bug fixes, patches,
 * or upgrades to the features, functionality or performance of the source
 * code ("Enhancements") to anyone; however, if you choose to make your
 * Enhancements available either publicly, or directly to Lawrence Berkeley
 * National Laboratory, without imposing a separate written license agreement
 * for such Enhancements, then you hereby grant the following license: a
 * non-exclusive, royalty-free perpetual license to install, use, modify,
 * prepare derivative works, incorporate into other computer software,
 * distribute, and sublicense such enhancements or derivative works thereof,
 * in binary and source code form.
 */

#pragma once

#ifdef EXCHCXX_HAS_CONFIG_H
#include <exchcxx/exchcxx_config.hpp>
#endif


namespace ExchCXX {

  using host_buffer_type       = double*;
  using const_host_buffer_type = const double*;

#ifdef EXCHCXX_ENABLE_DEVICE
  using device_buffer_type       = double*;
  using const_device_buffer_type = const double*;
#endif


}

#define NOTYPE

// LDA Parameters
#define TYPED_LDA_IPARAMS(INTT,BUFFER)      INTT N, BUFFER rho
#define TYPED_LDA_OPARAMS_EXC(BUFFER)       BUFFER eps
#define TYPED_LDA_OPARAMS_VXC(BUFFER)       BUFFER vxc
#define TYPED_LDA_OPARAMS_FXC(BUFFER)       BUFFER fxc
#define TYPED_LDA_OPARAMS_KXC(BUFFER)       BUFFER kxc
#define TYPED_LDA_OPARAMS_EXC_VXC(BUFFER) \
  TYPED_LDA_OPARAMS_EXC(BUFFER), TYPED_LDA_OPARAMS_VXC(BUFFER)
#define TYPED_LDA_OPARAMS_VXC_FXC(BUFFER) \
  TYPED_LDA_OPARAMS_VXC(BUFFER), TYPED_LDA_OPARAMS_FXC(BUFFER)

#define LDA_IPARAMS         TYPED_LDA_IPARAMS(int,const_host_buffer_type)
#define LDA_OPARAMS_EXC     TYPED_LDA_OPARAMS_EXC(host_buffer_type)
#define LDA_OPARAMS_VXC     TYPED_LDA_OPARAMS_VXC(host_buffer_type)
#define LDA_OPARAMS_FXC     TYPED_LDA_OPARAMS_FXC(host_buffer_type)
#define LDA_OPARAMS_KXC     TYPED_LDA_OPARAMS_KXC(host_buffer_type)
#define LDA_OPARAMS_EXC_VXC TYPED_LDA_OPARAMS_EXC_VXC(host_buffer_type)
#define LDA_OPARAMS_VXC_FXC TYPED_LDA_OPARAMS_VXC_FXC(host_buffer_type)

#define DEV_LDA_IPARAMS         TYPED_LDA_IPARAMS(int,const_device_buffer_type)
#define DEV_LDA_OPARAMS_EXC     TYPED_LDA_OPARAMS_EXC(device_buffer_type)
#define DEV_LDA_OPARAMS_VXC     TYPED_LDA_OPARAMS_VXC(device_buffer_type)
#define DEV_LDA_OPARAMS_FXC     TYPED_LDA_OPARAMS_FXC(device_buffer_type)
#define DEV_LDA_OPARAMS_KXC     TYPED_LDA_OPARAMS_KXC(device_buffer_type)
#define DEV_LDA_OPARAMS_EXC_VXC TYPED_LDA_OPARAMS_EXC_VXC(device_buffer_type)
#define DEV_LDA_OPARAMS_VXC_FXC TYPED_LDA_OPARAMS_VXC_FXC(device_buffer_type)

#define LDA_IPARAMS_NOTYPE         TYPED_LDA_IPARAMS(NOTYPE,NOTYPE)
#define LDA_OPARAMS_EXC_NOTYPE     TYPED_LDA_OPARAMS_EXC(NOTYPE)
#define LDA_OPARAMS_VXC_NOTYPE     TYPED_LDA_OPARAMS_VXC(NOTYPE)
#define LDA_OPARAMS_FXC_NOTYPE     TYPED_LDA_OPARAMS_FXC(NOTYPE)
#define LDA_OPARAMS_KXC_NOTYPE     TYPED_LDA_OPARAMS_KXC(NOTYPE)
#define LDA_OPARAMS_EXC_VXC_NOTYPE TYPED_LDA_OPARAMS_EXC_VXC(NOTYPE)
#define LDA_OPARAMS_VXC_FXC_NOTYPE TYPED_LDA_OPARAMS_VXC_FXC(NOTYPE)



// GGA Parameters
#define TYPED_GGA_IPARAMS(INTT,BUFFER)     INTT N, BUFFER rho, BUFFER sigma
#define TYPED_GGA_OPARAMS_EXC(BUFFER) BUFFER eps
#define TYPED_GGA_OPARAMS_VXC(BUFFER) BUFFER vrho,   BUFFER vsigma
#define TYPED_GGA_OPARAMS_FXC(BUFFER) \
  BUFFER v2rho2, BUFFER v2rhosigma, BUFFER v2sigma2
#define TYPED_GGA_OPARAMS_KXC(BUFFER) \
  BUFFER v3rho3, BUFFER v3rho2sigma, BUFFER v3rhosigma2, BUFFER v3sigma3 

#define TYPED_GGA_OPARAMS_EXC_VXC(BUFFER) \
  TYPED_GGA_OPARAMS_EXC(BUFFER), TYPED_GGA_OPARAMS_VXC(BUFFER)
#define TYPED_GGA_OPARAMS_VXC_FXC(BUFFER) \
  TYPED_GGA_OPARAMS_VXC(BUFFER), TYPED_GGA_OPARAMS_FXC(BUFFER)



#define GGA_IPARAMS         TYPED_GGA_IPARAMS(int,const_host_buffer_type)
#define GGA_OPARAMS_EXC     TYPED_GGA_OPARAMS_EXC(host_buffer_type)
#define GGA_OPARAMS_VXC     TYPED_GGA_OPARAMS_VXC(host_buffer_type)
#define GGA_OPARAMS_FXC     TYPED_GGA_OPARAMS_FXC(host_buffer_type)
#define GGA_OPARAMS_KXC     TYPED_GGA_OPARAMS_KXC(host_buffer_type)
#define GGA_OPARAMS_EXC_VXC TYPED_GGA_OPARAMS_EXC_VXC(host_buffer_type)
#define GGA_OPARAMS_VXC_FXC TYPED_GGA_OPARAMS_VXC_FXC(host_buffer_type)

#define DEV_GGA_IPARAMS         TYPED_GGA_IPARAMS(int,const_device_buffer_type)
#define DEV_GGA_OPARAMS_EXC     TYPED_GGA_OPARAMS_EXC(device_buffer_type)
#define DEV_GGA_OPARAMS_VXC     TYPED_GGA_OPARAMS_VXC(device_buffer_type)
#define DEV_GGA_OPARAMS_FXC     TYPED_GGA_OPARAMS_FXC(device_buffer_type)
#define DEV_GGA_OPARAMS_KXC     TYPED_GGA_OPARAMS_KXC(device_buffer_type)
#define DEV_GGA_OPARAMS_EXC_VXC TYPED_GGA_OPARAMS_EXC_VXC(device_buffer_type)
#define DEV_GGA_OPARAMS_VXC_FXC TYPED_GGA_OPARAMS_VXC_FXC(device_buffer_type)

#define GGA_IPARAMS_NOTYPE         TYPED_GGA_IPARAMS(NOTYPE,NOTYPE)
#define GGA_OPARAMS_EXC_NOTYPE     TYPED_GGA_OPARAMS_EXC(NOTYPE)
#define GGA_OPARAMS_VXC_NOTYPE     TYPED_GGA_OPARAMS_VXC(NOTYPE)
#define GGA_OPARAMS_FXC_NOTYPE     TYPED_GGA_OPARAMS_FXC(NOTYPE)
#define GGA_OPARAMS_KXC_NOTYPE     TYPED_GGA_OPARAMS_KXC(NOTYPE)
#define GGA_OPARAMS_EXC_VXC_NOTYPE TYPED_GGA_OPARAMS_EXC_VXC(NOTYPE)
#define GGA_OPARAMS_VXC_FXC_NOTYPE TYPED_GGA_OPARAMS_VXC_FXC(NOTYPE)


// MGGA Parameters
#define TYPED_MGGA_IPARAMS(INTT,BUFFER) \
  INTT N, BUFFER rho, BUFFER sigma, BUFFER lapl, BUFFER tau
#define TYPED_MGGA_OPARAMS_EXC(BUFFER) BUFFER eps
#define TYPED_MGGA_OPARAMS_VXC(BUFFER) \
  BUFFER vrho, BUFFER vsigma, BUFFER vlapl, BUFFER vtau
#define TYPED_MGGA_OPARAMS_FXC(BUFFER) \
  BUFFER v2rho2, BUFFER v2rhosigma, BUFFER v2rholapl, BUFFER v2rhotau, \
  BUFFER v2sigma2, BUFFER v2sigmalapl, BUFFER v2sigmatau, BUFFER v2lapl2, \
  BUFFER v2lapltau, BUFFER v2tau2

#define TYPED_MGGA_OPARAMS_EXC_VXC(BUFFER) \
  TYPED_MGGA_OPARAMS_EXC(BUFFER), TYPED_MGGA_OPARAMS_VXC(BUFFER)
#define TYPED_MGGA_OPARAMS_VXC_FXC(BUFFER) \
  TYPED_MGGA_OPARAMS_VXC(BUFFER), TYPED_MGGA_OPARAMS_FXC(BUFFER)

#define MGGA_IPARAMS         TYPED_MGGA_IPARAMS(int,const_host_buffer_type)
#define MGGA_OPARAMS_EXC     TYPED_MGGA_OPARAMS_EXC(host_buffer_type)
#define MGGA_OPARAMS_VXC     TYPED_MGGA_OPARAMS_VXC(host_buffer_type)
#define MGGA_OPARAMS_FXC     TYPED_MGGA_OPARAMS_FXC(host_buffer_type)
#define MGGA_OPARAMS_EXC_VXC TYPED_MGGA_OPARAMS_EXC_VXC(host_buffer_type)
#define MGGA_OPARAMS_VXC_FXC TYPED_MGGA_OPARAMS_VXC_FXC(host_buffer_type)

#define DEV_MGGA_IPARAMS         TYPED_MGGA_IPARAMS(int,const_device_buffer_type)
#define DEV_MGGA_OPARAMS_EXC     TYPED_MGGA_OPARAMS_EXC(device_buffer_type)
#define DEV_MGGA_OPARAMS_VXC     TYPED_MGGA_OPARAMS_VXC(device_buffer_type)
#define DEV_MGGA_OPARAMS_FXC     TYPED_MGGA_OPARAMS_FXC(device_buffer_type)
#define DEV_MGGA_OPARAMS_EXC_VXC TYPED_MGGA_OPARAMS_EXC_VXC(device_buffer_type)
#define DEV_MGGA_OPARAMS_VXC_FXC TYPED_MGGA_OPARAMS_VXC_FXC(device_buffer_type)

#define MGGA_IPARAMS_NOTYPE         TYPED_MGGA_IPARAMS(NOTYPE,NOTYPE)
#define MGGA_OPARAMS_EXC_NOTYPE     TYPED_MGGA_OPARAMS_EXC(NOTYPE)
#define MGGA_OPARAMS_VXC_NOTYPE     TYPED_MGGA_OPARAMS_VXC(NOTYPE)
#define MGGA_OPARAMS_FXC_NOTYPE     TYPED_MGGA_OPARAMS_FXC(NOTYPE)
#define MGGA_OPARAMS_EXC_VXC_NOTYPE TYPED_MGGA_OPARAMS_EXC_VXC(NOTYPE)
#define MGGA_OPARAMS_VXC_FXC_NOTYPE TYPED_MGGA_OPARAMS_VXC_FXC(NOTYPE)
