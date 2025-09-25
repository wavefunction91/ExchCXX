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

#include <exchcxx/impl/builtin/kernel_type.hpp>
#include <exchcxx/impl/builtin/util.hpp>
#include <exchcxx/impl/builtin/kernels.hpp>
#include <exchcxx/util/unused.hpp>


namespace ExchCXX {


// LDA interface
UNUSED_INTERFACE_GENERATOR( LDA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( LDA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( LDA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( LDA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )

// GGA interface
UNUSED_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( GGA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( GGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )

// MGGA interface
UNUSED_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc,     const )
UNUSED_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc, const )
UNUSED_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc,     const )
UNUSED_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc, const )


// INC interfaces
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( LDA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )

UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( GGA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )

UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc, const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_inc,     const )
UNUSED_INC_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc, const )


#ifdef EXCHCXX_ENABLE_DEVICE

// LDA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( LDA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )

// GGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( GGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )

// MGGA interface
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_device, const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_device,     const )
UNUSED_DEVICE_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_device, const )


// INC interfaces
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( LDA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( GGA,  VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC,     BuiltinKernel::eval_exc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, EXC_VXC, BuiltinKernel::eval_exc_vxc_inc_device, const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, FXC,     BuiltinKernel::eval_fxc_inc_device,     const )
UNUSED_DEVICE_INC_INTERFACE_GENERATOR( MGGA, VXC_FXC, BuiltinKernel::eval_vxc_fxc_inc_device, const )

#endif
}

