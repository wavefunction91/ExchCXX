/**
 * ExchCXX Copyright (c) 2020-2022, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from the U.S. Dept. of Energy). All rights reserved.
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


#include <exchcxx/impl/builtin/kernels/slater_exchange.hpp>
#include <exchcxx/impl/builtin/kernels/vwn3.hpp>
#include <exchcxx/impl/builtin/kernels/vwn_rpa.hpp>
#include <exchcxx/impl/builtin/kernels/pw91_lda.hpp>
#include <exchcxx/impl/builtin/kernels/pw91_lda_mod.hpp>
#include <exchcxx/impl/builtin/kernels/pw91_lda_rpa.hpp>
#include <exchcxx/impl/builtin/kernels/pz81.hpp>
#include <exchcxx/impl/builtin/kernels/pz81_mod.hpp>

#include <exchcxx/impl/builtin/kernels/b88.hpp>
#include <exchcxx/impl/builtin/kernels/lyp.hpp>
#include <exchcxx/impl/builtin/kernels/pbe_x.hpp>
#include <exchcxx/impl/builtin/kernels/rev_pbe_x.hpp>
#include <exchcxx/impl/builtin/kernels/pbe_c.hpp>

#include <exchcxx/impl/builtin/kernels/pbe0.hpp>
#include <exchcxx/impl/builtin/kernels/b3lyp.hpp>

#include <exchcxx/impl/builtin/kernels/scan_x.hpp>
#include <exchcxx/impl/builtin/kernels/scan_c.hpp>
#include <exchcxx/impl/builtin/kernels/r2scan_x.hpp>
#include <exchcxx/impl/builtin/kernels/r2scan_c.hpp>
#include <exchcxx/impl/builtin/kernels/ft98_x.hpp>
#include <exchcxx/impl/builtin/kernels/scanl_c.hpp>
#include <exchcxx/impl/builtin/kernels/scanl_x.hpp>
#include <exchcxx/impl/builtin/kernels/r2scanl_c.hpp>
#include <exchcxx/impl/builtin/kernels/r2scanl_x.hpp>
#include <exchcxx/impl/builtin/kernels/m06_2x_x.hpp>
#include <exchcxx/impl/builtin/kernels/m06_2x_c.hpp>

#include <exchcxx/impl/builtin/kernels/pc07_k.hpp>
#include <exchcxx/impl/builtin/kernels/pc07opt_k.hpp>

#include <exchcxx/impl/builtin/kernels/epc17_1.hpp>
#include <exchcxx/impl/builtin/kernels/epc17_2.hpp>
#include <exchcxx/impl/builtin/kernels/epc18_1.hpp>
#include <exchcxx/impl/builtin/kernels/epc18_2.hpp>

