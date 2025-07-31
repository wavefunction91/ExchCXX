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

#include <exchcxx/util/bidirectional_map.hpp>

namespace ExchCXX {

enum class Functional {
  SVWN3,
  SVWN5,
  BLYP,
  B3LYP,
  PBE,
  revPBE,
  PBE0,
  SCAN,
  R2SCAN,
  R2SCANL,
  M062X,
  PKZB,
  EPC17_1,
  EPC17_2,
  EPC18_1,
  EPC18_2,

  B97D,
  B97D3ZERO,
  CAMB3LYP,
  LDA,
  M06L,
  SCAN0,
  SPW92,
  TPSS,
  TPSSh,
  TPSS0,
  VWN3,
  VWN5,
  LRCwPBE,
  LRCwPBEh,
  BP86,
  HSE03,
  HSE06,
  revB3LYP,
  revPBE0,
  revTPSS,
  revTPSSh,
  PW91,
  mBEEF,
  B3PW91,
  O3LYP,
  OLYP,
  OPBE,
  MPW1K,
  RPBE,
  B88,
  MPW91,
  RSCAN,
  TUNEDCAMB3LYP,
  wB97,
  wB97X,
  wB97XD,
  wB97XD3,
  LCwPBE,
  X3LYP,
  XLYP,
  BHANDH,
  BMK,
  BP86VWN,
  PW86B95,
  PW86PBE,
  R2SCAN0,
  R2SCANh,
  R2SCAN50,
  M05,
  M06,
  M08HX,
  M08SO,
  M052X,
  M06SX,
  CF22D,
  SOGGA11X,
  M06HF,
  M11,
  MN12L,
  MN12SX,
  MN15,
  MN15L,
  revM06L,
};

extern BidirectionalMap<std::string, Functional> functional_map;

std::ostream &operator<<(std::ostream &out, Functional functional);

}
