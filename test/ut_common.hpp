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

#include "catch2/catch.hpp"
#include <exchcxx/exchcxx.hpp>
#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <random>

#include "reference_values.hpp"


enum class TestInterface {
  EXC,
  EXC_VXC,

  EXC_INC,
  EXC_VXC_INC
};

enum class EvalType {
  Regular,
  Small,
  Zero
};


static std::vector<ExchCXX::Kernel> lda_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5,
  ExchCXX::Kernel::PZ81,
  ExchCXX::Kernel::PZ81_MOD,
  ExchCXX::Kernel::PW91_LDA,
  ExchCXX::Kernel::PW91_LDA_MOD,
  ExchCXX::Kernel::PW91_LDA_RPA
};

static std::vector<ExchCXX::Kernel> gga_kernels = {
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,
  ExchCXX::Kernel::revPBE_X,
  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::B3LYP,
  ExchCXX::Kernel::PBE0
};

static std::vector<ExchCXX::Kernel> mgga_kernels = {
  ExchCXX::Kernel::SCAN_X,
  ExchCXX::Kernel::SCAN_C,
  ExchCXX::Kernel::SCANL_C,
  ExchCXX::Kernel::R2SCAN_X,
  ExchCXX::Kernel::R2SCAN_C,
  ExchCXX::Kernel::R2SCANL_X,
  ExchCXX::Kernel::R2SCANL_C,
  ExchCXX::Kernel::FT98_X,
  ExchCXX::Kernel::PC07_K,
  ExchCXX::Kernel::PC07OPT_K
};

static std::vector<ExchCXX::Kernel> builtin_supported_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5,
  ExchCXX::Kernel::PZ81,
  ExchCXX::Kernel::PZ81_MOD,
  ExchCXX::Kernel::PW91_LDA,
  ExchCXX::Kernel::PW91_LDA_MOD,
  ExchCXX::Kernel::PW91_LDA_RPA,

  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,
  ExchCXX::Kernel::revPBE_X,

  ExchCXX::Kernel::B3LYP,
  ExchCXX::Kernel::PBE0,

  ExchCXX::Kernel::SCAN_X, 
  ExchCXX::Kernel::SCAN_C,
  //ExchCXX::Kernel::SCANL_C,
  ExchCXX::Kernel::R2SCAN_X, 
  ExchCXX::Kernel::R2SCAN_C,
  ExchCXX::Kernel::FT98_X,

  ExchCXX::Kernel::PC07_K,
  ExchCXX::Kernel::PC07OPT_K
};

static constexpr std::array string_kernal_pairs = {
    std::pair("SlaterExchange", ExchCXX::Kernel::SlaterExchange),
    std::pair("PBE_X",ExchCXX::Kernel::PBE_X),
    std::pair("PBE_C", ExchCXX::Kernel::PBE_C),
    std::pair("SCAN_X",ExchCXX::Kernel::SCAN_X),
    std::pair("SCAN_C", ExchCXX::Kernel::SCAN_C),
    std::pair("SCANL_C", ExchCXX::Kernel::SCANL_C),
    std::pair("FT98_X",ExchCXX::Kernel::FT98_X),
    std::pair("PC07_K",ExchCXX::Kernel::PC07_K),
    std::pair("PC07OPT_K",ExchCXX::Kernel::PC07OPT_K),
    std::pair("R2SCANL_X",ExchCXX::Kernel::R2SCANL_X),
    std::pair("R2SCANL_C", ExchCXX::Kernel::R2SCANL_C),
    std::pair("R2SCAN_X",ExchCXX::Kernel::R2SCAN_X),
    std::pair("R2SCAN_C", ExchCXX::Kernel::R2SCAN_C),
    std::pair("LYP", ExchCXX::Kernel::LYP),
    std::pair("B3LYP", ExchCXX::Kernel::B3LYP),
    std::pair("PBE0", ExchCXX::Kernel::PBE0),
    std::pair("VWN3", ExchCXX::Kernel::VWN3),
    std::pair("VWN5", ExchCXX::Kernel::VWN5),
    std::pair("PZ81", ExchCXX::Kernel::PZ81),
    std::pair("PZ81_MOD", ExchCXX::Kernel::PZ81_MOD),
    std::pair("PW91_LDA", ExchCXX::Kernel::PW91_LDA),
    std::pair("PW91_LDA_MOD", ExchCXX::Kernel::PW91_LDA_MOD),
    std::pair("PW91_LDA_RPA", ExchCXX::Kernel::PW91_LDA_RPA),
    std::pair("B88", ExchCXX::Kernel::B88)
};

static constexpr std::array string_functional_pairs = {
    std::pair("SVWN3", ExchCXX::Functional::SVWN3),
    std::pair("SVWN5", ExchCXX::Functional::SVWN5),
    std::pair("BLYP", ExchCXX::Functional::BLYP),
    std::pair("B3LYP", ExchCXX::Functional::B3LYP),
    std::pair("PBE", ExchCXX::Functional::PBE),
    std::pair("SCAN", ExchCXX::Functional::SCAN),
    std::pair("R2SCANL", ExchCXX::Functional::R2SCANL),
    std::pair("PBE0", ExchCXX::Functional::PBE0)
};
