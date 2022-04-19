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

  ExchCXX::Kernel::B3LYP,
  ExchCXX::Kernel::PBE0
};

static constexpr std::array string_kernal_pairs = {
    std::pair("SlaterExchange", ExchCXX::Kernel::SlaterExchange),
    std::pair("PBE_X",ExchCXX::Kernel::PBE_X),
    std::pair("PBE_C", ExchCXX::Kernel::PBE_C),
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
    std::pair("PBE0", ExchCXX::Functional::PBE0)
};
