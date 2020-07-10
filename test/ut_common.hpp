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


static constexpr std::array lda_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5
};

static constexpr std::array gga_kernels = {
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,
  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::B3LYP,
  ExchCXX::Kernel::PBE0
};

static constexpr std::array builtin_supported_kernels = {
  ExchCXX::Kernel::SlaterExchange,
  ExchCXX::Kernel::VWN3,
  ExchCXX::Kernel::VWN5,

  ExchCXX::Kernel::B88,
  ExchCXX::Kernel::LYP,
  ExchCXX::Kernel::PBE_X,
  ExchCXX::Kernel::PBE_C,

  ExchCXX::Kernel::B3LYP,
  ExchCXX::Kernel::PBE0
};


