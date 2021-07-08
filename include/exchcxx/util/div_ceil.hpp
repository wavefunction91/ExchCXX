#pragma once
#include <type_traits>
#include <cstdlib>
#include <cstdint>

namespace ExchCXX  {
namespace util  {

template <typename Integral1, typename Integral2>
intmax_t div_ceil( Integral1 i, Integral2 j ) {

  intmax_t i_us = i;
  intmax_t j_us = j;

  auto d = std::div(i_us,j_us);
  return d.quot + !!d.rem;

};


}
}

