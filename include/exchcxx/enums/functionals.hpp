#pragma once

#include <exchcxx/bidirectional_map.hpp>

namespace ExchCXX {

enum class Functional {
  SVWN3,
  SVWN5,
  BLYP,
  B3LYP,
  PBE0,
};

extern BidirectionalMap<std::string, Functional> functional_map;

std::ostream &operator<<(std::ostream &out, Functional functional);

}
