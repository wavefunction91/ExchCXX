#pragma once
#include <map>

namespace ExchCXX {

template <typename Key, typename Val>
class BidirectionalMap{
  std::map<Key, Val> forward_map_;
  std::map<Val, Key> reverse_map_;

public:
  BidirectionalMap(std::map<Key, Val> map) : forward_map_(map){

    for (auto &&v : map) {
      reverse_map_.insert(std::make_pair(v.second, v.first));
    }

  }

  Val val(Key key) { return forward_map_.at(key); }

  Key key(Val val) { return reverse_map_.at(val); }
};

} // namespace ExchCXX