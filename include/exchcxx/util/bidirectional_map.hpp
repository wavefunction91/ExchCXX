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
      if (reverse_map_.find(v.second) != reverse_map_.end()){
        throw std::runtime_error("BidirectionalMap must have unique values");
      }
      reverse_map_.insert(std::make_pair(v.second, v.first));
    }

  }

  Val value(Key key) { return forward_map_.at(key); }

  Key key(Val val) { return reverse_map_.at(val); }

  bool key_exists(Key key){
    return forward_map_.find(key) != forward_map_.end();
  }

  bool value_exists(Val val){
    return reverse_map_.find(val) != reverse_map_.end();
  }

};

} // namespace ExchCXX