#pragma once

#include <iostream>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace limap {

struct CompareByDouble {
  bool operator()(const std::pair<size_t, double>& a,
                  const std::pair<size_t, double>& b) const {
    // If double values are the same, order by size_t value (ascending).
    if (a.second == b.second) {
      return a.first < b.first;
    }
    // Order primarily by the double value (descending).
    return a.second > b.second;
  }
};

class MaxValuePair {
 public:
  using ElemType = std::pair<size_t, double>;

 private:
  std::unordered_map<size_t, double>
      dataMap;  // Stores the key (size_t) and value (double) pairs for quick
                // access.
  std::set<ElemType, CompareByDouble>
      sortedSet;  // Set to keep elements ordered by double value in descending
                  // order.

 public:
  // Insert or update an element in the structure.
  void insertOrUpdate(size_t key, double value) {
    auto it = dataMap.find(key);
    if (it != dataMap.end()) {
      // If the element already exists, remove the old element from the set
      // before updating.
      if (!sortedSet.erase(ElemType(key, it->second))) {
        throw std::runtime_error("[insertOrUpdate] failed to erase SortedSet");
      }
      it->second = value;  // Update the value in the map.
    } else {
      // Insert the new key-value pair into the map if it doesn't exist.
      dataMap[key] = value;
    }
    // Insert the new or updated element into the sorted set.
    sortedSet.emplace(key, value);
  }

  void UpdateSortedSet(size_t key, double init_value, double value) {
    auto it = dataMap.find(key);
    if (it != dataMap.end()) {
      if (!sortedSet.erase(ElemType(key, init_value))) {
        throw std::runtime_error("[UpdateSortedSet] failed to erase SortedSet");
      }
      sortedSet.emplace(key, value);
    } else {
      throw std::runtime_error("dataMap.at(key) does not exist");
    }
  }

  // Find the maximum double value and its corresponding size_t key.
  ElemType findMaxValue() const {
    if (sortedSet.empty()) {
      throw std::runtime_error(
          "No elements present.");  // Exception if there are no elements.
    }
    return *sortedSet
                .begin();  // Return the element with the maximum double value.
  }

  // Erase an element by its key.
  void erase(size_t key) {
    auto it = dataMap.find(key);
    if (it != dataMap.end()) {
      // If the element exists, remove it from both the map and the set.
      if (!sortedSet.erase(ElemType(key, it->second))) {
        throw std::runtime_error("[erase] failed to erase SortedSet");
      }
      dataMap.erase(it);
    }
  }

  void eraseDataMap(size_t key) {
    auto it = dataMap.find(key);
    if (it != dataMap.end()) {
      dataMap.erase(it);
    }
  }

  void eraseSortedSet(size_t key, double value) {
    // If the element exists, remove it from both the map and the set.
    if (!sortedSet.erase(ElemType(key, value))) {
      throw std::runtime_error("[eraseSortedSet] failed to erase SortedSet");
    }
  }

  bool exist_key(size_t key) { return dataMap.find(key) != dataMap.end(); }

  double& get_value(size_t key) { return dataMap.at(key); }

  bool empty() { return dataMap.empty(); }

  bool check() {
    // check sortedSet and dataMap are consistent
    std::unordered_map<size_t, double> map1;
    for (const auto& pair : sortedSet) {
      map1.emplace(pair.first, pair.second);
    }
    if (map1 != dataMap) {
      return false;
    } else {
      return true;
    }
  }

  void clear() {
    dataMap.clear();
    sortedSet.clear();
  }
};

}  // namespace limap
