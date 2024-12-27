#pragma once

#include <vector>

namespace limap {

namespace math {

// Determine median value in vector. Returns NaN for empty vectors.
template <typename T>
double Median(const std::vector<T>& elems);

// Determine mean value in a vector.
template <typename T>
double Mean(const std::vector<T>& elems);

// Compute the n-th percentile in the given sequence.
template <typename T>
T Percentile(const std::vector<T>& elems, double p);

////////////////////////////////////////////////////////////////////////////////
// Implement
////////////////////////////////////////////////////////////////////////////////

template <typename T>
double Median(const std::vector<T>& elems) {
  CHECK(!elems.empty());

  const size_t mid_idx = elems.size() / 2;

  std::vector<T> ordered_elems = elems;
  std::nth_element(ordered_elems.begin(), ordered_elems.begin() + mid_idx,
                   ordered_elems.end());

  if (elems.size() % 2 == 0) {
    const T mid_element1 = ordered_elems[mid_idx];
    const T mid_element2 = *std::max_element(ordered_elems.begin(),
                                             ordered_elems.begin() + mid_idx);
    return (mid_element1 + mid_element2) / 2.0;
  } else {
    return ordered_elems[mid_idx];
  }
}

template <typename T>
double Mean(const std::vector<T>& elems) {
  CHECK(!elems.empty());
  double sum = 0;
  for (const auto el : elems) {
    sum += static_cast<double>(el);
  }
  return sum / elems.size();
}

template <typename T>
T Percentile(const std::vector<T>& elems, const double p) {
  CHECK(!elems.empty());
  CHECK_GE(p, 0);
  CHECK_LE(p, 100);

  const int idx = static_cast<int>(std::round(p / 100 * (elems.size() - 1)));
  const size_t percentile_idx =
      std::max(0, std::min(static_cast<int>(elems.size() - 1), idx));

  std::vector<T> ordered_elems = elems;
  std::nth_element(ordered_elems.begin(),
                   ordered_elems.begin() + percentile_idx, ordered_elems.end());

  return ordered_elems.at(percentile_idx);
}

}  // namespace math

}  // namespace limap
