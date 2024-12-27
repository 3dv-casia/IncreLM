#include "linemap/util/misc.h"

#include <iostream>

#include "linemap/util/types.h"

namespace limap {

void PrintHeading1(const std::string& heading) {
  std::cout << std::endl << std::string(78, '=') << std::endl;
  std::cout << heading << std::endl;
  std::cout << std::string(78, '=') << std::endl << std::endl;
}

void PrintHeading2(const std::string& heading) {
  std::cout << std::endl << heading << std::endl;
  std::cout << std::string(std::min<int>(heading.size(), 78), '-') << std::endl;
}

}  // namespace limap
