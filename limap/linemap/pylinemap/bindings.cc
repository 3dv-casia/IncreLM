#include "linemap/estimators/bindings.cc"

#include "linemap/mapping/bindings.cc"

#include "_limap/helpers.h"

namespace py = pybind11;

void bind_linemap_mapping(py::module&);
void bind_linemap_estimators(py::module&);

namespace limap {

void bind_linemap(py::module& m) {
  bind_linemap_mapping(m);
  bind_linemap_estimators(m);
}

}  // namespace limap
