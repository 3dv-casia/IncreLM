#include "linemap/estimators/line_triangulation_estimator.h"

#include <vector>

#include "_limap/helpers.h"
#include <Eigen/Core>
#include <RansacLib/ransac.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace py::literals;

namespace limap {

void bind_line_triangulation(py::module& m) {
  using namespace estimators::line_triangulation;

  py::class_<HybridLineTriangulationEstimatorOptions>(
      m, "HybridLineTriangulationEstimatorOptions")
      .def(py::init<>())
      .def_readwrite("ransac_options",
                     &HybridLineTriangulationEstimatorOptions::ransac_options)
      .def_readwrite(
          "cheirality_min_depth",
          &HybridLineTriangulationEstimatorOptions::cheirality_min_depth)
      .def_readwrite("solver_flags",
                     &HybridLineTriangulationEstimatorOptions::solver_flags)
      .def_readwrite("random", &HybridLineTriangulationEstimatorOptions::random)
      .def_readwrite("th_angle",
                     &HybridLineTriangulationEstimatorOptions::th_angle)
      .def_readwrite("th_perp",
                     &HybridLineTriangulationEstimatorOptions::th_perp)
      .def("Print", &HybridLineTriangulationEstimatorOptions::Print);
}

void bind_linemap_estimators(py::module& m) {
  using namespace estimators;

  bind_line_triangulation(m);
}

}  // namespace limap
