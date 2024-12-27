#include "linemap/mapping/line_mapper.h"

#include <vector>

#include "_limap/helpers.h"
#include <Eigen/Core>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

namespace limap {

void bind_linemap_line_mapper(py::module& m) {
  py::class_<LineMapper>(m, "LineMapper")
      .def(py::init<>())
      .def(py::init<const LineMapperConfig&, const ImageCollection&,
                    const std::unordered_map<image_t, std::vector<Line2d>>&,
                    const std::pair<V3D, V3D>&>())
      .def(py::init<py::dict, const ImageCollection&,
                    const std::unordered_map<image_t, std::vector<Line2d>>&,
                    const std::pair<V3D, V3D>&>())
      .def("SetPLBipartite2d", &LineMapper::SetPLBipartite2d)
      .def("SetSfMPoints", &LineMapper::SetSfMPoints)
      .def("SetPointTracks", &LineMapper::SetPointTracks)
      .def("SetVPResults", &LineMapper::SetVPResults)
      .def("SetHybridRansacOptions", &LineMapper::SetHybridRansacOptions)
      .def("Stop", &LineMapper::Stop)
      .def("ExistLineTrack", &LineMapper::ExistLineTrack)
      .def("SetLineTracks", &LineMapper::SetLineTracks)
      .def("Initialize", &LineMapper::Initialize)
      .def("MatchLines2dByEpipolarIoU", &LineMapper::MatchLines2dByEpipolarIoU)
      .def("LoadMatches", &LineMapper::LoadMatches)
      .def("ConstructTriNodes", &LineMapper::ConstructTriNodes)
      .def("ConstructTriEdges", &LineMapper::ConstructTriEdges)
      .def("FindAndAddBestLineTrack", &LineMapper::FindAndAddBestLineTrack)
      .def("FindAndAddBestLineTrackFastFindMax",
           &LineMapper::FindAndAddBestLineTrackFastFindMax)
      .def("ExtendTrack", &LineMapper::ExtendTrack)
      .def("MergeLineTracks", &LineMapper::MergeLineTracks)
      .def("FilterLineTracksByReprojection",
           &LineMapper::FilterLineTracksByReprojection)
      .def("FilterLineTracksBySensitivity",
           &LineMapper::FilterLineTracksBySensitivity)
      .def("FilterLineTracksByOverlap", &LineMapper::FilterLineTracksByOverlap)
      .def("FilterLineTracksByNumImages",
           &LineMapper::FilterLineTracksByNumImages)
      .def("GetLineTrackList", &LineMapper::GetLineTrackList)
      .def("GetLineTrackMap", &LineMapper::GetLineTrackMap)
      .def("CountMatches", &LineMapper::CountMatches)
      .def("GetAllHypotheses", &LineMapper::GetAllHypotheses);
}

void bind_linemap_mapping(py::module& m) { bind_linemap_line_mapper(m); }

}  // namespace limap
