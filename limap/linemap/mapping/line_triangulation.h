#pragma once

#include "linemap/util/types.h"

#include <unordered_map>
#include <vector>

#include "base/image_collection.h"
#include "base/line_linker.h"
#include "base/linebase.h"
#include "structures/pl_bipartite.h"

namespace limap {

void ComputePoint3dUncertainty(const PointTrack& point_track,
                               const ImageCollection& img_cols,
                               double* uncertainty, double var2d,
                               const std::string& u_method);

void FindSharedPoints3DId(const line2d_t& l1_idx, const line2d_t& l2_idx,
                          const structures::PL_Bipartite2d& bqt1,
                          const structures::PL_Bipartite2d& bqt2,
                          std::unordered_set<int>* point_ids);

// point3D_ids: the ID of 3D sfm points whose 2D observations lie on the 2D line
// segment.
void FindPoints3D(const line2d_t& line2d_idx,
                  const structures::PL_Bipartite2d& bqt,
                  std::unordered_set<int>* point3D_ids);

double ComputeTwo3DVectorAngle(const V3D& v1, const V3D& v2);

double ComputeTwo2DVectorAngle(const V2D& v1, const V2D& v2);

// method: ["point_PCA", "point_weighted_PCA", "direction_PCA",
//          "direction_weighted_PCA"]
//
// weights can be empty if the method is "point_PCA" or "direction_PCA"
std::pair<InfiniteLine3d, bool> FitInfiniteLine3dFromLines3d(
    const std::vector<Line3d>& lines, const std::vector<double>& weights,
    const std::string& method);

std::pair<Line3d, bool> GetLine3dWithInfLine3d(
    const InfiniteLine3d& line, const CameraView& view1,
    const CameraView& view2, const Line2d& line2d1, const Line2d& line2d2,
    const double th_angle, const double th_perp,
    const double cheirality_min_depth, const double var2d);

std::pair<InfiniteLine3d, bool> TriangulateInfLine3dFromTwoPlanes(
    const V4D& plane1, const V4D& plane2);

double ComputeDistance(const V3D& line2d, const V3D& point2d);

}  // namespace limap
