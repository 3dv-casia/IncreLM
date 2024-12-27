#pragma once

#include "linemap/util/types.h"

#include "base/image_collection.h"
#include "base/line_linker.h"
#include "base/linebase.h"
#include "base/linetrack.h"

namespace limap {

// re-build line track according to valid track elements
std::vector<LineTrack> ReBuildLineTracks(const LineTrack& linetrack,
                                         const std::vector<CameraView>& views,
                                         const std::vector<bool>& valid_lines2d,
                                         const int num_outliers,
                                         const int min_support_images,
                                         const double var2d,
                                         const bool filter_incontinuous_line3d);

void FindContinuesSegments(
    std::vector<std::tuple<double, double, size_t>>& segments,
    std::vector<std::vector<size_t>>& continues_segment_ids);

std::pair<Line3d, bool> MergeLines3d(
    const std::vector<Line3d>& lines3d, const std::vector<double>& weights,
    const std::unordered_set<image_line2d_t>& image_lines2d,
    const std::unordered_map<image_t, CameraView>& camviews,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const double th_angle, const double th_perp, const double var2d,
    const int num_outliers, const std::string& method);

Line3d ExtendLine3d(
    const Line3d& line3d,
    const std::unordered_set<image_line2d_t>& image_lines2d,
    const std::unordered_map<image_t, CameraView>& camviews,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const double var2d, const int num_outliers);

Line3d GetLine3dFromInfiniteLine3d(const InfiniteLine3d& inf_line3d,
                                   const std::vector<CameraView>& views,
                                   const std::vector<Line2d>& line2d_list,
                                   const double var2d, const int num_outliers);

bool TestDifferentInfLines3D(const Line2d& line2d1, const Line2d& line2d2,
                             const double overlap_th, const double angle_th,
                             const double perp_dist_th);

}  // namespace limap
