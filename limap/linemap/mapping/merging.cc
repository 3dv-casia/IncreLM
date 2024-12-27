#include "linemap/mapping/merging.h"

#include "linemap/mapping/line_triangulation.h"
#include "linemap/util/math.h"

#include <queue>

#include "merging/merging_utils.h"
#include "triangulation/functions.h"

namespace limap {

std::vector<LineTrack> ReBuildLineTracks(
    const LineTrack& linetrack, const std::vector<CameraView>& views,
    const std::vector<bool>& valid_lines2d, const int num_outliers,
    const int min_support_images, const double var2d,
    const bool filter_incontinuous_line3d) {
  const size_t n_lines = linetrack.image_id_list.size();
  THROW_CHECK_EQ(n_lines, valid_lines2d.size());

  if (std::isnan(linetrack.line.start[0]) ||
      std::isnan(linetrack.line.end[0])) {
    return {};
  }

  const InfiniteLine3d inf_line(linetrack.line);
  const V3D dir = inf_line.direction();
  const V3D p_ref = inf_line.point();

  // (tstart, tend, support id)
  std::vector<std::tuple<double, double, size_t>> segments;
  segments.reserve(n_lines);

  // (support id, (tstart, tend))
  std::unordered_map<size_t, std::pair<double, double>> points;
  points.reserve(n_lines);
  for (size_t i = 0; i < n_lines; ++i) {
    if (!valid_lines2d[i]) continue;

    const Line2d& line2d = linetrack.line2d_list[i];
    const auto& view = views[i];
    InfiniteLine2d inf_line2d_proj = inf_line.projection(view);
    // project the two 2D endpoints to the 2d projection of the infinite line
    V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
    V3D pstart_3d = inf_line.unprojection(pstart_2d, view);
    double tstart = (pstart_3d - p_ref).dot(dir);
    V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);
    V3D pend_3d = inf_line.unprojection(pend_2d, view);
    double tend = (pend_3d - p_ref).dot(dir);

    if (tstart < tend) {
      segments.emplace_back(tstart, tend, i);
    } else {
      segments.emplace_back(tend, tstart, i);
    }

    points.emplace(i, std::make_pair(tstart, tend));
  }

  // return {} if all lines2d are invalid
  if (segments.empty()) return {};

  std::vector<std::vector<size_t>> continues_segment_ids;
  FindContinuesSegments(segments, continues_segment_ids);
  THROW_CHECK_GT(continues_segment_ids.size(), 0);

  if (filter_incontinuous_line3d && continues_segment_ids.size() > 1) {
    return {};
  }

  std::vector<LineTrack> new_tracks;
  for (size_t track_idx = 0; track_idx < continues_segment_ids.size();
       track_idx++) {
    LineTrack track;
    std::vector<double> values;
    values.reserve(2 * continues_segment_ids[track_idx].size());
    std::vector<double> uncertainties;
    uncertainties.reserve(continues_segment_ids[track_idx].size());

    for (const auto& i : continues_segment_ids[track_idx]) {
      track.image_id_list.push_back(linetrack.image_id_list[i]);
      track.line_id_list.push_back(linetrack.line_id_list[i]);
      track.line2d_list.push_back(linetrack.line2d_list[i]);
      double tstart = points.at(i).first;
      double tend = points.at(i).second;

      values.push_back(tstart);
      values.push_back(tend);

      V3D median_endpoint = p_ref + dir * ((tstart + tend) / 2.0);
      double depth = views[i].pose.projdepth(median_endpoint);
      double u = views[i].cam.uncertainty(depth, var2d);
      THROW_CHECK_GT(u, 0);
      uncertainties.push_back(u);
    }

    if (values.empty()) continue;
    if (track.count_images() < min_support_images) continue;
    if (track.count_lines() <= num_outliers) {
      PRINT("warning: num_outliers is too large");
      continue;
    }

    std::sort(values.begin(), values.end());
    Line3d final_line;
    final_line.start = p_ref + dir * values[num_outliers];
    final_line.end = p_ref + dir * values[values.size() - 1 - num_outliers];
    // Important: compute uncertainty
    final_line.uncertainty = math::Median(uncertainties);
    THROW_CHECK_GT(final_line.uncertainty, 0);
    THROW_CHECK(!std::isnan(final_line.start[0]));
    THROW_CHECK(!std::isnan(final_line.end[0]));

    track.line = final_line;

    new_tracks.push_back(track);
  }

  return new_tracks;
}

void FindContinuesSegments(
    std::vector<std::tuple<double, double, size_t>>& segments,
    std::vector<std::vector<size_t>>& continues_segment_ids) {
  continues_segment_ids.clear();
  if (segments.empty()) return;

  std::sort(segments.begin(), segments.end(),
            [](const std::tuple<double, double, size_t>& tup1,
               const std::tuple<double, double, size_t>& tup2) {
              return std::get<0>(tup1) < std::get<0>(tup2);
            });

  // check start <= end
  for (const auto& tup : segments) {
    THROW_CHECK_LE(std::get<0>(tup), std::get<1>(tup));
  }

  double end = std::get<1>(segments[0]);
  continues_segment_ids.push_back(std::vector<size_t>());
  continues_segment_ids.back().push_back(std::get<2>(segments[0]));
  for (size_t i = 1; i < segments.size(); i++) {
    if (std::get<0>(segments[i]) < end) {
      // (might) extend end
      end = std::max(end, std::get<1>(segments[i]));
      // save segment id
      continues_segment_ids.back().push_back(std::get<2>(segments[i]));
    } else {
      // reset end
      end = std::get<1>(segments[i]);
      // save segment id
      continues_segment_ids.push_back(std::vector<size_t>());
      continues_segment_ids.back().push_back(std::get<2>(segments[i]));
    }
  }

  bool debug_mode = false;
  if (debug_mode) {
    if (continues_segment_ids.size() == 1) {
      THROW_CHECK_EQ(continues_segment_ids[0].size(), segments.size());
      for (size_t i = 0; i < segments.size(); i++) {
        THROW_CHECK_EQ(std::get<2>(segments[i]),
                       continues_segment_ids[0].at(i));
      }
    }
  }
}

std::pair<Line3d, bool> MergeLines3d(
    const std::vector<Line3d>& lines3d, const std::vector<double>& weights,
    const std::unordered_set<image_line2d_t>& image_lines2d,
    const std::unordered_map<image_t, CameraView>& camviews,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const double th_angle, const double th_perp, const double var2d,
    const int num_outliers, const std::string& method) {
  if (image_lines2d.size() <= num_outliers) {
    PRINT("warning: num_outliers is too large");
    return {Line3d(), false};
  }
  THROW_CHECK_GT(image_lines2d.size(), 0);

  const double squared_th_perp = th_perp * th_perp;

  //////////////////////////////////////////////////////////////////////////////
  // [1] compute the merged infinite 3D line
  //////////////////////////////////////////////////////////////////////////////

  // compute infinite 3D line using PCA

  auto res = FitInfiniteLine3dFromLines3d(lines3d, weights, method);
  if (!res.second) {
    PRINT("warning: failed to fit infinite 3D line");
    return {Line3d(), false};
  }
  const InfiniteLine3d& inf_line3d = res.first;

  //////////////////////////////////////////////////////////////////////////////
  // [2] check reprojection, positive depth and continuity for each 2D line
  //////////////////////////////////////////////////////////////////////////////

  // only accept merge if all track elements are inliers

  const V3D dir = inf_line3d.direction();
  const V3D p_ref = inf_line3d.point();

  std::vector<std::pair<double, double>> segments;  // (tstart, tend)
  segments.reserve(image_lines2d.size());
  std::vector<double> uncertainties;
  uncertainties.reserve(image_lines2d.size());

  for (const auto& image_line2d : image_lines2d) {
    const auto& view = camviews.at(image_line2d.first);
    const Line2d& line2d =
        all_lines2d.at(image_line2d.first)[image_line2d.second];

    InfiniteLine2d inf_line2d_proj = inf_line3d.projection(view);

    // project the two 2D endpoints to the 2d projection of the infinite line
    V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
    V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);

    // check reprojection
    double angle = ComputeTwo2DVectorAngle(inf_line2d_proj.direction(),
                                           line2d.direction());
    if (angle > th_angle) return {Line3d(), false};
    double squared_dist_start = (pstart_2d - line2d.start).squaredNorm();
    if (squared_dist_start > squared_th_perp) return {Line3d(), false};
    double squared_dist_end = (pend_2d - line2d.end).squaredNorm();
    if (squared_dist_end > squared_th_perp) return {Line3d(), false};

    // compute 3D endpoints of the 2D line
    V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
    double tstart = (pstart_3d - p_ref).dot(dir);
    pstart_3d = p_ref + dir * tstart;
    V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
    double tend = (pend_3d - p_ref).dot(dir);
    pend_3d = p_ref + dir * tend;

    // check positive depth
    double z_start = view.pose.projdepth(pstart_3d);
    double z_end = view.pose.projdepth(pend_3d);
    if (z_start < EPS || z_end < EPS) {
      return {Line3d(), false};
    }

    if (tstart < tend) {
      segments.emplace_back(tstart, tend);
    } else {
      segments.emplace_back(tend, tstart);
    }

    // compute uncertainty
    double z_median = (z_start + z_end) / 2.0;
    double u = view.cam.uncertainty(z_median, var2d);
    THROW_CHECK_GT(u, 0);
    uncertainties.push_back(u);
  }

  // check the back-projection of 2D lines along the 3D infinite line are
  // continuous
  std::sort(segments.begin(), segments.end(),
            [](const std::pair<double, double>& pair1,
               const std::pair<double, double>& pair2) {
              return pair1.first < pair2.first;
            });
  double end = segments[0].second;
  for (size_t i = 1; i < segments.size(); i++) {
    if (segments[i].first < end) {
      // (might) extend end
      end = std::max(end, segments[i].second);
    } else {
      PRINT(
          "warning: the back-projection of 2D lines along the 3D infinite "
          "line are not continuous");
      return {Line3d(), false};
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // [3] union the back-projection of 2D lines
  //////////////////////////////////////////////////////////////////////////////

  std::vector<double> values;
  values.reserve(2 * segments.size());

  for (const auto& pair : segments) {
    values.push_back(pair.first);
    values.push_back(pair.second);
  }
  std::sort(values.begin(), values.end());
  Line3d merged_line3d;
  merged_line3d.start = p_ref + dir * values[num_outliers];
  merged_line3d.end = p_ref + dir * values[values.size() - 1 - num_outliers];
  // Important: compute the uncertainty
  merged_line3d.uncertainty = math::Median(uncertainties);
  THROW_CHECK_GT(merged_line3d.uncertainty, 0);

  return {merged_line3d, true};
}

Line3d ExtendLine3d(
    const Line3d& line3d,
    const std::unordered_set<image_line2d_t>& image_lines2d,
    const std::unordered_map<image_t, CameraView>& camviews,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const double var2d, const int num_outliers) {
  THROW_CHECK_GT(image_lines2d.size(), num_outliers);

  InfiniteLine3d inf_line3d(line3d);
  const V3D dir = inf_line3d.direction();
  const V3D p_ref = inf_line3d.point();

  std::vector<double> uncertainties;
  uncertainties.reserve(image_lines2d.size());
  std::vector<double> values;
  values.reserve(2 * image_lines2d.size());

  for (const auto& image_line2d : image_lines2d) {
    const auto& view = camviews.at(image_line2d.first);
    const Line2d& line2d =
        all_lines2d.at(image_line2d.first)[image_line2d.second];

    InfiniteLine2d inf_line2d_proj = inf_line3d.projection(view);

    // project the two 2D endpoints to the 2d projection of the infinite line
    V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
    V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);

    // compute 3D endpoints of the 2D line
    V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
    double tstart = (pstart_3d - p_ref).dot(dir);
    pstart_3d = p_ref + dir * tstart;
    V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
    double tend = (pend_3d - p_ref).dot(dir);
    pend_3d = p_ref + dir * tend;

    values.push_back(tstart);
    values.push_back(tend);

    double z_start = view.pose.projdepth(pstart_3d);
    double z_end = view.pose.projdepth(pend_3d);
    THROW_CHECK_GT(z_start, EPS);
    THROW_CHECK_GT(z_end, EPS);

    // compute uncertainty
    double z_median = (z_start + z_end) / 2.0;
    double u = view.cam.uncertainty(z_median, var2d);
    THROW_CHECK_GT(u, 0);
    uncertainties.push_back(u);
  }

  std::sort(values.begin(), values.end());
  Line3d extended_line3d;
  extended_line3d.start = p_ref + dir * values[num_outliers];
  extended_line3d.end = p_ref + dir * values[values.size() - 1 - num_outliers];
  // Important: compute the uncertainty
  extended_line3d.uncertainty = math::Median(uncertainties);
  THROW_CHECK_GT(extended_line3d.uncertainty, 0);

  THROW_CHECK(!std::isnan(extended_line3d.start[0]))
  THROW_CHECK(!std::isnan(extended_line3d.end[0]))
  return extended_line3d;
}

Line3d GetLine3dFromInfiniteLine3d(const InfiniteLine3d& inf_line3d,
                                   const std::vector<CameraView>& views,
                                   const std::vector<Line2d>& line2d_list,
                                   const double var2d, const int num_outliers) {
  const size_t n_supports = views.size();
  THROW_CHECK_EQ(views.size(), line2d_list.size());
  THROW_CHECK_GT(n_supports, num_outliers);

  const V3D dir = inf_line3d.direction();
  const V3D p_ref = inf_line3d.point();

  std::vector<double> uncertainties;
  uncertainties.reserve(n_supports);
  std::vector<double> values;
  values.reserve(2 * n_supports);

  for (size_t i = 0; i < n_supports; i++) {
    const auto& view = views[i];
    const Line2d& line2d = line2d_list[i];

    InfiniteLine2d inf_line2d_proj = inf_line3d.projection(view);

    // project the two 2D endpoints to the 2d projection of the infinite line
    V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
    V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);

    // compute 3D endpoints of the 2D line
    V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
    double tstart = (pstart_3d - p_ref).dot(dir);
    pstart_3d = p_ref + dir * tstart;
    V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
    double tend = (pend_3d - p_ref).dot(dir);
    pend_3d = p_ref + dir * tend;

    values.push_back(tstart);
    values.push_back(tend);

    double z_start = view.pose.projdepth(pstart_3d);
    double z_end = view.pose.projdepth(pend_3d);

    // compute uncertainty
    double z_median = (z_start + z_end) / 2.0;
    double u = view.cam.uncertainty(z_median, var2d);
    // THROW_CHECK_GT(u, 0);
    uncertainties.push_back(u);
  }

  std::sort(values.begin(), values.end());
  Line3d extended_line3d;
  extended_line3d.start = p_ref + dir * values[num_outliers];
  extended_line3d.end = p_ref + dir * values[values.size() - 1 - num_outliers];
  // Important: compute the uncertainty
  extended_line3d.uncertainty = math::Median(uncertainties);
  return extended_line3d;
}

bool TestDifferentInfLines3D(const Line2d& line2d1, const Line2d& line2d2,
                             const double overlap_th, const double angle_th,
                             const double perp_dist_th) {
  double angle = compute_angle<Line2d>(line2d1, line2d2);
  if (angle > angle_th) return true;

  double overlap = compute_bioverlap<Line2d>(line2d1, line2d2);
  if (overlap > overlap_th) return true;

  // compute max-max perpendicular distance
  double perp_dist =
      std::max(compute_distance<Line2d>(line2d1, line2d2,
                                        LineDistType::PERPENDICULAR_ONEWAY),
               compute_distance<Line2d>(line2d2, line2d1,
                                        LineDistType::PERPENDICULAR_ONEWAY));

  if (perp_dist > perp_dist_th) return true;

  return false;
}

}  // namespace limap
