#include "linemap/mapping/line_triangulation.h"

#include "linemap/util/math.h"

#include "base/line_dists.h"
#include "merging/aggregator.h"
#include "triangulation/functions.h"

namespace limap {

void ComputePoint3dUncertainty(const PointTrack& point_track,
                               const ImageCollection& img_cols,
                               double* uncertainty, double var2d,
                               const std::string& u_method) {
  CHECK_NOTNULL(uncertainty);

  const V3D& xyz = point_track.p;
  const std::vector<int>& image_ids = point_track.image_id_list;
  const size_t num_images = image_ids.size();

  THROW_CHECK_GE(num_images, 1)

  std::vector<double> all_u;
  all_u.reserve(num_images);

  for (const auto& image_id : image_ids) {
    CameraView view = img_cols.camview(image_id);
    double depth = view.pose.projdepth(xyz);
    if (depth <= 0) {
      PRINT("warning: the depth of the SfM point <= 0");
    } else {
      all_u.push_back(view.cam.uncertainty(depth, var2d));
    }
  }

  if (u_method == "median") {
    *uncertainty = math::Median<double>(all_u);
  } else if (u_method == "average") {
    *uncertainty = math::Mean<double>(all_u);
  } else if (u_method == "min") {
    *uncertainty = *std::min_element(all_u.begin(), all_u.end());
  } else {
    throw std::invalid_argument("Error: Not Implemented.");
  }
}

double ComputeTwo3DVectorAngle(const V3D& v1, const V3D& v2) {
  return acos(std::min(std::abs(v1.normalized().dot(v2.normalized())), 1.0)) *
         180.0 / M_PI;
}

double ComputeTwo2DVectorAngle(const V2D& v1, const V2D& v2) {
  return acos(std::min(std::abs(v1.normalized().dot(v2.normalized())), 1.0)) *
         180.0 / M_PI;
}

std::pair<InfiniteLine3d, bool> TriangulateInfLine3dFromTwoPlanes(
    const V4D& plane1, const V4D& plane2) {
  // compute a 3D infinite line
  M4D L_dual = plane1 * plane2.transpose() - plane2 * plane1.transpose();
  V3D d = {-L_dual(1, 2), L_dual(0, 2), -L_dual(0, 1)};
  double d_norm = d.norm();
  V3D m = {L_dual(0, 3), -L_dual(3, 1), L_dual(2, 3)};
  if (d_norm == 0.0) {
    PRINT("warning: d_norm == 0.0");
    return std::make_pair(InfiniteLine3d(), false);
  }
  d = d / d_norm;
  m = m / d_norm;
  InfiniteLine3d inf_line3d(d, m, false);

  return std::make_pair(inf_line3d, true);
}

void FindSharedPoints3DId(const line2d_t& l1_idx, const line2d_t& l2_idx,
                          const structures::PL_Bipartite2d& bqt1,
                          const structures::PL_Bipartite2d& bqt2,
                          std::unordered_set<int>* point_ids) {
  CHECK_NOTNULL(point_ids);
  point_ids->clear();

  // Ids of 3D points of the 2D points on ref line.
  std::set<point3d_t> set1;

  for (const point2d_t& point_idx : bqt1.neighbor_points(l1_idx)) {
    // Point2d on ref line2d.
    auto p = bqt1.point(point_idx);
    set1.insert(p.point3D_id);
  }

  for (const point2d_t& point_idx : bqt2.neighbor_points(l2_idx)) {
    // Point2d on neighbor line2d.
    auto p = bqt2.point(point_idx);
    if (set1.find(p.point3D_id) != set1.end()) {
      // Save shared 3D points with respect to ref and neighbor line.
      point_ids->insert(p.point3D_id);
    }
  }
}

void FindPoints3D(const line2d_t& line2d_idx,
                  const structures::PL_Bipartite2d& bqt,
                  std::unordered_set<int>* point3D_ids) {
  CHECK_NOTNULL(point3D_ids);
  point3D_ids->clear();

  for (const point2d_t& point_idx : bqt.neighbor_points(line2d_idx)) {
    const Point2d& p = bqt.point(point_idx);
    point3D_ids->insert(p.point3D_id);
  }
}

std::pair<InfiniteLine3d, bool> FitInfiniteLine3dFromLines3d(
    const std::vector<Line3d>& lines, const std::vector<double>& weights,
    const std::string& method) {
  if (lines.empty()) return {InfiniteLine3d(), false};

  if (lines.size() == 1) {
    PRINT("warning: the input lines3d.size() == 1");
    if (std::isnan(lines[0].start[0]) || std::isnan(lines[0].end[0])) {
      PRINT("warning: the input line3d has nan value");
      return {InfiniteLine3d(), false};
    } else {
      return {InfiniteLine3d(lines[0]), true};
    }
  }

  const size_t num_lines = lines.size();

  V3D center(0.0, 0.0, 0.0);
  V3D direc(0.0, 0.0, 0.0);

  if (method == "point_PCA") {
    // compute center
    for (size_t i = 0; i < num_lines; ++i) {
      center += lines[i].start;
      center += lines[i].end;
    }
    center = center / (2 * num_lines);
    // compute direction
    Eigen::MatrixXd endpoints;
    endpoints.resize(num_lines * 2, 3);
    for (size_t i = 0; i < num_lines; ++i) {
      endpoints.row(2 * i) = lines[i].start - center;
      endpoints.row(2 * i + 1) = lines[i].end - center;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(endpoints, Eigen::ComputeThinV);
    direc = svd.matrixV().col(0).normalized();
  } else if (method == "point_weighted_PCA") {
    THROW_CHECK_EQ(lines.size(), weights.size());
    // compute center
    for (size_t i = 0; i < num_lines; ++i) {
      center += weights[i] * lines[i].start;
      center += weights[i] * lines[i].end;
    }
    center = center /
             (2.0 * std::accumulate(weights.begin(), weights.end(), 0.0) + EPS);
    // compute direction
    Eigen::MatrixXd endpoints;
    endpoints.resize(num_lines * 2, 3);
    for (size_t i = 0; i < num_lines; ++i) {
      endpoints.row(2 * i) = weights[i] * (lines[i].start - center);
      endpoints.row(2 * i + 1) = weights[i] * (lines[i].end - center);
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(endpoints, Eigen::ComputeThinV);
    direc = svd.matrixV().col(0).normalized();
  } else if (method == "direction_PCA") {
    // compute center
    for (size_t i = 0; i < num_lines; ++i) {
      center += lines[i].start;
      center += lines[i].end;
    }
    center = center / (2 * num_lines);
    // compute direction
    Eigen::MatrixXd directions;
    directions.resize(num_lines, 3);
    for (size_t i = 0; i < num_lines; ++i) {
      directions.row(i) = lines[i].direction().normalized();
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(directions, Eigen::ComputeThinV);
    direc = svd.matrixV().col(0).normalized();
  } else if (method == "direction_weighted_PCA") {
    THROW_CHECK_EQ(lines.size(), weights.size());
    // compute center
    for (size_t i = 0; i < num_lines; ++i) {
      center += weights[i] * lines[i].start;
      center += weights[i] * lines[i].end;
    }
    center = center /
             (2.0 * std::accumulate(weights.begin(), weights.end(), 0.0) + EPS);
    // compute direction
    Eigen::MatrixXd directions;
    directions.resize(num_lines, 3);
    for (size_t i = 0; i < num_lines; ++i) {
      directions.row(i) = weights[i] * lines[i].direction().normalized();
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(directions, Eigen::ComputeThinV);
    direc = svd.matrixV().col(0).normalized();
  } else {
    throw std::invalid_argument("Error: Not Implemented.");
  }

  if (std::isnan(center[0]) || std::isnan(direc[0])) {
    PRINT("warning: the infinite 3D line has nan values");
    return {InfiniteLine3d(), false};
  }

  InfiniteLine3d inf_line3d(center, direc);
  return {inf_line3d, true};
}

// take the union of line3d1 and line3d2 as the final line3d
std::pair<Line3d, bool> GetLine3dWithInfLine3d(
    const InfiniteLine3d& line, const CameraView& view1,
    const CameraView& view2, const Line2d& line2d1, const Line2d& line2d2,
    const double th_angle, const double th_perp,
    const double cheirality_min_depth, const double var2d) {
  const V3D dir = line.direction();
  const V3D p_ref = line.point();

  const auto C1 = view1.pose.center();
  const auto C2 = view2.pose.center();

  InfiniteLine2d inf_line2d_proj1 = line.projection(view1);

  // check reprojection
  double angle_1 = ComputeTwo2DVectorAngle(inf_line2d_proj1.direction(),
                                           line2d1.direction());
  if (angle_1 > th_angle) return std::make_pair(Line3d(), false);

  V2D pstart_2d_1 = inf_line2d_proj1.point_projection(line2d1.start);
  V3D ray_start_1 = view1.ray_direction(pstart_2d_1);
  InfiniteLine3d line_start_1(C1, ray_start_1);
  V3D pstart_3d_1 = line.project_from_infinite_line(line_start_1);
  double tstart_1 = (pstart_3d_1 - p_ref).dot(dir);

  // check reprojection
  double dist_start_1 = (pstart_2d_1 - line2d1.start).norm();
  if (dist_start_1 > th_perp) return std::make_pair(Line3d(), false);

  V2D pend_2d_1 = inf_line2d_proj1.point_projection(line2d1.end);
  V3D ray_end_1 = view1.ray_direction(pend_2d_1);
  InfiniteLine3d line_end_1(C1, ray_end_1);
  V3D pend_3d_1 = line.project_from_infinite_line(line_end_1);
  double tend_1 = (pend_3d_1 - p_ref).dot(dir);

  // check reprojection
  double dist_end_1 = (pend_2d_1 - line2d1.end).norm();
  if (dist_end_1 > th_perp) return std::make_pair(Line3d(), false);

  InfiniteLine2d inf_line2d_proj2 = line.projection(view2);

  // check reprojection
  double angle_2 = ComputeTwo2DVectorAngle(inf_line2d_proj2.direction(),
                                           line2d2.direction());
  if (angle_2 > th_angle) return std::make_pair(Line3d(), false);

  V2D pstart_2d_2 = inf_line2d_proj2.point_projection(line2d2.start);
  V3D ray_start_2 = view2.ray_direction(pstart_2d_2);
  InfiniteLine3d line_start_2(C2, ray_start_2);
  V3D pstart_3d_2 = line.project_from_infinite_line(line_start_2);
  double tstart_2 = (pstart_3d_2 - p_ref).dot(dir);

  // check reprojection
  double dist_start_2 = (pstart_2d_2 - line2d2.start).norm();
  if (dist_start_2 > th_perp) return std::make_pair(Line3d(), false);

  V2D pend_2d_2 = inf_line2d_proj2.point_projection(line2d2.end);
  V3D ray_end_2 = view2.ray_direction(pend_2d_2);
  InfiniteLine3d line_end_2(C2, ray_end_2);
  V3D pend_3d_2 = line.project_from_infinite_line(line_end_2);
  double tend_2 = (pend_3d_2 - p_ref).dot(dir);

  // check reprojection
  double dist_end_2 = (pend_2d_2 - line2d2.end).norm();
  if (dist_end_2 > th_perp) return std::make_pair(Line3d(), false);

  // get line3d1 and line3d2
  V3D line3d_start1 = p_ref + dir * tstart_1;
  V3D line3d_end1 = p_ref + dir * tend_1;
  V3D line3d_start2 = p_ref + dir * tstart_2;
  V3D line3d_end2 = p_ref + dir * tend_2;

  // check nan
  if (std::isnan(line3d_start1[0]) || std::isnan(line3d_end1[0]) ||
      std::isnan(line3d_start2[0]) || std::isnan(line3d_end2[0])) {
    return std::make_pair(Line3d(), false);
  }

  // check ill-posed condition
  auto ray_start1 = InfiniteLine3d(C1, view1.ray_direction(line2d1.start));
  auto ray_end1 = InfiniteLine3d(C1, view1.ray_direction(line2d1.end));
  auto ray_start2 = InfiniteLine3d(C2, view2.ray_direction(line2d2.start));
  auto ray_end2 = InfiniteLine3d(C2, view2.ray_direction(line2d2.end));
  double angle_start1 =
      acos(std::min(std::abs(ray_start1.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_start1 < 1.0) return std::make_pair(Line3d(), false);
  double angle_end1 =
      acos(std::min(std::abs(ray_end1.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_end1 < 1.0) return std::make_pair(Line3d(), false);
  double angle_start2 =
      acos(std::min(std::abs(ray_start2.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_start2 < 1.0) return std::make_pair(Line3d(), false);
  double angle_end2 =
      acos(std::min(std::abs(ray_end2.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_end2 < 1.0) return std::make_pair(Line3d(), false);

  // check 3D overlap
  std::vector<std::pair<double, double>> segments;
  if (tstart_1 < tend_1) {
    segments.emplace_back(tstart_1, tend_1);
  } else {
    segments.emplace_back(tend_1, tstart_1);
  }
  if (tstart_2 < tend_2) {
    segments.emplace_back(tstart_2, tend_2);
  } else {
    segments.emplace_back(tend_2, tstart_2);
  }
  std::sort(segments.begin(), segments.end(),
            [](const std::pair<double, double>& pair1,
               const std::pair<double, double>& pair2) {
              return pair1.first < pair2.first;
            });
  if (segments[1].first >= segments[0].second)
    return std::make_pair(Line3d(), false);

  // cheirality
  double z_start1 = view1.pose.projdepth(line3d_start1);
  double z_end1 = view1.pose.projdepth(line3d_end1);
  double z_start2 = view2.pose.projdepth(line3d_start2);
  double z_end2 = view2.pose.projdepth(line3d_end2);
  if (z_start1 < cheirality_min_depth || z_end1 < cheirality_min_depth ||
      z_start2 < cheirality_min_depth || z_end2 < cheirality_min_depth) {
    return std::make_pair(Line3d(), false);
  }

  // compute the union of line3d1 and line3d2
  std::vector<double> projections = {tstart_1, tend_1, tstart_2, tend_2};
  std::sort(projections.begin(), projections.end());
  Line3d line3d;
  line3d.start = p_ref + dir * projections[0];
  line3d.end = p_ref + dir * projections[3];

  // compute line uncertainty
  double z_median1 = (z_start1 + z_end1) / 2.0;
  double z_median2 = (z_start2 + z_end2) / 2.0;
  double u1 = view1.cam.uncertainty(z_median1, var2d);
  double u2 = view2.cam.uncertainty(z_median2, var2d);
  THROW_CHECK_GT(u1, 0.0);
  THROW_CHECK_GT(u2, 0.0);
  line3d.uncertainty = std::min(u1, u2);

  return std::make_pair(line3d, true);
}

double ComputeDistance(const V3D& line2d, const V3D& point2d) {
  V3D line = line2d / (line2d(2) + EPS);
  V3D point = point2d / (point2d(2) + EPS);

  double num = std::abs(line.dot(point));
  double den = std::sqrt(line(0) * line(0) + line(1) * line(1));

  return num / (den + EPS);
}

}  // namespace limap