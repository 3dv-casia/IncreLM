#include "linemap/estimators/line_triangulation_estimator.h"

#include "linemap/estimators/line_triangulation_hybrid_ransac.h"
#include "linemap/mapping/line_triangulation.h"

#include "base/line_dists.h"

namespace limap {

namespace estimators {

namespace line_triangulation {

std::pair<InfiniteLine3d, ransac_lib::HybridRansacStatistics>
LineTriangulation_Hybrid(const std::vector<V3D>& points,
                         const std::vector<double>& un_points,
                         const std::vector<V3D>& directions,
                         const Line2d& line2d1, const Line2d& line2d2,
                         const CameraView& view1, const CameraView& view2,
                         const HybridLineTriangulationEstimatorOptions& options,
                         const InfiniteLine3d& init_inf_line3d,
                         const bool init_valid) {
  ExtendedHybridLORansacOptions ransac_options = options.ransac_options;
  std::random_device rand_dev;
  if (options.random) ransac_options.random_seed_ = rand_dev();

  THROW_CHECK_EQ(ransac_options.data_type_weights_.size(), 2);
  THROW_CHECK_EQ(ransac_options.squared_inlier_thresholds_.size(), 2);

  LineTriangulationEstimator solver(points, un_points, directions, line2d1,
                                    line2d2, view1, view2, options.solver_flags,
                                    options.cheirality_min_depth,
                                    options.th_angle, options.th_perp);

  LineTriangulationHybridRansac<InfiniteLine3d, std::vector<InfiniteLine3d>,
                                LineTriangulationEstimator>
      hybrid_ransac;
  InfiniteLine3d best_model = init_inf_line3d;
  ransac_lib::HybridRansacStatistics ransac_stats;

  hybrid_ransac.EstimateModel(ransac_options, solver, &best_model,
                              &ransac_stats, init_valid);
  return std::make_pair(best_model, ransac_stats);
}

LineTriangulationEstimator::LineTriangulationEstimator(
    const std::vector<V3D>& points, const std::vector<double>& un_points,
    const std::vector<V3D>& directions, const Line2d& line2d1,
    const Line2d& line2d2, const CameraView& view1, const CameraView& view2,
    const std::vector<bool>& solver_flags, const double cheirality_min_depth,
    const double angle2d_th, const double dist2d_th)
    : points_(&points),
      un_points_(&un_points),
      num_points_(points.size()),
      directions_(&directions),
      num_directions_(directions.size()),
      line2d1_(&line2d1),
      line2d2_(&line2d2),
      view1_(&view1),
      view2_(&view2),
      solver_flags_(solver_flags),
      cheirality_min_depth_(cheirality_min_depth),
      th_angle_(angle2d_th),
      th_perp_(dist2d_th) {
  C1_ = view1_->pose.center();
  C2_ = view2_->pose.center();
  ray_start1_ = InfiniteLine3d(C1_, view1.ray_direction(line2d1.start));
  ray_end1_ = InfiniteLine3d(C1_, view1.ray_direction(line2d1.end));
  ray_start2_ = InfiniteLine3d(C2_, view2.ray_direction(line2d2.start));
  ray_end2_ = InfiniteLine3d(C2_, view2.ray_direction(line2d2.end));
}

int LineTriangulationEstimator::MinimalSolver(
    const std::vector<std::vector<int>>& sample, const int solver_idx,
    std::vector<InfiniteLine3d>* inf_lines3d) const {
  if (solver_idx == 0) {
    return TwoPointSolver(sample, inf_lines3d);
  } else if (solver_idx == 1) {
    return PointDirectionSolver(sample, inf_lines3d);
  }
}

bool LineTriangulationEstimator::CheckModel(const InfiniteLine3d& line) const {
  const V3D dir = line.direction();
  const V3D p_ref = line.point();

  InfiniteLine2d inf_line2d_proj1 = line.projection(*view1_);

  // check reprojection
  double angle_1 = ComputeTwo2DVectorAngle(inf_line2d_proj1.direction(),
                                           line2d1_->direction());
  if (angle_1 > th_angle_) return false;

  V2D pstart_2d_1 = inf_line2d_proj1.point_projection(line2d1_->start);
  V3D ray_start_1 = view1_->ray_direction(pstart_2d_1);
  InfiniteLine3d line_start_1(C1_, ray_start_1);
  V3D pstart_3d_1 = line.project_from_infinite_line(line_start_1);
  double tstart_1 = (pstart_3d_1 - p_ref).dot(dir);

  // check reprojection
  double dist_start_1 = (pstart_2d_1 - line2d1_->start).norm();
  if (dist_start_1 > th_perp_) return false;

  V2D pend_2d_1 = inf_line2d_proj1.point_projection(line2d1_->end);
  V3D ray_end_1 = view1_->ray_direction(pend_2d_1);
  InfiniteLine3d line_end_1(C1_, ray_end_1);
  V3D pend_3d_1 = line.project_from_infinite_line(line_end_1);
  double tend_1 = (pend_3d_1 - p_ref).dot(dir);

  // check reprojection
  double dist_end_1 = (pend_2d_1 - line2d1_->end).norm();
  if (dist_end_1 > th_perp_) return false;

  InfiniteLine2d inf_line2d_proj2 = line.projection(*view2_);

  // check reprojection
  double angle_2 = ComputeTwo2DVectorAngle(inf_line2d_proj2.direction(),
                                           line2d2_->direction());
  if (angle_2 > th_angle_) return false;

  V2D pstart_2d_2 = inf_line2d_proj2.point_projection(line2d2_->start);
  V3D ray_start_2 = view2_->ray_direction(pstart_2d_2);
  InfiniteLine3d line_start_2(C2_, ray_start_2);
  V3D pstart_3d_2 = line.project_from_infinite_line(line_start_2);
  double tstart_2 = (pstart_3d_2 - p_ref).dot(dir);

  // check reprojection
  double dist_start_2 = (pstart_2d_2 - line2d2_->start).norm();
  if (dist_start_2 > th_perp_) return false;

  V2D pend_2d_2 = inf_line2d_proj2.point_projection(line2d2_->end);
  V3D ray_end_2 = view2_->ray_direction(pend_2d_2);
  InfiniteLine3d line_end_2(C2_, ray_end_2);
  V3D pend_3d_2 = line.project_from_infinite_line(line_end_2);
  double tend_2 = (pend_3d_2 - p_ref).dot(dir);

  // check reprojection
  double dist_end_2 = (pend_2d_2 - line2d2_->end).norm();
  if (dist_end_2 > th_perp_) return false;

  // get line3d1 and line3d2
  V3D line3d_start1 = p_ref + dir * tstart_1;
  V3D line3d_end1 = p_ref + dir * tend_1;
  V3D line3d_start2 = p_ref + dir * tstart_2;
  V3D line3d_end2 = p_ref + dir * tend_2;

  // check nan
  if (std::isnan(line3d_start1[0]) || std::isnan(line3d_end1[0]) ||
      std::isnan(line3d_start2[0]) || std::isnan(line3d_end2[0])) {
    return false;
  }

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
  THROW_CHECK_LE(segments[0].first, segments[1].first);

  if (segments[1].first >= segments[0].second) return false;

  // check ill-posed condition
  double angle_start1 =
      acos(std::min(std::abs(ray_start1_.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_start1 < 1.0) return false;
  double angle_end1 =
      acos(std::min(std::abs(ray_end1_.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_end1 < 1.0) return false;
  double angle_start2 =
      acos(std::min(std::abs(ray_start2_.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_start2 < 1.0) return false;
  double angle_end2 =
      acos(std::min(std::abs(ray_end2_.d.dot(dir)), 1.0)) * 180.0 / M_PI;
  if (angle_end2 < 1.0) return false;

  // cheirality
  double z_start1 = view1_->pose.projdepth(line3d_start1);
  double z_end1 = view1_->pose.projdepth(line3d_end1);
  double z_start2 = view2_->pose.projdepth(line3d_start2);
  double z_end2 = view2_->pose.projdepth(line3d_end2);
  if (z_start1 < cheirality_min_depth_ || z_end1 < cheirality_min_depth_ ||
      z_start2 < cheirality_min_depth_ || z_end2 < cheirality_min_depth_) {
    return false;
  }

  return true;
}

int LineTriangulationEstimator::TwoPointSolver(
    const std::vector<std::vector<int>>& sample,
    std::vector<InfiniteLine3d>* inf_lines3d) const {
  THROW_CHECK_EQ(sample.size(), 2)
  THROW_CHECK_EQ(sample[0].size(), 2)
  THROW_CHECK_EQ(sample[1].size(), 0)
  CHECK_NOTNULL(inf_lines3d);
  inf_lines3d->clear();

  const auto& p1 = points_->at(sample[0][0]);
  const auto& p2 = points_->at(sample[0][1]);
  V3D d = (p2 - p1).normalized();
  if (abs(d.norm()) < EPS) return 0;

  inf_lines3d->push_back(InfiniteLine3d(p1, d));
  return 1;
}

int LineTriangulationEstimator::PointDirectionSolver(
    const std::vector<std::vector<int>>& sample,
    std::vector<InfiniteLine3d>* inf_lines3d) const {
  THROW_CHECK_EQ(sample.size(), 2)
  THROW_CHECK_EQ(sample[0].size(), 1)
  THROW_CHECK_EQ(sample[1].size(), 1)
  CHECK_NOTNULL(inf_lines3d);
  inf_lines3d->clear();

  const auto& p = points_->at(sample[0][0]);
  const auto& d = directions_->at(sample[1][0]);
  if (abs(d.norm()) < EPS) return 0;

  inf_lines3d->push_back(InfiniteLine3d(p, d));
  return 1;
}

double LineTriangulationEstimator::EvaluateModelOnPoint(
    const InfiniteLine3d& inf_line3d, int t, int i) const {
  if (t == 0) {
    const auto& p = points_->at(i);
    const double u = un_points_->at(i);
    THROW_CHECK_GT(u, 0);
    double si_dist = inf_line3d.point_distance(p) / u;
    return si_dist * si_dist;
  } else if (t == 1) {
    const auto& d = directions_->at(i);
    double angle =
        acos(std::min(std::abs(inf_line3d.direction().dot(d.normalized())),
                      1.0)) *
        180.0 / M_PI;
    return angle * angle;
  } else {
    std::runtime_error("type is error");
  }
}

}  // namespace line_triangulation

}  // namespace estimators

}  // namespace limap