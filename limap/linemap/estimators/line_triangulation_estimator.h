#pragma once

#include "base/camera_view.h"
#include "base/infinite_line.h"
#include "base/line_linker.h"
#include "base/linebase.h"
#include "estimators/extended_hybrid_ransac.h"
#include <RansacLib/ransac.h>

namespace limap {

namespace estimators {

namespace line_triangulation {

class HybridLineTriangulationEstimatorOptions {
 public:
  HybridLineTriangulationEstimatorOptions()
      : ransac_options(ExtendedHybridLORansacOptions()),
        cheirality_min_depth(0.0),
        solver_flags({true, true}),
        random(true),
        th_angle(5.0),
        th_perp(2.0) {}

  ExtendedHybridLORansacOptions ransac_options;
  double cheirality_min_depth = 0.0;
  // the first solver is two-points, the second solver is point-direction
  std::vector<bool> solver_flags = {true, true};
  bool random = true;
  // filter by reprojection
  double th_angle = 5.0;  // in degress
  double th_perp = 2.0;   // in pixles

  void Print() const {
    PRINT_VAR(ransac_options.success_probability_);
    PRINT_VAR(ransac_options.data_type_weights_.size());
    PRINT_VAR(ransac_options.data_type_weights_[0]);
    PRINT_VAR(ransac_options.data_type_weights_[1]);
    PRINT_VAR(ransac_options.final_least_squares_);
    PRINT_VAR(ransac_options.lo_starting_iterations_);
    PRINT_VAR(ransac_options.max_num_iterations_);
    PRINT_VAR(ransac_options.max_num_iterations_per_solver_);
    PRINT_VAR(ransac_options.min_num_iterations_);
    PRINT_VAR(ransac_options.min_sample_multiplicator_);
    PRINT_VAR(ransac_options.non_min_sample_multiplier_);
    PRINT_VAR(ransac_options.num_lo_steps_);
    PRINT_VAR(ransac_options.num_lsq_iterations_);
    PRINT_VAR(ransac_options.random_seed_);
    PRINT_VAR(ransac_options.squared_inlier_thresholds_.size());
    PRINT_VAR(ransac_options.squared_inlier_thresholds_[0]);
    PRINT_VAR(ransac_options.squared_inlier_thresholds_[1]);
    PRINT_VAR(cheirality_min_depth);
    PRINT_VAR(solver_flags.size());
    PRINT_VAR(solver_flags[0]);
    PRINT_VAR(solver_flags[1]);
    PRINT_VAR(random);
    PRINT_VAR(th_angle);
    PRINT_VAR(th_perp);
  }
};

class LineTriangulationEstimator {
 public:
  LineTriangulationEstimator(const std::vector<V3D>& points,
                             const std::vector<double>& un_points,
                             const std::vector<V3D>& directions,
                             const Line2d& line2d1, const Line2d& line2d2,
                             const CameraView& view1, const CameraView& view2,
                             const std::vector<bool>& solver_flags,
                             const double cheirality_min_depth,
                             const double angle2d_th, const double dist2d_th);

  inline int num_minimal_solvers() const { return 2; }

  void min_sample_sizes(std::vector<std::vector<int>>* min_sample_sizes) const {
    min_sample_sizes->resize(2);
    (*min_sample_sizes)[0] = std::vector<int>{2, 0};
    (*min_sample_sizes)[1] = std::vector<int>{1, 1};
  }

  inline int num_data_types() const { return 2; }

  inline void num_data(std::vector<int>* num_data) const {
    num_data->resize(2);
    (*num_data)[0] = num_points_;
    (*num_data)[1] = num_directions_;
  }

  // void solver_probabilities(std::vector<double>* solver_probabilities) const
  // {
  //   std::vector<std::vector<int>> sample_sizes;
  //   min_sample_sizes(&sample_sizes);
  //   solver_probabilities->resize(2);

  //   for (int i = 0; i < 2; i++) {
  //     if (!solver_flags_[i])
  //       solver_probabilities->at(i) = 0.0;
  //     else {
  //       solver_probabilities->at(i) =
  //           combination(num_points_, sample_sizes[i][0]) *
  //           combination(num_directions_, sample_sizes[i][1]);
  //     }
  //   }
  // }

  void solver_probabilities(std::vector<double>* solver_probabilities) const {
    std::vector<std::vector<int>> sample_sizes;
    min_sample_sizes(&sample_sizes);
    solver_probabilities->resize(2);

    for (int i = 0; i < 2; i++) {
      if (!solver_flags_[i])
        solver_probabilities->at(i) = 0.0;
      else {
        solver_probabilities->at(i) = 1.0;
      }
    }
  }

  int MinimalSolver(const std::vector<std::vector<int>>& sample,
                    const int solver_idx,
                    std::vector<InfiniteLine3d>* inf_lines3d) const;

  double EvaluateModelOnPoint(const InfiniteLine3d& inf_line3d, int t,
                              int i) const;

  bool CheckModel(const InfiniteLine3d& inf_line3d) const;

 protected:
  const std::vector<V3D>* points_;
  const std::vector<double>* un_points_;
  const int num_points_;
  const std::vector<V3D>* directions_;
  const int num_directions_;

  const Line2d* line2d1_;
  const Line2d* line2d2_;
  const CameraView* view1_;
  const CameraView* view2_;
  InfiniteLine3d ray_start1_;
  InfiniteLine3d ray_end1_;
  InfiniteLine3d ray_start2_;
  InfiniteLine3d ray_end2_;
  V3D C1_;
  V3D C2_;
  const double th_angle_;
  const double th_perp_;

  const std::vector<bool> solver_flags_;

  // cheirality and filtering options
  const double cheirality_min_depth_;

  // C_n^m
  unsigned long long combination(unsigned n, unsigned m) const {
    unsigned long long num = 1;
    unsigned long long denom = 1;
    for (int i = 0; i < m; i++) num *= n - i;
    for (int i = 0; i < m; i++) denom *= i + 1;
    return num / denom;
  }

  // tringulate a 3D infinite line using two 3D points
  int TwoPointSolver(const std::vector<std::vector<int>>& sample,
                     std::vector<InfiniteLine3d>* inf_lines3d) const;

  // tringulate a 3D infinite line using a 3D point and a 3D direction
  int PointDirectionSolver(const std::vector<std::vector<int>>& sample,
                           std::vector<InfiniteLine3d>* inf_lines3d) const;
};

std::pair<InfiniteLine3d, ransac_lib::HybridRansacStatistics>
LineTriangulation_Hybrid(const std::vector<V3D>& points,
                         const std::vector<double>& un_points,
                         const std::vector<V3D>& directions,
                         const Line2d& line2d1, const Line2d& line2d2,
                         const CameraView& view1, const CameraView& view2,
                         const HybridLineTriangulationEstimatorOptions& options,
                         const InfiniteLine3d& init_inf_line3d,
                         const bool init_valid);

}  // namespace line_triangulation

}  // namespace estimators

}  // namespace limap
