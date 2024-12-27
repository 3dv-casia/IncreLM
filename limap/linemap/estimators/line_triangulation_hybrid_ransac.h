
// Copyright (c) 2019, Torsten Sattler
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of the copyright holder nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// author: Torsten Sattler, torsten.sattler.de@googlemail.com

#pragma once

#include "linemap/util/types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include <RansacLib/hybrid_ransac.h>

namespace limap {

namespace estimators {

namespace line_triangulation {

using namespace ransac_lib;

// Our customized RANSAC based on HybridLocallyOptimizedMSAC from RansacLib
// [LINK]
// https://github.com/tsattler/RansacLib/blob/master/RansacLib/hybrid_ransac.h
template <class Model, class ModelVector, class HybridSolver,
          class Sampler = HybridUniformSampling<HybridSolver>>
class LineTriangulationHybridRansac : public HybridRansacBase {
 public:
  // Estimates a model using a given solver. Notice that the solver contains
  // all data.
  // Returns the number of inliers.
  int EstimateModel(const ExtendedHybridLORansacOptions& options,
                    const HybridSolver& solver, Model* best_model,
                    HybridRansacStatistics* statistics,
                    const bool init_valid) const {
    // Initializes all relevant variables.
    ResetStatistics(statistics);
    HybridRansacStatistics& stats = *statistics;

    const int kNumSolvers = solver.num_minimal_solvers();
    stats.num_iterations_per_solver.resize(kNumSolvers, 0);

    const int kNumDataTypes = solver.num_data_types();
    stats.inlier_ratios.resize(kNumDataTypes, 0.0);
    stats.inlier_indices.resize(kNumDataTypes);

    const std::vector<double>& kSqrInlierThresh =
        options.squared_inlier_thresholds_;

    std::vector<double> prior_probabilities;
    solver.solver_probabilities(&prior_probabilities);

    std::vector<std::vector<int>> min_sample_sizes;
    solver.min_sample_sizes(&min_sample_sizes);

    std::vector<int> num_data;
    solver.num_data(&num_data);

    uint32_t max_num_iterations =
        std::max(options.max_num_iterations_, options.min_num_iterations_);
    std::vector<uint32_t> max_num_iterations_per_solver(
        kNumSolvers, std::max(options.max_num_iterations_per_solver_,
                              options.min_num_iterations_));

    Model best_minimal_model;
    double best_min_model_score =
        std::numeric_limits<double>::max();  // number of outliers overall

    // score the initial best model if the initial best model is valid
    if (init_valid) {
      best_minimal_model = *best_model;
      ScoreModel(options, solver, *best_model, kSqrInlierThresh, kNumDataTypes,
                 num_data, &best_min_model_score);
      if (best_min_model_score < std::numeric_limits<double>::max()) {
        // Updates the best model.
        THROW_CHECK_EQ(stats.best_model_score,
                       std::numeric_limits<double>::max());
        // assuming that line-line solver type is -2
        UpdateBestModel(best_min_model_score, best_minimal_model, -2,
                        &(stats.best_model_score), best_model,
                        &(stats.best_solver_type));
        THROW_CHECK_LT(stats.best_model_score,
                       std::numeric_limits<double>::max());

        // Updates the number of RANSAC iterations for each solver as well
        // as the number of inliers and inlier ratios for each data type.
        UpdateRANSACTerminationCriteria(options, solver, *best_model,
                                        statistics,
                                        &max_num_iterations_per_solver);
      }
    }

    if (!VerifyData(min_sample_sizes, num_data, kNumSolvers, kNumDataTypes,
                    &prior_probabilities)) {
      if (init_valid) {  // if the initial best model is valid
        return stats.best_num_inliers;
      } else {
        return 0;
      }
    }

    Sampler sampler(options.random_seed_, solver);

    std::vector<std::vector<int>> minimal_sample(kNumDataTypes);
    ModelVector estimated_models;

    std::mt19937 rng;
    rng.seed(options.random_seed_);

    // Runs random sampling.
    for (stats.num_iterations_total = 0u;
         stats.num_iterations_total < max_num_iterations;
         ++stats.num_iterations_total) {
      const int kSolverType =
          SelectMinimalSolver(solver, prior_probabilities, stats,
                              options.min_num_iterations_, &rng);

      if (kSolverType == -1) {
        // Since no solver could be selected, we stop Hybrid RANSAC here.
        PRINT("kSolverType == -1");
        break;
      }

      stats.num_iterations_per_solver[kSolverType] += 1;

      sampler.Sample(min_sample_sizes[kSolverType], &minimal_sample);

      // MinimalSolver returns the number of estimated models.
      const int kNumEstimatedModels =
          solver.MinimalSolver(minimal_sample, kSolverType, &estimated_models);

      if (kNumEstimatedModels > 0) {
        // Finds the best model among all estimated models.
        double best_local_score = std::numeric_limits<double>::max();
        int best_local_model_id = 0;
        GetBestEstimatedModelId(options, solver, estimated_models,
                                kNumEstimatedModels, kSqrInlierThresh,
                                kNumDataTypes, num_data, &best_local_score,
                                &best_local_model_id);

        // Updates the best model found so far.
        if (best_local_score < best_min_model_score) {
          // New best model (estimated from inliers found.
          best_min_model_score = best_local_score;
          best_minimal_model = estimated_models[best_local_model_id];

          // Updates the best model.
          UpdateBestModel(best_min_model_score, best_minimal_model, kSolverType,
                          &(stats.best_model_score), best_model,
                          &(stats.best_solver_type));

          // Updates the number of RANSAC iterations for each solver as well
          // as the number of inliers and inlier ratios for each data type.
          UpdateRANSACTerminationCriteria(options, solver, *best_model,
                                          statistics,
                                          &max_num_iterations_per_solver);
        } else {
        }
      }

      // Terminate if the current solver reaches its maximum number of
      // iterations.
      if (stats.num_iterations_per_solver[kSolverType] >=
          max_num_iterations_per_solver[kSolverType]) {
        break;
      }
    }
    return stats.best_num_inliers;
  }

 protected:
  // Randomly selects a minimal solver. See Eq. 1 in Camposeco et al.
  int SelectMinimalSolver(const HybridSolver& solver,
                          const std::vector<double> prior_probabilities,
                          const HybridRansacStatistics& stats,
                          const uint32_t min_num_iterations,
                          std::mt19937* rng) const {
    double sum_probabilities = 0.0;
    const int kNumSolvers = static_cast<int>(prior_probabilities.size());
    std::vector<std::vector<int>> min_sample_sizes;
    solver.min_sample_sizes(&min_sample_sizes);
    std::vector<double> probabilities(kNumSolvers, 0);

    const int kNumDataTypes = solver.num_data_types();

    // There is a special case where all inlier ratios are 0. In this case, the
    // solvers should be sampled based on the priors.
    const double kSumInlierRatios = std::accumulate(
        stats.inlier_ratios.begin(), stats.inlier_ratios.end(), 0.0);

    if (kSumInlierRatios == 0.0) {
      for (int i = 0; i < kNumSolvers; ++i) {
        probabilities[i] = prior_probabilities[i];
        sum_probabilities += probabilities[i];
      }
    } else {
      for (int i = 0; i < kNumSolvers; ++i) {
        double num_iters =
            static_cast<double>(stats.num_iterations_per_solver[i]);
        if (num_iters > 0.0) {
          num_iters -= 1.0;
        }

        double all_inlier_prob = 1.0;
        for (int j = 0; j < kNumDataTypes; ++j) {
          all_inlier_prob *=
              std::pow(stats.inlier_ratios[j],
                       static_cast<double>(min_sample_sizes[i][j]));
        }

        if (num_iters < static_cast<double>(min_num_iterations)) {
          probabilities[i] = all_inlier_prob * prior_probabilities[i];
        } else {
          probabilities[i] = all_inlier_prob *
                             std::pow(1.0 - all_inlier_prob, num_iters) *
                             prior_probabilities[i];
        }
        sum_probabilities += probabilities[i];
      }
    }

    std::uniform_real_distribution<double> dist(0.0, sum_probabilities);

    const double kProb = dist(*rng);
    double current_prob = 0.0;
    for (int i = 0; i < kNumSolvers; ++i) {
      if (prior_probabilities[i] == 0.0) continue;

      current_prob += probabilities[i];
      if (kProb <= current_prob) return i;
    }
    return -1;
  }

  void GetBestEstimatedModelId(
      const ExtendedHybridLORansacOptions& options, const HybridSolver& solver,
      const ModelVector& models, const int num_models,
      const std::vector<double>& squared_inlier_thresholds,
      const int num_data_types, const std::vector<int> num_data,
      double* best_score, int* best_model_id) const {
    *best_score = std::numeric_limits<double>::max();
    *best_model_id = 0;

    for (int m = 0; m < num_models; ++m) {
      double score = std::numeric_limits<double>::max();
      ScoreModel(options, solver, models[m], squared_inlier_thresholds,
                 num_data_types, num_data, &score);

      if (score < *best_score) {
        *best_score = score;
        *best_model_id = m;
      }
    }
  }

  void ScoreModel(const ExtendedHybridLORansacOptions& options,
                  const HybridSolver& solver, const Model& model,
                  const std::vector<double>& squared_inlier_thresholds,
                  const int num_data_types, const std::vector<int> num_data,
                  double* score) const {
    if (!solver.CheckModel(model)) {
      *score = std::numeric_limits<double>::max();
      return;
    }

    *score = 0.0;

    for (int t = 0; t < num_data_types; ++t) {
      for (int i = 0; i < num_data[t]; ++i) {
        double squared_error = solver.EvaluateModelOnPoint(model, t, i);
        *score += ComputeScore(squared_error, squared_inlier_thresholds[t]) *
                  options.data_type_weights_[t];
      }
    }
  }

  // // MSAC (top-hat) scoring function.
  // inline double ComputeScore(const double squared_error,
  //                            const double squared_error_threshold) const {
  //   return std::min(squared_error, squared_error_threshold);
  // }

  // Inlier count scoring function.
  inline double ComputeScore(const double squared_error,
                             const double squared_error_threshold) const {
    if (squared_error < squared_error_threshold) {
      return 0.0;
    } else {
      return 1.0;
    }
  }

  int GetInliers(const HybridSolver& solver, const Model& model,
                 const std::vector<double>& squared_inlier_thresholds,
                 std::vector<std::vector<int>>* inliers) const {
    const int kNumDataTypes = solver.num_data_types();
    std::vector<int> num_data;
    solver.num_data(&num_data);
    if (inliers == nullptr) {
      int num_inliers = 0;
      for (int t = 0; t < kNumDataTypes; ++t) {
        for (int i = 0; i < num_data[i]; ++i) {
          double squared_error = solver.EvaluateModelOnPoint(model, t, i);
          if (squared_error < squared_inlier_thresholds[t]) {
            ++num_inliers;
          }
        }
      }
      return num_inliers;
    } else {
      inliers->clear();
      inliers->resize(kNumDataTypes);
      int num_inliers = 0;
      for (int t = 0; t < kNumDataTypes; ++t) {
        (*inliers)[t].clear();
        for (int i = 0; i < num_data[t]; ++i) {
          double squared_error = solver.EvaluateModelOnPoint(model, t, i);
          if (squared_error < squared_inlier_thresholds[t]) {
            ++num_inliers;
            (*inliers)[t].push_back(i);
          }
        }
      }
      return num_inliers;
    }
  }

  void UpdateRANSACTerminationCriteria(
      const ExtendedHybridLORansacOptions& options, const HybridSolver& solver,
      const Model& model, HybridRansacStatistics* statistics,
      std::vector<uint32_t>* max_iterations) const {
    statistics->best_num_inliers =
        GetInliers(solver, model, options.squared_inlier_thresholds_,
                   &(statistics->inlier_indices));

    const int kNumDataTypes = solver.num_data_types();
    std::vector<int> num_data;
    solver.num_data(&num_data);

    for (int d = 0; d < kNumDataTypes; ++d) {
      if (num_data[d] > 0) {
        statistics->inlier_ratios[d] =
            static_cast<double>(statistics->inlier_indices[d].size()) /
            static_cast<double>(num_data[d]);
      } else {
        statistics->inlier_ratios[d] = 0.0;
      }
    }

    std::vector<std::vector<int>> min_sample_sizes;
    solver.min_sample_sizes(&min_sample_sizes);
    const int kNumSolvers = solver.num_minimal_solvers();
    for (int s = 0; s < kNumSolvers; ++s) {
      (*max_iterations)[s] = utils::NumRequiredIterations(
          statistics->inlier_ratios, 1.0 - options.success_probability_,
          min_sample_sizes[s], options.min_num_iterations_,
          options.max_num_iterations_per_solver_);
    }
  }

  inline void UpdateBestModel(const double score_curr, const Model& m_curr,
                              const int solver_type, double* score_best,
                              Model* m_best, int* best_solver_type) const {
    if (score_curr < *score_best) {
      *score_best = score_curr;
      *m_best = m_curr;
      *best_solver_type = solver_type;
    }
  }

  // Determines whether enough data is available to run any of the minimal
  // solvers. Returns false otherwise. For those solvers that are not feasible
  // because not all data is available, the prior probability is set to 0.
  bool VerifyData(const std::vector<std::vector<int>>& min_sample_sizes,
                  const std::vector<int>& num_data, const int num_solvers,
                  const int num_data_types,
                  std::vector<double>* prior_probabilities) const {
    for (int i = 0; i < num_solvers; ++i) {
      for (int j = 0; j < num_data_types; ++j) {
        if (min_sample_sizes[i][j] > num_data[j] ||
            min_sample_sizes[i][j] < 0) {
          (*prior_probabilities)[i] = 0.0;
          break;
        }
      }
    }

    // No need to run HybridRANSAC if none of the solvers can be used.
    bool any_valid_solver = false;
    for (int i = 0; i < num_solvers; ++i) {
      if ((*prior_probabilities)[i] > 0.0) {
        any_valid_solver = true;
        break;
      }
    }

    if (!any_valid_solver) {
      return false;
    }

    return true;
  }
};

}  // namespace line_triangulation

}  // namespace estimators

}  // namespace limap
