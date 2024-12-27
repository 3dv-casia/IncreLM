#include "linemap/mapping/line_mapper.h"

#include "linemap/estimators/line_triangulation_estimator.h"
#include "linemap/mapping/line_triangulation.h"
#include "linemap/mapping/merging.h"
#include "linemap/util/math.h"
#include "linemap/util/misc.h"

#include <chrono>
#include <numeric>
#include <queue>
#include <unordered_set>

#include "base/line_dists.h"
#include "merging/aggregator.h"
#include "merging/merging_utils.h"
#include "progressbar.hpp"
#include "triangulation/functions.h"

namespace limap {

LineMapperConfig::LineMapperConfig(py::dict dict) {
  ASSIGN_PYDICT_ITEM(dict, debug_mode, bool)
  ASSIGN_PYDICT_ITEM(dict, add_halfpix, bool)
  ASSIGN_PYDICT_ITEM(dict, min_length_2d, double)
  ASSIGN_PYDICT_ITEM(dict, var2d, double)
  if (dict.contains("construct_3D_line_edges_linker2d_config"))
    construct_3D_line_edges_linker2d_config =
        LineLinker2dConfig(dict["construct_3D_line_edges_linker2d_config"]);
  if (dict.contains("construct_3D_line_edges_linker3d_config"))
    construct_3D_line_edges_linker3d_config =
        LineLinker3dConfig(dict["construct_3D_line_edges_linker3d_config"]);
  if (dict.contains("filter_by_reprojection_linker2d_config"))
    filter_by_reprojection_linker2d_config =
        LineLinker2dConfig(dict["filter_by_reprojection_linker2d_config"]);
  if (dict.contains("merge_tracks_linker2d_config"))
    merge_tracks_linker2d_config =
        LineLinker2dConfig(dict["merge_tracks_linker2d_config"]);
  if (dict.contains("merge_tracks_linker3d_config"))
    merge_tracks_linker3d_config =
        LineLinker3dConfig(dict["merge_tracks_linker3d_config"]);
  if (dict.contains("extend_tracks_linker2d_config"))
    extend_tracks_linker2d_config =
        LineLinker2dConfig(dict["extend_tracks_linker2d_config"]);
  if (dict.contains("extend_tracks_linker3d_config"))
    extend_tracks_linker3d_config =
        LineLinker3dConfig(dict["extend_tracks_linker3d_config"]);
  ASSIGN_PYDICT_ITEM(dict, max_num_one_to_many_matches, int)
  ASSIGN_PYDICT_ITEM(dict, use_range_filter, bool)
  ASSIGN_PYDICT_ITEM(dict, use_line_line_solver, bool)
  ASSIGN_PYDICT_ITEM(dict, line_tri_angle_threshold, double)
  ASSIGN_PYDICT_ITEM(dict, min_score, double)
  ASSIGN_PYDICT_ITEM(dict, num_outliers_aggregator_for_line_track_extension,
                     int)
  ASSIGN_PYDICT_ITEM(dict, extend_track_sensitivity_threshold, double)
  ASSIGN_PYDICT_ITEM(dict, require_tri_node, bool)
  ASSIGN_PYDICT_ITEM(dict, test_tri_node, bool)
  ASSIGN_PYDICT_ITEM(dict, num_outliers_aggregator_for_line_track_merging, int)
  ASSIGN_PYDICT_ITEM(dict, num_outliers_aggregator_for_reprojection_filtering,
                     int)
  ASSIGN_PYDICT_ITEM(dict, min_support_images, int)
  ASSIGN_PYDICT_ITEM(dict, filter_incontinuous_line3d, bool)
  ASSIGN_PYDICT_ITEM(dict, filtering_sensitivity_threshold, double)
  ASSIGN_PYDICT_ITEM(dict, th_depth, double)
  ASSIGN_PYDICT_ITEM(dict, th_sv_angular_3d, double)
  ASSIGN_PYDICT_ITEM(dict, th_sv_num_supports, int)
  ASSIGN_PYDICT_ITEM(dict, th_overlap, double)
  ASSIGN_PYDICT_ITEM(dict, th_overlap_num_supports, int)
  ASSIGN_PYDICT_ITEM(dict, use_unshared_points, bool)
  ASSIGN_PYDICT_ITEM(dict, max_point_line_reproject_distance, double)
  ASSIGN_PYDICT_ITEM(dict, use_collinearity2D, bool)
}

LineMapper::LineMapper(
    const LineMapperConfig& config, const ImageCollection& imagecols,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const std::pair<V3D, V3D>& ranges)
    : config_(config),
      imagecols_(imagecols),
      all_lines2d_(all_lines2d),
      ranges_(ranges),
      construct_3D_line_edges_linker_(
          config.construct_3D_line_edges_linker2d_config,
          config_.construct_3D_line_edges_linker3d_config),
      filter_by_reprojection_linker2d_(
          config.filter_by_reprojection_linker2d_config),
      merge_tracks_linker_(config.merge_tracks_linker2d_config,
                           config.merge_tracks_linker3d_config),
      extend_tracks_linker_(config.extend_tracks_linker2d_config,
                            config.extend_tracks_linker3d_config) {}

LineMapper::LineMapper(
    py::dict dict, const ImageCollection& imagecols,
    const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
    const std::pair<V3D, V3D>& ranges)
    : LineMapper(LineMapperConfig(dict), imagecols, all_lines2d, ranges) {}

void LineMapper::Initialize() {
  PreStartCheck();

  if (config_.add_halfpix) OffsetHalfPixel();

  // compute points uncertainty
  un_points_.clear();
  for (auto it = point_tracks_.begin(); it != point_tracks_.end(); it++) {
    point3d_t point_id = it->first;
    const auto& point_track = it->second;
    double uncertainty;
    ComputePoint3dUncertainty(point_track, imagecols_, &uncertainty, 1.0,
                              "median");
    un_points_.emplace(point_id, uncertainty);
  }
  // point tracks will not be used later
  point_tracks_.clear();

  // initialize empty containers
  size_t num_image = imagecols_.NumImages();
  const auto image_ids = imagecols_.get_img_ids();
  matches_.reserve(num_image);
  camviews_.clear();
  camviews_.reserve(num_image);
  num_lines2d_ = 0;
  for (image_t image_id : image_ids) {
    size_t num_lines2d = all_lines2d_.at(image_id).size();
    num_lines2d_ += num_lines2d;
    matches_.emplace(
        image_id,
        std::vector<
            std::unordered_map<image_t, std::unordered_map<line2d_t, tri_t>>>(
            num_lines2d));
    camviews_[image_id] = imagecols_.camview(image_id);
  }

  num_added_tracks_ = 0;
  tracks_.clear();
  stop_ = false;
  new_track_id_ = kInvalidLine3dId;

  // configure linker
  extend_tracks_linker_.linker_3d.config.set_to_spatial_merging();
  extend_tracks_linker_.linker_2d.config.set_to_default();
  merge_tracks_linker_.linker_3d.config.set_to_spatial_merging();

  // compute number of 2D LSs that has at least a associated 3D point / VP
  size_t num_point_association_line2d = 0;
  size_t num_vp_association_line2d = 0;
  for (image_t image_id : image_ids) {
    const auto& pl_bpt = all_pl_bpt2d_ptr_->at(image_id);
    size_t num_lines2d = all_lines2d_.at(image_id).size();
    const auto& vp_results = vp_results_ptr_->at(image_id);
    for (size_t line2d_idx = 0; line2d_idx < num_lines2d; line2d_idx++) {
      if (pl_bpt.ldegree(line2d_idx) > 0) num_point_association_line2d++;
      if (vp_results.HasVP(line2d_idx)) num_vp_association_line2d++;
    }
  }
  std::cout << "# 2D LSs with at least a associated 3D point / # all 2D LSs = "
            << num_point_association_line2d << " / " << num_lines2d_ << " = "
            << double(num_point_association_line2d) / double(num_lines2d_)
            << std::endl;
  std::cout << "# 2D LSs with a associated VP / # all 2D LSs = "
            << num_vp_association_line2d << " / " << num_lines2d_ << " = "
            << double(num_vp_association_line2d) / double(num_lines2d_)
            << std::endl;
}

void LineMapper::PreStartCheck() const {
  if (imagecols_.NumImages() == 0 || imagecols_.NumCameras() == 0) {
    throw std::runtime_error(
        "Error: imagecols_ is empty, please set it first.");
  }
  if (all_lines2d_.empty()) {
    throw std::runtime_error(
        "Error: all_lines2d_ is empty, please set it first.");
  }
  THROW_CHECK(imagecols_.IsUndistorted());
}

void LineMapper::OffsetHalfPixel() {
  std::vector<int> image_ids = imagecols_.get_img_ids();
  for (auto it = image_ids.begin(); it != image_ids.end(); ++it) {
    int img_id = *it;
    auto& lines = all_lines2d_.at(img_id);
    for (size_t line_id = 0; line_id < lines.size(); ++line_id) {
      auto& line = lines[line_id];
      line.start = line.start + V2D(0.5, 0.5);
      line.end = line.end + V2D(0.5, 0.5);
    }
  }
}

void LineMapper::MatchLines2dByEpipolarIoU(
    const std::unordered_map<image_t, std::vector<image_t>>& all_neighbors,
    const int n_matches, const double th_IoU) {
  std::vector<image_t> image_ids;
  image_ids.reserve(all_lines2d_.size());
  for (const auto& pair : all_lines2d_) {
    image_ids.emplace_back(pair.first);
  }

  // compute fundamental matrix
  std::unordered_map<std::pair<image_t, image_t>, M3D> Fs;
  for (int i = 0; i < image_ids.size(); i++) {
    image_t image_id = image_ids[i];
    const auto& view = camviews_.at(image_id);
    const auto& neighbors = all_neighbors.at(image_id);
    for (const auto& ng_image_id : neighbors) {
      const auto& ng_view = camviews_.at(ng_image_id);
      if (image_id == ng_image_id) continue;
      if (image_id < ng_image_id) {
        if (Fs.find(std::make_pair(image_id, ng_image_id)) != Fs.end()) {
          continue;
        }
        M3D F = triangulation::compute_fundamental_matrix(view, ng_view);
        Fs.emplace(std::make_pair(image_id, ng_image_id), F);
      } else {
        if (Fs.find(std::make_pair(ng_image_id, image_id)) != Fs.end()) {
          continue;
        }
        M3D F = triangulation::compute_fundamental_matrix(ng_view, view);
        Fs.emplace(std::make_pair(ng_image_id, image_id), F);
      }
    }
  }

  // matching
  progressbar bar(image_ids.size());
#pragma omp parallel for
  for (size_t i = 0; i < image_ids.size(); i++) {
    bar.update();
    const image_t image_id = image_ids[i];
    const auto& all_lines2d = all_lines2d_.at(image_id);
    const auto& neighbors = all_neighbors.at(image_id);

    std::vector<
        std::unordered_map<image_t, std::unordered_map<line2d_t, tri_t>>>
        all_matches(all_lines2d.size());

    for (const auto& ng_image_id : neighbors) {
      if (image_id == ng_image_id) continue;
      const auto& all_ng_lines2d = all_lines2d_.at(ng_image_id);

      M3D F;
      if (image_id < ng_image_id) {
        F = Fs.at(std::make_pair(image_id, ng_image_id));
      } else {
        F = (Fs.at(std::make_pair(ng_image_id, image_id))).transpose();
      }

      for (line2d_t line2d_idx = 0; line2d_idx < all_lines2d.size();
           line2d_idx++) {
        auto& matches = all_matches.at(line2d_idx);

        // init
        matches[ng_image_id] = std::unordered_map<line2d_t, tri_t>();
        auto& m = matches.at(ng_image_id);

        // compute epipolar line
        const auto& line2d = all_lines2d[line2d_idx];
        V3D coor_epline_start =
            (F * V3D(line2d.start[0], line2d.start[1], 1)).normalized();
        V3D coor_epline_end =
            (F * V3D(line2d.end[0], line2d.end[1], 1)).normalized();

        // initialize priority queue
        auto cmp = [](const std::pair<line2d_t, double>& a,
                      const std::pair<line2d_t, double>& b) {
          return a.second < b.second;
        };
        std::priority_queue<std::pair<line2d_t, double>,
                            std::vector<std::pair<line2d_t, double>>,
                            decltype(cmp)>
            pq(cmp);

        std::vector<line2d_t> valid_ng_lines2d;

        for (line2d_t ng_line2d_idx = 0; ng_line2d_idx < all_ng_lines2d.size();
             ng_line2d_idx++) {
          const auto& ng_line2d = all_ng_lines2d[ng_line2d_idx];
          // compute epipolar IoU
          V3D coor_ng_line = ng_line2d.coords();
          V3D homo_c_start = coor_ng_line.cross(coor_epline_start);
          V2D c_start = dehomogeneous(homo_c_start);
          V3D homo_c_end = coor_ng_line.cross(coor_epline_end);
          V2D c_end = dehomogeneous(homo_c_end);
          V2D ng_line2d_dir = ng_line2d.direction();
          double ng_line2d_length = ng_line2d.length();
          double c1 =
              (c_start - ng_line2d.start).dot(ng_line2d_dir) / ng_line2d_length;
          double c2 =
              (c_end - ng_line2d.start).dot(ng_line2d_dir) / ng_line2d_length;
          if (c1 > c2) std::swap(c1, c2);
          double IoU = (std::min(c2, 1.0) - std::max(c1, 0.0)) /
                       (std::max(c2, 1.0) - std::min(c1, 0.0));

          if (IoU >= th_IoU) {
            if (n_matches > 0) {  // use top K matching
              pq.push(std::make_pair(ng_line2d_idx, IoU));
            } else {
              valid_ng_lines2d.push_back(ng_line2d_idx);
            }
          }
        }

        // store matches
        if (n_matches > 0) {
          if (pq.empty()) {  // no match
            matches.erase(ng_image_id);
          } else {
            while (!pq.empty() && (m.size() < n_matches)) {
              m.emplace(pq.top().first, kInvalidTriIdx);
              pq.pop();
            }
          }
        } else {
          if (valid_ng_lines2d.empty()) {  // no match
            matches.erase(ng_image_id);
          } else {
            for (const auto& ng_line2d_idx : valid_ng_lines2d) {
              m.emplace(ng_line2d_idx, kInvalidTriIdx);
            }
          }
        }
      }
    }
#pragma omp critical
    { matches_.at(image_id) = all_matches; }
  }

  // make matches bidirectional
  for (const auto& pair1 : matches_) {
    image_t image_id1 = pair1.first;
    for (line2d_t line2d_idx1 = 0; line2d_idx1 < pair1.second.size();
         line2d_idx1++) {
      for (const auto& pair2 : pair1.second[line2d_idx1]) {
        image_t image_id2 = pair2.first;
        auto& matches_2 = matches_.at(image_id2);
        for (const auto& pair3 : pair2.second) {
          line2d_t line2d_idx2 = pair3.first;
          auto& m_2 = matches_2.at(line2d_idx2);
          if (m_2.find(image_id1) == m_2.end()) {
            m_2[image_id1] = std::unordered_map<line2d_t, tri_t>();
          }
          auto& m_22 = m_2.at(image_id1);
          if (m_22.find(line2d_idx1) == m_22.end()) {
            m_22.emplace(line2d_idx1, kInvalidTriIdx);
          }
        }
      }
    }
  }
}

void LineMapper::LoadMatches(
    const image_t& ref_image_id,
    const std::unordered_map<image_t, Eigen::MatrixXi>& init_matches) {
  const size_t num_ref_line2d = all_lines2d_.at(ref_image_id).size();
  if (num_ref_line2d == 0) return;
  auto& ref_matches = matches_.at(ref_image_id);

  // load matches and initialize tri_idx as kInvalidTriIdx
  for (const auto& init_match : init_matches) {
    const image_t ng_image_id = init_match.first;
    if (ref_image_id == ng_image_id) continue;
    auto& ng_matches = matches_.at(ng_image_id);

    const size_t num_ng_line2d = all_lines2d_.at(ng_image_id).size();

    const Eigen::MatrixXi& match_info = init_match.second;
    const size_t num_matches = match_info.rows();
    if (num_matches == 0) {
      continue;
    } else {
      THROW_CHECK_GT(num_matches, 0);
      THROW_CHECK_EQ(match_info.cols(), 2);
    }

    std::vector<size_t> num_one_to_many_matches(num_ref_line2d, 0);

    for (size_t i = 0; i < num_matches; i++) {
      const line2d_t ref_line2d_idx = match_info(i, 0);
      const line2d_t ng_line2d_idx = match_info(i, 1);
      if (ref_line2d_idx >= num_ref_line2d || ref_line2d_idx < 0 ||
          ng_line2d_idx >= num_ng_line2d || ng_line2d_idx < 0) {
        throw std::runtime_error(
            "Error: Out-of-index matches exist between image (img_id = " +
            std::to_string(ref_image_id) + ") and neighbor image (img_id = " +
            std::to_string(ng_image_id) + ").");
      }
      if (num_one_to_many_matches[ref_line2d_idx] >=
          config_.max_num_one_to_many_matches) {
        continue;
      }

      auto& matches = ref_matches[ref_line2d_idx][ng_image_id];
      matches.emplace(ng_line2d_idx, kInvalidTriIdx);
      ng_matches[ng_line2d_idx][ref_image_id].emplace(ref_line2d_idx,
                                                      kInvalidTriIdx);
      num_one_to_many_matches[ref_line2d_idx]++;
    }
  }
}

size_t LineMapper::CountMatches() {
  size_t num_matches = 0;
  for (const auto& pair1 : matches_) {
    for (const auto& pair2 : pair1.second) {
      for (const auto& pair3 : pair2) {
        num_matches += pair3.second.size();
      }
    }
  }
  num_matches = num_matches / 2;
  num_matches_ = num_matches;

  if (config_.debug_mode) {
    for (const auto& pair1 : matches_) {
      image_t image_id1 = pair1.first;
      for (line2d_t line2d_idx1 = 0;
           line2d_idx1 < all_lines2d_.at(image_id1).size(); line2d_idx1++) {
        const auto& pair2 = pair1.second.at(line2d_idx1);
        for (const auto& pair3 : pair2) {
          image_t image_id2 = pair3.first;
          for (const auto& pair4 : pair3.second) {
            const line2d_t line2d_idx2 = pair4.first;
            THROW_CHECK(
                matches_.at(image_id2).at(line2d_idx2).find(image_id1) !=
                matches_.at(image_id2).at(line2d_idx2).end());
            THROW_CHECK(
                matches_.at(image_id2)
                    .at(line2d_idx2)
                    .at(image_id1)
                    .find(line2d_idx1) !=
                matches_.at(image_id2).at(line2d_idx2).at(image_id1).end());
          }
        }
      }
    }
  }

  std::cout << "# line matches = " << num_matches_ << std::endl;
  return num_matches_;
}

void LineMapper::CreateCache() {
  PrintHeading1("create 2D line caches for extending tracks");

  simple_cache_to_image_line2d_.clear();
  image_line2d_to_simple_cache_.clear();
  simple_line2d_cache_list_.clear();
  simple_cache_to_image_line2d_.reserve(num_lines2d_);
  image_line2d_to_simple_cache_.reserve(num_lines2d_);
  simple_line2d_cache_list_.reserve(num_lines2d_);
  const auto image_ids = imagecols_.get_img_ids();
  size_t i = 0;
  for (size_t image_vec_idx = 0; image_vec_idx < image_ids.size();
       image_vec_idx++) {
    const image_t image_id = image_ids[image_vec_idx];
    auto& lines2d = all_lines2d_.at(image_id);
    const CameraView* view_ptr = &(camviews_.at(image_id));
    for (line2d_t line2d_idx = 0; line2d_idx < lines2d.size(); line2d_idx++) {
      SimpleLine2dCache line2d_cache;
      line2d_cache.line2d_ptr = &(lines2d[line2d_idx]);
      line2d_cache.view_ptr = view_ptr;

      simple_line2d_cache_list_.push_back(line2d_cache);
      image_line2d_to_simple_cache_[std::make_pair(image_id, line2d_idx)] = i++;
      simple_cache_to_image_line2d_.emplace_back(image_id, line2d_idx);
    }
  }
  THROW_CHECK_EQ(simple_line2d_cache_list_.size(), num_lines2d_);
  THROW_CHECK_EQ(image_line2d_to_simple_cache_.size(), num_lines2d_);
  THROW_CHECK_EQ(simple_cache_to_image_line2d_.size(), num_lines2d_);
  if (config_.debug_mode) {
    // test mapping
    for (size_t i = 0; i < simple_cache_to_image_line2d_.size(); i++) {
      image_t image_id = simple_cache_to_image_line2d_[i].first;
      line2d_t line2d_idx = simple_cache_to_image_line2d_[i].second;
      THROW_CHECK_EQ(image_line2d_to_simple_cache_.at(
                         std::make_pair(image_id, line2d_idx)),
                     i);
    }
  }

  progressbar bar(num_lines2d_);
#pragma omp parallel for
  for (size_t i = 0; i < num_lines2d_; i++) {
    bar.update();
    const image_t image_id = simple_cache_to_image_line2d_[i].first;
    const line2d_t line2d_idx = simple_cache_to_image_line2d_[i].second;
    auto& cache = simple_line2d_cache_list_[i];
    for (const auto& pair1 : matches_.at(image_id)[line2d_idx]) {
      const image_t ng_image_id = pair1.first;
      THROW_CHECK_NE(image_id, ng_image_id);
      for (const auto& pair2 : pair1.second) {
        const line2d_t ng_line2d_idx = pair2.first;
        const size_t cache_idx = image_line2d_to_simple_cache_.at(
            std::make_pair(ng_image_id, ng_line2d_idx));
        cache.ng_cache_tri_list.emplace_back(cache_idx, pair2.second);
      }
    }
    if (config_.debug_mode) {
      std::set<size_t> set1(cache.ng_lines2d_cache_idx.begin(),
                            cache.ng_lines2d_cache_idx.end());
      THROW_CHECK_EQ(set1.size(), cache.ng_lines2d_cache_idx.size());
    }
  }
}

void LineMapper::ConstructTriNodes() {
  tri_nodes_.clear();
  tri_nodes_.shrink_to_fit();

  // compute camera center
  std::unordered_map<image_t, V3D> centers;
  centers.reserve(camviews_.size());
  for (const auto& pair : camviews_) {
    centers[pair.first] = pair.second.pose.center();
  }

  // compute endpoints_inf_line3d, back-projected planes, vp directions
  std::unordered_map<image_t, std::vector<V4D>> backproj_planes;
  std::unordered_map<image_t, std::vector<V3D>> vp_directions;
  for (auto it = all_lines2d_.begin(); it != all_lines2d_.end(); it++) {
    const image_t image_id = it->first;
    const auto& lines2d = it->second;
    const auto& view = camviews_.at(image_id);
    const auto& vp_results = vp_results_ptr_->at(image_id);
    const V3D& C = centers.at(image_id);
    const size_t num_lines2d = lines2d.size();

    std::vector<V4D> planes;  // homogenous coordinates
    planes.resize(num_lines2d, V4D::Zero());

    std::vector<V3D> direcs;  // 3D direction of VP
    direcs.resize(num_lines2d, V3D::Zero());

#pragma omp parallel for
    for (size_t i = 0; i < num_lines2d; i++) {
      // compute endpoints_inf_line3d of 2D LS
      V3D ray_start = view.ray_direction(lines2d[i].start);
      V3D ray_end = view.ray_direction(lines2d[i].end);

      // compute back-projected plane of 2D LS
      V3D n = (ray_start.cross(ray_end)).normalized();
      double w = n.dot(C);
      planes[i] = {n[0], n[1], n[2], -w};

      // compute the 3D direction of VP
      if (vp_results.HasVP(i)) {
        direcs[i] =
            limap::triangulation::getDirectionFromVP(vp_results.GetVP(i), view);
      }
    }
    backproj_planes[image_id] = planes;
    vp_directions[image_id] = direcs;
  }

  // two-view line triangulation
  const auto image_ids = imagecols_.get_img_ids();
  size_t num_line_line = 0;
  size_t num_two_points = 0;
  size_t num_point_vp = 0;
  size_t num_degeneracy_line = 0;
  size_t num_use_two_points_solver_matches = 0;
  size_t num_use_point_vp_solver_matches = 0;
  size_t num_use_lineline_solver_matches = 0;
  progressbar bar(image_ids.size());
#pragma omp parallel for
  for (size_t image_vec_idx = 0; image_vec_idx < image_ids.size();
       image_vec_idx++) {
    bar.update();
    const image_t image_id1 = image_ids[image_vec_idx];
    const auto& all_matches = matches_.at(image_id1);
    const size_t num_lines1 = all_matches.size();
    const auto& all_lines2d1 = all_lines2d_.at(image_id1);
    const auto& view1 = camviews_.at(image_id1);
    const auto& planes1 = backproj_planes.at(image_id1);

    const auto& pl_bpt1 = all_pl_bpt2d_ptr_->at(image_id1);
    const auto& vp_directions1 = vp_directions.at(image_id1);

    // use buffer to speed up multithreading
    std::vector<std::tuple<Line3d, image_t, line2d_t, image_t, line2d_t>>
        tri_nodes_buffer;

    for (line2d_t line2d_idx1 = 0; line2d_idx1 < num_lines1; line2d_idx1++) {
      const Line2d& l1 = all_lines2d1[line2d_idx1];
      if (l1.length() <= config_.min_length_2d) continue;
      const V4D& plane1 = planes1.at(line2d_idx1);
      bool valid_direc1 = true;
      const V3D& direc1 = vp_directions1.at(line2d_idx1);
      if (direc1.isApprox(V3D::Zero())) valid_direc1 = false;

      std::unordered_set<int> ref_point3D_ids;
      if (config_.use_unshared_points) {
        // collect the ID of 3D sfm points whose 2D observations are on the
        // reference 2D line segment.
        FindPoints3D(line2d_idx1, pl_bpt1, &ref_point3D_ids);
      }

      for (const auto& matches : all_matches[line2d_idx1]) {
        const image_t image_id2 = matches.first;
        THROW_CHECK_NE(image_id1, image_id2);
        if (image_id1 < image_id2 && (image_id1 + image_id2) % 2 == 0) continue;
        if (image_id1 > image_id2 && (image_id1 + image_id2) % 2 == 1) continue;

        const auto& all_lines2d2 = all_lines2d_.at(image_id2);
        const auto& view2 = camviews_.at(image_id2);
        const auto& planes2 = backproj_planes.at(image_id2);

        const auto& pl_bpt2 = all_pl_bpt2d_ptr_->at(image_id2);
        const auto& vp_directions2 = vp_directions.at(image_id2);

        for (const auto& match : matches.second) {
          const line2d_t line2d_idx2 = match.first;
          const Line2d& l2 = all_lines2d2[line2d_idx2];
          if (l2.length() <= config_.min_length_2d) continue;
          bool valid_direc2 = true;
          const V3D& direc2 = vp_directions2.at(line2d_idx2);
          if (direc2.isApprox(V3D::Zero())) valid_direc2 = false;

          const V4D& plane2 = planes2.at(line2d_idx2);

          double plane_angle =
              ComputeTwo3DVectorAngle(plane1.head<3>(), plane2.head<3>());
          if (plane_angle < config_.line_tri_angle_threshold) {
#pragma omp atomic
            num_degeneracy_line += 1;
          }

          InfiniteLine3d algebraic_inf_line3d;
          bool algebraic_valid = false;

          if (config_.use_line_line_solver) {
            // triangulate
            auto res = TriangulateInfLine3dFromTwoPlanes(plane1, plane2);
            algebraic_inf_line3d = res.first;
            algebraic_valid = res.second;

#pragma omp atomic
            num_use_lineline_solver_matches++;
          }

          // collect shared 3D points
          std::unordered_set<int> shared_point_ids;
          FindSharedPoints3DId(line2d_idx1, line2d_idx2, pl_bpt1, pl_bpt2,
                               &shared_point_ids);

          std::unordered_set<int> valid_point3D_ids = shared_point_ids;
          shared_point_ids.clear();

          if (config_.use_unshared_points) {
            std::unordered_set<int> ng_point3D_ids;
            FindPoints3D(line2d_idx2, pl_bpt2, &ng_point3D_ids);

            // Filter point3D according to the distance between the reprojected
            // point2D to reference or neighboring 2D infinite line.
            for (const int& point3D_id : ng_point3D_ids) {
              if (point3D_id == -1) continue;
              if (valid_point3D_ids.find(point3D_id) != valid_point3D_ids.end())
                continue;
              const V3D& point3D =
                  sfm_points_ptr_->at(static_cast<point3d_t>(point3D_id));
              V2D point2D = view1.projection(point3D);
              V3D point2D_homo = homogeneous(point2D);
              double d = ComputeDistance(l1.coords(), point2D_homo);
              if (d < config_.max_point_line_reproject_distance) {
                valid_point3D_ids.insert(point3D_id);
              }
            }
            ng_point3D_ids.clear();
            for (const int& point3D_id : ref_point3D_ids) {
              if (point3D_id == -1) continue;
              if (valid_point3D_ids.find(point3D_id) != valid_point3D_ids.end())
                continue;
              const V3D& point3D =
                  sfm_points_ptr_->at(static_cast<point3d_t>(point3D_id));
              V2D ng_point2D = view2.projection(point3D);
              V3D ng_point2D_homo = homogeneous(ng_point2D);
              double d = ComputeDistance(l2.coords(), ng_point2D_homo);
              if (d < config_.max_point_line_reproject_distance) {
                valid_point3D_ids.insert(point3D_id);
              }
            }
          }

          std::vector<V3D> points3d;
          points3d.reserve(valid_point3D_ids.size());
          std::vector<double> un_points3d;
          un_points3d.reserve(valid_point3D_ids.size());
          for (const int& point3d_id : valid_point3D_ids) {
            points3d.push_back(
                sfm_points_ptr_->at(static_cast<point3d_t>(point3d_id)));
            un_points3d.push_back(
                un_points_.at(static_cast<point3d_t>(point3d_id)));
          }

          if (points3d.size() >= 2) {
#pragma omp atomic
            num_use_two_points_solver_matches++;
          }

          // collect vp directions
          std::vector<V3D> directions;
          if (valid_direc1) directions.push_back(direc1);
          if (valid_direc2) directions.push_back(direc2);

          if (directions.size() >= 1 && points3d.size() >= 1) {
#pragma omp atomic
            num_use_point_vp_solver_matches++;
          }

          // hybrid ransac for line triangulation
          auto inf_line3d_res =
              limap::estimators::line_triangulation::LineTriangulation_Hybrid(
                  points3d, un_points3d, directions, l1, l2, view1, view2,
                  config_.hybrid_ransac_options, algebraic_inf_line3d,
                  algebraic_valid);

          Line3d final_line3d;
          bool valid_final_line3d = false;

          if (inf_line3d_res.second.best_model_score <
              std::numeric_limits<double>::max()) {
            auto line3d_res = GetLine3dWithInfLine3d(
                inf_line3d_res.first, view1, view2, l1, l2,
                config_.hybrid_ransac_options.th_angle,
                config_.hybrid_ransac_options.th_perp,
                config_.hybrid_ransac_options.cheirality_min_depth,
                config_.var2d);
            THROW_CHECK(line3d_res.second);
            THROW_CHECK_GT(line3d_res.first.uncertainty, 0.0);

            if (config_.use_range_filter) {
              // filter the 3D line segment that are out of range
              if (!limap::triangulation::test_line_inside_ranges(
                      line3d_res.first, ranges_)) {
                continue;
              }
            }

            final_line3d = line3d_res.first;
            valid_final_line3d = true;

            if (inf_line3d_res.second.best_solver_type == -2) {
#pragma omp atomic
              num_line_line += 1;
            } else if (inf_line3d_res.second.best_solver_type == 0) {
#pragma omp atomic
              num_two_points += 1;
            } else if (inf_line3d_res.second.best_solver_type == 1) {
#pragma omp atomic
              num_point_vp += 1;
            } else {
              throw std::runtime_error("error: invalid best_solver_type");
            }
          }

          // add a tri node
          if (valid_final_line3d) {
            if (image_id1 < image_id2) {
              tri_nodes_buffer.emplace_back(
                  final_line3d, image_id1, line2d_idx1, image_id2, line2d_idx2);
            } else {
              tri_nodes_buffer.emplace_back(
                  final_line3d, image_id2, line2d_idx2, image_id1, line2d_idx1);
            }

            // // save hypotheses
            // #pragma omp critical
            //             {
            //               LineTrack track;
            //               track.line = final_line3d;
            //               track.image_id_list.push_back(image_id1);
            //               track.line_id_list.push_back(line2d_idx1);
            //               track.line2d_list.push_back(
            //                   all_lines2d_.at(image_id1).at(line2d_idx1));
            //               track.image_id_list.push_back(image_id2);
            //               track.line_id_list.push_back(line2d_idx2);
            //               track.line2d_list.push_back(
            //                   all_lines2d_.at(image_id2).at(line2d_idx2));
            //               all_hypotheses_.push_back(track);
            //             }
          }
        }
      }
    }
#pragma omp critical
    {
      // save
      tri_nodes_.insert(tri_nodes_.end(), tri_nodes_buffer.begin(),
                        tri_nodes_buffer.end());
    }
  }

  // clear
  un_points_.clear();
  vp_directions.clear();
  backproj_planes.clear();
  centers.clear();

  // save the mapping from the line match to tri_node
  for (tri_t tri_idx = 0; tri_idx < tri_nodes_.size(); tri_idx++) {
    image_t image_id1 = std::get<1>(tri_nodes_[tri_idx]);
    line2d_t line2d_idx1 = std::get<2>(tri_nodes_[tri_idx]);
    image_t image_id2 = std::get<3>(tri_nodes_[tri_idx]);
    line2d_t line2d_idx2 = std::get<4>(tri_nodes_[tri_idx]);
    THROW_CHECK_LT(image_id1, image_id2);

    if (config_.debug_mode) {
      THROW_CHECK_LT(line2d_idx1, matches_.at(image_id1).size());
      THROW_CHECK_EQ(
          matches_.at(image_id1).at(line2d_idx1).at(image_id2).at(line2d_idx2),
          kInvalidTriIdx);
      THROW_CHECK_LT(line2d_idx2, matches_.at(image_id2).size());
      THROW_CHECK_EQ(
          matches_.at(image_id2).at(line2d_idx2).at(image_id1).at(line2d_idx1),
          kInvalidTriIdx);
    }

    matches_.at(image_id1)[line2d_idx1].at(image_id2)[line2d_idx2] = tri_idx;
    matches_.at(image_id2)[line2d_idx2].at(image_id1)[line2d_idx1] = tri_idx;
  }

  // create cache
  CreateCache();

  PrintHeading2("[report] construct 3D LS hypothesis nodes");
  std::cout << "# 3D LS hypothesis nodes = " << tri_nodes_.size() << std::endl;
  std::cout << "# 3D LS hypothesis from line-line solver = " << num_line_line
            << std::endl;
  std::cout << "# 3D LS hypothesis from two-points solver = " << num_two_points
            << std::endl;
  std::cout << "# 3D LS hypothesis from point-vp solver = " << num_point_vp
            << std::endl;
  std::cout << "# use two-points solver matches / # all matches = "
            << num_use_two_points_solver_matches << " / " << num_matches_
            << " = "
            << double(num_use_two_points_solver_matches) / double(num_matches_)
            << std::endl;
  std::cout << "# use point-vp solver matches / # all matches = "
            << num_use_point_vp_solver_matches << " / " << num_matches_ << " = "
            << double(num_use_point_vp_solver_matches) / double(num_matches_)
            << std::endl;
  std::cout << "# use line-line solver matches / # all matches = "
            << num_use_lineline_solver_matches << " / " << num_matches_ << " = "
            << double(num_use_lineline_solver_matches) / double(num_matches_)
            << std::endl;
  std::cout << "# degenerate matches / # all matches = " << num_degeneracy_line
            << " / " << num_matches_ << " = "
            << double(num_degeneracy_line) / double(num_matches_) << std::endl;
}

void LineMapper::ConstructTriEdges() {
  // clear
  all_hypotheses_.clear();
  all_hypotheses_.shrink_to_fit();

  // configure linker
  LineLinker linker = construct_3D_line_edges_linker_;
  linker.linker_2d.config.set_to_default();
  linker.linker_3d.config.set_to_spatial_merging();
  const auto& linker2d = linker.linker_2d;
  const double score2d_th = linker.linker_2d.config.score_th;
  const double score3d_th = linker.linker_3d.config.score_th;
  const double score_th = std::min(score2d_th, score3d_th);

  // construct unidirectional edges between tri_nodes (tri_idx1 < tri_idx2)
  const size_t num_tri_nodes = tri_nodes_.size();
  tri_edges_.clear();
  tri_edges_.shrink_to_fit();
  tri_edges_.resize(num_tri_nodes);
  progressbar bar(num_tri_nodes);
#pragma omp parallel for
  for (tri_t tri_idx = 0; tri_idx < num_tri_nodes; tri_idx++) {
    bar.update();
    const auto& tri_node = tri_nodes_[tri_idx];
    const Line3d& line3d = std::get<0>(tri_node);
    const image_t image_id1 = std::get<1>(tri_node);
    const line2d_t line2d_idx1 = std::get<2>(tri_node);
    const image_t image_id2 = std::get<3>(tri_node);
    const line2d_t line2d_idx2 = std::get<4>(tri_node);
    const auto& matches1 = matches_.at(image_id1)[line2d_idx1];
    const auto& matches2 = matches_.at(image_id2)[line2d_idx2];
    const auto& view1 = camviews_.at(image_id1);
    const auto& view2 = camviews_.at(image_id2);
    const auto& line2d1 = all_lines2d_.at(image_id1)[line2d_idx1];
    const auto& line2d2 = all_lines2d_.at(image_id2)[line2d_idx2];
    auto& tri_edges = tri_edges_.at(tri_idx);

    for (auto it = matches1.begin(); it != matches1.end(); it++) {
      const image_t image_id3 = it->first;
      THROW_CHECK_NE(image_id3, image_id1);
      const auto& view3 = camviews_.at(image_id3);
      const auto& all_lines2d3 = all_lines2d_.at(image_id3);
      const Line2d proj_line2d_on_view3 = line3d.projection(view3);

      for (const auto& match : it->second) {
        const line2d_t line2d_idx3 = match.first;
        const tri_t ng_tri_idx = match.second;
        if (ng_tri_idx == kInvalidTriIdx) continue;
        THROW_CHECK_GE(ng_tri_idx, 0);
        THROW_CHECK_LT(ng_tri_idx, num_tri_nodes);
        if (tri_idx > ng_tri_idx) continue;
        if (tri_idx == ng_tri_idx) {
          THROW_CHECK_EQ(image_id3, image_id2);
          THROW_CHECK_EQ(line2d_idx3, line2d_idx2);
          continue;
        }

        // compute similarity score between line3d and ng_line3d in 3D and 2D
        const Line3d& ng_line3d = std::get<0>(tri_nodes_[ng_tri_idx]);
        double score3d = linker.compute_score_3d(line3d, ng_line3d);
        if (score3d == 0) continue;
        double score2d_12_to_3 = linker2d.compute_reprojection_score(
            all_lines2d3[line2d_idx3], proj_line2d_on_view3);
        if (score2d_12_to_3 == 0) continue;
        double score2d_13_to_2 = linker2d.compute_reprojection_score(
            line2d2, ng_line3d.projection(view2));
        if (score2d_13_to_2 == 0) continue;
        double score = std::min({score3d, score2d_12_to_3, score2d_13_to_2});
        THROW_CHECK_GE(score, score_th);

        // add an unidirectional edge and set weight as the similarity score
        THROW_CHECK_LT(tri_idx, ng_tri_idx);
        if (config_.debug_mode) {
          // test the edge (tri_idx ---> ng_tri_idx) was not added before
          THROW_CHECK(tri_edges.find(ng_tri_idx) == tri_edges.end());
        }
        tri_edges.emplace(ng_tri_idx, score);
      }
    }

    for (auto it = matches2.begin(); it != matches2.end(); it++) {
      const image_t image_id3 = it->first;
      THROW_CHECK_NE(image_id3, image_id2);
      const auto& view3 = camviews_.at(image_id3);

      const auto& all_lines2d3 = all_lines2d_.at(image_id3);
      const Line2d proj_line2d_on_view3 = line3d.projection(view3);

      for (const auto& match : it->second) {
        const line2d_t line2d_idx3 = match.first;
        const tri_t ng_tri_idx = match.second;
        if (ng_tri_idx == kInvalidTriIdx) continue;
        THROW_CHECK_GE(ng_tri_idx, 0);
        THROW_CHECK_LT(ng_tri_idx, num_tri_nodes);
        if (tri_idx > ng_tri_idx) continue;
        if (tri_idx == ng_tri_idx) {
          THROW_CHECK_EQ(image_id3, image_id1);
          THROW_CHECK_EQ(line2d_idx3, line2d_idx1);
          continue;
        }

        // compute similarity score between line3d and ng_line3d in 3D and 2D
        const Line3d& ng_line3d = std::get<0>(tri_nodes_[ng_tri_idx]);
        double score3d = linker.compute_score_3d(line3d, ng_line3d);
        if (score3d == 0) continue;
        double score2d_12_to_3 = linker2d.compute_reprojection_score(
            all_lines2d3[line2d_idx3], proj_line2d_on_view3);
        if (score2d_12_to_3 == 0) continue;
        double score2d_23_to_1 = linker2d.compute_reprojection_score(
            line2d1, ng_line3d.projection(view1));
        if (score2d_23_to_1 == 0) continue;
        double score = std::min({score3d, score2d_12_to_3, score2d_23_to_1});
        THROW_CHECK_GE(score, score_th);

        // add an unidirectional edge and set weight as the similarity score
        THROW_CHECK_LT(tri_idx, ng_tri_idx);
        if (config_.debug_mode) {
          // test the edge (tri_idx ---> ng_tri_idx) was not added before
          THROW_CHECK(tri_edges.find(ng_tri_idx) == tri_edges.end());
        }
        tri_edges.emplace(ng_tri_idx, score);
      }
    }
  }

  // make `tri_edges_` bidirectional
  for (tri_t tri_idx = 0; tri_idx < tri_edges_.size(); tri_idx++) {
    const auto& tri_edges = tri_edges_[tri_idx];
    for (const auto& edge : tri_edges) {
      const tri_t ng_tri_idx = edge.first;
      THROW_CHECK_NE(tri_idx, ng_tri_idx);
      if (tri_idx > ng_tri_idx) continue;
      if (config_.debug_mode) {
        // test the edge (ng_tri_idx ---> tri_idx) was not added before
        // note: tri_idx < ng_tri_idx
        THROW_CHECK(tri_edges_[ng_tri_idx].find(tri_idx) ==
                    tri_edges_[ng_tri_idx].end());
      }
      // add the edge (ng_tri_idx ---> tri_idx)
      tri_edges_[ng_tri_idx].emplace(tri_idx, edge.second);
    }
  }

  ////////////////////////////////
  // debug
  if (config_.debug_mode) {
    // test no self connected node in tri_edges_
    for (tri_t tri_idx = 0; tri_idx < tri_edges_.size(); tri_idx++) {
      const auto& ng_tris = tri_edges_[tri_idx];
      THROW_CHECK(ng_tris.find(tri_idx) == ng_tris.end());
    }
    // test both edge and score are bidirectional
    for (tri_t tri_idx = 0; tri_idx < tri_edges_.size(); tri_idx++) {
      const auto& ng_tris = tri_edges_.at(tri_idx);
      for (const auto& pair : ng_tris) {
        const auto& ng_ng_tris = tri_edges_.at(pair.first);
        // test edge is bidirectional
        THROW_CHECK(ng_ng_tris.find(tri_idx) != ng_ng_tris.end());
        // test score is bidirectional
        THROW_CHECK_EQ(ng_ng_tris.at(tri_idx), pair.second);
        THROW_CHECK_GE(pair.second, score_th);
      }
    }
  }
  ////////////////////////////////

  // init
  remained_nodes_ = 0;
  score_table_.clear();
  score_table_.reserve(tri_nodes_.size());
  flags_.clear();
  flags_.resize(tri_nodes_.size(), 1);
  THROW_CHECK_EQ(tri_edges_.size(), tri_nodes_.size());
  for (tri_t tri_idx = 0; tri_idx < tri_edges_.size(); tri_idx++) {
    double score = 0.0;
    for (const auto& pair : tri_edges_[tri_idx]) {
      score += pair.second;
    }
    if (score == 0) {
      flags_[tri_idx] = 0;
    } else {
      score_table_[tri_idx] = score;
      scores_.insertOrUpdate(tri_idx, score);
      remained_nodes_++;
    }
  }
  if (config_.debug_mode) {
    THROW_CHECK_EQ(flags_.size(), tri_nodes_.size());
    for (size_t i = 0; i < flags_.size(); i++) {
      if (flags_[i]) {
        THROW_CHECK(score_table_.find(i) != score_table_.end());
        THROW_CHECK(scores_.exist_key(i));
      } else if (!flags_[i]) {
        THROW_CHECK(score_table_.find(i) == score_table_.end());
        THROW_CHECK(!scores_.exist_key(i));
      } else {
        throw std::runtime_error("error");
      }
    }
    THROW_CHECK(scores_.check());
  }

  {
    size_t n_edges = 0;
    size_t valid_n_nodes = 0;
    for (tri_t tri_idx = 0; tri_idx < tri_edges_.size(); tri_idx++) {
      size_t n_e = tri_edges_[tri_idx].size();
      n_edges += n_e;
      if (n_e > 0) valid_n_nodes++;
    }
    n_edges /= 2;
    PrintHeading2("[report] construct edges");
    std::cout << "# all 3D line nodes = " << tri_nodes_.size() << std::endl;
    std::cout << "# nodes with at least 1 degree = " << valid_n_nodes
              << std::endl;
    std::cout << "# edges = " << n_edges << std::endl;
  }
}

line3d_t LineMapper::FindAndAddBestLineTrack() {
  const double score_th =
      std::max(construct_3D_line_edges_linker_.linker_3d.config.score_th,
               construct_3D_line_edges_linker_.linker_2d.config.score_th);

  size_t num_filtered_nodes = 0;

  // filter invalided nodes
  std::queue<tri_t> q;
  // in the first incremental iteration, new_track_id_ == kInvalidLine3dId
  if (new_track_id_ != kInvalidLine3dId) {
    const auto& new_track = tracks_.at(new_track_id_);
    for (size_t i = 0; i < new_track.image_id_list.size(); i++) {
      const image_t image_id = new_track.image_id_list[i];
      const line2d_t line2d_idx = new_track.line_id_list[i];
      const auto& matches = matches_.at(image_id).at(line2d_idx);
      for (const auto& pair1 : matches) {
        for (const auto& pair2 : pair1.second) {
          if (pair2.second != kInvalidTriIdx && flags_[pair2.second]) {
            // filter invalided node
            flags_[pair2.second] = 0;
            score_table_.erase(pair2.second);
            q.push(pair2.second);
            num_filtered_nodes++;
          }
        }
      }
    }
  }

  if (config_.debug_mode) {
    std::set<tri_t> set1;
    std::queue<tri_t> q2 = q;
    while (!q2.empty()) {
      set1.insert(q2.front());
      q2.pop();
    }
    THROW_CHECK_EQ(set1.size(), q.size());
    THROW_CHECK_EQ(flags_.size(), tri_nodes_.size());
    for (size_t i = 0; i < flags_.size(); i++) {
      if (flags_[i] == 1) {
        THROW_CHECK(score_table_.find(i) != score_table_.end());
      } else if (flags_[i] == 0) {
        THROW_CHECK(score_table_.find(i) == score_table_.end());
      } else {
        throw std::runtime_error("error");
      }
    }
  }

  // iteratively filter nodes and update score in score_table_
  while (!q.empty()) {
    const tri_t invalid_tri_idx = q.front();
    q.pop();
    const auto& ng_nodes = tri_edges_.at(invalid_tri_idx);
    for (const auto& ng_node : ng_nodes) {
      const tri_t ng_tri_idx = ng_node.first;
      if (!flags_.at(ng_tri_idx)) continue;
      auto& score = score_table_.at(ng_tri_idx);
      score -= ng_node.second;
      if (score < score_th) {  // ng_node is isolated now
        THROW_CHECK_LT(std::abs(score), EPS);
        flags_[ng_tri_idx] = 0;
        score_table_.erase(ng_tri_idx);
        q.push(ng_tri_idx);
        num_filtered_nodes++;
      }
    }
  }

  std::cout << "# filtered nodes = " << num_filtered_nodes << std::endl;
  remained_nodes_ -= num_filtered_nodes;
  std::cout << "# remained nodes = " << remained_nodes_ << std::endl;
  std::cout << "# all nodes = " << tri_nodes_.size() << std::endl;

  if (config_.debug_mode) {
    THROW_CHECK_EQ(flags_.size(), tri_nodes_.size());
    size_t n1 = 0;
    for (size_t i = 0; i < flags_.size(); i++) {
      if (flags_[i] == 1) {
        n1++;
        THROW_CHECK(score_table_.find(i) != score_table_.end());
      } else if (flags_[i] == 0) {
        THROW_CHECK(score_table_.find(i) == score_table_.end());
      } else {
        throw std::runtime_error("error");
      }
    }
    THROW_CHECK_EQ(n1, remained_nodes_);
  }

  // find the best track
  tri_t best_idx = kInvalidTriIdx;
  double max_score = -1.0;
  for (const auto& pair : score_table_) {
    if (pair.second > max_score) {
      max_score = pair.second;
      best_idx = pair.first;
    }
  }

  // stop incremental iteration
  if (best_idx == kInvalidTriIdx || max_score < config_.min_score) {
    stop_ = true;
    PrintHeading1("stop incremental iteration");

    // clear
    matches_.clear();
    tri_nodes_.clear();
    tri_nodes_.shrink_to_fit();
    tri_edges_.clear();
    tri_edges_.shrink_to_fit();
    simple_line2d_cache_list_.clear();
    simple_line2d_cache_list_.shrink_to_fit();
    image_line2d_to_simple_cache_.clear();
    simple_cache_to_image_line2d_.clear();
    simple_cache_to_image_line2d_.shrink_to_fit();
    score_table_.clear();
    scores_.clear();
    flags_.clear();
    flags_.shrink_to_fit();

    return kInvalidLine3dId;
  }

  // build the best track
  const auto& tri_node = tri_nodes_.at(best_idx);
  const Line3d& line3d = std::get<0>(tri_node);
  const image_t image_id1 = std::get<1>(tri_node);
  const line2d_t line2d_idx1 = std::get<2>(tri_node);
  const image_t image_id2 = std::get<3>(tri_node);
  const line2d_t line2d_idx2 = std::get<4>(tri_node);
  LineTrack track;
  track.line = line3d;
  track.image_id_list.push_back(image_id1);
  track.line_id_list.push_back(line2d_idx1);
  track.line2d_list.push_back(all_lines2d_.at(image_id1).at(line2d_idx1));
  track.image_id_list.push_back(image_id2);
  track.line_id_list.push_back(line2d_idx2);
  track.line2d_list.push_back(all_lines2d_.at(image_id2).at(line2d_idx2));

  // add the best track to the 3D LS map
  line3d_t track_id = AddTrack(track);
  new_track_id_ = track_id;

  return track_id;
}

line3d_t LineMapper::FindAndAddBestLineTrackFastFindMax() {
  const double score_th =
      std::max(construct_3D_line_edges_linker_.linker_3d.config.score_th,
               construct_3D_line_edges_linker_.linker_2d.config.score_th);

  size_t num_filtered_nodes = 0;

  // filter invalided nodes
  std::queue<tri_t> q;
  // in the first incremental iteration, new_track_id_ == kInvalidLine3dId
  if (new_track_id_ != kInvalidLine3dId) {
    const auto& new_track = tracks_.at(new_track_id_);
    for (size_t i = 0; i < new_track.image_id_list.size(); i++) {
      const image_t image_id = new_track.image_id_list[i];
      const line2d_t line2d_idx = new_track.line_id_list[i];
      const auto& matches = matches_.at(image_id).at(line2d_idx);
      for (const auto& pair1 : matches) {
        for (const auto& pair2 : pair1.second) {
          if (pair2.second != kInvalidTriIdx && flags_[pair2.second]) {
            // filter invalided node
            flags_[pair2.second] = 0;
            scores_.erase(pair2.second);
            q.push(pair2.second);
            num_filtered_nodes++;
          }
        }
      }
    }
  }

  if (config_.debug_mode) {
    std::set<tri_t> set1;
    std::queue<tri_t> q2 = q;
    while (!q2.empty()) {
      set1.insert(q2.front());
      q2.pop();
    }
    THROW_CHECK_EQ(set1.size(), q.size());
    THROW_CHECK_EQ(flags_.size(), tri_nodes_.size());
    for (size_t i = 0; i < flags_.size(); i++) {
      if (flags_[i] == 1) {
        THROW_CHECK(scores_.exist_key(i));
      } else if (flags_[i] == 0) {
        THROW_CHECK(!scores_.exist_key(i));
      } else {
        throw std::runtime_error("error");
      }
    }
    if (tracks_.size() % 500 == 0) {
      // scores_.check() is very slow
      THROW_CHECK(scores_.check());
    }
  }

  // iteratively filter nodes and update score in scores_.dataMap

  // (tri_idx, init score)
  std::unordered_map<size_t, double> remained_modified_nodes;
  std::vector<std::pair<size_t, double>> removed_nodes;
  while (!q.empty()) {
    const tri_t invalid_tri_idx = q.front();
    q.pop();
    const auto& ng_nodes = tri_edges_.at(invalid_tri_idx);
    for (const auto& ng_node : ng_nodes) {
      const tri_t ng_tri_idx = ng_node.first;
      if (!flags_.at(ng_tri_idx)) continue;
      double& score = scores_.get_value(ng_tri_idx);
      if (remained_modified_nodes.find(ng_tri_idx) ==
          remained_modified_nodes.end()) {
        remained_modified_nodes[ng_tri_idx] = score;  // initial score
      }
      score -= ng_node.second;

      if (config_.debug_mode) {
        THROW_CHECK_EQ(scores_.get_value(ng_tri_idx), score);
      }

      if (score < score_th) {  // ng_node is isolated now
        THROW_CHECK_LT(std::abs(score), EPS);
        flags_[ng_tri_idx] = 0;
        scores_.eraseDataMap(ng_tri_idx);
        removed_nodes.emplace_back(ng_tri_idx,
                                   remained_modified_nodes.at(ng_tri_idx));
        remained_modified_nodes.erase(ng_tri_idx);
        q.push(ng_tri_idx);
        num_filtered_nodes++;
      }
    }
  }

  std::cout << "# filtered nodes = " << num_filtered_nodes << std::endl;
  remained_nodes_ -= num_filtered_nodes;
  std::cout << "# remained nodes = " << remained_nodes_ << std::endl;
  std::cout << "# all_nodes = " << tri_nodes_.size() << std::endl;

  // remove keys in scores_.sortedSet
  for (const auto& pair : removed_nodes) {
    scores_.eraseSortedSet(pair.first, pair.second);
  }

  // update score in scores_.sortedSet
  for (const auto& pair : remained_modified_nodes) {
    double score = scores_.get_value(pair.first);
    scores_.UpdateSortedSet(pair.first, pair.second, score);
  }

  if (config_.debug_mode) {
    for (const auto& pair : remained_modified_nodes) {
      double score = scores_.get_value(pair.first);
      THROW_CHECK_GT(pair.second, score);
    }
    std::set<size_t> set1;
    for (const auto pair : removed_nodes) {
      set1.insert(pair.first);
    }
    THROW_CHECK_EQ(set1.size(), removed_nodes.size());
    THROW_CHECK_EQ(flags_.size(), tri_nodes_.size());
    size_t n1 = 0;
    for (size_t i = 0; i < flags_.size(); i++) {
      if (flags_[i] == 1) {
        n1++;
        THROW_CHECK(scores_.exist_key(i));
      } else if (flags_[i] == 0) {
        THROW_CHECK(!scores_.exist_key(i));
      } else {
        throw std::runtime_error("error");
      }
    }
    THROW_CHECK_EQ(n1, remained_nodes_);
    if (tracks_.size() % 500 == 0) {
      // note: scores_.check() is very slow
      THROW_CHECK(scores_.check());
      std::cout << "pass score check" << std::endl;
    }
  }

  // find the best track
  tri_t best_idx = kInvalidTriIdx;
  double max_score = -1.0;
  if (!scores_.empty()) {
    auto pair = scores_.findMaxValue();
    best_idx = pair.first;
    max_score = pair.second;
    THROW_CHECK_NE(best_idx, kInvalidTriIdx);
    THROW_CHECK_GE(max_score, score_th);
  }

  // stop incremental iteration
  if (best_idx == kInvalidTriIdx || max_score < config_.min_score) {
    stop_ = true;
    PrintHeading1("stop incremental iteration");

    // clear
    matches_.clear();
    tri_nodes_.clear();
    tri_nodes_.shrink_to_fit();
    tri_edges_.clear();
    tri_edges_.shrink_to_fit();
    simple_line2d_cache_list_.clear();
    simple_line2d_cache_list_.shrink_to_fit();
    image_line2d_to_simple_cache_.clear();
    simple_cache_to_image_line2d_.clear();
    simple_cache_to_image_line2d_.shrink_to_fit();
    score_table_.clear();
    scores_.clear();
    flags_.clear();
    flags_.shrink_to_fit();

    return kInvalidLine3dId;
  }

  // build the best track
  const auto& tri_node = tri_nodes_.at(best_idx);
  const Line3d& line3d = std::get<0>(tri_node);
  const image_t image_id1 = std::get<1>(tri_node);
  const line2d_t line2d_idx1 = std::get<2>(tri_node);
  const image_t image_id2 = std::get<3>(tri_node);
  const line2d_t line2d_idx2 = std::get<4>(tri_node);
  LineTrack track;
  track.line = line3d;
  track.image_id_list.push_back(image_id1);
  track.line_id_list.push_back(line2d_idx1);
  track.line2d_list.push_back(all_lines2d_.at(image_id1).at(line2d_idx1));
  track.image_id_list.push_back(image_id2);
  track.line_id_list.push_back(line2d_idx2);
  track.line2d_list.push_back(all_lines2d_.at(image_id2).at(line2d_idx2));

  // add the best track to the 3D LS map
  line3d_t track_id = AddTrack(track);
  new_track_id_ = track_id;

  return track_id;
}

void LineMapper::ExtendTrack(const line3d_t track_id) {
  PrintHeading1("extend the new track");

  // iteratively extend the new line track using untriangulated 2D LSs

  const auto& linker3d = extend_tracks_linker_.linker_3d;
  const auto& linker2d = extend_tracks_linker_.linker_2d;

  // init
  size_t n_prev_tracks = tracks_.size();
  size_t n_prev_tri_lines2d = 0;
  size_t n_new_tri_lines2d = 0;
  size_t n_prev_local_tri_lines2d = 0;
  if (config_.debug_mode) {
    for (size_t i = 0; i < simple_line2d_cache_list_.size(); i++) {
      if (simple_line2d_cache_list_[i].line2d_ptr->track_id !=
          kInvalidLine3dId) {
        n_prev_tri_lines2d++;
      }
    }
  }

  auto& track = tracks_.at(track_id);  // will be updated
  double extended_length = std::numeric_limits<double>::max();
  std::unordered_set<size_t> verified_caches;
  const double init_length = track.line.length();

  size_t iter = 0;
  while (true) {
    // extend the new line track using the matching graph

    // find untriangulated matched 2D LSs of track elements
    std::unordered_set<size_t> to_be_verified_cache_set;
    std::vector<std::pair<size_t, tri_t>> to_be_verified_cache_tri_list;

    std::unordered_set<size_t> invalid_cache_indices;
    for (size_t i = 0; i < track.image_id_list.size(); i++) {
      const image_t image_id = track.image_id_list[i];
      const line2d_t line2d_idx = track.line_id_list[i];
      const size_t cache_idx = image_line2d_to_simple_cache_.at(
          std::make_pair(image_id, line2d_idx));
      const auto& cache = simple_line2d_cache_list_[cache_idx];
      for (const auto& pair : cache.ng_cache_tri_list) {
        const size_t ng_cache_idx = pair.first;
        const tri_t tri_idx = pair.second;

        THROW_CHECK_NE(ng_cache_idx, cache_idx);
        if (config_.require_tri_node) {
          if (tri_idx == kInvalidTriIdx) continue;
        }

        const auto& ng_cache = simple_line2d_cache_list_[ng_cache_idx];

        // if ng 2D LS was triangulated, continue
        if (ng_cache.line2d_ptr->track_id != kInvalidLine3dId) continue;
        // if the 3D LS was extended and ng_cache_idx was verified before,
        // re-verify!
        if ((extended_length < 1e-6) &&
            (verified_caches.find(ng_cache_idx) != verified_caches.end())) {
          continue;
        }
        if (config_.use_collinearity2D) {
          if (invalid_cache_indices.find(ng_cache_idx) !=
              invalid_cache_indices.end())
            continue;

          bool flag = true;
          for (size_t j = 0; j < track.image_id_list.size(); j++) {
            const size_t inner_cache_idx = image_line2d_to_simple_cache_.at(
                std::make_pair(track.image_id_list[j], track.line_id_list[j]));
            const auto& inner_cache =
                simple_line2d_cache_list_[inner_cache_idx];
            if (inner_cache.view_ptr == ng_cache.view_ptr) {
              if (TestDifferentInfLines3D(*(inner_cache.line2d_ptr),
                                          *(ng_cache.line2d_ptr), 0.0, 2.0,
                                          2.0)) {
                invalid_cache_indices.insert(ng_cache_idx);
                flag = false;
                break;
              }
            }
          }
          if (!flag) continue;
        }

        if (config_.test_tri_node) {
          to_be_verified_cache_tri_list.push_back(pair);
        } else {
          to_be_verified_cache_set.insert(ng_cache_idx);
        }
        verified_caches.insert(ng_cache_idx);
      }
    }

    // note: `track` will be modified in each while loop
    const Line3d& line3d = track.line;

    // verify each 2D LS cache
    std::unordered_set<size_t> new_tri_caches;
    if (config_.test_tri_node) {
      for (size_t i = 0; i < to_be_verified_cache_tri_list.size(); i++) {
        const size_t cache_idx = to_be_verified_cache_tri_list[i].first;
        const size_t tri_idx = to_be_verified_cache_tri_list[i].second;
        const auto& cache = simple_line2d_cache_list_[cache_idx];
        const Line2d& line2d = *(cache.line2d_ptr);
        THROW_CHECK_EQ(line2d.track_id, kInvalidLine3dId);
        const auto& view = *(cache.view_ptr);

        // check reprojection
        Line2d proj_line2d = line3d.projection(view);
        double score2d =
            linker2d.compute_reprojection_score(line2d, proj_line2d);
        if (score2d == 0) continue;

        // check 3D proximity
        if (tri_idx != kInvalidTriIdx) {
          const Line3d& two_view_line3d = std::get<0>(tri_nodes_.at(tri_idx));
          double score3d = linker3d.compute_score(two_view_line3d, line3d);
          if (score3d == 0) continue;
        }

        // check sensitivity
        double sensitivity = line3d.sensitivity(view);
        if (sensitivity > config_.extend_track_sensitivity_threshold) {
          continue;
        }

        // check positive depth
        InfiniteLine3d inf_line3d(line3d);
        const V3D dir = inf_line3d.direction();
        const V3D p_ref = inf_line3d.point();
        InfiniteLine2d inf_line2d_proj = inf_line3d.projection(view);
        V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
        V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);
        // compute 3D endpoints of the 2D line
        V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
        double tstart = (pstart_3d - p_ref).dot(dir);
        pstart_3d = p_ref + dir * tstart;
        V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
        double tend = (pend_3d - p_ref).dot(dir);
        pend_3d = p_ref + dir * tend;
        double z_start = view.pose.projdepth(pstart_3d);
        double z_end = view.pose.projdepth(pend_3d);
        if (z_start < EPS || z_end < EPS) continue;

        new_tri_caches.insert(cache_idx);
      }
    } else {
      for (const auto& cache_idx : to_be_verified_cache_set) {
        const auto& cache = simple_line2d_cache_list_[cache_idx];
        const Line2d& line2d = *(cache.line2d_ptr);
        THROW_CHECK_EQ(line2d.track_id, kInvalidLine3dId);
        const auto& view = *(cache.view_ptr);

        // check reprojection
        Line2d proj_line2d = line3d.projection(view);
        double score2d =
            linker2d.compute_reprojection_score(line2d, proj_line2d);
        if (score2d == 0) continue;

        // check sensitivity
        double sensitivity = line3d.sensitivity(view);
        if (sensitivity > config_.extend_track_sensitivity_threshold) {
          continue;
        }

        // check positive depth
        InfiniteLine3d inf_line3d(line3d);
        const V3D dir = inf_line3d.direction();
        const V3D p_ref = inf_line3d.point();
        InfiniteLine2d inf_line2d_proj = inf_line3d.projection(view);
        V2D pstart_2d = inf_line2d_proj.point_projection(line2d.start);
        V2D pend_2d = inf_line2d_proj.point_projection(line2d.end);
        // compute 3D endpoints of the 2D line
        V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
        double tstart = (pstart_3d - p_ref).dot(dir);
        pstart_3d = p_ref + dir * tstart;
        V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
        double tend = (pend_3d - p_ref).dot(dir);
        pend_3d = p_ref + dir * tend;
        double z_start = view.pose.projdepth(pstart_3d);
        double z_end = view.pose.projdepth(pend_3d);
        if (z_start < EPS || z_end < EPS) continue;

        new_tri_caches.insert(cache_idx);
      }
    }

    for (const auto& cache_idx : new_tri_caches) {
      auto& cache = simple_line2d_cache_list_[cache_idx];
      // set the ng 2D LS as triangulated
      THROW_CHECK_EQ(cache.line2d_ptr->track_id, kInvalidLine3dId);
      cache.line2d_ptr->track_id = track_id;

      const image_t image_id = simple_cache_to_image_line2d_[cache_idx].first;
      const line2d_t line2d_idx =
          simple_cache_to_image_line2d_[cache_idx].second;

      track.image_id_list.push_back(image_id);
      track.line_id_list.push_back(line2d_idx);
      track.line2d_list.push_back(*(cache.line2d_ptr));
    }

    size_t n_new_tri_caches = new_tri_caches.size();
    n_new_tri_lines2d += n_new_tri_caches;
    if (config_.debug_mode) {
      if (iter == 0) {
        n_prev_local_tri_lines2d = n_prev_tri_lines2d;
      }
    }
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "track extension iter = " << iter++ << std::endl;
    std::cout << "# newly triangulated 2D lines = " << n_new_tri_caches
              << std::endl;
    if (config_.debug_mode) {
      std::cout << "# previous triangulated 2D lines = "
                << n_prev_local_tri_lines2d << std::endl;
      n_prev_local_tri_lines2d += n_new_tri_caches;
    }

    if (n_new_tri_caches == 0) break;

    // update 3D LS because track elements have been updated
    std::unordered_set<image_line2d_t> image_lines2d;
    for (size_t supp_idx = 0; supp_idx < track.image_id_list.size();
         supp_idx++) {
      image_lines2d.emplace(track.image_id_list[supp_idx],
                            track.line_id_list[supp_idx]);
    }
    Line3d extended_line3d = ExtendLine3d(
        line3d, image_lines2d, camviews_, all_lines2d_, config_.var2d,
        config_.num_outliers_aggregator_for_line_track_extension);
    extended_length = extended_line3d.length() - line3d.length();
    track.line = extended_line3d;

    ////////////////////////////////////
    // debug
    if (config_.debug_mode) {
      // test no repeat
      const auto& track = tracks_.at(track_id);
      std::vector<image_line2d_t> line2d_node_vec;
      for (size_t i = 0; i < track.image_id_list.size(); i++) {
        THROW_CHECK_EQ(track.image_id_list.size(), track.line_id_list.size());
        THROW_CHECK_EQ(track.image_id_list.size(), track.line2d_list.size());
        line2d_node_vec.emplace_back(track.image_id_list[i],
                                     track.line_id_list[i]);
      }
      std::set<image_line2d_t> line2d_node_set(line2d_node_vec.begin(),
                                               line2d_node_vec.end());
      THROW_CHECK_EQ(line2d_node_vec.size(), line2d_node_set.size());
    }
    ////////////////////////////////////
  }

  const double final_length = track.line.length();
  const double final_extended_length =
      final_length - init_length;  // may be negative due to outliers aggregator

  if (config_.debug_mode) {
    size_t n_tri_lines2d = 0;
    for (const auto& pair : all_lines2d_) {
      const image_t image_id = pair.first;
      for (line2d_t line2d_idx = 0; line2d_idx < pair.second.size();
           line2d_idx++) {
        if (pair.second[line2d_idx].track_id != kInvalidLine3dId) {
          n_tri_lines2d++;
        }
      }
    }
    THROW_CHECK_EQ(n_tri_lines2d, n_prev_tri_lines2d + n_new_tri_lines2d);

    // test 2D-3D correspondences
    const auto& track = tracks_.at(track_id);
    std::unordered_set<image_line2d_t> image_line2d_set;
    for (size_t i = 0; i < track.image_id_list.size(); i++) {
      image_line2d_set.emplace(track.image_id_list.at(i),
                               track.line_id_list.at(i));
    }
    THROW_CHECK_EQ(image_line2d_set.size(), track.image_id_list.size());
    for (const auto& pair : all_lines2d_) {
      for (size_t line2d_idx = 0; line2d_idx < pair.second.size();
           line2d_idx++) {
        const auto& line2d = pair.second[line2d_idx];
        if (line2d.track_id != kInvalidLine3dId) {
          // test the existence for the 3D line track
          THROW_CHECK(ExistLineTrack(line2d.track_id));
        }
      }
    }
  }

  THROW_CHECK_EQ(tracks_.size(), n_prev_tracks);

  PrintHeading2("[report] extend tracks");
  std::cout << "extended track_id = " << track_id << std::endl;
  std::cout << "extended length = " << final_extended_length << std::endl;
  if (config_.debug_mode) {
    std::cout << "# previous triangulated 2D lines = " << n_prev_tri_lines2d
              << std::endl;
    std::cout << "# current triangulated 2D lines / # all 2D lines = "
              << (n_prev_tri_lines2d + n_new_tri_lines2d) << " / "
              << num_lines2d_ << std::endl;
  }
  std::cout << "# total newly triangulated 2D lines = " << n_new_tri_lines2d
            << std::endl;
  std::cout << "# all detected 2D lines = " << num_lines2d_ << std::endl;
}

size_t LineMapper::MergeLineTracks() {
  const auto& linker3d = merge_tracks_linker_.linker_3d;
  const double th_angle2d = merge_tracks_linker_.linker_2d.config.th_angle;
  const double th_perp2d = merge_tracks_linker_.linker_2d.config.th_perp;

  std::vector<LineTrack> linetracks = GetLineTrackList();

  // delete all tracks
  auto track_ids = TrackIds();
  for (const auto track_id : track_ids) {
    DeleteTrack(track_id);
  }

  const size_t n_tracks = linetracks.size();

  // compute edges to remerge
  std::vector<std::tuple<size_t, size_t, double>>
      edges_tup_list;  // (track_id, track_id, score)
  std::vector<std::vector<std::tuple<size_t, size_t, double>>>
      edges_per_track_tup(n_tracks);  // (track_id, track_id, score)
  std::vector<int> active_ids;
  for (size_t i = 0; i < n_tracks; ++i) {
    if (linetracks[i].active) active_ids.push_back(i);
  }
  int n_active_ids = active_ids.size();
#pragma omp parallel for
  for (size_t k = 0; k < n_active_ids; ++k) {
    int i = active_ids[k];
    const Line3d& l1 = linetracks[i].line;
    for (size_t j = 0; j < n_tracks; ++j) {
      if (i == j) continue;
      if (n_active_ids == n_tracks) {
        if (i < j && (i + j) % 2 == 0) continue;
        if (i > j && (i + j) % 2 == 1) continue;
      }
      const Line3d& l2 = linetracks[j].line;
      // test similarity in 3D
      double score = linker3d.compute_score(l1, l2);
      if (score == 0) continue;
      if (i < j) {
        edges_per_track_tup[i].emplace_back(i, j, score);
      } else {
        edges_per_track_tup[i].emplace_back(j, i, score);
      }
    }
  }
  for (size_t i = 0; i < n_tracks; ++i) {
    edges_tup_list.insert(edges_tup_list.end(), edges_per_track_tup[i].begin(),
                          edges_per_track_tup[i].end());
  }
  edges_per_track_tup.clear();
  edges_per_track_tup.shrink_to_fit();
  std::map<std::pair<size_t, size_t>, double> edges_map;
  for (const auto& tup : edges_tup_list) {
    const auto pair = std::make_pair(std::get<0>(tup), std::get<1>(tup));
    if (edges_map.find(pair) == edges_map.end()) {
      edges_map.emplace(pair, std::get<2>(tup));
    } else {
      THROW_CHECK_EQ(edges_map.at(pair), std::get<2>(tup));
    }
  }
  edges_tup_list.clear();
  for (const auto& pair : edges_map) {
    edges_tup_list.emplace_back(pair.first.first, pair.first.second,
                                pair.second);
  }
  if (config_.debug_mode) {
    std::set<std::pair<size_t, size_t>> set1;
    for (const auto& tup : edges_tup_list) {
      set1.emplace(std::get<0>(tup), std::get<1>(tup));
    }
    THROW_CHECK_EQ(set1.size(), edges_tup_list.size());
  }

  // group connected components
  std::vector<int> parent_tracks(n_tracks, -1);
  std::vector<std::set<int>> tracks_in_group(n_tracks);
  std::vector<Line3d> avgline_in_track(n_tracks);
  std::vector<std::unordered_set<image_line2d_t>> lines2d_in_track(n_tracks);
#pragma omp parallel for
  for (size_t i = 0; i < n_tracks; ++i) {
    tracks_in_group[i].insert(i);
    const auto& track = linetracks[i];
    avgline_in_track[i] = track.line;
    for (size_t supp_idx = 0; supp_idx < track.image_id_list.size();
         supp_idx++) {
      lines2d_in_track[i].emplace(track.image_id_list[supp_idx],
                                  track.line_id_list[supp_idx]);
    }
  }

  // MST-like merging
  std::sort(edges_tup_list.begin(), edges_tup_list.end(),
            [](const std::tuple<size_t, size_t, double>& edge1,
               const std::tuple<size_t, size_t, double>& edge2) {
              return std::get<2>(edge1) > std::get<2>(edge2);
            });

  const size_t n_edges = edges_tup_list.size();
  for (size_t edge_idx = 0; edge_idx < n_edges; ++edge_idx) {
    const auto& edge = edges_tup_list[edge_idx];
    const size_t track_id1 = std::get<0>(edge);
    const size_t track_id2 = std::get<1>(edge);

    size_t root1 = union_find_get_root(track_id1, parent_tracks);
    size_t root2 = union_find_get_root(track_id2, parent_tracks);
    if (root1 != root2) {
      // try to merge the 3D LS of two groups
      const Line3d& line3d1 = avgline_in_track.at(root1);
      const Line3d& line3d2 = avgline_in_track.at(root2);

      // test1
      if (!linker3d.check_connection(line3d1, line3d2)) continue;

      const auto& image_lines2d1 = lines2d_in_track.at(root1);
      const auto& image_lines2d2 = lines2d_in_track.at(root2);

      // test2
      if (config_.use_collinearity2D) {
        bool flag = true;
        for (const auto& image_line2d1 : image_lines2d1) {
          const Line2d& line2d1 =
              all_lines2d_.at(image_line2d1.first)[image_line2d1.second];
          for (const auto& image_line2d2 : image_lines2d2) {
            if (image_line2d1.first == image_line2d2.first) {
              const Line2d& line2d2 =
                  all_lines2d_.at(image_line2d2.first)[image_line2d2.second];
              if (TestDifferentInfLines3D(line2d1, line2d2, 0.0, 2.0, 2.0)) {
                // cannot link
                flag = false;
                break;
              }
            }
          }
          if (!flag) break;
        }
        if (!flag) continue;
      }

      // test3
      std::vector<Line3d> lines3d = {line3d1, line3d2};
      std::vector<double> weights = {
          static_cast<double>(image_lines2d1.size()),
          static_cast<double>(image_lines2d2.size())};
      std::unordered_set<image_line2d_t> merged_image_lines2d = image_lines2d1;
      merged_image_lines2d.insert(image_lines2d2.begin(), image_lines2d2.end());
      THROW_CHECK_EQ(merged_image_lines2d.size(),
                     image_lines2d1.size() + image_lines2d2.size());
      auto res =
          MergeLines3d(lines3d, weights, merged_image_lines2d, camviews_,
                       all_lines2d_, th_angle2d, th_perp2d, config_.var2d,
                       config_.num_outliers_aggregator_for_line_track_merging,
                       "point_weighted_PCA");
      if (!res.second) continue;
      const Line3d merged_line3d = res.first;

      // union-find merging heuristic
      if (tracks_in_group[root1].size() < tracks_in_group[root2].size()) {
        parent_tracks[root1] = root2;
        // update tracks_in_group for root2
        tracks_in_group[root2].insert(tracks_in_group[root1].begin(),
                                      tracks_in_group[root1].end());
        tracks_in_group[root1].clear();
        lines2d_in_track[root2] = merged_image_lines2d;
        lines2d_in_track[root1].clear();
        avgline_in_track[root2] = merged_line3d;
      } else {
        parent_tracks[root2] = root1;
        // update tracks_in_group for root1
        tracks_in_group[root1].insert(tracks_in_group[root2].begin(),
                                      tracks_in_group[root2].end());
        tracks_in_group[root2].clear();
        lines2d_in_track[root1] = merged_image_lines2d;
        lines2d_in_track[root2].clear();
        avgline_in_track[root1] = merged_line3d;
      }
    }
  }

  if (config_.debug_mode) {
    for (size_t track_idx = 0; track_idx < n_tracks; ++track_idx) {
      if (parent_tracks[track_idx] == -1) {
        const auto& image_lines2d = lines2d_in_track[track_idx];
        std::unordered_set<image_line2d_t> image_lines2d_test;
        THROW_CHECK(!tracks_in_group[track_idx].empty());
        for (const auto& inner_track_idx : tracks_in_group[track_idx]) {
          const auto& track = linetracks[inner_track_idx];
          for (size_t i = 0; i < track.image_id_list.size(); ++i) {
            image_lines2d_test.emplace(track.image_id_list[i],
                                       track.line_id_list[i]);
          }
        }
        THROW_CHECK(image_lines2d == image_lines2d_test);
      } else {
        THROW_CHECK(tracks_in_group[track_idx].empty());
        THROW_CHECK(lines2d_in_track[track_idx].empty());
      }
    }
  }

  // collect roots
  std::vector<size_t> roots;
  for (size_t track_idx = 0; track_idx < n_tracks; ++track_idx) {
    if (parent_tracks[track_idx] == -1) {
      roots.push_back(track_idx);  // `track_idx` is root
    }
  }

  // create new tracks
  const size_t n_new_tracks = roots.size();
  std::vector<LineTrack> new_linetracks(n_new_tracks);
#pragma omp parallel for
  for (size_t new_track_idx = 0; new_track_idx < n_new_tracks;
       new_track_idx++) {
    auto& track = new_linetracks[new_track_idx];
    track.image_id_list.clear();
    track.line_id_list.clear();
    track.line2d_list.clear();
    const size_t root = roots[new_track_idx];
    for (const auto& image_line2d : lines2d_in_track.at(root)) {
      track.image_id_list.push_back(image_line2d.first);
      track.line_id_list.push_back(image_line2d.second);
      track.line2d_list.push_back(
          all_lines2d_.at(image_line2d.first).at(image_line2d.second));
    }
    track.line = avgline_in_track.at(root);
    if (tracks_in_group.at(root).size() == 1) {
      track.active = false;
    }
  }

  // add all new tracks
  for (const auto& track : new_linetracks) {
    AddTrack(track);
  }

  PrintHeading2("[report] merge tracks");
  std::cout << "# current tracks / # original tracks = "
            << new_linetracks.size() << " / " << linetracks.size() << std::endl;
  return linetracks.size() - new_linetracks.size();
}

line3d_t LineMapper::AddTrack(const LineTrack& track) {
  // add track
  line3d_t track_id = num_added_tracks_++;
  THROW_CHECK(!ExistLineTrack(track_id));
  tracks_.emplace(track_id, track);

  // set 2D LS as triangulated
  // TODO: Sync to track.line2d_list
  for (size_t i = 0; i < track.image_id_list.size(); i++) {
    auto& line2d =
        all_lines2d_.at(track.image_id_list[i])[track.line_id_list[i]];
    THROW_CHECK_EQ(line2d.track_id, kInvalidLine3dId);
    line2d.track_id = track_id;
  }

  return track_id;
}

void LineMapper::DeleteTrack(const line3d_t track_id) {
  THROW_CHECK(ExistLineTrack(track_id));
  const auto& track = tracks_.at(track_id);

  // reset 2D LS as untriangulated
  for (size_t i = 0; i < track.image_id_list.size(); i++) {
    auto& line2d =
        all_lines2d_.at(track.image_id_list[i])[track.line_id_list[i]];
    THROW_CHECK_EQ(line2d.track_id, track_id);
    line2d.track_id = kInvalidLine3dId;
  }

  // delete track
  tracks_.erase(track_id);
}

std::unordered_set<line3d_t> LineMapper::TrackIds() const {
  std::unordered_set<line3d_t> track_ids;
  track_ids.reserve(tracks_.size());

  for (const auto& track : tracks_) {
    track_ids.insert(track.first);
  }

  return track_ids;
}

void LineMapper::FilterLineTracksByReprojection() {
  PrintHeading1(
      "filter tracks by num images, NaN, reprojection, depth, sensitivity and "
      "ranges");
  LineLinker2d linker2d = filter_by_reprojection_linker2d_;
  linker2d.config.set_to_default();

  size_t old_num_tracks = tracks_.size();

  size_t n1_supports = 0;
  size_t n2_supports = 0;
  size_t n1_filtered_supports = 0;
  size_t n2_filtered_supports = 0;
  size_t n1_filtered_tracks = 0;
  size_t n2_filtered_tracks = 0;
  size_t n3_filtered_tracks = 0;

  auto track_ids = TrackIds();
  for (const line3d_t track_id : track_ids) {
    const auto& track = tracks_.at(track_id);

    const size_t num_supports = track.count_lines();
    n1_supports += num_supports;

    // [1] filter the track by num images

    if (track.count_images() < config_.min_support_images) {
      DeleteTrack(track_id);
      n1_filtered_supports += num_supports;
      n1_filtered_tracks++;
      continue;
    }

    // [2] filter the track by NaN

    if (std::isnan(track.line.start[0]) || std::isnan(track.line.end[0])) {
      DeleteTrack(track_id);
      n2_filtered_tracks++;
      continue;
    }

    // [3] filter track elements by reprojection, depth and sensitivity

    std::vector<CameraView> views;
    views.reserve(num_supports);
    for (size_t i = 0; i < track.image_id_list.size(); i++) {
      views.push_back(camviews_.at(track.image_id_list[i]));
    }

    std::vector<bool> results;
    limap::merging::CheckReprojectionAndDepth(
        linker2d, results, track, views, config_.th_depth,
        config_.filtering_sensitivity_threshold);
    THROW_CHECK_EQ(num_supports, results.size());
    size_t n_local_filtered_supports = 0;
    for (const auto& res : results) {
      if (!res) n_local_filtered_supports++;
    }
    n2_filtered_supports += n_local_filtered_supports;

    // [4] filter the track by ranges

    if (n_local_filtered_supports == 0) {
      if (config_.use_range_filter) {
        if (!limap::triangulation::test_line_inside_ranges(track.line,
                                                           ranges_)) {
          DeleteTrack(track_id);
          n3_filtered_tracks++;
        } else {
          // do not need filtering
          n2_supports += num_supports;
        }
      } else {
        // do not need filtering
        n2_supports += num_supports;
      }
    } else {
      auto new_tracks = ReBuildLineTracks(
          track, views, results,
          config_.num_outliers_aggregator_for_reprojection_filtering,
          config_.min_support_images, config_.var2d,
          config_.filter_incontinuous_line3d);

      DeleteTrack(track_id);

      for (const auto& new_track : new_tracks) {
        if (config_.use_range_filter) {
          if (!limap::triangulation::test_line_inside_ranges(new_track.line,
                                                             ranges_)) {
            n3_filtered_tracks++;
            continue;
          }
        }
        n2_supports += new_track.count_lines();
        AddTrack(new_track);
      }
    }
  }

  PrintHeading2(
      "[report] filter tracks by num images, nan, reprojection, depth, "
      "sensitivity and ranges");
  std::cout << "# filtered supports by num images = " << n1_filtered_supports
            << std::endl;
  std::cout << "# filtered tracks by num images = " << n1_filtered_tracks
            << std::endl;
  std::cout << "# filtered supports by reprojection, depth and sensitivity = "
            << n2_filtered_supports << std::endl;
  std::cout << "# filtered tracks by NaN = " << n2_filtered_tracks << std::endl;
  std::cout << "# filtered tracks by ranges = " << n3_filtered_tracks
            << std::endl;
  std::cout << "# previous supports = " << n1_supports << std::endl;
  std::cout << "# current supports = " << n2_supports << std::endl;
  std::cout << "# total filtered supports = " << n1_supports - n2_supports
            << std::endl;
  std::cout << "# previous tracks = " << old_num_tracks << std::endl;
  std::cout << "# current tracks = " << tracks_.size() << std::endl;
  std::cout << "# total filtered tracks = " << old_num_tracks - tracks_.size()
            << std::endl;
}

void LineMapper::FilterLineTracksBySensitivity() {
  PrintHeading1("filter tracks by sensitivity");

  size_t old_num_tracks = tracks_.size();

  size_t n1_supports = 0;
  size_t n2_supports = 0;

  std::vector<line3d_t> removed_tracks;
  for (const auto& pair : tracks_) {
    const line3d_t track_id = pair.first;
    const auto& track = pair.second;
    size_t num_supports = track.count_lines();
    n1_supports += num_supports;

    std::vector<bool> results;
    limap::merging::CheckSensitivity(results, track, imagecols_,
                                     config_.th_sv_angular_3d);
    THROW_CHECK_EQ(num_supports, results.size());

    std::set<int> support_images;
    for (size_t i = 0; i < num_supports; ++i) {
      if (results[i]) support_images.insert(track.image_id_list[i]);
    }

    int counter = support_images.size();
    if (counter < config_.th_sv_num_supports) {
      removed_tracks.push_back(track_id);
    } else {
      n2_supports += num_supports;
    }
  }

  // delete tracks
  for (line3d_t track_id : removed_tracks) {
    DeleteTrack(track_id);
  }

  PrintHeading2("[report] filter tracks by sensitivity");
  std::cout << "# previous supports = " << n1_supports << std::endl;
  std::cout << "# current supports = " << n2_supports << std::endl;
  std::cout << "# total filtered supports = " << n1_supports - n2_supports
            << std::endl;
  std::cout << "# previous tracks = " << old_num_tracks << std::endl;
  std::cout << "# current tracks = " << tracks_.size() << std::endl;
  std::cout << "# total filtered tracks = " << old_num_tracks - tracks_.size()
            << std::endl;
}

void LineMapper::FilterLineTracksByOverlap() {
  PrintHeading1("filter tracks by overlap");

  size_t old_num_tracks = tracks_.size();

  size_t n1_supports = 0;
  size_t n2_supports = 0;

  std::vector<line3d_t> removed_tracks;
  for (const auto& pair : tracks_) {
    const line3d_t track_id = pair.first;
    const auto& track = pair.second;
    size_t num_supports = track.count_lines();
    n1_supports += num_supports;

    std::set<int> support_images;
    for (size_t i = 0; i < num_supports; ++i) {
      Line2d line2d_proj =
          track.line.projection(camviews_.at(track.image_id_list[i]));
      const Line2d& line2d = track.line2d_list[i];
      double overlap = compute_overlap<Line2d>(line2d_proj, line2d);
      if (overlap >= config_.th_overlap)
        support_images.insert(track.image_id_list[i]);
    }

    int counter = support_images.size();
    if (counter < config_.th_overlap_num_supports) {
      removed_tracks.push_back(track_id);
    } else {
      n2_supports += num_supports;
    }
  }

  // delete tracks
  for (line3d_t track_id : removed_tracks) {
    DeleteTrack(track_id);
  }

  PrintHeading2("[report] filter tracks by overlap");
  std::cout << "# previous supports = " << n1_supports << std::endl;
  std::cout << "# current supports = " << n2_supports << std::endl;
  std::cout << "# total filtered supports = " << n1_supports - n2_supports
            << std::endl;
  std::cout << "# previous tracks = " << old_num_tracks << std::endl;
  std::cout << "# current tracks = " << tracks_.size() << std::endl;
  std::cout << "# total filtered tracks = " << old_num_tracks - tracks_.size()
            << std::endl;
}

void LineMapper::FilterLineTracksByNumImages(const int min_support_images) {
  PrintHeading1("filter tracks by num images");

  size_t old_num_tracks = tracks_.size();

  size_t n1_supports = 0;
  size_t n2_supports = 0;

  std::vector<line3d_t> removed_tracks;
  for (const auto& pair : tracks_) {
    const line3d_t track_id = pair.first;
    const auto& track = pair.second;
    size_t num_supports = track.count_lines();
    n1_supports += num_supports;

    if (track.count_images() < min_support_images) {
      removed_tracks.push_back(track_id);
    } else {
      n2_supports += num_supports;
    }
  }

  // delete tracks
  for (line3d_t track_id : removed_tracks) {
    DeleteTrack(track_id);
  }

  PrintHeading2("[report] filter tracks by num images");
  std::cout << "# previous supports = " << n1_supports << std::endl;
  std::cout << "# current supports = " << n2_supports << std::endl;
  std::cout << "# total filtered supports = " << n1_supports - n2_supports
            << std::endl;
  std::cout << "# previous tracks = " << old_num_tracks << std::endl;
  std::cout << "# current tracks = " << tracks_.size() << std::endl;
  std::cout << "# total filtered tracks = " << old_num_tracks - tracks_.size()
            << std::endl;
}

void LineMapper::SetLineTracks(const std::map<int, LineTrack>& tracks) {
  for (const auto& pair : tracks) {
    THROW_CHECK(ExistLineTrack(pair.first));
    tracks_.at(pair.first) = pair.second;
  }
}

std::vector<LineTrack> LineMapper::GetLineTrackList() const {
  std::vector<LineTrack> track_list;
  track_list.reserve(tracks_.size());

  for (const auto& pair : tracks_) {
    track_list.push_back(pair.second);
  }

  return track_list;
}

std::map<int, LineTrack> LineMapper::GetLineTrackMap() const {
  std::map<int, LineTrack> track_map;
  for (const auto& pair : tracks_) {
    track_map[pair.first] = pair.second;
  }

  return track_map;
}

std::vector<LineTrack> LineMapper::GetAllHypotheses() const {
  return all_hypotheses_;
}

}  // namespace limap