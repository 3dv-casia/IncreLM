#pragma once

#include "linemap/estimators/line_triangulation_estimator.h"
#include "linemap/mapping/functions.h"
#include "linemap/mapping/line_triangulation.h"
#include "linemap/util/types.h"

#include <map>
#include <unordered_map>
#include <vector>

#include "base/graph.h"
#include "base/image_collection.h"
#include "base/line_linker.h"
#include "base/linebase.h"
#include "structures/pl_bipartite.h"
#include "vplib/vpbase.h"

namespace limap {

struct LineMapperConfig {
 public:
  LineMapperConfig() {}
  LineMapperConfig(py::dict dict);

  // enable debug
  bool debug_mode = false;

  // add halfpix
  bool add_halfpix = false;

  // the min 2D line segment length
  double min_length_2d = 0.0;  // in pixels

  // uncertainty of line detection
  double var2d = 4.0;  // in pixels

  // line linker
  LineLinker2dConfig construct_3D_line_edges_linker2d_config;
  LineLinker3dConfig construct_3D_line_edges_linker3d_config;

  LineLinker2dConfig filter_by_reprojection_linker2d_config;

  LineLinker2dConfig merge_tracks_linker2d_config;
  LineLinker3dConfig merge_tracks_linker3d_config;

  LineLinker2dConfig extend_tracks_linker2d_config;
  LineLinker3dConfig extend_tracks_linker3d_config;

  // top $K$ matching
  int max_num_one_to_many_matches = 10;

  // whether to use range filter
  bool use_range_filter = true;

  // whether to use Line-Line solver
  bool use_line_line_solver = true;

  // min angle of two planes (only for statistics)
  double line_tri_angle_threshold = 1.0;  // in degrees

  // hybrid ransac options
  limap::estimators::line_triangulation::HybridLineTriangulationEstimatorOptions
      hybrid_ransac_options;

  // min score of the best node
  double min_score = 1.0;

  // extend tracks
  int num_outliers_aggregator_for_line_track_extension = 2;
  double extend_track_sensitivity_threshold = 50.0;
  bool require_tri_node = false;
  bool test_tri_node = false;

  // merging
  int num_outliers_aggregator_for_line_track_merging = 2;

  // filtering
  int num_outliers_aggregator_for_reprojection_filtering = 2;
  int min_support_images = 3;
  bool filter_incontinuous_line3d = true;
  double filtering_sensitivity_threshold = 90.0;  // in degrees
  double th_depth = EPS;
  double th_sv_angular_3d = 75.0;
  int th_sv_num_supports = 3;
  double th_overlap = 0.5;
  int th_overlap_num_supports = 3;

  // whether to use unshared sfm points for triangulation (inspired by CLMAP)
  bool use_unshared_points = false;
  double max_point_line_reproject_distance = 1.0;

  // whether to use 2D collinearity constraints (inspired by CLMAP)
  bool use_collinearity2D = false;
};

class LineMapper {
 public:
  // index for tri_nodes_
  typedef size_t tri_t;
  const tri_t kInvalidTriIdx = std::numeric_limits<tri_t>::max();

  struct SimpleLine2dCache {
    Line2d* line2d_ptr = nullptr;
    const CameraView* view_ptr = nullptr;
    std::vector<size_t> ng_lines2d_cache_idx;
    std::vector<std::pair<size_t, tri_t>> ng_cache_tri_list;
  };

  LineMapper() {}
  LineMapper(
      const LineMapperConfig& config, const ImageCollection& imagecols,
      const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
      const std::pair<V3D, V3D>& ranges);
  LineMapper(
      py::dict dict, const ImageCollection& imagecols,
      const std::unordered_map<image_t, std::vector<Line2d>>& all_lines2d,
      const std::pair<V3D, V3D>& ranges);

  inline void SetPLBipartite2d(
      const std::unordered_map<image_t, limap::structures::PL_Bipartite2d>&
          all_pl_bpt2d);
  inline void SetSfMPoints(
      const std::unordered_map<point3d_t, V3D>& sfm_points);
  inline void SetPointTracks(
      const std::unordered_map<point3d_t, PointTrack>& point_tracks);
  inline void SetVPResults(
      const std::unordered_map<image_t, limap::vplib::VPResult>& vpresults);
  inline void SetHybridRansacOptions(
      const limap::estimators::line_triangulation::
          HybridLineTriangulationEstimatorOptions& options);

  inline bool Stop() const;

  inline bool ExistLineTrack(const line3d_t& track_id) const;

  void SetLineTracks(const std::map<int, LineTrack>& tracks);

  void Initialize();

  // match 2D LS according to weak epipolar constraints
  void MatchLines2dByEpipolarIoU(
      const std::unordered_map<image_t, std::vector<image_t>>& all_neighbors,
      const int n_matches, const double th_IoU);

  // load one-to-many line matches for all 2D line segments and set the matching
  // as bi-directional
  void LoadMatches(
      const image_t& ref_image_id,
      const std::unordered_map<image_t, Eigen::MatrixXi>& init_ng_matches);

  size_t CountMatches();

  // construct tri nodes, i.e. triangulate 3D line segments from two views
  void ConstructTriNodes();

  // construct tri edges
  void ConstructTriEdges();

  // select the best tri node
  line3d_t FindAndAddBestLineTrack();
  line3d_t FindAndAddBestLineTrackFastFindMax();

  // extend a line track using untriangulated 2D line segments
  void ExtendTrack(const line3d_t track_id);

  // merge tracks
  size_t MergeLineTracks();

  // filter tracks
  void FilterLineTracksByReprojection();
  void FilterLineTracksBySensitivity();
  void FilterLineTracksByOverlap();
  void FilterLineTracksByNumImages(const int min_support_images);

  // outputs
  std::vector<LineTrack> GetLineTrackList() const;
  std::map<int, LineTrack> GetLineTrackMap() const;
  std::vector<LineTrack> GetAllHypotheses() const;

 private:
  void PreStartCheck() const;
  void OffsetHalfPixel();
  void CreateCache();

  // return all tracks id in `tracks_`
  std::unordered_set<line3d_t> TrackIds() const;

  // add a line track and return its track id
  line3d_t AddTrack(const LineTrack& track);

  // delete a line track
  void DeleteTrack(const line3d_t track_id);

  //////////////////////////////////////////////////////////////////////////////
  // Inputs
  //////////////////////////////////////////////////////////////////////////////

  // config of LineMapper
  LineMapperConfig config_;

  // linker
  LineLinker construct_3D_line_edges_linker_;
  LineLinker extend_tracks_linker_;
  LineLinker merge_tracks_linker_;
  const LineLinker2d filter_by_reprojection_linker2d_;

  // class that hold all cameras and images
  const ImageCollection imagecols_;

  // camviews
  std::unordered_map<image_t, CameraView> camviews_;

  // all Line2d for each image
  std::unordered_map<image_t, std::vector<Line2d>> all_lines2d_;

  // number of all detected 2D line segments.
  size_t num_lines2d_;

  // scene ranges
  const std::pair<V3D, V3D> ranges_;

  // 2D bipartite graph associating 2D point and 2D line segment
  std::shared_ptr<
      const std::unordered_map<image_t, limap::structures::PL_Bipartite2d>>
      all_pl_bpt2d_ptr_;

  // SfM points
  std::shared_ptr<const std::unordered_map<point3d_t, V3D>> sfm_points_ptr_;

  // SfM point tracks
  std::unordered_map<point3d_t, PointTrack> point_tracks_;

  // VP results for each image
  std::shared_ptr<const std::unordered_map<image_t, limap::vplib::VPResult>>
      vp_results_ptr_;

  //////////////////////////////////////////////////////////////////////////////
  // Outputs
  //////////////////////////////////////////////////////////////////////////////

  // matches_.at(image_id)[line2d_idx] = [ng_image_id -> [ng_line2d_idx ->
  // tri_idx]] (bidirectional)
  std::unordered_map<image_t,
                     std::vector<std::unordered_map<
                         image_t, std::unordered_map<line2d_t, tri_t>>>>
      matches_;

  // number of all line matches
  size_t num_matches_;

  // uncertainty of points3d, the keys are keys of point_tracks_ptr
  std::unordered_map<point3d_t, double> un_points_;

  // nodes of the 3D line segment hypothesis graph, which are all 3D line
  // segment hypotheses triangulated from two views
  std::vector<std::tuple<Line3d, image_t, line2d_t, image_t, line2d_t>>
      tri_nodes_;

  // edges of the 3D line segment hypothesis graph, which is bidirectional and
  // score >= min(linker_2d.config.score_th, linker_3d.config.score_th)
  std::vector<std::unordered_map<tri_t, double>> tri_edges_;

  // 3D line segment tracks
  std::unordered_map<line3d_t, LineTrack> tracks_;

  // total number of added line tracks, used to generate unique identifiers
  line3d_t num_added_tracks_;

  // cache
  std::vector<SimpleLine2dCache> simple_line2d_cache_list_;
  std::unordered_map<image_line2d_t, size_t> image_line2d_to_simple_cache_;
  std::vector<image_line2d_t> simple_cache_to_image_line2d_;

  // strength score for tri nodes, they will be updated during reconstruction
  std::unordered_map<tri_t, double> score_table_;
  MaxValuePair scores_;

  // 0: the tri node is filtered from graph, 1: otherwise.
  // the index of `flags_` is the index of tri nodes in `tri_nodes_`.
  std::vector<char> flags_;

  // each line track is a tri node
  std::vector<LineTrack> all_hypotheses_;

  // number of remained tri nodes during reconstruction
  size_t remained_nodes_;

  // new line track id
  line3d_t new_track_id_;

  // whether to stop incremental reconstruction
  bool stop_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

void LineMapper::SetPLBipartite2d(
    const std::unordered_map<image_t, limap::structures::PL_Bipartite2d>&
        all_pl_bpt2d) {
  all_pl_bpt2d_ptr_ = std::make_shared<
      const std::unordered_map<image_t, limap::structures::PL_Bipartite2d>>(
      all_pl_bpt2d);
}

void LineMapper::SetSfMPoints(
    const std::unordered_map<point3d_t, V3D>& sfm_points) {
  sfm_points_ptr_ =
      std::make_shared<const std::unordered_map<point3d_t, V3D>>(sfm_points);
}

void LineMapper::SetPointTracks(
    const std::unordered_map<point3d_t, PointTrack>& point_tracks) {
  point_tracks_ = point_tracks;
}

void LineMapper::SetVPResults(
    const std::unordered_map<image_t, limap::vplib::VPResult>& vp_results) {
  vp_results_ptr_ = std::make_shared<
      const std::unordered_map<image_t, limap::vplib::VPResult>>(vp_results);
}

void LineMapper::SetHybridRansacOptions(
    const limap::estimators::line_triangulation::
        HybridLineTriangulationEstimatorOptions& options) {
  config_.hybrid_ransac_options = options;
}

bool LineMapper::Stop() const { return stop_; }

bool LineMapper::ExistLineTrack(const line3d_t& track_id) const {
  return tracks_.find(track_id) != tracks_.end();
}

}  // namespace limap
