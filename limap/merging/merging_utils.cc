#include "merging/merging_utils.h"
#include "merging/aggregator.h"

#include "util/types.h"
#include "base/linebase.h"
#include "base/line_dists.h"

#include <third-party/progressbar.hpp>
#include <algorithm>

namespace limap
{

  namespace merging
  {

    std::vector<Line3d> SetUncertaintySegs3d(const std::vector<Line3d> &lines, const CameraView &view, const double var2d)
    {
      std::vector<Line3d> new_lines;
      for (size_t line_id = 0; line_id < lines.size(); ++line_id)
      {
        Line3d line = lines[line_id];
        line.set_uncertainty(line.computeUncertainty(view, var2d));
        new_lines.push_back(line);
      }
      return new_lines;
    }

    void CheckReprojection(std::vector<bool> &results,
                           const LineTrack &linetrack,
                           const ImageCollection &imagecols,
                           const double &th_angular2d, const double &th_perp2d)
    {
      results.clear();
      size_t num_supports = linetrack.count_lines();
      for (size_t i = 0; i < num_supports; ++i)
      {
        const Line2d &line2d = linetrack.line2d_list[i];
        Line2d line2d_projection = linetrack.line.projection(imagecols.camview(linetrack.image_id_list[i]));
        double angle = compute_angle<Line2d>(line2d, line2d_projection);
        if (angle > th_angular2d)
        {
          results.push_back(false);
          continue;
        }
        double dist_endperp = dist_endpoints_perpendicular_oneway<Line2d>(line2d, line2d_projection);
        if (dist_endperp > th_perp2d)
        {
          results.push_back(false);
          continue;
        }
        results.push_back(true);
      }
    }

    void CheckReprojection(const LineLinker2d &linker2d,
                           std::vector<bool> &results,
                           const LineTrack &linetrack,
                           const std::vector<CameraView> &views)
    {
      results.clear();
      size_t num_supports = linetrack.count_lines();

      for (size_t i = 0; i < num_supports; ++i)
      {
        const Line2d &line2d = linetrack.line2d_list[i];

        Line2d line2d_projection = linetrack.line.projection(views[i]);

        double score = linker2d.compute_reprojection_score(line2d, line2d_projection);
        if (score == 0)
        {
          results.push_back(false);
        }
        else
        {
          results.push_back(true);
        }
      }
    }

    void CheckReprojectionAndDepth(const LineLinker2d &linker2d,
                                   std::vector<bool> &results,
                                   const LineTrack &linetrack,
                                   const std::vector<CameraView> &views,
                                   const double th_depth,
                                   const double th_sensitivity)
    {
      results.clear();
      size_t num_supports = linetrack.count_lines();
      results.reserve(num_supports);

      const InfiniteLine3d inf_line3d(linetrack.line);
      const V3D dir = inf_line3d.direction();
      const V3D p_ref = inf_line3d.point();

      for (size_t i = 0; i < num_supports; ++i)
      {
        const auto &view = views[i];
        const Line2d &line2d = linetrack.line2d_list[i];

        Line2d line2d_projection = linetrack.line.projection(view);

        // check reprojection
        double score = linker2d.compute_reprojection_score(line2d, line2d_projection);
        if (score == 0)
        {
          results.push_back(false);
          continue;
        }

        // check positive depth
        InfiniteLine2d inf_line2d_projection(line2d_projection);
        V2D pstart_2d = inf_line2d_projection.point_projection(line2d.start);
        V2D pend_2d = inf_line2d_projection.point_projection(line2d.end);
        V3D pstart_3d = inf_line3d.unprojection(pstart_2d, view);
        double tstart = (pstart_3d - p_ref).dot(dir);
        pstart_3d = p_ref + dir * tstart;
        V3D pend_3d = inf_line3d.unprojection(pend_2d, view);
        double tend = (pend_3d - p_ref).dot(dir);
        pend_3d = p_ref + dir * tend;

        double z_start = view.pose.projdepth(pstart_3d);
        double z_end = view.pose.projdepth(pend_3d);
        if (z_start < th_depth || z_end < th_depth)
        {
          results.push_back(false);
          continue;
        }

        // check sensitivity
        double sensitivity = linetrack.line.sensitivity(view);
        if (sensitivity > th_sensitivity)
        {
          results.push_back(false);
          continue;
        }

        results.push_back(true);
      }
    }

    void FilterSupportingLines(std::vector<LineTrack> &new_linetracks,
                               const std::vector<LineTrack> &linetracks,
                               const ImageCollection &imagecols,
                               const double &th_angular2d, const double &th_perp2d,
                               const int num_outliers)
    {
      size_t num_tracks = linetracks.size();
      new_linetracks.clear();
      for (size_t track_id = 0; track_id < num_tracks; ++track_id)
      {
        const auto &track = linetracks[track_id];
        std::vector<bool> results;
        CheckReprojection(results, track, imagecols, th_angular2d, th_perp2d);
        size_t num_supports = track.count_lines();
        THROW_CHECK_EQ(num_supports, results.size());
        LineTrack newtrack;
        for (size_t support_id = 0; support_id < num_supports; ++support_id)
        {
          if (!results[support_id])
            continue;
          newtrack.node_id_list.push_back(track.node_id_list[support_id]);
          newtrack.image_id_list.push_back(track.image_id_list[support_id]);
          newtrack.line_id_list.push_back(track.line_id_list[support_id]);
          newtrack.line2d_list.push_back(track.line2d_list[support_id]);
          newtrack.line3d_list.push_back(track.line3d_list[support_id]);
          newtrack.score_list.push_back(track.score_list[support_id]);
        }
        if (newtrack.count_lines() == 0)
          continue;
        newtrack.line = Aggregator::aggregate_line3d_list(newtrack.line3d_list, newtrack.score_list, num_outliers);
        new_linetracks.push_back(newtrack);
      }
      // STDLOG(INFO) << "# tracks after filtering:" << " " << new_linetracks.size() << std::endl;
    }

    void CheckSensitivity(std::vector<bool> &results,
                          const LineTrack &linetrack,
                          const ImageCollection &imagecols,
                          const double &th_angular3d)
    {
      results.clear();
      size_t num_supports = linetrack.count_lines();
      for (size_t i = 0; i < num_supports; ++i)
      {
        // get line3d at the camera coordinate system
        const CameraView &view = imagecols.camview(linetrack.image_id_list[i]);

        // project and get camera ray corresponding to midpoint
        double sensitivity = linetrack.line.sensitivity(view);

        // compute angle
        if (sensitivity > th_angular3d)
          results.push_back(false);
        else
          results.push_back(true);
      }
    }

    void FilterTracksBySensitivity(std::vector<LineTrack> &new_linetracks,
                                   const std::vector<LineTrack> &linetracks,
                                   const ImageCollection &imagecols,
                                   const double &th_angular3d, const int &min_support_ns)
    {
      new_linetracks.clear();
      for (auto it = linetracks.begin(); it != linetracks.end(); ++it)
      {
        const auto &track = *it;
        size_t num_supports = track.count_lines();

        std::vector<bool> results;
        CheckSensitivity(results, track, imagecols, th_angular3d);
        std::set<int> support_images;
        for (size_t i = 0; i < num_supports; ++i)
        {
          if (results[i])
            support_images.insert(track.image_id_list[i]);
        }
        int counter = support_images.size();
        if (counter >= min_support_ns)
          new_linetracks.push_back(track);
      }
      std::cout << "filter line tracks by sensitivity (# after / # before): "
                << new_linetracks.size() << " / " << linetracks.size() << std::endl;
      // STDLOG(INFO) << "# tracks after filtering by sensitivity:" << " " << new_linetracks.size() << std::endl;
    }

    void FilterTracksByOverlap(std::vector<LineTrack> &new_linetracks,
                               const std::vector<LineTrack> &linetracks,
                               const ImageCollection &imagecols,
                               const double &th_overlap, const int &min_support_ns)
    {
      new_linetracks.clear();
      for (auto it = linetracks.begin(); it != linetracks.end(); ++it)
      {
        const auto &track = *it;
        size_t num_supports = track.count_lines();

        std::set<int> support_images;
        for (size_t i = 0; i < num_supports; ++i)
        {
          Line2d line2d_proj = track.line.projection(imagecols.camview(track.image_id_list[i]));
          const Line2d &line2d = track.line2d_list[i];
          double overlap = compute_overlap<Line2d>(line2d_proj, line2d);
          if (overlap >= th_overlap)
            support_images.insert(track.image_id_list[i]);
        }
        int counter = support_images.size();
        if (counter >= min_support_ns)
          new_linetracks.push_back(track);
      }
      std::cout << "filter line tracks by overlap (# after / # before): "
                << new_linetracks.size() << " / " << linetracks.size() << std::endl;
      // STDLOG(INFO) << "# tracks after filtering by overlap:" << " " << new_linetracks.size() << std::endl;
    }

  } // namespace merging

} // namespace limap
