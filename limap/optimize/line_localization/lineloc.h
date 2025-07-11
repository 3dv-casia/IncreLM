#ifndef LIMAP_OPTIMIZE_LINE_LOCALIZATION_LINELOC_H_
#define LIMAP_OPTIMIZE_LINE_LOCALIZATION_LINELOC_H_

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "_limap/helpers.h"

#include "util/types.h"

#include <ceres/ceres.h>
#include "base/linebase.h"
#include "optimize/line_localization/lineloc_config.h"

namespace py = pybind11;

namespace limap
{

  namespace optimize
  {

    namespace line_localization
    {

      class LineLocEngine
      {
      protected:
        LineLocConfig config_;
        std::vector<Line3d> l3ds;
        std::vector<std::vector<Line2d>> l2ds; // l3ds[i] <---> l2ds[i] (line2d track)
        V4D cam_kvec;
        Camera cam;
        CameraPose campose;

        // set up ceres problem
        void ParameterizeCamera();
        virtual void AddResiduals();

      public:
        LineLocEngine() {}
        LineLocEngine(const LineLocConfig &cfg) : config_(cfg)
        {
          // If the error of the line is measured in 3D, the error of the point is also measured in 3D
          if (config_.cost_function == E3DLineLineDist2 || config_.cost_function == E3DPlaneLineDist2)
            config_.points_3d_dist = true;
        }

        void Initialize(const std::vector<Line3d> &l3ds, const std::vector<std::vector<Line2d>> &l2ds, M3D K, M3D R, V3D T)
        {
          CHECK_EQ(l3ds.size(), l2ds.size());
          this->l3ds = l3ds;
          this->l2ds = l2ds;
          this->cam = Camera(colmap::PinholeCameraModel::model_id, K);
          this->campose = CameraPose(R, T);
        }
        void Initialize(const std::vector<Line3d> &l3ds, const std::vector<std::vector<Line2d>> &l2ds, V4D kvec, V4D qvec, V3D tvec)
        {
          CHECK_EQ(l3ds.size(), l2ds.size());
          this->l3ds = l3ds;
          this->l2ds = l2ds;
          M3D K = M3D::Zero();
          K(0, 0) = kvec(0);
          K(1, 1) = kvec(1);
          K(0, 2) = kvec(2);
          K(1, 2) = kvec(3);
          this->cam = Camera(colmap::PinholeCameraModel::model_id, K);
          this->campose = CameraPose(qvec, tvec);
        }
        void SetUp();
        bool Solve();

        // output
        M3D GetFinalR() const { return campose.R(); }
        V4D GetFinalQ() const { return campose.qvec; }
        V3D GetFinalT() const { return campose.T(); }
        double GetInitialCost() const { return std::sqrt(summary_.initial_cost / summary_.num_residuals_reduced); }
        double GetFinalCost() const { return std::sqrt(summary_.final_cost / summary_.num_residuals_reduced); }
        bool IsSolutionUsable() const { return summary_.IsSolutionUsable(); }

        // ceres
        std::unique_ptr<ceres::Problem> problem_;
        ceres::Solver::Summary summary_;
      };

      // Jointly optimize poses using points' and lines' 2D-3D correspondences, respectively.
      class JointLocEngine : public LineLocEngine
      {
      protected:
        std::vector<V3D> p3ds;
        std::vector<V2D> p2ds;
        void AddResiduals() override;

      public:
        JointLocEngine() {}
        JointLocEngine(const LineLocConfig &cfg) : LineLocEngine(cfg) {}

        void Initialize(const std::vector<Line3d> &l3ds, const std::vector<std::vector<Line2d>> &l2ds,
                        const std::vector<V3D> &p3ds, const std::vector<V2D> &p2ds, M3D K, M3D R, V3D T)
        {
          CHECK_EQ(p3ds.size(), p2ds.size());
          CHECK_EQ(l3ds.size(), l2ds.size());
          this->p3ds = p3ds;
          this->p2ds = p2ds;
          this->l3ds = l3ds;
          this->l2ds = l2ds;
          this->cam = Camera(colmap::PinholeCameraModel::model_id, K);
          this->campose = CameraPose(R, T);
        }
        void Initialize(const std::vector<Line3d> &l3ds, const std::vector<std::vector<Line2d>> &l2ds,
                        const std::vector<V3D> &p3ds, const std::vector<V2D> &p2ds, V4D kvec, V4D qvec, V3D tvec)
        {
          CHECK_EQ(p3ds.size(), p2ds.size());
          CHECK_EQ(l3ds.size(), l2ds.size());
          this->p3ds = p3ds;
          this->p2ds = p2ds;
          this->l3ds = l3ds;
          this->l2ds = l2ds;
          M3D K = M3D::Zero();
          K(0, 0) = kvec(0);
          K(1, 1) = kvec(1);
          K(0, 2) = kvec(2);
          K(1, 2) = kvec(3);
          this->cam = Camera(colmap::PinholeCameraModel::model_id, K);
          this->campose = CameraPose(qvec, tvec);
        }
      };

    } // namespace line_localization

  } // namespace optimize

} // namespace limap

#endif
