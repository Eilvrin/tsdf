#pragma once
#include <tsdf/TSDFGrid.h>

#include <pcl/common/common.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <fstream>
#include <limits>
#include <tuple>
#include <vector>

namespace tsdf {

struct TSDFCell {
  TSDFCell() : distances(0), weights(0), r(0), g(0), b(0) {}
  float distances;
  float weights;
  uint8_t r;
  uint8_t g;
  uint8_t b;
  uint8_t dummy;  // for better memory alignment
};

class TSDF {
  struct colorRGB {
    double r, g, b;
  };

 public:
  struct Options {
    float map_resolution = 0.01f;
    float truncation_dist = 0.02f;
    float min_weight = 0.0;
    bool integrate_color = false;
  };
  struct KinectParameters {
    float focalLength_x = 517.3f;
    float focalLength_y = 516.5f;
    float principalPointX = 318.6f;
    float principalPointY = 255.3f;
    float max_range = 5.0f;
  };

  TSDF(const Options& options)
      : options_(options),
        tsdfMap_(options.map_resolution,
                 std::ceil(options.truncation_dist / options.map_resolution)) {
    calcTableDeltaUV(step_size_range_, tableUV);
  }
  template <typename T>
  void addPointCloud(const typename pcl::PointCloud<T>::Ptr& cloud);

  tsdf::TSDFGrid<TSDFCell>& getGridMap() { return tsdfMap_; }
  void writeCSV(const std::string& filename);
  void setMapOffset(const Eigen::Vector3f& offset) {
    tsdfMap_.setOffset(offset);
  }

 protected:
  void calcTableDeltaUV(const float& step_size,
                        std::vector<std::tuple<float, int, int> >& tableUV);

  bool inUV(const Eigen::Vector2i& point_to_compare,
            const Eigen::Vector2i& point_center, const int& du,
            const int& dv) const;
  Eigen::Vector2i getProjectedPointIndex(const Eigen::Vector3f& point,
                                         const int width,
                                         const int height) const;
  Eigen::Vector3f cellCenterCartesian(const Eigen::Vector3f& point) const;
  Eigen::Vector3f transformPoint(const Eigen::Affine3f& transform,
                                 const Eigen::Vector3f& point) const;
  float calculateDistance(const Eigen::Vector3f& point_1,
                          const Eigen::Vector3f& point_2) const;

  template <typename U>
  void integrateColor(const Eigen::Vector3i& /*cell*/, const float /*w_old*/,
                      const float /*w_sum*/, const size_t /*point_idx*/,
                      const typename pcl::PointCloud<U>::Ptr& /*cloud*/) {}

 private:
  KinectParameters kinect_param_;
  Options options_;
  tsdf::TSDFGrid<TSDFCell> tsdfMap_;
  std::vector<std::tuple<float, int, int> > tableUV;
  const float step_size_range_ = 0.01f;
  const float weight_ = 1.0f;
  float normalizeAngle(float angle) {
    static constexpr double kTwoPI = 2.0 * M_PI;
    angle = std::fmod(angle, kTwoPI);
    if (angle < 0.0f) {
      angle += kTwoPI;
    }
    return angle;
  }
};

extern template void TSDF::addPointCloud<pcl::PointXYZRGB>(
    const typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr&);
extern template void TSDF::addPointCloud<pcl::PointXYZ>(
    const typename pcl::PointCloud<pcl::PointXYZ>::Ptr&);
}  // namespace tsdf
