#include <tsdf/TSDF.h>

namespace tsdf {

Eigen::Vector2i TSDF::getProjectedPointIndex(const Eigen::Vector3f& point,
                                             const int width,
                                             const int height) const {
  Eigen::Vector2i index;
  index(0) = floor((point.x() / point.z()) * kinect_param_.focalLength_x +
                   kinect_param_.principalPointX);
  index(1) = floor((point.y() / point.z()) * kinect_param_.focalLength_y +
                   kinect_param_.principalPointY);
  index(0) = std::min(width, index(0));
  index(0) = std::max(0, index(0));
  index(1) = std::min(height, index(1));
  index(1) = std::max(0, index(1));
  return index;
}

Eigen::Vector3f TSDF::transformPoint(const Eigen::Affine3f& transform,
                                     const Eigen::Vector3f& point) const {
  return transform * point;
}

float TSDF::calculateDistance(const Eigen::Vector3f& point_1,
                              const Eigen::Vector3f& point_2) const {
  Eigen::Vector3f diff = point_2 - point_1;
  return sqrt(diff(0) * diff(0) + diff(1) * diff(1) + diff(2) * diff(2));
}

Eigen::Vector3f TSDF::cellCenterCartesian(const Eigen::Vector3f& point) const {
  Eigen::Vector3f cell_center_cartesian;
  cell_center_cartesian << point(0) + options_.map_resolution / 2,
      point(1) + options_.map_resolution / 2,
      point(2) + options_.map_resolution / 2;
  return cell_center_cartesian;
}

void TSDF::calcTableDeltaUV(
    const float& step_size,
    std::vector<std::tuple<float, int, int> >& tableUV) {
  float max_range = kinect_param_.max_range + options_.truncation_dist;
  size_t max_range_idx = ceil(max_range / step_size);
  tableUV.resize(max_range_idx + 1);

  float half_cell = options_.map_resolution / 2;
  for (size_t i = 0; i < max_range_idx; i++) {
    float range = (i + 1) * step_size;
    int du = floor((half_cell / range) * kinect_param_.focalLength_x +
                   kinect_param_.principalPointX) -
             floor((-half_cell / range) * kinect_param_.focalLength_x +
                   kinect_param_.principalPointX);

    int dv = floor((half_cell / range) * kinect_param_.focalLength_y +
                   kinect_param_.principalPointY) -
             floor((-half_cell / range) * kinect_param_.focalLength_y +
                   kinect_param_.principalPointY);

    tableUV[i] = std::make_tuple(range, du, dv);
  }
}

bool TSDF::inUV(const Eigen::Vector2i& point_to_compare,
                const Eigen::Vector2i& point_center, const int& du,
                const int& dv) const {
  return ((point_to_compare(0) >= (point_center(0) - floor(du / 2))) &&
          (point_to_compare(1) >= (point_center(1) - floor(dv / 2))) &&
          (point_to_compare(0) <= (point_center(0) + floor(du / 2))) &&
          (point_to_compare(1) <= (point_center(1) + floor(dv / 2))));
}

void TSDF::writeCSV(const std::string& filename) {
  std::ofstream file(filename);
  file << "x, y, z, distances, weights\n";
  Eigen::Vector3i index;
  for (const auto& grid_x : tsdfMap_.getGrid()) {
    for (const auto& grid_y : grid_x.second) {
      for (const auto& grid_z : grid_y.second) {
        auto r = grid_x.first;
        auto c = grid_y.first;
        auto l = grid_z.first;
        index << r, c, l;
        file << r << "," << c << "," << l << "," << tsdfMap_.at(index).distances
             << "," << tsdfMap_.at(index).weights << "\n";
      }
    }
  }
}

template <typename T>
void TSDF::addPointCloud(const typename pcl::PointCloud<T>::Ptr& cloud) {
  Eigen::Vector3f camera_position = cloud->sensor_origin_.template head<3>();
  const Eigen::Affine3f global_T_local =
      Eigen::Translation3f(camera_position) * cloud->sensor_orientation_;
  const Eigen::Affine3f local_T_global = global_T_local.inverse();

  for (size_t point_idx = 0; point_idx < cloud->size(); point_idx++) {
    if (isnan(cloud->points[point_idx].x) ||
        isnan(cloud->points[point_idx].y) || isnan(cloud->points[point_idx].z))
      continue;

    // Get the central point
    Eigen::Vector3f point;
    point << cloud->points[point_idx].x, cloud->points[point_idx].y,
        cloud->points[point_idx].z;

    Eigen::Vector2i projected_point;
    projected_point =
        getProjectedPointIndex(point, cloud->width, cloud->height);
    // Distance from central point to camera
    float dist_camera_point =
        calculateDistance(point, Eigen::Vector3f(0, 0, 0));

    // Transform central point to global coordinates
    Eigen::Vector3f gl_point;
    gl_point = transformPoint(global_T_local, point);

    // Find voxel in global coordinates
    Eigen::Vector3i gl_point_cell = tsdfMap_.worldToGrid(gl_point);

    // Check neighborhood (in global coordinates)
    for (int zi = gl_point_cell(2) - tsdfMap_.getNumOfTruncCells();
         zi <= gl_point_cell(2) + tsdfMap_.getNumOfTruncCells(); zi++)
      for (int yi = gl_point_cell(1) - tsdfMap_.getNumOfTruncCells();
           yi <= gl_point_cell(1) + tsdfMap_.getNumOfTruncCells(); yi++)
        for (int xi = gl_point_cell(0) - tsdfMap_.getNumOfTruncCells();
             xi <= gl_point_cell(0) + tsdfMap_.getNumOfTruncCells(); xi++) {
          // Take a voxel in global xi,yi,zi and transform it to the world
          // from grid
          Eigen::Vector3i cell_idx;
          cell_idx << xi, yi, zi;
          Eigen::Vector3f gl_neigh_cell = tsdfMap_.gridToWorld(cell_idx);
          // Find center of the voxel in global coordinates
          gl_neigh_cell = cellCenterCartesian(gl_neigh_cell);
          // Transform center of the voxel to local coordinates
          Eigen::Vector3f l_neigh_cell =
              transformPoint(local_T_global, gl_neigh_cell);

          // Project local coordinate onto the projection plane
          Eigen::Vector2i l_projected_neigh;
          l_projected_neigh =
              getProjectedPointIndex(l_neigh_cell, cloud->width, cloud->height);

          int du =
              std::get<1>(tableUV[floor(l_neigh_cell(2) / step_size_range_)]);
          int dv =
              std::get<2>(tableUV[floor(l_neigh_cell(2) / step_size_range_)]);

          // Confirm that this reprojected cell is on the ray
          if (inUV(l_projected_neigh, projected_point, du, dv)) {
            // Distance from neighborhood point to camera
            float dist_camera_neigh =
                calculateDistance(l_neigh_cell, Eigen::Vector3f(0, 0, 0));
            // Distance from central point to neighborhood point
            float dist = calculateDistance(l_neigh_cell, point);
            int sign = 1;
            if (dist_camera_neigh > dist_camera_point) sign = -1;

            float w_old = tsdfMap_(cell_idx).weights;
            float w_sum = w_old + weight_;
            if (w_sum == 0) continue;

            if (options_.integrate_color) {
              integrateColor<T>(cell_idx, w_old, w_sum, point_idx, cloud);
            }

            float d_old = tsdfMap_(cell_idx).distances;

            tsdfMap_(cell_idx).distances =
                (d_old * w_old + sign * weight_ * dist) / w_sum;
            tsdfMap_(cell_idx).weights += weight_;
          }
        }
  }
}
template <>
void TSDF::integrateColor<pcl::PointXYZRGB>(
    const Eigen::Vector3i& cell, const float w_old, const float w_sum,
    const size_t point_idx,
    const pcl::PointCloud<pcl::PointXYZRGB>::Ptr& cloud) {
  uint8_t r_old = tsdfMap_.at(cell).r;
  uint8_t g_old = tsdfMap_.at(cell).g;
  uint8_t b_old = tsdfMap_.at(cell).b;
  tsdfMap_(cell).r = static_cast<uint8_t>(
      (w_old * r_old + weight_ * cloud->points[point_idx].r) / w_sum);
  tsdfMap_(cell).g = static_cast<uint8_t>(
      (w_old * g_old + weight_ * cloud->points[point_idx].g) / w_sum);
  tsdfMap_(cell).b = static_cast<uint8_t>(
      (w_old * b_old + weight_ * cloud->points[point_idx].b) / w_sum);
}

template void TSDF::addPointCloud<pcl::PointXYZRGB>(
    const typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr&);
template void TSDF::addPointCloud<pcl::PointXYZ>(
    const typename pcl::PointCloud<pcl::PointXYZ>::Ptr&);

}  // namespace tsdf
