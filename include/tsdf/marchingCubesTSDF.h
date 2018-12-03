#pragma once
#include <tsdf/TSDF.h>
#include <tsdf/TSDFGrid.h>

#include <pcl/common/common.h>
#include <pcl/surface/marching_cubes.h>
#include <pcl/surface/impl/marching_cubes.hpp>

#include <Eigen/Core>

namespace tsdf {
class MarchingCubesTSDF : public pcl::MarchingCubes<pcl::PointXYZ> {
 public:
  MarchingCubesTSDF() : pcl::MarchingCubes<pcl::PointXYZ>() {}

  void setTSDFData(tsdf::TSDFGrid<TSDFCell> &tsdfHashMap) {
    _dataTSDF = tsdfHashMap;
    setPercentageExtendGrid(0);
    setIsoLevel(0);
  }

  void performReconstruction(pcl::PolygonMesh &output);

  void setIntegrateColor(bool integrateColor) {
    _integrateColor = integrateColor;
  }

  void setMinWeight(float minWeight) { _minWeight = minWeight; }

 protected:
  void voxelizeData();
  bool getValidNeighborList1D(std::vector<float> &leaf,
                              Eigen::Vector3i &index3d);
  float getGridValue(Eigen::Vector3i pos);
  void createSurface(std::vector<float> &leaf_node, Eigen::Vector3i &index_3d,
                     pcl::PointCloud<pcl::PointXYZ> &cloud);

  tsdf::TSDFGrid<TSDFCell> _dataTSDF;
  bool _integrateColor = false;
  float _minWeight = 0.0f;
};
}  // namespace tsdf
