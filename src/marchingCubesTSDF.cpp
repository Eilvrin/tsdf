#include <tsdf/marchingCubesTSDF.h>
namespace tsdf {
void MarchingCubesTSDF::performReconstruction(pcl::PolygonMesh &output) {
  pcl::PointCloud<pcl::PointXYZ> cloud;
  pcl::PointCloud<pcl::PointXYZRGB> cloud_colored;

  for (const auto &grid_x : _dataTSDF.getGrid())
    for (const auto &grid_y : grid_x.second)
      for (const auto &grid_z : grid_y.second) {
        Eigen::Vector3i index_3d(grid_x.first, grid_y.first, grid_z.first);
        std::vector<float> leaf_node;
        if (!getValidNeighborList1D(leaf_node, index_3d)) continue;
        createSurface(leaf_node, index_3d, cloud);

        if (_integrateColor) {
          for (size_t i = cloud_colored.size(); i < cloud.size(); ++i) {
            pcl::PointXYZRGB pt_rgb;
            pt_rgb.getVector3fMap() = cloud.at(i).getVector3fMap();
            pt_rgb.r = _dataTSDF(index_3d).r;
            pt_rgb.g = _dataTSDF(index_3d).g;
            pt_rgb.b = _dataTSDF(index_3d).b;
            cloud_colored.push_back(pt_rgb);
          }
        }
      }
  if (_integrateColor)
    pcl::toPCLPointCloud2(cloud_colored, output.cloud);
  else
    pcl::toPCLPointCloud2(cloud, output.cloud);

  output.polygons.resize(cloud.size() / 3);
  for (size_t i = 0; i < output.polygons.size(); ++i) {
    pcl::Vertices v;
    v.vertices.resize(3);
    for (size_t j = 0; j < 3; ++j) v.vertices[j] = i * 3 + j;
    output.polygons[i] = v;
  }
}

bool MarchingCubesTSDF::getValidNeighborList1D(std::vector<float> &leaf,
                                               Eigen::Vector3i &index3d) {
  leaf = std::vector<float>(8, 0.0f);

  leaf[0] = getGridValue(index3d);
  if (pcl_isnan(leaf[0])) return (false);
  leaf[1] = getGridValue(index3d + Eigen::Vector3i(1, 0, 0));
  if (pcl_isnan(leaf[1])) return (false);
  leaf[2] = getGridValue(index3d + Eigen::Vector3i(1, 0, 1));
  if (pcl_isnan(leaf[2])) return (false);
  leaf[3] = getGridValue(index3d + Eigen::Vector3i(0, 0, 1));
  if (pcl_isnan(leaf[3])) return (false);
  leaf[4] = getGridValue(index3d + Eigen::Vector3i(0, 1, 0));
  if (pcl_isnan(leaf[4])) return (false);
  leaf[5] = getGridValue(index3d + Eigen::Vector3i(1, 1, 0));
  if (pcl_isnan(leaf[5])) return (false);
  leaf[6] = getGridValue(index3d + Eigen::Vector3i(1, 1, 1));
  if (pcl_isnan(leaf[6])) return (false);
  leaf[7] = getGridValue(index3d + Eigen::Vector3i(0, 1, 1));
  if (pcl_isnan(leaf[7])) return (false);
  return (true);
}

float MarchingCubesTSDF::getGridValue(Eigen::Vector3i pos) {
  if (_dataTSDF.hasKey(pos) && _dataTSDF(pos).weights > _minWeight) {
    return (_dataTSDF(pos).distances);
  } else
    return (std::numeric_limits<float>::quiet_NaN());
}

void MarchingCubesTSDF::createSurface(std::vector<float> &leaf_node,
                                      Eigen::Vector3i &index_3d,
                                      pcl::PointCloud<pcl::PointXYZ> &cloud) {
  int cubeindex = 0;
  Eigen::Vector3f vertex_list[12];
  if (leaf_node[0] < iso_level_) cubeindex |= 1;
  if (leaf_node[1] < iso_level_) cubeindex |= 2;
  if (leaf_node[2] < iso_level_) cubeindex |= 4;
  if (leaf_node[3] < iso_level_) cubeindex |= 8;
  if (leaf_node[4] < iso_level_) cubeindex |= 16;
  if (leaf_node[5] < iso_level_) cubeindex |= 32;
  if (leaf_node[6] < iso_level_) cubeindex |= 64;
  if (leaf_node[7] < iso_level_) cubeindex |= 128;

  if (pcl::edgeTable[cubeindex] == 0) return;

  std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > p;
  p.resize(8);
  p[0] = _dataTSDF.gridToWorld(index_3d);
  p[1] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(1, 0, 0));
  p[2] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(1, 0, 1));
  p[3] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(0, 0, 1));
  p[4] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(0, 1, 0));
  p[5] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(1, 1, 0));
  p[6] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(1, 1, 1));
  p[7] = _dataTSDF.gridToWorld(index_3d + Eigen::Vector3i(0, 1, 1));

  if (pcl::edgeTable[cubeindex] & 1)
    interpolateEdge(p[0], p[1], leaf_node[0], leaf_node[1], vertex_list[0]);
  if (pcl::edgeTable[cubeindex] & 2)
    interpolateEdge(p[1], p[2], leaf_node[1], leaf_node[2], vertex_list[1]);
  if (pcl::edgeTable[cubeindex] & 4)
    interpolateEdge(p[2], p[3], leaf_node[2], leaf_node[3], vertex_list[2]);
  if (pcl::edgeTable[cubeindex] & 8)
    interpolateEdge(p[3], p[0], leaf_node[3], leaf_node[0], vertex_list[3]);
  if (pcl::edgeTable[cubeindex] & 16)
    interpolateEdge(p[4], p[5], leaf_node[4], leaf_node[5], vertex_list[4]);
  if (pcl::edgeTable[cubeindex] & 32)
    interpolateEdge(p[5], p[6], leaf_node[5], leaf_node[6], vertex_list[5]);
  if (pcl::edgeTable[cubeindex] & 64)
    interpolateEdge(p[6], p[7], leaf_node[6], leaf_node[7], vertex_list[6]);
  if (pcl::edgeTable[cubeindex] & 128)
    interpolateEdge(p[7], p[4], leaf_node[7], leaf_node[4], vertex_list[7]);
  if (pcl::edgeTable[cubeindex] & 256)
    interpolateEdge(p[0], p[4], leaf_node[0], leaf_node[4], vertex_list[8]);
  if (pcl::edgeTable[cubeindex] & 512)
    interpolateEdge(p[1], p[5], leaf_node[1], leaf_node[5], vertex_list[9]);
  if (pcl::edgeTable[cubeindex] & 1024)
    interpolateEdge(p[2], p[6], leaf_node[2], leaf_node[6], vertex_list[10]);
  if (pcl::edgeTable[cubeindex] & 2048)
    interpolateEdge(p[3], p[7], leaf_node[3], leaf_node[7], vertex_list[11]);

  // Create the triangle
  for (int i = 0; pcl::triTable[cubeindex][i] != -1; i += 3) {
    pcl::PointXYZ p1, p2, p3;
    p1.x = vertex_list[pcl::triTable[cubeindex][i]][0];
    p1.y = vertex_list[pcl::triTable[cubeindex][i]][1];
    p1.z = vertex_list[pcl::triTable[cubeindex][i]][2];
    cloud.push_back(p1);
    p2.x = vertex_list[pcl::triTable[cubeindex][i + 1]][0];
    p2.y = vertex_list[pcl::triTable[cubeindex][i + 1]][1];
    p2.z = vertex_list[pcl::triTable[cubeindex][i + 1]][2];
    cloud.push_back(p2);
    p3.x = vertex_list[pcl::triTable[cubeindex][i + 2]][0];
    p3.y = vertex_list[pcl::triTable[cubeindex][i + 2]][1];
    p3.z = vertex_list[pcl::triTable[cubeindex][i + 2]][2];
    cloud.push_back(p3);
  }
}

void MarchingCubesTSDF::voxelizeData() {
  // Do nothing
}
}  // namespace tsdf
