#include <gtest/gtest.h>
#include <tsdf/TSDF.h>
#include <tsdf/TSDFGrid.h>

TEST(TestGridMap, TestWorldToGrid) {
  tsdf::TSDFGrid<tsdf::TSDFCell> grid_map(1, 2, Eigen::Vector3f(1.f, 1.f, 2.f));

  const Eigen::Vector3f point(3.5f, 2.1f, 4.9f);

  const Eigen::Vector3i grid_idx = grid_map.worldToGrid(point);

  EXPECT_TRUE(grid_idx == Eigen::Vector3i(2, 1, 2));
}
