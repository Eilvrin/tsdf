#include <pcl/common/common.h>
#include <pcl/common/time.h>
#include <pcl/console/parse.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

#include <tsdf/TSDF.h>
#include <tsdf/marchingCubesTSDF.h>

#include <Eigen/Core>
#include <iostream>
#include <string>
#include <vector>

// clang-format off
void printUsage(){
  std::cerr << "Usage: tsdf [options] " << std::endl;
  std::cerr << "-p:                 point cloud prefix" << std::endl;
  std::cerr << "-m  <int>:          maximum number of scans that should be processed" << std::endl;
  std::cerr << "-e  <int>:          process only every e-th scan" << std::endl;
  std::cerr << "-r  <float>:        resolution of TSDF map in meters (default 0.01)"  << std::endl;
  std::cerr << "-td <float>:        trunctation distance in meters (default 0.02)"  << std::endl;
  std::cerr << "-minw <float>:      minimum weight for marching cubes (default 0.0)"  << std::endl;
  std::cerr << "-o <output.ply>:    output mesh file name (default mesh.ply)"  << std::endl;
  std::cerr << "-c:                 integrate color (default false)"  << std::endl;
  std::cerr << "-csv :              write csv file (default false)"  << std::endl;
  std::cerr << "-ocsv <output.csv>: output csv file name (default tsdf.csv)"  << std::endl;
  std::cerr << "-h:                 this help" << std::endl;
}
// clang-format on

void parseTSDFOptions(int argc, char** argv, tsdf::TSDF::Options& options) {
  pcl::console::parse(argc, argv, "-r", options.map_resolution);
  pcl::console::parse(argc, argv, "-td", options.truncation_dist);
  pcl::console::parse(argc, argv, "-minw", options.min_weight);

  if (pcl::console::find_switch(argc, argv, "-c"))
    options.integrate_color = true;
}

int main(int argc, char** argv) {
  if (pcl::console::find_argument(argc, argv, "-h") >= 0) {
    printUsage();
    return EXIT_SUCCESS;
  }

  std::string pcPrefix = "ncloud";
  pcl::console::parse(argc, argv, "-p", pcPrefix);

  int maxStep = 1;
  pcl::console::parse(argc, argv, "-m", maxStep);

  int skipScansPerStep = 1;
  pcl::console::parse(argc, argv, "-e", skipScansPerStep);

  std::string outmeshFileName = "mesh.ply";
  pcl::console::parse(argc, argv, "-o", outmeshFileName);

  bool write_csv = false;
  if (pcl::console::find_switch(argc, argv, "-csv")) write_csv = true;

  std::string outFileName = "tsdf.csv";
  pcl::console::parse(argc, argv, "-ocsv", outFileName);

  tsdf::TSDF::Options tsdf_options;
  parseTSDFOptions(argc, argv, tsdf_options);

  tsdf::TSDF tsdf_algorithm(tsdf_options);

  double timing = pcl::getTime();

  for (int i = 0; i < maxStep; i = i + skipScansPerStep) {
    char cloudFileName[2048];
    sprintf(cloudFileName, "%s%05d.pcd", pcPrefix.c_str(), i);

    if (tsdf_options.integrate_color) {
      pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_colored(
          new pcl::PointCloud<pcl::PointXYZRGB>);
      pcl::io::loadPCDFile(cloudFileName, *cloud_colored);
      if (cloud_colored->size() == 0) continue;
      tsdf_algorithm.addPointCloud<pcl::PointXYZRGB>(cloud_colored);
    } else {
      pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(
          new pcl::PointCloud<pcl::PointXYZ>);
      pcl::io::loadPCDFile(cloudFileName, *cloud);
      if (cloud->size() == 0) continue;
      tsdf_algorithm.addPointCloud<pcl::PointXYZ>(cloud);
    }

    std::cout << "." << std::flush;
  }
  std::cerr << "TSDF took: " << (pcl::getTime() - timing) << "s." << std::endl;
  double time_mar = pcl::getTime();

  tsdf::MarchingCubesTSDF mcubes;
  mcubes.setTSDFData(tsdf_algorithm.getGridMap());
  pcl::PolygonMesh::Ptr mesh(new pcl::PolygonMesh);
  mcubes.setMinWeight(tsdf_options.min_weight);
  mcubes.setIntegrateColor(tsdf_options.integrate_color);
  mcubes.performReconstruction(*mesh);
  pcl::io::savePLYFile(outmeshFileName, *mesh);
  if (write_csv) tsdf_algorithm.writeCSV(outFileName);

  std::cerr << "Marching cubes took: " << (pcl::getTime() - time_mar) << "s."
            << std::endl;

  std::cout << std::endl;
  return 0;
}
