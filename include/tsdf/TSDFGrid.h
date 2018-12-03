#pragma once
#include <Eigen/Dense>

#include <cmath>
#include <unordered_map>

namespace tsdf {

template <typename CellType>
class TSDFGrid {
  using HashContainer = std::unordered_map<
      int, std::unordered_map<int, std::unordered_map<int, CellType>>>;

 public:
  TSDFGrid() = default;

  TSDFGrid(const float cell_size, const int number_of_truncated_cells,
           const Eigen::Vector3f offset = Eigen::Vector3f::Zero())
      : cell_size_(cell_size),
        offset_(offset),
        inverse_cell_size_(1.f / cell_size),
        number_of_truncated_cells_(number_of_truncated_cells) {}

  /**
   * @brief at
   * @param[in] Eigen::Vector3i containing x, y, z indixes of the cell
   * @return a reference to to the cell if succedes
   *
   * @throw std::out_of_range if no such element exists in a grid
   */
  const CellType& at(const Eigen::Vector3i& index) const {
    if (hasKey(index)) {
      return operator()(index);
    } else {
      throw std::out_of_range("Cell is not found in a HashGrid");
    }
  }
  const CellType& at(const int x, const int y, const int z) const {
    return at(Eigen::Vector3i(x, y, z));
  }

  /**
   * @brief operator () either gives an access to an existing element or creates
   * a new one othwerise
   * @param index[in] Eigen::Vector3i containing x, y, z coordinates of the cell
   * @return a reference to a cell
   */
  CellType& operator()(const Eigen::Vector3i& index) {
    return grid_[index(0)][index(1)][index(2)];
  }
  const CellType& operator()(const Eigen::Vector3i& index) const {
    return grid_.at(index(0)).at(index(1)).at(index(2));
  }
  CellType& operator()(const int x, const int y, const int z) const {
    return operator()(Eigen::Vector3i(x, y, z));
  }

  /**
   * @brief hasKey
   * @param[in] Eigen::Vector3i containing x, y, z indixes of the cell
   * @return true if the cell is in the hash grid, otherwise false
   */

  bool hasKey(const Eigen::Vector3i& index) const {
    // check x coordinate
    const auto itr_coord_x = grid_.find(index(0));
    if (itr_coord_x == grid_.cend()) return false;
    // check y coordinate
    const auto itr_coord_y = itr_coord_x->second.find(index(1));
    if (itr_coord_y == itr_coord_x->second.cend()) return false;
    // check z coordinate
    const auto itr_coord_z = itr_coord_y->second.find(index(2));
    if (itr_coord_z == itr_coord_y->second.cend()) return false;
    return true;
  }

  float gridCellSize() const { return cell_size_; }

  Eigen::Vector3f gridToWorld(const Eigen::Vector3i& index) const {
    return index.cast<float>() * cell_size_ + offset_;
  }

  Eigen::Vector3i worldToGrid(const Eigen::Vector3f& point) const {
    Eigen::Vector3f point_in_grid = (point - offset_) * inverse_cell_size_;
    return point_in_grid.cast<int>();
  }

  void setOffset(const Eigen::Vector3f& offset) { offset_ = offset; }

  int getNumOfTruncCells() const { return number_of_truncated_cells_; }

  const HashContainer& getGrid() const { return grid_; }

 protected:
 private:
  float cell_size_ = 1.f;  // cell size in meters
  Eigen::Vector3f offset_ = Eigen::Vector3f::Zero();
  float inverse_cell_size_ = 1.f;
  int number_of_truncated_cells_ = 2;
  HashContainer grid_;
};
}  // namespace tsdf
