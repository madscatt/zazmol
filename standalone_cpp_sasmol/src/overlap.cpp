#include "sasmol/overlap.hpp"

#include <stdexcept>

namespace sasmol {

bool has_overlap(const std::vector<Vec3>& first, const std::vector<Vec3>& second,
                 calc_type cutoff) {
  if (cutoff < 0.0) {
    throw std::invalid_argument("has_overlap requires a non-negative cutoff");
  }
  const auto cutoff_squared = cutoff * cutoff;
  for (const auto& a : first) {
    for (const auto& b : second) {
      const calc_type dx = static_cast<calc_type>(a.x) - b.x;
      const calc_type dy = static_cast<calc_type>(a.y) - b.y;
      const calc_type dz = static_cast<calc_type>(a.z) - b.z;
      if (dx * dx + dy * dy + dz * dz < cutoff_squared) {
        return true;
      }
    }
  }
  return false;
}

bool has_overlap(const Molecule& first, std::size_t first_frame,
                 const Molecule& second, std::size_t second_frame,
                 calc_type cutoff) {
  std::vector<Vec3> first_coordinates;
  std::vector<Vec3> second_coordinates;
  first_coordinates.reserve(first.natoms());
  second_coordinates.reserve(second.natoms());
  for (std::size_t atom = 0; atom < first.natoms(); ++atom) {
    first_coordinates.push_back(first.coordinate(first_frame, atom));
  }
  for (std::size_t atom = 0; atom < second.natoms(); ++atom) {
    second_coordinates.push_back(second.coordinate(second_frame, atom));
  }
  return has_overlap(first_coordinates, second_coordinates, cutoff);
}

}  // namespace sasmol
