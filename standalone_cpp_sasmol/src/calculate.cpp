#include "sasmol/calculate.hpp"

#include "sasmol/properties.hpp"

#include <limits>
#include <stdexcept>

namespace sasmol {

namespace {

std::vector<std::size_t> all_frame_indices(const Molecule& molecule) {
  std::vector<std::size_t> frames;
  frames.reserve(molecule.number_of_frames());
  for (std::size_t frame = 0; frame < molecule.number_of_frames(); ++frame) {
    frames.push_back(frame);
  }
  return frames;
}

}  // namespace

MassCalculationResult calculate_mass(Molecule& molecule) {
  if (molecule.element().size() != molecule.natoms()) {
    throw std::invalid_argument("calculate_mass requires one element per atom");
  }

  const auto& weights = amu();
  molecule.mass().assign(molecule.natoms(), calc_type{});
  MassCalculationResult result;

  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto found = weights.find(molecule.element()[atom]);
    if (found == weights.end()) {
      result.unknown_elements.push_back(molecule.element()[atom]);
      continue;
    }
    molecule.mass()[atom] = found->second;
    result.total_mass += found->second;
  }

  molecule.set_total_mass(result.total_mass);
  return result;
}

CoordinateBounds calculate_minimum_and_maximum(
    const Molecule& molecule, const std::vector<std::size_t>& frames) {
  if (molecule.natoms() == 0 || molecule.number_of_frames() == 0) {
    throw std::invalid_argument(
        "calculate_minimum_and_maximum requires atoms and frames");
  }

  const auto selected_frames = frames.empty() ? all_frame_indices(molecule) : frames;
  if (selected_frames.empty()) {
    throw std::invalid_argument("calculate_minimum_and_maximum requires frames");
  }

  CoordinateBounds bounds{
      {std::numeric_limits<coord_type>::max(),
       std::numeric_limits<coord_type>::max(),
       std::numeric_limits<coord_type>::max()},
      {std::numeric_limits<coord_type>::lowest(),
       std::numeric_limits<coord_type>::lowest(),
       std::numeric_limits<coord_type>::lowest()}};

  for (const auto frame : selected_frames) {
    if (frame >= molecule.number_of_frames()) {
      throw std::out_of_range("calculate_minimum_and_maximum frame is out of range");
    }
    for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
      const auto xyz = molecule.coordinate(frame, atom);
      if (xyz.x < bounds.minimum.x) bounds.minimum.x = xyz.x;
      if (xyz.y < bounds.minimum.y) bounds.minimum.y = xyz.y;
      if (xyz.z < bounds.minimum.z) bounds.minimum.z = xyz.z;
      if (xyz.x > bounds.maximum.x) bounds.maximum.x = xyz.x;
      if (xyz.y > bounds.maximum.y) bounds.maximum.y = xyz.y;
      if (xyz.z > bounds.maximum.z) bounds.maximum.z = xyz.z;
    }
  }

  return bounds;
}

CoordinateBounds calculate_minimum_and_maximum_all_steps(
    const Molecule& molecule) {
  return calculate_minimum_and_maximum(molecule);
}

CoordinateBounds calc_minmax_all_steps(const Molecule& molecule) {
  return calculate_minimum_and_maximum_all_steps(molecule);
}

}  // namespace sasmol
