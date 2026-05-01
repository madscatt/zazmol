#include "sasmol/calculate.hpp"

#include "sasmol/properties.hpp"

#include <cmath>
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

std::map<std::string, std::size_t> calculate_molecular_formula(
    Molecule& molecule) {
  if (molecule.element().size() != molecule.natoms()) {
    throw std::invalid_argument(
        "calculate_molecular_formula requires one element per atom");
  }

  std::map<std::string, std::size_t> formula;
  for (const auto& element : molecule.element()) {
    ++formula[element];
  }

  molecule.formula() = formula;
  return formula;
}

void calculate_residue_charge(Molecule& molecule) {
  if (molecule.resid().size() != molecule.natoms() ||
      molecule.atom_charge().size() != molecule.natoms()) {
    throw std::invalid_argument(
        "calculate_residue_charge requires resid and atom_charge per atom");
  }

  std::map<int, calc_type> charge_by_resid;
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    charge_by_resid[molecule.resid()[atom]] += molecule.atom_charge()[atom];
  }

  molecule.residue_charge().assign(molecule.natoms(), calc_type{});
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    molecule.residue_charge()[atom] = charge_by_resid[molecule.resid()[atom]];
  }
}

CalcVec3 calculate_center_of_mass(Molecule& molecule, std::size_t frame) {
  if (frame >= molecule.number_of_frames()) {
    throw std::out_of_range("calculate_center_of_mass frame is out of range");
  }
  if (molecule.natoms() == 0) {
    throw std::invalid_argument("calculate_center_of_mass requires atoms");
  }
  if (molecule.mass().size() != molecule.natoms() ||
      molecule.total_mass() <= calc_type{}) {
    const auto mass_result = calculate_mass(molecule);
    if (!mass_result.ok()) {
      throw std::invalid_argument(
          "calculate_center_of_mass requires known element masses");
    }
  }
  if (molecule.total_mass() <= calc_type{}) {
    throw std::invalid_argument("calculate_center_of_mass requires positive mass");
  }

  CalcVec3 center;
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto xyz = molecule.coordinate(frame, atom);
    const auto mass = molecule.mass()[atom];
    center.x += mass * static_cast<calc_type>(xyz.x);
    center.y += mass * static_cast<calc_type>(xyz.y);
    center.z += mass * static_cast<calc_type>(xyz.z);
  }

  center.x /= molecule.total_mass();
  center.y /= molecule.total_mass();
  center.z /= molecule.total_mass();
  return center;
}

calc_type calculate_radius_of_gyration(Molecule& molecule, std::size_t frame) {
  const auto center = calculate_center_of_mass(molecule, frame);
  calc_type rg_squared{};

  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto xyz = molecule.coordinate(frame, atom);
    const auto dx = static_cast<calc_type>(xyz.x) - center.x;
    const auto dy = static_cast<calc_type>(xyz.y) - center.y;
    const auto dz = static_cast<calc_type>(xyz.z) - center.z;
    rg_squared += dx * dx + dy * dy + dz * dz;
  }

  return std::sqrt(rg_squared / static_cast<calc_type>(molecule.natoms()));
}

calc_type calculate_root_mean_square_deviation(const Molecule& first,
                                               const Molecule& second) {
  if (first.natoms() == 0) {
    throw std::invalid_argument(
        "calculate_root_mean_square_deviation requires atoms");
  }
  if (first.natoms() != second.natoms() ||
      first.number_of_frames() != second.number_of_frames() ||
      first.coor().size() != second.coor().size()) {
    throw std::invalid_argument(
        "calculate_root_mean_square_deviation requires matching coordinates");
  }

  calc_type sum{};
  for (std::size_t index = 0; index < first.coor().size(); ++index) {
    const auto delta = static_cast<calc_type>(first.coor()[index]) -
                       static_cast<calc_type>(second.coor()[index]);
    sum += delta * delta;
  }

  return std::sqrt(sum / static_cast<calc_type>(first.natoms()));
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
