#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace sasmol {

struct SubsetResult {
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CoordinateSelection {
  std::vector<Vec3> coordinates;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

[[nodiscard]] CoordinateSelection get_coordinates_using_indices(
    const Molecule& molecule, std::size_t frame,
    const std::vector<std::size_t>& indices);

[[nodiscard]] SubsetResult set_coordinates_using_indices(
    Molecule& molecule, const Molecule& source, std::size_t frame,
    const std::vector<std::size_t>& indices);

[[nodiscard]] SubsetResult copy_molecule_using_indices(
    const Molecule& source, Molecule& destination,
    const std::vector<std::size_t>& indices, std::size_t frame);

[[nodiscard]] Molecule copied_molecule_using_indices(
    const Molecule& source, const std::vector<std::size_t>& indices,
    std::size_t frame);

}  // namespace sasmol
