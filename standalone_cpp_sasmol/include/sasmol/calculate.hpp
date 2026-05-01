#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <vector>

namespace sasmol {

struct CoordinateBounds {
  Vec3 minimum;
  Vec3 maximum;
};

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum(
    const Molecule& molecule, const std::vector<std::size_t>& frames = {});

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum_all_steps(
    const Molecule& molecule);

[[nodiscard]] CoordinateBounds calc_minmax_all_steps(const Molecule& molecule);

}  // namespace sasmol
