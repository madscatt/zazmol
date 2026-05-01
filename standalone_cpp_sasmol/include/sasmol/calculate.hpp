#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace sasmol {

struct CoordinateBounds {
  Vec3 minimum;
  Vec3 maximum;
};

struct CalcVec3 {
  calc_type x{};
  calc_type y{};
  calc_type z{};
};

struct MassCalculationResult {
  calc_type total_mass{};
  std::vector<std::string> unknown_elements;

  [[nodiscard]] bool ok() const noexcept { return unknown_elements.empty(); }
};

[[nodiscard]] MassCalculationResult calculate_mass(Molecule& molecule);

[[nodiscard]] std::map<std::string, std::size_t> calculate_molecular_formula(
    Molecule& molecule);

void calculate_residue_charge(Molecule& molecule);

[[nodiscard]] CalcVec3 calculate_center_of_mass(Molecule& molecule,
                                                std::size_t frame);

[[nodiscard]] calc_type calculate_radius_of_gyration(Molecule& molecule,
                                                     std::size_t frame);

[[nodiscard]] calc_type calculate_root_mean_square_deviation(
    const Molecule& first, const Molecule& second);

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum(
    const Molecule& molecule, const std::vector<std::size_t>& frames = {});

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum_all_steps(
    const Molecule& molecule);

[[nodiscard]] CoordinateBounds calc_minmax_all_steps(const Molecule& molecule);

}  // namespace sasmol
