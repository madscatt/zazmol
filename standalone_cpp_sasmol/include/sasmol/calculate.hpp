#pragma once

#include "sasmol/molecule.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
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

using CalcMatrix3 = std::array<std::array<calc_type, 3>, 3>;

struct PrincipalMomentsOfInertia {
  std::array<calc_type, 3> eigenvalues{};
  CalcMatrix3 eigenvectors{};
  CalcMatrix3 inertia{};
  bool singular{};
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

[[nodiscard]] PrincipalMomentsOfInertia calculate_principal_moments_of_inertia(
    Molecule& molecule, std::size_t frame);

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum(
    const Molecule& molecule, const std::vector<std::size_t>& frames = {});

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum_all_steps(
    const Molecule& molecule);

[[nodiscard]] CoordinateBounds calculate_minimum_and_maximum_all_steps(
    const std::filesystem::path& trajectory_filename);

[[nodiscard]] CoordinateBounds calc_minmax_all_steps(const Molecule& molecule);

[[nodiscard]] CoordinateBounds calc_minmax_all_steps(
    const std::filesystem::path& trajectory_filename);

}  // namespace sasmol
