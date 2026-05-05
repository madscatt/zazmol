#pragma once

#include "sasmol/calculate.hpp"
#include "sasmol/molecule.hpp"
#include "sasmol/subset.hpp"

#include <cstddef>
#include <vector>

namespace sasmol {

enum class Axis {
  x,
  y,
  z,
};

enum class CoordinateTransformConvention {
  column_vector,
  row_vector,
};

struct TranslationOptions {
  bool point{false};
};

struct Rotation {
  CalcMatrix3 matrix{};
  CoordinateTransformConvention convention{
      CoordinateTransformConvention::column_vector};
};

struct AlignmentPlan {
  std::vector<std::size_t> moving_basis_indices;
  std::vector<Vec3> centered_reference_basis;
  CalcVec3 reference_center_of_mass{};
};

[[nodiscard]] Rotation axis_rotation(Axis axis, calc_type theta);
[[nodiscard]] Rotation general_axis_rotation(calc_type theta, CalcVec3 unit_axis);
[[nodiscard]] Rotation euler_rotation(calc_type phi, calc_type theta,
                                      calc_type psi);

void apply_translation(CoordinateView coordinates, CalcVec3 delta);
void apply_rotation(CoordinateView coordinates, const Rotation& rotation);

[[nodiscard]] SubsetResult set_average_vdw(Molecule& molecule);

void translate(Molecule& molecule, std::size_t frame, CalcVec3 value,
               TranslationOptions options = {});
[[nodiscard]] Molecule translated(const Molecule& molecule, std::size_t frame,
                                  CalcVec3 value,
                                  TranslationOptions options = {});

void center(Molecule& molecule, std::size_t frame);
[[nodiscard]] Molecule centered(const Molecule& molecule, std::size_t frame);

void rotate(Molecule& molecule, std::size_t frame, Axis axis, calc_type theta);
[[nodiscard]] Molecule rotated(const Molecule& molecule, std::size_t frame,
                               Axis axis, calc_type theta);

void rotate_general_axis(Molecule& molecule, std::size_t frame, calc_type theta,
                         CalcVec3 unit_axis);
[[nodiscard]] Molecule rotated_general_axis(const Molecule& molecule,
                                            std::size_t frame, calc_type theta,
                                            CalcVec3 unit_axis);

void rotate_euler(Molecule& molecule, std::size_t frame, calc_type phi,
                  calc_type theta, calc_type psi);
[[nodiscard]] Molecule rotated_euler(const Molecule& molecule, std::size_t frame,
                                     calc_type phi, calc_type theta,
                                     calc_type psi);

void align_pmi_on_axis(Molecule& molecule, std::size_t frame,
                       std::size_t pmi_eigenvector, Axis alignment_axis);
[[nodiscard]] Molecule pmi_aligned_on_axis(const Molecule& molecule,
                                           std::size_t frame,
                                           std::size_t pmi_eigenvector,
                                           Axis alignment_axis);

void align_pmi_on_cardinal_axes(Molecule& molecule, std::size_t frame);
[[nodiscard]] Molecule pmi_aligned_on_cardinal_axes(const Molecule& molecule,
                                                    std::size_t frame);

[[nodiscard]] AlignmentPlan initialize_alignment(
    const Molecule& moving, const Molecule& reference,
    const std::vector<std::size_t>& moving_basis_indices,
    const std::vector<std::size_t>& reference_basis_indices, std::size_t frame);

void align(Molecule& moving, const AlignmentPlan& plan, std::size_t frame);
[[nodiscard]] Molecule aligned(const Molecule& moving, const AlignmentPlan& plan,
                               std::size_t frame);

}  // namespace sasmol
