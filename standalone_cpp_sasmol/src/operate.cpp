#include "sasmol/operate.hpp"

#include <cmath>

namespace sasmol {
namespace {

Vec3 transform_column_vector(Vec3 xyz, const CalcMatrix3& matrix) {
  return {static_cast<coord_type>(matrix[0][0] * xyz.x +
                                  matrix[0][1] * xyz.y +
                                  matrix[0][2] * xyz.z),
          static_cast<coord_type>(matrix[1][0] * xyz.x +
                                  matrix[1][1] * xyz.y +
                                  matrix[1][2] * xyz.z),
          static_cast<coord_type>(matrix[2][0] * xyz.x +
                                  matrix[2][1] * xyz.y +
                                  matrix[2][2] * xyz.z)};
}

Vec3 transform_row_vector(Vec3 xyz, const CalcMatrix3& matrix) {
  return {static_cast<coord_type>(xyz.x * matrix[0][0] +
                                  xyz.y * matrix[1][0] +
                                  xyz.z * matrix[2][0]),
          static_cast<coord_type>(xyz.x * matrix[0][1] +
                                  xyz.y * matrix[1][1] +
                                  xyz.z * matrix[2][1]),
          static_cast<coord_type>(xyz.x * matrix[0][2] +
                                  xyz.y * matrix[1][2] +
                                  xyz.z * matrix[2][2])};
}

}  // namespace

Rotation axis_rotation(Axis axis, calc_type theta) {
  const auto cs = std::cos(theta);
  const auto si = std::sin(theta);
  Rotation rotation;
  rotation.convention = CoordinateTransformConvention::column_vector;

  switch (axis) {
    case Axis::x:
      rotation.matrix = {{{1.0, 0.0, 0.0}, {0.0, cs, -si}, {0.0, si, cs}}};
      break;
    case Axis::y:
      rotation.matrix = {{{cs, 0.0, si}, {0.0, 1.0, 0.0}, {-si, 0.0, cs}}};
      break;
    case Axis::z:
      rotation.matrix = {{{cs, -si, 0.0}, {si, cs, 0.0}, {0.0, 0.0, 1.0}}};
      break;
  }

  return rotation;
}

Rotation general_axis_rotation(calc_type theta, CalcVec3 unit_axis) {
  const auto ux = unit_axis.x;
  const auto uy = unit_axis.y;
  const auto uz = unit_axis.z;
  const auto cs = std::cos(theta);
  const auto si = std::sin(theta);
  const auto one_minus_cs = 1.0 - cs;

  Rotation rotation;
  rotation.convention = CoordinateTransformConvention::row_vector;
  rotation.matrix = {{{cs + ux * ux * one_minus_cs,
                       ux * uy * one_minus_cs - uz * si,
                       ux * uz * one_minus_cs + uy * si},
                      {uy * ux * one_minus_cs + uz * si,
                       cs + uy * uy * one_minus_cs,
                       uy * uz * one_minus_cs - ux * si},
                      {uz * ux * one_minus_cs - uy * si,
                       uz * uy * one_minus_cs + ux * si,
                       cs + uz * uz * one_minus_cs}}};
  return rotation;
}

Rotation euler_rotation(calc_type phi, calc_type theta, calc_type psi) {
  Rotation rotation;
  rotation.convention = CoordinateTransformConvention::row_vector;
  rotation.matrix = {
      {{std::cos(theta) * std::cos(psi),
        std::cos(phi) * std::sin(psi) +
            std::sin(phi) * std::sin(theta) * std::cos(psi),
        std::sin(phi) * std::sin(psi) -
            std::cos(phi) * std::sin(theta) * std::cos(psi)},
       {-std::cos(theta) * std::sin(psi),
        std::cos(phi) * std::cos(psi) -
            std::sin(phi) * std::sin(theta) * std::sin(psi),
        std::sin(phi) * std::cos(psi) +
            std::cos(phi) * std::sin(theta) * std::sin(psi)},
       {std::sin(theta), -std::sin(phi) * std::cos(theta),
        std::cos(phi) * std::cos(theta)}}};
  return rotation;
}

void apply_translation(CoordinateView coordinates, CalcVec3 delta) {
  for (std::size_t atom = 0; atom < coordinates.natoms; ++atom) {
    const auto xyz = coordinates[atom];
    coordinates.set(atom, {static_cast<coord_type>(xyz.x + delta.x),
                           static_cast<coord_type>(xyz.y + delta.y),
                           static_cast<coord_type>(xyz.z + delta.z)});
  }
}

void apply_rotation(CoordinateView coordinates, const Rotation& rotation) {
  for (std::size_t atom = 0; atom < coordinates.natoms; ++atom) {
    const auto xyz = coordinates[atom];
    if (rotation.convention == CoordinateTransformConvention::column_vector) {
      coordinates.set(atom, transform_column_vector(xyz, rotation.matrix));
    } else {
      coordinates.set(atom, transform_row_vector(xyz, rotation.matrix));
    }
  }
}

void translate(Molecule& molecule, std::size_t frame, CalcVec3 value,
               TranslationOptions options) {
  if (options.point) {
    const auto com = calculate_center_of_mass(molecule, frame);
    apply_translation(molecule.coordinate_view(frame),
                      {-com.x + value.x, -com.y + value.y, -com.z + value.z});
  } else {
    apply_translation(molecule.coordinate_view(frame), value);
  }
}

Molecule translated(const Molecule& molecule, std::size_t frame, CalcVec3 value,
                    TranslationOptions options) {
  auto copy = molecule;
  translate(copy, frame, value, options);
  return copy;
}

void center(Molecule& molecule, std::size_t frame) {
  const auto com = calculate_center_of_mass(molecule, frame);
  apply_translation(molecule.coordinate_view(frame), {-com.x, -com.y, -com.z});
}

Molecule centered(const Molecule& molecule, std::size_t frame) {
  auto copy = molecule;
  center(copy, frame);
  return copy;
}

void rotate(Molecule& molecule, std::size_t frame, Axis axis, calc_type theta) {
  apply_rotation(molecule.coordinate_view(frame), axis_rotation(axis, theta));
}

Molecule rotated(const Molecule& molecule, std::size_t frame, Axis axis,
                 calc_type theta) {
  auto copy = molecule;
  rotate(copy, frame, axis, theta);
  return copy;
}

void rotate_general_axis(Molecule& molecule, std::size_t frame, calc_type theta,
                         CalcVec3 unit_axis) {
  apply_rotation(molecule.coordinate_view(frame),
                 general_axis_rotation(theta, unit_axis));
}

Molecule rotated_general_axis(const Molecule& molecule, std::size_t frame,
                              calc_type theta, CalcVec3 unit_axis) {
  auto copy = molecule;
  rotate_general_axis(copy, frame, theta, unit_axis);
  return copy;
}

void rotate_euler(Molecule& molecule, std::size_t frame, calc_type phi,
                  calc_type theta, calc_type psi) {
  apply_rotation(molecule.coordinate_view(frame), euler_rotation(phi, theta, psi));
}

Molecule rotated_euler(const Molecule& molecule, std::size_t frame,
                       calc_type phi, calc_type theta, calc_type psi) {
  auto copy = molecule;
  rotate_euler(copy, frame, phi, theta, psi);
  return copy;
}

}  // namespace sasmol
