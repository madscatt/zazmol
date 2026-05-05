#include "sasmol/operate.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

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

CalcVec3 axis_vector(Axis axis) {
  switch (axis) {
    case Axis::x:
      return {1.0, 0.0, 0.0};
    case Axis::y:
      return {0.0, 1.0, 0.0};
    case Axis::z:
      return {0.0, 0.0, 1.0};
  }
  return {};
}

const std::map<std::string, calc_type>& legacy_average_vdw_radii() {
  static const std::map<std::string, calc_type> radii{
      {"H", 0.928619230769}, {"C", 2.00249333333},
      {"N", 1.85},           {"O", 1.7392625},
      {"F", 1.7},            {"Ne", 1.5300},
      {"Na", 1.36375},       {"Mg", 1.18500},
      {"P", 2.15},           {"S", 2.000000},
      {"Cl", 2.27},          {"K", 1.76375},
      {"Ca", 1.367},         {"Fe", 0.650000},
      {"Zn", 1.09000},       {"Cs", 2.100},
      {"D", 0.928619230769}, {"1H", 0.928619230769}};
  return radii;
}

CalcVec3 matrix_column(const CalcMatrix3& matrix, std::size_t column) {
  return {matrix[0][column], matrix[1][column], matrix[2][column]};
}

calc_type dot(CalcVec3 lhs, CalcVec3 rhs) {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

CalcVec3 cross(CalcVec3 lhs, CalcVec3 rhs) {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

calc_type norm(CalcVec3 value) {
  return std::sqrt(dot(value, value));
}

CalcVec3 scaled(CalcVec3 value, calc_type scale) {
  return {value.x * scale, value.y * scale, value.z * scale};
}

CalcVec3 normalized(CalcVec3 value) {
  const auto length = norm(value);
  if (length <= 0.0 || !std::isfinite(length)) {
    throw std::invalid_argument("cannot normalize zero or non-finite vector");
  }
  return scaled(value, 1.0 / length);
}

std::array<std::size_t, 3> eigen_order_descending(
    const std::array<calc_type, 3>& values) {
  std::array<std::size_t, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return values[lhs] > values[rhs];
  });
  return order;
}

CalcMatrix3 identity_matrix3() {
  return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

calc_type max_abs_diagonal(const CalcMatrix3& matrix) {
  calc_type value{};
  for (std::size_t i = 0; i < 3; ++i) {
    value = std::max(value, std::fabs(matrix[i][i]));
  }
  return value;
}

std::pair<std::array<calc_type, 3>, CalcMatrix3> symmetric_eigen_decomposition(
    const CalcMatrix3& input) {
  auto matrix = input;
  auto eigenvectors = identity_matrix3();

  for (std::size_t iteration = 0; iteration < 60; ++iteration) {
    std::size_t p = 0;
    std::size_t q = 1;
    auto offdiag = std::fabs(matrix[p][q]);
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t column = row + 1; column < 3; ++column) {
        const auto value = std::fabs(matrix[row][column]);
        if (value > offdiag) {
          offdiag = value;
          p = row;
          q = column;
        }
      }
    }

    const auto scale = std::max<calc_type>(1.0, max_abs_diagonal(matrix));
    if (offdiag <= std::numeric_limits<calc_type>::epsilon() * scale) {
      break;
    }

    const auto app = matrix[p][p];
    const auto aqq = matrix[q][q];
    const auto apq = matrix[p][q];
    const auto angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
    const auto c = std::cos(angle);
    const auto s = std::sin(angle);

    matrix[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    matrix[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    matrix[p][q] = 0.0;
    matrix[q][p] = 0.0;

    for (std::size_t k = 0; k < 3; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const auto akp = matrix[k][p];
      const auto akq = matrix[k][q];
      matrix[k][p] = c * akp - s * akq;
      matrix[p][k] = matrix[k][p];
      matrix[k][q] = s * akp + c * akq;
      matrix[q][k] = matrix[k][q];
    }

    for (std::size_t row = 0; row < 3; ++row) {
      const auto vkp = eigenvectors[row][p];
      const auto vkq = eigenvectors[row][q];
      eigenvectors[row][p] = c * vkp - s * vkq;
      eigenvectors[row][q] = s * vkp + c * vkq;
    }
  }

  return {{matrix[0][0], matrix[1][1], matrix[2][2]}, eigenvectors};
}

CalcVec3 matrix_vector_multiply(const CalcMatrix3& matrix, CalcVec3 value) {
  return {matrix[0][0] * value.x + matrix[0][1] * value.y +
              matrix[0][2] * value.z,
          matrix[1][0] * value.x + matrix[1][1] * value.y +
              matrix[1][2] * value.z,
          matrix[2][0] * value.x + matrix[2][1] * value.y +
              matrix[2][2] * value.z};
}

Rotation alignment_rotation(const std::vector<Vec3>& reference_centered,
                            const std::vector<Vec3>& moving_centered) {
  if (reference_centered.size() != moving_centered.size() ||
      reference_centered.empty()) {
    throw std::invalid_argument(
        "alignment requires equal non-empty reference and moving bases");
  }

  CalcMatrix3 r{};
  for (std::size_t atom = 0; atom < reference_centered.size(); ++atom) {
    const auto x = reference_centered[atom];
    const auto y = moving_centered[atom];
    const calc_type ref[3] = {x.x, x.y, x.z};
    const calc_type mov[3] = {y.x, y.y, y.z};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        r[i][j] += mov[i] * ref[j];
      }
    }
  }

  CalcMatrix3 rtr{};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      for (std::size_t k = 0; k < 3; ++k) {
        rtr[i][j] += r[k][i] * r[k][j];
      }
    }
  }

  auto eigens = symmetric_eigen_decomposition(rtr);
  const auto order = eigen_order_descending(eigens.first);
  std::array<calc_type, 3> uk{};
  CalcMatrix3 ak{};
  for (std::size_t row = 0; row < 3; ++row) {
    uk[row] = eigens.first[order[row]];
    for (std::size_t column = 0; column < 3; ++column) {
      ak[row][column] = eigens.second[column][order[row]];
    }
  }
  const auto ak2 = cross({ak[0][0], ak[0][1], ak[0][2]},
                         {ak[1][0], ak[1][1], ak[1][2]});
  ak[2] = {ak2.x, ak2.y, ak2.z};

  auto rak0 = matrix_vector_multiply(r, {ak[0][0], ak[0][1], ak[0][2]});
  auto rak1 = matrix_vector_multiply(r, {ak[1][0], ak[1][1], ak[1][2]});
  auto urak0 = scaled(rak0, uk[0] == 0.0 ? 1.0e15 : 1.0 / std::sqrt(std::fabs(uk[0])));
  auto urak1 = scaled(rak1, uk[1] == 0.0 ? 1.0e15 : 1.0 / std::sqrt(std::fabs(uk[1])));
  auto urak2 = cross(urak0, urak1);

  CalcMatrix3 bk{{{urak0.x, urak0.y, urak0.z},
                  {urak1.x, urak1.y, urak1.z},
                  {urak2.x, urak2.y, urak2.z}}};

  Rotation rotation;
  rotation.convention = CoordinateTransformConvention::column_vector;
  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        rotation.matrix[j][i] += bk[k][i] * ak[k][j];
      }
    }
  }
  return rotation;
}

std::vector<Vec3> centered_coordinates(const Molecule& molecule,
                                       const std::vector<std::size_t>& indices,
                                       std::size_t frame, CalcVec3 center) {
  std::vector<Vec3> coordinates;
  coordinates.reserve(indices.size());
  for (const auto atom : indices) {
    const auto xyz = molecule.coordinate(frame, atom);
    coordinates.push_back({static_cast<coord_type>(xyz.x - center.x),
                           static_cast<coord_type>(xyz.y - center.y),
                           static_cast<coord_type>(xyz.z - center.z)});
  }
  return coordinates;
}

CalcVec3 selected_center_of_mass(const Molecule& molecule,
                                 const std::vector<std::size_t>& indices,
                                 std::size_t frame) {
  if (frame >= molecule.number_of_frames()) {
    throw std::out_of_range("selected center of mass frame is out of range");
  }
  if (indices.empty()) {
    throw std::invalid_argument("selected center of mass requires atoms");
  }

  auto mass_source = molecule;
  if (mass_source.mass().size() != mass_source.natoms() ||
      mass_source.total_mass() <= calc_type{}) {
    const auto mass_result = calculate_mass(mass_source);
    if (!mass_result.ok()) {
      throw std::invalid_argument(
          "selected center of mass requires known element masses");
    }
  }

  CalcVec3 center;
  calc_type total_mass{};
  for (const auto atom : indices) {
    if (atom >= mass_source.natoms()) {
      throw std::out_of_range("selected center of mass atom is out of range");
    }
    const auto xyz = mass_source.coordinate(frame, atom);
    const auto mass = mass_source.mass()[atom];
    center.x += mass * static_cast<calc_type>(xyz.x);
    center.y += mass * static_cast<calc_type>(xyz.y);
    center.z += mass * static_cast<calc_type>(xyz.z);
    total_mass += mass;
  }

  if (total_mass <= calc_type{}) {
    throw std::invalid_argument(
        "selected center of mass requires positive mass");
  }

  center.x /= total_mass;
  center.y /= total_mass;
  center.z /= total_mass;
  return center;
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

SubsetResult set_average_vdw(Molecule& molecule) {
  SubsetResult result;
  if (molecule.element().size() != molecule.natoms()) {
    result.errors.push_back(
        "set_average_vdw requires descriptor 'element' length to match natoms");
    return result;
  }

  const auto& radii = legacy_average_vdw_radii();
  std::vector<calc_type> atom_vdw;
  atom_vdw.reserve(molecule.natoms());
  for (const auto& element : molecule.element()) {
    const auto found = radii.find(element);
    atom_vdw.push_back(found == radii.end() ? calc_type{} : found->second);
  }
  molecule.atom_vdw() = std::move(atom_vdw);
  return result;
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

void align_pmi_on_axis(Molecule& molecule, std::size_t frame,
                       std::size_t pmi_eigenvector, Axis alignment_axis) {
  if (pmi_eigenvector >= 3) {
    throw std::out_of_range("PMI eigenvector index is out of range");
  }

  center(molecule, frame);
  const auto pmi = calculate_principal_moments_of_inertia(molecule, frame);
  if (pmi.singular) {
    throw std::invalid_argument("cannot align singular PMI tensor");
  }

  auto pmi_axis = normalized(matrix_column(pmi.eigenvectors, pmi_eigenvector));
  const auto target_axis = axis_vector(alignment_axis);
  auto rotation_axis = cross(target_axis, pmi_axis);
  const auto sine = norm(rotation_axis);
  const auto cosine = std::clamp(dot(pmi_axis, target_axis), -1.0, 1.0);

  if (sine > 1.0e-12) {
    rotation_axis = scaled(rotation_axis, 1.0 / sine);
  } else if (cosine > 0.0) {
    return;
  } else {
    const CalcVec3 bases[3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                               {0.0, 0.0, 1.0}};
    std::size_t basis_index = 0;
    auto smallest_component = std::fabs(pmi_axis.x);
    if (std::fabs(pmi_axis.y) < smallest_component) {
      basis_index = 1;
      smallest_component = std::fabs(pmi_axis.y);
    }
    if (std::fabs(pmi_axis.z) < smallest_component) {
      basis_index = 2;
    }
    rotation_axis = normalized(cross(pmi_axis, bases[basis_index]));
  }

  rotate_general_axis(molecule, frame, std::atan2(sine, cosine), rotation_axis);
}

Molecule pmi_aligned_on_axis(const Molecule& molecule, std::size_t frame,
                             std::size_t pmi_eigenvector,
                             Axis alignment_axis) {
  auto copy = molecule;
  align_pmi_on_axis(copy, frame, pmi_eigenvector, alignment_axis);
  return copy;
}

void align_pmi_on_cardinal_axes(Molecule& molecule, std::size_t frame) {
  align_pmi_on_axis(molecule, frame, 2, Axis::z);
  align_pmi_on_axis(molecule, frame, 1, Axis::y);
}

Molecule pmi_aligned_on_cardinal_axes(const Molecule& molecule,
                                      std::size_t frame) {
  auto copy = molecule;
  align_pmi_on_cardinal_axes(copy, frame);
  return copy;
}

AlignmentPlan initialize_alignment(
    const Molecule& moving, const Molecule& reference,
    const std::vector<std::size_t>& moving_basis_indices,
    const std::vector<std::size_t>& reference_basis_indices, std::size_t frame) {
  if (moving_basis_indices.size() != reference_basis_indices.size() ||
      moving_basis_indices.empty()) {
    throw std::invalid_argument(
        "alignment initialization requires matching non-empty basis indices");
  }
  if (frame >= moving.number_of_frames() || frame >= reference.number_of_frames()) {
    throw std::out_of_range("alignment frame is out of range");
  }

  AlignmentPlan plan;
  plan.moving_basis_indices = moving_basis_indices;
  const auto reference_center =
      selected_center_of_mass(reference, reference_basis_indices, frame);
  const auto moving_center =
      selected_center_of_mass(moving, moving_basis_indices, frame);
  plan.reference_center_of_mass = reference_center;
  plan.centered_reference_basis =
      centered_coordinates(reference, reference_basis_indices, frame,
                           reference_center);
  (void)centered_coordinates(moving, moving_basis_indices, frame, moving_center);
  return plan;
}

void align(Molecule& moving, const AlignmentPlan& plan, std::size_t frame) {
  if (plan.moving_basis_indices.empty() ||
      plan.centered_reference_basis.size() != plan.moving_basis_indices.size()) {
    throw std::invalid_argument("alignment plan is incomplete");
  }
  if (frame >= moving.number_of_frames()) {
    throw std::out_of_range("alignment frame is out of range");
  }

  const auto moving_center =
      selected_center_of_mass(moving, plan.moving_basis_indices, frame);
  const auto moving_centered =
      centered_coordinates(moving, plan.moving_basis_indices, frame,
                           moving_center);
  const auto rotation =
      alignment_rotation(plan.centered_reference_basis, moving_centered);

  auto view = moving.coordinate_view(frame);
  for (std::size_t atom = 0; atom < view.natoms; ++atom) {
    const auto xyz = view[atom];
    const Vec3 centered{static_cast<coord_type>(xyz.x - moving_center.x),
                        static_cast<coord_type>(xyz.y - moving_center.y),
                        static_cast<coord_type>(xyz.z - moving_center.z)};
    const auto rotated_xyz = transform_column_vector(centered, rotation.matrix);
    view.set(atom, {static_cast<coord_type>(rotated_xyz.x +
                                            plan.reference_center_of_mass.x),
                    static_cast<coord_type>(rotated_xyz.y +
                                            plan.reference_center_of_mass.y),
                    static_cast<coord_type>(rotated_xyz.z +
                                            plan.reference_center_of_mass.z)});
  }
}

Molecule aligned(const Molecule& moving, const AlignmentPlan& plan,
                 std::size_t frame) {
  auto copy = moving;
  align(copy, plan, frame);
  return copy;
}

}  // namespace sasmol
