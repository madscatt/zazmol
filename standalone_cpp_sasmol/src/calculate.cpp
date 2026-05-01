#include "sasmol/calculate.hpp"

#include "sasmol/properties.hpp"

#include <algorithm>
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

calc_type determinant(const CalcMatrix3& matrix) {
  return matrix[0][0] * (matrix[1][1] * matrix[2][2] -
                         matrix[1][2] * matrix[2][1]) -
         matrix[0][1] * (matrix[1][0] * matrix[2][2] -
                         matrix[1][2] * matrix[2][0]) +
         matrix[0][2] * (matrix[1][0] * matrix[2][1] -
                         matrix[1][1] * matrix[2][0]);
}

void sort_eigens_ascending(std::array<calc_type, 3>& eigenvalues,
                           CalcMatrix3& eigenvectors) {
  std::array<std::size_t, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return eigenvalues[lhs] < eigenvalues[rhs];
  });

  const auto original_values = eigenvalues;
  const auto original_vectors = eigenvectors;
  for (std::size_t column = 0; column < 3; ++column) {
    eigenvalues[column] = original_values[order[column]];
    for (std::size_t row = 0; row < 3; ++row) {
      eigenvectors[row][column] = original_vectors[row][order[column]];
    }
  }
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

  std::array<calc_type, 3> eigenvalues{matrix[0][0], matrix[1][1],
                                       matrix[2][2]};
  sort_eigens_ascending(eigenvalues, eigenvectors);
  return {eigenvalues, eigenvectors};
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

PrincipalMomentsOfInertia calculate_principal_moments_of_inertia(
    Molecule& molecule, std::size_t frame) {
  if (frame >= molecule.number_of_frames()) {
    throw std::out_of_range(
        "calculate_principal_moments_of_inertia frame is out of range");
  }
  if (molecule.natoms() == 0) {
    throw std::invalid_argument(
        "calculate_principal_moments_of_inertia requires atoms");
  }
  if (molecule.mass().size() != molecule.natoms() ||
      molecule.total_mass() <= calc_type{}) {
    const auto mass_result = calculate_mass(molecule);
    if (!mass_result.ok()) {
      throw std::invalid_argument(
          "calculate_principal_moments_of_inertia requires known element masses");
    }
  }

  const auto center = calculate_center_of_mass(molecule, frame);
  PrincipalMomentsOfInertia result;

  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto xyz = molecule.coordinate(frame, atom);
    const auto mass = molecule.mass()[atom];
    const auto x = static_cast<calc_type>(xyz.x) - center.x;
    const auto y = static_cast<calc_type>(xyz.y) - center.y;
    const auto z = static_cast<calc_type>(xyz.z) - center.z;
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(mass)) {
      throw std::invalid_argument(
          "calculate_principal_moments_of_inertia requires finite values");
    }

    result.inertia[0][0] += mass * (y * y + z * z);
    result.inertia[1][1] += mass * (x * x + z * z);
    result.inertia[2][2] += mass * (x * x + y * y);
    result.inertia[0][1] -= mass * x * y;
    result.inertia[0][2] -= mass * x * z;
    result.inertia[1][2] -= mass * y * z;
  }
  result.inertia[1][0] = result.inertia[0][1];
  result.inertia[2][0] = result.inertia[0][2];
  result.inertia[2][1] = result.inertia[1][2];

  const auto scale = std::max<calc_type>(1.0, max_abs_diagonal(result.inertia));
  const auto singular_tolerance =
      3.0 * std::numeric_limits<calc_type>::epsilon() * scale * scale * scale;
  result.singular = std::fabs(determinant(result.inertia)) <= singular_tolerance;
  if (!result.singular) {
    auto eigens = symmetric_eigen_decomposition(result.inertia);
    result.eigenvalues = eigens.first;
    result.eigenvectors = eigens.second;
  }
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
