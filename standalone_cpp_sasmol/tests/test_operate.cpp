#include "sasmol/file_io.hpp"
#include "sasmol/operate.hpp"
#include "sasmol/selection.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <vector>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-6) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

void assert_close_double(sasmol::calc_type actual, sasmol::calc_type expected,
                         double tolerance = 1.0e-6) {
  assert(std::fabs(actual - expected) < tolerance);
}

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_vec_close(sasmol::Vec3 actual, sasmol::Vec3 expected,
                      double tolerance = 1.0e-6) {
  assert_close(actual.x, expected.x, tolerance);
  assert_close(actual.y, expected.y, tolerance);
  assert_close(actual.z, expected.z, tolerance);
}

double axis_alignment(const sasmol::Molecule& mol, std::size_t frame,
                      std::size_t pmi_eigenvector,
                      std::array<double, 3> expected_axis) {
  auto copy = mol;
  const auto pmi = sasmol::calculate_principal_moments_of_inertia(copy, frame);
  assert(!pmi.singular);

  double dot{};
  double expected_norm{};
  for (std::size_t row = 0; row < 3; ++row) {
    dot += pmi.eigenvectors[row][pmi_eigenvector] * expected_axis[row];
    expected_norm += expected_axis[row] * expected_axis[row];
  }
  return std::fabs(dot) / std::sqrt(expected_norm);
}

double selected_rmsd(const sasmol::Molecule& first,
                     const sasmol::Molecule& second,
                     const std::vector<std::size_t>& indices,
                     std::size_t frame) {
  double sum{};
  for (const auto atom : indices) {
    const auto lhs = first.coordinate(frame, atom);
    const auto rhs = second.coordinate(frame, atom);
    const auto dx = static_cast<double>(lhs.x - rhs.x);
    const auto dy = static_cast<double>(lhs.y - rhs.y);
    const auto dz = static_cast<double>(lhs.z - rhs.z);
    sum += dx * dx + dy * dy + dz * dz;
  }
  return std::sqrt(sum / static_cast<double>(indices.size()));
}

std::vector<std::size_t> all_atom_indices(const sasmol::Molecule& molecule) {
  std::vector<std::size_t> indices;
  indices.reserve(molecule.natoms());
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    indices.push_back(atom);
  }
  return indices;
}

void test_translate_mutates_only_selected_frame() {
  sasmol::Molecule mol(1, 2);
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(1, 0, {4.0F, 5.0F, 6.0F});

  sasmol::translate(mol, 1, {10.0, -2.0, 0.5});

  assert_vec_close(mol.coordinate(0, 0), {1.0F, 2.0F, 3.0F});
  assert_vec_close(mol.coordinate(1, 0), {14.0F, 3.0F, 6.5F});
}

void test_translated_returns_copy_without_mutating_source() {
  sasmol::Molecule mol(1, 1);
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});

  const auto copy = sasmol::translated(mol, 0, {3.0, 4.0, 5.0});

  assert_vec_close(mol.coordinate(0, 0), {1.0F, 2.0F, 3.0F});
  assert_vec_close(copy.coordinate(0, 0), {4.0F, 6.0F, 8.0F});
}

void test_center_moves_com_to_origin() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  sasmol::center(mol, 0);
  const auto com = sasmol::calculate_center_of_mass(mol, 0);

  assert_close_double(com.x, 0.0, 1.0e-5);
  assert_close_double(com.y, 0.0, 1.0e-5);
  assert_close_double(com.z, 0.0, 1.0e-5);
}

void test_translate_point_moves_com_to_target() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  sasmol::translate(mol, 0, {3.0, 4.0, 5.0}, {.point = true});
  const auto com = sasmol::calculate_center_of_mass(mol, 0);

  assert_close_double(com.x, 3.0, 1.0e-5);
  assert_close_double(com.y, 4.0, 1.0e-5);
  assert_close_double(com.z, 5.0, 1.0e-5);
}

void test_axis_rotation_matches_python_column_vector_convention() {
  sasmol::Molecule mol(3, 1);
  mol.set_coordinate(0, 0, {0.0F, 1.0F, 0.0F});
  mol.set_coordinate(0, 1, {0.0F, 0.0F, 1.0F});
  mol.set_coordinate(0, 2, {1.0F, 0.0F, 0.0F});

  sasmol::rotate(mol, 0, sasmol::Axis::x, std::acos(-1.0) / 2.0);
  assert_vec_close(mol.coordinate(0, 0), {0.0F, 0.0F, 1.0F});
  assert_vec_close(mol.coordinate(0, 1), {0.0F, -1.0F, 0.0F});

  sasmol::rotate(mol, 0, sasmol::Axis::z, std::acos(-1.0) / 2.0);
  assert_vec_close(mol.coordinate(0, 2), {0.0F, 1.0F, 0.0F});
}

void test_general_axis_rotation_matches_python_row_vector_convention() {
  sasmol::Molecule mol(1, 1);
  mol.set_coordinate(0, 0, {1.0F, 0.0F, 0.0F});

  sasmol::rotate_general_axis(mol, 0, std::acos(-1.0) / 2.0,
                              {0.0, 0.0, 1.0});

  assert_vec_close(mol.coordinate(0, 0), {0.0F, -1.0F, 0.0F});
}

void test_euler_rotation_identity() {
  sasmol::Molecule mol(1, 1);
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});

  sasmol::rotate_euler(mol, 0, 0.0, 0.0, 0.0);

  assert_vec_close(mol.coordinate(0, 0), {1.0F, 2.0F, 3.0F});
}

void test_general_axis_preserves_python_non_unit_axis_behavior() {
  sasmol::Molecule mol(1, 1);
  mol.set_coordinate(0, 0, {73.944F, 41.799F, 41.652F});

  sasmol::rotate_general_axis(mol, 0, std::acos(-1.0) / 2.0,
                              {0.2, 1.3, -3.5});

  assert_vec_close(mol.coordinate(0, 0), {-215.775F, 167.484F, 356.058F},
                   0.001);
}

void test_align_pmi_on_axis_matches_python_alignment_contract() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());

  sasmol::align_pmi_on_axis(mol, 0, 2, sasmol::Axis::z);

  assert(axis_alignment(mol, 0, 2, {0.0, 0.0, 1.0}) > 0.9999);
}

void test_align_pmi_on_cardinal_axes_matches_python_contract() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());

  sasmol::align_pmi_on_cardinal_axes(mol, 0);

  assert(axis_alignment(mol, 0, 2, {0.0, 0.0, 1.0}) > 0.999);
  assert(axis_alignment(mol, 0, 1, {0.0, 1.0, 0.0}) > 0.999);
  assert(axis_alignment(mol, 0, 0, {1.0, 0.0, 0.0}) > 0.999);
}

void test_pmi_aligned_copy_does_not_mutate_source() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());
  const auto original = mol.coordinate(0, 0);

  const auto copy = sasmol::pmi_aligned_on_axis(mol, 0, 2, sasmol::Axis::z);

  assert_vec_close(mol.coordinate(0, 0), original);
  assert(axis_alignment(copy, 0, 2, {0.0, 0.0, 1.0}) > 0.9999);
}

void test_align_full_basis_rotated_fixture() {
  sasmol::PdbReader reader;
  sasmol::Molecule reference;
  sasmol::Molecule moving;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"),
                                reference);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("sasmol/operate", "1CRN-rot.pdb"),
                           moving);
  assert(status.ok());
  const auto indices = all_atom_indices(reference);

  const auto plan = sasmol::initialize_alignment(moving, reference, indices,
                                                 indices, 0);
  sasmol::align(moving, plan, 0);

  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(reference, moving), 0.0,
      0.02);
}

void test_align_full_basis_rotated_shifted_fixture() {
  sasmol::PdbReader reader;
  sasmol::Molecule reference;
  sasmol::Molecule moving;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"),
                                reference);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("sasmol/operate", "1CRN-rot-shift.pdb"),
                           moving);
  assert(status.ok());
  const auto indices = all_atom_indices(reference);

  const auto plan = sasmol::initialize_alignment(moving, reference, indices,
                                                 indices, 0);
  const auto aligned_copy = sasmol::aligned(moving, plan, 0);

  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(reference, aligned_copy), 0.0,
      0.02);
  assert(sasmol::calculate_root_mean_square_deviation(reference, moving) > 1.0);
}

void test_align_ca_subset_moves_whole_molecule_com() {
  sasmol::PdbReader reader;
  sasmol::Molecule reference;
  sasmol::Molecule moving;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"),
                                reference);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), moving);
  assert(status.ok());
  const auto selected = sasmol::select_indices(
      reference, "name[i] == \"CA\" and (resid[i] >= 20 and resid[i] <= 31)");
  assert(selected.ok());
  sasmol::rotate(moving, 0, sasmol::Axis::z, std::acos(-1.0) / 2.0);

  const auto plan = sasmol::initialize_alignment(
      moving, reference, selected.indices, selected.indices, 0);
  sasmol::align(moving, plan, 0);

  auto reference_for_com = reference;
  auto moving_for_com = moving;
  const auto expected_com =
      sasmol::calculate_center_of_mass(reference_for_com, 0);
  const auto result_com = sasmol::calculate_center_of_mass(moving_for_com, 0);
  assert_close_double(result_com.x, expected_com.x, 0.02);
  assert_close_double(result_com.y, expected_com.y, 0.02);
  assert_close_double(result_com.z, expected_com.z, 0.02);
}

void test_align_ca_subset_matches_reference_basis_coordinates() {
  sasmol::PdbReader reader;
  sasmol::Molecule reference;
  sasmol::Molecule moving;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"),
                                reference);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), moving);
  assert(status.ok());
  const auto selected = sasmol::select_indices(
      reference, "name[i] == \"CA\" and (resid[i] >= 20 and resid[i] <= 31)");
  assert(selected.ok());
  sasmol::rotate(moving, 0, sasmol::Axis::z, std::acos(-1.0) / 2.0);
  sasmol::translate(moving, 0, {15.0, -8.0, 3.0});

  const auto plan = sasmol::initialize_alignment(
      moving, reference, selected.indices, selected.indices, 0);
  sasmol::align(moving, plan, 0);

  assert_close_double(selected_rmsd(reference, moving, selected.indices, 0),
                      0.0, 0.02);
}

}  // namespace

int main() {
  test_translate_mutates_only_selected_frame();
  test_translated_returns_copy_without_mutating_source();
  test_center_moves_com_to_origin();
  test_translate_point_moves_com_to_target();
  test_axis_rotation_matches_python_column_vector_convention();
  test_general_axis_rotation_matches_python_row_vector_convention();
  test_euler_rotation_identity();
  test_general_axis_preserves_python_non_unit_axis_behavior();
  test_align_pmi_on_axis_matches_python_alignment_contract();
  test_align_pmi_on_cardinal_axes_matches_python_contract();
  test_pmi_aligned_copy_does_not_mutate_source();
  test_align_full_basis_rotated_fixture();
  test_align_full_basis_rotated_shifted_fixture();
  test_align_ca_subset_moves_whole_molecule_com();
  test_align_ca_subset_matches_reference_basis_coordinates();
  return 0;
}
