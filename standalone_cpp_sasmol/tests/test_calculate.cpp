#include "sasmol/calculate.hpp"
#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-6) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_bounds_close(const sasmol::CoordinateBounds& bounds,
                         sasmol::Vec3 expected_min,
                         sasmol::Vec3 expected_max) {
  assert_close(bounds.minimum.x, expected_min.x);
  assert_close(bounds.minimum.y, expected_min.y);
  assert_close(bounds.minimum.z, expected_min.z);
  assert_close(bounds.maximum.x, expected_max.x);
  assert_close(bounds.maximum.y, expected_max.y);
  assert_close(bounds.maximum.z, expected_max.z);
}

void test_calculate_minimum_and_maximum_all_loaded_frames() {
  sasmol::Molecule mol(2, 2);
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(0, 1, {-4.0F, 5.0F, 6.0F});
  mol.set_coordinate(1, 0, {7.0F, -8.0F, 9.0F});
  mol.set_coordinate(1, 1, {10.0F, 11.0F, -12.0F});

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {-4.0F, -8.0F, -12.0F},
                      {10.0F, 11.0F, 9.0F});
}

void test_calculate_minimum_and_maximum_selected_frames() {
  sasmol::Molecule mol(1, 3);
  mol.set_coordinate(0, 0, {-10.0F, -10.0F, -10.0F});
  mol.set_coordinate(1, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(2, 0, {4.0F, 5.0F, 6.0F});

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol, {1, 2});

  assert_bounds_close(bounds, {1.0F, 2.0F, 3.0F}, {4.0F, 5.0F, 6.0F});
}

void test_calculate_minimum_and_maximum_rejects_empty_molecule() {
  const sasmol::Molecule mol;
  bool threw = false;

  try {
    (void)sasmol::calculate_minimum_and_maximum(mol);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_minimum_and_maximum_rejects_bad_frame() {
  const sasmol::Molecule mol(1, 1);
  bool threw = false;

  try {
    (void)sasmol::calculate_minimum_and_maximum(mol, {1});
  } catch (const std::out_of_range&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_minimum_and_maximum_pdb_fixture_2aad() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {70.721F, 41.799F, 39.354F},
                      {79.712F, 46.273F, 43.910F});
}

void test_calculate_minimum_and_maximum_pdb_fixture_1crn() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {-3.097F, -0.516F, -7.422F},
                      {24.284F, 20.937F, 19.580F});
}

void test_calculate_minimum_and_maximum_all_steps_alias() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1ATM-1to2.pdb"), mol);
  assert(status.ok());

  const auto primary = sasmol::calculate_minimum_and_maximum_all_steps(mol);
  const auto alias = sasmol::calc_minmax_all_steps(mol);

  assert_bounds_close(primary, {73.944F, 38.799F, 41.652F},
                      {76.944F, 41.799F, 41.652F});
  assert_bounds_close(alias, primary.minimum, primary.maximum);
}

}  // namespace

int main() {
  test_calculate_minimum_and_maximum_all_loaded_frames();
  test_calculate_minimum_and_maximum_selected_frames();
  test_calculate_minimum_and_maximum_rejects_empty_molecule();
  test_calculate_minimum_and_maximum_rejects_bad_frame();
  test_calculate_minimum_and_maximum_pdb_fixture_2aad();
  test_calculate_minimum_and_maximum_pdb_fixture_1crn();
  test_calculate_minimum_and_maximum_all_steps_alias();
  return 0;
}
