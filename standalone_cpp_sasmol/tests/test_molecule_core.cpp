#include "sasmol/molecule.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected) {
  assert(std::fabs(static_cast<double>(actual - expected)) < 1.0e-6);
}

void test_defaults_and_integrity() {
  sasmol::Molecule mol(3, 2);

  assert(mol.natoms() == 3);
  assert(mol.number_of_frames() == 2);
  assert(mol.coor().size() == 18);
  assert(mol.record()[0] == "ATOM");
  assert(mol.index()[2] == 3);
  assert(mol.original_index()[1] == 2);
  assert(mol.occupancy()[0] == "1.00");
  assert(mol.beta()[0] == "0.00");

  const auto report = mol.check_integrity();
  assert(report.ok());
  assert(report.lengths.at("name") == 3);
  assert(report.lengths.at("coor") == 18);
}

void test_coordinates_are_frame_major_xyz_triplets() {
  sasmol::Molecule mol(2, 2);

  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(0, 1, {4.0F, 5.0F, 6.0F});
  mol.set_coordinate(1, 0, {7.0F, 8.0F, 9.0F});

  assert_close(mol.coor()[0], 1.0F);
  assert_close(mol.coor()[1], 2.0F);
  assert_close(mol.coor()[2], 3.0F);
  assert_close(mol.coor()[3], 4.0F);
  assert_close(mol.coor()[6], 7.0F);

  const auto value = mol.coordinate(1, 0);
  assert_close(value.x, 7.0F);
  assert_close(value.y, 8.0F);
  assert_close(value.z, 9.0F);
}

void test_coordinate_views_mutate_selected_frame() {
  sasmol::Molecule mol(2, 2);
  auto frame_one = mol.coordinate_view(1);

  frame_one.set(1, {10.0F, 11.0F, 12.0F});

  const auto frame_zero_value = mol.coordinate(0, 1);
  assert_close(frame_zero_value.x, 0.0F);
  assert_close(frame_zero_value.y, 0.0F);
  assert_close(frame_zero_value.z, 0.0F);

  const auto frame_one_value = mol.coordinate(1, 1);
  assert_close(frame_one_value.x, 10.0F);
  assert_close(frame_one_value.y, 11.0F);
  assert_close(frame_one_value.z, 12.0F);
}

void test_integrity_reports_descriptor_mismatch() {
  sasmol::Molecule mol(2, 1);
  mol.name().pop_back();

  const auto report = mol.check_integrity();

  assert(!report.ok());
  assert(report.lengths.at("name") == 1);
  assert(report.issues.size() == 1);
  assert(report.issues[0].field == "name");
  assert(report.issues[0].expected == 2);
  assert(report.issues[0].actual == 1);
}

void test_out_of_range_coordinates_throw() {
  sasmol::Molecule mol(1, 1);
  bool threw = false;

  try {
    (void)mol.coordinate(0, 1);
  } catch (const std::out_of_range&) {
    threw = true;
  }

  assert(threw);
}

}  // namespace

int main() {
  test_defaults_and_integrity();
  test_coordinates_are_frame_major_xyz_triplets();
  test_coordinate_views_mutate_selected_frame();
  test_integrity_reports_descriptor_mismatch();
  test_out_of_range_coordinates_throw();
  return 0;
}
