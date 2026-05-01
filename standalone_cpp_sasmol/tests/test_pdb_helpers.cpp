#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected) {
  assert(std::fabs(static_cast<double>(actual - expected)) < 1.0e-12);
}

void test_all_zero_axis_guard_only_nudges_first_atom() {
  sasmol::Molecule mol(2, 1);
  mol.set_coordinate(0, 0, {0.0F, 1.0F, 0.0F});
  mol.set_coordinate(0, 1, {0.0F, 2.0F, 0.0F});

  sasmol::PdbWriter writer;
  const auto status = writer.check_for_all_zero_columns(mol);

  assert(status.ok());
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 1.0e-10F);
  assert_close(xyz.y, 1.0F);
  assert_close(xyz.z, 1.0e-10F);

  xyz = mol.coordinate(0, 1);
  assert_close(xyz.x, 0.0F);
  assert_close(xyz.y, 2.0F);
  assert_close(xyz.z, 0.0F);
}

void test_nonzero_axes_are_untouched() {
  sasmol::Molecule mol(2, 1);
  mol.set_coordinate(0, 0, {1.0F, 0.0F, 0.0F});
  mol.set_coordinate(0, 1, {0.0F, 2.0F, 3.0F});

  sasmol::PdbWriter writer;
  const auto status = writer.check_for_all_zero_columns(mol);

  assert(status.ok());
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 1.0F);
  assert_close(xyz.y, 0.0F);
  assert_close(xyz.z, 0.0F);
}

void test_conect_lines_remap_original_to_current_indices() {
  sasmol::Molecule mol(3, 1);
  mol.original_index() = {30, 10, 20};
  mol.index() = {1, 2, 3};
  mol.conect()[0] = {10, 20};
  mol.conect()[1] = {30};

  sasmol::PdbWriter writer;
  const auto lines = writer.create_conect_pdb_lines(mol);

  assert((lines == std::vector<std::string>{
                       "CONECT    1    2    3",
                       "CONECT    2    1",
                   }));
}

}  // namespace

int main() {
  test_all_zero_axis_guard_only_nudges_first_atom();
  test_nonzero_axes_are_untouched();
  test_conect_lines_remap_original_to_current_indices();
  return 0;
}
