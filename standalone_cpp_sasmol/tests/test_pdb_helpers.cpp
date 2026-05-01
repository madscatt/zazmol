#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
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

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_scan(const std::filesystem::path& path, std::size_t natoms,
                 std::size_t nframes, sasmol::PdbFrameMode mode) {
  sasmol::PdbReader reader;
  sasmol::PdbFrameScan scan;

  const auto status = reader.scan_pdb_frames(path, scan);

  assert(status.ok());
  assert(scan.natoms == natoms);
  assert(scan.nframes == nframes);
  assert(scan.mode == mode);
}

void assert_scan_fails(const std::filesystem::path& path) {
  sasmol::PdbReader reader;
  sasmol::PdbFrameScan scan;

  const auto status = reader.scan_pdb_frames(path, scan);

  assert(!status.ok());
  assert(status.code == sasmol::IoCode::format_error);
}

void test_pdb_frame_scan_fixtures() {
  assert_scan(fixture_path("pdb_common", "1ATM.pdb"), 1, 1,
              sasmol::PdbFrameMode::end_records);
  assert_scan(fixture_path("pdb_common", "1ATM-1to2.pdb"), 1, 2,
              sasmol::PdbFrameMode::end_records);
  assert_scan(fixture_path("sasmol/file_io", "2AAD-1to3-END.pdb"), 15, 3,
              sasmol::PdbFrameMode::end_records);
  assert_scan(fixture_path("sasmol/file_io", "2AAD-1to3-MODEL.pdb"), 15, 3,
              sasmol::PdbFrameMode::model_records);
  assert_scan(fixture_path("sasmol/file_io", "1AA-NoEND.pdb"), 9, 1,
              sasmol::PdbFrameMode::single);
}

void test_pdb_frame_scan_failure_fixtures() {
  assert_scan_fails(
      fixture_path("sasmol/file_io", "2AAD-1to3-END_wrong_number_atoms.pdb"));
  assert_scan_fails(fixture_path(
      "sasmol/file_io", "2AAD-1to3-MODEL_wrong_number_atoms.pdb"));
  assert_scan_fails(fixture_path(
      "sasmol/file_io", "2AAD-1to3-MODEL_wrongnumber_mix_END.pdb"));
  assert_scan_fails(fixture_path(
      "sasmol/file_io", "2AAD-1to3-MODEL_mix_END_noterminating.pdb"));
  assert_scan_fails(fixture_path("pdb_common", "1PSI.pdb"));
}

}  // namespace

int main() {
  test_all_zero_axis_guard_only_nudges_first_atom();
  test_nonzero_axes_are_untouched();
  test_conect_lines_remap_original_to_current_indices();
  test_pdb_frame_scan_fixtures();
  test_pdb_frame_scan_failure_fixtures();
  return 0;
}
