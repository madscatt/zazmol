#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>

namespace {

std::filesystem::path fixture_path(const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / "dcd_common" / name;
}

void assert_header(const char* filename, std::size_t natoms,
                   std::size_t nframes, int charmm_flags) {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(fixture_path(filename));
  assert(status.ok());
  assert(reader.is_open());

  status = reader.read_header(header);
  assert(status.ok());
  assert(header.natoms == natoms);
  assert(header.nframes == nframes);
  assert(header.reverse_endian == false);
  assert(header.charmm_format == true);
  assert(header.has_unit_cell == true);
  assert(header.charmm_flags == charmm_flags);
  assert(header.namnf == 0);
  assert(header.nsavc == 1);

  status = reader.close_dcd_read();
  assert(status.ok());
}

void test_small_fixture_headers() {
  assert_header("1ATM.dcd", 1, 2, 5);
  assert_header("2AAD.dcd", 15, 3, 5);
  assert_header("rna-1to10.dcd", 10632, 10, 5);
}

void test_missing_file_is_file_error() {
  sasmol::DcdReader reader;
  const auto status = reader.open_dcd_read(fixture_path("missing.dcd"));

  assert(status.code == sasmol::IoCode::file_error);
  assert(!reader.is_open());
}

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 0.001) {
  assert(std::fabs(static_cast<double>(actual - expected)) <= tolerance);
}

void test_sequential_frame_reads_1atm() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());

  status = reader.read_next_frame(mol);
  assert(status.ok());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 2);
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);

  status = reader.read_next_frame(mol);
  assert(status.ok());
  xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, 73.944F);
  assert_close(xyz.y, 38.799F);
  assert_close(xyz.z, 41.652F);

  status = reader.read_next_frame(mol);
  assert(status.code == sasmol::IoCode::end_of_file);
}

void test_sequential_frame_reads_2aad_middle_and_final() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("2AAD.dcd"));
  assert(status.ok());

  status = reader.read_next_frame(mol);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());

  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 3);

  auto xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, -73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);

  xyz = mol.coordinate(2, 14);
  assert_close(xyz.x, 76.970F);
  assert_close(xyz.y, -46.273F);
  assert_close(xyz.z, 42.000F);
}

void test_sequential_frame_reads_rna_final_sample() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("rna-1to10.dcd"));
  assert(status.ok());

  for (int i = 0; i < 10; ++i) {
    status = reader.read_next_frame(mol);
    assert(status.ok());
  }

  assert(mol.natoms() == 10632);
  assert(mol.number_of_frames() == 10);

  const auto xyz = mol.coordinate(9, 10631);
  assert_close(xyz.x, -6.392F);
  assert_close(xyz.y, 14.348F);
  assert_close(xyz.z, 20.914F);
}

}  // namespace

int main() {
  test_small_fixture_headers();
  test_missing_file_is_file_error();
  test_sequential_frame_reads_1atm();
  test_sequential_frame_reads_2aad_middle_and_final();
  test_sequential_frame_reads_rna_final_sample();
  return 0;
}
