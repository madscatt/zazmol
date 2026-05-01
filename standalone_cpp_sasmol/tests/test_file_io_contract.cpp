#include "sasmol/file_io.hpp"

#include <cassert>
#include <filesystem>
#include <string>

namespace {

void test_pdb_contract_defaults_are_tolerant() {
  const sasmol::PdbReadOptions options;

  assert(options.tolerant);
  assert(options.resolve_elements);
  assert(options.preserve_conect);
  assert(options.apply_all_zero_coordinate_guard);
  assert(sasmol::PdbReader::tolerant_by_default());
}

void test_pdb_reader_reports_unsupported_for_multiframe_parser_slice() {
  sasmol::Molecule mol;
  sasmol::PdbReader reader;
  const auto fixture =
      std::filesystem::path(SASMOL_TEST_DATA_DIR) / "pdb_common" / "1ATM-1to2.pdb";

  const auto status = reader.read_pdb(fixture, mol);

  assert(!status.ok());
  assert(status.code == sasmol::IoCode::unsupported);
}

void test_dcd_contract_defaults_are_sequential() {
  const sasmol::DcdReadOptions options;

  assert(options.sequential);
  assert(!options.reopen_for_random_frame);
  assert(sasmol::DcdReader::sequential_by_default());
}

void test_dcd_reader_has_explicit_open_close_state() {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  const std::filesystem::path fixture =
      std::filesystem::path(SASMOL_TEST_DATA_DIR) / "dcd_common" / "1ATM.dcd";

  assert(!reader.is_open());
  auto status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::not_open);

  status = reader.open_dcd_read(fixture);
  assert(reader.is_open());
  assert(status.ok());

  status = reader.close_dcd_read();
  assert(status.ok());
  assert(!reader.is_open());
}

void test_dcd_random_frame_access_is_explicit_reopen_scan_contract() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;
  sasmol::DcdReadOptions options;
  const std::filesystem::path fixture =
      std::filesystem::path(SASMOL_TEST_DATA_DIR) / "dcd_common" / "1ATM.dcd";
  options.reopen_for_random_frame = true;

  const auto status = reader.read_single_dcd_step(fixture, 2, mol, options);

  assert(status.ok());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 1);
}

void test_dcd_writer_has_explicit_open_close_state() {
  sasmol::DcdWriter writer;
  sasmol::Molecule mol(1, 1);
  const auto output =
      std::filesystem::temp_directory_path() / "sasmol_contract_writer.dcd";

  auto status = writer.write_dcd_header(mol, 1);
  assert(status.code == sasmol::IoCode::not_open);

  status = writer.open_dcd_write(output);
  assert(writer.is_open());
  assert(status.ok());

  status = writer.close_dcd_write();
  assert(status.ok());
  assert(!writer.is_open());
  std::filesystem::remove(output);
}

}  // namespace

int main() {
  test_pdb_contract_defaults_are_tolerant();
  test_pdb_reader_reports_unsupported_for_multiframe_parser_slice();
  test_dcd_contract_defaults_are_sequential();
  test_dcd_reader_has_explicit_open_close_state();
  test_dcd_random_frame_access_is_explicit_reopen_scan_contract();
  test_dcd_writer_has_explicit_open_close_state();
  return 0;
}
