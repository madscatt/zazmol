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

void test_pdb_reader_reports_deferred_parser_without_mutating_molecule() {
  sasmol::Molecule mol(2, 1);
  sasmol::PdbReader reader;

  const auto status = reader.read_pdb("fixture.pdb", mol);

  assert(!status.ok());
  assert(status.code == sasmol::IoCode::not_implemented);
  assert(mol.natoms() == 2);
  assert(mol.number_of_frames() == 1);
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
  sasmol::Molecule mol(1, 1);
  sasmol::DcdReadOptions options;
  options.reopen_for_random_frame = true;

  const auto status =
      reader.read_single_dcd_step("fixture.dcd", 0, mol, options);

  assert(status.code == sasmol::IoCode::not_implemented);
  assert(status.message.find("reopen") != std::string::npos);
}

void test_dcd_writer_has_explicit_open_close_state() {
  sasmol::DcdWriter writer;
  sasmol::Molecule mol(1, 1);

  auto status = writer.write_dcd_header(mol, 1);
  assert(status.code == sasmol::IoCode::not_open);

  status = writer.open_dcd_write("fixture.dcd");
  assert(writer.is_open());
  assert(status.code == sasmol::IoCode::not_implemented);

  status = writer.close_dcd_write();
  assert(status.ok());
  assert(!writer.is_open());
}

}  // namespace

int main() {
  test_pdb_contract_defaults_are_tolerant();
  test_pdb_reader_reports_deferred_parser_without_mutating_molecule();
  test_dcd_contract_defaults_are_sequential();
  test_dcd_reader_has_explicit_open_close_state();
  test_dcd_random_frame_access_is_explicit_reopen_scan_contract();
  test_dcd_writer_has_explicit_open_close_state();
  return 0;
}
