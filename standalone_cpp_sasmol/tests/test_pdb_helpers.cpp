#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
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

std::string first_coordinate_line(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::string line;
  while (std::getline(input, line)) {
    if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
      return line;
    }
  }
  return {};
}

void test_parse_pdb_atom_record_uses_sasmol_field_names() {
  sasmol::PdbReader reader;
  sasmol::PdbAtomRecord atom;
  const auto line = first_coordinate_line(fixture_path("pdb_common", "2AAD.pdb"));

  const auto status = reader.parse_pdb_atom_record(line, atom);

  assert(status.ok());
  assert(atom.record == "ATOM");
  assert(atom.original_index == 53893);
  assert(atom.name == "N");
  assert(atom.loc == " ");
  assert(atom.resname == "ILE");
  assert(atom.chain == "N");
  assert(atom.resid == 515);
  assert(atom.original_resid == " 515");
  assert(atom.rescode == " ");
  assert_close(atom.coordinate.x, 73.944F);
  assert_close(atom.coordinate.y, 41.799F);
  assert_close(atom.coordinate.z, 41.652F);
  assert(atom.occupancy == "1.00");
  assert(atom.beta == "36.37");
  assert(atom.segname == "N");
  assert(atom.element == "N");
  assert(atom.charge == "  ");
}

void test_parse_pdb_atom_record_preserves_altloc_in_pdbscan_mode() {
  sasmol::PdbReader reader;
  sasmol::PdbReadOptions options;
  options.pdbscan = true;
  sasmol::PdbAtomRecord atom;
  std::string line(80, ' ');
  line.replace(0, 6, "ATOM  ");
  line.replace(6, 5, "    1");
  line.replace(12, 4, " CA ");
  line[16] = 'A';
  line.replace(17, 4, "ALA ");
  line[21] = 'A';
  line.replace(22, 4, "   1");
  line[26] = ' ';
  line.replace(30, 8, "  11.100");
  line.replace(38, 8, "  12.200");
  line.replace(46, 8, "  13.300");
  line.replace(76, 2, " C");

  const auto status = reader.parse_pdb_atom_record(line, atom, options);

  assert(status.ok());
  assert(atom.name == "CA");
  assert(atom.loc == "A");
  assert(atom.resname == "ALA");
  assert(atom.occupancy == "");
  assert(atom.beta == "");
  assert(atom.segname == "");
  assert(atom.element == "C");
}

void test_parse_pdb_atom_record_rejects_bad_required_numbers() {
  sasmol::PdbReader reader;
  sasmol::PdbAtomRecord atom;
  const std::string line =
      "ATOM      X  CA  ALA A   1      11.100  12.200  13.300  1.00  0.00";

  const auto status = reader.parse_pdb_atom_record(line, atom);

  assert(status.code == sasmol::IoCode::format_error);
}

void test_read_pdb_single_frame_1atm() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;

  const auto status = reader.read_pdb(fixture_path("pdb_common", "1ATM.pdb"), mol);

  assert(status.ok());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 1);
  assert(mol.record()[0] == "ATOM");
  assert(mol.original_index()[0] == 53893);
  assert(mol.index()[0] == 1);
  assert(mol.name()[0] == "N");
  assert(mol.resname()[0] == "ILE");
  assert(mol.chain()[0] == "N");
  assert(mol.resid()[0] == 515);
  const auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);
}

void test_read_pdb_single_frame_2aad_descriptors() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;

  const auto status = reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);

  assert(status.ok());
  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 1);
  assert(mol.name()[0] == "N");
  assert(mol.resname()[0] == "ILE");
  assert(mol.chain()[0] == "N");
  assert(mol.resid()[0] == 515);
  assert(mol.occupancy()[0] == "1.00");
  assert(mol.beta()[0] == "36.37");
  assert(mol.segname()[0] == "N");
  assert(mol.element()[0] == "N");
  assert(mol.charge()[0] == "  ");
  assert(mol.name()[14] == "CG2");
  const auto xyz = mol.coordinate(0, 14);
  assert_close(xyz.x, 76.970F);
  assert_close(xyz.y, 46.273F);
  assert_close(xyz.z, 42.000F);
}

void test_read_pdb_multi_frame_is_explicitly_deferred() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;

  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1ATM-1to2.pdb"), mol);

  assert(status.ok());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 2);
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);
  xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, 73.944F);
  assert_close(xyz.y, 38.799F);
  assert_close(xyz.z, 41.652F);
}

void test_read_pdb_end_separated_multiframe_coordinates() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;

  const auto status =
      reader.read_pdb(fixture_path("sasmol/file_io", "2AAD-1to3-END.pdb"), mol);

  assert(status.ok());
  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 3);
  assert(mol.name()[0] == "N");
  auto xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, -73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);
  xyz = mol.coordinate(2, 14);
  assert_close(xyz.x, 76.970F);
  assert_close(xyz.y, -46.273F);
  assert_close(xyz.z, 42.000F);
}

void test_read_pdb_model_multiframe_coordinates() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;

  const auto status =
      reader.read_pdb(fixture_path("sasmol/file_io", "2AAD-1to3-MODEL.pdb"), mol);

  assert(status.ok());
  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 3);
  auto xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, -73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);
}

}  // namespace

int main() {
  test_all_zero_axis_guard_only_nudges_first_atom();
  test_nonzero_axes_are_untouched();
  test_conect_lines_remap_original_to_current_indices();
  test_pdb_frame_scan_fixtures();
  test_pdb_frame_scan_failure_fixtures();
  test_parse_pdb_atom_record_uses_sasmol_field_names();
  test_parse_pdb_atom_record_preserves_altloc_in_pdbscan_mode();
  test_parse_pdb_atom_record_rejects_bad_required_numbers();
  test_read_pdb_single_frame_1atm();
  test_read_pdb_single_frame_2aad_descriptors();
  test_read_pdb_multi_frame_is_explicitly_deferred();
  test_read_pdb_end_separated_multiframe_coordinates();
  test_read_pdb_model_multiframe_coordinates();
  return 0;
}
