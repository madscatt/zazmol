#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-3) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

std::set<std::string> unique_values(const std::vector<std::string>& values) {
  return std::set<std::string>(values.begin(), values.end());
}

void write_moltype_mmcif(const std::filesystem::path& path) {
  const std::vector<std::tuple<int, std::string, std::string>> residues = {
      {1, "ALA", "P"},  {2, "ADE", "D"},  {3, "GUA", "D"},
      {4, "CYT", "D"},  {5, "THY", "D"},  {6, "URA", "R"},
      {7, "A", "R"},    {8, "DA", "D"},   {9, "TIP3", "W"},
      {10, "LIG", "L"}};

  std::ofstream output(path);
  output << "data_moltype_regression\n"
         << "loop_\n"
         << "_atom_site.group_PDB\n"
         << "_atom_site.id\n"
         << "_atom_site.type_symbol\n"
         << "_atom_site.label_atom_id\n"
         << "_atom_site.label_alt_id\n"
         << "_atom_site.label_comp_id\n"
         << "_atom_site.label_asym_id\n"
         << "_atom_site.label_entity_id\n"
         << "_atom_site.label_seq_id\n"
         << "_atom_site.pdbx_PDB_ins_code\n"
         << "_atom_site.Cartn_x\n"
         << "_atom_site.Cartn_y\n"
         << "_atom_site.Cartn_z\n"
         << "_atom_site.occupancy\n"
         << "_atom_site.B_iso_or_equiv\n"
         << "_atom_site.pdbx_formal_charge\n"
         << "_atom_site.auth_seq_id\n"
         << "_atom_site.auth_comp_id\n"
         << "_atom_site.auth_asym_id\n"
         << "_atom_site.auth_atom_id\n"
         << "_atom_site.pdbx_PDB_model_num\n";
  for (const auto& [atom_id, resname, chain] : residues) {
    output << "ATOM " << atom_id << " P P . " << resname << " " << chain
           << " 1 " << atom_id << " ? " << atom_id << ".0 "
           << atom_id + 1 << ".0 " << atom_id + 2
           << ".0 1.00 0.00 ? " << atom_id << " " << resname << " "
           << chain << " P 1\n";
  }
}

void test_1crn_mmcif_matches_existing_pdb_core_fields() {
  sasmol::PdbReader pdb_reader;
  sasmol::MmcifReader mmcif_reader;
  sasmol::Molecule pdb_molecule;
  sasmol::Molecule cif_molecule;

  auto status =
      pdb_reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"),
                          pdb_molecule);
  assert(status.ok());
  status = mmcif_reader.read_mmcif(fixture_path("mmcif_common", "1CRN.cif"),
                                   cif_molecule);
  assert(status.ok());

  assert(cif_molecule.natoms() == pdb_molecule.natoms());
  assert(cif_molecule.number_of_frames() == pdb_molecule.number_of_frames());
  assert(cif_molecule.index() == pdb_molecule.index());
  assert(cif_molecule.original_index() == pdb_molecule.original_index());
  assert(cif_molecule.name() == pdb_molecule.name());
  assert(cif_molecule.resname() == pdb_molecule.resname());
  assert(cif_molecule.chain() == pdb_molecule.chain());
  assert(cif_molecule.resid() == pdb_molecule.resid());
  assert(cif_molecule.original_resid() == pdb_molecule.original_resid());
  assert(cif_molecule.rescode() == pdb_molecule.rescode());
  assert(cif_molecule.segname() == pdb_molecule.segname());
  assert(cif_molecule.element() == pdb_molecule.element());
  for (std::size_t atom = 0; atom < cif_molecule.natoms(); ++atom) {
    const auto cif = cif_molecule.coordinate(0, atom);
    const auto pdb = pdb_molecule.coordinate(0, atom);
    assert_close(cif.x, pdb.x);
    assert_close(cif.y, pdb.y);
    assert_close(cif.z, pdb.z);
  }
}

void test_nucleic_acid_and_heterogen_fixtures_read() {
  sasmol::MmcifReader reader;
  sasmol::Molecule dna;
  sasmol::Molecule rna;
  sasmol::Molecule hemoglobin;

  auto status = reader.read_mmcif(fixture_path("mmcif_common", "1BNA.cif"), dna);
  assert(status.ok());
  status = reader.read_mmcif(fixture_path("mmcif_common", "1EHZ.cif"), rna);
  assert(status.ok());
  status =
      reader.read_mmcif(fixture_path("mmcif_common", "4HHB.cif"), hemoglobin);
  assert(status.ok());

  assert(dna.natoms() == 566);
  assert(unique_values(dna.moltype()).contains("dna"));
  assert(rna.natoms() == 1821);
  assert(unique_values(rna.moltype()).contains("rna"));
  assert(hemoglobin.natoms() == 4779);
  assert(unique_values(hemoglobin.moltype()).contains("protein"));
  assert(unique_values(hemoglobin.resname()).contains("HEM"));
}

void test_read_mmcif_classifies_overlap_resnames_as_nucleic() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_moltype_overlap.cif";
  write_moltype_mmcif(path);

  sasmol::MmcifReader reader;
  sasmol::Molecule molecule;
  const auto status = reader.read_mmcif(path, molecule);

  assert(status.ok());
  assert((molecule.moltype() == std::vector<std::string>{
                                  "protein", "nucleic", "nucleic", "nucleic",
                                  "dna", "rna", "rna", "dna", "water",
                                  "other"}));
  std::filesystem::remove(path);
}

void test_multimodel_nmr_fixture_groups_frames() {
  sasmol::MmcifReader reader;
  sasmol::Molecule molecule;

  const auto status =
      reader.read_mmcif(fixture_path("mmcif_common", "2K39.cif"), molecule);

  assert(status.ok());
  assert(molecule.natoms() == 1231);
  assert(molecule.number_of_frames() == 116);
  assert(molecule.coor().size() == 116 * 1231 * 3);
  assert(molecule.index()[0] == 1);
  assert(molecule.index()[1] == 2);
  assert(molecule.index()[2] == 3);
  assert(molecule.original_index()[0] == 1);
  assert(molecule.original_index()[1228] == 1229);

  const auto first_frame = molecule.coordinate(0, 0);
  const auto second_frame = molecule.coordinate(1, 0);
  assert_close(first_frame.x, 13.434F);
  assert_close(second_frame.x, 13.720F);
}

void test_large_multichain_fixture_reads_full_chain_ids() {
  sasmol::MmcifReader reader;
  sasmol::Molecule molecule;

  const auto status =
      reader.read_mmcif(fixture_path("mmcif_common", "1KP8.cif"), molecule);

  assert(status.ok());
  assert(molecule.natoms() == 57085);
  assert(molecule.number_of_frames() == 1);
  assert(unique_values(molecule.chain()).size() > 1);
  assert(molecule.segname() == molecule.chain());
}

void test_bad_identifier_does_not_mutate_existing_molecule() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_mmcif_id.cif";
  {
    std::ofstream output(path);
    output << "data_bad\n"
           << "loop_\n"
           << "_atom_site.group_PDB\n"
           << "_atom_site.id\n"
           << "_atom_site.type_symbol\n"
           << "_atom_site.label_atom_id\n"
           << "_atom_site.label_alt_id\n"
           << "_atom_site.label_comp_id\n"
           << "_atom_site.label_asym_id\n"
           << "_atom_site.label_seq_id\n"
           << "_atom_site.pdbx_PDB_ins_code\n"
           << "_atom_site.Cartn_x\n"
           << "_atom_site.Cartn_y\n"
           << "_atom_site.Cartn_z\n"
           << "_atom_site.occupancy\n"
           << "_atom_site.B_iso_or_equiv\n"
           << "_atom_site.pdbx_formal_charge\n"
           << "_atom_site.auth_seq_id\n"
           << "_atom_site.auth_comp_id\n"
           << "_atom_site.auth_asym_id\n"
           << "_atom_site.auth_atom_id\n"
           << "_atom_site.pdbx_PDB_model_num\n"
           << "ATOM BAD C CA . ALA A 1 ? 1.0 2.0 3.0 1.00 0.00 ? 1 ALA A CA 1\n";
  }

  sasmol::MmcifReader reader;
  sasmol::Molecule molecule(1, 1);
  molecule.original_index()[0] = 42;
  molecule.name()[0] = "KEEP";
  const auto status = reader.read_mmcif(path, molecule);

  assert(!status.ok());
  assert(status.code == sasmol::IoCode::format_error);
  assert(molecule.original_index()[0] == 42);
  assert(molecule.name()[0] == "KEEP");
  std::filesystem::remove(path);
}

}  // namespace

int main() {
  test_1crn_mmcif_matches_existing_pdb_core_fields();
  test_nucleic_acid_and_heterogen_fixtures_read();
  test_read_mmcif_classifies_overlap_resnames_as_nucleic();
  test_multimodel_nmr_fixture_groups_frames();
  test_large_multichain_fixture_reads_full_chain_ids();
  test_bad_identifier_does_not_mutate_existing_molecule();
  return 0;
}
