#include "sasmol/topology.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>

namespace {

std::filesystem::path topology_fixture(const char* filename) {
  return std::filesystem::path(SASMOL_TOPOLOGY_TEST_DATA_DIR) / filename;
}

void assert_close(sasmol::calc_type actual, sasmol::calc_type expected) {
  assert(std::fabs(actual - expected) < 1.0e-12);
}

void test_assign_charmm_types_sets_atom_aligned_values() {
  sasmol::Molecule mol(3, 1);

  const auto result = sasmol::assign_charmm_types(mol, {"NH1", "CT1", "C"});

  assert(result.ok());
  assert((mol.charmm_type() == std::vector<std::string>{"NH1", "CT1", "C"}));
}

void test_assign_charmm_types_rejects_length_mismatch_without_mutation() {
  sasmol::Molecule mol(2, 1);
  mol.charmm_type() = {"OLD", "TYPE"};

  const auto result = sasmol::assign_charmm_types(mol, {"NH1"});

  assert(!result.ok());
  assert((mol.charmm_type() == std::vector<std::string>{"OLD", "TYPE"}));
}

void test_assign_charmm_types_allows_empty_molecule_empty_types() {
  sasmol::Molecule mol;

  const auto result = sasmol::assign_charmm_types(mol, {});

  assert(result.ok());
  assert(mol.charmm_type().empty());
}

void test_assign_atom_charges_sets_atom_aligned_values() {
  sasmol::Molecule mol(3, 1);

  const auto result = sasmol::assign_atom_charges(mol, {0.1, -0.2, 0.3});

  assert(result.ok());
  assert(mol.atom_charge().size() == 3);
  assert_close(mol.atom_charge()[0], 0.1);
  assert_close(mol.atom_charge()[1], -0.2);
  assert_close(mol.atom_charge()[2], 0.3);
}

void test_assign_atom_charges_rejects_length_mismatch_without_mutation() {
  sasmol::Molecule mol(2, 1);
  mol.atom_charge() = {0.4, -0.4};

  const auto result = sasmol::assign_atom_charges(mol, {0.1});

  assert(!result.ok());
  assert(mol.atom_charge().size() == 2);
  assert_close(mol.atom_charge()[0], 0.4);
  assert_close(mol.atom_charge()[1], -0.4);
}

void test_assign_charmm_types_from_atom_table_validates_names() {
  sasmol::Molecule mol(3, 1);
  mol.name() = {"N", "CA", "C"};

  const auto result = sasmol::assign_charmm_types_from_atom_table(
      mol, {{"N", "NH1"}, {"CA", "CT1"}, {"C", "C"}});

  assert(result.ok());
  assert((mol.charmm_type() == std::vector<std::string>{"NH1", "CT1", "C"}));
}

void test_assign_charmm_types_from_atom_table_rejects_name_mismatch() {
  sasmol::Molecule mol(3, 1);
  mol.name() = {"N", "CA", "C"};
  mol.charmm_type() = {"OLD", "TYPE", "VALUES"};

  const auto result = sasmol::assign_charmm_types_from_atom_table(
      mol, {{"N", "NH1"}, {"CB", "CT1"}, {"C", "C"}});

  assert(!result.ok());
  assert((mol.charmm_type() ==
          std::vector<std::string>{"OLD", "TYPE", "VALUES"}));
}

void test_assign_charmm_types_from_atom_table_rejects_name_length_mismatch() {
  sasmol::Molecule mol(3, 1);
  mol.name().clear();
  mol.charmm_type() = {"OLD", "TYPE", "VALUES"};

  const auto result = sasmol::assign_charmm_types_from_atom_table(
      mol, {{"N", "NH1"}, {"CA", "CT1"}, {"C", "C"}});

  assert(!result.ok());
  assert((mol.charmm_type() ==
          std::vector<std::string>{"OLD", "TYPE", "VALUES"}));
}

void test_assign_charmm_types_and_atom_charges_from_atom_table() {
  sasmol::Molecule mol(3, 1);
  mol.name() = {"N", "CA", "C"};

  const auto result =
      sasmol::assign_charmm_types_and_atom_charges_from_atom_table(
          mol, {{"N", "NH1", -0.47}, {"CA", "CT1", 0.07}, {"C", "C", 0.51}});

  assert(result.ok());
  assert((mol.charmm_type() == std::vector<std::string>{"NH1", "CT1", "C"}));
  assert(mol.atom_charge().size() == 3);
  assert_close(mol.atom_charge()[0], -0.47);
  assert_close(mol.atom_charge()[1], 0.07);
  assert_close(mol.atom_charge()[2], 0.51);
}

void test_assign_charmm_types_and_atom_charges_rejects_name_mismatch() {
  sasmol::Molecule mol(3, 1);
  mol.name() = {"N", "CA", "C"};
  mol.charmm_type() = {"OLD", "TYPE", "VALUES"};
  mol.atom_charge() = {0.1, 0.2, 0.3};

  const auto result =
      sasmol::assign_charmm_types_and_atom_charges_from_atom_table(
          mol, {{"N", "NH1", -0.47}, {"CB", "CT1", 0.07}, {"C", "C", 0.51}});

  assert(!result.ok());
  assert((mol.charmm_type() ==
          std::vector<std::string>{"OLD", "TYPE", "VALUES"}));
  assert_close(mol.atom_charge()[0], 0.1);
  assert_close(mol.atom_charge()[1], 0.2);
  assert_close(mol.atom_charge()[2], 0.3);
}

void test_assign_charmm_types_and_atom_charges_rejects_length_mismatch() {
  sasmol::Molecule mol(3, 1);
  mol.name() = {"N", "CA", "C"};
  mol.charmm_type() = {"OLD", "TYPE", "VALUES"};
  mol.atom_charge() = {0.1, 0.2, 0.3};

  const auto result =
      sasmol::assign_charmm_types_and_atom_charges_from_atom_table(
          mol, {{"N", "NH1", -0.47}, {"CA", "CT1", 0.07}});

  assert(!result.ok());
  assert((mol.charmm_type() ==
          std::vector<std::string>{"OLD", "TYPE", "VALUES"}));
  assert_close(mol.atom_charge()[0], 0.1);
  assert_close(mol.atom_charge()[1], 0.2);
  assert_close(mol.atom_charge()[2], 0.3);
}

sasmol::CharmmResidueDefinition glycine_definition() {
  return {.resname = "GLY",
          .atoms = {{"N", "NH1", -0.47},
                    {"CA", "CT2", -0.02},
                    {"C", "C", 0.51},
                    {"O", "O", -0.51}}};
}

void test_validate_charmm_residue_atoms_accepts_exact_match_any_order() {
  const auto validation =
      sasmol::validate_charmm_residue_atoms({"CA", "O", "N", "C"},
                                            glycine_definition());

  assert(validation.ok());
  assert(validation.missing_atoms.empty());
  assert(validation.extra_atoms.empty());
  assert(validation.duplicate_molecule_atoms.empty());
  assert(validation.duplicate_topology_atoms.empty());
}

void test_validate_charmm_residue_atoms_reports_missing_atom() {
  const auto validation =
      sasmol::validate_charmm_residue_atoms({"N", "CA", "C"},
                                            glycine_definition());

  assert(!validation.ok());
  assert((validation.missing_atoms == std::vector<std::string>{"O"}));
  assert(validation.extra_atoms.empty());
}

void test_validate_charmm_residue_atoms_reports_extra_atom() {
  const auto validation =
      sasmol::validate_charmm_residue_atoms({"N", "CA", "C", "O", "CB"},
                                            glycine_definition());

  assert(!validation.ok());
  assert(validation.missing_atoms.empty());
  assert((validation.extra_atoms == std::vector<std::string>{"CB"}));
}

void test_validate_charmm_residue_atoms_reports_duplicate_molecule_atom() {
  const auto validation =
      sasmol::validate_charmm_residue_atoms({"N", "CA", "C", "O", "CA"},
                                            glycine_definition());

  assert(!validation.ok());
  assert((validation.duplicate_molecule_atoms ==
          std::vector<std::string>{"CA"}));
}

void test_validate_charmm_residue_atoms_reports_duplicate_topology_atom() {
  sasmol::CharmmResidueDefinition residue = glycine_definition();
  residue.atoms.push_back({"CA", "CT2", -0.02});

  const auto validation =
      sasmol::validate_charmm_residue_atoms({"N", "CA", "C", "O"}, residue);

  assert(!validation.ok());
  assert((validation.duplicate_topology_atoms ==
          std::vector<std::string>{"CA"}));
}

void test_parse_charmm_topology_globals_matches_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology_globals(
      topology_fixture("minimal_mass_only.rtf"));

  assert(result.ok());
  assert(result.topology.masses.size() == 2);
  assert(result.topology.masses[0].index == "1");
  assert(result.topology.masses[0].atom_type == "H");
  assert(result.topology.masses[0].mass == "1.00800");
  assert(result.topology.masses[1].index == "2");
  assert(result.topology.masses[1].atom_type == "C");
  assert(result.topology.masses[1].mass == "12.01100");
  assert((result.topology.declarations == std::vector<std::string>{"-C"}));
  assert((result.topology.defaults ==
          std::vector<std::string>{"FIRS", "NTER", "LAST", "CTER"}));
  assert((result.topology.auto_terms ==
          std::vector<std::string>{"ANGLES", "DIHE"}));
}

void test_parse_charmm_topology_globals_reports_malformed_records() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_global_topology.rtf";
  {
    std::ofstream output(path);
    output << "MASS 1 H\n";
    output << "DECL\n";
    output << "DEFA FIRS NTER LAST\n";
    output << "AUTO ANGLES\n";
  }

  const auto result = sasmol::parse_charmm_topology_globals(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 4);
  assert(result.topology.masses.empty());
  assert(result.topology.declarations.empty());
  assert(result.topology.defaults.empty());
  assert(result.topology.auto_terms.empty());
}

void test_parse_charmm_topology_residue_atoms_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_resi_atoms.rtf"));

  assert(result.ok());
  assert(result.topology.masses.size() == 4);
  assert(result.topology.entries.size() == 1);

  const auto& gly = result.topology.entries[0];
  assert(gly.kind == sasmol::CharmmTopologyEntryKind::Residue);
  assert(gly.name == "GLY");
  assert(gly.total_charge == "0.00");
  assert(gly.atoms.size() == 4);
  assert(gly.atoms[0].name == "N");
  assert(gly.atoms[0].charmm_type == "NH1");
  assert(gly.atoms[0].atom_charge == "-0.47");
  assert(gly.atoms[1].name == "CA");
  assert(gly.atoms[1].charmm_type == "CT2");
  assert(gly.atoms[1].atom_charge == "-0.02");
  assert(gly.atoms[2].name == "C");
  assert(gly.atoms[2].charmm_type == "C");
  assert(gly.atoms[2].atom_charge == "0.51");
  assert(gly.atoms[3].name == "O");
  assert(gly.atoms[3].charmm_type == "O");
  assert(gly.atoms[3].atom_charge == "-0.51");
}

void test_parse_charmm_topology_patch_atoms_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_pres_atoms_dele.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 2);

  const auto& gly = result.topology.entries[0];
  assert(gly.kind == sasmol::CharmmTopologyEntryKind::Residue);
  assert(gly.name == "GLY");
  assert(gly.total_charge == "0.00");
  assert(gly.atoms.size() == 4);

  const auto& nter = result.topology.entries[1];
  assert(nter.kind == sasmol::CharmmTopologyEntryKind::Patch);
  assert(nter.name == "NTER");
  assert(nter.total_charge == "1.00");
  assert(nter.atoms.size() == 4);
  assert(nter.atoms[0].name == "N");
  assert(nter.atoms[0].charmm_type == "NH3");
  assert(nter.atoms[0].atom_charge == "-0.30");
  assert(nter.atoms[1].name == "HT1");
  assert(nter.atoms[1].charmm_type == "HC");
  assert(nter.atoms[1].atom_charge == "0.33");
  assert(nter.atoms[2].name == "HT2");
  assert(nter.atoms[2].charmm_type == "HC");
  assert(nter.atoms[2].atom_charge == "0.33");
  assert(nter.atoms[3].name == "HT3");
  assert(nter.atoms[3].charmm_type == "HC");
  assert(nter.atoms[3].atom_charge == "0.33");
}

void test_parse_charmm_topology_reports_atom_before_residue_or_patch() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_atom_topology.rtf";
  {
    std::ofstream output(path);
    output << "ATOM N NH1 -0.47\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.empty());
}

}  // namespace

int main() {
  test_assign_charmm_types_sets_atom_aligned_values();
  test_assign_charmm_types_rejects_length_mismatch_without_mutation();
  test_assign_charmm_types_allows_empty_molecule_empty_types();
  test_assign_atom_charges_sets_atom_aligned_values();
  test_assign_atom_charges_rejects_length_mismatch_without_mutation();
  test_assign_charmm_types_from_atom_table_validates_names();
  test_assign_charmm_types_from_atom_table_rejects_name_mismatch();
  test_assign_charmm_types_from_atom_table_rejects_name_length_mismatch();
  test_assign_charmm_types_and_atom_charges_from_atom_table();
  test_assign_charmm_types_and_atom_charges_rejects_name_mismatch();
  test_assign_charmm_types_and_atom_charges_rejects_length_mismatch();
  test_validate_charmm_residue_atoms_accepts_exact_match_any_order();
  test_validate_charmm_residue_atoms_reports_missing_atom();
  test_validate_charmm_residue_atoms_reports_extra_atom();
  test_validate_charmm_residue_atoms_reports_duplicate_molecule_atom();
  test_validate_charmm_residue_atoms_reports_duplicate_topology_atom();
  test_parse_charmm_topology_globals_matches_python_oracle_fixture();
  test_parse_charmm_topology_globals_reports_malformed_records();
  test_parse_charmm_topology_residue_atoms_match_python_oracle_fixture();
  test_parse_charmm_topology_patch_atoms_match_python_oracle_fixture();
  test_parse_charmm_topology_reports_atom_before_residue_or_patch();
  return 0;
}
