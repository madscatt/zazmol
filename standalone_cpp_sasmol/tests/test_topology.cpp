#include "sasmol/topology.hpp"

#include "sasmol/file_io.hpp"

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

sasmol::Molecule fasta_test_molecule() {
  sasmol::Molecule molecule(7, 1);
  molecule.record() = {"ATOM", "ATOM", "ATOM", "ATOM", "ATOM", "ATOM",
                       "HETATM"};
  molecule.resname() = {"THR", "THR", "CYS", "CYS", "PRO", "ADE", "HOH"};
  molecule.resid() = {1, 1, 2, 2, 3, 4, 5};
  molecule.chain() = {"A", "A", "A", "A", "B", "B", "W"};
  molecule.segname() = {"SEG1", "SEG1", "SEG1", "SEG1", "SEG2", "SEG2",
                        "WAT"};
  return molecule;
}

sasmol::Molecule constraint_test_molecule() {
  sasmol::Molecule molecule(4, 1);
  molecule.record() = {"ATOM", "ATOM", "ATOM", "HETATM"};
  molecule.index() = {1, 2, 3, 4};
  molecule.name() = {"N", "H", "P", "O"};
  molecule.loc() = {" ", " ", " ", " "};
  molecule.resname() = {"ALA", "ALA", "ADE", "HOH"};
  molecule.chain() = {"A", "A", "B", "W"};
  molecule.resid() = {1, 1, 2, 3};
  molecule.rescode() = {" ", " ", " ", " "};
  molecule.occupancy() = {"0.25", "0.25", "0.25", "0.25"};
  molecule.beta() = {"0.50", "0.50", "0.50", "0.50"};
  molecule.segname() = {"SEG1", "SEG1", "SEG2", "WAT"};
  molecule.element() = {"N", "H", "P", "O"};
  molecule.charge() = {" ", " ", " ", " "};
  molecule.moltype() = {"protein", "protein", "nucleic", "water"};
  molecule.set_coordinate(0, 0, {1.0F, 0.0F, 0.0F});
  molecule.set_coordinate(0, 1, {0.0F, 1.0F, 0.0F});
  molecule.set_coordinate(0, 2, {0.0F, 0.0F, 1.0F});
  molecule.set_coordinate(0, 3, {1.0F, 1.0F, 1.0F});
  return molecule;
}

void test_create_fasta_default_sequence_and_in_place_string() {
  auto molecule = fasta_test_molecule();

  const auto result = sasmol::create_fasta(molecule);

  assert(result.ok());
  assert((result.sequence == std::vector<std::string>{"T", "C", "P", "A",
                                                      "X"}));
  assert(result.formatted.empty());

  const auto set_result = sasmol::create_fasta_in_place(molecule);
  assert(set_result.ok());
  assert(molecule.fasta() == "TCPAX");
}

void test_create_fasta_formatted_width_and_name() {
  auto molecule = fasta_test_molecule();
  sasmol::FastaOptions options;
  options.fasta_format = true;
  options.name = "demo";
  options.width = 3;

  const auto result = sasmol::create_fasta(molecule, options);

  assert(result.ok());
  assert(result.formatted == ">demo\nTCP\nAX\n");
}

void test_create_fasta_formatted_by_chain_and_segname() {
  auto molecule = fasta_test_molecule();
  sasmol::FastaOptions options;
  options.fasta_format = true;
  options.name = "demo";
  options.split_mode = sasmol::FastaSplitMode::by_chain;

  auto result = sasmol::create_fasta(molecule, options);

  assert(result.ok());
  assert(result.formatted.find(">demo chain:A\nTC\n") != std::string::npos);
  assert(result.formatted.find(">demo chain:B\nPA\n") != std::string::npos);
  assert(result.formatted.find(">demo chain:W\nX\n") != std::string::npos);

  options.split_mode = sasmol::FastaSplitMode::by_segname;
  result = sasmol::create_fasta(molecule, options);

  assert(result.ok());
  assert(result.formatted.find(">demo segname:SEG1\nTC\n") !=
         std::string::npos);
  assert(result.formatted.find(">demo segname:SEG2\nPA\n") !=
         std::string::npos);
  assert(result.formatted.find(">demo segname:WAT\nX\n") !=
         std::string::npos);
}

void test_create_fasta_excludes_hetatm_when_requested() {
  auto molecule = fasta_test_molecule();
  sasmol::FastaOptions options;
  options.fasta_format = true;
  options.exclude_hetatm = true;

  const auto result = sasmol::create_fasta(molecule, options);

  assert(result.ok());
  assert(result.sequence.back() == "X");
  assert(result.formatted == ">\nTCPA\n");
}

void test_create_fasta_rejects_unknown_residue_without_mutation() {
  auto molecule = fasta_test_molecule();
  molecule.resname()[0] = "BOG";
  molecule.fasta() = "old";

  const auto result = sasmol::create_fasta_in_place(molecule);

  assert(!result.ok());
  assert(molecule.fasta() == "old");
}

void test_renumber_default_updates_index_and_resid() {
  sasmol::Molecule molecule(5, 1);
  molecule.index() = {10, 11, 12, 13, 14};
  molecule.resid() = {20, 20, 25, 25, 99};

  const auto result = sasmol::renumber(molecule);

  assert(result.ok());
  assert((molecule.index() == std::vector<int>{1, 2, 3, 4, 5}));
  assert((molecule.resid() == std::vector<int>{1, 1, 2, 2, 3}));
}

void test_renumber_index_only_preserves_resid() {
  sasmol::Molecule molecule(3, 1);
  molecule.index() = {1, 2, 3};
  molecule.resid() = {7, 7, 8};
  sasmol::RenumberOptions options;
  options.index_start = 3;

  const auto result = sasmol::renumber(molecule, options);

  assert(result.ok());
  assert((molecule.index() == std::vector<int>{3, 4, 5}));
  assert((molecule.resid() == std::vector<int>{7, 7, 8}));
}

void test_renumber_resid_only_preserves_index() {
  sasmol::Molecule molecule(4, 1);
  molecule.index() = {4, 5, 6, 7};
  molecule.resid() = {100, 100, 101, 105};
  sasmol::RenumberOptions options;
  options.resid_start = 8;

  const auto result = sasmol::renumber(molecule, options);

  assert(result.ok());
  assert((molecule.index() == std::vector<int>{4, 5, 6, 7}));
  assert((molecule.resid() == std::vector<int>{8, 8, 9, 10}));
}

void test_renumber_index_and_resid_custom_starts() {
  sasmol::Molecule molecule(3, 1);
  molecule.index() = {1, 2, 3};
  molecule.resid() = {2, 3, 3};
  sasmol::RenumberOptions options;
  options.index_start = 223;
  options.resid_start = 18;

  const auto result = sasmol::renumber(molecule, options);

  assert(result.ok());
  assert((molecule.index() == std::vector<int>{223, 224, 225}));
  assert((molecule.resid() == std::vector<int>{18, 19, 19}));
}

void test_renumber_rejects_descriptor_mismatch_without_mutation() {
  sasmol::Molecule molecule(3, 1);
  molecule.index() = {1, 2, 3};
  molecule.resid() = {7, 8};

  const auto result = sasmol::renumber(molecule);

  assert(!result.ok());
  assert((molecule.index() == std::vector<int>{1, 2, 3}));
  assert((molecule.resid() == std::vector<int>{7, 8}));
}

void test_apply_constraint_descriptor_sets_heavy_beta() {
  auto molecule = constraint_test_molecule();

  const auto result =
      sasmol::apply_constraint_descriptor(molecule, sasmol::ConstraintBasis::heavy);

  assert(result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"1.00", "0.00", "1.00", "1.00"}));
}

void test_apply_constraint_descriptor_sets_basis_types() {
  auto molecule = constraint_test_molecule();

  auto result = sasmol::apply_constraint_descriptor(
      molecule, sasmol::ConstraintBasis::protein);
  assert(result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"1.00", "1.00", "0.00", "0.00"}));

  result =
      sasmol::apply_constraint_descriptor(molecule, sasmol::ConstraintBasis::nucleic);
  assert(result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"0.00", "0.00", "1.00", "0.00"}));

  result =
      sasmol::apply_constraint_descriptor(molecule, sasmol::ConstraintBasis::solute);
  assert(result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"1.00", "1.00", "1.00", "0.00"}));
}

void test_apply_constraint_descriptor_occupancy_and_reset_false() {
  auto molecule = constraint_test_molecule();
  sasmol::ConstraintPdbOptions options;
  options.field = sasmol::ConstraintField::occupancy;
  options.reset = false;

  const auto result = sasmol::apply_constraint_descriptor(
      molecule, sasmol::ConstraintBasis::nucleic, options);

  assert(result.ok());
  assert((molecule.occupancy() ==
          std::vector<std::string>{"0.25", "0.25", "1.00", "0.25"}));
  assert((molecule.beta() ==
          std::vector<std::string>{"0.50", "0.50", "0.50", "0.50"}));
}

void test_apply_constraint_descriptor_rejects_mismatch_without_mutation() {
  auto molecule = constraint_test_molecule();
  molecule.name().clear();

  const auto result = sasmol::apply_constraint_descriptor(
      molecule, sasmol::ConstraintBasis::backbone);

  assert(!result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"0.50", "0.50", "0.50", "0.50"}));
}

void test_make_constraint_pdb_writes_frame_zero_and_updates_descriptor() {
  auto molecule = constraint_test_molecule();
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_constraint_pdb_test.pdb";

  const auto result =
      sasmol::make_constraint_pdb(molecule, path, sasmol::ConstraintBasis::solute);

  assert(result.ok());
  assert((molecule.beta() ==
          std::vector<std::string>{"1.00", "1.00", "1.00", "0.00"}));

  sasmol::Molecule written;
  const auto status = sasmol::PdbReader{}.read_pdb(path, written);
  std::filesystem::remove(path);

  assert(status.ok());
  assert((written.beta() ==
          std::vector<std::string>{"1.00", "1.00", "1.00", "0.00"}));
}

void test_make_backbone_molecule_from_fasta_protein() {
  const auto result = sasmol::make_backbone_molecule_from_fasta(
      std::vector<std::string>{"T", "C", "P"}, sasmol::BackboneMoltype::protein);

  assert(result.ok());
  assert(result.molecule.natoms() == 3);
  assert((result.molecule.name() ==
          std::vector<std::string>{"CA", "CA", "CA"}));
  assert((result.molecule.resname() ==
          std::vector<std::string>{"THR", "CYS", "PRO"}));
  assert((result.molecule.index() == std::vector<int>{1, 2, 3}));
  assert((result.molecule.resid() == std::vector<int>{1, 2, 3}));
  assert(result.molecule.coordinate(0, 0).x == 0.0F);
}

void test_make_backbone_molecule_from_fasta_nucleic() {
  const auto result = sasmol::make_backbone_molecule_from_fasta(
      std::vector<std::string>{"G", "A", "U"}, sasmol::BackboneMoltype::nucleic);

  assert(result.ok());
  assert(result.molecule.natoms() == 3);
  assert((result.molecule.name() ==
          std::vector<std::string>{"O5'", "O5'", "O5'"}));
  assert((result.molecule.resname() ==
          std::vector<std::string>{"GUA", "ADE", "URA"}));
}

void test_make_backbone_molecule_from_fasta_uses_terminal_patches() {
  auto result = sasmol::make_backbone_molecule_from_fasta(
      std::vector<std::string>{"G", "A"}, sasmol::BackboneMoltype::protein);

  assert(result.ok());
  assert((result.molecule.resname() == std::vector<std::string>{"GLYP", "ALA"}));

  result = sasmol::make_backbone_molecule_from_fasta(
      std::vector<std::string>{"P", "A"}, sasmol::BackboneMoltype::protein);

  assert(result.ok());
  assert((result.molecule.resname() == std::vector<std::string>{"PROP", "ALA"}));
}

void test_make_backbone_molecule_from_formatted_fasta() {
  sasmol::Molecule source;
  source.fasta() = ">demo\nTC\nPA\n";

  const auto result = sasmol::make_backbone_molecule_from_fasta(
      source, sasmol::BackboneMoltype::protein);

  assert(result.ok());
  assert((result.molecule.resname() ==
          std::vector<std::string>{"THR", "CYS", "PRO", "ALA"}));
}

void test_make_backbone_molecule_rejects_unknown_residue() {
  const auto result = sasmol::make_backbone_molecule_from_fasta(
      std::vector<std::string>{"A", "X"}, sasmol::BackboneMoltype::protein);

  assert(!result.ok());
  assert(result.molecule.natoms() == 0);
}

void test_make_backbone_pdb_from_fasta_writes_pdb() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_backbone_fasta_test.pdb";

  const auto result = sasmol::make_backbone_pdb_from_fasta(
      std::vector<std::string>{"G", "A"}, path, sasmol::BackboneMoltype::protein);

  assert(result.ok());

  sasmol::Molecule written;
  const auto status = sasmol::PdbReader{}.read_pdb(path, written);
  std::filesystem::remove(path);

  assert(status.ok());
  assert(written.natoms() == 2);
  assert((written.name() == std::vector<std::string>{"CA", "CA"}));
  assert((written.resname() == std::vector<std::string>{"GLYP", "ALA"}));
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

void test_compare_list_ignore_order_matches_any_order() {
  assert(sasmol::compare_list_ignore_order({"N", "CA", "C"},
                                           {"C", "N", "CA"}));
}

void test_compare_list_ignore_order_rejects_length_mismatch() {
  assert(!sasmol::compare_list_ignore_order({"N", "CA", "C"}, {"N", "CA"}));
}

void test_compare_list_ignore_order_rejects_missing_item() {
  assert(
      !sasmol::compare_list_ignore_order({"N", "CA", "C"}, {"N", "CA", "O"}));
}

void test_compare_list_ignore_order_rejects_duplicate_in_second_list() {
  assert(!sasmol::compare_list_ignore_order({"N", "CA", "C"},
                                            {"N", "CA", "CA"}));
}

void test_compare_list_ignore_order_preserves_python_duplicate_first_list() {
  assert(sasmol::compare_list_ignore_order({"N", "N", "CA"},
                                           {"N", "CA", "C"}));
}

void test_setup_cys_patch_atoms_simple_removes_hg1() {
  const auto atoms = sasmol::setup_cys_patch_atoms_simple(
      {"N", "HN", "CA", "CB", "SG", "HG1", "C", "O"});

  assert((atoms == std::vector<std::string>{"N", "HN", "CA", "CB", "SG", "C",
                                            "O"}));
}

void test_setup_charmm_residue_atoms_builds_atom_name_lists() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "GLY",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"CA", "CT2", "-0.02"},
                 {"C", "C", "0.51"},
                 {"O", "O", "-0.51"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Patch,
       .name = "NTER",
       .total_charge = "1.00",
       .atoms = {{"N", "NH3", "-0.30"}, {"HT1", "HC", "0.33"}}});

  const auto result = sasmol::setup_charmm_residue_atoms(topology);

  assert(result.ok());
  assert(result.residue_atoms.size() == 2);
  assert((result.residue_atoms.at("GLY") ==
          std::vector<std::string>{"N", "CA", "C", "O"}));
  assert((result.residue_atoms.at("NTER") ==
          std::vector<std::string>{"N", "HT1"}));
}

void test_setup_charmm_residue_atoms_adds_disu_from_cys() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "CYS",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"CA", "CT1", "0.07"},
                 {"SG", "S", "-0.23"},
                 {"HG1", "HS", "0.16"},
                 {"C", "C", "0.51"}}});

  const auto result = sasmol::setup_charmm_residue_atoms(topology);

  assert(result.ok());
  assert((result.residue_atoms.at("CYS") ==
          std::vector<std::string>{"N", "CA", "SG", "HG1", "C"}));
  assert((result.residue_atoms.at("DISU") ==
          std::vector<std::string>{"N", "CA", "SG", "C"}));
}

void test_setup_charmm_residue_atoms_omits_disu_without_cys() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "GLY",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"}}});

  const auto result = sasmol::setup_charmm_residue_atoms(topology);

  assert(result.ok());
  assert(result.residue_atoms.find("DISU") == result.residue_atoms.end());
}

void test_patch_charmm_residue_atoms_applies_nter_atom_patch() {
  auto topology = sasmol::parse_charmm_topology(
                      topology_fixture("minimal_pres_atoms_dele.rtf"))
                      .topology;

  const auto result = sasmol::patch_charmm_residue_atoms(topology, "GLY", "NTER");

  assert(result.ok());
  assert(result.patched_entry.name == "GLY_NTER");
  assert(result.patched_entry.total_charge == "0.00");
  assert((result.atom_names ==
          std::vector<std::string>{"N", "HT1", "HT2", "HT3", "CA", "C",
                                   "O"}));
  assert(result.patched_entry.atoms[0].charmm_type == "NH3");
  assert(result.patched_entry.atoms[1].charmm_type == "HC");
}

void test_patch_charmm_residue_atoms_applies_cter_atom_patch() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "GLY",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"CA", "CT2", "-0.02"},
                 {"C", "C", "0.51"},
                 {"O", "O", "-0.51"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Patch,
       .name = "CTER",
       .total_charge = "-1.00",
       .atoms = {{"C", "CC", "0.34"}, {"OT1", "OC", "-0.67"}}});

  const auto result = sasmol::patch_charmm_residue_atoms(topology, "GLY", "CTER");

  assert(result.ok());
  assert(result.patched_entry.name == "GLY_CTER");
  assert((result.atom_names ==
          std::vector<std::string>{"N", "CA", "O", "C", "OT1"}));
  assert(result.patched_entry.atoms[3].charmm_type == "CC");
}

void test_patch_charmm_residue_atoms_reports_missing_residue() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Patch,
       .name = "NTER",
       .total_charge = "1.00",
       .atoms = {{"N", "NH3", "-0.30"}}});

  const auto result = sasmol::patch_charmm_residue_atoms(topology, "GLY", "NTER");

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.patched_entry.atoms.empty());
  assert(result.atom_names.empty());
}

void test_patch_charmm_residue_atoms_reports_missing_patch() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "GLY",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"}}});

  const auto result = sasmol::patch_charmm_residue_atoms(topology, "GLY", "NTER");

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.patched_entry.atoms.empty());
  assert(result.atom_names.empty());
}

sasmol::CharmmTopologyData residue_order_topology() {
  sasmol::CharmmTopologyData topology;
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "ALA",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"HN", "H", "0.31"},
                 {"CA", "CT1", "0.07"},
                 {"C", "C", "0.51"},
                 {"O", "O", "-0.51"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "CYS",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"CA", "CT1", "0.07"},
                 {"SG", "S", "-0.23"},
                 {"HG1", "HS", "0.16"},
                 {"C", "C", "0.51"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "HIS",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"},
                 {"HD1", "H", "0.30"},
                 {"HE2", "H", "0.30"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "HSE",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"}, {"HE2", "H", "0.30"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "HSD",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"}, {"HD1", "H", "0.30"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Residue,
       .name = "HSP",
       .total_charge = "0.00",
       .atoms = {{"N", "NH1", "-0.47"}, {"HD1", "H", "0.30"},
                 {"HE2", "H", "0.30"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Patch,
       .name = "NTER",
       .total_charge = "1.00",
       .atoms = {{"N", "NH3", "-0.30"},
                 {"HT1", "HC", "0.33"},
                 {"HT2", "HC", "0.33"},
                 {"HT3", "HC", "0.33"}},
       .deletes = {.atoms = {"HN"}}});
  topology.entries.push_back(
      {.kind = sasmol::CharmmTopologyEntryKind::Patch,
       .name = "CTER",
       .total_charge = "-1.00",
       .atoms = {{"C", "CC", "0.34"}, {"OT1", "OC", "-0.67"}}});
  return topology;
}

void test_choose_charmm_residue_atom_order_accepts_existing_match() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "ALA", 5, 1, 10,
      {"CA", "O", "N", "HN", "C"});

  assert(result.ok());
  assert(result.topology_residue_name == "ALA");
  assert((result.atom_order ==
          std::vector<std::string>{"N", "HN", "CA", "C", "O"}));
}

void test_choose_charmm_residue_atom_order_applies_nter_patch() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "ALA", 1, 1, 10,
      {"HT2", "CA", "N", "C", "O", "HT1", "HT3"});

  assert(result.ok());
  assert(result.topology_residue_name == "ALA_NTER");
  assert((result.atom_order ==
          std::vector<std::string>{"N", "HT1", "HT2", "HT3", "CA", "C",
                                   "O"}));
}

void test_choose_charmm_residue_atom_order_applies_cter_patch() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "ALA", 10, 1, 10,
      {"N", "HN", "CA", "O", "C", "OT1"});

  assert(result.ok());
  assert(result.topology_residue_name == "ALA_CTER");
  assert((result.atom_order ==
          std::vector<std::string>{"N", "HN", "CA", "O", "C", "OT1"}));
}

void test_choose_charmm_residue_atom_order_uses_disu_for_cys() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "CYS", 4, 1, 10,
      {"N", "CA", "SG", "C"});

  assert(result.ok());
  assert(result.topology_residue_name == "DISU");
  assert((result.atom_order == std::vector<std::string>{"N", "CA", "SG", "C"}));
}

void test_choose_charmm_residue_atom_order_uses_his_variant_fallback() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "HIS", 4, 1, 10, {"HD1", "N"});

  assert(result.ok());
  assert(result.topology_residue_name == "HSD");
  assert((result.atom_order == std::vector<std::string>{"N", "HD1"}));
}

void test_choose_charmm_residue_atom_order_reports_mismatch() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);

  const auto result = sasmol::choose_charmm_residue_atom_order(
      topology, residue_atoms.residue_atoms, "ALA", 5, 1, 10, {"N", "CA"});

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.atom_order.empty());
}

void test_plan_charmm_residue_reorder_indices_maps_topology_order() {
  const auto result = sasmol::plan_charmm_residue_reorder_indices(
      {"CA", "O", "N", "C"}, {"N", "CA", "C", "O"}, "GLY", 7);

  assert(result.ok());
  assert((result.observed_indices == std::vector<std::size_t>{2, 0, 3, 1}));
}

void test_plan_charmm_residue_reorder_indices_uses_first_duplicate_match() {
  const auto result = sasmol::plan_charmm_residue_reorder_indices(
      {"N", "CA", "N"}, {"N", "CA", "N"}, "DUP", 8);

  assert(result.ok());
  assert((result.observed_indices == std::vector<std::size_t>{0, 1, 0}));
}

void test_plan_charmm_residue_reorder_indices_reports_missing_atom() {
  const auto result = sasmol::plan_charmm_residue_reorder_indices(
      {"N", "CA", "C"}, {"N", "CA", "C", "O"}, "GLY", 7);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.observed_indices.empty());
}

void test_plan_charmm_residue_reorder_indices_reports_length_mismatch() {
  const auto result = sasmol::plan_charmm_residue_reorder_indices(
      {"N", "CA", "C", "O", "CB"}, {"N", "CA", "C", "O"}, "GLY", 7);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.observed_indices.empty());
}

void test_plan_charmm_molecule_reorder_plans_terminal_residues() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);
  sasmol::Molecule molecule(13, 1);
  molecule.segname() = {"A", "A", "A", "A", "A", "A", "A",
                        "A", "A", "A", "A", "A", "A"};
  molecule.resid() = {1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  molecule.resname() = {"ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA",
                        "ALA", "ALA", "ALA", "ALA", "ALA", "ALA"};
  molecule.name() = {"HT2", "CA", "N", "C", "O", "HT1", "HT3",
                     "N",   "HN", "CA", "O", "C",   "OT1"};

  const auto result = sasmol::plan_charmm_molecule_reorder(
      molecule, topology, residue_atoms.residue_atoms);

  assert(result.ok());
  assert(result.residues.size() == 2);
  assert(result.residues[0].topology_residue_name == "ALA_NTER");
  assert(result.residues[1].topology_residue_name == "ALA_CTER");
  assert((result.source_atom_indices ==
          std::vector<std::size_t>{2, 5, 0, 6, 1, 3, 4, 7, 8, 9, 10, 11,
                                   12}));
}

void test_plan_charmm_molecule_reorder_reports_mixed_residue_names() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);
  sasmol::Molecule molecule(2, 1);
  molecule.segname() = {"A", "A"};
  molecule.resid() = {1, 1};
  molecule.resname() = {"ALA", "GLY"};
  molecule.name() = {"N", "CA"};

  const auto result = sasmol::plan_charmm_molecule_reorder(
      molecule, topology, residue_atoms.residue_atoms);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.source_atom_indices.empty());
}

void test_plan_charmm_molecule_reorder_reports_descriptor_mismatch() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);
  sasmol::Molecule molecule(2, 1);
  molecule.name() = {"N", "CA"};
  molecule.resname() = {"ALA", "ALA"};
  molecule.resid() = {1, 1};
  molecule.segname().clear();

  const auto result = sasmol::plan_charmm_molecule_reorder(
      molecule, topology, residue_atoms.residue_atoms);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.source_atom_indices.empty());
}

void test_copy_reordered_charmm_molecule_reorders_copy_all_frames() {
  sasmol::Molecule molecule(3, 2);
  molecule.record() = {"ATOM", "HETATM", "ATOM"};
  molecule.index() = {10, 11, 12};
  molecule.original_index() = {100, 101, 102};
  molecule.original_resid() = {1, 1, 1};
  molecule.name() = {"CA", "O", "N"};
  molecule.loc() = {"A", "B", "C"};
  molecule.resname() = {"GLY", "GLY", "GLY"};
  molecule.chain() = {"A", "A", "A"};
  molecule.resid() = {7, 7, 7};
  molecule.rescode() = {"", "", ""};
  molecule.occupancy() = {"1.00", "0.50", "0.25"};
  molecule.beta() = {"10.0", "20.0", "30.0"};
  molecule.segname() = {"SEG", "SEG", "SEG"};
  molecule.element() = {"C", "O", "N"};
  molecule.charge() = {"0", "-1", "1"};
  molecule.atom_charge() = {0.1, -0.2, 0.3};
  molecule.atom_vdw() = {1.1, 1.2, 1.3};
  molecule.residue_flag() = {0, 1, 0};
  molecule.charmm_type() = {"CT2", "O", "NH1"};
  molecule.moltype() = {"protein", "protein", "protein"};
  molecule.mass() = {12.0, 16.0, 14.0};
  molecule.residue_charge() = {0.2, 0.2, 0.2};
  molecule.conect() = {{2}, {1, 3}, {2}};
  molecule.extra_string_descriptors()["tag"] = {"ca", "o", "n"};
  molecule.extra_int_descriptors()["rank"] = {1, 2, 3};
  molecule.extra_calc_descriptors()["score"] = {1.5, 2.5, 3.5};
  molecule.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  molecule.set_coordinate(0, 1, {4.0F, 5.0F, 6.0F});
  molecule.set_coordinate(0, 2, {7.0F, 8.0F, 9.0F});
  molecule.set_coordinate(1, 0, {11.0F, 12.0F, 13.0F});
  molecule.set_coordinate(1, 1, {14.0F, 15.0F, 16.0F});
  molecule.set_coordinate(1, 2, {17.0F, 18.0F, 19.0F});

  sasmol::CharmmMoleculeReorderPlan plan;
  plan.source_atom_indices = {2, 0, 1};

  const auto result = sasmol::copy_reordered_charmm_molecule(molecule, plan);

  assert(result.ok());
  assert((result.molecule.name() == std::vector<std::string>{"N", "CA", "O"}));
  assert((result.molecule.index() == std::vector<int>{12, 10, 11}));
  assert((result.molecule.charmm_type() ==
          std::vector<std::string>{"NH1", "CT2", "O"}));
  assert((result.molecule.extra_string_descriptors().at("tag") ==
          std::vector<std::string>{"n", "ca", "o"}));
  assert((result.molecule.extra_int_descriptors().at("rank") ==
          std::vector<int>{3, 1, 2}));
  assert_close(result.molecule.atom_charge()[0], 0.3);
  assert_close(result.molecule.extra_calc_descriptors().at("score")[1], 1.5);
  assert(result.molecule.coordinate(0, 0).x == 7.0F);
  assert(result.molecule.coordinate(0, 1).x == 1.0F);
  assert(result.molecule.coordinate(1, 0).x == 17.0F);
  assert(result.molecule.coordinate(1, 1).x == 11.0F);
  assert((molecule.name() == std::vector<std::string>{"CA", "O", "N"}));
}

void test_copy_reordered_charmm_molecule_rejects_bad_plan() {
  sasmol::Molecule molecule(2, 1);
  sasmol::CharmmMoleculeReorderPlan plan;
  plan.source_atom_indices = {1};

  const auto result = sasmol::copy_reordered_charmm_molecule(molecule, plan);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.molecule.natoms() == 0);
}

void test_copy_reordered_charmm_molecule_rejects_descriptor_mismatch() {
  sasmol::Molecule molecule(2, 1);
  molecule.extra_string_descriptors()["bad"] = {"only-one"};
  sasmol::CharmmMoleculeReorderPlan plan;
  plan.source_atom_indices = {1, 0};

  const auto result = sasmol::copy_reordered_charmm_molecule(molecule, plan);

  assert(!result.ok());
  assert(!result.errors.empty());
  assert(result.molecule.natoms() == 0);
}

void test_reorder_charmm_molecule_in_place_reorders_after_success() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);
  sasmol::Molecule molecule(13, 2);
  molecule.segname() = {"A", "A", "A", "A", "A", "A", "A",
                        "A", "A", "A", "A", "A", "A"};
  molecule.resid() = {1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  molecule.resname() = {"ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA",
                        "ALA", "ALA", "ALA", "ALA", "ALA", "ALA"};
  molecule.name() = {"HT2", "CA", "N", "C", "O", "HT1", "HT3",
                     "N",   "HN", "CA", "O", "C",   "OT1"};
  molecule.index() = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
  molecule.set_coordinate(0, 0, {1.0F, 0.0F, 0.0F});
  molecule.set_coordinate(0, 2, {3.0F, 0.0F, 0.0F});
  molecule.set_coordinate(0, 7, {8.0F, 0.0F, 0.0F});
  molecule.set_coordinate(1, 0, {11.0F, 0.0F, 0.0F});
  molecule.set_coordinate(1, 2, {13.0F, 0.0F, 0.0F});
  molecule.set_coordinate(1, 7, {18.0F, 0.0F, 0.0F});

  const auto result = sasmol::reorder_charmm_molecule_in_place(
      molecule, topology, residue_atoms.residue_atoms);

  assert(result.ok());
  assert((molecule.name() ==
          std::vector<std::string>{"N", "HT1", "HT2", "HT3", "CA", "C",
                                   "O", "N", "HN", "CA", "O", "C",
                                   "OT1"}));
  assert((molecule.index() ==
          std::vector<int>{3, 6, 1, 7, 2, 4, 5, 8, 9, 10, 11, 12, 13}));
  assert(molecule.coordinate(0, 0).x == 3.0F);
  assert(molecule.coordinate(0, 7).x == 8.0F);
  assert(molecule.coordinate(1, 0).x == 13.0F);
  assert(molecule.coordinate(1, 7).x == 18.0F);
}

void test_reorder_charmm_molecule_in_place_preserves_input_on_failure() {
  const auto topology = residue_order_topology();
  const auto residue_atoms = sasmol::setup_charmm_residue_atoms(topology);
  sasmol::Molecule molecule(2, 1);
  molecule.segname() = {"A", "A"};
  molecule.resid() = {1, 1};
  molecule.resname() = {"ALA", "ALA"};
  molecule.name() = {"N", "CA"};
  molecule.index() = {1, 2};

  const auto result = sasmol::reorder_charmm_molecule_in_place(
      molecule, topology, residue_atoms.residue_atoms);

  assert(!result.ok());
  assert((molecule.name() == std::vector<std::string>{"N", "CA"}));
  assert((molecule.index() == std::vector<int>{1, 2}));
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
  assert((nter.deletes.atoms == std::vector<std::string>{"HN"}));
  assert(nter.deletes.angles.size() == 1);
  assert((nter.deletes.angles[0] ==
          std::vector<std::string>{"HT1", "N", "CA"}));
}

void test_parse_charmm_topology_bond_pairs_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_resi_atoms_bonds.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 1);

  const auto& gly = result.topology.entries[0];
  assert(gly.name == "GLY");
  assert(gly.bonds.size() == 3);
  assert(gly.bonds[0].first == "N");
  assert(gly.bonds[0].second == "CA");
  assert(gly.bonds[1].first == "CA");
  assert(gly.bonds[1].second == "C");
  assert(gly.bonds[2].first == "C");
  assert(gly.bonds[2].second == "O");
  assert(gly.doubles.size() == 1);
  assert(gly.doubles[0].first == "C");
  assert(gly.doubles[0].second == "O");
}

void test_parse_charmm_topology_angle_triples_match_python_oracle_fixture() {
  const auto result =
      sasmol::parse_charmm_topology(topology_fixture("minimal_resi_angles.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 1);

  const auto& gly = result.topology.entries[0];
  assert(gly.name == "GLY");
  assert(gly.angles.size() == 2);
  assert(gly.angles[0].first == "N");
  assert(gly.angles[0].second == "CA");
  assert(gly.angles[0].third == "C");
  assert(gly.angles[1].first == "CA");
  assert(gly.angles[1].second == "C");
  assert(gly.angles[1].third == "O");
  assert(gly.thetas.size() == 1);
  assert(gly.thetas[0].first == "+C");
  assert(gly.thetas[0].second == "N");
  assert(gly.thetas[0].third == "CA");
}

void test_parse_charmm_topology_quad_terms_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_resi_four_terms.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 1);

  const auto& ala = result.topology.entries[0];
  assert(ala.name == "ALA");
  assert(ala.dihedrals.size() == 2);
  assert(ala.dihedrals[0].first == "N");
  assert(ala.dihedrals[0].second == "CA");
  assert(ala.dihedrals[0].third == "C");
  assert(ala.dihedrals[0].fourth == "O");
  assert(ala.dihedrals[1].first == "N");
  assert(ala.dihedrals[1].second == "CA");
  assert(ala.dihedrals[1].third == "CB");
  assert(ala.dihedrals[1].fourth == "C");
  assert(ala.impropers.size() == 1);
  assert(ala.impropers[0].first == "N");
  assert(ala.impropers[0].second == "-C");
  assert(ala.impropers[0].third == "CA");
  assert(ala.impropers[0].fourth == "HN");
  assert(ala.cmaps.size() == 1);
  assert(ala.cmaps[0].first == "-C");
  assert(ala.cmaps[0].second == "N");
  assert(ala.cmaps[0].third == "CA");
  assert(ala.cmaps[0].fourth == "C");
}

void test_parse_charmm_topology_donor_acceptor_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_resi_donor_acceptor.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 1);

  const auto& ser = result.topology.entries[0];
  assert(ser.name == "SER");
  assert((ser.donors == std::vector<std::string>{"HN", "N", "HG1", "OG"}));
  assert((ser.acceptors == std::vector<std::string>{"O", "OG"}));
}

void test_parse_charmm_topology_ic_records_match_python_oracle_fixture() {
  const auto result =
      sasmol::parse_charmm_topology(topology_fixture("minimal_resi_ic.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 1);

  const auto& gly = result.topology.entries[0];
  assert(gly.name == "GLY");
  assert(gly.internal_coordinates.size() == 2);
  assert((gly.internal_coordinates[0].fields ==
          std::vector<std::string>{"-C", "N", "CA", "C", "1.3551",
                                   "126.4900", "180.0000", "114.4400",
                                   "1.5390"}));
  assert((gly.internal_coordinates[1].fields ==
          std::vector<std::string>{"N", "CA", "C", "O", "1.4592",
                                   "114.4400", "180.0000", "120.9900",
                                   "1.2310"}));
}

void test_parse_charmm_topology_ic_preserves_short_python_slice() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_short_ic_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "IC N CA C\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].internal_coordinates.size() == 1);
  assert((result.topology.entries[0].internal_coordinates[0].fields ==
          std::vector<std::string>{"N", "CA", "C"}));
}

void test_parse_charmm_topology_dele_records_match_python_oracle_fixture() {
  const auto result = sasmol::parse_charmm_topology(
      topology_fixture("minimal_pres_atoms_dele.rtf"));

  assert(result.ok());
  assert(result.topology.entries.size() == 2);

  const auto& gly = result.topology.entries[0];
  assert(gly.deletes.atoms.empty());
  assert(gly.deletes.angles.empty());

  const auto& nter = result.topology.entries[1];
  assert((nter.deletes.atoms == std::vector<std::string>{"HN"}));
  assert(nter.deletes.angles.size() == 1);
  assert((nter.deletes.angles[0] ==
          std::vector<std::string>{"HT1", "N", "CA"}));
}

void test_parse_charmm_topology_dele_atom_keeps_python_first_atom_only() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_dele_atom_extra_topology.rtf";
  {
    std::ofstream output(path);
    output << "PRES NTER 1.00\n";
    output << "DELE ATOM HN EXTRA_IGNORED\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert((result.topology.entries[0].deletes.atoms ==
          std::vector<std::string>{"HN"}));
}

void test_parse_charmm_topology_dele_unknown_type_is_ignored() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_dele_unknown_topology.rtf";
  {
    std::ofstream output(path);
    output << "PRES NTER 1.00\n";
    output << "DELE BOND HT1 N\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].deletes.atoms.empty());
  assert(result.topology.entries[0].deletes.angles.empty());
}

void test_parse_charmm_topology_bond_pairs_stop_at_inline_comment() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bond_comment_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "BOND N CA ! CA C ignored\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].bonds.size() == 1);
  assert(result.topology.entries[0].bonds[0].first == "N");
  assert(result.topology.entries[0].bonds[0].second == "CA");
}

void test_parse_charmm_topology_angle_triples_stop_at_inline_comment() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_angle_comment_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "ANGL N CA C ! CA C O ignored\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].angles.size() == 1);
  assert(result.topology.entries[0].angles[0].first == "N");
  assert(result.topology.entries[0].angles[0].second == "CA");
  assert(result.topology.entries[0].angles[0].third == "C");
}

void test_parse_charmm_topology_quad_terms_stop_at_inline_comment() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_quad_comment_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "DIHE N CA C O ! N CA CB C ignored\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].dihedrals.size() == 1);
  assert(result.topology.entries[0].dihedrals[0].first == "N");
  assert(result.topology.entries[0].dihedrals[0].second == "CA");
  assert(result.topology.entries[0].dihedrals[0].third == "C");
  assert(result.topology.entries[0].dihedrals[0].fourth == "O");
}

void test_parse_charmm_topology_donor_acceptor_stop_at_inline_comment() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_donor_comment_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI SER 0.00\n";
    output << "DONO HN N ! HG1 OG ignored\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(result.ok());
  assert(result.topology.entries.size() == 1);
  assert((result.topology.entries[0].donors ==
          std::vector<std::string>{"HN", "N"}));
}

void test_parse_charmm_topology_reports_malformed_bond_pair() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_bond_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "BOND N !\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].bonds.empty());
}

void test_parse_charmm_topology_reports_incomplete_bond_pair() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_incomplete_bond_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "BOND N CA C\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].bonds.size() == 1);
  assert(result.topology.entries[0].bonds[0].first == "N");
  assert(result.topology.entries[0].bonds[0].second == "CA");
}

void test_parse_charmm_topology_reports_malformed_angle_triple() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_angle_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "ANGL N @ C\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].angles.empty());
}

void test_parse_charmm_topology_reports_incomplete_angle_triple() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_incomplete_angle_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "THET N CA C CA\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].thetas.size() == 1);
  assert(result.topology.entries[0].thetas[0].first == "N");
  assert(result.topology.entries[0].thetas[0].second == "CA");
  assert(result.topology.entries[0].thetas[0].third == "C");
}

void test_parse_charmm_topology_reports_malformed_quad_term() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_bad_quad_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "IMPR N CA @ O\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].impropers.empty());
}

void test_parse_charmm_topology_reports_incomplete_quad_term() {
  const auto path =
      std::filesystem::temp_directory_path() / "sasmol_incomplete_quad_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI GLY 0.00\n";
    output << "CMAP N CA C O CA\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert(result.topology.entries[0].cmaps.size() == 1);
  assert(result.topology.entries[0].cmaps[0].first == "N");
  assert(result.topology.entries[0].cmaps[0].second == "CA");
  assert(result.topology.entries[0].cmaps[0].third == "C");
  assert(result.topology.entries[0].cmaps[0].fourth == "O");
}

void test_parse_charmm_topology_reports_malformed_donor_acceptor_token() {
  const auto path = std::filesystem::temp_directory_path() /
                    "sasmol_bad_donor_acceptor_topology.rtf";
  {
    std::ofstream output(path);
    output << "RESI SER 0.00\n";
    output << "ACCE O @\n";
  }

  const auto result = sasmol::parse_charmm_topology(path);
  std::filesystem::remove(path);

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(result.topology.entries.size() == 1);
  assert((result.topology.entries[0].acceptors ==
          std::vector<std::string>{"O"}));
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
  test_create_fasta_default_sequence_and_in_place_string();
  test_create_fasta_formatted_width_and_name();
  test_create_fasta_formatted_by_chain_and_segname();
  test_create_fasta_excludes_hetatm_when_requested();
  test_create_fasta_rejects_unknown_residue_without_mutation();
  test_renumber_default_updates_index_and_resid();
  test_renumber_index_only_preserves_resid();
  test_renumber_resid_only_preserves_index();
  test_renumber_index_and_resid_custom_starts();
  test_renumber_rejects_descriptor_mismatch_without_mutation();
  test_apply_constraint_descriptor_sets_heavy_beta();
  test_apply_constraint_descriptor_sets_basis_types();
  test_apply_constraint_descriptor_occupancy_and_reset_false();
  test_apply_constraint_descriptor_rejects_mismatch_without_mutation();
  test_make_constraint_pdb_writes_frame_zero_and_updates_descriptor();
  test_make_backbone_molecule_from_fasta_protein();
  test_make_backbone_molecule_from_fasta_nucleic();
  test_make_backbone_molecule_from_fasta_uses_terminal_patches();
  test_make_backbone_molecule_from_formatted_fasta();
  test_make_backbone_molecule_rejects_unknown_residue();
  test_make_backbone_pdb_from_fasta_writes_pdb();
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
  test_compare_list_ignore_order_matches_any_order();
  test_compare_list_ignore_order_rejects_length_mismatch();
  test_compare_list_ignore_order_rejects_missing_item();
  test_compare_list_ignore_order_rejects_duplicate_in_second_list();
  test_compare_list_ignore_order_preserves_python_duplicate_first_list();
  test_setup_cys_patch_atoms_simple_removes_hg1();
  test_setup_charmm_residue_atoms_builds_atom_name_lists();
  test_setup_charmm_residue_atoms_adds_disu_from_cys();
  test_setup_charmm_residue_atoms_omits_disu_without_cys();
  test_patch_charmm_residue_atoms_applies_nter_atom_patch();
  test_patch_charmm_residue_atoms_applies_cter_atom_patch();
  test_patch_charmm_residue_atoms_reports_missing_residue();
  test_patch_charmm_residue_atoms_reports_missing_patch();
  test_choose_charmm_residue_atom_order_accepts_existing_match();
  test_choose_charmm_residue_atom_order_applies_nter_patch();
  test_choose_charmm_residue_atom_order_applies_cter_patch();
  test_choose_charmm_residue_atom_order_uses_disu_for_cys();
  test_choose_charmm_residue_atom_order_uses_his_variant_fallback();
  test_choose_charmm_residue_atom_order_reports_mismatch();
  test_plan_charmm_residue_reorder_indices_maps_topology_order();
  test_plan_charmm_residue_reorder_indices_uses_first_duplicate_match();
  test_plan_charmm_residue_reorder_indices_reports_missing_atom();
  test_plan_charmm_residue_reorder_indices_reports_length_mismatch();
  test_plan_charmm_molecule_reorder_plans_terminal_residues();
  test_plan_charmm_molecule_reorder_reports_mixed_residue_names();
  test_plan_charmm_molecule_reorder_reports_descriptor_mismatch();
  test_copy_reordered_charmm_molecule_reorders_copy_all_frames();
  test_copy_reordered_charmm_molecule_rejects_bad_plan();
  test_copy_reordered_charmm_molecule_rejects_descriptor_mismatch();
  test_reorder_charmm_molecule_in_place_reorders_after_success();
  test_reorder_charmm_molecule_in_place_preserves_input_on_failure();
  test_parse_charmm_topology_globals_matches_python_oracle_fixture();
  test_parse_charmm_topology_globals_reports_malformed_records();
  test_parse_charmm_topology_residue_atoms_match_python_oracle_fixture();
  test_parse_charmm_topology_patch_atoms_match_python_oracle_fixture();
  test_parse_charmm_topology_bond_pairs_match_python_oracle_fixture();
  test_parse_charmm_topology_angle_triples_match_python_oracle_fixture();
  test_parse_charmm_topology_quad_terms_match_python_oracle_fixture();
  test_parse_charmm_topology_donor_acceptor_match_python_oracle_fixture();
  test_parse_charmm_topology_ic_records_match_python_oracle_fixture();
  test_parse_charmm_topology_ic_preserves_short_python_slice();
  test_parse_charmm_topology_dele_records_match_python_oracle_fixture();
  test_parse_charmm_topology_dele_atom_keeps_python_first_atom_only();
  test_parse_charmm_topology_dele_unknown_type_is_ignored();
  test_parse_charmm_topology_bond_pairs_stop_at_inline_comment();
  test_parse_charmm_topology_angle_triples_stop_at_inline_comment();
  test_parse_charmm_topology_quad_terms_stop_at_inline_comment();
  test_parse_charmm_topology_donor_acceptor_stop_at_inline_comment();
  test_parse_charmm_topology_reports_malformed_bond_pair();
  test_parse_charmm_topology_reports_incomplete_bond_pair();
  test_parse_charmm_topology_reports_malformed_angle_triple();
  test_parse_charmm_topology_reports_incomplete_angle_triple();
  test_parse_charmm_topology_reports_malformed_quad_term();
  test_parse_charmm_topology_reports_incomplete_quad_term();
  test_parse_charmm_topology_reports_malformed_donor_acceptor_token();
  test_parse_charmm_topology_reports_atom_before_residue_or_patch();
  return 0;
}
