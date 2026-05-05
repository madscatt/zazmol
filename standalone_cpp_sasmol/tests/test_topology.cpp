#include "sasmol/topology.hpp"

#include <cassert>
#include <cmath>

namespace {

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
  return 0;
}
