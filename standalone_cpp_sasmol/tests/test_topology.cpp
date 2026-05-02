#include "sasmol/topology.hpp"

#include <cassert>

namespace {

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

}  // namespace

int main() {
  test_assign_charmm_types_sets_atom_aligned_values();
  test_assign_charmm_types_rejects_length_mismatch_without_mutation();
  test_assign_charmm_types_allows_empty_molecule_empty_types();
  test_assign_charmm_types_from_atom_table_validates_names();
  test_assign_charmm_types_from_atom_table_rejects_name_mismatch();
  test_assign_charmm_types_from_atom_table_rejects_name_length_mismatch();
  return 0;
}
