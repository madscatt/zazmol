#include "sasmol/file_io.hpp"
#include "sasmol/selection.hpp"
#include "sasmol/subset.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-6) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

void assert_vec_close(sasmol::Vec3 actual, sasmol::Vec3 expected,
                      double tolerance = 1.0e-6) {
  assert_close(actual.x, expected.x, tolerance);
  assert_close(actual.y, expected.y, tolerance);
  assert_close(actual.z, expected.z, tolerance);
}

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

sasmol::Molecule read_fixture(const char* name) {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status = reader.read_pdb(fixture_path("pdb_common", name), mol);
  assert(status.ok());
  return mol;
}

void test_get_coordinates_using_indices() {
  const auto mol = read_fixture("2AAD.pdb");
  const std::vector<std::size_t> indices{0, 4, 14};

  const auto result = sasmol::get_coordinates_using_indices(mol, 0, indices);

  assert(result.ok());
  assert(result.coordinates.size() == 3);
  assert_vec_close(result.coordinates[0], mol.coordinate(0, 0));
  assert_vec_close(result.coordinates[1], mol.coordinate(0, 4));
  assert_vec_close(result.coordinates[2], mol.coordinate(0, 14));
}

void test_get_coordinates_rejects_bad_index() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result = sasmol::get_coordinates_using_indices(mol, 0, {0, 99});

  assert(!result.ok());
  assert(result.coordinates.empty());
}

void test_get_indices_from_mask() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result =
      sasmol::get_indices_from_mask(mol, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1});

  assert(result.ok());
  assert((result.indices == std::vector<std::size_t>{0, 4, 14}));
}

void test_get_coordinates_using_mask_delegates_to_indices() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result =
      sasmol::get_coordinates_using_mask(mol, 0,
                                         {1, 0, 0, 0, 1, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 1});

  assert(result.ok());
  assert(result.coordinates.size() == 3);
  assert_vec_close(result.coordinates[1], mol.coordinate(0, 4));
}

void test_select_mask_flows_into_coordinate_subset() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto selected = sasmol::select_mask(mol, "name[i] == \"CA\"");
  assert(selected.ok());
  const auto result = sasmol::get_coordinates_using_mask(mol, 0, selected.mask);

  assert(result.ok());
  assert(result.coordinates.size() == 2);
  assert_vec_close(result.coordinates[0], mol.coordinate(0, 1));
  assert_vec_close(result.coordinates[1], mol.coordinate(0, 9));
}

void test_named_basis_mask_flows_into_copy_subset() {
  auto source = read_fixture("2AAD.pdb");
  source.name()[0] = "H1";
  source.name()[4] = "HA";
  sasmol::Molecule subset;

  const auto selected = sasmol::select_named_basis_mask(source, "heavy");
  assert(selected.ok());
  const auto result =
      sasmol::copy_molecule_using_mask(source, subset, selected.mask, 0);

  assert(result.ok());
  assert(subset.natoms() == source.natoms() - 2);
  for (const auto& name : subset.name()) {
    assert(name.empty() || name.front() != 'H');
  }
}

void test_mask_rejects_bad_shape_and_values() {
  const auto mol = read_fixture("2AAD.pdb");

  auto result = sasmol::get_indices_from_mask(mol, {1, 0});
  assert(!result.ok());
  assert(result.indices.empty());

  result = sasmol::get_indices_from_mask(
      mol, {1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1});
  assert(!result.ok());
  assert(result.indices.empty());
}

void test_get_dihedral_subset_mask_protein() {
  sasmol::Molecule mol(18, 1);
  mol.name() = {"N",  "CA", "C",  "O",  "CB", "CG",
                "N",  "CA", "C",  "O",  "CB", "CG",
                "N",  "CA", "C",  "O",  "CB", "CG"};
  mol.resid() = {1, 1, 1, 1, 1, 1, 2, 2, 2,
                 2, 2, 2, 3, 3, 3, 3, 3, 3};

  const auto result = sasmol::get_dihedral_subset_mask(mol, {2}, 0);

  assert(result.ok());
  assert(result.masks.size() == 1);
  assert((result.masks[0] ==
          std::vector<int>{0, 0, 1, 0, 0, 0, 1, 1, 1,
                           0, 0, 0, 1, 0, 0, 0, 0, 0}));
}

void test_get_dihedral_subset_mask_rna() {
  sasmol::Molecule mol(18, 1);
  mol.name() = {"P", "O5'", "C5'", "C4'", "C3'", "O3'",
                "P", "O5'", "C5'", "C4'", "C3'", "O3'",
                "P", "O5'", "C5'", "C4'", "C3'", "O3'"};
  mol.resid() = {1, 1, 1, 1, 1, 1, 2, 2, 2,
                 2, 2, 2, 3, 3, 3, 3, 3, 3};

  const auto result = sasmol::get_dihedral_subset_mask(mol, {2}, 1);

  assert(result.ok());
  assert(result.masks.size() == 1);
  assert((result.masks[0] ==
          std::vector<int>{0, 0, 0, 0, 0, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 0, 0, 0, 0}));
}

void test_get_dihedral_subset_mask_rejects_bad_inputs() {
  sasmol::Molecule mol(1, 1);
  mol.name().clear();

  auto result = sasmol::get_dihedral_subset_mask(mol, {1}, 0);
  assert(!result.ok());
  assert(result.masks.empty());

  mol.name() = {"N"};
  result = sasmol::get_dihedral_subset_mask(mol, {1}, 99);
  assert(!result.ok());
  assert(result.masks.empty());
}

void test_copy_molecule_using_indices_preserves_descriptors() {
  const auto source = read_fixture("1CRN.pdb");
  const auto selected = sasmol::select_indices(
      source, "name[i] == \"CA\" and (resid[i] >= 20 and resid[i] <= 31)");
  assert(selected.ok());
  sasmol::Molecule subset;

  const auto result =
      sasmol::copy_molecule_using_indices(source, subset, selected.indices, 0);

  assert(result.ok());
  assert(subset.natoms() == selected.indices.size());
  assert(subset.number_of_frames() == 1);
  for (std::size_t atom = 0; atom < selected.indices.size(); ++atom) {
    const auto source_atom = selected.indices[atom];
    assert(subset.name()[atom] == source.name()[source_atom]);
    assert(subset.resid()[atom] == source.resid()[source_atom]);
    assert(subset.resname()[atom] == source.resname()[source_atom]);
    assert_vec_close(subset.coordinate(0, atom), source.coordinate(0, source_atom));
  }
}

void test_copy_molecule_using_mask_preserves_descriptors() {
  const auto source = read_fixture("2AAD.pdb");
  sasmol::Molecule subset;

  const auto result = sasmol::copy_molecule_using_mask(
      source, subset, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 0);

  assert(result.ok());
  assert(subset.natoms() == 3);
  assert(subset.name()[0] == source.name()[0]);
  assert(subset.name()[1] == source.name()[4]);
  assert(subset.name()[2] == source.name()[14]);
}

void test_copy_molecule_using_indices_preserves_extended_descriptors() {
  auto source = read_fixture("2AAD.pdb");
  source.residue_flag()[0] = 1;
  source.residue_flag()[4] = 1;
  source.charmm_type().assign(source.natoms(), "");
  source.charmm_type()[0] = "CT1";
  source.charmm_type()[4] = "NH1";
  source.extra_string_descriptors()["trial_name"] =
      std::vector<std::string>(source.natoms(), "");
  source.extra_string_descriptors()["trial_name"][4] = "keep";
  source.extra_int_descriptors()["trial_int"] =
      std::vector<int>(source.natoms(), 0);
  source.extra_int_descriptors()["trial_int"][4] = 7;
  source.extra_calc_descriptors()["trial_calc"] =
      std::vector<sasmol::calc_type>(source.natoms(), 0.0);
  source.extra_calc_descriptors()["trial_calc"][4] = 3.5;
  sasmol::Molecule subset;

  const auto result =
      sasmol::copy_molecule_using_indices(source, subset, {0, 4}, 0);

  assert(result.ok());
  assert((subset.residue_flag() == std::vector<int>{1, 1}));
  assert((subset.charmm_type() == std::vector<std::string>{"CT1", "NH1"}));
  assert((subset.extra_string_descriptors().at("trial_name") ==
          std::vector<std::string>{"", "keep"}));
  assert((subset.extra_int_descriptors().at("trial_int") ==
          std::vector<int>{0, 7}));
  assert(std::fabs(subset.extra_calc_descriptors().at("trial_calc")[1] - 3.5) <
         1.0e-12);
}

void test_copied_molecule_using_indices_returns_value() {
  const auto source = read_fixture("2AAD.pdb");

  const auto subset = sasmol::copied_molecule_using_indices(source, {0, 1}, 0);

  assert(subset.natoms() == 2);
  assert(subset.index()[0] == source.index()[0]);
  assert(subset.index()[1] == source.index()[1]);
}

void test_copied_molecule_using_mask_returns_value() {
  const auto source = read_fixture("2AAD.pdb");

  const auto subset = sasmol::copied_molecule_using_mask(
      source, {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 0);

  assert(subset.natoms() == 3);
  assert(subset.name()[0] == source.name()[0]);
  assert(subset.name()[1] == source.name()[4]);
  assert(subset.name()[2] == source.name()[14]);
  assert_vec_close(subset.coordinate(0, 1), source.coordinate(0, 4));
}

void test_duplicate_molecule_returns_deep_value_copies() {
  auto source = read_fixture("2AAD.pdb");
  source.formula() = {{"C", 2}, {"N", 1}};
  source.fasta() = "AG";
  source.unitcell() = {1.0, 2.0, 3.0, 90.0, 90.0, 120.0};

  auto duplicates = sasmol::duplicate_molecule(source, 2);

  assert(duplicates.size() == 2);
  assert(duplicates[0].natoms() == source.natoms());
  assert(duplicates[0].number_of_frames() == source.number_of_frames());
  assert(duplicates[0].name() == source.name());
  assert(duplicates[0].resid() == source.resid());
  assert(duplicates[0].coor() == source.coor());
  assert(duplicates[0].conect() == source.conect());
  assert(duplicates[0].formula() == source.formula());
  assert(duplicates[0].fasta() == source.fasta());
  assert(duplicates[0].unitcell() == source.unitcell());

  duplicates[0].name()[0] = "XX";
  duplicates[0].set_coordinate(0, 0, {10.0F, 20.0F, 30.0F});
  duplicates[0].formula()["C"] = 99;

  assert(source.name()[0] != "XX");
  assert(duplicates[1].name()[0] == source.name()[0]);
  assert_vec_close(source.coordinate(0, 0), duplicates[1].coordinate(0, 0));
  assert(source.formula().at("C") == 2);
  assert(duplicates[1].formula().at("C") == 2);
}

void test_duplicate_molecule_allows_zero_duplicates() {
  const auto source = read_fixture("2AAD.pdb");

  const auto duplicates = sasmol::duplicate_molecule(source, 0);

  assert(duplicates.empty());
}

void test_merge_two_molecules_combines_core_descriptors_and_coordinates() {
  const auto mol1 = read_fixture("1ATM.pdb");
  const auto mol2 = read_fixture("2AAD.pdb");
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert(merged.natoms() == mol1.natoms() + mol2.natoms());
  assert(merged.number_of_frames() == 1);
  assert(merged.name()[0] == mol1.name()[0]);
  assert(merged.name()[mol1.natoms()] == mol2.name()[0]);
  assert(merged.resid()[0] == mol1.resid()[0]);
  assert(merged.resid()[mol1.natoms()] == mol2.resid()[0]);
  assert_vec_close(merged.coordinate(0, 0), mol1.coordinate(0, 0));
  assert_vec_close(merged.coordinate(0, mol1.natoms()), mol2.coordinate(0, 0));
}

void test_merge_two_molecules_combines_extended_descriptors() {
  sasmol::Molecule mol1(1, 1);
  sasmol::Molecule mol2(1, 1);
  mol1.residue_flag()[0] = 1;
  mol2.residue_flag()[0] = 0;
  mol1.charmm_type() = {"CT1"};
  mol2.charmm_type() = {"NH1"};
  mol1.extra_string_descriptors()["trial"] = {"left"};
  mol2.extra_string_descriptors()["trial"] = {"right"};
  mol1.extra_int_descriptors()["trial"] = {1};
  mol2.extra_int_descriptors()["trial"] = {2};
  mol1.extra_calc_descriptors()["trial"] = {1.5};
  mol2.extra_calc_descriptors()["trial"] = {2.5};
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert((merged.residue_flag() == std::vector<int>{1, 0}));
  assert((merged.charmm_type() == std::vector<std::string>{"CT1", "NH1"}));
  assert((merged.extra_string_descriptors().at("trial") ==
          std::vector<std::string>{"left", "right"}));
  assert((merged.extra_int_descriptors().at("trial") ==
          std::vector<int>{1, 2}));
  assert(std::fabs(merged.extra_calc_descriptors().at("trial")[0] - 1.5) <
         1.0e-12);
  assert(std::fabs(merged.extra_calc_descriptors().at("trial")[1] - 2.5) <
         1.0e-12);
}

void test_merge_two_molecules_reports_skipped_extra_descriptor() {
  sasmol::Molecule mol1(1, 1);
  sasmol::Molecule mol2(1, 1);
  mol1.extra_string_descriptors()["trial"] = {"left"};
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(
      mol1, mol2, merged, {.report_skipped_descriptors = true});

  assert(!result.ok());
  assert(merged.natoms() == 2);
  assert(merged.extra_string_descriptors().empty());
}

void test_merge_two_molecules_regenerates_second_indices() {
  auto mol1 = read_fixture("1ATM.pdb");
  const auto mol2 = read_fixture("2AAD.pdb");
  mol1.index()[0] = 20;
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert(merged.index()[0] == 20);
  for (std::size_t atom = 0; atom < mol2.natoms(); ++atom) {
    assert(merged.index()[mol1.natoms() + atom] == 21 + static_cast<int>(atom));
  }
}

void test_merge_two_molecules_uses_frame_zero_only() {
  const auto mol1 = read_fixture("1ATM.pdb");
  const auto mol2 = read_fixture("2AAD-1to3.pdb");
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert(merged.number_of_frames() == 1);
  assert_vec_close(merged.coordinate(0, mol1.natoms()), mol2.coordinate(0, 0));
}

void test_merge_two_molecules_rejects_empty_first_molecule() {
  const sasmol::Molecule mol1;
  const auto mol2 = read_fixture("1ATM.pdb");
  sasmol::Molecule merged(1, 1);
  merged.name()[0] = "KEEP";

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(!result.ok());
  assert(merged.natoms() == 1);
  assert(merged.name()[0] == "KEEP");
}

void test_merge_two_molecules_allows_empty_second_molecule() {
  auto mol1 = read_fixture("1ATM.pdb");
  mol1.atom_charge()[0] = 0.25;
  const sasmol::Molecule mol2;
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert(merged.natoms() == mol1.natoms());
  assert(merged.name() == mol1.name());
  assert(merged.atom_charge() == mol1.atom_charge());
}

void test_merge_two_molecules_rejects_bad_coordinates_before_mutation() {
  auto mol1 = read_fixture("1ATM.pdb");
  const auto mol2 = read_fixture("2AAD.pdb");
  mol1.coor().pop_back();
  sasmol::Molecule merged(1, 1);
  merged.name()[0] = "KEEP";

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(!result.ok());
  assert(merged.natoms() == 1);
  assert(merged.name()[0] == "KEEP");
}

void test_merge_two_molecules_reports_skipped_optional_numeric_descriptor() {
  auto mol1 = read_fixture("1ATM.pdb");
  auto mol2 = read_fixture("1ATM.pdb");
  mol2.atom_charge().clear();
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(
      mol1, mol2, merged, {.report_skipped_descriptors = true});

  assert(!result.ok());
  assert(result.errors.size() == 1);
  assert(merged.natoms() == 2);
  assert(merged.atom_charge().empty());
}

void test_merge_two_molecules_copies_conect_without_aliasing() {
  sasmol::Molecule mol1(1, 1);
  sasmol::Molecule mol2(1, 1);
  mol1.name()[0] = "A";
  mol2.name()[0] = "B";
  mol1.original_index()[0] = 10;
  mol2.original_index()[0] = 20;
  mol1.conect()[0] = {20};
  mol2.conect()[0] = {10};
  sasmol::Molecule merged;

  const auto result = sasmol::merge_two_molecules(mol1, mol2, merged);

  assert(result.ok());
  assert((merged.conect()[0] == std::vector<int>{20}));
  assert((merged.conect()[1] == std::vector<int>{10}));
  mol1.conect()[0].clear();
  assert((merged.conect()[0] == std::vector<int>{20}));
}

void test_set_coordinates_using_mask_replaces_selected_atoms_only() {
  sasmol::Molecule target(3, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  target.set_coordinate(0, 2, {3.0F, 3.0F, 3.0F});
  sasmol::Molecule source(2, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});
  source.set_coordinate(0, 1, {20.0F, 21.0F, 22.0F});

  const auto result =
      sasmol::set_coordinates_using_mask(target, source, 0, {1, 0, 1});

  assert(result.ok());
  assert_vec_close(target.coordinate(0, 0), {10.0F, 11.0F, 12.0F});
  assert_vec_close(target.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
  assert_vec_close(target.coordinate(0, 2), {20.0F, 21.0F, 22.0F});
}

void test_set_coordinates_using_mask_rejects_before_mutation() {
  sasmol::Molecule target(2, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  sasmol::Molecule source(1, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});

  const auto result =
      sasmol::set_coordinates_using_mask(target, source, 0, {1, 2});

  assert(!result.ok());
  assert_vec_close(target.coordinate(0, 0), {1.0F, 1.0F, 1.0F});
  assert_vec_close(target.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
}

void test_set_coordinates_using_indices_replaces_selected_atoms_only() {
  sasmol::Molecule target(3, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  target.set_coordinate(0, 2, {3.0F, 3.0F, 3.0F});
  sasmol::Molecule source(2, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});
  source.set_coordinate(0, 1, {20.0F, 21.0F, 22.0F});

  const auto result =
      sasmol::set_coordinates_using_indices(target, source, 0, {0, 2});

  assert(result.ok());
  assert_vec_close(target.coordinate(0, 0), {10.0F, 11.0F, 12.0F});
  assert_vec_close(target.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
  assert_vec_close(target.coordinate(0, 2), {20.0F, 21.0F, 22.0F});
}

void test_set_coordinates_rejects_shape_mismatch() {
  sasmol::Molecule target(3, 1);
  sasmol::Molecule source(1, 1);

  const auto result =
      sasmol::set_coordinates_using_indices(target, source, 0, {0, 2});

  assert(!result.ok());
}

void test_with_coordinates_using_indices_returns_modified_copy() {
  sasmol::Molecule target(3, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  target.set_coordinate(0, 2, {3.0F, 3.0F, 3.0F});
  sasmol::Molecule source(2, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});
  source.set_coordinate(0, 1, {20.0F, 21.0F, 22.0F});

  const auto updated =
      sasmol::with_coordinates_using_indices(target, source, 0, {0, 2});

  assert_vec_close(updated.coordinate(0, 0), {10.0F, 11.0F, 12.0F});
  assert_vec_close(updated.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
  assert_vec_close(updated.coordinate(0, 2), {20.0F, 21.0F, 22.0F});
  assert_vec_close(target.coordinate(0, 0), {1.0F, 1.0F, 1.0F});
  assert_vec_close(target.coordinate(0, 2), {3.0F, 3.0F, 3.0F});
}

void test_with_coordinates_using_mask_returns_modified_copy() {
  sasmol::Molecule target(3, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  target.set_coordinate(0, 2, {3.0F, 3.0F, 3.0F});
  sasmol::Molecule source(2, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});
  source.set_coordinate(0, 1, {20.0F, 21.0F, 22.0F});

  const auto updated =
      sasmol::with_coordinates_using_mask(target, source, 0, {1, 0, 1});

  assert_vec_close(updated.coordinate(0, 0), {10.0F, 11.0F, 12.0F});
  assert_vec_close(updated.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
  assert_vec_close(updated.coordinate(0, 2), {20.0F, 21.0F, 22.0F});
  assert_vec_close(target.coordinate(0, 0), {1.0F, 1.0F, 1.0F});
  assert_vec_close(target.coordinate(0, 2), {3.0F, 3.0F, 3.0F});
}

void test_with_coordinates_rejects_bad_mask_without_mutating_target() {
  sasmol::Molecule target(2, 1);
  target.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});
  target.set_coordinate(0, 1, {2.0F, 2.0F, 2.0F});
  sasmol::Molecule source(1, 1);
  source.set_coordinate(0, 0, {10.0F, 11.0F, 12.0F});

  bool threw = false;
  try {
    (void)sasmol::with_coordinates_using_mask(target, source, 0, {1, 2});
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
  assert_vec_close(target.coordinate(0, 0), {1.0F, 1.0F, 1.0F});
  assert_vec_close(target.coordinate(0, 1), {2.0F, 2.0F, 2.0F});
}

void test_copy_conect_filters_to_selected_atoms() {
  sasmol::Molecule source(3, 1);
  source.index() = {1, 2, 3};
  source.original_index() = {10, 20, 30};
  source.conect()[0] = {20, 30};
  source.conect()[1] = {10};
  source.conect()[2] = {10};
  sasmol::Molecule subset;

  const auto result =
      sasmol::copy_molecule_using_indices(source, subset, {0, 1}, 0);

  assert(result.ok());
  assert((subset.conect()[0] == std::vector<int>{20}));
  assert((subset.conect()[1] == std::vector<int>{10}));
}

void test_get_and_set_string_descriptor_using_indices() {
  auto mol = read_fixture("2AAD.pdb");

  auto values = sasmol::get_string_descriptor_using_indices(
      mol, sasmol::StringDescriptor::beta, {0, 4, 14});
  assert(values.ok());
  assert((values.values ==
          std::vector<std::string>{mol.beta()[0], mol.beta()[4],
                                   mol.beta()[14]}));

  const auto result = sasmol::set_string_descriptor_using_indices(
      mol, sasmol::StringDescriptor::beta, {0, 14}, "7.50");

  assert(result.ok());
  assert(mol.beta()[0] == "7.50");
  assert(mol.beta()[4] == values.values[1]);
  assert(mol.beta()[14] == "7.50");
}

void test_get_and_set_string_descriptor_using_mask() {
  auto mol = read_fixture("2AAD.pdb");

  const auto result = sasmol::set_string_descriptor_using_mask(
      mol, sasmol::StringDescriptor::segname,
      {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, "TEST");

  assert(result.ok());
  auto values = sasmol::get_string_descriptor_using_mask(
      mol, sasmol::StringDescriptor::segname,
      {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1});
  assert(values.ok());
  assert((values.values == std::vector<std::string>{"TEST", "TEST", "TEST"}));
}

void test_get_and_set_int_descriptor_using_indices() {
  auto mol = read_fixture("2AAD.pdb");
  const auto unchanged_resid = mol.resid()[4];

  const auto result = sasmol::set_int_descriptor_using_indices(
      mol, sasmol::IntDescriptor::resid, {0, 14}, 99);

  assert(result.ok());
  auto values = sasmol::get_int_descriptor_using_indices(
      mol, sasmol::IntDescriptor::resid, {0, 4, 14});
  assert(values.ok());
  assert((values.values == std::vector<int>{99, unchanged_resid, 99}));
}

void test_get_and_set_calc_descriptor_using_mask() {
  sasmol::Molecule mol(3, 1);
  mol.atom_charge() = {0.1, 0.2, 0.3};

  const auto result = sasmol::set_calc_descriptor_using_mask(
      mol, sasmol::CalcDescriptor::atom_charge, {1, 0, 1}, -1.5);

  assert(result.ok());
  auto values = sasmol::get_calc_descriptor_using_indices(
      mol, sasmol::CalcDescriptor::atom_charge, {0, 1, 2});
  assert(values.ok());
  assert(std::fabs(values.values[0] + 1.5) < 1.0e-12);
  assert(std::fabs(values.values[1] - 0.2) < 1.0e-12);
  assert(std::fabs(values.values[2] + 1.5) < 1.0e-12);
}

void test_descriptor_set_rejects_bad_mask_before_mutation() {
  auto mol = read_fixture("2AAD.pdb");
  const auto original = mol.beta()[0];

  const auto result = sasmol::set_string_descriptor_using_mask(
      mol, sasmol::StringDescriptor::beta,
      {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "9.99");

  assert(!result.ok());
  assert(mol.beta()[0] == original);
}

void test_descriptor_get_rejects_descriptor_length_mismatch() {
  auto mol = read_fixture("2AAD.pdb");
  mol.beta().clear();

  const auto result = sasmol::get_string_descriptor_using_indices(
      mol, sasmol::StringDescriptor::beta, {0});

  assert(!result.ok());
  assert(result.values.empty());
}

void test_apply_biomt_transforms_assembles_in_transform_order() {
  sasmol::Molecule source(1, 1);
  source.name()[0] = "N";
  source.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});

  const std::vector<sasmol::BiomtTransform> transforms{
      {{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
       {0.0, 0.0, 0.0}},
      {{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
       {10.0, -1.0, 0.5}},
  };

  sasmol::Molecule transformed;
  const auto result =
      sasmol::apply_biomt_transforms(source, 0, transforms, transformed);

  assert(result.ok());
  assert(transformed.number_of_frames() == 1);
  assert(transformed.natoms() == 2);
  assert_vec_close(transformed.coordinate(0, 0), {1.0F, 2.0F, 3.0F});
  assert_vec_close(transformed.coordinate(0, 1), {11.0F, 1.0F, 3.5F});
}

void test_biomt_transformed_returns_value_without_mutating_source() {
  sasmol::Molecule source(2, 1);
  source.name() = {"A", "B"};
  source.set_coordinate(0, 0, {0.0F, 0.0F, 0.0F});
  source.set_coordinate(0, 1, {2.0F, 3.0F, 4.0F});
  const auto original = source.coor();

  const std::vector<sasmol::BiomtTransform> transforms{
      {{{{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}},
       {5.0, 0.0, 0.0}},
  };

  const auto transformed = sasmol::biomt_transformed(source, 0, transforms);

  assert(transformed.natoms() == 2);
  assert_vec_close(transformed.coordinate(0, 0), {5.0F, 0.0F, 0.0F});
  assert_vec_close(transformed.coordinate(0, 1), {2.0F, 2.0F, 4.0F});
  assert(source.coor() == original);
}

void test_apply_biomt_transforms_from_metadata_uses_recorded_transforms() {
  sasmol::Molecule source(1, 1);
  source.name()[0] = "N";
  source.set_coordinate(0, 0, {2.0F, 0.0F, 0.0F});

  sasmol::BiomtRecord record;
  record.rot.push_back({{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}});
  record.trans.push_back({0.0, 0.0, 0.0});
  record.rot.push_back({{{-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}});
  record.trans.push_back({1.0, 0.0, 0.0});
  source.biomt()[1] = record;

  sasmol::Molecule transformed;
  const auto result =
      sasmol::apply_biomt_transforms_from_metadata(source, 0, 1, transformed);

  assert(result.ok());
  assert(transformed.natoms() == 2);
  assert_vec_close(transformed.coordinate(0, 0), {2.0F, 0.0F, 0.0F});
  assert_vec_close(transformed.coordinate(0, 1), {-1.0F, 0.0F, 0.0F});
}

void test_apply_biomt_transforms_rejects_invalid_inputs_without_mutation() {
  sasmol::Molecule source(1, 1);
  source.name()[0] = "N";
  source.set_coordinate(0, 0, {1.0F, 1.0F, 1.0F});

  sasmol::Molecule transformed(1, 1);
  transformed.name()[0] = "KEEP";
  transformed.set_coordinate(0, 0, {9.0F, 9.0F, 9.0F});

  auto result = sasmol::apply_biomt_transforms(source, 0, {}, transformed);
  assert(!result.ok());
  assert(transformed.natoms() == 1);
  assert(transformed.name()[0] == "KEEP");
  assert_vec_close(transformed.coordinate(0, 0), {9.0F, 9.0F, 9.0F});

  const std::vector<sasmol::BiomtTransform> non_finite{
      {{{{1.0, 0.0, 0.0}, {0.0, NAN, 0.0}, {0.0, 0.0, 1.0}}},
       {0.0, 0.0, 0.0}},
  };
  result = sasmol::apply_biomt_transforms(source, 0, non_finite, transformed);
  assert(!result.ok());
  assert(transformed.name()[0] == "KEEP");

  result = sasmol::apply_biomt_transforms_from_metadata(source, 0, 9, transformed);
  assert(!result.ok());
  assert(transformed.name()[0] == "KEEP");
}

}  // namespace

int main() {
  test_get_coordinates_using_indices();
  test_get_coordinates_rejects_bad_index();
  test_get_indices_from_mask();
  test_get_coordinates_using_mask_delegates_to_indices();
  test_select_mask_flows_into_coordinate_subset();
  test_named_basis_mask_flows_into_copy_subset();
  test_mask_rejects_bad_shape_and_values();
  test_get_dihedral_subset_mask_protein();
  test_get_dihedral_subset_mask_rna();
  test_get_dihedral_subset_mask_rejects_bad_inputs();
  test_copy_molecule_using_indices_preserves_descriptors();
  test_copy_molecule_using_mask_preserves_descriptors();
  test_copy_molecule_using_indices_preserves_extended_descriptors();
  test_copied_molecule_using_indices_returns_value();
  test_copied_molecule_using_mask_returns_value();
  test_duplicate_molecule_returns_deep_value_copies();
  test_duplicate_molecule_allows_zero_duplicates();
  test_merge_two_molecules_combines_core_descriptors_and_coordinates();
  test_merge_two_molecules_combines_extended_descriptors();
  test_merge_two_molecules_reports_skipped_extra_descriptor();
  test_merge_two_molecules_regenerates_second_indices();
  test_merge_two_molecules_uses_frame_zero_only();
  test_merge_two_molecules_rejects_empty_first_molecule();
  test_merge_two_molecules_allows_empty_second_molecule();
  test_merge_two_molecules_rejects_bad_coordinates_before_mutation();
  test_merge_two_molecules_reports_skipped_optional_numeric_descriptor();
  test_merge_two_molecules_copies_conect_without_aliasing();
  test_set_coordinates_using_mask_replaces_selected_atoms_only();
  test_set_coordinates_using_mask_rejects_before_mutation();
  test_set_coordinates_using_indices_replaces_selected_atoms_only();
  test_set_coordinates_rejects_shape_mismatch();
  test_with_coordinates_using_indices_returns_modified_copy();
  test_with_coordinates_using_mask_returns_modified_copy();
  test_with_coordinates_rejects_bad_mask_without_mutating_target();
  test_copy_conect_filters_to_selected_atoms();
  test_get_and_set_string_descriptor_using_indices();
  test_get_and_set_string_descriptor_using_mask();
  test_get_and_set_int_descriptor_using_indices();
  test_get_and_set_calc_descriptor_using_mask();
  test_descriptor_set_rejects_bad_mask_before_mutation();
  test_descriptor_get_rejects_descriptor_length_mismatch();
  test_apply_biomt_transforms_assembles_in_transform_order();
  test_biomt_transformed_returns_value_without_mutating_source();
  test_apply_biomt_transforms_from_metadata_uses_recorded_transforms();
  test_apply_biomt_transforms_rejects_invalid_inputs_without_mutation();
  return 0;
}
