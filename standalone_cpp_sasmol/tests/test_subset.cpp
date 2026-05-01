#include "sasmol/file_io.hpp"
#include "sasmol/selection.hpp"
#include "sasmol/subset.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>

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

void test_copied_molecule_using_indices_returns_value() {
  const auto source = read_fixture("2AAD.pdb");

  const auto subset = sasmol::copied_molecule_using_indices(source, {0, 1}, 0);

  assert(subset.natoms() == 2);
  assert(subset.index()[0] == source.index()[0]);
  assert(subset.index()[1] == source.index()[1]);
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

}  // namespace

int main() {
  test_get_coordinates_using_indices();
  test_get_coordinates_rejects_bad_index();
  test_get_indices_from_mask();
  test_get_coordinates_using_mask_delegates_to_indices();
  test_mask_rejects_bad_shape_and_values();
  test_copy_molecule_using_indices_preserves_descriptors();
  test_copy_molecule_using_mask_preserves_descriptors();
  test_copied_molecule_using_indices_returns_value();
  test_set_coordinates_using_mask_replaces_selected_atoms_only();
  test_set_coordinates_using_mask_rejects_before_mutation();
  test_set_coordinates_using_indices_replaces_selected_atoms_only();
  test_set_coordinates_rejects_shape_mismatch();
  test_copy_conect_filters_to_selected_atoms();
  return 0;
}
