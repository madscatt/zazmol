#include "sasmol/file_io.hpp"
#include "sasmol/selection.hpp"

#include <cassert>
#include <filesystem>

namespace {

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

void test_indices_all() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result = sasmol::indices_all(mol);

  assert(result.ok());
  assert(result.indices.size() == mol.natoms());
  assert(result.indices.front() == 0);
  assert(result.indices.back() == mol.natoms() - 1);
}

void test_explicit_helpers() {
  const auto mol = read_fixture("1CRN.pdb");

  auto result = sasmol::indices_by_name(mol, "CA");
  assert(result.ok());
  assert(result.indices.size() == 46);

  result = sasmol::indices_by_resid_range(mol, 20, 31);
  assert(result.ok());
  assert(!result.indices.empty());
  for (const auto atom : result.indices) {
    assert(mol.resid()[atom] >= 20);
    assert(mol.resid()[atom] <= 31);
  }
}

void test_select_all_keyword() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result = sasmol::select_indices(mol, "all");

  assert(result.ok());
  assert(result.indices.size() == 15);
}

void test_safe_expression_matches_align_basis_example() {
  const auto mol = read_fixture("1CRN.pdb");

  const auto result = sasmol::select_indices(
      mol, "name[i] == \"CA\" and (resid[i] >= 20 and resid[i] <= 31)");

  assert(result.ok());
  assert(result.indices.size() == 12);
  for (const auto atom : result.indices) {
    assert(mol.name()[atom] == "CA");
    assert(mol.resid()[atom] >= 20);
    assert(mol.resid()[atom] <= 31);
  }
}

void test_safe_expression_supports_or_and_not_equal() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto result = sasmol::select_indices(
      mol, "element[i] == \"N\" or name[i] != \"CA\"");

  assert(result.ok());
  assert(result.indices.size() > 2);
}

void test_safe_expression_supports_unary_not() {
  const auto mol = read_fixture("2AAD.pdb");

  const auto all = sasmol::select_indices(mol, "not name[i] == None");
  assert(all.ok());
  assert(all.indices.size() == mol.natoms());

  const auto not_ca = sasmol::select_indices(mol, "not name[i] == \"CA\"");
  assert(not_ca.ok());
  for (const auto atom : not_ca.indices) {
    assert(mol.name()[atom] != "CA");
  }
}

void test_safe_expression_supports_string_offset_for_heavy_atoms() {
  sasmol::Molecule mol(4);
  mol.name() = {"H1", "CA", "N", ""};

  const auto result = sasmol::select_indices(mol, "not name[i][0] == \"H\"");

  assert(result.ok());
  assert((result.indices == std::vector<std::size_t>{1, 2, 3}));
}

void test_safe_expression_supports_additional_python_mask_descriptors() {
  sasmol::Molecule mol(3);
  mol.loc() = {" ", "A", " "};
  mol.residue_flag() = {0, 1, 1};

  auto result = sasmol::select_indices(mol, "loc[i] == \" \"");
  assert(result.ok());
  assert((result.indices == std::vector<std::size_t>{0, 2}));

  result = sasmol::select_indices(mol, "residue_flag[i] == True");
  assert(result.ok());
  assert((result.indices == std::vector<std::size_t>{1, 2}));
}

void test_unsupported_expression_fails_without_indices() {
  const auto mol = read_fixture("1CRN.pdb");

  const auto result =
      sasmol::select_indices(mol, "name[i].startswith(\"C\")");

  assert(!result.ok());
  assert(result.indices.empty());
}

void test_unsupported_descriptor_fails_without_indices() {
  const auto mol = read_fixture("1CRN.pdb");

  const auto result = sasmol::select_indices(mol, "foobar[i] == \"CA\"");

  assert(!result.ok());
  assert(result.indices.empty());
}

void test_no_match_is_error_without_indices() {
  const auto mol = read_fixture("1CRN.pdb");

  const auto result = sasmol::select_indices(mol, "name[i] == \"NOPE\"");

  assert(!result.ok());
  assert(result.indices.empty());
}

void test_descriptor_length_mismatch_fails_without_indices() {
  auto mol = read_fixture("1CRN.pdb");
  mol.name().clear();

  const auto result = sasmol::select_indices(mol, "name[i] == \"CA\"");

  assert(!result.ok());
  assert(result.indices.empty());
}

}  // namespace

int main() {
  test_indices_all();
  test_explicit_helpers();
  test_select_all_keyword();
  test_safe_expression_matches_align_basis_example();
  test_safe_expression_supports_or_and_not_equal();
  test_safe_expression_supports_unary_not();
  test_safe_expression_supports_string_offset_for_heavy_atoms();
  test_safe_expression_supports_additional_python_mask_descriptors();
  test_unsupported_expression_fails_without_indices();
  test_unsupported_descriptor_fails_without_indices();
  test_no_match_is_error_without_indices();
  test_descriptor_length_mismatch_fails_without_indices();
  return 0;
}
