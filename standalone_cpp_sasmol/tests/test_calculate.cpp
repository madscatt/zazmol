#include "sasmol/calculate.hpp"
#include "sasmol/file_io.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-6) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

void assert_close_double(sasmol::calc_type actual, sasmol::calc_type expected,
                         double tolerance = 1.0e-6) {
  assert(std::fabs(actual - expected) < tolerance);
}

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_bounds_close(const sasmol::CoordinateBounds& bounds,
                         sasmol::Vec3 expected_min,
                         sasmol::Vec3 expected_max) {
  assert_close(bounds.minimum.x, expected_min.x);
  assert_close(bounds.minimum.y, expected_min.y);
  assert_close(bounds.minimum.z, expected_min.z);
  assert_close(bounds.maximum.x, expected_max.x);
  assert_close(bounds.maximum.y, expected_max.y);
  assert_close(bounds.maximum.z, expected_max.z);
}

void assert_calc_vec_close(const sasmol::CalcVec3& actual,
                           sasmol::CalcVec3 expected,
                           double tolerance = 1.0e-5) {
  assert_close_double(actual.x, expected.x, tolerance);
  assert_close_double(actual.y, expected.y, tolerance);
  assert_close_double(actual.z, expected.z, tolerance);
}

void assert_matrix_close(const sasmol::CalcMatrix3& actual,
                         const sasmol::CalcMatrix3& expected,
                         double tolerance = 1.0e-6) {
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t column = 0; column < 3; ++column) {
      assert_close_double(actual[row][column], expected[row][column], tolerance);
    }
  }
}

double absolute_eigenvector_dot(const sasmol::CalcMatrix3& eigenvectors,
                                std::size_t column,
                                std::array<double, 3> expected) {
  double dot{};
  for (std::size_t row = 0; row < 3; ++row) {
    dot += eigenvectors[row][column] * expected[row];
  }
  return std::fabs(dot);
}

void assert_pmi_eigenpairs_are_orthonormal(
    const sasmol::PrincipalMomentsOfInertia& pmi,
    double tolerance = 1.0e-6) {
  for (std::size_t column = 0; column < 3; ++column) {
    double norm{};
    for (std::size_t row = 0; row < 3; ++row) {
      norm += pmi.eigenvectors[row][column] * pmi.eigenvectors[row][column];
    }
    assert_close_double(norm, 1.0, tolerance);

    for (std::size_t other = column + 1; other < 3; ++other) {
      double dot{};
      for (std::size_t row = 0; row < 3; ++row) {
        dot += pmi.eigenvectors[row][column] * pmi.eigenvectors[row][other];
      }
      assert_close_double(dot, 0.0, tolerance);
    }

    for (std::size_t row = 0; row < 3; ++row) {
      double left{};
      for (std::size_t k = 0; k < 3; ++k) {
        left += pmi.inertia[row][k] * pmi.eigenvectors[k][column];
      }
      const auto right = pmi.eigenvalues[column] * pmi.eigenvectors[row][column];
      assert_close_double(left, right, tolerance);
    }
  }
}

void test_calculate_minimum_and_maximum_all_loaded_frames() {
  sasmol::Molecule mol(2, 2);
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(0, 1, {-4.0F, 5.0F, 6.0F});
  mol.set_coordinate(1, 0, {7.0F, -8.0F, 9.0F});
  mol.set_coordinate(1, 1, {10.0F, 11.0F, -12.0F});

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {-4.0F, -8.0F, -12.0F},
                      {10.0F, 11.0F, 9.0F});
}

void test_calculate_principal_moments_of_inertia_synthetic() {
  sasmol::Molecule mol(3, 1);
  mol.element() = {"C", "C", "C"};
  mol.set_coordinate(0, 0, {1.0F, 0.0F, 0.0F});
  mol.set_coordinate(0, 1, {0.0F, 1.0F, 0.0F});
  mol.set_coordinate(0, 2, {0.0F, 0.0F, 1.0F});

  const auto result = sasmol::calculate_principal_moments_of_inertia(mol, 0);

  const auto mass = 12.01078;
  const sasmol::CalcMatrix3 expected_inertia{
      {{mass * 4.0 / 3.0, mass / 3.0, mass / 3.0},
       {mass / 3.0, mass * 4.0 / 3.0, mass / 3.0},
       {mass / 3.0, mass / 3.0, mass * 4.0 / 3.0}}};
  assert(!result.singular);
  assert_matrix_close(result.inertia, expected_inertia, 1.0e-6);
  assert_close_double(result.eigenvalues[0], mass, 1.0e-6);
  assert_close_double(result.eigenvalues[1], mass, 1.0e-6);
  assert_close_double(result.eigenvalues[2], 2.0 * mass, 1.0e-6);
  assert_pmi_eigenpairs_are_orthonormal(result);
}

void test_calculate_principal_moments_of_inertia_2aad_fixture() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  const auto result = sasmol::calculate_principal_moments_of_inertia(mol, 0);

  const sasmol::CalcMatrix3 expected_inertia{
      {{589.533746311651, -64.328461565902, -439.375385704833},
       {-64.328461565902, 1532.135608483722, -65.398994295632},
       {-439.375385704833, -65.398994295632, 1407.583009456943}}};
  assert(!result.singular);
  assert_matrix_close(result.inertia, expected_inertia, 1.0e-6);
  assert_close_double(result.eigenvalues[0], 391.910570545354, 1.0e-6);
  assert_close_double(result.eigenvalues[1], 1523.060334881807, 1.0e-6);
  assert_close_double(result.eigenvalues[2], 1614.281458831898, 1.0e-6);
  assert(result.eigenvalues[0] <= result.eigenvalues[1]);
  assert(result.eigenvalues[1] <= result.eigenvalues[2]);
  assert_pmi_eigenpairs_are_orthonormal(result);
  assert(absolute_eigenvector_dot(result.eigenvectors, 0,
                                  {-0.91349702, -0.07447765,
                                   -0.39997034}) > 0.999999);
  assert(absolute_eigenvector_dot(result.eigenvectors, 1,
                                  {0.22717101, -0.90894624,
                                   -0.34958556}) > 0.999999);
  assert(absolute_eigenvector_dot(result.eigenvectors, 2,
                                  {0.33751523, 0.41020703,
                                   -0.84723885}) > 0.999999);
}

void test_calculate_principal_moments_of_inertia_reports_singular_tensor() {
  sasmol::Molecule mol(1, 1);
  mol.element()[0] = "C";
  mol.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});

  const auto result = sasmol::calculate_principal_moments_of_inertia(mol, 0);

  assert(result.singular);
  assert_matrix_close(result.inertia, {{{0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0}}});
}

void test_calculate_principal_moments_of_inertia_rejects_bad_frame() {
  sasmol::Molecule mol(1, 1);
  mol.element()[0] = "C";
  bool threw = false;

  try {
    (void)sasmol::calculate_principal_moments_of_inertia(mol, 1);
  } catch (const std::out_of_range&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_mass_2aad_fixture() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  const auto result = sasmol::calculate_mass(mol);

  const std::vector<sasmol::calc_type> expected_mass = {
      14.00672, 12.01078, 12.01078, 15.99943, 12.01078,
      12.01078, 12.01078, 12.01078, 14.00672, 12.01078,
      12.01078, 15.99943, 12.01078, 15.99943, 12.01078};
  assert(result.ok());
  assert(mol.mass().size() == expected_mass.size());
  for (std::size_t atom = 0; atom < expected_mass.size(); ++atom) {
    assert_close_double(mol.mass()[atom], expected_mass[atom]);
  }
  assert_close_double(result.total_mass, 196.11953);
  assert_close_double(mol.total_mass(), result.total_mass);
}

void test_calculate_mass_rna_and_1crn_fixture_totals() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "rna.pdb"), mol);
  assert(status.ok());

  auto result = sasmol::calculate_mass(mol);
  assert(result.ok());
  assert_close_double(result.total_mass, 106197.087, 0.001);

  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());
  result = sasmol::calculate_mass(mol);
  assert(result.ok());
  assert_close_double(result.total_mass, 4412.904, 0.001);
}

void test_calculate_mass_reports_unknown_elements() {
  sasmol::Molecule mol(2, 1);
  mol.element()[0] = "C";
  mol.element()[1] = "XX";

  const auto result = sasmol::calculate_mass(mol);

  assert(!result.ok());
  assert((result.unknown_elements == std::vector<std::string>{"XX"}));
  assert_close_double(result.total_mass, 12.01078);
  assert_close_double(mol.mass()[0], 12.01078);
  assert_close_double(mol.mass()[1], 0.0);
  assert_close_double(mol.total_mass(), 12.01078);
}

void test_calculate_mass_rejects_element_mismatch() {
  sasmol::Molecule mol(1, 1);
  mol.element().clear();
  bool threw = false;

  try {
    (void)sasmol::calculate_mass(mol);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_molecular_formula_fixtures() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1ATM.pdb"), mol);
  assert(status.ok());

  auto formula = sasmol::calculate_molecular_formula(mol);
  assert((formula == std::map<std::string, std::size_t>{{"N", 1}}));
  assert(mol.formula() == formula);

  status = reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());
  formula = sasmol::calculate_molecular_formula(mol);
  assert(formula.at("N") == 2);
  assert(formula.at("O") == 3);
  assert(formula.at("C") == 10);
  assert(mol.formula() == formula);
}

void test_calculate_molecular_formula_rejects_element_mismatch() {
  sasmol::Molecule mol(1, 1);
  mol.element().clear();
  bool threw = false;

  try {
    (void)sasmol::calculate_molecular_formula(mol);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_residue_charge_two_residues() {
  sasmol::Molecule mol(4, 1);
  mol.resid() = {1, 1, 2, 2};
  mol.atom_charge() = {0.1, 0.2, -0.4, 0.3};

  sasmol::calculate_residue_charge(mol);

  assert(mol.residue_charge().size() == 4);
  assert_close_double(mol.residue_charge()[0], 0.3);
  assert_close_double(mol.residue_charge()[1], 0.3);
  assert_close_double(mol.residue_charge()[2], -0.1);
  assert_close_double(mol.residue_charge()[3], -0.1);
}

void test_calculate_residue_charge_rejects_descriptor_mismatch() {
  sasmol::Molecule mol(1, 1);
  mol.atom_charge().clear();
  bool threw = false;

  try {
    sasmol::calculate_residue_charge(mol);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_center_of_mass_autocalculates_mass() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());
  mol.set_total_mass(0.0);

  const auto center = sasmol::calculate_center_of_mass(mol, 0);

  assert_calc_vec_close(center, {75.68045, 43.70790, 41.27621}, 1.0e-5);
  assert(mol.total_mass() > 0.0);
}

void test_calculate_center_of_mass_fixture_totals() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "rna.pdb"), mol);
  assert(status.ok());
  auto center = sasmol::calculate_center_of_mass(mol, 0);
  assert_calc_vec_close(center, {-8.033, 4.352, 9.231}, 0.001);

  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());
  center = sasmol::calculate_center_of_mass(mol, 0);
  assert_calc_vec_close(center, {9.30026, 9.77488, 6.97776}, 1.0e-5);
}

void test_calculate_center_of_mass_rejects_unknown_mass() {
  sasmol::Molecule mol(1, 1);
  mol.element()[0] = "XX";
  bool threw = false;

  try {
    (void)sasmol::calculate_center_of_mass(mol, 0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_center_of_mass_rejects_bad_frame() {
  sasmol::Molecule mol(1, 1);
  mol.element()[0] = "C";
  bool threw = false;

  try {
    (void)sasmol::calculate_center_of_mass(mol, 1);
  } catch (const std::out_of_range&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_radius_of_gyration_fixtures() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_pdb(fixture_path("pdb_common", "1ATM.pdb"), mol);
  assert(status.ok());
  assert_close_double(sasmol::calculate_radius_of_gyration(mol, 0), 0.0);

  status = reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());
  assert_close_double(sasmol::calculate_radius_of_gyration(mol, 0), 2.998744,
                      1.0e-5);

  status = reader.read_pdb(fixture_path("pdb_common", "rna.pdb"), mol);
  assert(status.ok());
  assert_close_double(sasmol::calculate_radius_of_gyration(mol, 0),
                      65.62305792313865, 0.001);

  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());
  assert_close_double(sasmol::calculate_radius_of_gyration(mol, 0),
                      9.666350273306351, 0.001);
}

void test_calculate_radius_of_gyration_reuses_center_errors() {
  sasmol::Molecule mol(1, 1);
  mol.element()[0] = "XX";
  bool threw = false;

  try {
    (void)sasmol::calculate_radius_of_gyration(mol, 0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_root_mean_square_deviation_synthetic() {
  sasmol::Molecule first(1, 1);
  sasmol::Molecule second(1, 1);
  first.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  second.set_coordinate(0, 0, {4.0F, 5.0F, 6.0F});

  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(first, second),
      3.0 * std::sqrt(3.0));

  first.resize(2, 1);
  second.resize(2, 1);
  first.set_coordinate(0, 0, {7.0F, 8.0F, 9.0F});
  first.set_coordinate(0, 1, {1.0F, 3.0F, 5.0F});
  second.set_coordinate(0, 0, {12.0F, 53.0F, 67.0F});
  second.set_coordinate(0, 1, {76.0F, 87.0F, 96.0F});

  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(first, second),
      114.83901775964473);
}

void test_calculate_root_mean_square_deviation_fixtures() {
  sasmol::PdbReader reader;
  sasmol::Molecule first;
  sasmol::Molecule second;

  auto status = reader.read_pdb(fixture_path("pdb_common", "1ATM.pdb"), first);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("pdb_common", "1ATM.pdb"), second);
  assert(status.ok());
  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(first, second), 0.0);

  status = reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), first);
  assert(status.ok());
  status = reader.read_pdb(fixture_path("sasmol/calculate", "1CRN-rot.pdb"),
                           second);
  assert(status.ok());
  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(first, second), 29.008,
      0.001);

  status = reader.read_pdb(fixture_path("sasmol/calculate", "1CRN-rot-shift.pdb"),
                           second);
  assert(status.ok());
  assert_close_double(
      sasmol::calculate_root_mean_square_deviation(first, second), 19.831,
      0.001);
}

void test_calculate_root_mean_square_deviation_rejects_shape_mismatch() {
  const sasmol::Molecule first(1, 1);
  const sasmol::Molecule second(2, 1);
  bool threw = false;

  try {
    (void)sasmol::calculate_root_mean_square_deviation(first, second);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_minimum_and_maximum_selected_frames() {
  sasmol::Molecule mol(1, 3);
  mol.set_coordinate(0, 0, {-10.0F, -10.0F, -10.0F});
  mol.set_coordinate(1, 0, {1.0F, 2.0F, 3.0F});
  mol.set_coordinate(2, 0, {4.0F, 5.0F, 6.0F});

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol, {1, 2});

  assert_bounds_close(bounds, {1.0F, 2.0F, 3.0F}, {4.0F, 5.0F, 6.0F});
}

void test_calculate_minimum_and_maximum_rejects_empty_molecule() {
  const sasmol::Molecule mol;
  bool threw = false;

  try {
    (void)sasmol::calculate_minimum_and_maximum(mol);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_minimum_and_maximum_rejects_bad_frame() {
  const sasmol::Molecule mol(1, 1);
  bool threw = false;

  try {
    (void)sasmol::calculate_minimum_and_maximum(mol, {1});
  } catch (const std::out_of_range&) {
    threw = true;
  }

  assert(threw);
}

void test_calculate_minimum_and_maximum_pdb_fixture_2aad() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), mol);
  assert(status.ok());

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {70.721F, 41.799F, 39.354F},
                      {79.712F, 46.273F, 43.910F});
}

void test_calculate_minimum_and_maximum_pdb_fixture_1crn() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), mol);
  assert(status.ok());

  const auto bounds = sasmol::calculate_minimum_and_maximum(mol);

  assert_bounds_close(bounds, {-3.097F, -0.516F, -7.422F},
                      {24.284F, 20.937F, 19.580F});
}

void test_calculate_minimum_and_maximum_all_steps_alias() {
  sasmol::PdbReader reader;
  sasmol::Molecule mol;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1ATM-1to2.pdb"), mol);
  assert(status.ok());

  const auto primary = sasmol::calculate_minimum_and_maximum_all_steps(mol);
  const auto alias = sasmol::calc_minmax_all_steps(mol);

  assert_bounds_close(primary, {73.944F, 38.799F, 41.652F},
                      {76.944F, 41.799F, 41.652F});
  assert_bounds_close(alias, primary.minimum, primary.maximum);
}

void test_calculate_minimum_and_maximum_all_steps_streams_dcd() {
  const auto bounds = sasmol::calculate_minimum_and_maximum_all_steps(
      fixture_path("dcd_common", "1ATM.dcd"));
  const auto alias =
      sasmol::calc_minmax_all_steps(fixture_path("dcd_common", "1ATM.dcd"));

  assert_bounds_close(bounds, {73.944F, 38.799F, 41.652F},
                      {76.944F, 41.799F, 41.652F});
  assert_bounds_close(alias, bounds.minimum, bounds.maximum);
}

void test_calculate_minimum_and_maximum_dcd_matches_loaded_trajectory() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;
  const auto path = fixture_path("dcd_common", "2AAD.dcd");
  const auto status = reader.read_dcd(path, mol);
  assert(status.ok());

  const auto loaded = sasmol::calculate_minimum_and_maximum(mol);
  const auto streamed = sasmol::calculate_minimum_and_maximum_all_steps(path);

  assert_bounds_close(streamed, loaded.minimum, loaded.maximum);
}

}  // namespace

int main() {
  test_calculate_mass_2aad_fixture();
  test_calculate_mass_rna_and_1crn_fixture_totals();
  test_calculate_mass_reports_unknown_elements();
  test_calculate_mass_rejects_element_mismatch();
  test_calculate_molecular_formula_fixtures();
  test_calculate_molecular_formula_rejects_element_mismatch();
  test_calculate_residue_charge_two_residues();
  test_calculate_residue_charge_rejects_descriptor_mismatch();
  test_calculate_center_of_mass_autocalculates_mass();
  test_calculate_center_of_mass_fixture_totals();
  test_calculate_center_of_mass_rejects_unknown_mass();
  test_calculate_center_of_mass_rejects_bad_frame();
  test_calculate_radius_of_gyration_fixtures();
  test_calculate_radius_of_gyration_reuses_center_errors();
  test_calculate_root_mean_square_deviation_synthetic();
  test_calculate_root_mean_square_deviation_fixtures();
  test_calculate_root_mean_square_deviation_rejects_shape_mismatch();
  test_calculate_principal_moments_of_inertia_synthetic();
  test_calculate_principal_moments_of_inertia_2aad_fixture();
  test_calculate_principal_moments_of_inertia_reports_singular_tensor();
  test_calculate_principal_moments_of_inertia_rejects_bad_frame();
  test_calculate_minimum_and_maximum_all_loaded_frames();
  test_calculate_minimum_and_maximum_selected_frames();
  test_calculate_minimum_and_maximum_rejects_empty_molecule();
  test_calculate_minimum_and_maximum_rejects_bad_frame();
  test_calculate_minimum_and_maximum_pdb_fixture_2aad();
  test_calculate_minimum_and_maximum_pdb_fixture_1crn();
  test_calculate_minimum_and_maximum_all_steps_alias();
  test_calculate_minimum_and_maximum_all_steps_streams_dcd();
  test_calculate_minimum_and_maximum_dcd_matches_loaded_trajectory();
  return 0;
}
