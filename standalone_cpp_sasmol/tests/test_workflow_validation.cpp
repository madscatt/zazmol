#include "sasmol/calculate.hpp"
#include "sasmol/file_io.hpp"
#include "sasmol/selection.hpp"
#include "sasmol/subset.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <vector>

namespace {

std::filesystem::path fixture_path(const char* area, const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / area / name;
}

void assert_close(float actual, float expected, double tolerance = 0.001) {
  assert(std::fabs(static_cast<double>(actual - expected)) <= tolerance);
}

void test_sassie_selection_copy_dcd_round_trip_workflow() {
  sasmol::PdbReader pdb_reader;
  sasmol::Molecule source;
  auto status = pdb_reader.read_pdb(fixture_path("pdb_common", "2AAD.pdb"), source);
  assert(status.ok());

  const auto selected = sasmol::select_sassie_basis_mask(
      source, "name N or name CA or name C or name O");
  assert(selected.ok());

  sasmol::Molecule subset;
  const auto copy_status =
      sasmol::copy_molecule_using_mask(source, subset, selected.mask, 0);
  assert(copy_status.ok());
  assert(subset.natoms() == 8);
  for (const auto& name : subset.name()) {
    assert(name == "N" || name == "CA" || name == "C" || name == "O");
  }

  const auto output =
      std::filesystem::temp_directory_path() / "sasmol_workflow_subset.dcd";
  sasmol::DcdWriter writer;
  status = writer.write_dcd(output, subset);
  assert(status.ok());

  sasmol::DcdReader dcd_reader;
  sasmol::Molecule round_trip;
  status = dcd_reader.read_dcd(output, round_trip);
  assert(status.ok());
  assert(round_trip.natoms() == subset.natoms());
  assert(round_trip.number_of_frames() == 1);
  for (std::size_t atom = 0; atom < subset.natoms(); ++atom) {
    const auto expected = subset.coordinate(0, atom);
    const auto actual = round_trip.coordinate(0, atom);
    assert_close(actual.x, expected.x);
    assert_close(actual.y, expected.y);
    assert_close(actual.z, expected.z);
  }

  std::filesystem::remove(output);
}

void test_dcd_streaming_bounds_match_loaded_workflow() {
  const auto trajectory = fixture_path("dcd_common", "2AAD.dcd");

  sasmol::DcdReader reader;
  sasmol::Molecule loaded;
  const auto status = reader.read_dcd(trajectory, loaded);
  assert(status.ok());

  const auto loaded_bounds = sasmol::calculate_minimum_and_maximum(loaded);
  const auto streamed_bounds =
      sasmol::calculate_minimum_and_maximum_all_steps(trajectory);

  assert_close(streamed_bounds.minimum.x, loaded_bounds.minimum.x);
  assert_close(streamed_bounds.minimum.y, loaded_bounds.minimum.y);
  assert_close(streamed_bounds.minimum.z, loaded_bounds.minimum.z);
  assert_close(streamed_bounds.maximum.x, loaded_bounds.maximum.x);
  assert_close(streamed_bounds.maximum.y, loaded_bounds.maximum.y);
  assert_close(streamed_bounds.maximum.z, loaded_bounds.maximum.z);
}

void test_pdb_selection_calculation_workflow() {
  sasmol::PdbReader reader;
  sasmol::Molecule source;
  const auto status =
      reader.read_pdb(fixture_path("pdb_common", "1CRN.pdb"), source);
  assert(status.ok());

  const auto selected = sasmol::select_sassie_basis_mask(
      source, "resid >= 20 and resid <= 31 and name CA");
  assert(selected.ok());

  sasmol::Molecule subset;
  const auto copy_status =
      sasmol::copy_molecule_using_mask(source, subset, selected.mask, 0);
  assert(copy_status.ok());
  assert(subset.natoms() == 12);

  const auto mass = sasmol::calculate_mass(subset);
  assert(mass.ok());
  const auto center = sasmol::calculate_center_of_mass(subset, 0);
  const auto rg = sasmol::calculate_radius_of_gyration(subset, 0);

  assert(center.x > 0.0);
  assert(center.y > 0.0);
  assert(center.z > 0.0);
  assert(rg > 0.0);
}

}  // namespace

int main() {
  test_sassie_selection_copy_dcd_round_trip_workflow();
  test_dcd_streaming_bounds_match_loaded_workflow();
  test_pdb_selection_calculation_workflow();
  return 0;
}
