#include "sasmol/properties.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>

namespace {

std::filesystem::path fixture_path(const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / "sasmol" /
         "properties" / name;
}

std::vector<sasmol::calc_type> parse_value_list(std::string values) {
  for (auto& character : values) {
    if (character == '[' || character == ']' || character == ',') {
      character = ' ';
    }
  }

  std::istringstream stream(values);
  std::vector<sasmol::calc_type> parsed;
  sasmol::calc_type value{};
  while (stream >> value) {
    parsed.push_back(value);
  }
  return parsed;
}

sasmol::PropertyTable read_property_table(const char* name) {
  std::ifstream input(fixture_path(name));
  assert(input);

  sasmol::PropertyTable table;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream stream(line);
    std::string key;
    stream >> key;
    std::string values;
    std::getline(stream, values);
    table[key] = parse_value_list(values);
  }
  return table;
}

sasmol::OptionalPropertyTable read_optional_property_table(const char* name) {
  std::ifstream input(fixture_path(name));
  assert(input);

  sasmol::OptionalPropertyTable table;
  std::string key;
  std::string value;
  while (input >> key >> value) {
    if (value == "None") {
      table[key] = std::nullopt;
    } else {
      table[key] = std::stod(value);
    }
  }
  return table;
}

void assert_close(sasmol::calc_type actual, sasmol::calc_type expected) {
  assert(std::fabs(actual - expected) < 1.0e-12);
}

void assert_table_matches_fixture(const sasmol::PropertyTable& actual,
                                  const char* fixture_name) {
  const auto expected = read_property_table(fixture_name);
  assert(actual.size() == expected.size());
  for (const auto& [key, expected_values] : expected) {
    const auto found = actual.find(key);
    assert(found != actual.end());
    assert(found->second.size() == expected_values.size());
    for (std::size_t i = 0; i < expected_values.size(); ++i) {
      assert_close(found->second[i], expected_values[i]);
    }
  }
}

void assert_optional_table_matches_fixture(
    const sasmol::OptionalPropertyTable& actual, const char* fixture_name) {
  const auto expected = read_optional_property_table(fixture_name);
  assert(actual.size() == expected.size());
  for (const auto& [key, expected_value] : expected) {
    const auto found = actual.find(key);
    assert(found != actual.end());
    assert(found->second.has_value() == expected_value.has_value());
    if (expected_value.has_value()) {
      assert_close(*found->second, *expected_value);
    }
  }
}

void test_property_tables_match_python_fixtures() {
  assert_table_matches_fixture(sasmol::amino_acid_sld(), "amino_acid_sld.txt");
  assert_table_matches_fixture(sasmol::element_scattering_lengths(),
                               "element_scattering_lengths.txt");
  assert_table_matches_fixture(sasmol::nucleotide_scattering_lengths(),
                               "nucleotide_scattering_lengths.txt");
  assert_table_matches_fixture(sasmol::dna_scattering_lengths(),
                               "dna_scattering_lengths.txt");
  assert_table_matches_fixture(sasmol::rna_scattering_lengths(),
                               "rna_scattering_lengths.txt");
  assert_table_matches_fixture(sasmol::protein_scattering_lengths(),
                               "protein_scattering_lengths.txt");
  assert_optional_table_matches_fixture(sasmol::van_der_waals_radii(),
                                        "van_der_waals_radii.txt");
}

void test_aliases_return_canonical_tables() {
  assert(&sasmol::element_sl() == &sasmol::element_scattering_lengths());
  assert(&sasmol::nucleotide_sl() ==
         &sasmol::nucleotide_scattering_lengths());
  assert(&sasmol::dna_sl() == &sasmol::dna_scattering_lengths());
  assert(&sasmol::rna_sl() == &sasmol::rna_scattering_lengths());
  assert(&sasmol::protein_sl() == &sasmol::protein_scattering_lengths());
  assert(&sasmol::van_der_Waals_radii() == &sasmol::van_der_waals_radii());
}

void test_vdw_preserves_python_none_values() {
  const auto& vdw = sasmol::van_der_waals_radii();
  assert(vdw.at("CE") == std::nullopt);
  assert(vdw.at("YB") == std::nullopt);
  assert_close(*vdw.at("H"), 1.2);
}

}  // namespace

int main() {
  test_property_tables_match_python_fixtures();
  test_aliases_return_canonical_tables();
  test_vdw_preserves_python_none_values();
  return 0;
}
