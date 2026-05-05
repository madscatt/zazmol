#include <array>
#include <cassert>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct FixtureBiomtRecord {
  std::vector<std::string> subdivs;
  std::string auth_bio_unit;
  std::string soft_bio_unit;
  std::vector<std::array<std::array<double, 3>, 3>> rot;
  std::vector<std::array<double, 3>> trans;
};

using FixtureBiomtMap = std::map<int, FixtureBiomtRecord>;

FixtureBiomtRecord empty_record() {
  return FixtureBiomtRecord{};
}

std::vector<std::string> split_commas(const std::string& text) {
  std::vector<std::string> tokens;
  std::stringstream ss(text);
  std::string token;
  while (std::getline(ss, token, ',')) {
    const auto first = token.find_first_not_of(' ');
    if (first == std::string::npos) {
      continue;
    }
    const auto last = token.find_last_not_of(' ');
    tokens.push_back(token.substr(first, last - first + 1));
  }
  return tokens;
}

FixtureBiomtMap parse_biomt_fixture_contract(
    const std::vector<std::string>& header_lines) {
  FixtureBiomtMap biomt;
  int current_biomol_no = -1;
  int current_transform_id = -1;
  std::array<std::array<double, 3>, 3> current_rot{};
  std::array<double, 3> current_trans{};
  int next_row = 1;

  for (const auto& line : header_lines) {
    if (line.rfind("REMARK 350", 0) != 0) {
      continue;
    }
    if (line.size() < 11) {
      continue;
    }

    auto text = line.substr(10);
    const auto first = text.find_first_not_of(' ');
    if (first == std::string::npos) {
      continue;
    }
    text = text.substr(first);

    if (text.rfind("BIOMOLECULE:", 0) == 0) {
      current_biomol_no = -1;
      current_transform_id = -1;
      next_row = 1;

      const auto colon = text.find(':');
      if (colon == std::string::npos) {
        continue;
      }
      const auto candidates = split_commas(text.substr(colon + 1));
      for (const auto& token : candidates) {
        try {
          current_biomol_no = std::stoi(token);
          break;
        } catch (...) {
          continue;
        }
      }
      if (current_biomol_no != -1 && biomt.find(current_biomol_no) == biomt.end()) {
        biomt[current_biomol_no] = empty_record();
      }
      continue;
    }

    if (current_biomol_no == -1) {
      continue;
    }

    auto& record = biomt[current_biomol_no];

    if (text.rfind("AUTHOR DETERMINED BIOLOGICAL UNIT:", 0) == 0) {
      const auto colon = text.find(':');
      if (colon != std::string::npos) {
        record.auth_bio_unit = text.substr(colon + 1);
        const auto trim = record.auth_bio_unit.find_first_not_of(' ');
        record.auth_bio_unit = trim == std::string::npos ? "" : record.auth_bio_unit.substr(trim);
      }
      continue;
    }

    if (text.rfind("SOFTWARE DETERMINED", 0) == 0) {
      const auto colon = text.find(':');
      if (colon != std::string::npos) {
        record.soft_bio_unit = text.substr(colon + 1);
        const auto trim = record.soft_bio_unit.find_first_not_of(' ');
        record.soft_bio_unit = trim == std::string::npos ? "" : record.soft_bio_unit.substr(trim);
      }
      continue;
    }

    if (text.rfind("APPLY THE FOLLOWING TO CHAINS:", 0) == 0 ||
        text.rfind("AND CHAINS:", 0) == 0) {
      const auto colon = text.find(':');
      if (colon != std::string::npos) {
        for (const auto& chain : split_commas(text.substr(colon + 1))) {
          bool seen = false;
          for (const auto& existing : record.subdivs) {
            if (existing == chain) {
              seen = true;
              break;
            }
          }
          if (!seen) {
            record.subdivs.push_back(chain);
          }
        }
      }
      continue;
    }

    if (text.rfind("BIOMT", 0) != 0) {
      continue;
    }

    std::stringstream cols(text);
    std::string row_token;
    std::string transform_id_token;
    std::string c1;
    std::string c2;
    std::string c3;
    std::string t;
    cols >> row_token >> transform_id_token >> c1 >> c2 >> c3 >> t;
    if (row_token.empty() || transform_id_token.empty() || c1.empty() || c2.empty() ||
        c3.empty() || t.empty()) {
      continue;
    }
    if (row_token.size() < 6) {
      continue;
    }

    int row = 0;
    int transform_id = -1;
    double r0 = 0.0;
    double r1 = 0.0;
    double r2 = 0.0;
    double trans = 0.0;
    try {
      row = std::stoi(row_token.substr(row_token.size() - 1));
      transform_id = std::stoi(transform_id_token);
      r0 = std::stod(c1);
      r1 = std::stod(c2);
      r2 = std::stod(c3);
      trans = std::stod(t);
    } catch (...) {
      continue;
    }

    if (row < 1 || row > 3) {
      continue;
    }

    if (row == 1) {
      current_transform_id = transform_id;
      current_rot = {};
      current_trans = {};
      next_row = 1;
    }

    if (current_transform_id == -1) {
      continue;
    }
    if (transform_id != current_transform_id) {
      if (row != 1) {
        continue;
      }
      current_transform_id = transform_id;
      current_rot = {};
      current_trans = {};
      next_row = 1;
    }
    if (row != next_row) {
      if (row != 1) {
        continue;
      }
      current_transform_id = transform_id;
      current_rot = {};
      current_trans = {};
      next_row = 1;
    }

    current_rot[static_cast<std::size_t>(row - 1)] = {r0, r1, r2};
    current_trans[static_cast<std::size_t>(row - 1)] = trans;

    if (row == 3) {
      record.rot.push_back(current_rot);
      record.trans.push_back(current_trans);
      current_transform_id = -1;
      next_row = 1;
    } else {
      next_row = row + 1;
    }
  }

  return biomt;
}

void test_empty_when_no_biomt_records() {
  const std::vector<std::string> header{
      "HEADER    TEST",
      "ATOM      1  N   ALA A   1      10.000  10.000  10.000",
  };
  const auto biomt = parse_biomt_fixture_contract(header);
  assert(biomt.empty());
}

void test_single_operator_metadata_record() {
  const std::vector<std::string> header{
      "REMARK 350 BIOMOLECULE: 1",
      "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC",
      "REMARK 350 APPLY THE FOLLOWING TO CHAINS: N",
      "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
      "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
      "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
  };
  const auto biomt = parse_biomt_fixture_contract(header);

  assert(biomt.size() == 1);
  const auto& record = biomt.at(1);
  assert(record.subdivs.size() == 1);
  assert(record.subdivs[0] == "N");
  assert(record.auth_bio_unit == "MONOMERIC");
  assert(record.rot.size() == 1);
  assert(record.trans.size() == 1);
  assert(record.rot[0][0][0] == 1.0);
  assert(record.rot[0][1][1] == 1.0);
  assert(record.rot[0][2][2] == 1.0);
}

void test_chain_continuation_and_multiple_operators() {
  const std::vector<std::string> header{
      "REMARK 350 BIOMOLECULE: 1",
      "REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC",
      "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B",
      "REMARK 350                    AND CHAINS: C",
      "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
      "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
      "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
      "REMARK 350   BIOMT1   2  0.000000 -1.000000  0.000000       61.15700",
      "REMARK 350   BIOMT2   2 -1.000000  0.000000  0.000000       61.15700",
      "REMARK 350   BIOMT3   2  0.000000  0.000000 -1.000000       61.53850",
  };
  const auto biomt = parse_biomt_fixture_contract(header);

  assert(biomt.size() == 1);
  const auto& record = biomt.at(1);
  assert(record.soft_bio_unit == "DIMERIC");
  assert((record.subdivs == std::vector<std::string>{"A", "B", "C"}));
  assert(record.rot.size() == 2);
  assert(record.trans.size() == 2);
  assert(record.trans[1][0] == 61.157);
  assert(record.trans[1][2] == 61.5385);
}

void test_incomplete_triplet_is_ignored() {
  const std::vector<std::string> header{
      "REMARK 350 BIOMOLECULE: 1",
      "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A",
      "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
      "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
      "REMARK 350 BIOMOLECULE: 2",
      "REMARK 350 APPLY THE FOLLOWING TO CHAINS: B",
      "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000",
      "REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000",
      "REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000",
  };
  const auto biomt = parse_biomt_fixture_contract(header);

  assert(biomt.size() == 2);
  assert(biomt.at(1).rot.empty());
  assert(biomt.at(1).trans.empty());
  assert(biomt.at(2).rot.size() == 1);
  assert(biomt.at(2).trans.size() == 1);
}

}  // namespace

int main() {
  test_empty_when_no_biomt_records();
  test_single_operator_metadata_record();
  test_chain_continuation_and_multiple_operators();
  test_incomplete_triplet_is_ignored();
  return 0;
}
