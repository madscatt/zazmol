#include "sasmol/file_io.hpp"

#include <cassert>
#include <string>
#include <vector>

namespace {

void test_empty_when_no_biomt_records() {
  const std::vector<std::string> header{
      "HEADER    TEST",
      "ATOM      1  N   ALA A   1      10.000  10.000  10.000",
  };
  const auto biomt = sasmol::parse_biomt_header_records(header);
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
  const auto biomt = sasmol::parse_biomt_header_records(header);

  assert(biomt.size() == 1);
  const auto& record = biomt.at(1);
  assert(record.subdivs.size() == 1);
  assert(record.subdivs[0] == "N");
  assert(record.auth_bio_unit == "MONOMERIC");
  assert(record.rot.size() == 1);
  assert(record.trans.size() == 1);
  assert(record.rot[0][0][0] == static_cast<sasmol::calc_type>(1.0));
  assert(record.rot[0][1][1] == static_cast<sasmol::calc_type>(1.0));
  assert(record.rot[0][2][2] == static_cast<sasmol::calc_type>(1.0));
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
  const auto biomt = sasmol::parse_biomt_header_records(header);

  assert(biomt.size() == 1);
  const auto& record = biomt.at(1);
  assert(record.soft_bio_unit == "DIMERIC");
  assert((record.subdivs == std::vector<std::string>{"A", "B", "C"}));
  assert(record.rot.size() == 2);
  assert(record.trans.size() == 2);
  assert(record.trans[1][0] == static_cast<sasmol::calc_type>(61.157));
  assert(record.trans[1][2] == static_cast<sasmol::calc_type>(61.5385));
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
  const auto biomt = sasmol::parse_biomt_header_records(header);

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
