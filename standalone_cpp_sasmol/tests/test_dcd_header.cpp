#include "sasmol/file_io.hpp"

#include <cassert>
#include <filesystem>

namespace {

std::filesystem::path fixture_path(const char* name) {
  return std::filesystem::path(SASMOL_TEST_DATA_DIR) / "dcd_common" / name;
}

void assert_header(const char* filename, std::size_t natoms,
                   std::size_t nframes, int charmm_flags) {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(fixture_path(filename));
  assert(status.ok());
  assert(reader.is_open());

  status = reader.read_header(header);
  assert(status.ok());
  assert(header.natoms == natoms);
  assert(header.nframes == nframes);
  assert(header.reverse_endian == false);
  assert(header.charmm_format == true);
  assert(header.has_unit_cell == true);
  assert(header.charmm_flags == charmm_flags);
  assert(header.namnf == 0);
  assert(header.nsavc == 1);

  status = reader.close_dcd_read();
  assert(status.ok());
}

void test_small_fixture_headers() {
  assert_header("1ATM.dcd", 1, 2, 5);
  assert_header("2AAD.dcd", 15, 3, 5);
  assert_header("rna-1to10.dcd", 10632, 10, 5);
}

void test_missing_file_is_file_error() {
  sasmol::DcdReader reader;
  const auto status = reader.open_dcd_read(fixture_path("missing.dcd"));

  assert(status.code == sasmol::IoCode::file_error);
  assert(!reader.is_open());
}

}  // namespace

int main() {
  test_small_fixture_headers();
  test_missing_file_is_file_error();
  return 0;
}
