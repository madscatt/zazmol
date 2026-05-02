#include "sasmol/file_io.hpp"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

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

void write_bytes(const std::filesystem::path& path,
                 const std::vector<char>& bytes) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  assert(out.good());
}

void append_i32(std::vector<char>& bytes, std::int32_t value) {
  const auto* raw = reinterpret_cast<const char*>(&value);
  bytes.insert(bytes.end(), raw, raw + sizeof(value));
}

void append_f64(std::vector<char>& bytes, double value) {
  const auto* raw = reinterpret_cast<const char*>(&value);
  bytes.insert(bytes.end(), raw, raw + sizeof(value));
}

std::vector<char> dcd_header_prefix(std::int32_t nframes,
                                    std::int32_t nsavc = 1,
                                    std::int32_t namnf = 0) {
  std::vector<char> bytes;
  append_i32(bytes, 84);
  bytes.insert(bytes.end(), {'C', 'O', 'R', 'D'});
  append_i32(bytes, nframes);
  append_i32(bytes, 0);
  append_i32(bytes, nsavc);
  for (int i = 0; i < 5; ++i) {
    append_i32(bytes, 0);
  }
  append_i32(bytes, namnf);
  append_f64(bytes, 1.0);
  for (int i = 0; i < 9; ++i) {
    append_i32(bytes, 0);
  }
  append_i32(bytes, 84);
  assert(bytes.size() == 92);
  return bytes;
}

void append_valid_title_and_atom_count(std::vector<char>& bytes,
                                       std::int32_t natoms) {
  append_i32(bytes, 84);
  append_i32(bytes, 1);
  bytes.insert(bytes.end(), 80, ' ');
  append_i32(bytes, 84);
  append_i32(bytes, 4);
  append_i32(bytes, natoms);
  append_i32(bytes, 4);
}

std::filesystem::path temp_dcd_path(const char* name) {
  const auto tag = std::filesystem::path(SASMOL_TEST_BUILD_DIR).filename();
  return std::filesystem::temp_directory_path() /
         ("sasmol_" + tag.string() + "_" + name);
}

void test_truncated_header_returns_status() {
  const auto truncated = temp_dcd_path("sasmol_truncated_header.dcd");
  write_bytes(truncated, {'n', 'o', 't', 'd', 'c', 'd'});

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(truncated);
  assert(status.ok());
  status = reader.read_header(header);
  assert(!status.ok());
  assert(status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(truncated);
}

void test_bad_cord_magic_returns_status() {
  const auto bad_magic = temp_dcd_path("sasmol_bad_cord.dcd");
  auto bytes = std::vector<char>(4 + 84 + 4, '\0');
  const std::int32_t marker = 84;
  std::memcpy(bytes.data(), &marker, sizeof(marker));
  bytes[4] = 'B';
  bytes[5] = 'A';
  bytes[6] = 'D';
  bytes[7] = '!';
  std::memcpy(bytes.data() + 4 + 84, &marker, sizeof(marker));
  write_bytes(bad_magic, bytes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(bad_magic);
  assert(status.ok());
  status = reader.read_header(header);
  assert(!status.ok());
  assert(status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(bad_magic);
}

void test_invalid_title_block_size_returns_status() {
  const auto bad_title = temp_dcd_path("sasmol_bad_title_size.dcd");
  auto bytes = dcd_header_prefix(1);
  append_i32(bytes, 5);
  write_bytes(bad_title, bytes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(bad_title);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(bad_title);
}

void test_invalid_title_count_returns_status() {
  const auto bad_title = temp_dcd_path("sasmol_bad_title_count.dcd");
  auto bytes = dcd_header_prefix(1);
  append_i32(bytes, 84);
  append_i32(bytes, 2);
  write_bytes(bad_title, bytes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(bad_title);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(bad_title);
}

void test_negative_header_count_returns_status() {
  const auto bad_count = temp_dcd_path("sasmol_negative_header_count.dcd");
  const auto bytes = dcd_header_prefix(-1);
  write_bytes(bad_count, bytes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(bad_count);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(bad_count);
}

void test_oversized_atom_count_returns_unsupported_before_allocation() {
  const auto oversized = temp_dcd_path("sasmol_oversized_atom_count.dcd");
  auto bytes = dcd_header_prefix(1);
  append_valid_title_and_atom_count(bytes, std::numeric_limits<std::int32_t>::max());
  write_bytes(oversized, bytes);

  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(oversized);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.code == sasmol::IoCode::unsupported);

  std::filesystem::remove(oversized);
}

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 0.001) {
  assert(std::fabs(static_cast<double>(actual - expected)) <= tolerance);
}

void assert_close_double(double actual, double expected, double tolerance) {
  assert(std::fabs(actual - expected) <= tolerance);
}

double coordinate_sum(const sasmol::Molecule& mol) {
  double total = 0.0;
  for (const auto value : mol.coor()) {
    total += static_cast<double>(value);
  }
  return total;
}

double coordinate_sum(const std::vector<sasmol::Vec3>& coordinates) {
  double total = 0.0;
  for (const auto& xyz : coordinates) {
    total += static_cast<double>(xyz.x);
    total += static_cast<double>(xyz.y);
    total += static_cast<double>(xyz.z);
  }
  return total;
}

void test_sequential_frame_reads_1atm() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());

  status = reader.read_next_frame(mol);
  assert(status.ok());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 2);
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);

  status = reader.read_next_frame(mol);
  assert(status.ok());
  xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, 73.944F);
  assert_close(xyz.y, 38.799F);
  assert_close(xyz.z, 41.652F);

  status = reader.read_next_frame(mol);
  assert(status.code == sasmol::IoCode::end_of_file);

  assert_close_double(coordinate_sum(mol), 314.790, 0.001);
}

void test_sequential_frame_reads_2aad_middle_and_final() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("2AAD.dcd"));
  assert(status.ok());

  status = reader.read_next_frame(mol);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());

  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 3);

  auto xyz = mol.coordinate(1, 0);
  assert_close(xyz.x, -73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);

  xyz = mol.coordinate(2, 14);
  assert_close(xyz.x, 76.970F);
  assert_close(xyz.y, -46.273F);
  assert_close(xyz.z, 42.000F);

  assert_close_double(coordinate_sum(mol), 3644.294, 0.01);
}

void test_sequential_frame_reads_rna_final_sample() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("rna-1to10.dcd"));
  assert(status.ok());

  for (int i = 0; i < 10; ++i) {
    status = reader.read_next_frame(mol);
    assert(status.ok());
  }

  assert(mol.natoms() == 10632);
  assert(mol.number_of_frames() == 10);

  const auto xyz = mol.coordinate(9, 10631);
  assert_close(xyz.x, -6.392F);
  assert_close(xyz.y, 14.348F);
  assert_close(xyz.z, 20.914F);
  assert_close_double(coordinate_sum(mol), -430804.378, 0.1);
}

void test_streaming_coordinates_reuses_caller_buffer_without_accumulation() {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  std::vector<sasmol::Vec3> coordinates;
  coordinates.reserve(100);
  const auto reserved_capacity = coordinates.capacity();

  auto status = reader.open_dcd_read(fixture_path("2AAD.dcd"));
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.ok());
  assert(header.natoms == 15);
  assert(header.nframes == 3);

  status = reader.read_next_frame_coordinates(coordinates);
  assert(status.ok());
  assert(coordinates.size() == 15);
  assert(coordinates.capacity() == reserved_capacity);
  const auto* first_buffer = coordinates.data();
  assert_close(coordinates[0].x, 73.944F);
  assert_close(coordinates[0].y, 41.799F);
  assert_close(coordinates[0].z, 41.652F);

  status = reader.read_next_frame_coordinates(coordinates);
  assert(status.ok());
  assert(coordinates.size() == 15);
  assert(coordinates.capacity() == reserved_capacity);
  assert(coordinates.data() == first_buffer);
  assert_close(coordinates[0].x, -73.944F);
  assert_close(coordinates[0].y, 41.799F);
  assert_close(coordinates[0].z, 41.652F);

  status = reader.close_dcd_read();
  assert(status.ok());
}

void test_streaming_coordinates_supports_running_reduction() {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  std::vector<sasmol::Vec3> coordinates;
  double total = 0.0;
  std::size_t frames_read = 0;

  auto status = reader.open_dcd_read(fixture_path("2AAD.dcd"));
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.ok());

  while (true) {
    status = reader.read_next_frame_coordinates(coordinates);
    if (status.code == sasmol::IoCode::end_of_file) {
      break;
    }
    assert(status.ok());
    assert(coordinates.size() == header.natoms);
    total += coordinate_sum(coordinates);
    ++frames_read;
  }

  assert(frames_read == header.nframes);
  assert_close_double(total, 3644.294, 0.01);

  status = reader.close_dcd_read();
  assert(status.ok());
}

void test_single_step_reopen_scan_uses_one_based_frames() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.read_single_dcd_step(fixture_path("2AAD.dcd"), 2, mol);
  assert(status.ok());
  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 1);

  const auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, -73.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);
  assert_close_double(coordinate_sum(mol), 140.846, 0.01);
}

void test_single_step_rejects_zero_frame_number() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  const auto status =
      reader.read_single_dcd_step(fixture_path("1ATM.dcd"), 0, mol);

  assert(status.code == sasmol::IoCode::format_error);
}

void test_truncated_frame_returns_status() {
  const auto source = fixture_path("1ATM.dcd");
  const auto truncated = temp_dcd_path("sasmol_truncated_frame.dcd");
  std::filesystem::copy_file(source, truncated,
                             std::filesystem::copy_options::overwrite_existing);
  const auto original_size = std::filesystem::file_size(truncated);
  std::filesystem::resize_file(truncated, original_size - 8);

  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(truncated);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(!status.ok());
  assert(status.code == sasmol::IoCode::end_of_file ||
         status.code == sasmol::IoCode::format_error);

  std::filesystem::remove(truncated);
}

void test_whole_trajectory_read_1atm() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  const auto status = reader.read_dcd(fixture_path("1ATM.dcd"), mol);

  assert(status.ok());
  assert(!reader.is_open());
  assert(mol.natoms() == 1);
  assert(mol.number_of_frames() == 2);
  assert_close_double(coordinate_sum(mol), 314.790, 0.001);
}

void test_whole_trajectory_read_2aad() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  const auto status = reader.read_dcd(fixture_path("2AAD.dcd"), mol);

  assert(status.ok());
  assert(!reader.is_open());
  assert(mol.natoms() == 15);
  assert(mol.number_of_frames() == 3);
  const auto xyz = mol.coordinate(2, 14);
  assert_close(xyz.x, 76.970F);
  assert_close(xyz.y, -46.273F);
  assert_close(xyz.z, 42.000F);
  assert_close_double(coordinate_sum(mol), 3644.294, 0.01);
}

void test_whole_trajectory_read_rna() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  const auto status = reader.read_dcd(fixture_path("rna-1to10.dcd"), mol);

  assert(status.ok());
  assert(!reader.is_open());
  assert(mol.natoms() == 10632);
  assert(mol.number_of_frames() == 10);
  const auto xyz = mol.coordinate(9, 10631);
  assert_close(xyz.x, -6.392F);
  assert_close(xyz.y, 14.348F);
  assert_close(xyz.z, 20.914F);
  assert_close_double(coordinate_sum(mol), -430804.378, 0.1);
}

void write_round_trip_fixture(const char* source_name, const char* temp_name,
                              double expected_sum, double sum_tolerance) {
  sasmol::DcdReader reader;
  sasmol::Molecule source;
  auto status = reader.read_dcd(fixture_path(source_name), source);
  assert(status.ok());

  const auto output = temp_dcd_path(temp_name);
  sasmol::DcdWriter writer;
  status = writer.open_dcd_write(output);
  assert(status.ok());
  status = writer.write_dcd_header(source, source.number_of_frames());
  assert(status.ok());
  for (std::size_t frame = 0; frame < source.number_of_frames(); ++frame) {
    status = writer.write_dcd_step(source, frame, frame + 1);
    assert(status.ok());
  }
  status = writer.close_dcd_write();
  assert(status.ok());

  sasmol::Molecule round_trip;
  status = reader.read_dcd(output, round_trip);
  assert(status.ok());
  assert(round_trip.natoms() == source.natoms());
  assert(round_trip.number_of_frames() == source.number_of_frames());
  assert_close_double(coordinate_sum(round_trip), expected_sum, sum_tolerance);

  std::filesystem::remove(output);
}

void test_writer_round_trips_1atm() {
  write_round_trip_fixture("1ATM.dcd", "sasmol_cpp_1atm_roundtrip.dcd",
                           314.790, 0.001);
}

void test_writer_round_trips_2aad() {
  write_round_trip_fixture("2AAD.dcd", "sasmol_cpp_2aad_roundtrip.dcd",
                           3644.294, 0.01);
}

void test_writer_convenience_round_trips_rna() {
  sasmol::DcdReader reader;
  sasmol::Molecule source;
  auto status = reader.read_dcd(fixture_path("rna-1to10.dcd"), source);
  assert(status.ok());

  const auto output = temp_dcd_path("sasmol_cpp_rna_roundtrip.dcd");
  sasmol::DcdWriter writer;
  status = writer.write_dcd(output, source);
  assert(status.ok());
  assert(!writer.is_open());

  sasmol::Molecule round_trip;
  status = reader.read_dcd(output, round_trip);
  assert(status.ok());
  assert(round_trip.natoms() == 10632);
  assert(round_trip.number_of_frames() == 10);
  const auto xyz = round_trip.coordinate(9, 10631);
  assert_close(xyz.x, -6.392F);
  assert_close(xyz.y, 14.348F);
  assert_close(xyz.z, 20.914F);
  assert_close_double(coordinate_sum(round_trip), -430804.378, 0.1);

  std::filesystem::remove(output);
}

void test_read_dcd_failure_closes_reader() {
  const auto source = fixture_path("1ATM.dcd");
  const auto truncated = temp_dcd_path("sasmol_truncated_read_dcd.dcd");
  std::filesystem::copy_file(source, truncated,
                             std::filesystem::copy_options::overwrite_existing);
  const auto original_size = std::filesystem::file_size(truncated);
  std::filesystem::resize_file(truncated, original_size - 8);

  sasmol::DcdReader reader;
  sasmol::Molecule mol;
  const auto status = reader.read_dcd(truncated, mol);

  assert(!status.ok());
  assert(!reader.is_open());

  std::filesystem::remove(truncated);
}

void test_single_step_past_end_closes_reader() {
  sasmol::DcdReader reader;
  sasmol::Molecule mol;

  const auto status =
      reader.read_single_dcd_step(fixture_path("1ATM.dcd"), 3, mol);

  assert(status.code == sasmol::IoCode::end_of_file);
  assert(!reader.is_open());
}

void test_writer_rejects_unit_cell_option() {
  sasmol::DcdWriter writer;
  sasmol::DcdWriteOptions options;
  options.include_unit_cell = true;

  const auto status =
      writer.open_dcd_write(temp_dcd_path("sasmol_unit_cell_write.dcd"), options);

  assert(status.code == sasmol::IoCode::unsupported);
  assert(!writer.is_open());
}

void test_writer_rejects_out_of_range_frame() {
  sasmol::DcdWriter writer;
  sasmol::Molecule mol(1, 1);
  const auto output = temp_dcd_path("sasmol_out_of_range_write.dcd");

  auto status = writer.open_dcd_write(output);
  assert(status.ok());
  status = writer.write_dcd_header(mol, 1);
  assert(status.ok());
  status = writer.write_dcd_step(mol, 1, 1);
  assert(status.code == sasmol::IoCode::format_error);
  status = writer.close_dcd_write();
  assert(status.ok());

  std::filesystem::remove(output);
}

void test_repeated_open_close_is_safe() {
  sasmol::DcdReader reader;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());
  status = reader.open_dcd_read(fixture_path("2AAD.dcd"));
  assert(status.ok());

  sasmol::DcdHeader header;
  status = reader.read_header(header);
  assert(status.ok());
  assert(header.natoms == 15);

  status = reader.close_dcd_read();
  assert(status.ok());
  status = reader.close_dcd_read();
  assert(status.ok());
}

void test_read_header_resets_sequential_position() {
  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  sasmol::Molecule mol;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());
  auto xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);

  status = reader.read_header(header);
  assert(status.ok());
  status = reader.read_next_frame(mol);
  assert(status.ok());
  xyz = mol.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);
}

void test_read_and_write_after_close_return_not_open() {
  sasmol::DcdReader reader;
  sasmol::DcdWriter writer;
  sasmol::Molecule mol(1, 1);
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());
  status = reader.close_dcd_read();
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::not_open);
  status = reader.read_next_frame(mol);
  assert(status.code == sasmol::IoCode::not_open);

  const auto output = temp_dcd_path("sasmol_after_close_write.dcd");
  status = writer.open_dcd_write(output);
  assert(status.ok());
  status = writer.close_dcd_write();
  assert(status.ok());
  status = writer.write_dcd_header(mol, 1);
  assert(status.code == sasmol::IoCode::not_open);
  status = writer.write_dcd_step(mol, 0, 1);
  assert(status.code == sasmol::IoCode::not_open);

  std::filesystem::remove(output);
}

}  // namespace

int main() {
  test_small_fixture_headers();
  test_missing_file_is_file_error();
  test_truncated_header_returns_status();
  test_bad_cord_magic_returns_status();
  test_invalid_title_block_size_returns_status();
  test_invalid_title_count_returns_status();
  test_negative_header_count_returns_status();
  test_oversized_atom_count_returns_unsupported_before_allocation();
  test_sequential_frame_reads_1atm();
  test_sequential_frame_reads_2aad_middle_and_final();
  test_sequential_frame_reads_rna_final_sample();
  test_streaming_coordinates_reuses_caller_buffer_without_accumulation();
  test_streaming_coordinates_supports_running_reduction();
  test_single_step_reopen_scan_uses_one_based_frames();
  test_single_step_rejects_zero_frame_number();
  test_truncated_frame_returns_status();
  test_whole_trajectory_read_1atm();
  test_whole_trajectory_read_2aad();
  test_whole_trajectory_read_rna();
  test_writer_round_trips_1atm();
  test_writer_round_trips_2aad();
  test_writer_convenience_round_trips_rna();
  test_read_dcd_failure_closes_reader();
  test_single_step_past_end_closes_reader();
  test_writer_rejects_unit_cell_option();
  test_writer_rejects_out_of_range_frame();
  test_repeated_open_close_is_safe();
  test_read_header_resets_sequential_position();
  test_read_and_write_after_close_return_not_open();
  return 0;
}
