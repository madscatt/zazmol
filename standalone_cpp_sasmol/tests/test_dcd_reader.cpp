#include "sasmol/calculate.hpp"
#include "sasmol/file_io.hpp"

#include <algorithm>
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

void overwrite_i32(const std::filesystem::path& path, std::streamoff offset,
                   std::int32_t value) {
  std::fstream stream(path, std::ios::binary | std::ios::in | std::ios::out);
  stream.seekp(offset, std::ios::beg);
  stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
  assert(stream.good());
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

void test_fixed_free_atom_header_is_explicitly_unsupported() {
  const auto fixed_free = temp_dcd_path("sasmol_fixed_free_header.dcd");
  auto bytes = dcd_header_prefix(1, 1, 1);
  append_valid_title_and_atom_count(bytes, 5);
  write_bytes(fixed_free, bytes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;

  auto status = reader.open_dcd_read(fixed_free);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.code == sasmol::IoCode::unsupported);

  std::filesystem::remove(fixed_free);
}

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 0.001) {
  assert(std::fabs(static_cast<double>(actual - expected)) <= tolerance);
}

void assert_close_double(double actual, double expected, double tolerance) {
  assert(std::fabs(actual - expected) <= tolerance);
}

sasmol::Vec3 generated_coordinate(std::size_t frame, std::size_t atom) {
  const auto frame_value = static_cast<float>(frame);
  const auto atom_value = static_cast<float>(atom);
  const auto shell = static_cast<float>(atom % 7);
  return {frame_value * 0.125F + atom_value * 0.5F,
          -frame_value * 0.25F + shell * 1.25F,
          frame_value * 0.375F - atom_value * 0.125F};
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

double expected_generated_sum(std::size_t nframes, std::size_t natoms) {
  double total = 0.0;
  for (std::size_t frame = 0; frame < nframes; ++frame) {
    for (std::size_t atom = 0; atom < natoms; ++atom) {
      const auto xyz = generated_coordinate(frame, atom);
      total += static_cast<double>(xyz.x);
      total += static_cast<double>(xyz.y);
      total += static_cast<double>(xyz.z);
    }
  }
  return total;
}

double expected_generated_sum(std::size_t start, std::size_t end,
                              std::size_t natoms) {
  double total = 0.0;
  for (std::size_t frame = start; frame < end; ++frame) {
    for (std::size_t atom = 0; atom < natoms; ++atom) {
      const auto xyz = generated_coordinate(frame, atom);
      total += static_cast<double>(xyz.x);
      total += static_cast<double>(xyz.y);
      total += static_cast<double>(xyz.z);
    }
  }
  return total;
}

sasmol::CoordinateBounds expected_generated_bounds(std::size_t nframes,
                                                   std::size_t natoms) {
  sasmol::CoordinateBounds bounds{
      generated_coordinate(0, 0), generated_coordinate(0, 0)};
  for (std::size_t frame = 0; frame < nframes; ++frame) {
    for (std::size_t atom = 0; atom < natoms; ++atom) {
      const auto xyz = generated_coordinate(frame, atom);
      bounds.minimum.x = std::min(bounds.minimum.x, xyz.x);
      bounds.minimum.y = std::min(bounds.minimum.y, xyz.y);
      bounds.minimum.z = std::min(bounds.minimum.z, xyz.z);
      bounds.maximum.x = std::max(bounds.maximum.x, xyz.x);
      bounds.maximum.y = std::max(bounds.maximum.y, xyz.y);
      bounds.maximum.z = std::max(bounds.maximum.z, xyz.z);
    }
  }
  return bounds;
}

std::filesystem::path generate_streaming_dcd(const char* name,
                                             std::size_t natoms,
                                             std::size_t nframes) {
  const auto output = temp_dcd_path(name);
  sasmol::Molecule frame(natoms, 1);
  sasmol::DcdWriter writer;

  auto status = writer.open_dcd_write(output);
  assert(status.ok());
  status = writer.write_dcd_header(frame, nframes);
  assert(status.ok());
  for (std::size_t frame_index = 0; frame_index < nframes; ++frame_index) {
    for (std::size_t atom = 0; atom < natoms; ++atom) {
      frame.set_coordinate(0, atom, generated_coordinate(frame_index, atom));
    }
    status = writer.write_dcd_step(frame, 0, frame_index + 1);
    assert(status.ok());
  }
  status = writer.close_dcd_write();
  assert(status.ok());

  return output;
}

void assert_bounds_close(const sasmol::CoordinateBounds& actual,
                         const sasmol::CoordinateBounds& expected,
                         double tolerance = 0.001) {
  assert_close(actual.minimum.x, expected.minimum.x, tolerance);
  assert_close(actual.minimum.y, expected.minimum.y, tolerance);
  assert_close(actual.minimum.z, expected.minimum.z, tolerance);
  assert_close(actual.maximum.x, expected.maximum.x, tolerance);
  assert_close(actual.maximum.y, expected.maximum.y, tolerance);
  assert_close(actual.maximum.z, expected.maximum.z, tolerance);
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

void test_streaming_coordinates_auto_reads_header() {
  sasmol::DcdReader reader;
  std::vector<sasmol::Vec3> coordinates;

  auto status = reader.open_dcd_read(fixture_path("1ATM.dcd"));
  assert(status.ok());
  status = reader.read_next_frame_coordinates(coordinates);
  assert(status.ok());

  assert(coordinates.size() == 1);
  assert_close(coordinates[0].x, 76.944F);
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

void test_generated_dcd_streams_without_static_fixture() {
  const std::size_t natoms = 37;
  const std::size_t nframes = 257;
  const auto generated =
      generate_streaming_dcd("sasmol_generated_streaming.dcd", natoms, nframes);

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  std::vector<sasmol::Vec3> coordinates;
  double total = 0.0;
  std::size_t frames_read = 0;

  auto status = reader.open_dcd_read(generated);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.ok());
  assert(header.natoms == natoms);
  assert(header.nframes == nframes);

  while (true) {
    status = reader.read_next_frame_coordinates(coordinates);
    if (status.code == sasmol::IoCode::end_of_file) {
      break;
    }
    assert(status.ok());
    assert(coordinates.size() == natoms);
    if (frames_read == 128) {
      const auto expected = generated_coordinate(frames_read, 17);
      assert_close(coordinates[17].x, expected.x);
      assert_close(coordinates[17].y, expected.y);
      assert_close(coordinates[17].z, expected.z);
    }
    total += coordinate_sum(coordinates);
    ++frames_read;
  }

  assert(frames_read == nframes);
  assert_close_double(total, expected_generated_sum(nframes, natoms), 0.01);
  status = reader.close_dcd_read();
  assert(status.ok());

  std::filesystem::remove(generated);
}

void test_generated_dcd_whole_read_and_minmax() {
  const std::size_t natoms = 37;
  const std::size_t nframes = 257;
  const auto generated =
      generate_streaming_dcd("sasmol_generated_whole_read.dcd", natoms, nframes);

  sasmol::DcdReader reader;
  sasmol::Molecule mol;
  auto status = reader.read_dcd(generated, mol);
  assert(status.ok());
  assert(mol.natoms() == natoms);
  assert(mol.number_of_frames() == nframes);
  assert_close_double(coordinate_sum(mol), expected_generated_sum(nframes, natoms),
                      0.01);

  const auto expected_bounds = expected_generated_bounds(nframes, natoms);
  const auto streamed_bounds =
      sasmol::calculate_minimum_and_maximum_all_steps(generated);
  assert_bounds_close(streamed_bounds, expected_bounds);

  std::filesystem::remove(generated);
}

void test_generated_dcd_single_step_scans_to_late_frame() {
  const std::size_t natoms = 37;
  const std::size_t nframes = 257;
  const std::size_t one_based_frame = 251;
  const auto generated =
      generate_streaming_dcd("sasmol_generated_single_step.dcd", natoms, nframes);

  sasmol::DcdReader reader;
  sasmol::Molecule mol;
  const auto status = reader.read_single_dcd_step(generated, one_based_frame, mol);

  assert(status.ok());
  assert(mol.natoms() == natoms);
  assert(mol.number_of_frames() == 1);
  const auto expected = generated_coordinate(one_based_frame - 1, natoms - 1);
  const auto actual = mol.coordinate(0, natoms - 1);
  assert_close(actual.x, expected.x);
  assert_close(actual.y, expected.y);
  assert_close(actual.z, expected.z);

  std::filesystem::remove(generated);
}

void test_generated_dcd_bad_frame_marker_returns_status() {
  const std::size_t natoms = 37;
  const std::size_t nframes = 3;
  const auto generated =
      generate_streaming_dcd("sasmol_generated_bad_marker.dcd", natoms, nframes);
  const std::streamoff first_x_marker_offset = 276;
  overwrite_i32(generated, first_x_marker_offset,
                static_cast<std::int32_t>(natoms * sizeof(float) + 4));

  sasmol::DcdReader reader;
  sasmol::DcdHeader header;
  std::vector<sasmol::Vec3> coordinates;

  auto status = reader.open_dcd_read(generated);
  assert(status.ok());
  status = reader.read_header(header);
  assert(status.ok());
  status = reader.read_next_frame_coordinates(coordinates);
  assert(status.code == sasmol::IoCode::format_error);

  status = reader.close_dcd_read();
  assert(status.ok());
  std::filesystem::remove(generated);
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

void test_write_dcd_frames_writes_first_and_last_fixture_frames() {
  sasmol::DcdReader reader;
  sasmol::Molecule source;
  auto status = reader.read_dcd(fixture_path("1ATM.dcd"), source);
  assert(status.ok());

  const auto first_output = temp_dcd_path("sasmol_cpp_1atm_first_frame.dcd");
  sasmol::DcdWriter writer;
  status = writer.write_dcd_frames(first_output, source, 0, 1);
  assert(status.ok());

  sasmol::Molecule first;
  status = reader.read_dcd(first_output, first);
  assert(status.ok());
  assert(first.natoms() == 1);
  assert(first.number_of_frames() == 1);
  auto xyz = first.coordinate(0, 0);
  assert_close(xyz.x, 76.944F);
  assert_close(xyz.y, 41.799F);
  assert_close(xyz.z, 41.652F);

  const auto last_output = temp_dcd_path("sasmol_cpp_1atm_last_frame.dcd");
  status = writer.write_dcd_frames(last_output, source, 1, 2);
  assert(status.ok());

  sasmol::Molecule last;
  status = reader.read_dcd(last_output, last);
  assert(status.ok());
  assert(last.natoms() == 1);
  assert(last.number_of_frames() == 1);
  xyz = last.coordinate(0, 0);
  assert_close(xyz.x, 73.944F);
  assert_close(xyz.y, 38.799F);
  assert_close(xyz.z, 41.652F);

  std::filesystem::remove(first_output);
  std::filesystem::remove(last_output);
}

void test_write_dcd_frames_writes_generated_middle_range() {
  const std::size_t natoms = 37;
  const std::size_t nframes = 257;
  const std::size_t start = 120;
  const std::size_t end = 128;
  const auto source_path =
      generate_streaming_dcd("sasmol_generated_range_source.dcd", natoms,
                             nframes);

  sasmol::DcdReader reader;
  sasmol::Molecule source;
  auto status = reader.read_dcd(source_path, source);
  assert(status.ok());

  const auto output = temp_dcd_path("sasmol_generated_range_output.dcd");
  sasmol::DcdWriter writer;
  status = writer.write_dcd_frames(output, source, start, end);
  assert(status.ok());

  sasmol::Molecule range;
  status = reader.read_dcd(output, range);
  assert(status.ok());
  assert(range.natoms() == natoms);
  assert(range.number_of_frames() == end - start);
  assert_close_double(coordinate_sum(range),
                      expected_generated_sum(start, end, natoms), 0.01);
  const auto expected = generated_coordinate(start + 3, 17);
  const auto actual = range.coordinate(3, 17);
  assert_close(actual.x, expected.x);
  assert_close(actual.y, expected.y);
  assert_close(actual.z, expected.z);

  std::filesystem::remove(source_path);
  std::filesystem::remove(output);
}

void test_write_dcd_frames_rejects_invalid_ranges() {
  sasmol::Molecule mol(1, 2);
  sasmol::DcdWriter writer;
  const auto output = temp_dcd_path("sasmol_invalid_range.dcd");

  auto status = writer.write_dcd_frames(output, mol, 1, 1);
  assert(status.code == sasmol::IoCode::format_error);
  assert(!writer.is_open());
  status = writer.write_dcd_frames(output, mol, 0, 3);
  assert(status.code == sasmol::IoCode::format_error);
  assert(!writer.is_open());
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
  test_fixed_free_atom_header_is_explicitly_unsupported();
  test_sequential_frame_reads_1atm();
  test_sequential_frame_reads_2aad_middle_and_final();
  test_sequential_frame_reads_rna_final_sample();
  test_streaming_coordinates_reuses_caller_buffer_without_accumulation();
  test_streaming_coordinates_auto_reads_header();
  test_streaming_coordinates_supports_running_reduction();
  test_generated_dcd_streams_without_static_fixture();
  test_generated_dcd_whole_read_and_minmax();
  test_generated_dcd_single_step_scans_to_late_frame();
  test_generated_dcd_bad_frame_marker_returns_status();
  test_single_step_reopen_scan_uses_one_based_frames();
  test_single_step_rejects_zero_frame_number();
  test_truncated_frame_returns_status();
  test_whole_trajectory_read_1atm();
  test_whole_trajectory_read_2aad();
  test_whole_trajectory_read_rna();
  test_writer_round_trips_1atm();
  test_writer_round_trips_2aad();
  test_writer_convenience_round_trips_rna();
  test_write_dcd_frames_writes_first_and_last_fixture_frames();
  test_write_dcd_frames_writes_generated_middle_range();
  test_write_dcd_frames_rejects_invalid_ranges();
  test_read_dcd_failure_closes_reader();
  test_single_step_past_end_closes_reader();
  test_writer_rejects_unit_cell_option();
  test_writer_rejects_out_of_range_frame();
  test_repeated_open_close_is_safe();
  test_read_header_resets_sequential_position();
  test_read_and_write_after_close_return_not_open();
  return 0;
}
