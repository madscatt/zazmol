#include "sasmol/file_io.hpp"

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <ios>
#include <vector>
#include <utility>

namespace sasmol {

namespace {

constexpr int dcd_is_charmm = 0x01;
constexpr int dcd_has_4dims = 0x02;
constexpr int dcd_has_extra_block = 0x04;

std::uint32_t byte_swap32(std::uint32_t value) {
  return ((value & 0x000000FFU) << 24U) | ((value & 0x0000FF00U) << 8U) |
         ((value & 0x00FF0000U) >> 8U) | ((value & 0xFF000000U) >> 24U);
}

std::uint64_t byte_swap64(std::uint64_t value) {
  return ((value & 0x00000000000000FFULL) << 56U) |
         ((value & 0x000000000000FF00ULL) << 40U) |
         ((value & 0x0000000000FF0000ULL) << 24U) |
         ((value & 0x00000000FF000000ULL) << 8U) |
         ((value & 0x000000FF00000000ULL) >> 8U) |
         ((value & 0x0000FF0000000000ULL) >> 24U) |
         ((value & 0x00FF000000000000ULL) >> 40U) |
         ((value & 0xFF00000000000000ULL) >> 56U);
}

template <typename T>
bool read_exact(std::istream& stream, T* value) {
  stream.read(reinterpret_cast<char*>(value), sizeof(T));
  return static_cast<bool>(stream);
}

bool read_bytes(std::istream& stream, char* data, std::streamsize count) {
  stream.read(data, count);
  return static_cast<bool>(stream);
}

int read_int32_from_bytes(const char* data, bool reverse_endian) {
  std::uint32_t raw{};
  std::memcpy(&raw, data, sizeof(raw));
  if (reverse_endian) {
    raw = byte_swap32(raw);
  }
  std::int32_t value{};
  std::memcpy(&value, &raw, sizeof(value));
  return value;
}

float read_float_from_bytes(const char* data, bool reverse_endian) {
  std::uint32_t raw{};
  std::memcpy(&raw, data, sizeof(raw));
  if (reverse_endian) {
    raw = byte_swap32(raw);
  }
  float value{};
  std::memcpy(&value, &raw, sizeof(value));
  return value;
}

double read_double_from_bytes(const char* data, bool reverse_endian) {
  std::uint64_t raw{};
  std::memcpy(&raw, data, sizeof(raw));
  if (reverse_endian) {
    raw = byte_swap64(raw);
  }
  double value{};
  std::memcpy(&value, &raw, sizeof(value));
  return value;
}

IoStatus read_record_marker(std::istream& stream, int expected,
                            bool reverse_endian,
                            const std::string& context) {
  std::uint32_t raw{};
  if (!read_exact(stream, &raw)) {
    return {IoCode::end_of_file, "Unexpected EOF while " + context + "."};
  }
  if (reverse_endian) {
    raw = byte_swap32(raw);
  }
  std::int32_t marker{};
  std::memcpy(&marker, &raw, sizeof(marker));
  if (marker != expected) {
    return {IoCode::format_error, "Unexpected DCD record marker while " +
                                      context + "."};
  }
  return IoStatus::success();
}

IoStatus skip_bytes(std::istream& stream, std::streamoff count,
                    const std::string& context) {
  stream.seekg(count, std::ios::cur);
  if (!stream) {
    return {IoCode::end_of_file, "Unexpected EOF while " + context + "."};
  }
  return IoStatus::success();
}

IoStatus read_float_block(std::istream& stream, std::vector<float>& values,
                          bool reverse_endian, const std::string& axis) {
  const int expected_size = static_cast<int>(values.size() * sizeof(float));
  auto status = read_record_marker(stream, expected_size, reverse_endian,
                                   "reading " + axis + " block marker");
  if (!status) {
    return status;
  }

  if (!read_bytes(stream, reinterpret_cast<char*>(values.data()),
                  static_cast<std::streamsize>(expected_size))) {
    return {IoCode::end_of_file, "Unexpected EOF while reading DCD " + axis +
                                " coordinate block."};
  }

  if (reverse_endian) {
    for (float& value : values) {
      std::uint32_t raw{};
      std::memcpy(&raw, &value, sizeof(raw));
      raw = byte_swap32(raw);
      std::memcpy(&value, &raw, sizeof(value));
    }
  }

  return read_record_marker(stream, expected_size, reverse_endian,
                            "reading trailing " + axis + " block marker");
}

}  // namespace

IoStatus IoStatus::success() { return {}; }

IoStatus IoStatus::not_implemented(std::string message) {
  return {IoCode::not_implemented, std::move(message)};
}

IoStatus PdbReader::read_pdb(const std::filesystem::path& filename,
                             Molecule& molecule,
                             const PdbReadOptions& options) const {
  (void)filename;
  (void)molecule;
  (void)options;
  return IoStatus::not_implemented(
      "PDB parsing is intentionally deferred until the tolerance contract is "
      "fully reviewed against Python zazmol fixtures.");
}

IoStatus PdbWriter::write_pdb(const std::filesystem::path& filename,
                              const Molecule& molecule,
                              const PdbWriteOptions& options) const {
  (void)filename;
  (void)molecule;
  (void)options;
  return IoStatus::not_implemented(
      "PDB writing is intentionally deferred until descriptor and formatting "
      "parity are reviewed against Python zazmol fixtures.");
}

IoStatus DcdReader::open_dcd_read(const std::filesystem::path& filename,
                                  const DcdReadOptions& options) {
  (void)close_dcd_read();
  filename_ = filename;
  options_ = options;
  header_ = {};
  header_read_ = false;
  current_frame_ = 0;

  stream_.open(filename_, std::ios::binary);
  if (!stream_) {
    filename_.clear();
    return {IoCode::file_error, "Failed to open DCD file: " +
                                    filename.string()};
  }

  open_ = true;
  return IoStatus::success();
}

IoStatus DcdReader::read_header(DcdHeader& header) {
  if (!open_) {
    return {IoCode::not_open, "DCD reader is not open."};
  }

  stream_.clear();
  stream_.seekg(0, std::ios::beg);
  if (!stream_) {
    return {IoCode::file_error, "Failed to seek to DCD header."};
  }

  std::uint32_t raw_marker{};
  if (!read_exact(stream_, &raw_marker)) {
    return {IoCode::format_error, "Failed to read initial DCD record marker."};
  }

  std::int32_t marker{};
  std::memcpy(&marker, &raw_marker, sizeof(marker));
  bool reverse_endian = false;
  if (marker != 84) {
    raw_marker = byte_swap32(raw_marker);
    std::memcpy(&marker, &raw_marker, sizeof(marker));
    if (marker != 84) {
      return {IoCode::format_error, "Invalid initial DCD record marker."};
    }
    reverse_endian = true;
  }

  std::array<char, 84> header_block{};
  if (!read_bytes(stream_, header_block.data(), header_block.size())) {
    return {IoCode::format_error, "Failed to read DCD 84-byte header block."};
  }

  if (header_block[0] != 'C' || header_block[1] != 'O' ||
      header_block[2] != 'R' || header_block[3] != 'D') {
    return {IoCode::format_error, "DCD header does not contain CORD magic."};
  }

  auto status = read_record_marker(stream_, 84, reverse_endian,
                                   "reading trailing header marker");
  if (!status) {
    return status;
  }

  const int charmm_version =
      read_int32_from_bytes(header_block.data() + 80, reverse_endian);
  int charmm_flags = 0;
  if (charmm_version != 0) {
    charmm_flags = dcd_is_charmm;
    if (read_int32_from_bytes(header_block.data() + 44, reverse_endian) == 1) {
      charmm_flags |= dcd_has_extra_block;
    }
    if (read_int32_from_bytes(header_block.data() + 48, reverse_endian) == 1) {
      charmm_flags |= dcd_has_4dims;
    }
  }

  DcdHeader parsed;
  parsed.nframes = static_cast<std::size_t>(
      read_int32_from_bytes(header_block.data() + 4, reverse_endian));
  parsed.istart = read_int32_from_bytes(header_block.data() + 8, reverse_endian);
  parsed.nsavc = read_int32_from_bytes(header_block.data() + 12, reverse_endian);
  parsed.namnf = read_int32_from_bytes(header_block.data() + 36, reverse_endian);
  parsed.delta = (charmm_flags & dcd_is_charmm)
                     ? static_cast<double>(
                           read_float_from_bytes(header_block.data() + 40,
                                                 reverse_endian))
                     : read_double_from_bytes(header_block.data() + 40,
                                              reverse_endian);
  parsed.reverse_endian = reverse_endian;
  parsed.has_unit_cell = (charmm_flags & dcd_has_extra_block) != 0;
  parsed.charmm_format = (charmm_flags & dcd_is_charmm) != 0;
  parsed.charmm_flags = charmm_flags;

  std::uint32_t raw_title_size{};
  if (!read_exact(stream_, &raw_title_size)) {
    return {IoCode::format_error, "Failed to read DCD title block size."};
  }
  if (reverse_endian) {
    raw_title_size = byte_swap32(raw_title_size);
  }
  std::int32_t title_size{};
  std::memcpy(&title_size, &raw_title_size, sizeof(title_size));
  if (title_size < 4 || ((title_size - 4) % 80) != 0) {
    return {IoCode::format_error, "Invalid DCD title block size."};
  }

  std::uint32_t raw_ntitle{};
  if (!read_exact(stream_, &raw_ntitle)) {
    return {IoCode::format_error, "Failed to read DCD title count."};
  }
  if (reverse_endian) {
    raw_ntitle = byte_swap32(raw_ntitle);
  }
  std::int32_t ntitle{};
  std::memcpy(&ntitle, &raw_ntitle, sizeof(ntitle));
  if (ntitle < 0 || 4 + (ntitle * 80) != title_size) {
    return {IoCode::format_error, "Invalid DCD title count."};
  }

  status = skip_bytes(stream_, static_cast<std::streamoff>(ntitle) * 80,
                      "skipping DCD title strings");
  if (!status) {
    return status;
  }

  status =
      read_record_marker(stream_, title_size, reverse_endian,
                         "reading trailing DCD title block marker");
  if (!status) {
    return status;
  }

  status = read_record_marker(stream_, 4, reverse_endian,
                              "reading DCD atom-count block marker");
  if (!status) {
    return status;
  }

  std::uint32_t raw_natoms{};
  if (!read_exact(stream_, &raw_natoms)) {
    return {IoCode::format_error, "Failed to read DCD atom count."};
  }
  if (reverse_endian) {
    raw_natoms = byte_swap32(raw_natoms);
  }
  std::int32_t natoms{};
  std::memcpy(&natoms, &raw_natoms, sizeof(natoms));
  if (natoms < 0) {
    return {IoCode::format_error, "Invalid negative DCD atom count."};
  }
  parsed.natoms = static_cast<std::size_t>(natoms);

  status = read_record_marker(stream_, 4, reverse_endian,
                              "reading trailing DCD atom-count marker");
  if (!status) {
    return status;
  }

  if (parsed.namnf != 0) {
    return {IoCode::unsupported,
            "DCD fixed/free atom index tables are detected but not yet "
            "implemented."};
  }

  header_ = parsed;
  header_read_ = true;
  current_frame_ = 0;
  header = parsed;
  return IoStatus::success();
}

IoStatus DcdReader::read_next_frame(Molecule& molecule) {
  if (!open_) {
    return {IoCode::not_open, "DCD reader is not open."};
  }
  if (!header_read_) {
    DcdHeader header;
    auto status = read_header(header);
    if (!status) {
      return status;
    }
  }
  if (current_frame_ >= header_.nframes) {
    return {IoCode::end_of_file, "No more DCD frames are available."};
  }
  if (header_.namnf != 0) {
    return {IoCode::unsupported,
            "DCD fixed/free atom frame reading is not yet implemented."};
  }

  const bool reverse_endian = header_.reverse_endian;
  if (header_.has_unit_cell) {
    std::uint32_t raw_size{};
    if (!read_exact(stream_, &raw_size)) {
      return {IoCode::end_of_file,
              "Unexpected EOF while reading DCD unit-cell block marker."};
    }
    if (reverse_endian) {
      raw_size = byte_swap32(raw_size);
    }
    std::int32_t block_size{};
    std::memcpy(&block_size, &raw_size, sizeof(block_size));
    if (block_size < 0) {
      return {IoCode::format_error, "Invalid negative DCD unit-cell block."};
    }
    auto status = skip_bytes(stream_, block_size, "skipping DCD unit-cell block");
    if (!status) {
      return status;
    }
    status = read_record_marker(stream_, block_size, reverse_endian,
                                "reading trailing DCD unit-cell block marker");
    if (!status) {
      return status;
    }
  }

  std::vector<float> x(header_.natoms);
  std::vector<float> y(header_.natoms);
  std::vector<float> z(header_.natoms);

  auto status = read_float_block(stream_, x, reverse_endian, "X");
  if (!status) {
    return status;
  }
  status = read_float_block(stream_, y, reverse_endian, "Y");
  if (!status) {
    return status;
  }
  status = read_float_block(stream_, z, reverse_endian, "Z");
  if (!status) {
    return status;
  }

  if (header_.charmm_flags & dcd_has_4dims) {
    std::uint32_t raw_size{};
    if (!read_exact(stream_, &raw_size)) {
      return {IoCode::end_of_file,
              "Unexpected EOF while reading DCD 4D block marker."};
    }
    if (reverse_endian) {
      raw_size = byte_swap32(raw_size);
    }
    std::int32_t block_size{};
    std::memcpy(&block_size, &raw_size, sizeof(block_size));
    if (block_size < 0) {
      return {IoCode::format_error, "Invalid negative DCD 4D block."};
    }
    status = skip_bytes(stream_, block_size, "skipping DCD 4D block");
    if (!status) {
      return status;
    }
    status = read_record_marker(stream_, block_size, reverse_endian,
                                "reading trailing DCD 4D block marker");
    if (!status) {
      return status;
    }
  }

  if (molecule.natoms() != header_.natoms ||
      molecule.number_of_frames() != header_.nframes) {
    molecule.resize(header_.natoms, header_.nframes);
  }

  for (std::size_t atom = 0; atom < header_.natoms; ++atom) {
    molecule.set_coordinate(current_frame_, atom, {x[atom], y[atom], z[atom]});
  }

  ++current_frame_;
  return IoStatus::success();
}

IoStatus DcdReader::close_dcd_read() {
  if (stream_.is_open()) {
    stream_.close();
  }
  open_ = false;
  header_read_ = false;
  header_ = {};
  current_frame_ = 0;
  filename_.clear();
  return IoStatus::success();
}

IoStatus DcdReader::read_single_dcd_step(const std::filesystem::path& filename,
                                         std::size_t frame, Molecule& molecule,
                                         const DcdReadOptions& options) {
  if (frame == 0) {
    return {IoCode::format_error,
            "DCD single-step reads use one-based frame numbers."};
  }

  DcdReader reader;
  auto status = reader.open_dcd_read(filename, options);
  if (!status) {
    return status;
  }

  Molecule trajectory;
  for (std::size_t step = 0; step < frame; ++step) {
    status = reader.read_next_frame(trajectory);
    if (!status) {
      (void)reader.close_dcd_read();
      return status;
    }
  }

  molecule.resize(trajectory.natoms(), 1);
  for (std::size_t atom = 0; atom < trajectory.natoms(); ++atom) {
    molecule.set_coordinate(0, atom, trajectory.coordinate(frame - 1, atom));
  }

  return reader.close_dcd_read();
}

IoStatus DcdReader::read_dcd(const std::filesystem::path& filename,
                             Molecule& molecule,
                             const DcdReadOptions& options) {
  auto status = open_dcd_read(filename, options);
  if (!status) {
    return status;
  }

  DcdHeader header;
  status = read_header(header);
  if (!status) {
    (void)close_dcd_read();
    return status;
  }

  molecule.resize(header.natoms, header.nframes);
  for (std::size_t frame = 0; frame < header.nframes; ++frame) {
    status = read_next_frame(molecule);
    if (!status) {
      (void)close_dcd_read();
      return status;
    }
  }

  return close_dcd_read();
}

IoStatus DcdWriter::open_dcd_write(const std::filesystem::path& filename,
                                   const DcdWriteOptions& options) {
  filename_ = filename;
  options_ = options;
  open_ = true;
  return IoStatus::not_implemented(
      "DCD writing is intentionally deferred until Python/C-extension write "
      "parity is captured in tests.");
}

IoStatus DcdWriter::write_dcd_header(const Molecule& molecule,
                                     std::size_t nframes) {
  (void)molecule;
  (void)nframes;
  if (!open_) {
    return {IoCode::not_open, "DCD writer is not open."};
  }
  return IoStatus::not_implemented(
      "DCD header writing is intentionally deferred until legacy header fields "
      "are documented.");
}

IoStatus DcdWriter::write_dcd_step(const Molecule& molecule, std::size_t frame,
                                   std::size_t step) {
  (void)molecule;
  (void)frame;
  (void)step;
  if (!open_) {
    return {IoCode::not_open, "DCD writer is not open."};
  }
  return IoStatus::not_implemented(
      "DCD frame writing is intentionally deferred until coordinate layout and "
      "round-trip parity are locked.");
}

IoStatus DcdWriter::close_dcd_write() {
  open_ = false;
  filename_.clear();
  return IoStatus::success();
}

}  // namespace sasmol
