#include "sasmol/file_io.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <ios>
#include <limits>
#include <map>
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
  if (values.size() >
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) /
          sizeof(float)) {
    return {IoCode::unsupported,
            "DCD coordinate block is too large for int32 record markers."};
  }
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

template <typename T>
IoStatus write_exact(std::ostream& stream, const T& value,
                     const std::string& context) {
  stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
  if (!stream) {
    return {IoCode::file_error, "Failed while " + context + "."};
  }
  return IoStatus::success();
}

IoStatus write_bytes(std::ostream& stream, const char* data,
                     std::streamsize count, const std::string& context) {
  stream.write(data, count);
  if (!stream) {
    return {IoCode::file_error, "Failed while " + context + "."};
  }
  return IoStatus::success();
}

std::array<char, 80> padded_title(const std::string& text) {
  std::array<char, 80> title{};
  title.fill(' ');
  const std::size_t count = std::min(text.size(), title.size());
  std::memcpy(title.data(), text.data(), count);
  return title;
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

IoStatus PdbWriter::check_for_all_zero_columns(Molecule& molecule,
                                               std::size_t frame) const {
  if (frame >= molecule.number_of_frames()) {
    return {IoCode::format_error, "PDB coordinate frame is out of range."};
  }
  if (molecule.natoms() == 0) {
    return {IoCode::format_error, "PDB coordinate guard requires atoms."};
  }

  constexpr calc_type small = 1.0e-10;
  calc_type x_sum{};
  calc_type y_sum{};
  calc_type z_sum{};

  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto xyz = molecule.coordinate(frame, atom);
    x_sum += static_cast<calc_type>(xyz.x) * static_cast<calc_type>(xyz.x);
    y_sum += static_cast<calc_type>(xyz.y) * static_cast<calc_type>(xyz.y);
    z_sum += static_cast<calc_type>(xyz.z) * static_cast<calc_type>(xyz.z);
  }

  auto first = molecule.coordinate(frame, 0);
  if (x_sum < small) {
    first.x = static_cast<coord_type>(small);
  }
  if (y_sum < small) {
    first.y = static_cast<coord_type>(small);
  }
  if (z_sum < small) {
    first.z = static_cast<coord_type>(small);
  }
  molecule.set_coordinate(frame, 0, first);
  return IoStatus::success();
}

std::vector<std::string> PdbWriter::create_conect_pdb_lines(
    const Molecule& molecule) const {
  std::map<int, int> original_to_current;
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    original_to_current[molecule.original_index()[atom]] = molecule.index()[atom];
  }

  std::map<int, std::vector<int>> remapped;
  const auto& conect = molecule.conect();
  for (std::size_t atom = 0; atom < molecule.natoms() && atom < conect.size();
       ++atom) {
    if (conect[atom].empty()) {
      continue;
    }

    const auto base_it = original_to_current.find(molecule.original_index()[atom]);
    if (base_it == original_to_current.end()) {
      continue;
    }

    auto& linked = remapped[base_it->second];
    for (const int original_linked : conect[atom]) {
      const auto linked_it = original_to_current.find(original_linked);
      if (linked_it != original_to_current.end()) {
        linked.push_back(linked_it->second);
      }
    }
  }

  std::vector<std::string> lines;
  for (const auto& [base, linked] : remapped) {
    std::string line = "CONECT";
    auto append_index = [&line](int value) {
      const std::string text = std::to_string(value);
      line.append(5 - std::min<std::size_t>(5, text.size()), ' ');
      line += text;
    };
    append_index(base);
    for (const int value : linked) {
      append_index(value);
    }
    lines.push_back(line);
  }
  return lines;
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

  const int nframes =
      read_int32_from_bytes(header_block.data() + 4, reverse_endian);
  const int istart = read_int32_from_bytes(header_block.data() + 8, reverse_endian);
  const int nsavc = read_int32_from_bytes(header_block.data() + 12, reverse_endian);
  const int namnf = read_int32_from_bytes(header_block.data() + 36, reverse_endian);
  if (nframes < 0 || nsavc < 0 || namnf < 0) {
    return {IoCode::format_error, "Invalid negative DCD header count."};
  }

  DcdHeader parsed;
  parsed.nframes = static_cast<std::size_t>(nframes);
  parsed.istart = istart;
  parsed.nsavc = nsavc;
  parsed.namnf = namnf;
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
  if (header_.natoms >
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) /
          sizeof(float)) {
    return {IoCode::unsupported,
            "DCD coordinate block is too large for int32 record markers."};
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
  (void)close_dcd_write();
  if (options.include_unit_cell) {
    return {IoCode::unsupported,
            "DCD writing with unit-cell blocks is not yet implemented."};
  }
  filename_ = filename;
  options_ = options;
  stream_.open(filename_, std::ios::binary | std::ios::trunc);
  if (!stream_) {
    filename_.clear();
    return {IoCode::file_error, "Failed to open DCD file for writing: " +
                                    filename.string()};
  }
  open_ = true;
  return IoStatus::success();
}

IoStatus DcdWriter::write_dcd_header(const Molecule& molecule,
                                     std::size_t nframes) {
  if (!open_) {
    return {IoCode::not_open, "DCD writer is not open."};
  }
  if (molecule.natoms() >
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) ||
      nframes >
          static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
    return {IoCode::unsupported, "DCD writer supports int32-sized headers."};
  }

  const std::int32_t marker84 = 84;
  const std::int32_t marker164 = 164;
  const std::int32_t marker4 = 4;
  const std::int32_t zero = 0;
  const std::int32_t one = 1;
  const std::int32_t two = 2;
  const std::int32_t natoms = static_cast<std::int32_t>(molecule.natoms());
  const std::int32_t frame_count = static_cast<std::int32_t>(nframes);
  const double delta = 1.0;

  auto status = write_exact(stream_, marker84, "writing DCD header marker");
  if (!status) return status;
  status = write_bytes(stream_, "CORD", 4, "writing DCD CORD magic");
  if (!status) return status;
  status = write_exact(stream_, frame_count, "writing DCD NSET");
  if (!status) return status;
  status = write_exact(stream_, zero, "writing DCD ISTART");
  if (!status) return status;
  status = write_exact(stream_, one, "writing DCD NSAVC");
  if (!status) return status;
  for (int i = 0; i < 6; ++i) {
    status = write_exact(stream_, zero, "writing DCD reserved header integer");
    if (!status) return status;
  }
  status = write_exact(stream_, delta, "writing DCD DELTA");
  if (!status) return status;
  for (int i = 0; i < 9; ++i) {
    status = write_exact(stream_, zero, "writing DCD reserved header integer");
    if (!status) return status;
  }
  status = write_exact(stream_, marker84, "writing trailing DCD header marker");
  if (!status) return status;

  status = write_exact(stream_, marker164, "writing DCD title block marker");
  if (!status) return status;
  status = write_exact(stream_, two, "writing DCD title count");
  if (!status) return status;
  const auto title1 = padded_title("REMARKS FILENAME=A.DCD :: SASSIE");
  const auto title2 =
      padded_title("REMARKS CREATED BY STANDALONE C++ SASMOL");
  status =
      write_bytes(stream_, title1.data(), title1.size(), "writing DCD title 1");
  if (!status) return status;
  status =
      write_bytes(stream_, title2.data(), title2.size(), "writing DCD title 2");
  if (!status) return status;
  status = write_exact(stream_, marker164,
                       "writing trailing DCD title block marker");
  if (!status) return status;

  status = write_exact(stream_, marker4, "writing DCD atom-count marker");
  if (!status) return status;
  status = write_exact(stream_, natoms, "writing DCD atom count");
  if (!status) return status;
  return write_exact(stream_, marker4, "writing trailing DCD atom-count marker");
}

IoStatus DcdWriter::write_dcd_step(const Molecule& molecule, std::size_t frame,
                                   std::size_t step) {
  if (!open_) {
    return {IoCode::not_open, "DCD writer is not open."};
  }
  if (frame >= molecule.number_of_frames()) {
    return {IoCode::format_error, "DCD frame index is out of range."};
  }
  if (molecule.natoms() >
      static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()) /
          sizeof(float)) {
    return {IoCode::unsupported, "DCD writer supports int32-sized atom counts."};
  }

  const auto natoms = molecule.natoms();
  const std::int32_t block_size =
      static_cast<std::int32_t>(natoms * sizeof(float));
  std::vector<float> x(natoms);
  std::vector<float> y(natoms);
  std::vector<float> z(natoms);

  for (std::size_t atom = 0; atom < natoms; ++atom) {
    const auto xyz = molecule.coordinate(frame, atom);
    x[atom] = xyz.x;
    y[atom] = xyz.y;
    z[atom] = xyz.z;
  }

  auto write_axis = [&](const std::vector<float>& axis,
                        const std::string& label) -> IoStatus {
    auto status = write_exact(stream_, block_size,
                              "writing DCD " + label + " block marker");
    if (!status) return status;
    status = write_bytes(stream_, reinterpret_cast<const char*>(axis.data()),
                         block_size, "writing DCD " + label + " coordinates");
    if (!status) return status;
    return write_exact(stream_, block_size,
                       "writing trailing DCD " + label + " block marker");
  };

  auto status = write_axis(x, "X");
  if (!status) return status;
  status = write_axis(y, "Y");
  if (!status) return status;
  status = write_axis(z, "Z");
  if (!status) return status;

  if (step > 0 && step <= static_cast<std::size_t>(INT32_MAX)) {
    const std::int32_t written_step = static_cast<std::int32_t>(step);
    stream_.seekp(8, std::ios::beg);
    status = write_exact(stream_, written_step, "updating DCD frame count");
    if (!status) return status;
    stream_.seekp(20, std::ios::beg);
    status = write_exact(stream_, written_step, "updating DCD step count");
    if (!status) return status;
    stream_.seekp(0, std::ios::end);
  }

  return IoStatus::success();
}

IoStatus DcdWriter::close_dcd_write() {
  if (stream_.is_open()) {
    stream_.close();
  }
  open_ = false;
  filename_.clear();
  return IoStatus::success();
}

IoStatus DcdWriter::write_dcd(const std::filesystem::path& filename,
                              const Molecule& molecule,
                              const DcdWriteOptions& options) {
  auto status = open_dcd_write(filename, options);
  if (!status) {
    return status;
  }

  status = write_dcd_header(molecule, molecule.number_of_frames());
  if (!status) {
    (void)close_dcd_write();
    return status;
  }

  for (std::size_t frame = 0; frame < molecule.number_of_frames(); ++frame) {
    status = write_dcd_step(molecule, frame, frame + 1);
    if (!status) {
      (void)close_dcd_write();
      return status;
    }
  }

  return close_dcd_write();
}

}  // namespace sasmol
