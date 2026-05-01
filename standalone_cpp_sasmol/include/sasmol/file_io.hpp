#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>

namespace sasmol {

enum class IoCode {
  ok,
  not_open,
  end_of_file,
  file_error,
  format_error,
  unsupported,
  not_implemented,
};

struct IoStatus {
  IoCode code{IoCode::ok};
  std::string message;

  [[nodiscard]] bool ok() const noexcept { return code == IoCode::ok; }
  [[nodiscard]] explicit operator bool() const noexcept { return ok(); }

  [[nodiscard]] static IoStatus success();
  [[nodiscard]] static IoStatus not_implemented(std::string message);
};

struct PdbReadOptions {
  bool tolerant{true};
  bool resolve_elements{true};
  bool preserve_conect{true};
  bool apply_all_zero_coordinate_guard{true};
};

struct PdbWriteOptions {
  std::size_t frame{0};
  bool write_all_frames{false};
  bool include_conect{true};
};

class PdbReader {
 public:
  [[nodiscard]] IoStatus read_pdb(const std::filesystem::path& filename,
                                  Molecule& molecule,
                                  const PdbReadOptions& options = {}) const;

  [[nodiscard]] static constexpr bool tolerant_by_default() noexcept {
    return true;
  }
};

class PdbWriter {
 public:
  [[nodiscard]] IoStatus write_pdb(const std::filesystem::path& filename,
                                   const Molecule& molecule,
                                   const PdbWriteOptions& options = {}) const;
};

struct DcdHeader {
  std::size_t natoms{};
  std::size_t nframes{};
  int istart{};
  int nsavc{};
  double delta{};
  int namnf{};
  bool reverse_endian{};
  bool has_unit_cell{};
  bool charmm_format{};
  int charmm_flags{};
};

struct DcdReadOptions {
  bool sequential{true};
  bool reopen_for_random_frame{false};
};

struct DcdWriteOptions {
  bool include_unit_cell{false};
};

class DcdReader {
 public:
  DcdReader() = default;

  [[nodiscard]] IoStatus open_dcd_read(const std::filesystem::path& filename,
                                       const DcdReadOptions& options = {});
  [[nodiscard]] IoStatus read_header(DcdHeader& header);
  [[nodiscard]] IoStatus read_next_frame(Molecule& molecule);
  [[nodiscard]] IoStatus close_dcd_read();

  [[nodiscard]] IoStatus read_single_dcd_step(
      const std::filesystem::path& filename, std::size_t frame,
      Molecule& molecule, const DcdReadOptions& options = {});

  [[nodiscard]] bool is_open() const noexcept { return open_; }
  [[nodiscard]] static constexpr bool sequential_by_default() noexcept {
    return true;
  }

 private:
  std::filesystem::path filename_;
  DcdReadOptions options_;
  std::ifstream stream_;
  DcdHeader header_;
  std::size_t current_frame_{};
  bool header_read_{false};
  bool open_{false};
};

class DcdWriter {
 public:
  DcdWriter() = default;

  [[nodiscard]] IoStatus open_dcd_write(const std::filesystem::path& filename,
                                        const DcdWriteOptions& options = {});
  [[nodiscard]] IoStatus write_dcd_header(const Molecule& molecule,
                                          std::size_t nframes);
  [[nodiscard]] IoStatus write_dcd_step(const Molecule& molecule,
                                        std::size_t frame, std::size_t step);
  [[nodiscard]] IoStatus close_dcd_write();

  [[nodiscard]] bool is_open() const noexcept { return open_; }

 private:
  std::filesystem::path filename_;
  DcdWriteOptions options_;
  bool open_{false};
};

}  // namespace sasmol
