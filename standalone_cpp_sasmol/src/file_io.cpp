#include "sasmol/file_io.hpp"

#include <utility>

namespace sasmol {

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
  filename_ = filename;
  options_ = options;
  open_ = true;
  return IoStatus::not_implemented(
      "DCD opening is intentionally deferred until the C-extension binary "
      "contract is captured in parity tests.");
}

IoStatus DcdReader::read_header(DcdHeader& header) {
  (void)header;
  if (!open_) {
    return {IoCode::not_open, "DCD reader is not open."};
  }
  return IoStatus::not_implemented(
      "DCD header parsing is intentionally deferred until CHARMM/NAMD/unit-cell "
      "variants are captured in parity tests.");
}

IoStatus DcdReader::read_next_frame(Molecule& molecule) {
  (void)molecule;
  if (!open_) {
    return {IoCode::not_open, "DCD reader is not open."};
  }
  return IoStatus::not_implemented(
      "DCD frame reading is intentionally deferred; the primary API is "
      "sequential read_next_frame().");
}

IoStatus DcdReader::close_dcd_read() {
  open_ = false;
  filename_.clear();
  return IoStatus::success();
}

IoStatus DcdReader::read_single_dcd_step(const std::filesystem::path& filename,
                                         std::size_t frame, Molecule& molecule,
                                         const DcdReadOptions& options) {
  (void)filename;
  (void)frame;
  (void)molecule;
  (void)options;
  return IoStatus::not_implemented(
      "Random-frame DCD access will be implemented, if needed, as explicit "
      "reopen/scan behavior rather than hidden seeking.");
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
