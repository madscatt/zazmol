#include "sasmol_direct_usage.h"

#include "sasmol/file_io.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <exception>
#include <limits>
#include <string>

namespace {

void copy_message(char* output, std::size_t output_size,
                  const std::string& message) {
  if (output == nullptr || output_size == 0) {
    return;
  }
  const auto length = std::min(output_size - 1, message.size());
  std::memcpy(output, message.data(), length);
  output[length] = '\0';
}

bool fits_int(std::size_t value) {
  return value <= static_cast<std::size_t>(std::numeric_limits<int>::max());
}

}  // namespace

extern "C" int sasmol_pdb_to_dcd(const char* pdb_path, const char* dcd_path,
                                  char* error_message,
                                  std::size_t error_message_size,
                                  int* atom_count, int* frame_count) {
  if (atom_count != nullptr) {
    *atom_count = 0;
  }
  if (frame_count != nullptr) {
    *frame_count = 0;
  }
  copy_message(error_message, error_message_size, "");

  if (pdb_path == nullptr || dcd_path == nullptr) {
    copy_message(error_message, error_message_size,
                 "pdb_path and dcd_path must be non-null");
    return 0;
  }

  try {
    sasmol::PdbReader pdb_reader;
    sasmol::Molecule molecule;
    auto status = pdb_reader.read_pdb(pdb_path, molecule);
    if (!status) {
      copy_message(error_message, error_message_size,
                   "read_pdb failed: " + status.message);
      return 0;
    }

    sasmol::DcdWriter dcd_writer;
    status = dcd_writer.write_dcd(dcd_path, molecule);
    if (!status) {
      copy_message(error_message, error_message_size,
                   "write_dcd failed: " + status.message);
      return 0;
    }

    if (!fits_int(molecule.natoms()) || !fits_int(molecule.number_of_frames())) {
      copy_message(error_message, error_message_size,
                   "molecule counts exceed the C example's int output range");
      return 0;
    }

    if (atom_count != nullptr) {
      *atom_count = static_cast<int>(molecule.natoms());
    }
    if (frame_count != nullptr) {
      *frame_count = static_cast<int>(molecule.number_of_frames());
    }
    return 1;
  } catch (const std::exception& exc) {
    copy_message(error_message, error_message_size, exc.what());
    return 0;
  } catch (...) {
    copy_message(error_message, error_message_size, "unknown C++ exception");
    return 0;
  }
}
