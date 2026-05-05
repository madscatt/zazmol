#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/types.hpp"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace sasmol {

enum class ViewCode {
  ok,
  invalid_argument,
  unsupported,
  transport_error,
};

struct ViewStatus {
  ViewCode code{ViewCode::ok};
  std::string message;

  [[nodiscard]] bool ok() const noexcept { return code == ViewCode::ok; }
};

struct VmdCoordinateArrays {
  std::vector<coord_type> x;
  std::vector<coord_type> y;
  std::vector<coord_type> z;
};

struct VmdCoordinateArrayResult {
  VmdCoordinateArrays coordinates;
  ViewStatus status;

  [[nodiscard]] bool ok() const noexcept { return status.ok(); }
};

using VmdCoordinateSender =
    int (*)(std::span<const coord_type> x, std::span<const coord_type> y,
            std::span<const coord_type> z, int port, bool clear_socket);

[[nodiscard]] VmdCoordinateArrayResult prepare_vmd_coordinate_arrays(
    const Molecule& molecule, std::size_t frame = 0);

[[nodiscard]] ViewStatus send_coordinate_arrays_to_vmd(
    std::span<const coord_type> x, std::span<const coord_type> y,
    std::span<const coord_type> z, int port, bool clear_socket);

[[nodiscard]] ViewStatus send_coordinate_arrays_to_vmd(
    std::span<const coord_type> x, std::span<const coord_type> y,
    std::span<const coord_type> z, int port, bool clear_socket,
    VmdCoordinateSender sender);

[[nodiscard]] ViewStatus send_coordinates_to_vmd(const Molecule& molecule,
                                                 int port, bool clear_socket,
                                                 std::size_t frame = 0);

[[nodiscard]] ViewStatus send_coordinates_to_vmd(const Molecule& molecule,
                                                 int port, bool clear_socket,
                                                 std::size_t frame,
                                                 VmdCoordinateSender sender);

}  // namespace sasmol
