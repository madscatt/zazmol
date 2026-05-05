#include "sasmol/view.hpp"

#include <limits>

namespace sasmol {
namespace {

#ifdef SASMOL_ENABLE_VMD_ADAPTER
extern "C" int send_coordinates_to_vmd(int n_atoms, float* x, float* y,
                                       float* z, int port,
                                       int flag_clear_sock);

int legacy_vmd_sender(std::span<const coord_type> x,
                      std::span<const coord_type> y,
                      std::span<const coord_type> z, int port,
                      bool clear_socket) {
  return send_coordinates_to_vmd(
      static_cast<int>(x.size()), const_cast<coord_type*>(x.data()),
      const_cast<coord_type*>(y.data()), const_cast<coord_type*>(z.data()),
      port, clear_socket ? 1 : 0);
}
#endif

ViewStatus validate_coordinate_arrays(std::span<const coord_type> x,
                                      std::span<const coord_type> y,
                                      std::span<const coord_type> z) {
  if (x.size() != y.size() || x.size() != z.size()) {
    return {ViewCode::invalid_argument,
            "VMD coordinate arrays must have matching lengths."};
  }
  if (x.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    return {ViewCode::unsupported,
            "VMD coordinate arrays exceed int-sized atom counts."};
  }
  return {};
}

}  // namespace

VmdCoordinateArrayResult prepare_vmd_coordinate_arrays(
    const Molecule& molecule, std::size_t frame) {
  VmdCoordinateArrayResult result;
  if (frame >= molecule.number_of_frames()) {
    result.status = {ViewCode::invalid_argument,
                     "VMD frame index is out of range."};
    return result;
  }
  const auto expected_size = molecule.natoms() * molecule.number_of_frames() * 3;
  if (molecule.coor().size() != expected_size) {
    result.status = {ViewCode::invalid_argument,
                     "Molecule coordinate storage does not match shape."};
    return result;
  }

  result.coordinates.x.resize(molecule.natoms());
  result.coordinates.y.resize(molecule.natoms());
  result.coordinates.z.resize(molecule.natoms());
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto coordinate = molecule.coordinate(frame, atom);
    result.coordinates.x[atom] = coordinate.x;
    result.coordinates.y[atom] = coordinate.y;
    result.coordinates.z[atom] = coordinate.z;
  }
  return result;
}

ViewStatus send_coordinate_arrays_to_vmd(std::span<const coord_type> x,
                                         std::span<const coord_type> y,
                                         std::span<const coord_type> z,
                                         int port, bool clear_socket) {
#ifdef SASMOL_ENABLE_VMD_ADAPTER
  return send_coordinate_arrays_to_vmd(x, y, z, port, clear_socket,
                                       legacy_vmd_sender);
#else
  (void)port;
  (void)clear_socket;
  const auto validation = validate_coordinate_arrays(x, y, z);
  if (!validation.ok()) {
    return validation;
  }
  return {ViewCode::unsupported,
          "VMD adapter is disabled in this standalone C++ build."};
#endif
}

ViewStatus send_coordinate_arrays_to_vmd(std::span<const coord_type> x,
                                         std::span<const coord_type> y,
                                         std::span<const coord_type> z,
                                         int port, bool clear_socket,
                                         VmdCoordinateSender sender) {
  const auto validation = validate_coordinate_arrays(x, y, z);
  if (!validation.ok()) {
    return validation;
  }
  if (sender == nullptr) {
    return {ViewCode::invalid_argument, "VMD coordinate sender is null."};
  }
  const auto result = sender(x, y, z, port, clear_socket);
  if (result != 0) {
    return {ViewCode::transport_error,
            "VMD coordinate sender returned a nonzero status."};
  }
  return {};
}

ViewStatus send_coordinates_to_vmd(const Molecule& molecule, int port,
                                   bool clear_socket, std::size_t frame) {
  const auto arrays = prepare_vmd_coordinate_arrays(molecule, frame);
  if (!arrays.ok()) {
    return arrays.status;
  }
  return send_coordinate_arrays_to_vmd(arrays.coordinates.x,
                                       arrays.coordinates.y,
                                       arrays.coordinates.z, port,
                                       clear_socket);
}

ViewStatus send_coordinates_to_vmd(const Molecule& molecule, int port,
                                   bool clear_socket, std::size_t frame,
                                   VmdCoordinateSender sender) {
  const auto arrays = prepare_vmd_coordinate_arrays(molecule, frame);
  if (!arrays.ok()) {
    return arrays.status;
  }
  return send_coordinate_arrays_to_vmd(arrays.coordinates.x,
                                       arrays.coordinates.y,
                                       arrays.coordinates.z, port,
                                       clear_socket, sender);
}

}  // namespace sasmol
