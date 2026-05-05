#include "sasmol/view.hpp"

#include <cassert>
#include <span>
#include <vector>

namespace {

struct CapturedSend {
  std::vector<sasmol::coord_type> x;
  std::vector<sasmol::coord_type> y;
  std::vector<sasmol::coord_type> z;
  int port{};
  bool clear_socket{};
  int calls{};
};

CapturedSend captured;

int capture_sender(std::span<const sasmol::coord_type> x,
                   std::span<const sasmol::coord_type> y,
                   std::span<const sasmol::coord_type> z, int port,
                   bool clear_socket) {
  captured.x.assign(x.begin(), x.end());
  captured.y.assign(y.begin(), y.end());
  captured.z.assign(z.begin(), z.end());
  captured.port = port;
  captured.clear_socket = clear_socket;
  ++captured.calls;
  return 0;
}

int failing_sender(std::span<const sasmol::coord_type>,
                   std::span<const sasmol::coord_type>,
                   std::span<const sasmol::coord_type>, int, bool) {
  return 7;
}

void test_prepare_vmd_coordinate_arrays_uses_requested_frame() {
  sasmol::Molecule molecule(2, 2);
  molecule.set_coordinate(0, 0, {1.0F, 2.0F, 3.0F});
  molecule.set_coordinate(0, 1, {4.0F, 5.0F, 6.0F});
  molecule.set_coordinate(1, 0, {7.0F, 8.0F, 9.0F});
  molecule.set_coordinate(1, 1, {10.0F, 11.0F, 12.0F});

  const auto result = sasmol::prepare_vmd_coordinate_arrays(molecule, 1);

  assert(result.ok());
  assert((result.coordinates.x == std::vector<sasmol::coord_type>{7.0F, 10.0F}));
  assert((result.coordinates.y == std::vector<sasmol::coord_type>{8.0F, 11.0F}));
  assert((result.coordinates.z == std::vector<sasmol::coord_type>{9.0F, 12.0F}));
}

void test_prepare_vmd_coordinate_arrays_rejects_bad_frame() {
  const sasmol::Molecule molecule(1, 1);

  const auto result = sasmol::prepare_vmd_coordinate_arrays(molecule, 1);

  assert(!result.ok());
  assert(result.status.code == sasmol::ViewCode::invalid_argument);
  assert(result.coordinates.x.empty());
}

void test_send_coordinates_to_vmd_passes_python_style_arrays() {
  captured = {};
  sasmol::Molecule molecule(2, 1);
  molecule.set_coordinate(0, 0, {1.25F, 2.5F, 3.75F});
  molecule.set_coordinate(0, 1, {4.0F, 5.5F, 6.75F});

  const auto status =
      sasmol::send_coordinates_to_vmd(molecule, 1085, false, 0, capture_sender);

  assert(status.ok());
  assert(captured.calls == 1);
  assert(captured.port == 1085);
  assert(!captured.clear_socket);
  assert((captured.x == std::vector<sasmol::coord_type>{1.25F, 4.0F}));
  assert((captured.y == std::vector<sasmol::coord_type>{2.5F, 5.5F}));
  assert((captured.z == std::vector<sasmol::coord_type>{3.75F, 6.75F}));
}

void test_send_coordinate_arrays_rejects_mismatched_lengths() {
  const std::vector<sasmol::coord_type> x{1.0F};
  const std::vector<sasmol::coord_type> y{2.0F, 3.0F};
  const std::vector<sasmol::coord_type> z{4.0F};

  const auto status =
      sasmol::send_coordinate_arrays_to_vmd(x, y, z, 1085, false, capture_sender);

  assert(status.code == sasmol::ViewCode::invalid_argument);
}

void test_send_coordinate_arrays_reports_transport_failure() {
  const std::vector<sasmol::coord_type> x{1.0F};
  const std::vector<sasmol::coord_type> y{2.0F};
  const std::vector<sasmol::coord_type> z{3.0F};

  const auto status =
      sasmol::send_coordinate_arrays_to_vmd(x, y, z, 1085, false, failing_sender);

  assert(status.code == sasmol::ViewCode::transport_error);
}

void test_default_vmd_sender_is_optional_when_disabled() {
#ifndef SASMOL_ENABLE_VMD_ADAPTER
  const std::vector<sasmol::coord_type> x{1.0F};
  const std::vector<sasmol::coord_type> y{2.0F};
  const std::vector<sasmol::coord_type> z{3.0F};

  const auto status = sasmol::send_coordinate_arrays_to_vmd(x, y, z, 1085, false);

  assert(status.code == sasmol::ViewCode::unsupported);
#endif
}

}  // namespace

int main() {
  test_prepare_vmd_coordinate_arrays_uses_requested_frame();
  test_prepare_vmd_coordinate_arrays_rejects_bad_frame();
  test_send_coordinates_to_vmd_passes_python_style_arrays();
  test_send_coordinate_arrays_rejects_mismatched_lengths();
  test_send_coordinate_arrays_reports_transport_failure();
  test_default_vmd_sender_is_optional_when_disabled();
  return 0;
}
