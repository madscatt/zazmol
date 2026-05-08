#include "sasmol/molecule.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <vector>

namespace {

using Mat3 = std::array<std::array<sasmol::calc_type, 3>, 3>;

struct FixtureTransform {
  Mat3 rotation;
  std::array<sasmol::calc_type, 3> translation;
};

void assert_close(sasmol::coord_type actual, sasmol::coord_type expected,
                  double tolerance = 1.0e-6) {
  assert(std::fabs(static_cast<double>(actual - expected)) < tolerance);
}

void assert_vec_close(sasmol::Vec3 actual, sasmol::Vec3 expected,
                      double tolerance = 1.0e-6) {
  assert_close(actual.x, expected.x, tolerance);
  assert_close(actual.y, expected.y, tolerance);
  assert_close(actual.z, expected.z, tolerance);
}

sasmol::Vec3 apply_transform(const FixtureTransform& transform, sasmol::Vec3 xyz) {
  const auto x = static_cast<sasmol::calc_type>(xyz.x);
  const auto y = static_cast<sasmol::calc_type>(xyz.y);
  const auto z = static_cast<sasmol::calc_type>(xyz.z);

  const auto rx = transform.rotation[0][0] * x + transform.rotation[0][1] * y +
                  transform.rotation[0][2] * z + transform.translation[0];
  const auto ry = transform.rotation[1][0] * x + transform.rotation[1][1] * y +
                  transform.rotation[1][2] * z + transform.translation[1];
  const auto rz = transform.rotation[2][0] * x + transform.rotation[2][1] * y +
                  transform.rotation[2][2] * z + transform.translation[2];

  return {static_cast<sasmol::coord_type>(rx), static_cast<sasmol::coord_type>(ry),
          static_cast<sasmol::coord_type>(rz)};
}

std::vector<sasmol::Vec3> transformed_copy(
    const std::vector<sasmol::Vec3>& source_coordinates,
    const FixtureTransform& transform) {
  std::vector<sasmol::Vec3> output;
  output.reserve(source_coordinates.size());
  for (const auto xyz : source_coordinates) {
    output.push_back(apply_transform(transform, xyz));
  }
  return output;
}

void test_biomt_fixture_transform_order_semantics() {
  const std::vector<sasmol::Vec3> source{
      {1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
  };

  const FixtureTransform translate_x10{
      .rotation = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
      .translation = {10.0, 0.0, 0.0},
  };
  const FixtureTransform rotate_z_90_translate_y5{
      .rotation = {{{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}},
      .translation = {0.0, 5.0, 0.0},
  };

  const auto first = transformed_copy(source, translate_x10);
  const auto second = transformed_copy(source, rotate_z_90_translate_y5);

  assert_vec_close(first[0], {11.0f, 0.0f, 0.0f});
  assert_vec_close(first[1], {10.0f, 1.0f, 0.0f});
  assert_vec_close(first[2], {10.0f, 0.0f, 1.0f});

  assert_vec_close(second[0], {0.0f, 6.0f, 0.0f});
  assert_vec_close(second[1], {-1.0f, 5.0f, 0.0f});
  assert_vec_close(second[2], {0.0f, 5.0f, 1.0f});
}

void test_biomt_fixture_frame_target_semantics() {
  sasmol::Molecule molecule(2, 2);
  molecule.set_coordinate(0, 0, {1.0f, 0.0f, 0.0f});
  molecule.set_coordinate(0, 1, {0.0f, 1.0f, 0.0f});
  molecule.set_coordinate(1, 0, {3.0f, 0.0f, 0.0f});
  molecule.set_coordinate(1, 1, {0.0f, 3.0f, 0.0f});

  const FixtureTransform translate_x2{
      .rotation = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}},
      .translation = {2.0, 0.0, 0.0},
  };

  const std::vector<sasmol::Vec3> frame_one{
      molecule.coordinate(1, 0),
      molecule.coordinate(1, 1),
  };
  const auto transformed = transformed_copy(frame_one, translate_x2);

  assert_vec_close(transformed[0], {5.0f, 0.0f, 0.0f});
  assert_vec_close(transformed[1], {2.0f, 3.0f, 0.0f});
}

void test_biomt_fixture_source_coordinates_unchanged() {
  const std::vector<sasmol::Vec3> source{
      {2.0f, 2.0f, 0.0f},
      {3.0f, -1.0f, 1.0f},
  };
  const auto source_copy = source;

  const FixtureTransform rotate_z_180{
      .rotation = {{{-1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}}},
      .translation = {0.0, 0.0, 0.0},
  };

  const auto transformed = transformed_copy(source, rotate_z_180);

  assert_vec_close(transformed[0], {-2.0f, -2.0f, 0.0f});
  assert_vec_close(transformed[1], {-3.0f, 1.0f, 1.0f});
  assert_vec_close(source[0], source_copy[0]);
  assert_vec_close(source[1], source_copy[1]);
}

}  // namespace

int main() {
  test_biomt_fixture_transform_order_semantics();
  test_biomt_fixture_frame_target_semantics();
  test_biomt_fixture_source_coordinates_unchanged();
  return 0;
}
