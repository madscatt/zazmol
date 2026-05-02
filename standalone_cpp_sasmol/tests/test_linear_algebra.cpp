#include "sasmol/linear_algebra.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

void assert_close(double actual, double expected, double tolerance = 1.0e-6) {
  assert(std::fabs(actual - expected) <= tolerance);
}

void test_cross_product_matches_python_example() {
  const auto result =
      sasmol::cross_product({1.0, 2.0, 3.0}, {-1.0, 6.0, 8.0});

  assert_close(result.x, -2.0);
  assert_close(result.y, -11.0);
  assert_close(result.z, 8.0);
}

void test_matrix_multiply_matches_extension_contract() {
  const sasmol::Matrix a{{5.0, 3.0, 1.0}, {2.0, 3.0, 5.0}};
  const sasmol::Matrix b{{2.0}, {-4.0}, {8.0}};

  const auto result = sasmol::matrix_multiply(a, b);

  assert(result.size() == 2);
  assert(result[0].size() == 1);
  assert_close(result[0][0], 6.0);
  assert_close(result[1][0], 32.0);
}

void test_matrix_multiply_rejects_bad_shapes() {
  bool threw = false;
  try {
    (void)sasmol::matrix_multiply({{1.0, 2.0}}, {{1.0, 2.0}});
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

void test_vector_helpers_and_angles() {
  const auto difference =
      sasmol::vec_sub({1.0, 2.0, 3.0}, {-1.0, 6.0, 8.0});
  assert_close(difference.x, 2.0);
  assert_close(difference.y, -4.0);
  assert_close(difference.z, -5.0);

  const auto scaled = sasmol::vec_scale(2.0, {-1.0, 6.0, 8.0});
  assert_close(scaled.x, -2.0);
  assert_close(scaled.y, 12.0);
  assert_close(scaled.z, 16.0);

  const auto signed_angle = sasmol::signed_angle(
      {1.0, 2.0, 3.0}, {-1.0, 6.0, 8.0}, {-4.0, -1.0, 4.0});
  assert_close(signed_angle, 21.444512921997863, 1.0e-10);

  const auto dihedral = sasmol::dihedral_angle(
      {1.0, 2.0, 3.0}, {-1.0, 6.0, 8.0}, {-4.0, -1.0, 4.0},
      {-3.0, -41.0, 3.0});
  assert_close(dihedral, 85.950635659264, 1.0e-10);

  const auto angle = sasmol::calculate_angle({1.0, 0.0, 0.0},
                                             {0.0, 0.0, 0.0},
                                             {0.0, 1.0, 0.0});
  assert_close(angle, 1.5707963267948966, 1.0e-12);
}

}  // namespace

int main() {
  test_cross_product_matches_python_example();
  test_matrix_multiply_matches_extension_contract();
  test_matrix_multiply_rejects_bad_shapes();
  test_vector_helpers_and_angles();
  return 0;
}
