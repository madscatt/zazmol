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
  assert(sasmol::cmp(3.0, 2.0) == 1);
  assert(sasmol::cmp(2.0, 3.0) == -1);
  assert(sasmol::cmp(3.0, 3.0) == 0);

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

void assert_matrix_close(const sasmol::Matrix& actual,
                         const sasmol::Matrix& expected,
                         double tolerance = 1.0e-2) {
  assert(actual.size() == expected.size());
  for (std::size_t row = 0; row < expected.size(); ++row) {
    assert(actual[row].size() == expected[row].size());
    for (std::size_t column = 0; column < expected[row].size(); ++column) {
      assert_close(actual[row][column], expected[row][column], tolerance);
    }
  }
}

void test_find_u_matches_python_mathematica_examples() {
  assert_matrix_close(
      sasmol::find_u({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
                     {{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}),
      {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});

  assert_matrix_close(
      sasmol::find_u({{1.0, 1.0, 1.0}, {2.0, 1.0, 1.0}},
                     {{1.0, 1.0, 1.0}, {1.0, 2.0, 1.0}}),
      {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}});

  assert_matrix_close(
      sasmol::find_u({{2.920, -2.367, 1.693},
                      {-0.770, -0.827, -0.417},
                      {-2.150, 3.193, -1.277}},
                     {{1.663, -1.170, 3.567},
                      {-1.197, -1.460, -0.523},
                      {-0.467, 2.630, -3.043}}),
      {{0.902737, -0.0539463, 0.426797},
       {0.23798, 0.8891, -0.390981},
       {-0.358373, 0.454522, 0.815462}});

  assert_matrix_close(
      sasmol::find_u({{0.357, -12.123, 2.098},
                      {1.209, 10.209, -50.082},
                      {-1.098, 3.572, 2.982},
                      {1.231, -1.230, 0.589},
                      {12.398, -30.289, 19.482},
                      {12.123, 0.980, 19.309}},
                     {{90.380, 12.987, 0.392},
                      {3.219, 83.390, 0.028},
                      {0.002, 10.298, -18.820},
                      {12.879, -10.298, 0.987},
                      {0.986, 12.984, 0.367},
                      {12.359, -12.402, 1.298}}),
      {{0.121253, 0.025345, 0.992298},
       {-0.992602, 0.00937959, 0.12105},
       {-0.00623933, -0.999635, 0.0262948}});
}

void test_find_u_rejects_bad_shapes() {
  bool threw = false;
  try {
    (void)sasmol::find_u({{1.0, 2.0}}, {{1.0, 2.0, 3.0}});
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);

  threw = false;
  try {
    (void)sasmol::find_u({{1.0, 2.0, 3.0}},
                         {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}});
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  test_cross_product_matches_python_example();
  test_matrix_multiply_matches_extension_contract();
  test_matrix_multiply_rejects_bad_shapes();
  test_vector_helpers_and_angles();
  test_find_u_matches_python_mathematica_examples();
  test_find_u_rejects_bad_shapes();
  return 0;
}
