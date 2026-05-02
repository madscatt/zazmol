#include "sasmol/linear_algebra.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sasmol {
namespace {

constexpr calc_type pi = 3.141592653589793238462643383279502884;

calc_type dot(CalcVec3 a, CalcVec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

calc_type norm_squared(CalcVec3 value) {
  return dot(value, value);
}

int compare(calc_type a, calc_type b) {
  return (a > b) - (a < b);
}

}  // namespace

CalcVec3 cross_product(CalcVec3 a, CalcVec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

Matrix matrix_multiply(const Matrix& a, const Matrix& b) {
  if (a.empty() || b.empty() || b.front().empty()) {
    throw std::invalid_argument("matrix_multiply requires non-empty matrices");
  }
  const auto inner = a.front().size();
  if (inner == 0 || b.size() != inner) {
    throw std::invalid_argument("matrix_multiply received incompatible matrices");
  }
  for (const auto& row : a) {
    if (row.size() != inner) {
      throw std::invalid_argument("matrix_multiply requires rectangular matrix a");
    }
  }
  const auto columns = b.front().size();
  for (const auto& row : b) {
    if (row.size() != columns) {
      throw std::invalid_argument("matrix_multiply requires rectangular matrix b");
    }
  }

  Matrix result(a.size(), std::vector<calc_type>(columns, 0.0));
  for (std::size_t row = 0; row < a.size(); ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      for (std::size_t k = 0; k < inner; ++k) {
        result[row][column] += a[row][k] * b[k][column];
      }
    }
  }
  return result;
}

CalcVec3 vec_sub(CalcVec3 b, CalcVec3 c) {
  return {b.x - c.x, b.y - c.y, b.z - c.z};
}

CalcVec3 vec_scale(calc_type scale, CalcVec3 value) {
  return {scale * value.x, scale * value.y, scale * value.z};
}

calc_type signed_angle(CalcVec3 a, CalcVec3 b, CalcVec3 c) {
  const auto ada = norm_squared(a);
  const auto bdb = norm_squared(b);
  if (ada * bdb <= 0.0) {
    return 180.0;
  }

  const auto argument = dot(a, b) / std::sqrt(ada * bdb);
  if (!std::isfinite(argument)) {
    return 180.0;
  }
  const auto clamped = std::clamp(argument, -1.0, 1.0);
  const auto angle = (180.0 / pi) * std::acos(clamped);
  const auto cp = cross_product(a, b);
  const auto sign = compare(dot(cp, c), 0.0);
  return static_cast<calc_type>(sign) * angle;
}

calc_type dihedral_angle(CalcVec3 a1, CalcVec3 a2, CalcVec3 a3,
                         CalcVec3 a4) {
  const auto r1 = vec_sub(a1, a2);
  const auto r2 = vec_sub(a3, a2);
  const auto r3 = vec_scale(-1.0, r2);
  const auto r4 = vec_sub(a4, a3);
  const auto n1 = cross_product(r1, r2);
  const auto n2 = cross_product(r3, r4);
  return signed_angle(n1, n2, r2);
}

calc_type calculate_angle(CalcVec3 a, CalcVec3 b, CalcVec3 c) {
  const auto u = vec_sub(a, b);
  const auto v = vec_sub(c, b);
  const auto denominator = std::sqrt(norm_squared(u) * norm_squared(v));
  if (denominator <= 0.0 || !std::isfinite(denominator)) {
    throw std::invalid_argument("calculate_angle requires nonzero vectors");
  }
  const auto cosine = std::clamp(dot(u, v) / denominator, -1.0, 1.0);
  return std::acos(cosine);
}

}  // namespace sasmol
