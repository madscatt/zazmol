#include "sasmol/linear_algebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

using Matrix3 = std::array<std::array<calc_type, 3>, 3>;

Matrix3 identity_matrix3() {
  return {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
}

calc_type max_abs_diagonal(const Matrix3& matrix) {
  calc_type value{};
  for (std::size_t i = 0; i < 3; ++i) {
    value = std::max(value, std::fabs(matrix[i][i]));
  }
  return value;
}

std::pair<std::array<calc_type, 3>, Matrix3> symmetric_eigen_decomposition(
    const Matrix3& input) {
  auto matrix = input;
  auto eigenvectors = identity_matrix3();

  for (std::size_t iteration = 0; iteration < 60; ++iteration) {
    std::size_t p = 0;
    std::size_t q = 1;
    auto offdiag = std::fabs(matrix[p][q]);
    for (std::size_t row = 0; row < 3; ++row) {
      for (std::size_t column = row + 1; column < 3; ++column) {
        const auto value = std::fabs(matrix[row][column]);
        if (value > offdiag) {
          offdiag = value;
          p = row;
          q = column;
        }
      }
    }

    const auto scale = std::max<calc_type>(1.0, max_abs_diagonal(matrix));
    if (offdiag <= std::numeric_limits<calc_type>::epsilon() * scale) {
      break;
    }

    const auto app = matrix[p][p];
    const auto aqq = matrix[q][q];
    const auto apq = matrix[p][q];
    const auto tau = (aqq - app) / (2.0 * apq);
    const auto t = (tau >= 0.0 ? 1.0 : -1.0) /
                   (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
    const auto c = 1.0 / std::sqrt(1.0 + t * t);
    const auto s = t * c;

    for (std::size_t k = 0; k < 3; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const auto mkp = matrix[k][p];
      const auto mkq = matrix[k][q];
      matrix[k][p] = matrix[p][k] = c * mkp - s * mkq;
      matrix[k][q] = matrix[q][k] = s * mkp + c * mkq;
    }

    matrix[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    matrix[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    matrix[p][q] = matrix[q][p] = 0.0;

    for (std::size_t k = 0; k < 3; ++k) {
      const auto vkp = eigenvectors[k][p];
      const auto vkq = eigenvectors[k][q];
      eigenvectors[k][p] = c * vkp - s * vkq;
      eigenvectors[k][q] = s * vkp + c * vkq;
    }
  }

  return {{{matrix[0][0], matrix[1][1], matrix[2][2]}}, eigenvectors};
}

void validate_coordinate_matrix(const Matrix& value, const char* label) {
  if (value.empty()) {
    throw std::invalid_argument(std::string(label) + " must not be empty");
  }
  for (const auto& row : value) {
    if (row.size() != 3) {
      throw std::invalid_argument(std::string(label) +
                                  " must have three columns");
    }
  }
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

Matrix find_u(const Matrix& x, const Matrix& y) {
  validate_coordinate_matrix(x, "find_u x");
  validate_coordinate_matrix(y, "find_u y");
  if (x.size() != y.size()) {
    throw std::invalid_argument("find_u x and y must have matching rows");
  }

  Matrix3 r{};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      for (std::size_t n = 0; n < x.size(); ++n) {
        r[i][j] += y[n][i] * x[n][j];
      }
    }
  }

  Matrix3 rtr{};
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      for (std::size_t k = 0; k < 3; ++k) {
        rtr[i][j] += r[k][i] * r[k][j];
      }
    }
  }

  auto [eigenvalues, eigenvectors_by_column] =
      symmetric_eigen_decomposition(rtr);
  std::array<std::size_t, 3> order{0, 1, 2};
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return eigenvalues[lhs] > eigenvalues[rhs];
  });

  Matrix3 ak{};
  std::array<calc_type, 3> uk{};
  for (std::size_t row = 0; row < 3; ++row) {
    uk[row] = eigenvalues[order[row]];
    for (std::size_t column = 0; column < 3; ++column) {
      ak[row][column] = eigenvectors_by_column[column][order[row]];
    }
  }

  const CalcVec3 cross_row = cross_product({ak[0][0], ak[0][1], ak[0][2]},
                                           {ak[1][0], ak[1][1], ak[1][2]});
  ak[2] = {cross_row.x, cross_row.y, cross_row.z};

  Matrix3 bk{};
  for (std::size_t basis = 0; basis < 2; ++basis) {
    std::array<calc_type, 3> rak{};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        rak[i] += r[i][k] * ak[basis][k];
      }
    }
    const auto scale = (uk[basis] == 0.0)
                           ? calc_type{1.0e15}
                           : 1.0 / std::sqrt(std::fabs(uk[basis]));
    for (std::size_t i = 0; i < 3; ++i) {
      bk[basis][i] = scale * rak[i];
    }
  }

  const CalcVec3 bk_cross =
      cross_product({bk[0][0], bk[0][1], bk[0][2]},
                    {bk[1][0], bk[1][1], bk[1][2]});
  bk[2] = {bk_cross.x, bk_cross.y, bk_cross.z};

  Matrix u(3, std::vector<calc_type>(3, 0.0));
  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k) {
        u[j][i] += bk[k][i] * ak[k][j];
      }
    }
  }
  return u;
}

CalcVec3 vec_sub(CalcVec3 b, CalcVec3 c) {
  return {b.x - c.x, b.y - c.y, b.z - c.z};
}

CalcVec3 vec_scale(calc_type scale, CalcVec3 value) {
  return {scale * value.x, scale * value.y, scale * value.z};
}

int cmp(calc_type a, calc_type b) {
  return compare(a, b);
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
