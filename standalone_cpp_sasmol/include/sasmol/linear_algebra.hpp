#pragma once

#include "sasmol/calculate.hpp"
#include "sasmol/types.hpp"

#include <vector>

namespace sasmol {

using Matrix = std::vector<std::vector<calc_type>>;

[[nodiscard]] CalcVec3 cross_product(CalcVec3 a, CalcVec3 b);
[[nodiscard]] Matrix matrix_multiply(const Matrix& a, const Matrix& b);
[[nodiscard]] Matrix find_u(const Matrix& x, const Matrix& y);
[[nodiscard]] CalcVec3 vec_sub(CalcVec3 b, CalcVec3 c);
[[nodiscard]] CalcVec3 vec_scale(calc_type scale, CalcVec3 value);
[[nodiscard]] int cmp(calc_type a, calc_type b);
[[nodiscard]] calc_type signed_angle(CalcVec3 a, CalcVec3 b, CalcVec3 c);
[[nodiscard]] calc_type dihedral_angle(CalcVec3 a1, CalcVec3 a2, CalcVec3 a3,
                                       CalcVec3 a4);
[[nodiscard]] calc_type calculate_angle(CalcVec3 a, CalcVec3 b, CalcVec3 c);

}  // namespace sasmol
