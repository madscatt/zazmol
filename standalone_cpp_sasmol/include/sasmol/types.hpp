#pragma once

#include <cstddef>
#include <span>

namespace sasmol {

using coord_type = float;
using calc_type = double;

struct Vec3 {
  coord_type x{};
  coord_type y{};
  coord_type z{};
};

struct ConstCoordinateView {
  std::span<const coord_type> xyz;
  std::size_t natoms{};

  [[nodiscard]] Vec3 operator[](std::size_t atom) const;
};

struct CoordinateView {
  std::span<coord_type> xyz;
  std::size_t natoms{};

  [[nodiscard]] Vec3 operator[](std::size_t atom) const;
  void set(std::size_t atom, Vec3 value);
};

}  // namespace sasmol
