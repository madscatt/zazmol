#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/types.hpp"

#include <cstddef>
#include <vector>

namespace sasmol {

[[nodiscard]] bool has_overlap(const std::vector<Vec3>& first,
                               const std::vector<Vec3>& second,
                               calc_type cutoff);
[[nodiscard]] bool has_overlap(const Molecule& first, std::size_t first_frame,
                               const Molecule& second,
                               std::size_t second_frame, calc_type cutoff);

}  // namespace sasmol
