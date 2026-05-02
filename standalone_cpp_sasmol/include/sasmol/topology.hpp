#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/subset.hpp"

#include <string>
#include <vector>

namespace sasmol {

[[nodiscard]] SubsetResult assign_charmm_types(
    Molecule& molecule, const std::vector<std::string>& types);

}  // namespace sasmol
