#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/subset.hpp"

#include <string>
#include <vector>

namespace sasmol {

struct CharmmTypeAssignment {
  std::string atom_name;
  std::string charmm_type;
};

[[nodiscard]] SubsetResult assign_charmm_types(
    Molecule& molecule, const std::vector<std::string>& types);

[[nodiscard]] SubsetResult assign_charmm_types_from_atom_table(
    Molecule& molecule, const std::vector<CharmmTypeAssignment>& assignments);

}  // namespace sasmol
