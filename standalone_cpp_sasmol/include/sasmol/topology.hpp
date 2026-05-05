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

struct CharmmTypeChargeAssignment {
  std::string atom_name;
  std::string charmm_type;
  calc_type atom_charge{};
};

[[nodiscard]] SubsetResult assign_charmm_types(
    Molecule& molecule, const std::vector<std::string>& types);

[[nodiscard]] SubsetResult assign_atom_charges(
    Molecule& molecule, const std::vector<calc_type>& charges);

[[nodiscard]] SubsetResult assign_charmm_types_from_atom_table(
    Molecule& molecule, const std::vector<CharmmTypeAssignment>& assignments);

[[nodiscard]] SubsetResult assign_charmm_types_and_atom_charges_from_atom_table(
    Molecule& molecule,
    const std::vector<CharmmTypeChargeAssignment>& assignments);

}  // namespace sasmol
