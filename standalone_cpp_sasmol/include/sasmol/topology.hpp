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

struct CharmmAtomDefinition {
  std::string name;
  std::string charmm_type;
  calc_type atom_charge{};
};

struct CharmmResidueDefinition {
  std::string resname;
  std::vector<CharmmAtomDefinition> atoms;
};

struct CharmmResidueValidation {
  std::vector<std::string> missing_atoms;
  std::vector<std::string> extra_atoms;
  std::vector<std::string> duplicate_molecule_atoms;
  std::vector<std::string> duplicate_topology_atoms;

  [[nodiscard]] bool ok() const noexcept {
    return missing_atoms.empty() && extra_atoms.empty() &&
           duplicate_molecule_atoms.empty() && duplicate_topology_atoms.empty();
  }
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

[[nodiscard]] CharmmResidueValidation validate_charmm_residue_atoms(
    const std::vector<std::string>& molecule_atom_names,
    const CharmmResidueDefinition& residue);

}  // namespace sasmol
