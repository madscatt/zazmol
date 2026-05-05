#include "sasmol/topology.hpp"

namespace sasmol {

SubsetResult assign_charmm_types(Molecule& molecule,
                                 const std::vector<std::string>& types) {
  SubsetResult result;
  if (types.size() != molecule.natoms()) {
    result.errors.push_back("charmm_type length does not match natoms");
    return result;
  }
  molecule.charmm_type() = types;
  return result;
}

SubsetResult assign_atom_charges(Molecule& molecule,
                                 const std::vector<calc_type>& charges) {
  SubsetResult result;
  if (charges.size() != molecule.natoms()) {
    result.errors.push_back("atom_charge length does not match natoms");
    return result;
  }
  molecule.atom_charge() = charges;
  return result;
}

SubsetResult assign_charmm_types_from_atom_table(
    Molecule& molecule, const std::vector<CharmmTypeAssignment>& assignments) {
  SubsetResult result;
  if (assignments.size() != molecule.natoms()) {
    result.errors.push_back("CHARMM atom table length does not match natoms");
    return result;
  }
  if (molecule.name().size() != molecule.natoms()) {
    result.errors.push_back("atom name length does not match natoms");
    return result;
  }

  std::vector<std::string> types;
  types.reserve(assignments.size());
  for (std::size_t atom = 0; atom < assignments.size(); ++atom) {
    if (assignments[atom].atom_name != molecule.name()[atom]) {
      result.errors.push_back("CHARMM atom table name mismatch at atom " +
                              std::to_string(atom));
      return result;
    }
    types.push_back(assignments[atom].charmm_type);
  }

  molecule.charmm_type() = types;
  return result;
}

SubsetResult assign_charmm_types_and_atom_charges_from_atom_table(
    Molecule& molecule,
    const std::vector<CharmmTypeChargeAssignment>& assignments) {
  SubsetResult result;
  if (assignments.size() != molecule.natoms()) {
    result.errors.push_back("CHARMM atom table length does not match natoms");
    return result;
  }
  if (molecule.name().size() != molecule.natoms()) {
    result.errors.push_back("atom name length does not match natoms");
    return result;
  }

  std::vector<std::string> types;
  std::vector<calc_type> charges;
  types.reserve(assignments.size());
  charges.reserve(assignments.size());
  for (std::size_t atom = 0; atom < assignments.size(); ++atom) {
    if (assignments[atom].atom_name != molecule.name()[atom]) {
      result.errors.push_back("CHARMM atom table name mismatch at atom " +
                              std::to_string(atom));
      return result;
    }
    types.push_back(assignments[atom].charmm_type);
    charges.push_back(assignments[atom].atom_charge);
  }

  molecule.charmm_type() = std::move(types);
  molecule.atom_charge() = std::move(charges);
  return result;
}

}  // namespace sasmol
