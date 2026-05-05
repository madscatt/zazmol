#include "sasmol/topology.hpp"

#include <map>

namespace sasmol {

namespace {

std::map<std::string, int> count_atom_names(const std::vector<std::string>& names) {
  std::map<std::string, int> counts;
  for (const auto& name : names) {
    ++counts[name];
  }
  return counts;
}

std::map<std::string, int> count_topology_atom_names(
    const std::vector<CharmmAtomDefinition>& atoms) {
  std::map<std::string, int> counts;
  for (const auto& atom : atoms) {
    ++counts[atom.name];
  }
  return counts;
}

void append_duplicates(const std::map<std::string, int>& counts,
                       std::vector<std::string>& duplicates) {
  for (const auto& [name, count] : counts) {
    if (count > 1) {
      duplicates.push_back(name);
    }
  }
}

}  // namespace

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

CharmmResidueValidation validate_charmm_residue_atoms(
    const std::vector<std::string>& molecule_atom_names,
    const CharmmResidueDefinition& residue) {
  CharmmResidueValidation validation;
  const auto molecule_counts = count_atom_names(molecule_atom_names);
  const auto topology_counts = count_topology_atom_names(residue.atoms);

  append_duplicates(molecule_counts, validation.duplicate_molecule_atoms);
  append_duplicates(topology_counts, validation.duplicate_topology_atoms);

  for (const auto& [name, count] : topology_counts) {
    (void)count;
    if (molecule_counts.find(name) == molecule_counts.end()) {
      validation.missing_atoms.push_back(name);
    }
  }
  for (const auto& [name, count] : molecule_counts) {
    (void)count;
    if (topology_counts.find(name) == topology_counts.end()) {
      validation.extra_atoms.push_back(name);
    }
  }

  return validation;
}

}  // namespace sasmol
