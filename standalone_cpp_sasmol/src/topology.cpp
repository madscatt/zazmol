#include "sasmol/topology.hpp"

#include <cstddef>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>

namespace sasmol {

namespace {

std::string line_context(std::size_t line_number, const std::string& message) {
  return "line " + std::to_string(line_number) + ": " + message;
}

std::vector<std::string> split_tokens(const std::string& line) {
  std::istringstream stream(line);
  std::vector<std::string> tokens;
  std::string token;
  while (stream >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

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

CharmmTopologyParseResult parse_charmm_topology_impl(
    const std::filesystem::path& filename, bool include_entries) {
  CharmmTopologyParseResult result;
  std::ifstream input(filename);
  if (!input) {
    result.errors.push_back("failed to open CHARMM topology file: " +
                            filename.string());
    return result;
  }

  std::string line;
  std::size_t line_number = 0;
  CharmmTopologyEntry* current_entry = nullptr;
  while (std::getline(input, line)) {
    ++line_number;
    const auto tokens = split_tokens(line);
    if (tokens.empty() || tokens.front().starts_with("!")) {
      continue;
    }

    const auto& record = tokens.front();
    if (record == "MASS") {
      if (tokens.size() < 4) {
        result.errors.push_back(line_context(
            line_number, "MASS record requires index, atom type, and mass"));
        continue;
      }
      result.topology.masses.push_back({tokens[1], tokens[2], tokens[3]});
    } else if (record == "DECL") {
      if (tokens.size() < 2) {
        result.errors.push_back(
            line_context(line_number, "DECL record requires one token"));
        continue;
      }
      result.topology.declarations.push_back(tokens[1]);
    } else if (record == "DEFA") {
      if (tokens.size() < 5) {
        result.errors.push_back(
            line_context(line_number, "DEFA record requires four tokens"));
        continue;
      }
      result.topology.defaults.insert(result.topology.defaults.end(),
                                      tokens.begin() + 1, tokens.begin() + 5);
    } else if (record == "AUTO") {
      if (tokens.size() < 3) {
        result.errors.push_back(
            line_context(line_number, "AUTO record requires two tokens"));
        continue;
      }
      result.topology.auto_terms.insert(result.topology.auto_terms.end(),
                                        tokens.begin() + 1,
                                        tokens.begin() + 3);
    } else if (include_entries && (record == "RESI" || record == "PRES")) {
      if (tokens.size() < 3) {
        result.errors.push_back(line_context(
            line_number, record + " record requires name and total charge"));
        current_entry = nullptr;
        continue;
      }
      result.topology.entries.push_back(
          {.kind = record == "RESI" ? CharmmTopologyEntryKind::Residue
                                    : CharmmTopologyEntryKind::Patch,
           .name = tokens[1],
           .total_charge = tokens[2],
           .atoms = {}});
      current_entry = &result.topology.entries.back();
    } else if (include_entries && record == "ATOM") {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, "ATOM record appeared before RESI or PRES"));
        continue;
      }
      if (tokens.size() < 4) {
        result.errors.push_back(line_context(
            line_number,
            "ATOM record requires atom name, CHARMM type, and charge"));
        continue;
      }
      current_entry->atoms.push_back({tokens[1], tokens[2], tokens[3]});
    }
  }

  return result;
}

}  // namespace

CharmmTopologyParseResult parse_charmm_topology_globals(
    const std::filesystem::path& filename) {
  return parse_charmm_topology_impl(filename, false);
}

CharmmTopologyParseResult parse_charmm_topology(
    const std::filesystem::path& filename) {
  return parse_charmm_topology_impl(filename, true);
}

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
