#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/subset.hpp"

#include <filesystem>
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

struct CharmmMassRecord {
  std::string index;
  std::string atom_type;
  std::string mass;
};

enum class CharmmTopologyEntryKind {
  Residue,
  Patch,
};

struct CharmmTopologyAtomRecord {
  std::string name;
  std::string charmm_type;
  std::string atom_charge;
};

struct CharmmTopologyPairRecord {
  std::string first;
  std::string second;
};

struct CharmmTopologyTripleRecord {
  std::string first;
  std::string second;
  std::string third;
};

struct CharmmTopologyQuadRecord {
  std::string first;
  std::string second;
  std::string third;
  std::string fourth;
};

struct CharmmTopologyInternalCoordinateRecord {
  std::vector<std::string> fields;
};

struct CharmmTopologyEntry {
  CharmmTopologyEntryKind kind{};
  std::string name;
  std::string total_charge;
  std::vector<CharmmTopologyAtomRecord> atoms;
  std::vector<CharmmTopologyPairRecord> bonds;
  std::vector<CharmmTopologyPairRecord> doubles;
  std::vector<CharmmTopologyTripleRecord> angles;
  std::vector<CharmmTopologyTripleRecord> thetas;
  std::vector<CharmmTopologyQuadRecord> dihedrals;
  std::vector<CharmmTopologyQuadRecord> impropers;
  std::vector<CharmmTopologyQuadRecord> cmaps;
  std::vector<std::string> donors;
  std::vector<std::string> acceptors;
  std::vector<CharmmTopologyInternalCoordinateRecord> internal_coordinates;
};

struct CharmmTopologyData {
  std::vector<CharmmMassRecord> masses;
  std::vector<std::string> declarations;
  std::vector<std::string> defaults;
  std::vector<std::string> auto_terms;
  std::vector<CharmmTopologyEntry> entries;
};

struct CharmmTopologyParseResult {
  CharmmTopologyData topology;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
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

[[nodiscard]] CharmmTopologyParseResult parse_charmm_topology(
    const std::filesystem::path& filename);

[[nodiscard]] CharmmTopologyParseResult parse_charmm_topology_globals(
    const std::filesystem::path& filename);

}  // namespace sasmol
