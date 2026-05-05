#pragma once

#include "sasmol/molecule.hpp"
#include "sasmol/subset.hpp"

#include <cstddef>
#include <filesystem>
#include <map>
#include <optional>
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

struct CharmmTopologyDeleteRecords {
  std::vector<std::string> atoms;
  std::vector<std::vector<std::string>> angles;
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
  CharmmTopologyDeleteRecords deletes;
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

struct CharmmResidueAtomListResult {
  std::map<std::string, std::vector<std::string>> residue_atoms;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CharmmPatchResult {
  CharmmTopologyEntry patched_entry;
  std::vector<std::string> atom_names;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CharmmResidueOrderResult {
  std::string topology_residue_name;
  std::vector<std::string> atom_order;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CharmmResidueReorderPlan {
  std::vector<std::size_t> observed_indices;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CharmmResidueReorderSelection {
  std::string segname;
  int resid{};
  std::string resname;
  std::string topology_residue_name;
  std::vector<std::size_t> atom_indices;
  std::vector<std::string> observed_atom_names;
  std::vector<std::string> topology_atom_order;
  std::vector<std::size_t> observed_indices;
  std::vector<std::size_t> source_atom_indices;
};

struct CharmmMoleculeReorderPlan {
  std::vector<CharmmResidueReorderSelection> residues;
  std::vector<std::size_t> source_atom_indices;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CharmmReorderedMoleculeResult {
  Molecule molecule;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

enum class FastaSplitMode {
  none,
  by_chain,
  by_segname,
};

enum class ConstraintBasis {
  backbone,
  heavy,
  protein,
  nucleic,
  solute,
};

enum class ConstraintField {
  beta,
  occupancy,
};

struct FastaOptions {
  bool fasta_format{false};
  bool exclude_hetatm{false};
  FastaSplitMode split_mode{FastaSplitMode::none};
  std::string name;
  std::size_t width{80};
};

struct FastaResult {
  std::vector<std::string> sequence;
  std::string formatted;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct RenumberOptions {
  std::optional<int> index_start;
  std::optional<int> resid_start;
};

struct ConstraintPdbOptions {
  ConstraintField field{ConstraintField::beta};
  bool reset{true};
};

[[nodiscard]] FastaResult create_fasta(const Molecule& molecule,
                                       const FastaOptions& options = {});

[[nodiscard]] SubsetResult create_fasta_in_place(
    Molecule& molecule,
    const FastaOptions& options = {});

[[nodiscard]] SubsetResult renumber(Molecule& molecule,
                                    const RenumberOptions& options = {});

[[nodiscard]] SubsetResult apply_constraint_descriptor(
    Molecule& molecule,
    ConstraintBasis basis,
    const ConstraintPdbOptions& options = {});

[[nodiscard]] SubsetResult make_constraint_pdb(
    Molecule& molecule,
    const std::filesystem::path& filename,
    ConstraintBasis basis,
    const ConstraintPdbOptions& options = {});

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

[[nodiscard]] bool compare_list_ignore_order(
    const std::vector<std::string>& first,
    const std::vector<std::string>& second);

[[nodiscard]] std::vector<std::string> setup_cys_patch_atoms_simple(
    const std::vector<std::string>& cys_atom_names);

[[nodiscard]] CharmmResidueAtomListResult setup_charmm_residue_atoms(
    const CharmmTopologyData& topology);

[[nodiscard]] CharmmPatchResult patch_charmm_residue_atoms(
    const CharmmTopologyData& topology,
    const std::string& residue,
    const std::string& patch);

[[nodiscard]] CharmmResidueOrderResult choose_charmm_residue_atom_order(
    const CharmmTopologyData& topology,
    const std::map<std::string, std::vector<std::string>>& residue_atoms,
    const std::string& residue_name,
    int residue_id,
    int segment_n_terminal_residue_id,
    int segment_c_terminal_residue_id,
    const std::vector<std::string>& observed_atom_names);

[[nodiscard]] CharmmResidueReorderPlan plan_charmm_residue_reorder_indices(
    const std::vector<std::string>& observed_atom_names,
    const std::vector<std::string>& topology_atom_order,
    const std::string& residue_name,
    int residue_id);

[[nodiscard]] CharmmMoleculeReorderPlan plan_charmm_molecule_reorder(
    const Molecule& molecule,
    const CharmmTopologyData& topology,
    const std::map<std::string, std::vector<std::string>>& residue_atoms);

[[nodiscard]] CharmmReorderedMoleculeResult copy_reordered_charmm_molecule(
    const Molecule& molecule,
    const CharmmMoleculeReorderPlan& plan);

[[nodiscard]] SubsetResult reorder_charmm_molecule_in_place(
    Molecule& molecule,
    const CharmmTopologyData& topology,
    const std::map<std::string, std::vector<std::string>>& residue_atoms);

[[nodiscard]] CharmmTopologyParseResult parse_charmm_topology(
    const std::filesystem::path& filename);

[[nodiscard]] CharmmTopologyParseResult parse_charmm_topology_globals(
    const std::filesystem::path& filename);

}  // namespace sasmol
