#include "sasmol/topology.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <iterator>
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

bool is_charmm_atom_reference_token(const std::string& token) {
  if (token.empty()) {
    return false;
  }
  const auto first = static_cast<unsigned char>(token.front());
  return std::isalnum(first) || token.front() == '+' || token.front() == '-';
}

void parse_pair_record(const std::vector<std::string>& tokens,
                       const std::string& record_name,
                       std::size_t line_number,
                       std::vector<CharmmTopologyPairRecord>& records,
                       std::vector<std::string>& errors) {
  std::size_t token = 1;
  while (token + 1 < tokens.size() && !tokens[token].starts_with("!")) {
    if (!is_charmm_atom_reference_token(tokens[token]) ||
        !is_charmm_atom_reference_token(tokens[token + 1])) {
      errors.push_back(line_context(line_number,
                                    record_name + " record has invalid pair"));
      return;
    }
    records.push_back({tokens[token], tokens[token + 1]});
    token += 2;
  }
  if (token < tokens.size() && !tokens[token].starts_with("!")) {
    errors.push_back(line_context(line_number,
                                  record_name + " record has incomplete pair"));
  }
}

void parse_single_record(const std::vector<std::string>& tokens,
                         const std::string& record_name,
                         std::size_t line_number,
                         std::vector<std::string>& records,
                         std::vector<std::string>& errors) {
  std::size_t token = 1;
  while (token < tokens.size() && !tokens[token].starts_with("!")) {
    if (!is_charmm_atom_reference_token(tokens[token])) {
      errors.push_back(line_context(line_number,
                                    record_name + " record has invalid token"));
      return;
    }
    records.push_back(tokens[token]);
    ++token;
  }
}

void parse_triple_record(const std::vector<std::string>& tokens,
                         const std::string& record_name,
                         std::size_t line_number,
                         std::vector<CharmmTopologyTripleRecord>& records,
                         std::vector<std::string>& errors) {
  std::size_t token = 1;
  while (token + 2 < tokens.size() && !tokens[token].starts_with("!")) {
    if (!is_charmm_atom_reference_token(tokens[token]) ||
        !is_charmm_atom_reference_token(tokens[token + 1]) ||
        !is_charmm_atom_reference_token(tokens[token + 2])) {
      errors.push_back(line_context(line_number,
                                    record_name + " record has invalid triple"));
      return;
    }
    records.push_back({tokens[token], tokens[token + 1], tokens[token + 2]});
    token += 3;
  }
  if (token < tokens.size() && !tokens[token].starts_with("!")) {
    errors.push_back(line_context(
        line_number, record_name + " record has incomplete triple"));
  }
}

void parse_quad_record(const std::vector<std::string>& tokens,
                       const std::string& record_name,
                       std::size_t line_number,
                       std::vector<CharmmTopologyQuadRecord>& records,
                       std::vector<std::string>& errors) {
  std::size_t token = 1;
  while (token + 3 < tokens.size() && !tokens[token].starts_with("!")) {
    if (!is_charmm_atom_reference_token(tokens[token]) ||
        !is_charmm_atom_reference_token(tokens[token + 1]) ||
        !is_charmm_atom_reference_token(tokens[token + 2]) ||
        !is_charmm_atom_reference_token(tokens[token + 3])) {
      errors.push_back(line_context(line_number,
                                    record_name + " record has invalid quad"));
      return;
    }
    records.push_back({tokens[token], tokens[token + 1], tokens[token + 2],
                       tokens[token + 3]});
    token += 4;
  }
  if (token < tokens.size() && !tokens[token].starts_with("!")) {
    errors.push_back(
        line_context(line_number, record_name + " record has incomplete quad"));
  }
}

void parse_internal_coordinate_record(
    const std::vector<std::string>& tokens,
    std::vector<CharmmTopologyInternalCoordinateRecord>& records) {
  auto end = tokens.end();
  if (tokens.size() > 10) {
    end = tokens.begin() + 10;
  }
  records.push_back({std::vector<std::string>(tokens.begin() + 1, end)});
}

void parse_delete_record(const std::vector<std::string>& tokens,
                         std::size_t line_number,
                         CharmmTopologyDeleteRecords& records,
                         std::vector<std::string>& errors) {
  if (tokens.size() < 2) {
    errors.push_back(
        line_context(line_number, "DELE record requires a delete type"));
    return;
  }

  const auto delete_type = tokens[1].substr(0, 4);
  if (delete_type == "ATOM") {
    if (tokens.size() < 3) {
      errors.push_back(
          line_context(line_number, "DELE ATOM record requires atom name"));
      return;
    }
    records.atoms.push_back(tokens[2]);
  } else if (delete_type == "ANGL") {
    records.angles.emplace_back(tokens.begin() + 2, tokens.end());
  }
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

template <typename T>
void require_atom_aligned(const Molecule& molecule,
                          const std::vector<T>& values,
                          const std::string& field,
                          std::vector<std::string>& errors) {
  if (values.size() != molecule.natoms()) {
    errors.push_back(field + " length does not match natoms");
  }
}

template <typename T>
void require_optional_atom_aligned(const Molecule& molecule,
                                   const std::vector<T>& values,
                                   const std::string& field,
                                   std::vector<std::string>& errors) {
  if (!values.empty() && values.size() != molecule.natoms()) {
    errors.push_back(field + " length does not match natoms");
  }
}

template <typename T>
void copy_required_reordered_vector(const std::vector<T>& source,
                                    std::vector<T>& destination,
                                    const std::vector<std::size_t>& order) {
  for (std::size_t atom = 0; atom < order.size(); ++atom) {
    destination[atom] = source[order[atom]];
  }
}

template <typename T>
void copy_optional_reordered_vector(const std::vector<T>& source,
                                    std::vector<T>& destination,
                                    const std::vector<std::size_t>& order) {
  if (source.empty()) {
    destination.clear();
    return;
  }
  copy_required_reordered_vector(source, destination, order);
}

template <typename T>
void copy_reordered_descriptor_map(
    const std::map<std::string, std::vector<T>>& source,
    std::map<std::string, std::vector<T>>& destination,
    const std::vector<std::size_t>& order) {
  destination.clear();
  for (const auto& [name, values] : source) {
    auto& output = destination[name];
    output.resize(order.size());
    for (std::size_t atom = 0; atom < order.size(); ++atom) {
      output[atom] = values[order[atom]];
    }
  }
}

const CharmmTopologyEntry* find_topology_entry(
    const CharmmTopologyData& topology, const std::string& name) {
  for (const auto& entry : topology.entries) {
    if (entry.name == name) {
      return &entry;
    }
  }
  return nullptr;
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
           .atoms = {},
           .bonds = {},
           .doubles = {},
           .angles = {},
           .thetas = {},
           .dihedrals = {},
           .impropers = {},
           .cmaps = {},
           .donors = {},
           .acceptors = {},
           .internal_coordinates = {},
           .deletes = {}});
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
    } else if (include_entries && (record.starts_with("BOND") ||
                                   record.starts_with("DOUB"))) {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, record + " record appeared before RESI or PRES"));
        continue;
      }
      auto& records = record.starts_with("BOND") ? current_entry->bonds
                                                 : current_entry->doubles;
      parse_pair_record(tokens, record.starts_with("BOND") ? "BOND" : "DOUB",
                        line_number, records, result.errors);
    } else if (include_entries && (record.starts_with("ANGL") ||
                                   record.starts_with("THET"))) {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, record + " record appeared before RESI or PRES"));
        continue;
      }
      auto& records = record.starts_with("ANGL") ? current_entry->angles
                                                 : current_entry->thetas;
      parse_triple_record(tokens, record.starts_with("ANGL") ? "ANGL" : "THET",
                          line_number, records, result.errors);
    } else if (include_entries && (record.starts_with("DIHE") ||
                                   record.starts_with("IMPR") ||
                                   record.starts_with("CMAP"))) {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, record + " record appeared before RESI or PRES"));
        continue;
      }
      auto* records = &current_entry->cmaps;
      std::string record_name = "CMAP";
      if (record.starts_with("DIHE")) {
        records = &current_entry->dihedrals;
        record_name = "DIHE";
      } else if (record.starts_with("IMPR")) {
        records = &current_entry->impropers;
        record_name = "IMPR";
      }
      parse_quad_record(tokens, record_name, line_number, *records,
                        result.errors);
    } else if (include_entries && (record.starts_with("DONO") ||
                                   record.starts_with("ACCE"))) {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, record + " record appeared before RESI or PRES"));
        continue;
      }
      auto& records = record.starts_with("DONO") ? current_entry->donors
                                                 : current_entry->acceptors;
      parse_single_record(tokens, record.starts_with("DONO") ? "DONO" : "ACCE",
                          line_number, records, result.errors);
    } else if (include_entries && record == "IC") {
      if (current_entry == nullptr) {
        result.errors.push_back(
            line_context(line_number, "IC record appeared before RESI or PRES"));
        continue;
      }
      parse_internal_coordinate_record(tokens,
                                       current_entry->internal_coordinates);
    } else if (include_entries && record.starts_with("DELE")) {
      if (current_entry == nullptr) {
        result.errors.push_back(line_context(
            line_number, "DELE record appeared before RESI or PRES"));
        continue;
      }
      parse_delete_record(tokens, line_number, current_entry->deletes,
                          result.errors);
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

bool compare_list_ignore_order(const std::vector<std::string>& first,
                               const std::vector<std::string>& second) {
  if (first.size() != second.size()) {
    return false;
  }
  for (const auto& item : first) {
    if (std::count(second.begin(), second.end(), item) != 1) {
      return false;
    }
  }
  return true;
}

std::vector<std::string> setup_cys_patch_atoms_simple(
    const std::vector<std::string>& cys_atom_names) {
  std::vector<std::string> atoms;
  atoms.reserve(cys_atom_names.size());
  for (const auto& atom_name : cys_atom_names) {
    if (atom_name != "HG1") {
      atoms.push_back(atom_name);
    }
  }
  return atoms;
}

CharmmResidueAtomListResult setup_charmm_residue_atoms(
    const CharmmTopologyData& topology) {
  CharmmResidueAtomListResult result;
  for (const auto& entry : topology.entries) {
    if (entry.atoms.empty()) {
      continue;
    }
    auto& atom_names = result.residue_atoms[entry.name];
    atom_names.reserve(entry.atoms.size());
    for (const auto& atom : entry.atoms) {
      atom_names.push_back(atom.name);
    }
  }

  const auto cys = result.residue_atoms.find("CYS");
  if (cys != result.residue_atoms.end()) {
    result.residue_atoms["DISU"] = setup_cys_patch_atoms_simple(cys->second);
  }
  return result;
}

CharmmPatchResult patch_charmm_residue_atoms(const CharmmTopologyData& topology,
                                             const std::string& residue,
                                             const std::string& patch) {
  CharmmPatchResult result;
  const auto* residue_entry = find_topology_entry(topology, residue);
  if (residue_entry == nullptr) {
    result.errors.push_back("Residue " + residue +
                            " not found in the topology information list "
                            "during patching!");
    return result;
  }
  const auto* patch_entry = find_topology_entry(topology, patch);
  if (patch_entry == nullptr) {
    result.errors.push_back("Patch " + patch +
                            " not found in the topology information list "
                            "during patching!");
    return result;
  }

  result.patched_entry = *residue_entry;
  result.patched_entry.name = residue + "_" + patch;

  for (const auto& atom_delete : patch_entry->deletes.atoms) {
    for (auto atom = result.patched_entry.atoms.begin();
         atom != result.patched_entry.atoms.end();) {
      if (atom->name == atom_delete) {
        atom = result.patched_entry.atoms.erase(atom);
      } else {
        ++atom;
      }
    }
  }

  std::size_t num_added = 0;
  for (const auto& patch_atom : patch_entry->atoms) {
    for (auto atom = result.patched_entry.atoms.begin();
         atom != result.patched_entry.atoms.end();) {
      if (atom->name == patch_atom.name) {
        atom = result.patched_entry.atoms.erase(atom);
      } else {
        ++atom;
      }
    }

    if (patch == "NTER" || patch == "GLYP" || patch == "PROP") {
      const auto insert_at =
          std::min(num_added, result.patched_entry.atoms.size());
      result.patched_entry.atoms.insert(result.patched_entry.atoms.begin() +
                                            static_cast<std::ptrdiff_t>(insert_at),
                                        patch_atom);
    } else if (patch == "CTER") {
      result.patched_entry.atoms.push_back(patch_atom);
    }
    ++num_added;
  }

  result.atom_names.reserve(result.patched_entry.atoms.size());
  for (const auto& atom : result.patched_entry.atoms) {
    result.atom_names.push_back(atom.name);
  }
  return result;
}

CharmmResidueOrderResult choose_charmm_residue_atom_order(
    const CharmmTopologyData& topology,
    const std::map<std::string, std::vector<std::string>>& residue_atoms,
    const std::string& residue_name,
    int residue_id,
    int segment_n_terminal_residue_id,
    int segment_c_terminal_residue_id,
    const std::vector<std::string>& observed_atom_names) {
  CharmmResidueOrderResult result;
  auto working_topology = topology;
  auto working_residue_atoms = residue_atoms;
  auto topology_residue_name = residue_name;

  auto atom_order_for = [&](const std::string& name)
      -> const std::vector<std::string>* {
    const auto found = working_residue_atoms.find(name);
    if (found == working_residue_atoms.end()) {
      return nullptr;
    }
    return &found->second;
  };

  auto apply_patch = [&](const std::string& patch_name) -> bool {
    auto patch_result =
        patch_charmm_residue_atoms(working_topology, topology_residue_name,
                                   patch_name);
    if (!patch_result.ok()) {
      result.errors.insert(result.errors.end(), patch_result.errors.begin(),
                           patch_result.errors.end());
      return false;
    }
    topology_residue_name += "_" + patch_name;
    working_topology.entries.push_back(std::move(patch_result.patched_entry));
    working_residue_atoms[topology_residue_name] =
        std::move(patch_result.atom_names);
    return true;
  };

  auto* atom_order = atom_order_for(topology_residue_name);
  if (atom_order == nullptr) {
    result.errors.push_back("For residue: " + topology_residue_name +
                            "\nthe atom names doesn't match those in the "
                            "charmm topology file!\n[]\n");
    return result;
  }

  if (!compare_list_ignore_order(*atom_order, observed_atom_names)) {
    if (topology_residue_name == "CYS") {
      topology_residue_name = "DISU";
      atom_order = atom_order_for(topology_residue_name);
    }
    if (topology_residue_name == "HIS") {
      topology_residue_name = "HSE";
      atom_order = atom_order_for(topology_residue_name);
      if (atom_order == nullptr ||
          !compare_list_ignore_order(*atom_order, observed_atom_names)) {
        topology_residue_name = "HSD";
        atom_order = atom_order_for(topology_residue_name);
        if (atom_order == nullptr ||
            !compare_list_ignore_order(*atom_order, observed_atom_names)) {
          topology_residue_name = "HSP";
          atom_order = atom_order_for(topology_residue_name);
        }
      }
    }

    if (residue_id == segment_n_terminal_residue_id) {
      std::string patch = "NTER";
      if (topology_residue_name == "GLY") {
        patch = "GLYP";
      } else if (topology_residue_name == "PRO" ||
                 topology_residue_name == "PROP") {
        patch = "PROP";
      }
      if (!apply_patch(patch)) {
        return result;
      }
      atom_order = atom_order_for(topology_residue_name);
    }
    if (residue_id == segment_c_terminal_residue_id) {
      if (!apply_patch("CTER")) {
        return result;
      }
      atom_order = atom_order_for(topology_residue_name);
    }
    if (atom_order == nullptr ||
        !compare_list_ignore_order(*atom_order, observed_atom_names)) {
      result.errors.push_back("For residue: " + topology_residue_name +
                              "\nthe atom names doesn't match those in the "
                              "charmm topology file!");
      return result;
    }
  }

  atom_order = atom_order_for(topology_residue_name);
  if (atom_order == nullptr) {
    result.errors.push_back("For residue: " + topology_residue_name +
                            "\nthe atom names doesn't match those in the "
                            "charmm topology file!");
    return result;
  }
  result.topology_residue_name = topology_residue_name;
  result.atom_order = *atom_order;
  return result;
}

CharmmResidueReorderPlan plan_charmm_residue_reorder_indices(
    const std::vector<std::string>& observed_atom_names,
    const std::vector<std::string>& topology_atom_order,
    const std::string& residue_name,
    int residue_id) {
  CharmmResidueReorderPlan result;
  result.observed_indices.reserve(topology_atom_order.size());
  for (const auto& atom_name : topology_atom_order) {
    const auto found =
        std::find(observed_atom_names.begin(), observed_atom_names.end(),
                  atom_name);
    if (found == observed_atom_names.end()) {
      result.errors.push_back("Atom " + atom_name +
                              " not found while planning CHARMM order for "
                              "resname: " +
                              residue_name + "\nand resid: " +
                              std::to_string(residue_id) + "!");
      result.observed_indices.clear();
      return result;
    }
    result.observed_indices.push_back(static_cast<std::size_t>(
        std::distance(observed_atom_names.begin(), found)));
  }

  if (observed_atom_names.size() != result.observed_indices.size()) {
    result.errors.push_back(
        "Number of atoms doesnt match that in charmm topology file for "
        "\nresname: " +
        residue_name + "\nand resid: " + std::to_string(residue_id) + "!");
    result.observed_indices.clear();
    return result;
  }
  return result;
}

CharmmMoleculeReorderPlan plan_charmm_molecule_reorder(
    const Molecule& molecule,
    const CharmmTopologyData& topology,
    const std::map<std::string, std::vector<std::string>>& residue_atoms) {
  CharmmMoleculeReorderPlan result;
  require_atom_aligned(molecule, molecule.name(), "name", result.errors);
  require_atom_aligned(molecule, molecule.resname(), "resname", result.errors);
  require_atom_aligned(molecule, molecule.resid(), "resid", result.errors);
  require_atom_aligned(molecule, molecule.segname(), "segname", result.errors);
  if (!result.errors.empty()) {
    return result;
  }

  struct ResidueGroup {
    std::string segname;
    int resid{};
    std::string resname;
    std::vector<std::size_t> atom_indices;
    std::vector<std::string> atom_names;
  };

  struct SegmentGroup {
    std::string segname;
    std::vector<ResidueGroup> residues;
  };

  std::vector<SegmentGroup> segments;
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    const auto& segname = molecule.segname()[atom];
    const auto resid = molecule.resid()[atom];
    const auto& resname = molecule.resname()[atom];

    auto segment = std::find_if(
        segments.begin(), segments.end(), [&](const SegmentGroup& candidate) {
          return candidate.segname == segname;
        });
    if (segment == segments.end()) {
      segments.push_back({.segname = segname, .residues = {}});
      segment = std::prev(segments.end());
    }

    auto residue = std::find_if(
        segment->residues.begin(), segment->residues.end(),
        [&](const ResidueGroup& candidate) { return candidate.resid == resid; });
    if (residue == segment->residues.end()) {
      segment->residues.push_back({.segname = segname,
                                   .resid = resid,
                                   .resname = resname,
                                   .atom_indices = {},
                                   .atom_names = {}});
      residue = std::prev(segment->residues.end());
    } else if (residue->resname != resname) {
      result.errors.push_back("Residue " + std::to_string(resid) +
                              " in segment " + segname +
                              " has mixed residue names");
      return result;
    }

    residue->atom_indices.push_back(atom);
    residue->atom_names.push_back(molecule.name()[atom]);
  }

  result.source_atom_indices.reserve(molecule.natoms());
  for (const auto& segment : segments) {
    if (segment.residues.empty()) {
      continue;
    }
    const auto n_terminal_resid = segment.residues.front().resid;
    const auto c_terminal_resid = segment.residues.back().resid;

    for (const auto& residue : segment.residues) {
      auto order = choose_charmm_residue_atom_order(
          topology, residue_atoms, residue.resname, residue.resid,
          n_terminal_resid, c_terminal_resid, residue.atom_names);
      if (!order.ok()) {
        result.errors.insert(result.errors.end(), order.errors.begin(),
                             order.errors.end());
        result.source_atom_indices.clear();
        result.residues.clear();
        return result;
      }

      auto residue_plan = plan_charmm_residue_reorder_indices(
          residue.atom_names, order.atom_order, order.topology_residue_name,
          residue.resid);
      if (!residue_plan.ok()) {
        result.errors.insert(result.errors.end(), residue_plan.errors.begin(),
                             residue_plan.errors.end());
        result.source_atom_indices.clear();
        result.residues.clear();
        return result;
      }

      CharmmResidueReorderSelection selection;
      selection.segname = residue.segname;
      selection.resid = residue.resid;
      selection.resname = residue.resname;
      selection.topology_residue_name = order.topology_residue_name;
      selection.atom_indices = residue.atom_indices;
      selection.observed_atom_names = residue.atom_names;
      selection.topology_atom_order = order.atom_order;
      selection.observed_indices = residue_plan.observed_indices;
      selection.source_atom_indices.reserve(residue_plan.observed_indices.size());
      for (const auto observed_index : residue_plan.observed_indices) {
        const auto source_atom = residue.atom_indices[observed_index];
        selection.source_atom_indices.push_back(source_atom);
        result.source_atom_indices.push_back(source_atom);
      }
      result.residues.push_back(std::move(selection));
    }
  }

  if (result.source_atom_indices.size() != molecule.natoms()) {
    result.errors.push_back("Planned CHARMM reorder atom count does not match "
                            "natoms");
    result.source_atom_indices.clear();
    result.residues.clear();
  }
  return result;
}

CharmmReorderedMoleculeResult copy_reordered_charmm_molecule(
    const Molecule& molecule,
    const CharmmMoleculeReorderPlan& plan) {
  CharmmReorderedMoleculeResult result;
  if (!plan.ok()) {
    result.errors.insert(result.errors.end(), plan.errors.begin(),
                         plan.errors.end());
    return result;
  }
  if (plan.source_atom_indices.size() != molecule.natoms()) {
    result.errors.push_back("CHARMM reorder plan length does not match natoms");
    return result;
  }
  for (const auto source_atom : plan.source_atom_indices) {
    if (source_atom >= molecule.natoms()) {
      result.errors.push_back("CHARMM reorder plan contains out-of-range atom "
                              "index");
      return result;
    }
  }

  auto integrity = molecule.check_integrity();
  if (!integrity.ok()) {
    for (const auto& issue : integrity.issues) {
      result.errors.push_back(issue.field + " length does not match natoms");
    }
    return result;
  }

  for (const auto& [name, values] : molecule.extra_string_descriptors()) {
    require_optional_atom_aligned(molecule, values,
                                  "extra_string_descriptors." + name,
                                  result.errors);
  }
  for (const auto& [name, values] : molecule.extra_int_descriptors()) {
    require_optional_atom_aligned(molecule, values,
                                  "extra_int_descriptors." + name,
                                  result.errors);
  }
  for (const auto& [name, values] : molecule.extra_calc_descriptors()) {
    require_optional_atom_aligned(molecule, values,
                                  "extra_calc_descriptors." + name,
                                  result.errors);
  }
  if (!result.errors.empty()) {
    return result;
  }

  result.molecule = molecule;
  const auto& order = plan.source_atom_indices;
  copy_required_reordered_vector(molecule.record(), result.molecule.record(),
                                 order);
  copy_required_reordered_vector(molecule.index(), result.molecule.index(),
                                 order);
  copy_required_reordered_vector(molecule.original_index(),
                                 result.molecule.original_index(), order);
  copy_required_reordered_vector(molecule.original_resid(),
                                 result.molecule.original_resid(), order);
  copy_required_reordered_vector(molecule.name(), result.molecule.name(), order);
  copy_required_reordered_vector(molecule.loc(), result.molecule.loc(), order);
  copy_required_reordered_vector(molecule.resname(), result.molecule.resname(),
                                 order);
  copy_required_reordered_vector(molecule.chain(), result.molecule.chain(),
                                 order);
  copy_required_reordered_vector(molecule.resid(), result.molecule.resid(),
                                 order);
  copy_required_reordered_vector(molecule.rescode(), result.molecule.rescode(),
                                 order);
  copy_required_reordered_vector(molecule.occupancy(),
                                 result.molecule.occupancy(), order);
  copy_required_reordered_vector(molecule.beta(), result.molecule.beta(), order);
  copy_required_reordered_vector(molecule.segname(), result.molecule.segname(),
                                 order);
  copy_required_reordered_vector(molecule.element(), result.molecule.element(),
                                 order);
  copy_required_reordered_vector(molecule.charge(), result.molecule.charge(),
                                 order);
  copy_required_reordered_vector(molecule.atom_charge(),
                                 result.molecule.atom_charge(), order);
  copy_required_reordered_vector(molecule.atom_vdw(), result.molecule.atom_vdw(),
                                 order);
  copy_required_reordered_vector(molecule.residue_flag(),
                                 result.molecule.residue_flag(), order);
  copy_optional_reordered_vector(molecule.charmm_type(),
                                 result.molecule.charmm_type(), order);
  copy_required_reordered_vector(molecule.moltype(), result.molecule.moltype(),
                                 order);
  copy_required_reordered_vector(molecule.mass(), result.molecule.mass(), order);
  copy_required_reordered_vector(molecule.residue_charge(),
                                 result.molecule.residue_charge(), order);
  copy_required_reordered_vector(molecule.conect(), result.molecule.conect(),
                                 order);

  copy_reordered_descriptor_map(molecule.extra_string_descriptors(),
                                result.molecule.extra_string_descriptors(),
                                order);
  copy_reordered_descriptor_map(molecule.extra_int_descriptors(),
                                result.molecule.extra_int_descriptors(), order);
  copy_reordered_descriptor_map(molecule.extra_calc_descriptors(),
                                result.molecule.extra_calc_descriptors(),
                                order);

  for (std::size_t frame = 0; frame < molecule.number_of_frames(); ++frame) {
    for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
      result.molecule.set_coordinate(
          frame, atom, molecule.coordinate(frame, order[atom]));
    }
  }

  return result;
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
