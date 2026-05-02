#include "sasmol/subset.hpp"

#include <stdexcept>
#include <set>

namespace sasmol {
namespace {

std::vector<std::string> validate_indices(const Molecule& molecule,
                                          const std::vector<std::size_t>& indices,
                                          std::size_t frame) {
  std::vector<std::string> errors;
  if (frame >= molecule.number_of_frames()) {
    errors.push_back("subset frame is out of range");
    return errors;
  }
  if (indices.empty()) {
    errors.push_back("subset indices are empty");
    return errors;
  }
  for (const auto atom : indices) {
    if (atom >= molecule.natoms()) {
      errors.push_back("subset atom index is out of range");
      return errors;
    }
  }
  return errors;
}

std::vector<std::string> validate_atom_indices(
    const Molecule& molecule, const std::vector<std::size_t>& indices) {
  std::vector<std::string> errors;
  if (indices.empty()) {
    errors.push_back("subset indices are empty");
    return errors;
  }
  for (const auto atom : indices) {
    if (atom >= molecule.natoms()) {
      errors.push_back("subset atom index is out of range");
      return errors;
    }
  }
  return errors;
}

template <typename T>
void copy_selected_vector(const std::vector<T>& source, std::vector<T>& destination,
                          const std::vector<std::size_t>& indices) {
  destination.clear();
  if (source.empty()) {
    return;
  }
  if (source.size() < indices.size()) {
    throw std::invalid_argument("descriptor vector is smaller than selection");
  }
  destination.reserve(indices.size());
  for (const auto atom : indices) {
    if (atom >= source.size()) {
      throw std::invalid_argument("descriptor vector length does not cover selection");
    }
    destination.push_back(source[atom]);
  }
}

std::vector<std::string>& string_descriptor(Molecule& molecule,
                                            StringDescriptor descriptor) {
  switch (descriptor) {
    case StringDescriptor::record:
      return molecule.record();
    case StringDescriptor::name:
      return molecule.name();
    case StringDescriptor::loc:
      return molecule.loc();
    case StringDescriptor::resname:
      return molecule.resname();
    case StringDescriptor::chain:
      return molecule.chain();
    case StringDescriptor::rescode:
      return molecule.rescode();
    case StringDescriptor::occupancy:
      return molecule.occupancy();
    case StringDescriptor::beta:
      return molecule.beta();
    case StringDescriptor::segname:
      return molecule.segname();
    case StringDescriptor::element:
      return molecule.element();
    case StringDescriptor::charge:
      return molecule.charge();
    case StringDescriptor::charmm_type:
      return molecule.charmm_type();
    case StringDescriptor::moltype:
      return molecule.moltype();
  }
  return molecule.name();
}

const std::vector<std::string>& string_descriptor(
    const Molecule& molecule, StringDescriptor descriptor) {
  switch (descriptor) {
    case StringDescriptor::record:
      return molecule.record();
    case StringDescriptor::name:
      return molecule.name();
    case StringDescriptor::loc:
      return molecule.loc();
    case StringDescriptor::resname:
      return molecule.resname();
    case StringDescriptor::chain:
      return molecule.chain();
    case StringDescriptor::rescode:
      return molecule.rescode();
    case StringDescriptor::occupancy:
      return molecule.occupancy();
    case StringDescriptor::beta:
      return molecule.beta();
    case StringDescriptor::segname:
      return molecule.segname();
    case StringDescriptor::element:
      return molecule.element();
    case StringDescriptor::charge:
      return molecule.charge();
    case StringDescriptor::charmm_type:
      return molecule.charmm_type();
    case StringDescriptor::moltype:
      return molecule.moltype();
  }
  return molecule.name();
}

std::vector<int>& int_descriptor(Molecule& molecule, IntDescriptor descriptor) {
  switch (descriptor) {
    case IntDescriptor::index:
      return molecule.index();
    case IntDescriptor::original_index:
      return molecule.original_index();
    case IntDescriptor::original_resid:
      return molecule.original_resid();
    case IntDescriptor::resid:
      return molecule.resid();
    case IntDescriptor::residue_flag:
      return molecule.residue_flag();
  }
  return molecule.index();
}

const std::vector<int>& int_descriptor(const Molecule& molecule,
                                       IntDescriptor descriptor) {
  switch (descriptor) {
    case IntDescriptor::index:
      return molecule.index();
    case IntDescriptor::original_index:
      return molecule.original_index();
    case IntDescriptor::original_resid:
      return molecule.original_resid();
    case IntDescriptor::resid:
      return molecule.resid();
    case IntDescriptor::residue_flag:
      return molecule.residue_flag();
  }
  return molecule.index();
}

std::vector<calc_type>& calc_descriptor(Molecule& molecule,
                                        CalcDescriptor descriptor) {
  switch (descriptor) {
    case CalcDescriptor::atom_charge:
      return molecule.atom_charge();
    case CalcDescriptor::atom_vdw:
      return molecule.atom_vdw();
    case CalcDescriptor::mass:
      return molecule.mass();
  }
  return molecule.mass();
}

const std::vector<calc_type>& calc_descriptor(const Molecule& molecule,
                                              CalcDescriptor descriptor) {
  switch (descriptor) {
    case CalcDescriptor::atom_charge:
      return molecule.atom_charge();
    case CalcDescriptor::atom_vdw:
      return molecule.atom_vdw();
    case CalcDescriptor::mass:
      return molecule.mass();
  }
  return molecule.mass();
}

template <typename T>
std::vector<std::string> validate_descriptor_vector(
    const Molecule& molecule, const std::vector<T>& descriptor) {
  if (descriptor.size() != molecule.natoms()) {
    return {"descriptor length does not match natoms"};
  }
  return {};
}

template <typename T>
void append_vector(std::vector<T>& destination, const std::vector<T>& source) {
  destination.insert(destination.end(), source.begin(), source.end());
}

template <typename T>
void validate_required_vector(const Molecule& molecule,
                              const std::vector<T>& values,
                              const std::string& label,
                              std::vector<std::string>& errors) {
  if (values.size() != molecule.natoms()) {
    errors.push_back(label + " length does not match natoms");
  }
}

template <typename T>
void merge_required_vector(const std::vector<T>& left,
                           const std::vector<T>& right,
                           std::vector<T>& destination) {
  destination.clear();
  destination.reserve(left.size() + right.size());
  append_vector(destination, left);
  append_vector(destination, right);
}

template <typename T>
void merge_optional_vector(const std::vector<T>& left,
                           std::size_t left_natoms,
                           const std::vector<T>& right,
                           std::size_t right_natoms,
                           std::vector<T>& destination,
                           const std::string& label,
                           MergeOptions options,
                           SubsetResult& result) {
  destination.clear();
  const bool left_present = left_natoms == 0 || !left.empty();
  const bool right_present = right_natoms == 0 || !right.empty();
  if (!left_present && !right_present) {
    return;
  }
  if (left_present != right_present) {
    if (options.report_skipped_descriptors) {
      result.errors.push_back("skipped descriptor " + label +
                              " missing from one input molecule");
    }
    return;
  }
  destination.reserve(left.size() + right.size());
  append_vector(destination, left);
  append_vector(destination, right);
}

template <typename T>
void copy_selected_descriptor_map(
    const std::map<std::string, std::vector<T>>& source,
    std::map<std::string, std::vector<T>>& destination,
    const std::vector<std::size_t>& indices) {
  destination.clear();
  for (const auto& [name, values] : source) {
    std::vector<T> selected;
    copy_selected_vector(values, selected, indices);
    destination[name] = selected;
  }
}

template <typename T>
void validate_descriptor_map(const Molecule& molecule,
                             const std::map<std::string, std::vector<T>>& map,
                             const std::string& label,
                             std::vector<std::string>& errors) {
  for (const auto& [name, values] : map) {
    validate_required_vector(molecule, values, label + "." + name, errors);
  }
}

template <typename T>
void merge_descriptor_map(const std::map<std::string, std::vector<T>>& left,
                          std::size_t left_natoms,
                          const std::map<std::string, std::vector<T>>& right,
                          std::size_t right_natoms,
                          std::map<std::string, std::vector<T>>& destination,
                          const std::string& label,
                          MergeOptions options,
                          SubsetResult& result) {
  destination.clear();
  for (const auto& [name, values] : left) {
    const auto found = right.find(name);
    if (found == right.end()) {
      if (right_natoms == 0) {
        destination[name] = values;
        continue;
      }
      if (options.report_skipped_descriptors) {
        result.errors.push_back("skipped descriptor " + label + "." + name +
                                " missing from mol2");
      }
      continue;
    }
    std::vector<T> merged;
    merge_optional_vector(values, left_natoms, found->second, right_natoms,
                          merged, label + "." + name, options, result);
    if (!merged.empty() || (left_natoms == 0 && right_natoms == 0)) {
      destination[name] = merged;
    }
  }
  for (const auto& [name, values] : right) {
    if (!left.contains(name) && options.report_skipped_descriptors) {
      result.errors.push_back("skipped descriptor " + label + "." + name +
                              " missing from mol1");
    }
  }
}

template <typename T>
std::vector<T> get_descriptor_values(const std::vector<T>& descriptor,
                                     const std::vector<std::size_t>& indices) {
  std::vector<T> values;
  values.reserve(indices.size());
  for (const auto atom : indices) {
    values.push_back(descriptor[atom]);
  }
  return values;
}

template <typename T>
SubsetResult set_descriptor_values(std::vector<T>& descriptor,
                                   const Molecule& molecule,
                                   const std::vector<std::size_t>& indices,
                                   const T& value) {
  SubsetResult result;
  result.errors = validate_atom_indices(molecule, indices);
  if (!result.ok()) {
    return result;
  }
  result.errors = validate_descriptor_vector(molecule, descriptor);
  if (!result.ok()) {
    return result;
  }
  for (const auto atom : indices) {
    descriptor[atom] = value;
  }
  return result;
}

void copy_selected_conect(const Molecule& source, Molecule& destination,
                          const std::vector<std::size_t>& indices) {
  destination.conect().assign(indices.size(), {});
  if (source.original_index().size() != source.natoms() ||
      source.conect().empty()) {
    return;
  }

  std::set<int> selected_original_indices;
  for (const auto atom : indices) {
    selected_original_indices.insert(source.original_index()[atom]);
  }

  for (std::size_t selected = 0; selected < indices.size(); ++selected) {
    const auto source_atom = indices[selected];
    if (source_atom >= source.conect().size()) {
      continue;
    }
    for (const auto original_linked : source.conect()[source_atom]) {
      if (selected_original_indices.contains(original_linked)) {
        destination.conect()[selected].push_back(original_linked);
      }
    }
  }
}

std::vector<std::string> validate_merge_source(const Molecule& molecule,
                                               const std::string& label,
                                               bool require_atoms) {
  std::vector<std::string> errors;
  if (require_atoms && molecule.natoms() == 0) {
    errors.push_back(label + " has no atoms to merge");
    return errors;
  }
  if (molecule.natoms() == 0) {
    return errors;
  }
  if (molecule.number_of_frames() == 0) {
    errors.push_back(label + " has no coordinate frames");
    return errors;
  }
  const auto expected_coor = molecule.natoms() * molecule.number_of_frames() * 3;
  if (molecule.coor().size() != expected_coor) {
    errors.push_back(label + " coordinate storage does not match shape");
    return errors;
  }

  validate_required_vector(molecule, molecule.record(), label + " record",
                           errors);
  validate_required_vector(molecule, molecule.index(), label + " index",
                           errors);
  validate_required_vector(molecule, molecule.original_index(),
                           label + " original_index", errors);
  validate_required_vector(molecule, molecule.original_resid(),
                           label + " original_resid", errors);
  validate_required_vector(molecule, molecule.name(), label + " name", errors);
  validate_required_vector(molecule, molecule.loc(), label + " loc", errors);
  validate_required_vector(molecule, molecule.resname(), label + " resname",
                           errors);
  validate_required_vector(molecule, molecule.chain(), label + " chain", errors);
  validate_required_vector(molecule, molecule.resid(), label + " resid", errors);
  validate_required_vector(molecule, molecule.rescode(), label + " rescode",
                           errors);
  validate_required_vector(molecule, molecule.occupancy(),
                           label + " occupancy", errors);
  validate_required_vector(molecule, molecule.beta(), label + " beta", errors);
  validate_required_vector(molecule, molecule.segname(), label + " segname",
                           errors);
  validate_required_vector(molecule, molecule.element(), label + " element",
                           errors);
  validate_required_vector(molecule, molecule.charge(), label + " charge",
                           errors);
  validate_required_vector(molecule, molecule.moltype(), label + " moltype",
                           errors);
  validate_required_vector(molecule, molecule.mass(), label + " mass", errors);
  if (!molecule.atom_charge().empty()) {
    validate_required_vector(molecule, molecule.atom_charge(),
                             label + " atom_charge", errors);
  }
  if (!molecule.atom_vdw().empty()) {
    validate_required_vector(molecule, molecule.atom_vdw(),
                             label + " atom_vdw", errors);
  }
  validate_required_vector(molecule, molecule.residue_flag(),
                           label + " residue_flag", errors);
  if (!molecule.charmm_type().empty()) {
    validate_required_vector(molecule, molecule.charmm_type(),
                             label + " charmm_type", errors);
  }
  if (!molecule.residue_charge().empty()) {
    validate_required_vector(molecule, molecule.residue_charge(),
                             label + " residue_charge", errors);
  }
  validate_descriptor_map(molecule, molecule.extra_string_descriptors(),
                          label + " extra_string_descriptors", errors);
  validate_descriptor_map(molecule, molecule.extra_int_descriptors(),
                          label + " extra_int_descriptors", errors);
  validate_descriptor_map(molecule, molecule.extra_calc_descriptors(),
                          label + " extra_calc_descriptors", errors);
  return errors;
}

}  // namespace

IndexSelection get_indices_from_mask(const Molecule& molecule,
                                     const std::vector<int>& mask) {
  IndexSelection result;
  if (mask.size() != molecule.natoms()) {
    result.errors.push_back("mask length does not match natoms");
    return result;
  }
  for (std::size_t atom = 0; atom < mask.size(); ++atom) {
    if (mask[atom] == 1) {
      result.indices.push_back(atom);
    } else if (mask[atom] != 0) {
      result.indices.clear();
      result.errors.push_back("mask values must be 0 or 1");
      return result;
    }
  }
  if (result.indices.empty()) {
    result.errors.push_back("mask selects no atoms");
  }
  return result;
}

MaskMatrixSelection get_dihedral_subset_mask(
    const Molecule& molecule, const std::vector<int>& flexible_residues,
    int molecule_type) {
  MaskMatrixSelection result;
  if (molecule.name().size() != molecule.natoms() ||
      molecule.resid().size() != molecule.natoms()) {
    result.errors.push_back("name and resid descriptors must match natoms");
    return result;
  }
  if (molecule_type != 0 && molecule_type != 1) {
    result.errors.push_back("dihedral molecule type must be 0 protein or 1 RNA");
    return result;
  }

  result.masks.assign(flexible_residues.size(),
                      std::vector<int>(molecule.natoms(), 0));
  for (std::size_t row = 0; row < flexible_residues.size(); ++row) {
    const auto q0 = flexible_residues[row];
    for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
      const auto resid = molecule.resid()[atom];
      const auto& name = molecule.name()[atom];
      bool include = false;
      if (molecule_type == 0) {
        include = (resid == q0 - 1 && name == "C") ||
                  (resid == q0 && (name == "N" || name == "CA" ||
                                   name == "C")) ||
                  (resid == q0 + 1 && name == "N");
      } else {
        include = (resid == q0 - 1 && name == "O3'") ||
                  (resid == q0 &&
                   (name == "P" || name == "O5'" || name == "C5'" ||
                    name == "C4'" || name == "C3'" || name == "O3'")) ||
                  (resid == q0 + 1 && (name == "P" || name == "O5'"));
      }
      result.masks[row][atom] = include ? 1 : 0;
    }
  }
  return result;
}

CoordinateSelection get_coordinates_using_indices(
    const Molecule& molecule, std::size_t frame,
    const std::vector<std::size_t>& indices) {
  CoordinateSelection result;
  result.errors = validate_indices(molecule, indices, frame);
  if (!result.ok()) {
    return result;
  }

  result.coordinates.reserve(indices.size());
  for (const auto atom : indices) {
    result.coordinates.push_back(molecule.coordinate(frame, atom));
  }
  return result;
}

CoordinateSelection get_coordinates_using_mask(
    const Molecule& molecule, std::size_t frame, const std::vector<int>& mask) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return get_coordinates_using_indices(molecule, frame, selected.indices);
}

SubsetResult set_coordinates_using_indices(
    Molecule& molecule, const Molecule& source, std::size_t frame,
    const std::vector<std::size_t>& indices) {
  SubsetResult result;
  result.errors = validate_indices(molecule, indices, frame);
  if (!result.ok()) {
    return result;
  }
  if (source.number_of_frames() <= frame) {
    result.errors.push_back("source frame is out of range");
    return result;
  }
  if (source.natoms() != indices.size()) {
    result.errors.push_back(
        "source atom count must match selected destination atom count");
    return result;
  }

  for (std::size_t selected = 0; selected < indices.size(); ++selected) {
    molecule.set_coordinate(frame, indices[selected],
                            source.coordinate(frame, selected));
  }
  return result;
}

SubsetResult set_coordinates_using_mask(
    Molecule& molecule, const Molecule& source, std::size_t frame,
    const std::vector<int>& mask) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {selected.errors};
  }
  return set_coordinates_using_indices(molecule, source, frame, selected.indices);
}

Molecule with_coordinates_using_indices(
    const Molecule& target, const Molecule& source, std::size_t frame,
    const std::vector<std::size_t>& indices) {
  Molecule copy = target;
  const auto result = set_coordinates_using_indices(copy, source, frame, indices);
  if (!result.ok()) {
    throw std::invalid_argument(result.errors.front());
  }
  return copy;
}

Molecule with_coordinates_using_mask(const Molecule& target,
                                     const Molecule& source, std::size_t frame,
                                     const std::vector<int>& mask) {
  Molecule copy = target;
  const auto result = set_coordinates_using_mask(copy, source, frame, mask);
  if (!result.ok()) {
    throw std::invalid_argument(result.errors.front());
  }
  return copy;
}

SubsetResult copy_molecule_using_indices(
    const Molecule& source, Molecule& destination,
    const std::vector<std::size_t>& indices, std::size_t frame) {
  SubsetResult result;
  result.errors = validate_indices(source, indices, frame);
  if (!result.ok()) {
    return result;
  }

  try {
    destination.resize(indices.size(), 1);
    copy_selected_vector(source.record(), destination.record(), indices);
    copy_selected_vector(source.index(), destination.index(), indices);
    copy_selected_vector(source.original_index(), destination.original_index(),
                         indices);
    copy_selected_vector(source.original_resid(), destination.original_resid(),
                         indices);
    copy_selected_vector(source.name(), destination.name(), indices);
    copy_selected_vector(source.loc(), destination.loc(), indices);
    copy_selected_vector(source.resname(), destination.resname(), indices);
    copy_selected_vector(source.chain(), destination.chain(), indices);
    copy_selected_vector(source.resid(), destination.resid(), indices);
    copy_selected_vector(source.rescode(), destination.rescode(), indices);
    copy_selected_vector(source.occupancy(), destination.occupancy(), indices);
    copy_selected_vector(source.beta(), destination.beta(), indices);
    copy_selected_vector(source.segname(), destination.segname(), indices);
    copy_selected_vector(source.element(), destination.element(), indices);
    copy_selected_vector(source.charge(), destination.charge(), indices);
    copy_selected_vector(source.atom_charge(), destination.atom_charge(), indices);
    copy_selected_vector(source.atom_vdw(), destination.atom_vdw(), indices);
    copy_selected_vector(source.residue_flag(), destination.residue_flag(),
                         indices);
    if (!source.charmm_type().empty()) {
      copy_selected_vector(source.charmm_type(), destination.charmm_type(),
                           indices);
    } else {
      destination.charmm_type().clear();
    }
    copy_selected_vector(source.moltype(), destination.moltype(), indices);
    copy_selected_vector(source.mass(), destination.mass(), indices);
    copy_selected_descriptor_map(source.extra_string_descriptors(),
                                 destination.extra_string_descriptors(),
                                 indices);
    copy_selected_descriptor_map(source.extra_int_descriptors(),
                                 destination.extra_int_descriptors(), indices);
    copy_selected_descriptor_map(source.extra_calc_descriptors(),
                                 destination.extra_calc_descriptors(), indices);

    for (std::size_t selected = 0; selected < indices.size(); ++selected) {
      destination.set_coordinate(0, selected,
                                 source.coordinate(frame, indices[selected]));
    }
    copy_selected_conect(source, destination, indices);
  } catch (const std::exception& error) {
    result.errors.push_back(error.what());
  }

  return result;
}

SubsetResult copy_molecule_using_mask(const Molecule& source,
                                      Molecule& destination,
                                      const std::vector<int>& mask,
                                      std::size_t frame) {
  auto selected = get_indices_from_mask(source, mask);
  if (!selected.ok()) {
    return {selected.errors};
  }
  return copy_molecule_using_indices(source, destination, selected.indices, frame);
}

Molecule copied_molecule_using_indices(const Molecule& source,
                                       const std::vector<std::size_t>& indices,
                                       std::size_t frame) {
  Molecule destination;
  const auto result = copy_molecule_using_indices(source, destination, indices,
                                                  frame);
  if (!result.ok()) {
    throw std::invalid_argument(result.errors.front());
  }
  return destination;
}

Molecule copied_molecule_using_mask(const Molecule& source,
                                    const std::vector<int>& mask,
                                    std::size_t frame) {
  Molecule destination;
  const auto result = copy_molecule_using_mask(source, destination, mask, frame);
  if (!result.ok()) {
    throw std::invalid_argument(result.errors.front());
  }
  return destination;
}

std::vector<Molecule> duplicate_molecule(const Molecule& molecule,
                                         std::size_t number_of_duplicates) {
  std::vector<Molecule> duplicates;
  duplicates.reserve(number_of_duplicates);
  for (std::size_t i = 0; i < number_of_duplicates; ++i) {
    duplicates.push_back(molecule);
  }
  return duplicates;
}

SubsetResult merge_two_molecules(const Molecule& mol1, const Molecule& mol2,
                                 Molecule& merged, MergeOptions options) {
  SubsetResult result;
  auto errors = validate_merge_source(mol1, "mol1", true);
  append_vector(result.errors, errors);
  errors = validate_merge_source(mol2, "mol2", false);
  append_vector(result.errors, errors);
  if (!result.ok()) {
    return result;
  }

  const auto natoms1 = mol1.natoms();
  const auto natoms2 = mol2.natoms();
  const auto natoms = natoms1 + natoms2;

  Molecule candidate(natoms, 1);

  merge_required_vector(mol1.record(), mol2.record(), candidate.record());
  merge_required_vector(mol1.original_index(), mol2.original_index(),
                        candidate.original_index());
  merge_required_vector(mol1.original_resid(), mol2.original_resid(),
                        candidate.original_resid());
  merge_required_vector(mol1.name(), mol2.name(), candidate.name());
  merge_required_vector(mol1.loc(), mol2.loc(), candidate.loc());
  merge_required_vector(mol1.resname(), mol2.resname(), candidate.resname());
  merge_required_vector(mol1.chain(), mol2.chain(), candidate.chain());
  merge_required_vector(mol1.resid(), mol2.resid(), candidate.resid());
  merge_required_vector(mol1.rescode(), mol2.rescode(), candidate.rescode());
  merge_required_vector(mol1.occupancy(), mol2.occupancy(),
                        candidate.occupancy());
  merge_required_vector(mol1.beta(), mol2.beta(), candidate.beta());
  merge_required_vector(mol1.segname(), mol2.segname(), candidate.segname());
  merge_required_vector(mol1.element(), mol2.element(), candidate.element());
  merge_required_vector(mol1.charge(), mol2.charge(), candidate.charge());
  merge_required_vector(mol1.residue_flag(), mol2.residue_flag(),
                        candidate.residue_flag());
  merge_required_vector(mol1.moltype(), mol2.moltype(), candidate.moltype());
  merge_required_vector(mol1.mass(), mol2.mass(), candidate.mass());

  candidate.index().clear();
  candidate.index().reserve(natoms);
  append_vector(candidate.index(), mol1.index());
  const auto last_index = mol1.index().back();
  for (std::size_t atom = 0; atom < natoms2; ++atom) {
    candidate.index().push_back(last_index + static_cast<int>(atom) + 1);
  }

  merge_optional_vector(mol1.atom_charge(), natoms1, mol2.atom_charge(), natoms2,
                        candidate.atom_charge(), "atom_charge", options,
                        result);
  merge_optional_vector(mol1.atom_vdw(), natoms1, mol2.atom_vdw(), natoms2,
                        candidate.atom_vdw(), "atom_vdw", options, result);
  merge_optional_vector(mol1.charmm_type(), natoms1, mol2.charmm_type(), natoms2,
                        candidate.charmm_type(), "charmm_type", options,
                        result);
  merge_optional_vector(mol1.residue_charge(), natoms1, mol2.residue_charge(),
                        natoms2,
                        candidate.residue_charge(), "residue_charge", options,
                        result);
  merge_descriptor_map(mol1.extra_string_descriptors(), natoms1,
                       mol2.extra_string_descriptors(), natoms2,
                       candidate.extra_string_descriptors(),
                       "extra_string_descriptors", options, result);
  merge_descriptor_map(mol1.extra_int_descriptors(), natoms1,
                       mol2.extra_int_descriptors(), natoms2,
                       candidate.extra_int_descriptors(),
                       "extra_int_descriptors", options, result);
  merge_descriptor_map(mol1.extra_calc_descriptors(), natoms1,
                       mol2.extra_calc_descriptors(), natoms2,
                       candidate.extra_calc_descriptors(),
                       "extra_calc_descriptors", options, result);

  for (std::size_t atom = 0; atom < natoms1; ++atom) {
    candidate.set_coordinate(0, atom, mol1.coordinate(0, atom));
  }
  for (std::size_t atom = 0; atom < natoms2; ++atom) {
    candidate.set_coordinate(0, natoms1 + atom, mol2.coordinate(0, atom));
  }

  candidate.conect().assign(natoms, {});
  for (std::size_t atom = 0; atom < natoms1 && atom < mol1.conect().size();
       ++atom) {
    candidate.conect()[atom] = mol1.conect()[atom];
  }
  for (std::size_t atom = 0; atom < natoms2 && atom < mol2.conect().size();
       ++atom) {
    candidate.conect()[natoms1 + atom] = mol2.conect()[atom];
  }

  candidate.unitcell() = mol1.unitcell();
  if (!mol1.fasta().empty() || !mol2.fasta().empty()) {
    candidate.fasta() = mol1.fasta() + mol2.fasta();
  }
  if (mol1.total_mass() != calc_type{} && mol2.total_mass() != calc_type{}) {
    candidate.set_total_mass(mol1.total_mass() + mol2.total_mass());
  }

  merged = candidate;
  return result;
}

Molecule merged_two_molecules(const Molecule& mol1, const Molecule& mol2,
                              MergeOptions options) {
  Molecule merged;
  const auto result = merge_two_molecules(mol1, mol2, merged, options);
  if (!result.ok()) {
    throw std::invalid_argument(result.errors.front());
  }
  return merged;
}

StringSelection get_string_descriptor_using_indices(
    const Molecule& molecule, StringDescriptor descriptor,
    const std::vector<std::size_t>& indices) {
  StringSelection result;
  result.errors = validate_atom_indices(molecule, indices);
  if (!result.ok()) {
    return result;
  }
  const auto& values = string_descriptor(molecule, descriptor);
  result.errors = validate_descriptor_vector(molecule, values);
  if (!result.ok()) {
    return result;
  }
  result.values = get_descriptor_values(values, indices);
  return result;
}

StringSelection get_string_descriptor_using_mask(
    const Molecule& molecule, StringDescriptor descriptor,
    const std::vector<int>& mask) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return get_string_descriptor_using_indices(molecule, descriptor,
                                             selected.indices);
}

SubsetResult set_string_descriptor_using_indices(
    Molecule& molecule, StringDescriptor descriptor,
    const std::vector<std::size_t>& indices, const std::string& value) {
  return set_descriptor_values(string_descriptor(molecule, descriptor), molecule,
                               indices, value);
}

SubsetResult set_string_descriptor_using_mask(
    Molecule& molecule, StringDescriptor descriptor, const std::vector<int>& mask,
    const std::string& value) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {selected.errors};
  }
  return set_string_descriptor_using_indices(molecule, descriptor,
                                             selected.indices, value);
}

IntSelection get_int_descriptor_using_indices(
    const Molecule& molecule, IntDescriptor descriptor,
    const std::vector<std::size_t>& indices) {
  IntSelection result;
  result.errors = validate_atom_indices(molecule, indices);
  if (!result.ok()) {
    return result;
  }
  const auto& values = int_descriptor(molecule, descriptor);
  result.errors = validate_descriptor_vector(molecule, values);
  if (!result.ok()) {
    return result;
  }
  result.values = get_descriptor_values(values, indices);
  return result;
}

IntSelection get_int_descriptor_using_mask(const Molecule& molecule,
                                           IntDescriptor descriptor,
                                           const std::vector<int>& mask) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return get_int_descriptor_using_indices(molecule, descriptor, selected.indices);
}

SubsetResult set_int_descriptor_using_indices(
    Molecule& molecule, IntDescriptor descriptor,
    const std::vector<std::size_t>& indices, int value) {
  return set_descriptor_values(int_descriptor(molecule, descriptor), molecule,
                               indices, value);
}

SubsetResult set_int_descriptor_using_mask(Molecule& molecule,
                                           IntDescriptor descriptor,
                                           const std::vector<int>& mask,
                                           int value) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {selected.errors};
  }
  return set_int_descriptor_using_indices(molecule, descriptor, selected.indices,
                                          value);
}

CalcSelection get_calc_descriptor_using_indices(
    const Molecule& molecule, CalcDescriptor descriptor,
    const std::vector<std::size_t>& indices) {
  CalcSelection result;
  result.errors = validate_atom_indices(molecule, indices);
  if (!result.ok()) {
    return result;
  }
  const auto& values = calc_descriptor(molecule, descriptor);
  result.errors = validate_descriptor_vector(molecule, values);
  if (!result.ok()) {
    return result;
  }
  result.values = get_descriptor_values(values, indices);
  return result;
}

CalcSelection get_calc_descriptor_using_mask(const Molecule& molecule,
                                             CalcDescriptor descriptor,
                                             const std::vector<int>& mask) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return get_calc_descriptor_using_indices(molecule, descriptor, selected.indices);
}

SubsetResult set_calc_descriptor_using_indices(
    Molecule& molecule, CalcDescriptor descriptor,
    const std::vector<std::size_t>& indices, calc_type value) {
  return set_descriptor_values(calc_descriptor(molecule, descriptor), molecule,
                               indices, value);
}

SubsetResult set_calc_descriptor_using_mask(Molecule& molecule,
                                            CalcDescriptor descriptor,
                                            const std::vector<int>& mask,
                                            calc_type value) {
  auto selected = get_indices_from_mask(molecule, mask);
  if (!selected.ok()) {
    return {selected.errors};
  }
  return set_calc_descriptor_using_indices(molecule, descriptor, selected.indices,
                                           value);
}

}  // namespace sasmol
