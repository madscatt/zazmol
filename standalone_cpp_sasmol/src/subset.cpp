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
    copy_selected_vector(source.moltype(), destination.moltype(), indices);
    copy_selected_vector(source.mass(), destination.mass(), indices);

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

std::vector<Molecule> duplicate_molecule(const Molecule& molecule,
                                         std::size_t number_of_duplicates) {
  std::vector<Molecule> duplicates;
  duplicates.reserve(number_of_duplicates);
  for (std::size_t i = 0; i < number_of_duplicates; ++i) {
    duplicates.push_back(molecule);
  }
  return duplicates;
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
