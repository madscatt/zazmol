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

}  // namespace sasmol
