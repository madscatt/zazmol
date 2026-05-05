#pragma once

#include "sasmol/molecule.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace sasmol {

enum class StringDescriptor {
  record,
  name,
  loc,
  resname,
  chain,
  rescode,
  occupancy,
  beta,
  segname,
  element,
  charge,
  charmm_type,
  moltype,
};

enum class IntDescriptor {
  index,
  original_index,
  original_resid,
  resid,
  residue_flag,
};

enum class CalcDescriptor {
  atom_charge,
  atom_vdw,
  mass,
};

struct SubsetResult {
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CoordinateSelection {
  std::vector<Vec3> coordinates;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct IndexSelection {
  std::vector<std::size_t> indices;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct StringSelection {
  std::vector<std::string> values;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct IntSelection {
  std::vector<int> values;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct CalcSelection {
  std::vector<calc_type> values;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct MaskMatrixSelection {
  std::vector<std::vector<int>> masks;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct MergeOptions {
  bool report_skipped_descriptors{false};
};

struct BiomtTransform {
  std::array<std::array<calc_type, 3>, 3> rotation{};
  std::array<calc_type, 3> translation{};
};

[[nodiscard]] IndexSelection get_indices_from_mask(
    const Molecule& molecule, const std::vector<int>& mask);

[[nodiscard]] MaskMatrixSelection get_dihedral_subset_mask(
    const Molecule& molecule, const std::vector<int>& flexible_residues,
    int molecule_type);

[[nodiscard]] CoordinateSelection get_coordinates_using_indices(
    const Molecule& molecule, std::size_t frame,
    const std::vector<std::size_t>& indices);

[[nodiscard]] CoordinateSelection get_coordinates_using_mask(
    const Molecule& molecule, std::size_t frame, const std::vector<int>& mask);

[[nodiscard]] SubsetResult set_coordinates_using_indices(
    Molecule& molecule, const Molecule& source, std::size_t frame,
    const std::vector<std::size_t>& indices);

[[nodiscard]] SubsetResult set_coordinates_using_mask(
    Molecule& molecule, const Molecule& source, std::size_t frame,
    const std::vector<int>& mask);

[[nodiscard]] Molecule with_coordinates_using_indices(
    const Molecule& target, const Molecule& source, std::size_t frame,
    const std::vector<std::size_t>& indices);

[[nodiscard]] Molecule with_coordinates_using_mask(
    const Molecule& target, const Molecule& source, std::size_t frame,
    const std::vector<int>& mask);

[[nodiscard]] SubsetResult copy_molecule_using_indices(
    const Molecule& source, Molecule& destination,
    const std::vector<std::size_t>& indices, std::size_t frame);

[[nodiscard]] SubsetResult copy_molecule_using_mask(
    const Molecule& source, Molecule& destination, const std::vector<int>& mask,
    std::size_t frame);

[[nodiscard]] Molecule copied_molecule_using_indices(
    const Molecule& source, const std::vector<std::size_t>& indices,
    std::size_t frame);

[[nodiscard]] Molecule copied_molecule_using_mask(
    const Molecule& source, const std::vector<int>& mask, std::size_t frame);

[[nodiscard]] std::vector<Molecule> duplicate_molecule(
    const Molecule& molecule, std::size_t number_of_duplicates);

[[nodiscard]] SubsetResult merge_two_molecules(
    const Molecule& mol1, const Molecule& mol2, Molecule& merged,
    MergeOptions options = {});

[[nodiscard]] Molecule merged_two_molecules(
    const Molecule& mol1, const Molecule& mol2, MergeOptions options = {});

[[nodiscard]] SubsetResult apply_biomt_transforms(
    const Molecule& source, std::size_t frame,
    const std::vector<BiomtTransform>& transforms, Molecule& transformed);

[[nodiscard]] Molecule biomt_transformed(
    const Molecule& source, std::size_t frame,
    const std::vector<BiomtTransform>& transforms);

[[nodiscard]] SubsetResult apply_biomt_transforms_from_metadata(
    const Molecule& source, std::size_t frame, int biomol_id,
    Molecule& transformed);

[[nodiscard]] Molecule biomt_transformed_from_metadata(
    const Molecule& source, std::size_t frame, int biomol_id);

[[nodiscard]] StringSelection get_string_descriptor_using_indices(
    const Molecule& molecule, StringDescriptor descriptor,
    const std::vector<std::size_t>& indices);
[[nodiscard]] StringSelection get_string_descriptor_using_mask(
    const Molecule& molecule, StringDescriptor descriptor,
    const std::vector<int>& mask);
[[nodiscard]] SubsetResult set_string_descriptor_using_indices(
    Molecule& molecule, StringDescriptor descriptor,
    const std::vector<std::size_t>& indices, const std::string& value);
[[nodiscard]] SubsetResult set_string_descriptor_using_mask(
    Molecule& molecule, StringDescriptor descriptor, const std::vector<int>& mask,
    const std::string& value);

[[nodiscard]] IntSelection get_int_descriptor_using_indices(
    const Molecule& molecule, IntDescriptor descriptor,
    const std::vector<std::size_t>& indices);
[[nodiscard]] IntSelection get_int_descriptor_using_mask(
    const Molecule& molecule, IntDescriptor descriptor,
    const std::vector<int>& mask);
[[nodiscard]] SubsetResult set_int_descriptor_using_indices(
    Molecule& molecule, IntDescriptor descriptor,
    const std::vector<std::size_t>& indices, int value);
[[nodiscard]] SubsetResult set_int_descriptor_using_mask(
    Molecule& molecule, IntDescriptor descriptor, const std::vector<int>& mask,
    int value);

[[nodiscard]] CalcSelection get_calc_descriptor_using_indices(
    const Molecule& molecule, CalcDescriptor descriptor,
    const std::vector<std::size_t>& indices);
[[nodiscard]] CalcSelection get_calc_descriptor_using_mask(
    const Molecule& molecule, CalcDescriptor descriptor,
    const std::vector<int>& mask);
[[nodiscard]] SubsetResult set_calc_descriptor_using_indices(
    Molecule& molecule, CalcDescriptor descriptor,
    const std::vector<std::size_t>& indices, calc_type value);
[[nodiscard]] SubsetResult set_calc_descriptor_using_mask(
    Molecule& molecule, CalcDescriptor descriptor, const std::vector<int>& mask,
    calc_type value);

}  // namespace sasmol
