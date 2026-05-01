#pragma once

#include "sasmol/molecule.hpp"

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
  moltype,
};

enum class IntDescriptor {
  index,
  original_index,
  original_resid,
  resid,
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

[[nodiscard]] IndexSelection get_indices_from_mask(
    const Molecule& molecule, const std::vector<int>& mask);

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

[[nodiscard]] SubsetResult copy_molecule_using_indices(
    const Molecule& source, Molecule& destination,
    const std::vector<std::size_t>& indices, std::size_t frame);

[[nodiscard]] SubsetResult copy_molecule_using_mask(
    const Molecule& source, Molecule& destination, const std::vector<int>& mask,
    std::size_t frame);

[[nodiscard]] Molecule copied_molecule_using_indices(
    const Molecule& source, const std::vector<std::size_t>& indices,
    std::size_t frame);

[[nodiscard]] std::vector<Molecule> duplicate_molecule(
    const Molecule& molecule, std::size_t number_of_duplicates);

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
