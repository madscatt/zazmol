#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace sasmol {

struct SelectionResult {
  std::vector<std::size_t> indices;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

[[nodiscard]] SelectionResult indices_all(const Molecule& molecule);
[[nodiscard]] SelectionResult indices_by_name(const Molecule& molecule,
                                              const std::string& name);
[[nodiscard]] SelectionResult indices_by_resname(const Molecule& molecule,
                                                 const std::string& resname);
[[nodiscard]] SelectionResult indices_by_resid_range(const Molecule& molecule,
                                                     int first_resid,
                                                     int last_resid);

[[nodiscard]] SelectionResult select_indices(const Molecule& molecule,
                                             const std::string& expression);

}  // namespace sasmol
