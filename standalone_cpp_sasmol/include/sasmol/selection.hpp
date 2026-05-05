#pragma once

#include "sasmol/molecule.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace sasmol {

struct SelectionResult {
  std::vector<std::size_t> indices;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct MaskSelectionResult {
  std::vector<int> mask;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

enum class SassieBasisContext {
  generic,
  protein,
  nucleic,
  nucleic_overlap,
};

struct BasisExpressionResult {
  std::string expression;
  std::vector<std::string> errors;

  [[nodiscard]] bool ok() const noexcept { return errors.empty(); }
};

struct SegmentBasisExpressionResult {
  std::vector<std::string> expressions;
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

[[nodiscard]] std::optional<std::string> basis_expression(
    const std::string& basis_name);
[[nodiscard]] BasisExpressionResult sassie_basis_expression(
    const std::string& basis, SassieBasisContext context = SassieBasisContext::generic);
[[nodiscard]] SegmentBasisExpressionResult sassie_segment_basis_expressions(
    const std::string& basis, const std::vector<std::string>& segnames,
    const std::vector<SassieBasisContext>& contexts = {});
[[nodiscard]] SelectionResult select_named_basis(const Molecule& molecule,
                                                 const std::string& basis_name);
[[nodiscard]] SelectionResult select_sassie_basis(
    const Molecule& molecule, const std::string& basis,
    SassieBasisContext context = SassieBasisContext::generic);
[[nodiscard]] SelectionResult select_indices(const Molecule& molecule,
                                             const std::string& expression);
[[nodiscard]] MaskSelectionResult mask_from_indices(
    const Molecule& molecule, const std::vector<std::size_t>& indices);
[[nodiscard]] MaskSelectionResult select_mask(const Molecule& molecule,
                                              const std::string& expression);
[[nodiscard]] MaskSelectionResult select_named_basis_mask(
    const Molecule& molecule, const std::string& basis_name);
[[nodiscard]] MaskSelectionResult select_sassie_basis_mask(
    const Molecule& molecule, const std::string& basis,
    SassieBasisContext context = SassieBasisContext::generic);

}  // namespace sasmol
