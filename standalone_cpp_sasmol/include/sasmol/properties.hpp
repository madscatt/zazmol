#pragma once

#include "sasmol/types.hpp"

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace sasmol {

using PropertyTable = std::map<std::string, std::vector<calc_type>>;
using OptionalPropertyTable = std::map<std::string, std::optional<calc_type>>;

[[nodiscard]] const std::map<std::string, calc_type>& amu();
[[nodiscard]] const PropertyTable& amino_acid_sld();
[[nodiscard]] const PropertyTable& element_scattering_lengths();
[[nodiscard]] const PropertyTable& element_sl();
[[nodiscard]] const PropertyTable& nucleotide_scattering_lengths();
[[nodiscard]] const PropertyTable& nucleotide_sl();
[[nodiscard]] const PropertyTable& dna_scattering_lengths();
[[nodiscard]] const PropertyTable& dna_sl();
[[nodiscard]] const PropertyTable& rna_scattering_lengths();
[[nodiscard]] const PropertyTable& rna_sl();
[[nodiscard]] const PropertyTable& protein_scattering_lengths();
[[nodiscard]] const PropertyTable& protein_sl();
[[nodiscard]] const OptionalPropertyTable& van_der_waals_radii();
[[nodiscard]] const OptionalPropertyTable& van_der_Waals_radii();

}  // namespace sasmol
