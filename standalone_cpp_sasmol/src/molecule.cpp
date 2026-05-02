#include "sasmol/molecule.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <tuple>

namespace sasmol {

namespace {

void require_index(std::size_t value, std::size_t limit,
                   const char* label) {
  if (value >= limit) {
    throw std::out_of_range(std::string(label) + " is out of range");
  }
}

template <typename T>
void resize_descriptor(std::vector<T>& values, std::size_t natoms,
                       const T& default_value) {
  values.assign(natoms, default_value);
}

template <typename T>
T value_or(const std::vector<T>& values, std::size_t index,
           const T& default_value) {
  if (index < values.size()) {
    return values[index];
  }
  return default_value;
}

void append_unique(std::vector<std::string>& values,
                   const std::string& value) {
  if (std::find(values.begin(), values.end(), value) == values.end()) {
    values.push_back(value);
  }
}

bool contains(const std::set<std::string>& values, const std::string& value) {
  return values.find(value) != values.end();
}

std::string join(const std::vector<std::string>& values,
                 const std::string& separator) {
  std::string result;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      result += separator;
    }
    result += values[i];
  }
  return result;
}

}  // namespace

Vec3 ConstCoordinateView::operator[](std::size_t atom) const {
  require_index(atom, natoms, "atom");
  const std::size_t base = atom * 3;
  return {xyz[base], xyz[base + 1], xyz[base + 2]};
}

Vec3 CoordinateView::operator[](std::size_t atom) const {
  require_index(atom, natoms, "atom");
  const std::size_t base = atom * 3;
  return {xyz[base], xyz[base + 1], xyz[base + 2]};
}

void CoordinateView::set(std::size_t atom, Vec3 value) {
  require_index(atom, natoms, "atom");
  const std::size_t base = atom * 3;
  xyz[base] = value.x;
  xyz[base + 1] = value.y;
  xyz[base + 2] = value.z;
}

Molecule::Molecule(std::size_t natoms, std::size_t number_of_frames) {
  resize(natoms, number_of_frames);
}

void Molecule::resize(std::size_t natoms, std::size_t number_of_frames) {
  natoms_ = natoms;
  number_of_frames_ = number_of_frames;
  coor_.assign(natoms_ * number_of_frames_ * 3, coord_type{});

  resize_descriptor(record_, natoms_, std::string{"ATOM"});
  index_.resize(natoms_);
  original_index_.resize(natoms_);
  original_resid_.assign(natoms_, 1);
  resize_descriptor(name_, natoms_, std::string{});
  resize_descriptor(loc_, natoms_, std::string{});
  resize_descriptor(resname_, natoms_, std::string{});
  resize_descriptor(chain_, natoms_, std::string{});
  resid_.assign(natoms_, 1);
  resize_descriptor(rescode_, natoms_, std::string{});
  resize_descriptor(occupancy_, natoms_, std::string{"1.00"});
  resize_descriptor(beta_, natoms_, std::string{"0.00"});
  resize_descriptor(segname_, natoms_, std::string{});
  resize_descriptor(element_, natoms_, std::string{});
  resize_descriptor(charge_, natoms_, std::string{});
  resize_descriptor(atom_charge_, natoms_, calc_type{});
  resize_descriptor(atom_vdw_, natoms_, calc_type{});
  resize_descriptor(residue_flag_, natoms_, 0);
  charmm_type_.clear();
  resize_descriptor(moltype_, natoms_, std::string{});
  resize_descriptor(mass_, natoms_, calc_type{});
  resize_descriptor(residue_charge_, natoms_, calc_type{});
  conect_.assign(natoms_, {});

  for (std::size_t i = 0; i < natoms_; ++i) {
    const int one_based = static_cast<int>(i + 1);
    index_[i] = one_based;
    original_index_[i] = one_based;
  }

  total_mass_ = calc_type{};
  formula_.clear();
  fasta_.clear();
  unitcell_.fill(calc_type{});
  extra_string_descriptors_.clear();
  extra_int_descriptors_.clear();
  extra_calc_descriptors_.clear();
}

std::size_t Molecule::offset(std::size_t frame, std::size_t atom) const {
  require_index(frame, number_of_frames_, "frame");
  require_index(atom, natoms_, "atom");
  return ((frame * natoms_) + atom) * 3;
}

Vec3 Molecule::coordinate(std::size_t frame, std::size_t atom) const {
  const std::size_t base = offset(frame, atom);
  return {coor_[base], coor_[base + 1], coor_[base + 2]};
}

void Molecule::set_coordinate(std::size_t frame, std::size_t atom,
                              Vec3 value) {
  const std::size_t base = offset(frame, atom);
  coor_[base] = value.x;
  coor_[base + 1] = value.y;
  coor_[base + 2] = value.z;
}

CoordinateView Molecule::coordinate_view(std::size_t frame) {
  require_index(frame, number_of_frames_, "frame");
  const std::size_t base = frame * natoms_ * 3;
  return {std::span<coord_type>{coor_}.subspan(base, natoms_ * 3), natoms_};
}

ConstCoordinateView Molecule::coordinate_view(std::size_t frame) const {
  require_index(frame, number_of_frames_, "frame");
  const std::size_t base = frame * natoms_ * 3;
  return {std::span<const coord_type>{coor_}.subspan(base, natoms_ * 3),
          natoms_};
}

IntegrityReport Molecule::check_integrity(bool fast_check) const {
  IntegrityReport report;
  const std::size_t expected_coor = natoms_ * number_of_frames_ * 3;

  report.lengths["coor"] = coor_.size();
  if (coor_.size() != expected_coor) {
    report.issues.push_back({"coor", expected_coor, coor_.size()});
    if (fast_check) {
      return report;
    }
  }

  auto add = [&](const std::string& field, const auto& values) {
    report.lengths[field] = values.size();
    if (values.size() != natoms_) {
      report.issues.push_back({field, natoms_, values.size()});
    }
  };

  add("record", record_);
  add("index", index_);
  add("original_index", original_index_);
  add("original_resid", original_resid_);
  add("name", name_);
  add("loc", loc_);
  add("resname", resname_);
  add("chain", chain_);
  add("resid", resid_);
  add("rescode", rescode_);
  add("occupancy", occupancy_);
  add("beta", beta_);
  add("segname", segname_);
  add("element", element_);
  add("charge", charge_);
  add("atom_charge", atom_charge_);
  add("atom_vdw", atom_vdw_);
  add("residue_flag", residue_flag_);
  if (!charmm_type_.empty()) {
    add("charmm_type", charmm_type_);
  }
  add("moltype", moltype_);
  add("mass", mass_);
  add("residue_charge", residue_charge_);
  add("conect", conect_);

  for (const auto& [name, values] : extra_string_descriptors_) {
    add("extra_string_descriptors." + name, values);
  }
  for (const auto& [name, values] : extra_int_descriptors_) {
    add("extra_int_descriptors." + name, values);
  }
  for (const auto& [name, values] : extra_calc_descriptors_) {
    add("extra_calc_descriptors." + name, values);
  }

  if (fast_check && !report.issues.empty()) {
    report.issues.resize(1);
  }

  return report;
}

MoltypeReport Molecule::moltype_by_segname_report() const {
  static const std::set<std::string> dna = {
      "NUSA", "NUSG", "NUSC", "NUSU", "DA", "DG", "DC", "DT",
      "ADE",  "GUA",  "CYT",  "THY"};
  static const std::set<std::string> rna = {
      "RNUS", "RNUA", "RUUG", "RNUC", "A", "C", "G", "U",
      "ADE",  "CYT",  "GUA",  "URA"};
  static const std::set<std::string> rna_atom_names = {"O2'", "O2*"};

  const std::string empty;
  MoltypeReport report;
  std::map<std::string, std::set<std::tuple<int, std::string, std::string>>>
      residues_by_segment;

  for (std::size_t atom = 0; atom < natoms_; ++atom) {
    const auto& seg = value_or(segname_, atom, empty);
    const auto& resn = value_or(resname_, atom, empty);
    const auto& atom_name = value_or(name_, atom, empty);
    const auto& assigned_moltype = value_or(moltype_, atom, empty);
    const int residue_id = value_or(resid_, atom, 0);
    const auto& insertion_code = value_or(rescode_, atom, empty);

    auto& segment = report.segments[seg];
    if (segment.segname.empty()) {
      segment.segname = seg;
    }
    ++segment.atom_count;
    append_unique(segment.assigned_moltypes, assigned_moltype);
    append_unique(segment.resnames, resn);
    residues_by_segment[seg].insert({residue_id, insertion_code, resn});

    if (contains(dna, resn) && contains(rna, resn)) {
      append_unique(segment.ambiguous_resnames, resn);
    }
    if (contains(dna, resn) && !contains(rna, resn)) {
      append_unique(segment.dna_resname_evidence, resn);
    }
    if (contains(rna, resn) && !contains(dna, resn)) {
      append_unique(segment.rna_resname_evidence, resn);
    }
    if (contains(rna_atom_names, atom_name)) {
      append_unique(segment.rna_atom_evidence, atom_name);
    }
  }

  bool has_mixed = false;
  bool has_ambiguous = false;
  bool has_all_other = false;

  for (auto& [seg, segment] : report.segments) {
    segment.residue_count = residues_by_segment[seg].size();
    std::vector<std::string> assigned;
    for (const auto& moltype_value : segment.assigned_moltypes) {
      if (!moltype_value.empty()) {
        assigned.push_back(moltype_value);
      }
    }

    const bool has_ambiguous_nucleic = !segment.ambiguous_resnames.empty();
    const bool has_specific_nucleic_evidence =
        !segment.dna_resname_evidence.empty() ||
        !segment.rna_resname_evidence.empty() ||
        !segment.rna_atom_evidence.empty();

    if (assigned.size() > 1) {
      segment.status = "mixed";
      has_mixed = true;
      const auto moltype_list = join(assigned, ", ");
      segment.evidence.push_back("multiple assigned moltypes: " + moltype_list);
      report.messages.push_back(
          "Segment " + seg +
          " contains multiple assigned moltypes: " + moltype_list + ".");
    } else if (assigned.size() == 1 && assigned.front() == "other") {
      segment.status = "all_other";
      has_all_other = true;
      segment.evidence.push_back("all atoms are assigned moltype other");
      report.messages.push_back(
          "Segment " + seg + " contains only moltype other assignments.");
    } else if (has_ambiguous_nucleic && !has_specific_nucleic_evidence) {
      segment.status = "ambiguous_nucleic";
      has_ambiguous = true;
      const auto resname_list = join(segment.ambiguous_resnames, ", ");
      segment.evidence.push_back(
          "DNA/RNA-overlap residue names without DNA- or RNA-specific "
          "evidence: " +
          resname_list);
      report.messages.push_back(
          "Segment " + seg +
          " contains DNA/RNA-overlap residue names without DNA- or "
          "RNA-specific evidence: " +
          resname_list + ".");
    } else {
      segment.status = "clean";
      if (!segment.dna_resname_evidence.empty()) {
        segment.evidence.push_back("DNA-specific residue names: " +
                                   join(segment.dna_resname_evidence, ", "));
      }
      if (!segment.rna_resname_evidence.empty()) {
        segment.evidence.push_back("RNA-specific residue names: " +
                                   join(segment.rna_resname_evidence, ", "));
      }
      if (!segment.rna_atom_evidence.empty()) {
        segment.evidence.push_back("RNA atom-name evidence: " +
                                   join(segment.rna_atom_evidence, ", "));
      }
    }
  }

  if (has_mixed) {
    report.overall_status = "mixed_by_segname";
  } else if (has_ambiguous) {
    report.overall_status = "ambiguous_nucleic";
  } else if (has_all_other) {
    report.overall_status = "unknown";
  }

  return report;
}

}  // namespace sasmol
