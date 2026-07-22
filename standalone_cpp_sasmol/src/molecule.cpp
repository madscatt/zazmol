#include "sasmol/molecule.hpp"

#include <algorithm>
#include <cctype>
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

void append_unique(std::vector<std::string>& values,
                   const std::string& value) {
  if (std::find(values.begin(), values.end(), value) == values.end()) {
    values.push_back(value);
  }
}

bool contains(const std::set<std::string>& values, const std::string& value) {
  return values.find(value) != values.end();
}

std::string trim_copy(std::string value) {
  auto not_space = [](unsigned char c) { return !std::isspace(c); };
  value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
  value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(),
              value.end());
  return value;
}

struct MoltypeGroupKey {
  std::string source;
  std::string value;

  [[nodiscard]] std::string key() const {
    if (source == "chain") {
      return "chain:" + value;
    }
    return value;
  }
};

MoltypeGroupKey effective_moltype_group(const std::string& segname,
                                        const std::string& chain) {
  const auto trimmed_segname = trim_copy(segname);
  if (!trimmed_segname.empty()) {
    return {"segname", trimmed_segname};
  }

  const auto trimmed_chain = trim_copy(chain);
  if (!trimmed_chain.empty()) {
    return {"chain", trimmed_chain};
  }

  return {"none", ""};
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
  biomt_.clear();
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

MoltypeReport Molecule::classify_nucleic_moltypes() {
  // Coordinate-only counterpart of the Python residue-first classifier.  This
  // compact first pass deliberately ignores pre-existing moltype values.
  static const std::set<std::string> deoxy = {
      "DA", "DC", "DG", "DT", "DI", "DU", "NUSA", "NUSG", "NUSC", "NUSU"};
  static const std::set<std::string> ribo = {"RNUS", "RNUA", "RUUG", "RNUC"};
  static const std::set<std::string> ambiguous = {
      "A", "C", "G", "T", "U", "ADE", "CYT", "GUA", "THY", "URA"};
  static const std::set<std::string> o2 = {"O2'", "O2*"};
  static const std::array<std::set<std::string>, 5> ring = {{{"C1'", "C1*"},
      {"C2'", "C2*"}, {"C3'", "C3*"}, {"C4'", "C4*"}, {"O4'", "O4*"}}};

  struct Residue {
    std::vector<std::size_t> atoms;
    std::string seg;
    std::string chain;
    int resid{};
    std::string rescode;
    std::string name;
    std::set<std::string> names;
    std::set<std::string> records;
    std::vector<std::string> o2prime_atoms;
    std::vector<std::string> conflict_reasons;
    std::string identity;
    std::string evidence;
    bool complete_deoxy_sugar{};
  };
  std::vector<Residue> residues;
  std::string previous;
  for (std::size_t i = 0; i < natoms_; ++i) {
    const auto& rn = resname_[i];
    if (!contains(deoxy, rn) && !contains(ribo, rn) && !contains(ambiguous, rn)) {
      previous.clear(); continue;
    }
    const std::string signature = segname_[i] + "\x1f" + chain_[i] + "\x1f" +
        std::to_string(resid_[i]) + "\x1f" + rescode_[i] + "\x1f" + rn;
    if (signature != previous) {
      residues.push_back({{}, segname_[i], chain_[i], resid_[i], rescode_[i],
                          rn, {}, {}, {}, {}, "", "", false});
      previous = signature;
    }
    residues.back().atoms.push_back(i);
    residues.back().names.insert(trim_copy(name_[i]));
    const auto record_class = trim_copy(record_[i]);
    if (!record_class.empty()) {
      residues.back().records.insert(record_class);
    }
  }
  std::map<std::string, std::vector<Residue*>> groups;
  for (auto& residue : residues) {
    const bool has_o2 = std::any_of(o2.begin(), o2.end(), [&](const auto& a) { return residue.names.contains(a); });
    bool complete = !has_o2;
    for (const auto& aliases : ring) {
      bool found = false; for (const auto& alias : aliases) found = found || residue.names.contains(alias);
      complete = complete && found;
    }
    residue.complete_deoxy_sugar = complete;
    for (const auto& alias : o2) {
      if (residue.names.contains(alias)) {
        residue.o2prime_atoms.push_back(alias);
      }
    }
    if (contains(deoxy, residue.name) && has_o2) {
      residue.identity = "conflict";
      residue.evidence = "conflict";
      residue.conflict_reasons.push_back(
          "explicit deoxy name with O2prime evidence");
    } else if (contains(ribo, residue.name) && complete) {
      residue.identity = "conflict";
      residue.evidence = "conflict";
      residue.conflict_reasons.push_back(
          "explicit ribonucleotide name with complete deoxy sugar");
    } else if (has_o2) {
      residue.identity = "rna";
      residue.evidence = "o2prime";
    } else if (contains(deoxy, residue.name)) {
      residue.identity = "dna";
      residue.evidence = "explicit_deoxy_name";
    } else if (complete) {
      residue.identity = "dna";
      residue.evidence = "complete_deoxy_sugar";
    } else if (contains(ribo, residue.name)) {
      residue.identity = "rna";
      residue.evidence = "explicit_ribonucleotide_name";
    } else {
      residue.identity = "ambiguous";
      residue.evidence = "insufficient_coordinate_evidence";
    }
    const auto group = effective_moltype_group(residue.seg, residue.chain);
    const auto key = group.source == "none" ? "unassigned:" + std::to_string(1 + &residue - residues.data()) : group.key();
    groups[key].push_back(&residue);
  }
  MoltypeReport report;
  for (const auto& [key, group] : groups) {
    bool dna = false, rna = false, conflict = false;
    for (const auto* residue : group) { dna = dna || residue->identity == "dna"; rna = rna || residue->identity == "rna"; conflict = conflict || residue->identity == "conflict"; }
    const std::string canonical = (conflict || (dna && rna)) ? "nucleic" : dna ? "dna" : rna ? "rna" : "nucleic";
    for (const auto* residue : group) for (const auto atom : residue->atoms) moltype_[atom] = canonical;
    auto& segment = report.segments[key];
    const auto group_info = effective_moltype_group(group.front()->seg, group.front()->chain);
    segment.segname = group.front()->seg;
    segment.chain = group.front()->chain;
    segment.grouping_source = group_info.source;
    segment.grouping_value = group_info.value;
    segment.canonical_moltype = canonical;
    segment.identity = conflict || (dna && rna) ? "conflict" : dna ? "resolved_dna" : rna ? "resolved_rna" : "ambiguous";
    segment.status = segment.identity;
    segment.assigned_moltypes = {canonical};
    segment.decision = "assigned recognized nucleic atoms as " + canonical;
    if (segment.identity == "ambiguous") {
      segment.evidence.push_back("insufficient coordinate evidence to resolve DNA or RNA");
    }
    for (const auto* residue : group) {
      ++segment.residue_count;
      segment.atom_count += residue->atoms.size();
      append_unique(segment.chains, residue->chain);
      append_unique(segment.resnames, residue->name);
      if (contains(deoxy, residue->name)) append_unique(segment.dna_resname_evidence, residue->name);
      if (contains(ribo, residue->name)) append_unique(segment.rna_resname_evidence, residue->name);
      if (contains(ambiguous, residue->name) && residue->identity == "ambiguous") append_unique(segment.ambiguous_resnames, residue->name);
      for (const auto& atom_name : residue->o2prime_atoms) {
        append_unique(segment.rna_atom_evidence, atom_name);
      }
      if (residue->complete_deoxy_sugar) {
        append_unique(segment.complete_deoxy_sugar_evidence, residue->name);
      }
      for (const auto& record_class : residue->records) {
        append_unique(segment.record_classes, record_class);
      }
      MoltypeResidueDecision decision;
      decision.resname = residue->name;
      decision.resid = residue->resid;
      decision.rescode = residue->rescode;
      decision.chain = residue->chain;
      decision.identity = residue->identity;
      decision.evidence = residue->evidence;
      decision.o2prime_atoms = residue->o2prime_atoms;
      decision.complete_deoxy_sugar = residue->complete_deoxy_sugar;
      decision.conflict_reasons = residue->conflict_reasons;
      segment.residue_decisions.push_back(decision);
      if (residue->identity == "conflict") {
        segment.conflicting_residues.push_back(decision);
      }
    }
    if (segment.identity == "conflict") {
      report.overall_status = "conflict";
      report.messages.push_back("Segment " + key + ": conflicting DNA and RNA coordinate evidence.");
    } else if (segment.identity == "ambiguous" && report.overall_status != "conflict") {
      report.overall_status = "ambiguous";
      report.messages.push_back("Segment " + key + ": insufficient coordinate evidence to resolve DNA or RNA.");
    }
  }
  return report;
}

MoltypeReport Molecule::moltype_by_segname_report() const {
  Molecule classified = *this;
  return classified.classify_nucleic_moltypes();

}

}  // namespace sasmol
