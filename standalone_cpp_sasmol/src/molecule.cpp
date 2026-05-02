#include "sasmol/molecule.hpp"

#include <stdexcept>

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

}  // namespace sasmol
