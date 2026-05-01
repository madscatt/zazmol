#pragma once

#include "sasmol/types.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace sasmol {

struct IntegrityIssue {
  std::string field;
  std::size_t expected{};
  std::size_t actual{};
};

struct IntegrityReport {
  std::map<std::string, std::size_t> lengths;
  std::vector<IntegrityIssue> issues;

  [[nodiscard]] bool ok() const noexcept { return issues.empty(); }
};

class Molecule {
 public:
  Molecule() = default;
  explicit Molecule(std::size_t natoms, std::size_t number_of_frames = 1);

  void resize(std::size_t natoms, std::size_t number_of_frames = 1);

  [[nodiscard]] std::size_t natoms() const noexcept { return natoms_; }
  [[nodiscard]] std::size_t number_of_frames() const noexcept {
    return number_of_frames_;
  }

  [[nodiscard]] std::vector<std::string>& record() noexcept { return record_; }
  [[nodiscard]] const std::vector<std::string>& record() const noexcept {
    return record_;
  }
  [[nodiscard]] std::vector<int>& index() noexcept { return index_; }
  [[nodiscard]] const std::vector<int>& index() const noexcept {
    return index_;
  }
  [[nodiscard]] std::vector<int>& original_index() noexcept {
    return original_index_;
  }
  [[nodiscard]] const std::vector<int>& original_index() const noexcept {
    return original_index_;
  }
  [[nodiscard]] std::vector<int>& original_resid() noexcept {
    return original_resid_;
  }
  [[nodiscard]] const std::vector<int>& original_resid() const noexcept {
    return original_resid_;
  }
  [[nodiscard]] std::vector<std::string>& name() noexcept { return name_; }
  [[nodiscard]] const std::vector<std::string>& name() const noexcept {
    return name_;
  }
  [[nodiscard]] std::vector<std::string>& loc() noexcept { return loc_; }
  [[nodiscard]] const std::vector<std::string>& loc() const noexcept {
    return loc_;
  }
  [[nodiscard]] std::vector<std::string>& resname() noexcept {
    return resname_;
  }
  [[nodiscard]] const std::vector<std::string>& resname() const noexcept {
    return resname_;
  }
  [[nodiscard]] std::vector<std::string>& chain() noexcept { return chain_; }
  [[nodiscard]] const std::vector<std::string>& chain() const noexcept {
    return chain_;
  }
  [[nodiscard]] std::vector<int>& resid() noexcept { return resid_; }
  [[nodiscard]] const std::vector<int>& resid() const noexcept {
    return resid_;
  }
  [[nodiscard]] std::vector<std::string>& rescode() noexcept {
    return rescode_;
  }
  [[nodiscard]] const std::vector<std::string>& rescode() const noexcept {
    return rescode_;
  }
  [[nodiscard]] std::vector<std::string>& occupancy() noexcept {
    return occupancy_;
  }
  [[nodiscard]] const std::vector<std::string>& occupancy() const noexcept {
    return occupancy_;
  }
  [[nodiscard]] std::vector<std::string>& beta() noexcept { return beta_; }
  [[nodiscard]] const std::vector<std::string>& beta() const noexcept {
    return beta_;
  }
  [[nodiscard]] std::vector<std::string>& segname() noexcept {
    return segname_;
  }
  [[nodiscard]] const std::vector<std::string>& segname() const noexcept {
    return segname_;
  }
  [[nodiscard]] std::vector<std::string>& element() noexcept {
    return element_;
  }
  [[nodiscard]] const std::vector<std::string>& element() const noexcept {
    return element_;
  }
  [[nodiscard]] std::vector<std::string>& charge() noexcept { return charge_; }
  [[nodiscard]] const std::vector<std::string>& charge() const noexcept {
    return charge_;
  }
  [[nodiscard]] std::vector<calc_type>& atom_charge() noexcept {
    return atom_charge_;
  }
  [[nodiscard]] const std::vector<calc_type>& atom_charge() const noexcept {
    return atom_charge_;
  }
  [[nodiscard]] std::vector<calc_type>& atom_vdw() noexcept {
    return atom_vdw_;
  }
  [[nodiscard]] const std::vector<calc_type>& atom_vdw() const noexcept {
    return atom_vdw_;
  }
  [[nodiscard]] std::vector<std::string>& moltype() noexcept {
    return moltype_;
  }
  [[nodiscard]] const std::vector<std::string>& moltype() const noexcept {
    return moltype_;
  }
  [[nodiscard]] std::vector<calc_type>& mass() noexcept { return mass_; }
  [[nodiscard]] const std::vector<calc_type>& mass() const noexcept {
    return mass_;
  }
  [[nodiscard]] calc_type total_mass() const noexcept { return total_mass_; }
  void set_total_mass(calc_type value) noexcept { total_mass_ = value; }

  [[nodiscard]] std::map<std::string, std::size_t>& formula() noexcept {
    return formula_;
  }
  [[nodiscard]] const std::map<std::string, std::size_t>& formula()
      const noexcept {
    return formula_;
  }
  [[nodiscard]] std::vector<calc_type>& residue_charge() noexcept {
    return residue_charge_;
  }
  [[nodiscard]] const std::vector<calc_type>& residue_charge() const noexcept {
    return residue_charge_;
  }
  [[nodiscard]] std::string& fasta() noexcept { return fasta_; }
  [[nodiscard]] const std::string& fasta() const noexcept { return fasta_; }
  [[nodiscard]] std::array<calc_type, 6>& unitcell() noexcept {
    return unitcell_;
  }
  [[nodiscard]] const std::array<calc_type, 6>& unitcell() const noexcept {
    return unitcell_;
  }
  [[nodiscard]] std::vector<std::vector<int>>& conect() noexcept {
    return conect_;
  }
  [[nodiscard]] const std::vector<std::vector<int>>& conect() const noexcept {
    return conect_;
  }

  [[nodiscard]] std::vector<coord_type>& coor() noexcept { return coor_; }
  [[nodiscard]] const std::vector<coord_type>& coor() const noexcept {
    return coor_;
  }

  [[nodiscard]] Vec3 coordinate(std::size_t frame, std::size_t atom) const;
  void set_coordinate(std::size_t frame, std::size_t atom, Vec3 value);
  [[nodiscard]] CoordinateView coordinate_view(std::size_t frame);
  [[nodiscard]] ConstCoordinateView coordinate_view(std::size_t frame) const;

  [[nodiscard]] IntegrityReport check_integrity(bool fast_check = false) const;

 private:
  [[nodiscard]] std::size_t offset(std::size_t frame,
                                   std::size_t atom) const;

  std::size_t natoms_{};
  std::size_t number_of_frames_{};
  std::vector<coord_type> coor_;

  std::vector<std::string> record_;
  std::vector<int> index_;
  std::vector<int> original_index_;
  std::vector<int> original_resid_;
  std::vector<std::string> name_;
  std::vector<std::string> loc_;
  std::vector<std::string> resname_;
  std::vector<std::string> chain_;
  std::vector<int> resid_;
  std::vector<std::string> rescode_;
  std::vector<std::string> occupancy_;
  std::vector<std::string> beta_;
  std::vector<std::string> segname_;
  std::vector<std::string> element_;
  std::vector<std::string> charge_;
  std::vector<calc_type> atom_charge_;
  std::vector<calc_type> atom_vdw_;
  std::vector<std::string> moltype_;
  std::vector<calc_type> mass_;

  calc_type total_mass_{};
  std::map<std::string, std::size_t> formula_;
  std::vector<calc_type> residue_charge_;
  std::string fasta_;
  std::array<calc_type, 6> unitcell_{};
  std::vector<std::vector<int>> conect_;
};

}  // namespace sasmol
