#include "sasmol/overlap.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>

namespace {

void test_vector_overlap_detects_close_points() {
  const std::vector<sasmol::Vec3> first{{0.0F, 0.0F, 0.0F},
                                        {10.0F, 0.0F, 0.0F}};
  const std::vector<sasmol::Vec3> second{{2.0F, 0.0F, 0.0F},
                                         {10.4F, 0.0F, 0.0F}};

  assert(sasmol::has_overlap(first, second, 0.5));
  assert(!sasmol::has_overlap(first, second, 0.3));
}

void test_molecule_overlap_uses_requested_frames() {
  sasmol::Molecule first(1, 2);
  sasmol::Molecule second(1, 2);
  first.set_coordinate(0, 0, {0.0F, 0.0F, 0.0F});
  first.set_coordinate(1, 0, {100.0F, 0.0F, 0.0F});
  second.set_coordinate(0, 0, {100.0F, 0.0F, 0.0F});
  second.set_coordinate(1, 0, {0.25F, 0.0F, 0.0F});

  assert(sasmol::has_overlap(first, 0, second, 1, 0.5));
  assert(!sasmol::has_overlap(first, 1, second, 1, 0.5));
}

void test_negative_cutoff_is_rejected() {
  bool threw = false;
  try {
    (void)sasmol::has_overlap(std::vector<sasmol::Vec3>{{0.0F, 0.0F, 0.0F}},
                              std::vector<sasmol::Vec3>{}, -1.0);
  } catch (const std::invalid_argument&) {
    threw = true;
  }

  assert(threw);
}

}  // namespace

int main() {
  test_vector_overlap_detects_close_points();
  test_molecule_overlap_uses_requested_frames();
  test_negative_cutoff_is_rejected();
  return 0;
}
