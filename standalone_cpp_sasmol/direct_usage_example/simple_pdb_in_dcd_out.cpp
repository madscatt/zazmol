#include "sasmol/file_io.hpp"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: simple_pdb_in_dcd_out_cpp input.pdb output.dcd\n";
    return EXIT_FAILURE;
  }

  sasmol::PdbReader pdb_reader;
  sasmol::Molecule molecule;
  auto status = pdb_reader.read_pdb(argv[1], molecule);
  if (!status) {
    std::cerr << "read_pdb failed: " << status.message << '\n';
    return EXIT_FAILURE;
  }

  sasmol::DcdWriter dcd_writer;
  status = dcd_writer.write_dcd(argv[2], molecule);
  if (!status) {
    std::cerr << "write_dcd failed: " << status.message << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "wrote " << argv[2] << " from " << argv[1]
            << " atoms=" << molecule.natoms()
            << " frames=" << molecule.number_of_frames() << '\n';
  return EXIT_SUCCESS;
}
