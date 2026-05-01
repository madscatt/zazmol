#include "sasmol/file_io.hpp"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: pdb_roundtrip_writer input.pdb output.pdb\n";
    return EXIT_FAILURE;
  }

  sasmol::PdbReader reader;
  sasmol::Molecule molecule;
  auto status = reader.read_pdb(argv[1], molecule);
  if (!status) {
    std::cerr << status.message << '\n';
    return EXIT_FAILURE;
  }

  sasmol::PdbWriter writer;
  sasmol::PdbWriteOptions options;
  options.write_all_frames = molecule.number_of_frames() > 1;
  status = writer.write_pdb(argv[2], molecule, options);
  if (!status) {
    std::cerr << status.message << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "wrote " << argv[2] << " frames=" << molecule.number_of_frames()
            << " atoms=" << molecule.natoms() << '\n';
  return EXIT_SUCCESS;
}
