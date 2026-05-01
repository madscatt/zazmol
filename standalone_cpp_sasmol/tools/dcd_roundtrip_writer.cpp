#include "sasmol/file_io.hpp"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: dcd_roundtrip_writer input.dcd output.dcd\n";
    return EXIT_FAILURE;
  }

  sasmol::DcdReader reader;
  sasmol::Molecule molecule;
  auto status = reader.read_dcd(argv[1], molecule);
  if (!status) {
    std::cerr << status.message << '\n';
    return EXIT_FAILURE;
  }

  sasmol::DcdWriter writer;
  status = writer.write_dcd(argv[2], molecule);
  if (!status) {
    std::cerr << status.message << '\n';
    return EXIT_FAILURE;
  }

  std::cout << "wrote " << argv[2] << " frames=" << molecule.number_of_frames()
            << " atoms=" << molecule.natoms() << '\n';
  return EXIT_SUCCESS;
}
