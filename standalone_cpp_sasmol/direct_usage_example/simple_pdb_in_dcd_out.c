#include "sasmol_direct_usage.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: simple_pdb_in_dcd_out_c input.pdb output.dcd\n");
    return EXIT_FAILURE;
  }

  char error_message[512];
  int atom_count = 0;
  int frame_count = 0;
  const int ok = sasmol_pdb_to_dcd(argv[1], argv[2], error_message,
                                   sizeof(error_message), &atom_count,
                                   &frame_count);
  if (!ok) {
    fprintf(stderr, "sasmol_pdb_to_dcd failed: %s\n", error_message);
    return EXIT_FAILURE;
  }

  printf("wrote %s from %s atoms=%d frames=%d\n", argv[2], argv[1],
         atom_count, frame_count);
  return EXIT_SUCCESS;
}
