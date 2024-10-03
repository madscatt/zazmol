#include <stdio.h>
#include <stdlib.h>

#define DCD_BADFORMAT -1

int read_dcdstep(FILE * fd, int N, float *X, float *Y, float *Z, int num_fixed,
        int first, int reverseEndian, int charmm) 
{
  int ret_val;          /*  Return value from read          */
  int input_integer;    /*  Integer buffer space            */
  int i;                /*  Loop counter              */

  // Debug: Print file pointer position before reading the frame
  long pos_before = ftell(fd);
  printf("File pointer position before reading frame: %ld\n", pos_before);

  if ((num_fixed == 0) || first) {
    // Read frame data
    ret_val = fread(&input_integer, sizeof(int), 1, fd);
    if (ret_val != 1) {
      printf("Error reading frame data\n");
      return DCD_BADFORMAT;
    }

    // Debug: Print file pointer position after reading the frame
    long pos_after = ftell(fd);
    printf("File pointer position after reading frame: %ld\n", pos_after);

    // Read X, Y, Z coordinates
    ret_val = fread(X, sizeof(float), N, fd);
    if (ret_val != N) {
      printf("Error reading X coordinates\n");
      return DCD_BADFORMAT;
    }
    ret_val = fread(Y, sizeof(float), N, fd);
    if (ret_val != N) {
      printf("Error reading Y coordinates\n");
      return DCD_BADFORMAT;
    }
    ret_val = fread(Z, sizeof(float), N, fd);
    if (ret_val != N) {
      printf("Error reading Z coordinates\n");
      return DCD_BADFORMAT;
    }
  }

  return 0;
}

int main() {
  FILE *fp = fopen("h.dcd", "rb");
  if (!fp) {
    printf("Failed to open file\n");
    return 1;
  }

  // Read header (example, adjust as needed)
  int header[10];
  fread(header, sizeof(int), 10, fp);

  // Debug: Print file pointer position after reading the header
  long pos_header = ftell(fp);
  printf("File pointer position after reading header: %ld\n", pos_header);

  // Rewind file pointer if necessary
  fseek(fp, pos_header, SEEK_SET);

  // Read first 5000 frames
  int N = 6730; // Example number of atoms
  float *X = (float *)malloc(N * sizeof(float));
  float *Y = (float *)malloc(N * sizeof(float));
  float *Z = (float *)malloc(N * sizeof(float));
  int num_fixed = 0;
  int first = 1;
  int reverseEndian = 0;
  int charmm = 0;

  for (int frame = 0; frame < 5000; frame++) {
    int result = read_dcdstep(fp, N, X, Y, Z, num_fixed, first, reverseEndian, charmm);
    if (result != 0) {
      printf("Failed to read DCD step at frame %d\n", frame);
      break;
    }
    first = 0; // Only the first frame is considered the first
  }

  fclose(fp);
  free(X);
  free(Y);
  free(Z);

  return 0;
}