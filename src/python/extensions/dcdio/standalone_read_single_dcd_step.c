#include <stdio.h>
#include <stdlib.h>
#include "dcdio.h"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <dcd file> <frame number>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    int frame = atoi(argv[2]);

    FILE *infile = open_dcd_read(filename);
    if (!infile) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        return 1;
    }

    int nnatoms, nsets, istart, nsavc, nfixed, reverseEndian, charmm;
    double delta;
    int *freeind = NULL;
    float *tx = NULL, *ty = NULL, *tz = NULL;

    if (read_dcdheader(infile, &nnatoms, &nsets, &istart, &nsavc, &delta, &nfixed, &reverseEndian, &charmm) != 0) {
        fprintf(stderr, "Failed to read DCD header\n");
        close_dcd_read(infile);
        return 1;
    }

    tx = (float *)malloc(nnatoms * sizeof(float));
    ty = (float *)malloc(nnatoms * sizeof(float));
    tz = (float *)malloc(nnatoms * sizeof(float));
    if (!tx || !ty || !tz) {
        fprintf(stderr, "Failed to allocate memory\n");
        free(tx);
        free(ty);
        free(tz);
        close_dcd_read(infile);
        return 1;
    }

    // Loop over frames until the desired frame is reached
    for (int i = 0; i < frame; i++) {
        if (read_dcdstep(infile, nnatoms, tx, ty, tz, 0, 0, reverseEndian, charmm) != 0) {
            fprintf(stderr, "Failed to read DCD step\n");
            free(tx);
            free(ty);
            free(tz);
            close_dcd_read(infile);
            return 1;
        }
    }

    // Read the desired frame
    if (read_dcdstep(infile, nnatoms, tx, ty, tz, 0, 0, reverseEndian, charmm) != 0) {
        fprintf(stderr, "Failed to read DCD step\n");
        free(tx);
        free(ty);
        free(tz);
        close_dcd_read(infile);
        return 1;
    }

    // Print the coordinates for verification
    for (int i = 0; i < nnatoms; i++) {
        printf("Atom %d: x = %f, y = %f, z = %f\n", i, tx[i], ty[i], tz[i]);
    }

    free(tx);
    free(ty);
    free(tz);
    close_dcd_read(infile);
    return 0;
}