#include <stdio.h>
#include <stdlib.h>

#define DCD_BADFORMAT -6

int reverseFourByteWord(int input) {
    return ((input >> 24) & 0x000000FF) |
           ((input >> 8)  & 0x0000FF00) |
           ((input << 8)  & 0x00FF0000) |
           ((input << 24) & 0xFF000000);
}

int read_dcdheader(FILE *fd) {
    int input_integer;
    size_t ret_val;

    // Ensure the file pointer is at the beginning of the file
    fseek(fd, 0, SEEK_SET);

    // Read the magic number from the file header
    ret_val = fread(&input_integer, sizeof(int), 1, fd);
    if (ret_val != 1) {
        fprintf(stderr, "Error reading first int from DCD file\n");
        fflush(stderr);
        return DCD_BADFORMAT;
    }
    fprintf(stderr, "read_dcdheader: input_integer = %d\n", input_integer);
    fflush(stderr);

    /* Check magic number in file header and determine byte order*/
    if (input_integer != 84) {
        // Reverse the byte order
        input_integer = reverseFourByteWord(input_integer);

        if (input_integer != 84) {
            fprintf(stderr, "Invalid magic number in file header: %d\n", input_integer);
            fflush(stderr);
            return DCD_BADFORMAT;
        }
    }

    fprintf(stderr, "Valid DCD file header\n");
    fflush(stderr);
    return 0;
}

int main() {
    FILE *file = fopen("h200.dcd", "rb");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }

    int result = read_dcdheader(file);
    if (result != 0) {
        fprintf(stderr, "Error reading DCD header: Failed to read DCD header\n");
    } else {
        fprintf(stderr, "DCD header read successfully\n");
    }

    fclose(file);
    return 0;
}