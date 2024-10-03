import struct

DCD_BADFORMAT = -6

def reverse_four_byte_word(input_integer):
    return ((input_integer >> 24) & 0x000000FF) | \
           ((input_integer >> 8)  & 0x0000FF00) | \
           ((input_integer << 8)  & 0x00FF0000) | \
           ((input_integer << 24) & 0xFF000000)

def read_dcdheader(file_path):
    try:
        with open(file_path, 'rb') as fd:
            # Ensure the file pointer is at the beginning of the file
            fd.seek(0)

            # Read the magic number from the file header
            input_integer_bytes = fd.read(4)
            if len(input_integer_bytes) != 4:
                print("Error reading first int from DCD file")
                return DCD_BADFORMAT

            # Unpack the integer in little-endian format
            input_integer = struct.unpack('<i', input_integer_bytes)[0]
            print(f"read_dcdheader: input_integer = {input_integer}")

            # Check magic number in file header and determine byte order
            if input_integer != 84:
                # Reverse the byte order
                input_integer = reverse_four_byte_word(input_integer)

                if input_integer != 84:
                    print(f"Invalid magic number in file header: {input_integer}")
                    return DCD_BADFORMAT

            # Rewind the file pointer to the beginning of the file
            fd.seek(0)

            # Continue with the rest of the header reading process
            # (Add additional header reading logic here as needed)

            print("DCD header read successfully")
            return 0  # Return success code or other relevant information

    except IOError as e:
        print(f"Failed to open file: {e}")
        return DCD_BADFORMAT

if __name__ == "__main__":
    file_path = "h200.dcd"
    result = read_dcdheader(file_path)
    if result == DCD_BADFORMAT:
        print("Failed to read DCD header")
    else:
        print("DCD header read successfully")
