import dcdio_module as dcdio
import sys
import os

def main():
    # Redirect stderr to a file
    stderr_file = "stderr_output.txt"
    sys.stderr = open(stderr_file, "w")

    try:
        # Open the DCD file
        file_capsule = dcdio.open_dcd_read("h200.dcd")

        # Read the DCD header
        readheaderresult, pyfd, nnatoms, nset, istart, nsavc, namnf, delta, reverseEndian, charmm = dcdio.read_dcdheader(file_capsule)
        print("DCD Header Information:")
        print(f"Number of atoms: {nnatoms}")
        print(f"Number of sets: {nset}")
        print(f"Starting timestep: {istart}")
        print(f"Timesteps between saves: {nsavc}")
        print(f"Number of atoms with fixed coordinates: {namnf}")
        print(f"Time step size: {delta}")
        print(f"Reverse endian: {reverseEndian}")
        print(f"CHARMM version: {charmm}")

    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
        # Close the redirected stderr
        sys.stderr.close()
        sys.stderr = sys.__stderr__

        # Print the contents of the stderr file
        if os.path.exists(stderr_file):
            with open(stderr_file, "r") as f:
                stderr_contents = f.read()
                if stderr_contents:
                    print("Captured stderr output from C code:")
                    print(stderr_contents)
            os.remove(stderr_file)


        # Determine and print the endianness
        if reverseEndian == 0:
            print("Endianness: Little-endian")
        else:
            print("Endianness: Big-endian")


if __name__ == "__main__":
    main()
