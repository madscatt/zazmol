import dcdio_module as dcdio
import numpy as np

def main():
    # Open the DCD file
    try:
        file_capsule = dcdio.open_dcd_read("hiv1_gag_200_frames.dcd")
    except IOError as e:
        print(f"Failed to open DCD file: {e}")
        return

    # Read the DCD header
    try:
        result, file_capsule, nnatoms, nset, istart, nsavc, namnf, delta, reverseEndian, charmm = dcdio.read_dcdheader(file_capsule)
        if result != 0:
            print("Failed to read DCD header")
            return

        print("DCD Header Information:")
        print(f"Number of atoms: {nnatoms}")
        print(f"Number of sets: {nset}")
        print(f"Starting timestep: {istart}")
        print(f"Timesteps between saves: {nsavc}")
        print(f"Number of atoms with fixed coordinates: {namnf}")
        print(f"Time step size: {delta}")
        print(f"Reverse endian: {reverseEndian}")
        print(f"CHARMM version: {charmm}")

        # Determine and print the endianness
        if reverseEndian == 0:
            print("Endianness: Little-endian")
        else:
            print("Endianness: Big-endian")

        # Prepare arrays to hold the coordinates
        x_array = np.zeros(nnatoms, dtype=np.float32)
        y_array = np.zeros(nnatoms, dtype=np.float32)
        z_array = np.zeros(nnatoms, dtype=np.float32)

        # Read all 200 steps
        for step in range(nset):
            print('CALLING DCDIO.READ_DCDSTEP FOR STEP = ', step)
            #result = dcdio.read_dcdstep(file_capsule, nnatoms, namnf, step, reverseEndian, charmm, x_array, y_array, z_array)
            #result = dcdio.read_dcdstep(file_capsule, nnatoms, x_array, y_array, z_array, namnf,  reverseEndian, charmm)
            #result = read_dcdstep(fp, natoms, (float*)PyArray_DATA(x_array), (float*)PyArray_DATA(y_array), (float*)PyArray_DATA(z_array), num_fixed, first, reverseEndian, charmm)
            result = dcdio.read_dcdstep(file_capsule, nnatoms, x_array, y_array, z_array, namnf, step, reverseEndian, charmm)

            if result != 0:
                print(f"Failed to read DCD step {step}")
                return

            # Process the coordinates as needed
            print(f"Step {step + 1}:")
            print(f"x: {x_array}")
            print(f"y: {y_array}")
            print(f"z: {z_array}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
