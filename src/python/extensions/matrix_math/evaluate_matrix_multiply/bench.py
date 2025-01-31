import numpy as np
import time
import matrix_math
import matplotlib.pyplot as plt
import ctypes
from sasmol.linear_algebra import matrix_multiply as linear_algebra_matrix_multiply  # Import the C extension method

# Load the shared library
lib = ctypes.CDLL('./matrix_multiply.so')

# Define the argument and return types for the C function
lib.matrix_multiply_c.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                  ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float),
                                  ctypes.POINTER(ctypes.c_float)]

# Function to benchmark and compare results
def benchmark_and_compare(size):
    print(f"Benchmarking for size: {size}x{size}")

    # Create random matrices of the given size
    a = np.random.rand(size, size).astype(np.float32)
    b = np.random.rand(size, size).astype(np.float32)

    # Benchmark the first implementation (NumPy's dot)
    start_time = time.time()
    result1 = np.dot(a, b)
    end_time = time.time()
    time1 = end_time - start_time
    print(f"NumPy dot time: {time1:.6f} seconds")

    # Benchmark the second implementation (C extension)
    start_time = time.time()
    result2 = matrix_math.matrix_multiply(a, b)
    end_time = time.time()
    time2 = end_time - start_time
    print(f"C extension time: {time2:.6f} seconds")

    # Benchmark the third implementation (C function via ctypes)
    result3 = np.zeros((size, size), dtype=np.float32)
    start_time = time.time()
    lib.matrix_multiply_c(size, size, size,
                          a.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                          b.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                          result3.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
    end_time = time.time()
    time3 = end_time - start_time
    print(f"C function via ctypes time: {time3:.6f} seconds")

    # Benchmark the fourth implementation (linear_algebra matrix_multiply)
    start_time = time.time()
    error, result4 = linear_algebra_matrix_multiply(a, b)
    end_time = time.time()
    time4 = end_time - start_time
    print(f"linear_algebra matrix_multiply time: {time4:.6f} seconds")

    # Compare the numerical results
    if np.allclose(result1, result2, atol=1e-6) and np.allclose(result1, result3, atol=1e-6) and np.allclose(result1, result4, atol=1e-6):
        print("Results are numerically close.")
    else:
        print("Results differ.")

    # Calculate and compare the trace of the matrices
    trace1 = np.trace(result1)
    trace2 = np.trace(result2)
    trace3 = np.trace(result3)
    trace4 = np.trace(result4)
    print(f"Trace comparison: NumPy={trace1}, C extension={trace2}, C function via ctypes={trace3}, linear_algebra={trace4}")

    # Calculate and compare the Frobenius norm of the matrices
    norm1 = np.linalg.norm(result1)
    norm2 = np.linalg.norm(result2)
    norm3 = np.linalg.norm(result3)
    norm4 = np.linalg.norm(result4)
    print(f"Frobenius norm comparison: NumPy={norm1}, C extension={norm2}, C function via ctypes={norm3}, linear_algebra={norm4}")

    # Plot the numerical differences
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle(f'Matrix Multiplication Results for {size}x{size}')
    im = axs[0].imshow(result1, cmap='viridis')
    axs[0].set_title('NumPy dot result')
    fig.colorbar(im, ax=axs[0])
    im = axs[1].imshow(result2, cmap='viridis')
    axs[1].set_title('C extension result')
    fig.colorbar(im, ax=axs[1])
    im = axs[2].imshow(result3, cmap='viridis')
    axs[2].set_title('C function via ctypes result')
    fig.colorbar(im, ax=axs[2])
    im = axs[3].imshow(result4, cmap='viridis')
    axs[3].set_title('linear_algebra result')
    fig.colorbar(im, ax=axs[3])
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'matrix_results_{size}x{size}.png')
    plt.close()

    # Comment out the parts that plot the differences
    # fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    # fig.suptitle(f'Differences in Matrix Multiplication Results for {size}x{size}')
    # im = axs[0].imshow(result1 - result2, cmap='viridis')
    # axs[0].set_title('Difference (NumPy - C extension)')
    # fig.colorbar(im, ax=axs[0])
    # im = axs[1].imshow(result1 - result3, cmap='viridis')
    # axs[1].set_title('Difference (NumPy - C function via ctypes)')
    # fig.colorbar(im, ax=axs[1])
    # im = axs[2].imshow(result1 - result4, cmap='viridis')
    # axs[2].set_title('Difference (NumPy - linear_algebra)')
    # fig.colorbar(im, ax=axs[2])
    # plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.show()

    print("NumPy dot result:\n", result1)
    print("C extension result:\n", result2)
    print("C function via ctypes result:\n", result3)
    print("linear_algebra matrix_multiply result:\n", result4)

    return time1, time2, time3, time4, trace1, trace2, trace3, trace4, norm1, norm2, norm3, norm4

# List of matrix sizes to benchmark
sizes = [3, 10, 100, 500, 1000]

# Lists to store timing results
numpy_times = []
c_extension_times = []
c_function_times = []
linear_algebra_times = []

# Lists to store trace and norm results
trace_results = []
norm_results = []

# Run the benchmark for each size
for size in sizes:
    time1, time2, time3, time4, trace1, trace2, trace3, trace4, norm1, norm2, norm3, norm4 = benchmark_and_compare(size)
    numpy_times.append(time1)
    c_extension_times.append(time2)
    c_function_times.append(time3)
    linear_algebra_times.append(time4)
    trace_results.append((trace1, trace2, trace3, trace4))
    norm_results.append((norm1, norm2, norm3, norm4))

# Plot the timing results
plt.figure(figsize=(10, 6))
plt.plot(sizes, numpy_times, label='NumPy dot', marker='o')
plt.plot(sizes, c_extension_times, label='C extension', marker='o')
plt.plot(sizes, c_function_times, label='C function via ctypes', marker='o')
plt.plot(sizes, linear_algebra_times, label='linear_algebra matrix_multiply', marker='o')
plt.xlabel('Matrix Size (NxN)')
plt.ylabel('Time (seconds)')
plt.title('Matrix Multiplication Benchmark')
plt.legend()
plt.grid(True)
plt.savefig('timings_mm.png')
plt.close()

# Summarize trace and norm comparisons
print("\nTrace Comparison Summary:")
for i, size in enumerate(sizes):
    trace1, trace2, trace3, trace4 = trace_results[i]
    print(f"Size {size}x{size}: NumPy={trace1}, C extension={trace2}, C function via ctypes={trace3}, linear_algebra={trace4}")

print("\nFrobenius Norm Comparison Summary:")
for i, size in enumerate(sizes):
    norm1, norm2, norm3, norm4 = norm_results[i]
    print(f"Size {size}x{size}: NumPy={norm1}, C extension={norm2}, C function via ctypes={norm3}, linear_algebra={norm4}")