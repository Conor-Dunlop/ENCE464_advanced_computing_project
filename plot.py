import numpy as np
import matplotlib.pyplot as plt
import sys
import os

    
def threads():
    
    thread_vals = np.array([1, 2, 3, 4, 6, 8, 10, 12, 14, 16])
    times = np.empty(thread_vals.size)

    call_path = os.path.relpath('./build/poisson')

    for i in range(thread_vals.size):
        res = os.popen(f"{call_path} -n 301 -i 100 -t {thread_vals[i]} -r 10").read()
        time = float(res.splitlines()[-1].split()[-1])
        times[i] = time
        print(f"Completed run with {thread_vals[i]} threads")

    fig, ax = plt.subplots()
    ax.plot(thread_vals, times)
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Execution time (s)")
    ax.grid()
    plt.ylim(bottom=0)

    plt.show()

def time_exec(name, size, runs):
    call_path = os.path.relpath(f'./build/{name}')
    res = os.popen(f"{call_path} -n {size} -i 100 -t 8 -r {runs}").read()
    time = float(res.splitlines()[-1].split()[-1])
    print(f"Completed {name} run with size {size} ")
    return time

def simd():
    
    sizes = np.array([101, 201, 301, 401, 501, 601, 701, 801, 901])
    runs = np.array([10, 10, 10, 10, 8, 6, 4, 2, 1])

    times_scalar = np.empty(sizes.size)
    times_avx2 = np.empty(sizes.size)
    times_avx512 = np.empty(sizes.size)

    names = ['poisson_scalar', 'poisson_avx2', 'poisson_avx512']
    times = {names[0]: np.empty(sizes.size), names[1]: np.empty(sizes.size), names[2]: np.empty(sizes.size)}

    for name in names:
        for i in range(sizes.size):
            times[name][i] = time_exec(name, sizes[i], runs[i])

    fig, ax = plt.subplots()
    ax.plot(sizes, times[names[0]], label='Scalar')
    ax.plot(sizes, times[names[1]], label='AVX2')
    ax.plot(sizes, times[names[2]], label='AVX512')
    ax.set_xlabel("Size (N)")
    ax.set_ylabel("Execution time (s)")
    ax.grid()
    plt.ylim(bottom=0)

    plt.show()

def main():
    #threads()
    simd()

if __name__ == "__main__":
    sys.exit(int(main() or 0))