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
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Execution time (s)")
    ax.plot(thread_vals, times)

    plt.show()

def main():
    threads()

if __name__ == "__main__":
    sys.exit(int(main() or 0))