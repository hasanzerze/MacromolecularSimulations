import random
import math
import numpy as np
import matplotlib.pyplot as plt

# Plot fonts for book-quality figures
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
})

def generate_saw(N_max, max_attempts=10000):
    """Attempt to generate a single SAW of length N_max (number of steps = N_max).
    Returns list of positions (x,y) of length N_max+1 on success, or None on failure.
    """
    attempts = 0
    while attempts < max_attempts:
        attempts += 1
        visited = set()
        x = y = 0
        visited.add((x,y))
        positions = [(x,y)]
        success = True
        for step in range(N_max):
            neighbors = [(x+1,y),(x-1,y),(x,y+1),(x,y-1)]
            free = [nb for nb in neighbors if nb not in visited]
            if not free:
                success = False
                break
            nx, ny = random.choice(free)
            x, y = nx, ny
            visited.add((x,y))
            positions.append((x,y))
        if success:
            return positions
    return None

def compute_mean_r2(N_max=200, n_trials=800, max_attempts_per_walk=3000):
    """Generate many SAWs and compute <r^2_n> for n=1..N_max."""
    sum_r2 = np.zeros(N_max+1, dtype=float)
    counts = np.zeros(N_max+1, dtype=int)
    trials_done = 0
    while trials_done < n_trials:
        walk = generate_saw(N_max, max_attempts=max_attempts_per_walk)
        if walk is None:
            continue
        trials_done += 1
        for n in range(1, N_max+1):
            x,y = walk[n]
            r2 = x*x + y*y
            sum_r2[n] += r2
            counts[n] += 1
        if trials_done % (max(1, n_trials//10)) == 0:
            print(f"Completed {trials_done}/{n_trials} SAW samples")
    n_vals = np.arange(1, N_max+1)
    mean_r2 = np.zeros(N_max, dtype=float)
    for n in n_vals:
        mean_r2[n-1] = sum_r2[n] / counts[n]
    return n_vals, mean_r2

if __name__ == "__main__":
    # Parameters (adjust for speed/accuracy)
    N_max = 200        # maximum contour length
    n_trials = 800     # number of successful SAWs to average
    max_attempts_per_walk = 3000

    n_vals, mean_r2 = compute_mean_r2(N_max=N_max, n_trials=n_trials, max_attempts_per_walk=max_attempts_per_walk)

    # Fit a power law <r^2> = A * n^alpha
    fit_min = 4
    fit_max = N_max
    log_n = np.log(n_vals[fit_min-1:fit_max])
    log_r2 = np.log(mean_r2[fit_min-1:fit_max])
    coeffs = np.polyfit(log_n, log_r2, 1)
    alpha = coeffs[0]
    A = np.exp(coeffs[1])
    print(f"Fitted power-law: <r^2> = {A:.3f} * n^{alpha:.3f} (fit range n={fit_min}..{fit_max})")

    # Log-log plot
    plt.figure(figsize=(7,5))
    plt.loglog(n_vals, mean_r2, 'ko', markersize=4, label=r'SARW (MC)')
    n_fine = np.linspace(1, N_max, 200)
    plt.loglog(n_fine, A * n_fine**alpha, 'k-', linewidth=2, label=f'Fit: n^{alpha:.2f}')
    # Plot ideal RW ~ n scaled to match at n=10 for visibility
    scale_point = 10
    scale = mean_r2[scale_point-1] / (scale_point)
    plt.loglog(n_fine, scale * n_fine, 'k:', linewidth=2, label='Ideal RW ~ n')
    plt.xlabel('Contour length n')
    plt.ylabel(r'$\langle r_n^2 \rangle$ (lattice units$^2$)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.tight_layout()
#    plt.show()
    plt.savefig("r2_vs_n.png")

    # Local exponent plot
    ln_n = np.log(n_vals)
    ln_r2 = np.log(mean_r2)
    local_exp = np.gradient(ln_r2, ln_n)
    plt.figure(figsize=(7,5))
    plt.plot(n_vals, local_exp, marker='.', markersize=4, linestyle='-')
    plt.axhline(y=1.0, linestyle=':', linewidth=1.5, label='RW exponent = 1')
    plt.axhline(y=1.5, linestyle='--', linewidth=1.5, label='SAW theory (2ν) = 1.5')
    plt.xlabel('Contour length n')
    plt.ylabel('Local exponent $d\ln\langle r^2\\rangle / d\ln n$')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

