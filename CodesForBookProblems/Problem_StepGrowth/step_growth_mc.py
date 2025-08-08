import random
import numpy as np
import matplotlib.pyplot as plt

def run_one_trial(N0, record_every=1, seed=None):
    """Run one Monte Carlo trial of linear step-growth polymerization.
    Returns arrays of (r_list, p_list, Mn_list, Mw_list) recorded after every `record_every` reactions.
    """
    if seed is not None:
        random.seed(seed)
    chains = [1]*N0  # chain lengths (number of monomers)
    total_monomers = N0
    r = 0
    r_list = []
    p_list = []
    Mn_list = []
    Mw_list = []
    max_reactions = N0 - 1
    while r < max_reactions:
        # pick two distinct molecules uniformly at random
        i, j = random.sample(range(len(chains)), 2)
        Li = chains[i]
        Lj = chains[j]
        # create new chain
        newL = Li + Lj
        # remove the two (remove by index; ensure j>i)
        if i > j:
            i, j = j, i
        # Remove j first, then i
        chains.pop(j)
        chains.pop(i)
        chains.append(newL)
        r += 1
        if r % record_every == 0 or r == max_reactions:
            N = len(chains)
            lengths = np.array(chains, dtype=float)
            Mn = total_monomers / N  # in monomer units
            Mw = np.sum(lengths**2) / np.sum(lengths)
            p = r / N0
            r_list.append(r)
            p_list.append(p)
            Mn_list.append(Mn)
            Mw_list.append(Mw)
    return np.array(r_list), np.array(p_list), np.array(Mn_list), np.array(Mw_list)

def ensemble_average(N0, trials=200, record_every=1):
    # we will record at each r from 1..N0-1; to align averaging we will store by reaction count index
    max_r = N0 - 1
    Mn_sum = np.zeros(max_r)
    Mw_sum = np.zeros(max_r)
    for t in range(trials):
        r_list, p_list, Mn_list, Mw_list = run_one_trial(N0, record_every=1)
        # r_list goes 1..max_r; store into sums (index r-1)
        for idx, r in enumerate(r_list):
            Mn_sum[r-1] += Mn_list[idx]
            Mw_sum[r-1] += Mw_list[idx]
    Mn_avg = Mn_sum / trials
    Mw_avg = Mw_sum / trials
    r_axis = np.arange(1, max_r+1)
    p_axis = r_axis / N0
    return r_axis, p_axis, Mn_avg, Mw_avg

if __name__ == "__main__":
    N0 = 1000       # starting number of monomers (adjust for speed / statistics)
    trials = 200    # ensemble average runs
    r_axis, p_axis, Mn_avg, Mw_avg = ensemble_average(N0, trials=trials)

    # Convert DP to molecular weights assuming M0 = 100 g/mol (example)
    M0 = 100.0
    Mn_avg_MW = Mn_avg * M0
    Mw_avg_MW = Mw_avg * M0

    # Downsample by taking every 'skip' points
    skip = 20
    p_axis_ds   = p_axis[::skip]
    Mn_avg_ds   = Mn_avg_MW[::skip]
    Mw_avg_ds   = Mw_avg_MW[::skip]

    # ---- Plot settings ----
    plt.rcParams.update({
        "font.size": 16,          # default text size
        "axes.labelsize": 18,     # x and y labels
        "xtick.labelsize": 16,    # x tick labels
        "ytick.labelsize": 16,    # y tick labels
        "legend.fontsize": 14,    # legend font
    })

    plt.figure(figsize=(6,4))
    plt.plot(p_axis_ds, Mn_avg_ds, 'ko', label=r'$M_n$ (MC average)')
    plt.plot(p_axis_ds, Mw_avg_ds, 'k^', label=r'$M_w$ (MC average)')
    # plot Flory predictions for comparison
    p_theory = p_axis
    Mn_theory = M0/(1 - p_theory)
    Mw_theory = M0*(1 + p_theory)/(1 - p_theory)
    plt.plot(p_theory, Mn_theory, 'k-', label=r'$M_n$ (Flory)')
    plt.plot(p_theory, Mw_theory, 'k--', label=r'$M_w$ (Flory)')
    plt.xlabel('Conversion p')
    plt.ylabel('Molecular weight (g/mol)')
    plt.legend()
#    plt.title('Step-growth polymerization: $M_n$ and $M_w$ vs conversion p')
    plt.tight_layout()
    plt.yscale('log')
#    plt.show()
    plt.savefig('MnMw_vs_p.png')

    # PDI vs p
    PDI_mc = Mw_avg / Mn_avg
    PDI_mc_ds = Mw_avg_ds / Mn_avg_ds
    PDI_theory = 1 + p_axis
    plt.figure(figsize=(6,4))
    plt.plot(p_axis_ds, PDI_mc_ds, 'ko', label='PDI (MC)')
    plt.plot(p_axis, PDI_theory, 'k-', label='PDI (Flory: 1+p)')
    plt.xlabel('Conversion p')
    plt.ylabel('PDI = $M_w/M_n$')
    plt.legend()
    plt.xlim(0,0.9)
    plt.tight_layout()
 #   plt.show()
    plt.savefig('PDI_vs_p.png')


