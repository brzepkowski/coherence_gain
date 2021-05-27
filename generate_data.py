import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_no_normalization import calc_g_av

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    t_limit = 0.005
    t_step = 0.001
    tau_step = 0.001

    ts = np.arange(t_step, t_limit + t_step, t_step) # List of all "t"s

    results = []
    for t in ts:
        results_t = []
        taus = list(np.arange(0, t, tau_step)) # List of all "tau"s
        for tau in taus:
            print("t: ", format(t, '.3f'), " / ", format(ts[-1], '.3f'), " || tau: ", format(tau, '.3f'), " / ", format(taus[-1], '.3f'), end="\r")
            g_av = calc_g_av(t, tau, T)
            results_t.append(g_av)
        results.append([t, taus, results_t])
        # plt.plot(taus, results_t, ".-")
        # plt.title(r'$t\ =\ ' + str(t) + '$')
        # plt.ylabel(r'$g_{av}$')
        # plt.xlabel(r'$\tau$')
        # plt.show()

    filename = 'data_T=' + str(T) + '.pkl'
    with open(filename, 'wb') as file:
        pickle.dump(results, file)

if __name__ == "__main__":
    main()
