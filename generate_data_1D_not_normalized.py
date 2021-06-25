import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_not_normalized import calc_g_av

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    t = 2
    tau_step = 0.0001

    taus = list(np.arange(0, 0.01, tau_step)) # List of all "tau"s
    results = []

    for tau in taus:
        print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'), end="\r")
        _, _, _, final_result = calc_g_av(t, tau, T)
        results.append(final_result)
    plt.plot(taus, results, ".-")
    plt.title(r'$t\ =\ ' + str(t) + '$')
    plt.ylabel(r'$g_{av}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    # plt.show()
    filename = 'T=' + str(T) + '_t=' + str(t) + '.pdf'
    plt.savefig(filename)

    # filename = 'data_T=' + str(T) + '.pkl'
    # with open(filename, 'wb') as file:
    #     pickle.dump(results, file)

    print("t: ", format(t, '.5f'), " || tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'))

if __name__ == "__main__":
    main()
