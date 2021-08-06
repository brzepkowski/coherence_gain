import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from coherence_gain_t_inf import calc_g_av

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide only the temperature T [K] as an argument.")
        sys.exit()

    T = float(sys.argv[1])

    tau_step = 0.0001

    taus = list(np.arange(6, 6.005, tau_step)) # List of all "tau"s
    results_D_plus = []
    results_D_minus = []
    results_D_t = []
    results_p_plus = []
    results_p_minus = []
    results_g_av = []

    for tau in taus:
        print("tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'), end="\r")
        D_plus, D_minus, D_t, p_plus, p_minus, g_av = calc_g_av(tau, T)
        results_D_plus.append(D_plus)
        results_D_minus.append(D_minus)
        results_D_t.append(D_t)

        results_p_plus.append(p_plus)
        results_p_minus.append(p_minus)

        results_g_av.append(g_av)
    print("tau: ", format(tau, '.5f'), " / ", format(taus[-1], '.5f'))

    plt.plot(taus, results_D_plus, ".-", label=r'$D_+^{\infty}$')
    plt.plot(taus, results_D_minus, ".-", label=r'$D_-^{\infty}$')
    plt.plot(taus, results_D_t, ".-", label=r'$D_t^{\infty}$')
    plt.title(r'$t\ \rightarrow \infty$')
    plt.ylabel(r'$D^{\infty}$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'Ds_T=' + str(T) + '.pdf'
    plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    plt.plot(taus, results_p_plus, ".-", label=r'$p_+$')
    plt.plot(taus, results_p_minus, ".-", label=r'$p_-$')
    plt.title(r'$t\ \rightarrow \infty$')
    plt.ylabel(r'$p$')
    plt.xlabel(r'$\tau$')
    plt.legend()
    plt.grid()
    filename = 'probabilities_T=' + str(T) + '.pdf'
    plt.savefig(filename)
    plt.show()

    # Clear the plot
    plt.clf()

    plt.plot(taus, results_g_av, ".-")
    plt.title(r'$t\ \rightarrow \infty$')
    plt.ylabel(r'$g_{av}^{\infty}$')
    plt.xlabel(r'$\tau$')
    plt.grid()
    filename = 'g_av_T=' + str(T) + '.pdf'
    plt.savefig(filename)
    plt.show()

if __name__ == "__main__":
    main()
