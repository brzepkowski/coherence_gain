import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from numpy import absolute
import coherence_gain_t_finite
import coherence_gain_t_inf

def main():
    if len(sys.argv) != 3:
        print("Wrong number of parameters! Please provide them in the following fashion:")
        print("- temperature T [K],")
        print("- time of measurement t [ps] (it's the upper bound).")
        sys.exit()

    T = float(sys.argv[1])
    t = float(sys.argv[2])

    t_step = 1

    ts = list(np.arange(0, t, t_step)) # List of all "t"s
    results_D_t_finite = []
    results_D_t_inf = []

    _, _, D_t_inf, _, _, _ = coherence_gain_t_inf.calc_g_av(0, T) # Here tau doesn't play any role, because we are extracting only D(t)

    for t_intermediate in ts:
        print("t: ", format(t_intermediate, '.5f'), " / ", format(ts[-1], '.5f'))
        _, _, _, _, _, _, _, D_t_finite, _, _, _ = coherence_gain_t_finite.calc_g_av(t_intermediate, 0, T) # Here tau doesn't play any role, because we are extracting only D_t.

        results_D_t_finite.append(D_t_finite)
        results_D_t_inf.append(D_t_inf)

    plt.plot(ts, results_D_t_finite, ".-", label=r'$D_t$')
    plt.plot(ts, results_D_t_inf, ".-", label=r'$D_t^\infty$')
    plt.title(r'$t\ finite\ vs.\ infinite$')
    plt.xlabel(r'$t$')
    plt.grid()
    plt.legend()
    filename = 'T=' + str(T) + '_t_finite_vs_infinite.pdf'
    plt.savefig(filename)
    plt.show()

if __name__ == "__main__":
    main()
