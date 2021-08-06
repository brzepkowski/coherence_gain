import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
from numpy import absolute
from coherence_gain_components import calc_W_one_part

def main():
    if len(sys.argv) != 2:
        print("Wrong number of parameters! Please provide them in the following fashion:")
        print("- time of measurement t [ps] (it's the upper bound).")
        sys.exit()

    # T = float(sys.argv[1])
    Ts = [0.00001, 10, 100, 300]
    t = float(sys.argv[1])

    t_step = 0.5

    ts = list(np.arange(0, t, t_step)) # List of all "t"s
    results_W = [[] for _ in ts]

    for T in Ts:
        print("#### T: ", T, " ####")
        index = 0
        for t_intermediate in ts:
            print("t: ", format(t_intermediate, '.5f'), " / ", format(ts[-1], '.5f'))
            W_basic = calc_W_one_part(t_intermediate, T) # Here tau doesn't play any role, because we are extracting only W_basic
            results_W[index].append(absolute(W_basic))
            # print("absolute(W_basic): ", absolute(W_basic))
            index += 1

    plt.plot(ts, results_W, ".-")
    plt.legend(["0 K","10 K","100 K","300 K"], loc=1)
    # plt.title(r'$t\ finite\ vs.\ infinite$')
    plt.xlabel(r'$t$')
    plt.grid()
    filename = 't_finite_multiple_temps.pdf'
    plt.savefig(filename)
    plt.show()

if __name__ == "__main__":
    main()
